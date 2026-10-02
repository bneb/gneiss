use nalgebra::{Matrix3, SMatrix, SVector, Vector3};

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::coords::ecef_to_llh;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::signal::get_wavelength;
use gneiss_core::time::GpsTime;

/// Doppler-derived receiver velocity solution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DopplerVelocitySolution {
    /// 3D ECEF velocity in m/s
    pub vel_ecef: Vector3<f64>,
    /// Receiver clock drift in m/s (c * d(dt)/dt)
    pub clk_drift_m_s: f64,
    /// 3x3 covariance matrix of velocity in (m/s)^2
    pub cov: Matrix3<f64>,
    /// Number of satellites used in solution
    pub n_sats: usize,
    /// Velocity Dilution of Precision
    pub vdop: f64,
}

#[derive(Clone, Copy)]
struct SingleDopplerMeas {
    #[allow(dead_code)]
    sat: SatelliteId,
    e_los: Vector3<f64>,
    y_range_rate: f64,
    weight: f64,
}

fn compute_up_vector(pos: &Vector3<f64>) -> Vector3<f64> {
    let llh = ecef_to_llh(*pos);
    let cos_lat = llh.x.cos();
    let sin_lat = llh.x.sin();
    let cos_lon = llh.y.cos();
    let sin_lon = llh.y.sin();
    Vector3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)
}

fn extract_satellite_doppler(
    sat_obs: &gneiss_core::obs::SatObs,
    ephems: &[Ephemeris],
    t_rx: GpsTime,
    approx_pos: &Vector3<f64>,
    up: &Vector3<f64>,
) -> Option<SingleDopplerMeas> {
    let (doppler_hz, band) = sat_obs
        .get_doppler(1)
        .map(|d| (d, 1))
        .or_else(|| sat_obs.get_doppler(2).map(|d| (d, 2)))
        .or_else(|| sat_obs.get_doppler(5).map(|d| (d, 5)))?;

    // GLONASS (FDMA / leap-second time offset) and QZSS require dedicated handling; use robust CDMA constellations
    if sat_obs.sat.constellation == Constellation::Glonass || sat_obs.sat.constellation == Constellation::Qzss {
        return None;
    }

    let eph = crate::swfg::engine::epoch::select_best_ephemeris(ephems, sat_obs.sat, t_rx)?;
    let lambda = get_wavelength(sat_obs.sat, band, eph.freq_num());
    if lambda <= 0.0 {
        return None;
    }

    let (sat_p0, _, _, _) = eph.position(t_rx);
    let tau = (sat_p0 - approx_pos).norm() / SPEED_OF_LIGHT_M_S;
    let t_tx = t_rx - tau;
    let (sat_pos, _, _, sat_drift) = eph.position(t_tx);
    // Keplerian analytical velocity omits harmonic perturbations (~10-50 m/s error).
    // Compute satellite velocity via exact central difference of position.
    let dt = 0.05;
    let (p_plus, _, _, _) = eph.position(t_tx + dt);
    let (p_minus, _, _, _) = eph.position(t_tx - dt);
    let sat_vel = (p_plus - p_minus) / (2.0 * dt);

    // Rotate satellite position and velocity from transmission ECEF to reception ECEF frame (Sagnac effect)
    let theta = 7.292_115_146_7e-5 * tau;
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let sat_pos = Vector3::new(
        sat_pos.x * cos_t + sat_pos.y * sin_t,
        -sat_pos.x * sin_t + sat_pos.y * cos_t,
        sat_pos.z,
    );
    let sat_vel = Vector3::new(
        sat_vel.x * cos_t + sat_vel.y * sin_t,
        -sat_vel.x * sin_t + sat_vel.y * cos_t,
        sat_vel.z,
    );

    let los = sat_pos - approx_pos;
    let dist = los.norm();
    if dist < 1e6 {
        return None;
    }
    let e_los = los / dist;
    let sin_el = e_los.dot(up);
    if sin_el < 0.1736 {
        // Mask elevation < 10 degrees
        return None;
    }

    let c_sat_drift = SPEED_OF_LIGHT_M_S * sat_drift;
    let rho_meas = -lambda * doppler_hz;
    let rho_sat = e_los.dot(&sat_vel) - c_sat_drift;
    let y_range_rate = rho_meas - rho_sat;
    let weight = (sin_el * sin_el).clamp(0.05, 1.0);

    Some(SingleDopplerMeas {
        sat: sat_obs.sat,
        e_los,
        y_range_rate,
        weight,
    })
}

type Matrix6<T = f64> = SMatrix<T, 6, 6>;
type Vector6<T = f64> = SVector<T, 6>;

#[derive(Clone, Copy)]
struct ConstellationMap {
    gps: Option<usize>,
    gal: Option<usize>,
    bds: Option<usize>,
    n_params: usize,
}

impl ConstellationMap {
    fn new(meas: &[SingleDopplerMeas]) -> Self {
        let mut n = 3;
        let has_gps = meas.iter().any(|m| matches!(m.sat.constellation, Constellation::Gps | Constellation::Qzss));
        let has_gal = meas.iter().any(|m| m.sat.constellation == Constellation::Galileo);
        let has_bds = meas.iter().any(|m| m.sat.constellation == Constellation::Beidou);

        let gps = if has_gps { let col = n; n += 1; Some(col) } else { None };
        let gal = if has_gal { let col = n; n += 1; Some(col) } else { None };
        let bds = if has_bds { let col = n; n += 1; Some(col) } else { None };

        Self { gps, gal, bds, n_params: n }
    }

    fn col_for(&self, c: Constellation) -> Option<usize> {
        match c {
            Constellation::Gps | Constellation::Qzss => self.gps,
            Constellation::Galileo => self.gal,
            Constellation::Beidou => self.bds,
            _ => None,
        }
    }
}

fn solve_wls(meas: &[SingleDopplerMeas], cmap: &ConstellationMap) -> Option<(Vector6<f64>, Matrix6<f64>)> {
    let p = cmap.n_params;
    if meas.len() < p {
        return None;
    }
    let mut normal = Matrix6::zeros();
    let mut rhs = Vector6::zeros();

    for m in meas {
        let clk_col = cmap.col_for(m.sat.constellation)?;
        let mut h = Vector6::zeros();
        h[0] = -m.e_los.x;
        h[1] = -m.e_los.y;
        h[2] = -m.e_los.z;
        h[clk_col] = 1.0;

        let hw = h * m.weight;
        normal += hw * h.transpose();
        rhs += hw * m.y_range_rate;
    }

    for i in p..6 {
        normal[(i, i)] = 1.0;
    }

    let inv = normal.try_inverse()?;
    let sol = inv * rhs;
    Some((sol, inv))
}

fn compute_residuals(meas: &[SingleDopplerMeas], x: &Vector6<f64>, cmap: &ConstellationMap) -> Vec<f64> {
    meas.iter()
        .map(|m| {
            let clk_col = cmap.col_for(m.sat.constellation).unwrap_or(3);
            let pred = -m.e_los.x * x[0] - m.e_los.y * x[1] - m.e_los.z * x[2] + x[clk_col];
            m.y_range_rate - pred
        })
        .collect()
}

/// Nominal carrier phase Doppler standard deviation (m/s).
const NOMINAL_DOPPLER_SIGMA_M_S: f64 = 0.05;

fn filter_outliers(
    meas: &[SingleDopplerMeas],
    x: &Vector6<f64>,
    inv: &Matrix6<f64>,
    cmap: &ConstellationMap,
) -> Option<Vec<SingleDopplerMeas>> {
    let res = compute_residuals(meas, x, cmap);
    let mut max_idx = 0;
    let mut max_stat = 0.0;
    for (i, m) in meas.iter().enumerate() {
        let clk_col = match cmap.col_for(m.sat.constellation) {
            Some(c) => c,
            None => continue,
        };
        let mut h = Vector6::zeros();
        h[0] = -m.e_los.x;
        h[1] = -m.e_los.y;
        h[2] = -m.e_los.z;
        h[clk_col] = 1.0;

        let h_ii = m.weight * (h.transpose() * inv * h)[(0, 0)];
        let q_vv = ((1.0 - h_ii).max(0.01) / m.weight).max(1e-4);
        let sigma_r = NOMINAL_DOPPLER_SIGMA_M_S * q_vv.sqrt();
        let stat = res[i].abs() / sigma_r;
        if stat > max_stat {
            max_stat = stat;
            max_idx = i;
        }
    }
    if max_stat > 3.0 && meas.len() > cmap.n_params {
        let mut filtered = meas.to_vec();
        filtered.remove(max_idx);
        Some(filtered)
    } else {
        None
    }
}

/// Estimate instantaneous receiver 3D velocity and clock drift from multi-GNSS Doppler observations.
pub fn estimate_doppler_velocity(
    epoch: &EpochObs,
    ephems: &[Ephemeris],
    approx_pos: Vector3<f64>,
) -> Option<DopplerVelocitySolution> {
    let up = compute_up_vector(&approx_pos);
    let mut meas: Vec<SingleDopplerMeas> = epoch
        .satellites
        .iter()
        .filter_map(|s| extract_satellite_doppler(s, ephems, epoch.time, &approx_pos, &up))
        .collect();

    let mut cmap = ConstellationMap::new(&meas);
    if meas.len() < cmap.n_params {
        return None;
    }

    let (mut sol, mut cov6) = solve_wls(&meas, &cmap)?;
    if let Some(cleaned) = filter_outliers(&meas, &sol, &cov6, &cmap) {
        let c2map = ConstellationMap::new(&cleaned);
        if let Some((s2, c2)) = solve_wls(&cleaned, &c2map) {
            sol = s2;
            cov6 = c2;
            meas = cleaned;
            cmap = c2map;
        }
    }

    let res = compute_residuals(&meas, &sol, &cmap);
    let dof = (meas.len() as f64 - cmap.n_params as f64).max(1.0);
    let chi2: f64 = meas
        .iter()
        .zip(&res)
        .map(|(m, r)| m.weight * r * r)
        .sum::<f64>()
        / dof;
    let sigma_scale = chi2.clamp(0.01, 10.0);

    let mut cov3 = Matrix3::zeros();
    for r in 0..3 {
        for c in 0..3 {
            cov3[(r, c)] = cov6[(r, c)] * sigma_scale;
        }
    }
    // Floor minimal velocity variance to physically realistic measurement noise (2 cm/s)^2
    for i in 0..3 {
        if cov3[(i, i)] < 0.0004 {
            cov3[(i, i)] = 0.0004;
        }
    }

    let vdop = (cov6[(0, 0)] + cov6[(1, 1)] + cov6[(2, 2)]).sqrt();

    Some(DopplerVelocitySolution {
        vel_ecef: sol.fixed_rows::<3>(0).into_owned(),
        clk_drift_m_s: sol[3],
        cov: cov3,
        n_sats: meas.len(),
        vdop,
    })
}

#[cfg(test)]
mod tests;
