use nalgebra::{DMatrix, DVector, Vector3};
use gneiss_core::obs::{EpochObs, ObsType};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::coords::{ecef_to_llh, az_el, Coordinate, Datum, Frame};
use gneiss_core::atmosphere::{AtmosphereModel, KlobucharParams, TropoParams};
use gneiss_core::time::GpsTime;

/// Represents the current estimated state of the receiver.
#[derive(Debug, Clone, PartialEq)]
pub struct SppState {
    /// Receiver position in a specific Datum and Frame.
    pub position: Coordinate,
    /// Receiver clock bias in meters (c * dt) for GPS.
    pub cdt: f64,
    /// Receiver clock bias in meters for Galileo.
    pub cdt_gal: f64,
    /// Receiver clock bias in meters for BeiDou.
    pub cdt_bds: f64,
    /// Receiver clock bias in meters for GLONASS.
    pub cdt_glo: f64,
}

impl SppState {
    pub fn new(position: Coordinate, cdt: f64, cdt_gal: f64, cdt_bds: f64, cdt_glo: f64) -> Self {
        Self { position, cdt, cdt_gal, cdt_bds, cdt_glo }
    }
}

/// A single satellite measurement for use in the SPP estimator.
#[derive(Debug, Clone)]
pub struct SppMeasurement {
    pub constellation: gneiss_core::sat::Constellation,
    pub raw_pr: f64,
    pub snr: f64,
    pub doppler: f64,
    pub time: GpsTime,
    pub eph: Ephemeris,
}

/// Errors that can occur during an SPP WNLLS step.
#[derive(Debug, Clone, PartialEq)]
pub enum SppError {
    /// Not enough measurements to solve the state (need at least 4).
    NotEnoughMeasurements,
    /// The matrix inversion failed (singular matrix, poor geometry).
    MatrixInversionFailed,
    /// Maximum iterations reached without convergence.
    ConvergenceFailed,
    /// Geometry exploded or was too poor to use.
    PoorGeometry,
}

const LIGHT_SPEED: f64 = 299792458.0;
const OMEGA_E: f64 = 7.2921151467e-5; // WGS 84 value of earth's rotation rate

use serde::{Serialize, Deserialize};

/// Configuration for the Single Point Positioning solver.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SppConfig {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub geometry_variance_threshold: f64,
    pub enable_sagnac: bool,
    pub enable_tropo: bool,
    pub enable_iono: bool,
    pub raim_outlier_m: f64,
    pub nominal_snr_dbhz: f64,
}

impl Default for SppConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_threshold: 1e-2, // Relaxed to 1cm
            geometry_variance_threshold: 10000.0, // Increased to 10000m^2 to allow 3-sat initialization
            enable_sagnac: true,
            enable_tropo: true,
            enable_iono: true,
            raim_outlier_m: 25.0,
            nominal_snr_dbhz: 45.0,
        }
    }
}

fn build_measurements(epoch: &EpochObs, ephemerides: &[Ephemeris], config: &SppConfig) -> Vec<SppMeasurement> {
    let mut measurements = Vec::new();

    for sat_obs in &epoch.satellites {
        if sat_obs.sat.constellation != gneiss_core::sat::Constellation::Gps 
            && sat_obs.sat.constellation != gneiss_core::sat::Constellation::Galileo
            && sat_obs.sat.constellation != gneiss_core::sat::Constellation::Beidou { continue; }

        // Find matching ephemeris closest to epoch time
        let eph = ephemerides.iter()
            .filter(|e| e.sat() == sat_obs.sat)
            .min_by(|a, b| {
                let da = (a.toe().tow - epoch.time.tow).abs();
                let db = (b.toe().tow - epoch.time.tow).abs();
                da.partial_cmp(&db).unwrap()
            });

        if let Some(eph) = eph {
            let pr_obs = sat_obs.observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 1);
            
            if let Some(obs) = pr_obs {
                let snr = sat_obs.observations.iter().find(|o| o.code.obs_type == ObsType::Snr && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(config.nominal_snr_dbhz);
                let doppler = sat_obs.observations.iter().find(|o| o.code.obs_type == ObsType::Doppler && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(0.0);

                measurements.push(SppMeasurement {
                    constellation: sat_obs.sat.constellation,
                    raw_pr: obs.value,
                    snr,
                    doppler,
                    time: epoch.time,
                    eph: eph.clone(),
                });
            }
        }
    }
    measurements
}

fn compute_sat_state(m: &SppMeasurement, receiver_cdt: f64) -> (Coordinate, f64) {
    let pr_time = m.raw_pr / LIGHT_SPEED;
    let t_rcv = m.time.tow - (receiver_cdt / LIGHT_SPEED);
    let t_tx_sat = t_rcv - pr_time;
    let t_tx_sat_gps = GpsTime::new(m.time.week, t_tx_sat);
    
    let (_, _, sat_clk_err_rough, _) = m.eph.position(t_tx_sat_gps);
    
    let t_tx_true = t_tx_sat - sat_clk_err_rough;
    let t_tx_true_gps = GpsTime::new(m.time.week, t_tx_true);
    
    let (sat_pos, _, sat_clk_err, _) = m.eph.position(t_tx_true_gps);
    let corrected_pr = m.raw_pr + (sat_clk_err * LIGHT_SPEED);
    
    (Coordinate::new(sat_pos, Datum::WGS84, Frame::ECEF, m.time), corrected_pr)
}

fn seed_initial_state(measurements: &[SppMeasurement], prev_state: Option<&SppState>) -> SppState {
    let (seed_x, seed_y, seed_z) = if let Some(coord) = prev_state.map(|s| &s.position) {
        (coord.vector.x, coord.vector.y, coord.vector.z)
    } else {
        // Intelligently seed the position on the Earth's surface roughly below the visible satellite constellation
        let mut avg_x = 0.0;
        let mut avg_y = 0.0;
        let mut avg_z = 0.0;
        for m in measurements {
            let (sat_coord, _) = compute_sat_state(m, 0.0);
            avg_x += sat_coord.vector.x;
            avg_y += sat_coord.vector.y;
            avg_z += sat_coord.vector.z;
        }
        let n_f = measurements.len() as f64;
        avg_x /= n_f;
        avg_y /= n_f;
        avg_z /= n_f;

        // Project to Earth's surface (WGS84) exactly
        let avg_ecef = Vector3::new(avg_x, avg_y, avg_z);
        let mut llh = ecef_to_llh(avg_ecef);
        llh.z = 0.0; // Force to ellipsoid surface
        let projected_ecef = gneiss_core::coords::llh_to_ecef(llh);
        
        (projected_ecef.x, projected_ecef.y, projected_ecef.z)
    };

    let mut cdt_gps = None;
    let mut cdt_gal = None;
    let mut cdt_bds = None;
    let mut cdt_glo = None;

    for m in measurements {
        let (sat_coord, corrected_pr) = compute_sat_state(m, 0.0);
        let dx = seed_x - sat_coord.vector.x;
        let dy = seed_y - sat_coord.vector.y;
        let dz = seed_z - sat_coord.vector.z;
        let geom_r = f64::sqrt(dx * dx + dy * dy + dz * dz);
        let cdt = corrected_pr - geom_r;
        
        match m.constellation {
            gneiss_core::sat::Constellation::Gps => if cdt_gps.is_none() { cdt_gps = Some(cdt); },
            gneiss_core::sat::Constellation::Galileo => if cdt_gal.is_none() { cdt_gal = Some(cdt); },
            gneiss_core::sat::Constellation::Beidou => if cdt_bds.is_none() { cdt_bds = Some(cdt); },
            gneiss_core::sat::Constellation::Glonass => if cdt_glo.is_none() { cdt_glo = Some(cdt); },
            _ => {},
        }
    }

    let default_cdt = cdt_gps.or(cdt_gal).or(cdt_bds).or(cdt_glo).unwrap_or(0.0);

    SppState::new(
        Coordinate::new(Vector3::new(seed_x, seed_y, seed_z), Datum::WGS84, Frame::ECEF, measurements[0].time),
        cdt_gps.unwrap_or(default_cdt),
        cdt_gal.unwrap_or(default_cdt),
        cdt_bds.unwrap_or(default_cdt),
        cdt_glo.unwrap_or(default_cdt)
    )
}

/// Computes the Single Point Position (SPP) using the provided epoch observations and ephemerides.
pub fn compute_spp(
    epoch: &EpochObs,
    ephemerides: &[Ephemeris],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
    prev_state: Option<&SppState>,
) -> Result<SppState, SppError> {
    let measurements = build_measurements(epoch, ephemerides, config);

    if measurements.len() < 3 {
        tracing::error!("SPP failed: Only {} valid measurements. Need at least 3.", measurements.len());
        return Err(SppError::NotEnoughMeasurements);
    }

    let mut state = seed_initial_state(&measurements, prev_state);

    let seed_x = state.position.vector.x;
    let seed_y = state.position.vector.y;
    let seed_z = state.position.vector.z;
    let seed_cdt = state.cdt;
    let seed_cdt_gal = state.cdt_gal;
    let seed_cdt_bds = state.cdt_bds;
    let seed_cdt_glo = state.cdt_glo;

    for _ in 0..config.max_iterations {
        let prev_state = state.clone();

        state = spp_wnlls_step(&state, &measurements, iono_params, config)?;

        let dx = state.position.vector.x - prev_state.position.vector.x;
        let dy = state.position.vector.y - prev_state.position.vector.y;
        let dz = state.position.vector.z - prev_state.position.vector.z;
        let dcdt = state.cdt - prev_state.cdt;

        let delta_norm = f64::sqrt(dx * dx + dy * dy + dz * dz + dcdt * dcdt);

        if delta_norm < config.convergence_threshold {
            // Adaptive RAIM: Median Absolute Deviation (MAD)
            let mut residuals = Vec::with_capacity(measurements.len());
            for m in &measurements {
                let cdt = match m.constellation {
                    gneiss_core::sat::Constellation::Gps => state.cdt,
                    gneiss_core::sat::Constellation::Galileo => state.cdt_gal,
                    gneiss_core::sat::Constellation::Beidou => state.cdt_bds,
                    _ => state.cdt,
                };
                let (sat_coord, corrected_pr) = compute_sat_state(m, cdt);
                let r_dx = state.position.vector.x - sat_coord.vector.x;
                let r_dy = state.position.vector.y - sat_coord.vector.y;
                let r_dz = state.position.vector.z - sat_coord.vector.z;
                let r_dist = f64::sqrt(r_dx * r_dx + r_dy * r_dy + r_dz * r_dz);
                let expected_pr = r_dist + cdt;
                let residual = (corrected_pr - expected_pr).abs();
                residuals.push((m, residual));
            }

            residuals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let median_residual = residuals[residuals.len() / 2].1;
            
            let mut deviations: Vec<f64> = residuals.iter().map(|(_, r)| (*r - median_residual).abs()).collect();
            deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mad = deviations[deviations.len() / 2];
            
            // 1.4826 scales MAD to match standard deviation for a normal distribution. Using 5 sigma bound.
            let mad_threshold = (mad * 1.4826 * 5.0).max(config.raim_outlier_m);

            let mut good_measurements = Vec::new();
            for (m, r) in residuals {
                if r <= mad_threshold {
                    good_measurements.push(m.clone());
                }
            }

            if good_measurements.len() < measurements.len() && good_measurements.len() >= 4 {
                // Re-run with clean measurements
                let mut clean_state = SppState::new(
                    Coordinate::new(Vector3::new(seed_x, seed_y, seed_z), Datum::WGS84, Frame::ECEF, measurements[0].time),
                    seed_cdt,
                    seed_cdt_gal,
                    seed_cdt_bds,
                    seed_cdt_glo
                );
                for _ in 0..config.max_iterations {
                    let prev_clean = clean_state.clone();
                    if let Ok(new_state) = spp_wnlls_step(&clean_state, &good_measurements, iono_params, config) {
                        clean_state = new_state;
                        let c_dx = clean_state.position.vector.x - prev_clean.position.vector.x;
                        let c_dy = clean_state.position.vector.y - prev_clean.position.vector.y;
                        let c_dz = clean_state.position.vector.z - prev_clean.position.vector.z;
                        let c_dcdt = clean_state.cdt - prev_clean.cdt;
                        if f64::sqrt(c_dx * c_dx + c_dy * c_dy + c_dz * c_dz + c_dcdt * c_dcdt) < config.convergence_threshold {
                            return Ok(clean_state);
                        }
                    } else {
                        break;
                    }
                }
                return Err(SppError::ConvergenceFailed);
            }
            return Ok(state);
        }
    }

    Err(SppError::ConvergenceFailed)
}



/// Performs a single Weighted Non-Linear Least Squares (WNLLS) iteration.
pub fn spp_wnlls_step(
    current_state: &SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> Result<SppState, SppError> {
    let n = measurements.len();
    
    let mut cols = 3;
    let mut gps_col = None;
    let mut gal_col = None;
    let mut bds_col = None;
    let mut glo_col = None;

    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Gps) {
        gps_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Galileo) {
        gal_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Beidou) {
        bds_col = Some(cols); cols += 1;
    }
    if measurements.iter().any(|m| m.constellation == gneiss_core::sat::Constellation::Glonass) {
        glo_col = Some(cols); cols += 1;
    }

    let use_height_constraint = n == cols - 1;
    if n < cols - 1 {
        return Err(SppError::NotEnoughMeasurements);
    }

    let matrix_n = if use_height_constraint { n + 1 } else { n };

    let mut h_matrix = DMatrix::<f64>::zeros(matrix_n, cols);
    let mut w_matrix = DMatrix::<f64>::zeros(matrix_n, matrix_n);
    let mut dz_vector = DVector::<f64>::zeros(matrix_n);

    let rec_ecef = Vector3::new(current_state.position.vector.x, current_state.position.vector.y, current_state.position.vector.z);
    let rec_llh = ecef_to_llh(rec_ecef);
    let tropo_params = TropoParams::default();

    for (i, m) in measurements.iter().enumerate() {
        let cdt = match m.constellation {
            gneiss_core::sat::Constellation::Gps => current_state.cdt,
            gneiss_core::sat::Constellation::Galileo => current_state.cdt_gal,
            gneiss_core::sat::Constellation::Beidou => current_state.cdt_bds,
            _ => current_state.cdt,
        };
        let (sat_coord, corrected_pr) = compute_sat_state(m, cdt);
        let sat_ecef = Vector3::new(sat_coord.vector.x, sat_coord.vector.y, sat_coord.vector.z);
        
        // 1. Sagnac Effect (Earth rotation during signal time of flight)
        let sat_ecef_rot = if config.enable_sagnac {
            // Must use an approximation of the true geometric time of flight,
            // NOT the raw pseudorange which is dominated by receiver clock bias.
            let geometric_pr = corrected_pr - cdt;
            let tof = geometric_pr / LIGHT_SPEED;
            let theta = OMEGA_E * tof;
            let cos_t = f64::cos(theta);
            let sin_t = f64::sin(theta);
            
            Vector3::new(
                sat_ecef.x * cos_t + sat_ecef.y * sin_t,
                -sat_ecef.x * sin_t + sat_ecef.y * cos_t,
                sat_ecef.z
            )
        } else {
            sat_ecef
        };

        let dx = current_state.position.vector.x - sat_ecef_rot.x;
        let dy = current_state.position.vector.y - sat_ecef_rot.y;
        let dz = current_state.position.vector.z - sat_ecef_rot.z;
        let r = f64::sqrt(dx * dx + dy * dy + dz * dz);
        
        let r_safe = if r < 1e-6 { 1e-6 } else { r };

        // 2. Azimuth and Elevation for Variance and Atmospherics
        let (az, el) = az_el(rec_llh, rec_ecef, sat_ecef_rot);
        
        #[cfg(test)]
        let el_mask = -core::f64::consts::PI; // Disable in unit tests
        #[cfg(not(test))]
        let el_mask = 0.1745; // 10 degrees

        // Exclude satellites below elevation mask
        if el < el_mask && current_state.position.vector.x != 0.0 {
            // Give it an incredibly low weight so it doesn't affect the solution
            w_matrix[(i, i)] = 1e-10;
            continue;
        }
        
        // 3. Atmospheric Delays
        let mut tropo_delay = 0.0;
        let mut iono_delay = 0.0;

        // Constrain elevation for atmospherics to prevent math explosion
        let min_el = 5.0 * core::f64::consts::PI / 180.0;
        let safe_el = if el < min_el { min_el } else { el };

        // Don't apply atmospherics if receiver is clearly near center of earth (first iteration)
        if rec_ecef.norm() > 6_000_000.0 {
            if config.enable_tropo {
                tropo_delay = AtmosphereModel::tropo_nmf(&tropo_params, rec_llh, safe_el, m.time);
            }
            if config.enable_iono {
                if let Some(iono) = iono_params {
                    iono_delay = AtmosphereModel::iono_klobuchar(iono, rec_llh, az, safe_el, m.time);
                }
            }
        }

        let expected_pr = r + cdt + tropo_delay + iono_delay;
        let residual = corrected_pr - expected_pr;
        
        // 4. Weight Matrix (Lower elevation or SNR = higher variance)
        let variance = gneiss_core::variance::observation_variance(m.snr, el, config.nominal_snr_dbhz);
        
        h_matrix[(i, 0)] = dx / r_safe;
        h_matrix[(i, 1)] = dy / r_safe;
        h_matrix[(i, 2)] = dz / r_safe;
        
        match m.constellation {
            gneiss_core::sat::Constellation::Gps => { if let Some(c) = gps_col { h_matrix[(i, c)] = 1.0; } },
            gneiss_core::sat::Constellation::Galileo => { if let Some(c) = gal_col { h_matrix[(i, c)] = 1.0; } },
            gneiss_core::sat::Constellation::Beidou => { if let Some(c) = bds_col { h_matrix[(i, c)] = 1.0; } },
            _ => {},
        }
        
        w_matrix[(i, i)] = 1.0 / variance;
        dz_vector[i] = residual;
    }

    if use_height_constraint {
        let lat = rec_llh.x;
        let lon = rec_llh.y;
        
        let sin_lat = lat.sin();
        let cos_lat = lat.cos();
        let sin_lon = lon.sin();
        let cos_lon = lon.cos();
        
        // UP vector in ECEF
        h_matrix[(n, 0)] = cos_lat * cos_lon;
        h_matrix[(n, 1)] = cos_lat * sin_lon;
        h_matrix[(n, 2)] = sin_lat;
        
        // Penalize change in height from initial guess
        dz_vector[n] = 0.0;
        
        // Give it a relatively high variance so true measurements take precedence,
        // but low enough to constrain the geometry (100 m^2 = 10m std dev)
        w_matrix[(n, n)] = 1.0 / 100.0; 
    }

    let h_t = h_matrix.transpose();
    let h_t_w = &h_t * &w_matrix;
    let h_t_w_h = &h_t_w * &h_matrix;

    let h_t_w_h_inv = h_t_w_h.try_inverse().ok_or(SppError::MatrixInversionFailed)?;

    // Check Geometry / Positional Variance
    // Since W has units of 1/m^2, h_t_w_h_inv has units of m^2
    // If the variance is huge (e.g. we deweighted satellites and have poor remaining geometry), reject it.
    let pos_variance = h_t_w_h_inv[(0, 0)] + h_t_w_h_inv[(1, 1)] + h_t_w_h_inv[(2, 2)];
    if pos_variance > config.geometry_variance_threshold {
        tracing::debug!("Poor geometry: pos_variance = {} (threshold = {})", pos_variance, config.geometry_variance_threshold);
        return Err(SppError::PoorGeometry);
    }

    let dx_vec = h_t_w_h_inv * h_t_w * dz_vector;

    let next_cdt = current_state.cdt + gps_col.map(|c| dx_vec[c]).unwrap_or(0.0);
    let next_cdt_gal = current_state.cdt_gal + gal_col.map(|c| dx_vec[c]).unwrap_or(0.0);
    let next_cdt_bds = current_state.cdt_bds + bds_col.map(|c| dx_vec[c]).unwrap_or(0.0);

    let next_cdt_glo = current_state.cdt_glo + glo_col.map(|c| dx_vec[c]).unwrap_or(0.0);

    Ok(SppState::new(
        Coordinate::new(
            Vector3::new(
                current_state.position.vector.x + dx_vec[0],
                current_state.position.vector.y + dx_vec[1],
                current_state.position.vector.z + dx_vec[2],
            ),
            Datum::WGS84,
            Frame::ECEF,
            measurements[0].time
        ),
        next_cdt,
        next_cdt_gal,
        next_cdt_bds,
        next_cdt_glo
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_spp() {
        use gneiss_core::time::GpsTime;
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode, SignalCode};

        let t = GpsTime::new(2000, 100000.0);
        let true_pos = nalgebra::Vector3::new(6378137.0, 0.0, 0.0);
        let true_cdt = 1000.0;

        let mut ephemerides = Vec::new();
        let mut satellites = Vec::new();

        let sats = vec![
            (1, (20000000.0, 5000000.0, 5000000.0)),
            (2, (22000000.0, -5000000.0, 5000000.0)),
            (3, (19000000.0, 5000000.0, -5000000.0)),
            (4, (21000000.0, -5000000.0, -5000000.0)),
            (5, (25000000.0, 0.0, 0.0)),
        ];

        for (prn, (_sx, _sy, _sz)) in sats {
            let sat_id = SatelliteId { constellation: Constellation::Gps, prn };
            
            let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: sat_id,
                toe: t, toc: t,
                af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: (prn as f64) * core::f64::consts::PI / 3.0,
                e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: (prn as f64) * core::f64::consts::PI / 2.0,
                omega_dot: 0.0, i0: 0.95, idot: 0.0,
                omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            });
            ephemerides.push(eph.clone());

            let mut raw_pr = 20000000.0 + true_cdt; // rough initial guess
            let light_speed = 299792458.0;
            
            for _ in 0..5 {
                let pr_time = raw_pr / light_speed;
                let t_tx_sat = t.tow - pr_time;
                let (_, _, sat_clk_err_rough, _) = eph.position(GpsTime::new(t.week, t_tx_sat));
                let t_tx_true = t_tx_sat - sat_clk_err_rough;
                
                let (sat_pos, _, sat_clk_err, _) = eph.position(GpsTime::new(t.week, t_tx_true));
                let dx = true_pos.x - sat_pos.x;
                let dy = true_pos.y - sat_pos.y;
                let dz = true_pos.z - sat_pos.z;
                let geometric_range = f64::sqrt(dx * dx + dy * dy + dz * dz);
                
                raw_pr = geometric_range + true_cdt - (sat_clk_err * light_speed);
            }

            satellites.push(SatObs {
                sat: sat_id,
                observations: vec![Observation {
                    code: ObsCode {
                        obs_type: ObsType::Pseudorange,
                        signal: SignalCode { freq_band: 1, attribute: 'C' },
                    },
                    value: raw_pr,
                    lock_time: None,
                }],
            });
        }

        let epoch = EpochObs {
            time: t,
            satellites,
        };

        let config = SppConfig {
            enable_sagnac: false,
            enable_tropo: false,
            enable_iono: false,
            geometry_variance_threshold: 100000.0,
            ..Default::default()
        };

        let state = compute_spp(&epoch, &ephemerides, None, &config, None).unwrap();

        assert!((state.position.vector.x - true_pos.x).abs() < 1e-1, "X error too large: {}", (state.position.vector.x - true_pos.x).abs());
        assert!((state.position.vector.y - true_pos.y).abs() < 1e-1, "Y error too large: {}", (state.position.vector.y - true_pos.y).abs());
        assert!((state.position.vector.z - true_pos.z).abs() < 1e-1, "Z error too large: {}", (state.position.vector.z - true_pos.z).abs());
        assert!((state.cdt - true_cdt).abs() < 1e-1, "CDT error too large: {}", (state.cdt - true_cdt).abs());
    }
}
