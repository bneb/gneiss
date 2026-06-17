use nalgebra::{DMatrix, DVector, Vector3};
use gneiss_core::obs::EpochObs;
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

    pub fn get_cdt(&self, constellation: gneiss_core::sat::Constellation) -> f64 {
        match constellation {
            gneiss_core::sat::Constellation::Galileo => self.cdt_gal,
            gneiss_core::sat::Constellation::Beidou => self.cdt_bds,
            gneiss_core::sat::Constellation::Glonass => self.cdt_glo,
            _ => self.cdt,
        }
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
    pub is_iono_free: bool,
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

const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
const OMEGA_E: f64 = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S; // WGS 84 value of earth's rotation rate
const MIN_WEIGHT: f64 = 1e-10;
const HEIGHT_CONSTRAINT_VAR: f64 = 100.0;
const MIN_EARTH_RADIUS_M: f64 = 6_000_000.0;
const MIN_ATMOSPHERE_ELEVATION_RAD: f64 = 5.0 * core::f64::consts::PI / 180.0;

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
    pub snr_a: f64,
    pub snr_b: f64,
    pub min_measurements_init: usize,
    pub raim_mad_multiplier: f64,
    pub elevation_mask_rad: f64,
}

impl Default for SppConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_threshold: 1e-4,
            geometry_variance_threshold: 10000.0,
            enable_sagnac: true,
            enable_tropo: true,
            enable_iono: true,
            raim_outlier_m: 25.0,
            snr_a: 1.0,
            snr_b: 150.0,
            min_measurements_init: 3,
            raim_mad_multiplier: 7.413, // 1.4826 * 5.0 sigma
            elevation_mask_rad: 0.261799, // 15 degrees
        }
    }
}

fn build_single_measurement(sat_obs: &gneiss_core::obs::SatObs, ephemerides: &[Ephemeris], epoch_time: GpsTime) -> Option<SppMeasurement> {
    let eph = ephemerides.iter()
        .filter(|e| e.sat() == sat_obs.sat)
        .min_by(|a, b| (a.toe().tow - epoch_time.tow).abs().partial_cmp(&(b.toe().tow - epoch_time.tow).abs()).unwrap())?;

    let (f1, mut f2) = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num());
    if f2 == 0.0 { f2 = f1; }

    let p1_opt = match sat_obs.sat.constellation {
        gneiss_core::sat::Constellation::Beidou => sat_obs.get_observable(2),
        _ => sat_obs.get_observable(1),
    };
    let p2_opt = match sat_obs.sat.constellation {
        gneiss_core::sat::Constellation::Galileo => sat_obs.get_observable(7).or(sat_obs.get_observable(5)),
        gneiss_core::sat::Constellation::Beidou => sat_obs.get_observable(7).or(sat_obs.get_observable(6)),
        _ => sat_obs.get_observable(2),
    };

    let (raw_pr, is_iono_free) = if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
        let f1_sq = f1 * f1; let f2_sq = f2 * f2;
        ((f1_sq * p1 - f2_sq * p2) / (f1_sq - f2_sq), true)
    } else {
        (p1_opt?, false)
    };

    Some(SppMeasurement {
        constellation: sat_obs.sat.constellation,
        raw_pr,
        snr: sat_obs.get_snr(1).unwrap_or(45) as f64,
        doppler: sat_obs.get_doppler(1).unwrap_or(0.0),
        time: epoch_time,
        eph: eph.clone(),
        is_iono_free,
    })
}

pub fn build_measurements(epoch: &EpochObs, ephemerides: &[Ephemeris], _config: &SppConfig) -> Vec<SppMeasurement> {
    epoch.satellites.iter().filter_map(|s| build_single_measurement(s, ephemerides, epoch.time)).collect()
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

fn compute_seed_position(measurements: &[SppMeasurement]) -> (f64, f64, f64) {
    let mut avg = Vector3::zeros();
    for m in measurements {
        let (sat_coord, _) = compute_sat_state(m, 0.0);
        avg += sat_coord.vector;
    }
    avg /= measurements.len() as f64;
    let mut llh = ecef_to_llh(avg);
    llh.z = 0.0;
    let proj = gneiss_core::coords::llh_to_ecef(llh);
    (proj.x, proj.y, proj.z)
}

fn compute_seed_clocks(measurements: &[SppMeasurement], sx: f64, sy: f64, sz: f64) -> (f64, f64, f64, f64) {
    let mut cdt_gps = None; let mut cdt_gal = None; let mut cdt_bds = None; let mut cdt_glo = None;
    for m in measurements {
        let (sat_coord, corrected_pr) = compute_sat_state(m, 0.0);
        let dx = sx - sat_coord.vector.x; let dy = sy - sat_coord.vector.y; let dz = sz - sat_coord.vector.z;
        let cdt = corrected_pr - f64::sqrt(dx * dx + dy * dy + dz * dz);
        match m.constellation {
            gneiss_core::sat::Constellation::Gps => if cdt_gps.is_none() { cdt_gps = Some(cdt); },
            gneiss_core::sat::Constellation::Galileo => if cdt_gal.is_none() { cdt_gal = Some(cdt); },
            gneiss_core::sat::Constellation::Beidou => if cdt_bds.is_none() { cdt_bds = Some(cdt); },
            gneiss_core::sat::Constellation::Glonass => if cdt_glo.is_none() { cdt_glo = Some(cdt); },
            _ => {},
        }
        tracing::debug!("SPP seed: SAT={}, raw_pr={:.3}, cdt={:.3}", m.eph.sat(), m.raw_pr, cdt);
    }
    let default_cdt = cdt_gps.or(cdt_gal).or(cdt_bds).or(cdt_glo).unwrap_or(0.0);
    (cdt_gps.unwrap_or(default_cdt), cdt_gal.unwrap_or(default_cdt), cdt_bds.unwrap_or(default_cdt), cdt_glo.unwrap_or(default_cdt))
}

fn seed_initial_state(measurements: &[SppMeasurement], prev_state: Option<&SppState>) -> SppState {
    let (seed_x, seed_y, seed_z) = if let Some(coord) = prev_state.map(|s| &s.position) {
        (coord.vector.x, coord.vector.y, coord.vector.z)
    } else {
        compute_seed_position(measurements)
    };
    let (cdt, cdt_gal, cdt_bds, cdt_glo) = compute_seed_clocks(measurements, seed_x, seed_y, seed_z);
    SppState::new(
        Coordinate::new(Vector3::new(seed_x, seed_y, seed_z), Datum::WGS84, Frame::ECEF, measurements[0].time),
        cdt, cdt_gal, cdt_bds, cdt_glo
    )
}

pub fn compute_spp(
    epoch: &EpochObs, ephemerides: &[Ephemeris], iono_params: Option<&KlobucharParams>, config: &SppConfig, prev_state: Option<&SppState>,
) -> Result<SppState, SppError> {
    let measurements = build_measurements(epoch, ephemerides, config);
    if measurements.len() < config.min_measurements_init {
        tracing::error!("SPP failed: Only {} valid measurements. Need at least {}.", measurements.len(), config.min_measurements_init);
        return Err(SppError::NotEnoughMeasurements);
    }
    let seed_state = seed_initial_state(&measurements, prev_state);
    let state = solve_spp_iteratively(seed_state.clone(), &measurements, iono_params, config)?;
    apply_raim(state, seed_state, &measurements, iono_params, config)
}

fn solve_spp_iteratively(
    mut state: SppState, measurements: &[SppMeasurement], iono_params: Option<&KlobucharParams>, config: &SppConfig,
) -> Result<SppState, SppError> {
    for _ in 0..config.max_iterations {
        let prev_state = state.clone();
        state = spp_wnlls_step(&state, measurements, iono_params, config)?;
        let dx = state.position.vector.x - prev_state.position.vector.x;
        let dy = state.position.vector.y - prev_state.position.vector.y;
        let dz = state.position.vector.z - prev_state.position.vector.z;
        let dcdt = state.cdt - prev_state.cdt;
        if f64::sqrt(dx * dx + dy * dy + dz * dz + dcdt * dcdt) < config.convergence_threshold {
            return Ok(state);
        }
    }
    Err(SppError::ConvergenceFailed)
}

fn apply_raim(
    state: SppState, seed_state: SppState, measurements: &[SppMeasurement], iono_params: Option<&KlobucharParams>, config: &SppConfig,
) -> Result<SppState, SppError> {
    let good_measurements = filter_raim_outliers(&state, measurements, iono_params, config);
    if good_measurements.len() < measurements.len() && good_measurements.len() >= 4 {
        return solve_spp_iteratively(seed_state, &good_measurements, iono_params, config);
    }
    Ok(state)
}

fn filter_raim_outliers(state: &SppState, measurements: &[SppMeasurement], iono_params: Option<&KlobucharParams>, config: &SppConfig) -> Vec<SppMeasurement> {
    let rec_ecef = state.position.vector;
    let rec_llh = ecef_to_llh(rec_ecef);

    let mut residuals: Vec<(&SppMeasurement, f64)> = measurements.iter().map(|m| {
        let (_, _, _, _, residual, _) = compute_measurement_residuals(state, m, rec_ecef, rec_llh, iono_params, config);
        (m, residual.abs())
    }).collect();

    residuals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let median_residual = residuals[residuals.len() / 2].1;
    let mut deviations: Vec<f64> = residuals.iter().map(|(_, r)| (*r - median_residual).abs()).collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad_threshold = (deviations[deviations.len() / 2] * config.raim_mad_multiplier).max(config.raim_outlier_m);

    residuals.into_iter().filter(|(_, r)| *r <= mad_threshold).map(|(m, _)| m.clone()).collect()
}



/// Performs a single Weighted Non-Linear Least Squares (WNLLS) iteration.
pub fn spp_wnlls_step(
    current_state: &SppState, measurements: &[SppMeasurement], iono_params: Option<&KlobucharParams>, config: &SppConfig,
) -> Result<SppState, SppError> {
    let (mut h_matrix, mut w_matrix, mut dz_vector, cols, clocks) = build_design_matrix(current_state, measurements, iono_params, config)?;
    let n = measurements.len();
    if n == cols - 1 {
        let rec_llh = ecef_to_llh(Vector3::new(current_state.position.vector.x, current_state.position.vector.y, current_state.position.vector.z));
        apply_height_constraint(rec_llh, n, &mut h_matrix, &mut w_matrix, &mut dz_vector);
    }

    let h_t = h_matrix.transpose();
    let h_t_w = &h_t * &w_matrix;
    let h_t_w_h_inv_opt = (&h_t_w * &h_matrix).try_inverse();
    let h_t_w_h_inv = match h_t_w_h_inv_opt {
        Some(inv) => inv,
        None => {
            tracing::debug!("SPP MatrixInversionFailed. Measurements: {}", measurements.len());
            return Err(SppError::MatrixInversionFailed);
        }
    };

    let trace = h_t_w_h_inv[(0, 0)] + h_t_w_h_inv[(1, 1)] + h_t_w_h_inv[(2, 2)];
    if trace > config.geometry_variance_threshold { 
        let valid_sats = w_matrix.diagonal().iter().filter(|&w| *w > MIN_WEIGHT).count();
        tracing::debug!("SPP PoorGeometry. Trace: {:.3}, threshold: {:.3}, valid_sats: {}, total_sats: {}", trace, config.geometry_variance_threshold, valid_sats, measurements.len());
        return Err(SppError::PoorGeometry); 
    }

    let dx_vec = h_t_w_h_inv * h_t_w * dz_vector;
    tracing::debug!("SPP iter: dx_vec=[{:.2}, {:.2}, {:.2}, {:.2}], trace={:.2}", dx_vec[0], dx_vec[1], dx_vec[2], dx_vec[3], trace);

    
    Ok(SppState::new(
        Coordinate::new(
            Vector3::new(current_state.position.vector.x + dx_vec[0], current_state.position.vector.y + dx_vec[1], current_state.position.vector.z + dx_vec[2]), 
            Datum::WGS84, Frame::ECEF, measurements[0].time
        ),
        current_state.cdt + clocks.0.map(|c| dx_vec[c]).unwrap_or(0.0),
        current_state.cdt_gal + clocks.1.map(|c| dx_vec[c]).unwrap_or(0.0),
        current_state.cdt_bds + clocks.2.map(|c| dx_vec[c]).unwrap_or(0.0),
        current_state.cdt_glo + clocks.3.map(|c| dx_vec[c]).unwrap_or(0.0)
    ))
}

struct ClockCols(Option<usize>, Option<usize>, Option<usize>, Option<usize>);

fn find_clock_cols(measurements: &[SppMeasurement]) -> (usize, ClockCols) {
    let mut cols = 3;
    let has = |c| measurements.iter().any(|m| m.constellation == c);
    let gps_col = if has(gneiss_core::sat::Constellation::Gps) || has(gneiss_core::sat::Constellation::Qzss) { cols += 1; Some(cols - 1) } else { None };
    let gal_col = if has(gneiss_core::sat::Constellation::Galileo) { cols += 1; Some(cols - 1) } else { None };
    let bds_col = if has(gneiss_core::sat::Constellation::Beidou) { cols += 1; Some(cols - 1) } else { None };
    let glo_col = if has(gneiss_core::sat::Constellation::Glonass) { cols += 1; Some(cols - 1) } else { None };
    (cols, ClockCols(gps_col, gal_col, bds_col, glo_col))
}

fn build_design_matrix(
    state: &SppState, measurements: &[SppMeasurement], iono_params: Option<&KlobucharParams>, config: &SppConfig,
) -> Result<(DMatrix<f64>, DMatrix<f64>, DVector<f64>, usize, ClockCols), SppError> {
    let (cols, clocks) = find_clock_cols(measurements);
    let n = measurements.len();
    if n < cols - 1 { return Err(SppError::NotEnoughMeasurements); }
    let matrix_n = if n == cols - 1 { n + 1 } else { n };

    let mut h_matrix = DMatrix::<f64>::zeros(matrix_n, cols);
    let mut w_matrix = DMatrix::<f64>::zeros(matrix_n, matrix_n);
    let mut dz_vector = DVector::<f64>::zeros(matrix_n);
    let rec_ecef = state.position.vector;
    let rec_llh = ecef_to_llh(rec_ecef);

    for (i, m) in measurements.iter().enumerate() {
        let (dx, dy, dz, r, residual, el) = compute_measurement_residuals(state, m, rec_ecef, rec_llh, iono_params, config);
        if el < config.elevation_mask_rad && rec_ecef.x != 0.0 { w_matrix[(i, i)] = MIN_WEIGHT; continue; }

        h_matrix[(i, 0)] = dx / r; h_matrix[(i, 1)] = dy / r; h_matrix[(i, 2)] = dz / r;
        if let Some(c) = match m.constellation {
            gneiss_core::sat::Constellation::Gps | gneiss_core::sat::Constellation::Qzss => clocks.0,
            gneiss_core::sat::Constellation::Galileo => clocks.1,
            gneiss_core::sat::Constellation::Beidou => clocks.2,
            gneiss_core::sat::Constellation::Glonass => clocks.3,
            _ => None,
        } { h_matrix[(i, c)] = 1.0; }
        
        w_matrix[(i, i)] = 1.0 / gneiss_core::variance::observation_variance(m.snr, el, config.snr_a, config.snr_b);
        dz_vector[i] = residual;
    }
    
    Ok((h_matrix, w_matrix, dz_vector, cols, clocks))
}

fn compute_measurement_residuals(
    current_state: &SppState, m: &SppMeasurement, rec_ecef: Vector3<f64>, rec_llh: Vector3<f64>,
    iono_params: Option<&KlobucharParams>, config: &SppConfig,
) -> (f64, f64, f64, f64, f64, f64) {
    let cdt = current_state.get_cdt(m.constellation);
    let (sat_coord, corrected_pr) = compute_sat_state(m, cdt);
    
    let sat_ecef = if config.enable_sagnac {
        compute_sagnac_correction(sat_coord.vector, corrected_pr - cdt)
    } else { sat_coord.vector };

    let dx = rec_ecef.x - sat_ecef.x;
    let dy = rec_ecef.y - sat_ecef.y;
    let dz = rec_ecef.z - sat_ecef.z;
    let r = f64::sqrt(dx * dx + dy * dy + dz * dz).max(1e-6);

    let (az, el) = az_el(rec_llh, rec_ecef, sat_ecef);
    let (tropo, mut iono) = compute_atmospheric_delays(rec_ecef, rec_llh, az, el, m.time, iono_params, config);
    if m.is_iono_free { iono = 0.0; }

    (dx, dy, dz, r, corrected_pr - (r + cdt + tropo + iono), el)
}

fn compute_sagnac_correction(sat_ecef: Vector3<f64>, geometric_pr: f64) -> Vector3<f64> {
    let tof = geometric_pr / LIGHT_SPEED;
    let theta = OMEGA_E * tof;
    let cos_t = f64::cos(theta);
    let sin_t = f64::sin(theta);
    Vector3::new(sat_ecef.x * cos_t + sat_ecef.y * sin_t, -sat_ecef.x * sin_t + sat_ecef.y * cos_t, sat_ecef.z)
}

fn compute_atmospheric_delays(
    rec_ecef: Vector3<f64>, rec_llh: Vector3<f64>, az: f64, el: f64, time: gneiss_core::time::GpsTime,
    iono_params: Option<&KlobucharParams>, config: &SppConfig
) -> (f64, f64) {
    if rec_ecef.norm() <= MIN_EARTH_RADIUS_M { return (0.0, 0.0); }
    let safe_el = el.max(MIN_ATMOSPHERE_ELEVATION_RAD);
    
    let tropo = if config.enable_tropo {
        AtmosphereModel::tropo_nmf(&TropoParams::default(), rec_llh, safe_el, time)
    } else { 0.0 };
    
    let iono = if config.enable_iono {
        if let Some(iono) = iono_params {
            AtmosphereModel::iono_klobuchar(iono, rec_llh, az, safe_el, time)
        } else { 0.0 }
    } else { 0.0 };
    
    (tropo, iono)
}
fn apply_height_constraint(
    rec_llh: Vector3<f64>,
    row_idx: usize,
    h_matrix: &mut DMatrix<f64>,
    w_matrix: &mut DMatrix<f64>,
    dz_vector: &mut DVector<f64>
) {
    let lat = rec_llh.x;
    let lon = rec_llh.y;
    
    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_lon = lon.sin();
    let cos_lon = lon.cos();
    
    // UP vector in ECEF
    h_matrix[(row_idx, 0)] = cos_lat * cos_lon;
    h_matrix[(row_idx, 1)] = cos_lat * sin_lon;
    h_matrix[(row_idx, 2)] = sin_lat;
    
    // Penalize change in height from initial guess
    dz_vector[row_idx] = 0.0;
    
    // Give it a relatively high variance so true measurements take precedence,
    // but low enough to constrain the geometry (100 m^2 = 10m std dev)
    w_matrix[(row_idx, row_idx)] = 1.0 / HEIGHT_CONSTRAINT_VAR; 
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::obs::ObsType;

    #[test]
    fn test_compute_spp() {
        use gneiss_core::time::GpsTime;
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode, SignalCode};

        let t = GpsTime::new(2000, 100000.0);
        let true_pos = nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
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
            let light_speed = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
            
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
                    lli: None,
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
            elevation_mask_rad: -core::f64::consts::PI, // Bypass elevation mask for tests
            ..Default::default()
        };

        let state = compute_spp(&epoch, &ephemerides, None, &config, None).unwrap();

        assert!((state.position.vector.x - true_pos.x).abs() < 1e-1, "X error too large: {}", (state.position.vector.x - true_pos.x).abs());
        assert!((state.position.vector.y - true_pos.y).abs() < 1e-1, "Y error too large: {}", (state.position.vector.y - true_pos.y).abs());
        assert!((state.position.vector.z - true_pos.z).abs() < 1e-1, "Z error too large: {}", (state.position.vector.z - true_pos.z).abs());
        assert!((state.cdt - true_cdt).abs() < 1e-1, "CDT error too large: {}", (state.cdt - true_cdt).abs());
    }
    #[test]
    fn test_spp_not_enough_measurements() {
        use gneiss_core::time::GpsTime;
        use gneiss_core::obs::EpochObs;

        let t = GpsTime::new(2000, 100000.0);
        let epoch = EpochObs {
            time: t,
            satellites: vec![],
        };
        let config = SppConfig {
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };
        let res = compute_spp(&epoch, &[], None, &config, None);
        assert_eq!(res.unwrap_err(), SppError::NotEnoughMeasurements);
    }

    #[test]
    fn test_spp_wnlls_step_not_enough_measurements() {
        use gneiss_core::time::GpsTime;
        use gneiss_core::sat::{Constellation, SatelliteId};

        let t = GpsTime::new(2000, 100000.0);
        let m1 = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            }),
            is_iono_free: false,
        };

        let config = SppConfig {
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };
        let state = SppState::new(Coordinate::new(nalgebra::Vector3::zeros(), Datum::WGS84, Frame::ECEF, t), 0.0, 0.0, 0.0, 0.0);
        
        let res = spp_wnlls_step(&state, &[m1], None, &config);
        assert_eq!(res.unwrap_err(), SppError::NotEnoughMeasurements);
    }

    #[test]
    fn test_build_measurements() {
        use gneiss_core::time::GpsTime;
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode, SignalCode, ObsType};
        
        let t = GpsTime::new(2000, 100000.0);
        let t_eph1 = GpsTime::new(2000, 100010.0); // da = 10
        let t_eph2 = GpsTime::new(2000, 99980.0); // db = 20. eph1 is closer!
        
        let sat_gps = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_glo = SatelliteId { constellation: Constellation::Glonass, prn: 1 };
        let sat_gal = SatelliteId { constellation: Constellation::Galileo, prn: 2 };
        let sat_bds = SatelliteId { constellation: Constellation::Beidou, prn: 3 };

        let mut epoch = EpochObs { time: t, satellites: vec![] };
        
        // GLONASS omitted for this test

        // Add GPS with no pseudorange (should be skipped by pr check)
        epoch.satellites.push(SatObs {
            sat: sat_gps,
            observations: vec![Observation {
                code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                value: 100000.0, lock_time: None, lli: None,
            }],
        });
        
        // Add Galileo with pseudorange, snr, doppler, and correct freq band
        epoch.satellites.push(SatObs {
            sat: sat_gal,
            observations: vec![
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                    value: 21000000.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Snr, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                    value: 42.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Doppler, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                    value: 1500.0, lock_time: None, lli: None,
                },
                // Wrong freq band pseudorange just to test
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 2, attribute: 'C' } },
                    value: 22000000.0, lock_time: None, lli: None,
                },
            ],
        });
        
        // Beidou omitted for this test

        let eph_gps = gneiss_core::ephemeris::GpsEphemeris {
            sat: sat_gps, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        };
        
        let eph_gal_1 = gneiss_core::ephemeris::GalileoEphemeris {
            sat: sat_gal, toe: t_eph1, toc: t_eph1, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0, omega: 0.0, bgd_e1_e5a: 0.0, iod_nav: 1,
        };
        
        let eph_gal_2 = gneiss_core::ephemeris::GalileoEphemeris {
            sat: sat_gal, toe: t_eph2, toc: t_eph2, af0: 1.0, af1: 0.0, af2: 0.0, // closer in time!
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0, omega: 0.0, bgd_e1_e5a: 0.0, iod_nav: 2,
        };
        
        let eph_bds = gneiss_core::ephemeris::BeidouEphemeris {
            sat: sat_bds, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0, omega: 0.0, tgd1: 0.0, aode: 1, aodc: 1,
        };
        
        let eph_glo = gneiss_core::ephemeris::GlonassEphemeris {
            sat: sat_glo, toe: t, freq_num: 1, tau_n: 0.0, gamma_n: 0.0, delta_tau_n: 0.0,
            x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, vz: 0.0, ax: 0.0, ay: 0.0, az: 0.0,
        };

        let ephemerides = vec![
            Ephemeris::Gps(eph_gps),
            Ephemeris::Galileo(eph_gal_1),
            Ephemeris::Galileo(eph_gal_2.clone()),
            Ephemeris::Beidou(eph_bds),
            Ephemeris::Glonass(eph_glo),
        ];
        
        let config = SppConfig::default();
        let measurements = build_measurements(&epoch, &ephemerides, &config);
        
        assert_eq!(measurements.len(), 1, "Should have exactly 1 measurement (Galileo)");
        
        // Galileo check
        assert_eq!(measurements[0].snr, 42.0);
        assert_eq!(measurements[0].doppler, 1500.0);
        match &measurements[0].eph {
            Ephemeris::Galileo(g) => assert_eq!(g.iod_nav, 1), // eph1 is now closer
            _ => panic!("Wrong ephemeris type"),
        }
        assert!(!measurements[0].is_iono_free); // missing L2
    }
}
