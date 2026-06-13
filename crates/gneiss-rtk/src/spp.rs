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
    pub nominal_snr_dbhz: f64,
    pub min_measurements_init: usize,
    pub raim_mad_multiplier: f64,
    pub elevation_mask_rad: f64,
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
            min_measurements_init: 3,
            raim_mad_multiplier: 7.413, // 1.4826 * 5.0 sigma
            elevation_mask_rad: 0.1745, // ~10 degrees
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
        tracing::debug!("SPP seed: SAT={}, raw_pr={:.3}, dt_sat_m={:.3}, geom_r={:.3}, cdt={:.3}", m.eph.sat(), m.raw_pr, (corrected_pr - m.raw_pr), geom_r, cdt);
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

    if measurements.len() < config.min_measurements_init {
        tracing::error!("SPP failed: Only {} valid measurements. Need at least {}.", measurements.len(), config.min_measurements_init);
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
            
            let mad_threshold = (mad * config.raim_mad_multiplier).max(config.raim_outlier_m);

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
    let h_t_w_h_inv = (&h_t_w * &h_matrix).try_inverse().ok_or(SppError::MatrixInversionFailed)?;

    if (h_t_w_h_inv[(0, 0)] + h_t_w_h_inv[(1, 1)] + h_t_w_h_inv[(2, 2)]) > config.geometry_variance_threshold { 
        return Err(SppError::PoorGeometry); 
    }

    let dx_vec = h_t_w_h_inv * h_t_w * dz_vector;
    
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
    let gps_col = if has(gneiss_core::sat::Constellation::Gps) { cols += 1; Some(cols - 1) } else { None };
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

    let rec_ecef = Vector3::new(state.position.vector.x, state.position.vector.y, state.position.vector.z);
    let rec_llh = ecef_to_llh(rec_ecef);

    for (i, m) in measurements.iter().enumerate() {
        let (dx, dy, dz, r, residual, el) = compute_measurement_residuals(state, m, rec_ecef, rec_llh, iono_params, config);
        
        let el_mask = config.elevation_mask_rad;

        if el < el_mask && state.position.vector.x != 0.0 { w_matrix[(i, i)] = MIN_WEIGHT; continue; }

        h_matrix[(i, 0)] = dx / r; h_matrix[(i, 1)] = dy / r; h_matrix[(i, 2)] = dz / r;
        if m.constellation == gneiss_core::sat::Constellation::Gps { if let Some(c) = clocks.0 { h_matrix[(i, c)] = 1.0; } }
        else if m.constellation == gneiss_core::sat::Constellation::Galileo { if let Some(c) = clocks.1 { h_matrix[(i, c)] = 1.0; } }
        else if m.constellation == gneiss_core::sat::Constellation::Beidou { if let Some(c) = clocks.2 { h_matrix[(i, c)] = 1.0; } }
        
        w_matrix[(i, i)] = 1.0 / gneiss_core::variance::observation_variance(m.snr, el, config.nominal_snr_dbhz);
        dz_vector[i] = residual;
    }
    
    Ok((h_matrix, w_matrix, dz_vector, cols, clocks))
}

fn compute_measurement_residuals(
    current_state: &SppState,
    m: &SppMeasurement,
    rec_ecef: Vector3<f64>,
    rec_llh: Vector3<f64>,
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> (f64, f64, f64, f64, f64, f64) {
    let cdt = match m.constellation {
        gneiss_core::sat::Constellation::Gps => current_state.cdt,
        gneiss_core::sat::Constellation::Galileo => current_state.cdt_gal,
        gneiss_core::sat::Constellation::Beidou => current_state.cdt_bds,
        _ => current_state.cdt,
    };
    
    let (sat_coord, corrected_pr) = compute_sat_state(m, cdt);
    let sat_ecef = Vector3::new(sat_coord.vector.x, sat_coord.vector.y, sat_coord.vector.z);
    
    let sat_ecef_rot = if config.enable_sagnac {
        compute_sagnac_correction(sat_ecef, corrected_pr - cdt)
    } else { sat_ecef };

    let dx = current_state.position.vector.x - sat_ecef_rot.x;
    let dy = current_state.position.vector.y - sat_ecef_rot.y;
    let dz = current_state.position.vector.z - sat_ecef_rot.z;
    let r = f64::sqrt(dx * dx + dy * dy + dz * dz).max(1e-6);

    let (az, el) = az_el(rec_llh, rec_ecef, sat_ecef_rot);
    let (tropo_delay, iono_delay) = compute_atmospheric_delays(rec_ecef, rec_llh, az, el, m.time, iono_params, config);

    let expected_pr = r + cdt + tropo_delay + iono_delay;
    (dx, dy, dz, r, corrected_pr - expected_pr, el)
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
        };

        let config = SppConfig {
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };
        let state = SppState::new(Coordinate::new(nalgebra::Vector3::zeros(), Datum::WGS84, Frame::ECEF, t), 0.0, 0.0, 0.0, 0.0);
        
        let res = spp_wnlls_step(&state, &[m1], None, &config);
        assert_eq!(res.unwrap_err(), SppError::NotEnoughMeasurements);
    }

}
