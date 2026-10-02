//! WNLLS solver iteration, design matrix building, and atmospheric corrections for SPP.

use nalgebra::{DMatrix, DVector, Vector3};
use gneiss_core::atmosphere::{AtmosphereModel, KlobucharParams, TropoParams};
use gneiss_core::coords::{az_el, ecef_to_llh, Coordinate, Datum, Frame};
use super::measurements::compute_sat_state;
use super::{
    SppConfig, SppError, SppMeasurement, SppState, HEIGHT_CONSTRAINT_VAR, LIGHT_SPEED,
    MIN_ATMOSPHERE_ELEVATION_RAD, MIN_EARTH_RADIUS_M, MIN_WEIGHT, OMEGA_E,
};

pub(crate) struct ClockCols(pub Option<usize>, pub Option<usize>, pub Option<usize>, pub Option<usize>);

pub fn solve_spp_iteratively(
    mut state: SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
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

/// Performs a single Weighted Non-Linear Least Squares (WNLLS) iteration.
pub fn spp_wnlls_step(
    current_state: &SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> Result<SppState, SppError> {
    let (mut h_matrix, mut w_matrix, mut dz_vector, cols, clocks) =
        build_design_matrix(current_state, measurements, iono_params, config)?;
    let n = measurements.len();
    if n == cols - 1 {
        let rec_llh = ecef_to_llh(Vector3::new(
            current_state.position.vector.x,
            current_state.position.vector.y,
            current_state.position.vector.z,
        ));
        apply_height_constraint(rec_llh, n, &mut h_matrix, &mut w_matrix, &mut dz_vector);
    }

    let h_t = h_matrix.transpose();
    let h_t_w = &h_t * &w_matrix;
    let h_t_w_h_inv_opt = (&h_t_w * &h_matrix).try_inverse();
    let h_t_w_h_inv = match h_t_w_h_inv_opt {
        Some(inv) => inv,
        None => {
            return Err(SppError::MatrixInversionFailed);
        }
    };

    let trace = h_t_w_h_inv[(0, 0)] + h_t_w_h_inv[(1, 1)] + h_t_w_h_inv[(2, 2)];
    if trace > config.geometry_variance_threshold {
        return Err(SppError::PoorGeometry);
    }

    let dx_vec = h_t_w_h_inv * h_t_w * dz_vector;

    Ok(SppState::new(
        Coordinate::new(
            Vector3::new(
                current_state.position.vector.x + dx_vec[0],
                current_state.position.vector.y + dx_vec[1],
                current_state.position.vector.z + dx_vec[2],
            ),
            Datum::WGS84,
            Frame::ECEF,
            measurements[0].time,
        ),
        current_state.cdt + clocks.0.map(|c| dx_vec[c]).unwrap_or(0.0),
        current_state.cdt_gal + clocks.1.map(|c| dx_vec[c]).unwrap_or(0.0),
        current_state.cdt_bds + clocks.2.map(|c| dx_vec[c]).unwrap_or(0.0),
        current_state.cdt_glo + clocks.3.map(|c| dx_vec[c]).unwrap_or(0.0),
    ))
}

pub(crate) fn find_clock_cols(measurements: &[SppMeasurement]) -> (usize, ClockCols) {
    let mut cols = 3;
    let has = |c| measurements.iter().any(|m| m.constellation == c);
    let gps_col = if has(gneiss_core::sat::Constellation::Gps)
        || has(gneiss_core::sat::Constellation::Qzss)
    {
        cols += 1;
        Some(cols - 1)
    } else {
        None
    };
    let gal_col = if has(gneiss_core::sat::Constellation::Galileo) {
        cols += 1;
        Some(cols - 1)
    } else {
        None
    };
    let bds_col = if has(gneiss_core::sat::Constellation::Beidou) {
        cols += 1;
        Some(cols - 1)
    } else {
        None
    };
    let glo_col = if has(gneiss_core::sat::Constellation::Glonass) {
        cols += 1;
        Some(cols - 1)
    } else {
        None
    };
    (cols, ClockCols(gps_col, gal_col, bds_col, glo_col))
}

#[allow(clippy::type_complexity)]
pub(crate) fn build_design_matrix(
    state: &SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> Result<(DMatrix<f64>, DMatrix<f64>, DVector<f64>, usize, ClockCols), SppError> {
    let (cols, clocks) = find_clock_cols(measurements);
    let n = measurements.len();
    if n < cols - 1 {
        return Err(SppError::NotEnoughMeasurements);
    }
    let matrix_n = if n == cols - 1 { n + 1 } else { n };

    let mut h_matrix = DMatrix::<f64>::zeros(matrix_n, cols);
    let mut w_matrix = DMatrix::<f64>::zeros(matrix_n, matrix_n);
    let mut dz_vector = DVector::<f64>::zeros(matrix_n);
    let rec_ecef = state.position.vector;
    let rec_llh = ecef_to_llh(rec_ecef);

    for (i, m) in measurements.iter().enumerate() {
        let (dx, dy, dz, r, residual, el) =
            compute_measurement_residuals(state, m, rec_ecef, rec_llh, iono_params, config);

        if el < config.elevation_mask_rad && rec_ecef.x != 0.0 {
            w_matrix[(i, i)] = MIN_WEIGHT;
            continue;
        }

        h_matrix[(i, 0)] = dx / r;
        h_matrix[(i, 1)] = dy / r;
        h_matrix[(i, 2)] = dz / r;
        if let Some(c) = match m.constellation {
            gneiss_core::sat::Constellation::Gps | gneiss_core::sat::Constellation::Qzss => {
                clocks.0
            }
            gneiss_core::sat::Constellation::Galileo => clocks.1,
            gneiss_core::sat::Constellation::Beidou => clocks.2,
            gneiss_core::sat::Constellation::Glonass => clocks.3,
            _ => None,
        } {
            h_matrix[(i, c)] = 1.0;
        }

        w_matrix[(i, i)] = 1.0
            / gneiss_core::variance::observation_variance(m.snr, el, config.snr_a, config.snr_b);
        dz_vector[i] = residual;
    }

    Ok((h_matrix, w_matrix, dz_vector, cols, clocks))
}

pub(crate) fn compute_measurement_residuals(
    current_state: &SppState,
    m: &SppMeasurement,
    rec_ecef: Vector3<f64>,
    rec_llh: Vector3<f64>,
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> (f64, f64, f64, f64, f64, f64) {
    let cdt = current_state.get_cdt(m.constellation);
    let (sat_coord, corrected_pr) = compute_sat_state(m, cdt);

    let sat_ecef = if config.enable_sagnac {
        compute_sagnac_correction(sat_coord.vector, corrected_pr - cdt)
    } else {
        sat_coord.vector
    };

    let dx = rec_ecef.x - sat_ecef.x;
    let dy = rec_ecef.y - sat_ecef.y;
    let dz = rec_ecef.z - sat_ecef.z;
    let r = f64::sqrt(dx * dx + dy * dy + dz * dz).max(1e-6);

    let (az, el) = az_el(rec_llh, rec_ecef, sat_ecef);
    let (tropo, mut iono) =
        compute_atmospheric_delays(rec_ecef, rec_llh, az, el, m.time, iono_params, config);
    if m.is_iono_free {
        iono = 0.0;
    }

    (dx, dy, dz, r, corrected_pr - (r + cdt + tropo + iono), el)
}

pub(crate) fn compute_sagnac_correction(sat_ecef: Vector3<f64>, geometric_pr: f64) -> Vector3<f64> {
    let tof = geometric_pr / LIGHT_SPEED;
    let theta = OMEGA_E * tof;
    let cos_t = f64::cos(theta);
    let sin_t = f64::sin(theta);
    Vector3::new(
        sat_ecef.x * cos_t + sat_ecef.y * sin_t,
        -sat_ecef.x * sin_t + sat_ecef.y * cos_t,
        sat_ecef.z,
    )
}

pub(crate) fn compute_atmospheric_delays(
    rec_ecef: Vector3<f64>,
    rec_llh: Vector3<f64>,
    az: f64,
    el: f64,
    time: gneiss_core::time::GpsTime,
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> (f64, f64) {
    if rec_ecef.norm() <= MIN_EARTH_RADIUS_M {
        return (0.0, 0.0);
    }
    let safe_el = el.max(MIN_ATMOSPHERE_ELEVATION_RAD);

    let tropo = if config.enable_tropo {
        AtmosphereModel::tropo_nmf(&TropoParams::default(), rec_llh, safe_el, time)
    } else {
        0.0
    };

    let iono = if config.enable_iono {
        if let Some(iono) = iono_params {
            AtmosphereModel::iono_klobuchar(iono, rec_llh, az, safe_el, time)
        } else {
            0.0
        }
    } else {
        0.0
    };

    (tropo, iono)
}

pub(crate) fn apply_height_constraint(
    rec_llh: Vector3<f64>,
    row_idx: usize,
    h_matrix: &mut DMatrix<f64>,
    w_matrix: &mut DMatrix<f64>,
    dz_vector: &mut DVector<f64>,
) {
    let lat = rec_llh.x;
    let lon = rec_llh.y;

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_lon = lon.sin();
    let cos_lon = lon.cos();

    h_matrix[(row_idx, 0)] = cos_lat * cos_lon;
    h_matrix[(row_idx, 1)] = cos_lat * sin_lon;
    h_matrix[(row_idx, 2)] = sin_lat;

    dz_vector[row_idx] = 0.0;
    w_matrix[(row_idx, row_idx)] = 1.0 / HEIGHT_CONSTRAINT_VAR;
}
