//! RAIM fault detection, outlier exclusion, and initial state seeding for SPP.

use gneiss_core::coords::{ecef_to_llh, Coordinate, Datum, Frame};
use gneiss_core::atmosphere::KlobucharParams;
use nalgebra::Vector3;
use super::measurements::compute_sat_state;
use super::solver::{compute_measurement_residuals, solve_spp_iteratively};
use super::{SppConfig, SppError, SppMeasurement, SppState};

pub(crate) fn compute_seed_position(measurements: &[SppMeasurement]) -> (f64, f64, f64) {
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

pub(crate) fn compute_seed_clocks(
    measurements: &[SppMeasurement],
    sx: f64,
    sy: f64,
    sz: f64,
) -> (f64, f64, f64, f64) {
    let mut cdt_gps = None;
    let mut cdt_gal = None;
    let mut cdt_bds = None;
    let mut cdt_glo = None;
    for m in measurements {
        let (sat_coord, corrected_pr) = compute_sat_state(m, 0.0);
        let dx = sx - sat_coord.vector.x;
        let dy = sy - sat_coord.vector.y;
        let dz = sz - sat_coord.vector.z;
        let cdt = corrected_pr - f64::sqrt(dx * dx + dy * dy + dz * dz);
        match m.constellation {
            gneiss_core::sat::Constellation::Gps => {
                if cdt_gps.is_none() {
                    cdt_gps = Some(cdt);
                }
            }
            gneiss_core::sat::Constellation::Galileo => {
                if cdt_gal.is_none() {
                    cdt_gal = Some(cdt);
                }
            }
            gneiss_core::sat::Constellation::Beidou => {
                if cdt_bds.is_none() {
                    cdt_bds = Some(cdt);
                }
            }
            gneiss_core::sat::Constellation::Glonass if cdt_glo.is_none() => {
                cdt_glo = Some(cdt);
            }
            _ => {}
        }
    }
    let default_cdt = cdt_gps.or(cdt_gal).or(cdt_bds).or(cdt_glo).unwrap_or(0.0);
    (
        cdt_gps.unwrap_or(default_cdt),
        cdt_gal.unwrap_or(default_cdt),
        cdt_bds.unwrap_or(default_cdt),
        cdt_glo.unwrap_or(default_cdt),
    )
}

pub(crate) fn seed_initial_state(measurements: &[SppMeasurement], prev_state: Option<&SppState>) -> SppState {
    let (seed_x, seed_y, seed_z) = if let Some(coord) = prev_state.map(|s| &s.position) {
        (coord.vector.x, coord.vector.y, coord.vector.z)
    } else {
        compute_seed_position(measurements)
    };
    let (cdt, cdt_gal, cdt_bds, cdt_glo) =
        compute_seed_clocks(measurements, seed_x, seed_y, seed_z);
    SppState::new(
        Coordinate::new(
            Vector3::new(seed_x, seed_y, seed_z),
            Datum::WGS84,
            Frame::ECEF,
            measurements[0].time,
        ),
        cdt,
        cdt_gal,
        cdt_bds,
        cdt_glo,
    )
}

pub(crate) fn apply_raim(
    state: SppState,
    seed_state: SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> Result<SppState, SppError> {
    let good_measurements = filter_raim_outliers(&state, measurements, iono_params, config);
    if good_measurements.len() < measurements.len() && good_measurements.len() >= 4 {
        return solve_spp_iteratively(seed_state, &good_measurements, iono_params, config);
    }
    Ok(state)
}

pub(crate) fn filter_raim_outliers(
    state: &SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> Vec<SppMeasurement> {
    let rec_ecef = state.position.vector;
    let rec_llh = ecef_to_llh(rec_ecef);

    let mut residuals: Vec<(&SppMeasurement, f64)> = measurements
        .iter()
        .map(|m| {
            let (_, _, _, _, residual, _) =
                compute_measurement_residuals(state, m, rec_ecef, rec_llh, iono_params, config);
            (m, residual.abs())
        })
        .collect();

    residuals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let median_residual = residuals[residuals.len() / 2].1;
    let mut deviations: Vec<f64> = residuals
        .iter()
        .map(|(_, r)| (*r - median_residual).abs())
        .collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad_threshold =
        (deviations[deviations.len() / 2] * config.raim_mad_multiplier).max(config.raim_outlier_m);

    residuals
        .into_iter()
        .filter(|(_, r)| *r <= mad_threshold)
        .map(|(m, _)| m.clone())
        .collect()
}
