//! Iterated Extended Kalman Filter (IEKF) Double-Difference Measurement Update.

pub mod kalman;
pub mod robust;
pub mod system;

#[cfg(test)]
mod tests;

use nalgebra::Vector3;
use crate::estimators::rtk_iekf::state::{DoubleDiffKey, RtkState};
use self::kalman::{apply_joseph_update, compute_iekf_step};
use self::system::build_measurement_system;

pub use self::robust::{
    if_residual_outliers, phase_innovation_outliers, robust_inflate, update_zwd_scalar,
    GRAD_INIT_VAR_M2, GRAD_MIN_SIN_EL, GRAD_RW_M2_PER_S, MAX_ZWD_STEP_M,
    PHASE_INNOVATION_GATE_CYCLES, ROBUST_INNOVATION_THRESHOLD, ZWD_INIT_VAR_M2, ZWD_RW_M2_PER_S,
};
pub use self::system::compute_tropo_dd;

/// Double-difference observation for a satellite pair on a single frequency band.
#[derive(Debug, Clone)]
pub struct DoubleDiffMeasurement {
    pub key: DoubleDiffKey,
    pub dd_pr_m: f64,
    pub dd_cp_cycles: Option<f64>,
    pub sat_pos: Vector3<f64>,
    pub ref_pos: Vector3<f64>,
    pub base_pos: Vector3<f64>,
    pub lambda: f64,
    pub pr_var_m2: f64,
    pub cp_var_cycles2: f64,
    pub dm_wet_rov: f64,
    pub dgrad_n_rov: f64,
    pub dgrad_e_rov: f64,
    pub tide_dd_m: f64,
    pub dd_pcv_m: f64,
}

/// DD phase observation with differential receiver-PCV removed.
pub fn pcv_corrected_cp(m: &DoubleDiffMeasurement) -> Option<f64> {
    m.dd_cp_cycles.map(|cp| cp - m.dd_pcv_m / m.lambda)
}

/// Perform Iterated Extended Kalman Filter (IEKF) update with double-differenced measurements.
pub fn iekf_update(state: &mut RtkState, measurements: &[DoubleDiffMeasurement]) -> Result<f64, String> {
    iekf_update_gated(state, measurements, 1.0)
}

/// IEKF update with a profile-dependent robust-gate scale.
pub fn iekf_update_gated(
    state: &mut RtkState,
    measurements: &[DoubleDiffMeasurement],
    innov_gate_scale: f64,
) -> Result<f64, String> {
    if measurements.is_empty() {
        return Ok(0.0);
    }

    let p0 = state.cov.clone();
    let x0 = state.to_dvector();
    let mut x = x0.clone();
    let mut last_rms = 0.0;

    for _iter in 0..4 {
        let (h, y, r) = build_measurement_system(state, &x, measurements, innov_gate_scale);
        if y.is_empty() {
            return Ok(0.0);
        }
        let delta_x = match compute_iekf_step(&p0, &h, &y, &r, &x, &x0) {
            Some(dx) => dx,
            None => return Err("Innovation covariance matrix inversion failed".to_string()),
        };
        x = &x0 + &delta_x;
        last_rms = delta_x.rows(0, 3).norm();
        if last_rms < 1e-4 {
            break;
        }
    }

    apply_joseph_update(state, &p0, &x, measurements, innov_gate_scale);
    state.update_from_dvector(&x);
    Ok(last_rms)
}
