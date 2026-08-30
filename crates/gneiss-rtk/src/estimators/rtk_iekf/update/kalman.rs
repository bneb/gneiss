//! Iterated Extended Kalman Filter (IEKF) step computation and Joseph covariance update.

use nalgebra::{DMatrix, DVector};
use crate::estimators::rtk_iekf::state::RtkState;
use super::system::build_measurement_system;
use super::DoubleDiffMeasurement;

/// Compute a single IEKF state innovation step via fast Cholesky inverse.
pub fn compute_iekf_step(
    p0: &DMatrix<f64>,
    h: &DMatrix<f64>,
    y: &DVector<f64>,
    r: &DMatrix<f64>,
    x: &DVector<f64>,
    x0: &DVector<f64>,
) -> Option<DVector<f64>> {
    let ht = h.transpose();
    let hp = h * p0;
    let s = &hp * &ht + r;
    let dx_iter = x - x0;
    let innov = y + h * dx_iter;

    let s_inv = s.clone().cholesky().map(|c| c.inverse()).or_else(|| s.try_inverse())?;
    let k = p0 * &ht * &s_inv;
    Some(&k * innov)
}

/// Apply stabilized Joseph-form covariance update.
pub fn apply_joseph_update(
    state: &mut RtkState,
    p0: &DMatrix<f64>,
    x: &DVector<f64>,
    measurements: &[DoubleDiffMeasurement],
    innov_gate_scale: f64,
) {
    let (h_final, y_final, r_final) = build_measurement_system(state, x, measurements, innov_gate_scale);
    if y_final.is_empty() {
        return;
    }
    let ht = h_final.transpose();
    let hp = &h_final * p0;
    let s = &hp * &ht + &r_final;
    let s_inv = s.clone().cholesky().map(|c| c.inverse()).or_else(|| s.try_inverse());

    if let Some(s_inv_mat) = s_inv {
        let k = p0 * &ht * &s_inv_mat;
        let i_kh = DMatrix::identity(state.dim(), state.dim()) - &k * &h_final;
        state.cov = &i_kh * p0 * i_kh.transpose() + &k * &r_final * k.transpose();
    }
}
