//! State & Covariance Time Propagation for DD-IEKF.

use nalgebra::DMatrix;
use gneiss_core::time::GpsTime;
use super::state::RtkState;

/// Predict state and covariance to target GPS time.
/// Returns the state transition matrix F for smoothing.
pub fn predict_state(state: &mut RtkState, target_time: GpsTime, q_accel: f64) -> DMatrix<f64> {
    predict_state_gated(state, target_time, q_accel, false)
}

/// Predict with `reverse_safe` selecting corrected backward-time handling.
///
/// The legacy path builds Q from the signed dt: a negative dt makes the
/// position/velocity process-noise DIAGONALS negative (only the cross term
/// was sign-flipped), so every step of a long reverse pass bleeds the
/// covariance toward non-positive-definiteness until inversions explode.
/// With `reverse_safe` the noise magnitudes use |dt| — variances always
/// grow; the time direction is carried by F alone. Opt-in to keep legacy
/// datasets bit-identical.
pub fn predict_state_gated(
    state: &mut RtkState,
    target_time: GpsTime,
    q_accel: f64,
    reverse_safe: bool,
) -> DMatrix<f64> {
    let dt = target_time.tow - state.time.tow;
    let dim = state.dim();
    let f = build_transition_matrix(dim, dt);
    let q = build_process_noise(state, dim, dt, q_accel, reverse_safe && dt < 0.0);

    // Propagate state vector
    state.pos_ecef += state.vel_ecef * dt;
    state.time = target_time;

    // Propagate covariance: P = F * P * F^T + Q
    state.cov = &f * &state.cov * f.transpose() + q;

    f
}

/// Build state transition matrix F for dimension (6 + m).
fn build_transition_matrix(dim: usize, dt: f64) -> DMatrix<f64> {
    let mut f = DMatrix::identity(dim, dim);
    for i in 0..3 {
        f[(i, i + 3)] = dt;
    }
    f
}

/// Build process noise matrix Q for constant-velocity kinematic model.
fn build_process_noise(state: &RtkState, dim: usize, dt: f64, q_accel: f64, reverse_safe: bool) -> DMatrix<f64> {
    let mut q = DMatrix::zeros(dim, dim);
    let dt_abs = dt.abs();
    let dt2 = dt_abs * dt_abs;
    let dt3 = dt2 * dt_abs;
    // Backward propagation: magnitudes must mirror forward time exactly.
    let signed_dt3 = if reverse_safe { dt3 } else { dt3 * dt.signum() };
    let signed_dt = if reverse_safe { dt_abs } else { dt_abs * dt.signum() };

    let q_pos = signed_dt3 / 3.0 * q_accel;
    let q_pos_vel = dt2 / 2.0 * q_accel;
    let q_vel = signed_dt * q_accel;

    for i in 0..3 {
        q[(i, i)] = q_pos;
        q[(i, i + 3)] = if dt >= 0.0 || reverse_safe { q_pos_vel } else { -q_pos_vel };
        q[(i + 3, i)] = if dt >= 0.0 || reverse_safe { q_pos_vel } else { -q_pos_vel };
        q[(i + 3, i + 3)] = q_vel;
    }

    let off = state.amb_offset();
    let n_amb = state.ambiguities.len();
    for i in off..off + n_amb {
        q[(i, i)] = 1e-7 * dt_abs;
    }
    if let Some(zi) = state.zwd_idx() {
        q[(zi, zi)] = super::update::ZWD_RW_M2_PER_S * dt_abs;
    }
    if let Some((gn, ge)) = state.grad_idx() {
        q[(gn, gn)] = super::update::GRAD_RW_M2_PER_S * dt_abs;
        q[(ge, ge)] = super::update::GRAD_RW_M2_PER_S * dt_abs;
    }

    q
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_prediction_forward_and_backward() {
        let mut state = RtkState::new(Vector3::new(10.0, 20.0, 30.0), GpsTime::new(2000, 100.0));
        state.vel_ecef = Vector3::new(1.0, 2.0, 3.0);

        let fwd_time = GpsTime::new(2000, 101.0);
        let f = predict_state(&mut state, fwd_time, 1.0);
        assert_eq!(state.pos_ecef, Vector3::new(11.0, 22.0, 33.0));
        assert!(state.cov[(0, 0)] > 0.0);
        assert_eq!(f[(0, 3)], 1.0);

        let bwd_time = GpsTime::new(2000, 100.0);
        let f_bwd = predict_state(&mut state, bwd_time, 1.0);
        assert_eq!(state.pos_ecef, Vector3::new(10.0, 20.0, 30.0));
        assert_eq!(f_bwd[(0, 3)], -1.0);
    }

    #[test]
    fn test_reverse_safe_keeps_process_noise_positive() {
        // Legacy path: negative dt poisons the position variance directly.
        let mut legacy = DMatrix::<f64>::zeros(6, 6);
        let q_legacy = build_process_noise(&RtkState::new(Vector3::zeros(), GpsTime::new(2000,0.0)), 6, -30.0, 1.0, false);
        legacy += q_legacy;
        assert!(legacy[(0, 0)] < 0.0, "legacy behavior must stay observable");
        // Reverse-safe path: same |dt|, positive variances throughout. The
        // time direction lives in F (negative dt off-diagonal), so the Q
        // cross term must stay positive — flipping it here would
        // double-negate through F P F^T.
        let q_safe = build_process_noise(&RtkState::new(Vector3::zeros(), GpsTime::new(2000,0.0)), 6, -30.0, 1.0, true);
        assert!(q_safe[(0, 0)] > 0.0 && q_safe[(3, 3)] > 0.0);
        assert!(q_safe[(0, 3)] > 0.0);
        // Magnitudes match an equivalent forward step exactly.
        let q_fwd = build_process_noise(&RtkState::new(Vector3::zeros(), GpsTime::new(2000,0.0)), 6, 30.0, 1.0, true);
        assert!((q_safe[(0, 0)] - q_fwd[(0, 0)]).abs() < 1e-12);
        assert!((q_safe[(0, 3)] - q_fwd[(0, 3)]).abs() < 1e-12);
    }
}
