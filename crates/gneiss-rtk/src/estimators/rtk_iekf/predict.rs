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
    // Per-pair iono residual: slow random walk.
    const IONO_RW_M2_PER_S: f64 = 1e-8;
    let io_off = state.iono_offset();
    for i in 0..state.ionos.len() {
        let idx = io_off + i;
        if idx < q.nrows() {
            q[(idx, idx)] += IONO_RW_M2_PER_S * dt_abs;
        }
    }
    // Satellite-mapped slant iono: slow random walk.
    const SAT_IONO_RW_M2_PER_S: f64 = 1e-8;
    let so_off = state.sat_iono_offset();
    for i in 0..state.sat_ionos.len() {
        let idx = so_off + i;
        if idx < q.nrows() {
            q[(idx, idx)] += SAT_IONO_RW_M2_PER_S * dt_abs;
        }
    }

    q
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use crate::post_process::dynamics::{ProcessingDynamics, STATIC_Q_ACCEL, KINEMATIC_Q_ACCEL};

    fn q_for(q_accel: f64, dt: f64) -> DMatrix<f64> {
        let state = RtkState::new(Vector3::zeros(), GpsTime::new(2000, 0.0));
        build_process_noise(&state, 6, dt, q_accel, false)
    }

    #[test]
    fn env_unset_profile_selects_static_legacy_q_bitwise() {
        // Env unset -> Static -> the exact legacy constants, so every
        // Q element produced through profile resolution is bitwise
        // identical to the historical hard-coded path.
        let dynamics = ProcessingDynamics::from_env_value(None);
        assert_eq!(dynamics, ProcessingDynamics::Static);
        let via_profile = q_for(dynamics.q_accel_default(), 30.0);
        let legacy_binary = q_for(1e-6_f64, 30.0);
        for i in 0..6 {
            for j in 0..6 {
                assert_eq!(
                    via_profile[(i, j)].to_bits(),
                    legacy_binary[(i, j)].to_bits(),
                    "Q[{i}][{j}] must be bit-identical when env unset"
                );
            }
        }
    }

    #[test]
    fn kinematic_q_position_at_least_100x_static_and_velocity_nonzero() {
        let dt = 30.0;
        let q_static = q_for(STATIC_Q_ACCEL, dt);
        let q_kin = q_for(KINEMATIC_Q_ACCEL, dt);
        for i in 0..3 {
            assert!(
                q_kin[(i, i)] >= 100.0 * q_static[(i, i)],
                "position process noise axis {i}: kin {} vs static {}",
                q_kin[(i, i)], q_static[(i, i)]
            );
            // Velocity states carry live random-walk noise under both
            // profiles, but only kinematic makes them meaningful.
            assert!(q_static[(i + 3, i + 3)] > 0.0);
            assert!(q_kin[(i + 3, i + 3)] > 100.0 * q_static[(i + 3, i + 3)]);
            assert!((q_kin[(i + 3, i + 3)] - dt * KINEMATIC_Q_ACCEL).abs() < 1e-12,
                "velocity Q must equal dt*q_accel");
        }
        // Documented magnitude: dt^3/3 * q at 30 s epochs ~ 55 m std.
        let expected_pos_q = dt * dt * dt / 3.0 * KINEMATIC_Q_ACCEL;
        assert!((q_kin[(0, 0)] - expected_pos_q).abs() < 1e-9);
    }

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
