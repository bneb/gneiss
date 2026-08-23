//! State & Covariance Time Propagation for DD-IEKF.

use nalgebra::DMatrix;
use gneiss_core::time::GpsTime;
use super::state::RtkState;

/// Predict state and covariance to target GPS time.
/// Returns the state transition matrix F for smoothing.
pub fn predict_state(state: &mut RtkState, target_time: GpsTime, q_accel: f64) -> DMatrix<f64> {
    let dt = target_time.tow - state.time.tow;
    let dim = state.dim();
    let f = build_transition_matrix(dim, dt);
    let q = build_process_noise(dim, dt, q_accel);

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
fn build_process_noise(dim: usize, dt: f64, q_accel: f64) -> DMatrix<f64> {
    let mut q = DMatrix::zeros(dim, dim);
    let dt_abs = dt.abs();
    let dt2 = dt_abs * dt_abs;
    let dt3 = dt2 * dt_abs;

    let q_pos = dt3 / 3.0 * q_accel;
    let q_pos_vel = dt2 / 2.0 * q_accel;
    let q_vel = dt_abs * q_accel;

    for i in 0..3 {
        q[(i, i)] = q_pos;
        q[(i, i + 3)] = if dt >= 0.0 { q_pos_vel } else { -q_pos_vel };
        q[(i + 3, i)] = if dt >= 0.0 { q_pos_vel } else { -q_pos_vel };
        q[(i + 3, i + 3)] = q_vel;
    }

    // Small process noise on carrier phase ambiguities to prevent singular covariance
    for i in 6..dim {
        q[(i, i)] = 1e-7 * dt_abs;
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
}
