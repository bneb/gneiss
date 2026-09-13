use nalgebra::{Matrix2, Matrix2x3, Matrix3, SMatrix, Vector2, Vector3};

use super::types::{skew_symmetric, EngineError, EskfState};
use super::update::{apply_error_injection, joseph_form_update};

pub type Matrix2x15<T = f64> = SMatrix<T, 2, 15>;
pub type Matrix3x15<T = f64> = SMatrix<T, 3, 15>;

fn build_nhc_system(
    state: &EskfState,
    _lever_arm: &Vector3<f64>,
) -> (Vector2<f64>, Matrix2x15<f64>) {
    let r_e2b = state.attitude.to_rotation_matrix().into_inner().transpose();
    let v_b = r_e2b * state.vel_ecef;
    // Lateral (Y) and vertical (Z) body velocity constraint
    let y = Vector2::new(-v_b.y, -v_b.z);

    let e23 = Matrix2x3::new(
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    );
    let h_vel = e23 * r_e2b;
    let mut h = Matrix2x15::zeros();
    h.fixed_view_mut::<2, 3>(0, 3).copy_from(&h_vel);

    if state.vel_ecef.norm() > 0.5 {
        let v_skew = skew_symmetric(&state.vel_ecef);
        let h_att = e23 * r_e2b * v_skew;
        h.fixed_view_mut::<2, 3>(0, 6).copy_from(&h_att);
    }
    (y, h)
}

pub fn update_nhc(
    state: &mut EskfState,
    lever_arm: &Vector3<f64>,
    r_nhc: &Matrix2<f64>,
) -> Result<(), EngineError> {
    let (y, h) = build_nhc_system(state, lever_arm);
    let s = h * state.cov * h.transpose() + r_nhc;
    let s_inv = s.try_inverse().ok_or(EngineError::InversionError)?;

    let k = state.cov * h.transpose() * s_inv;
    let mut dx = k * y;
    dx.fixed_rows_mut::<6>(9).fill(0.0);

    apply_error_injection(state, &dx);
    state.cov = joseph_form_update(&state.cov, &h, &k, r_nhc);
    Ok(())
}

fn build_zupt_system(state: &EskfState) -> (Vector3<f64>, Matrix3x15<f64>) {
    let y = -state.vel_ecef;
    let mut h = Matrix3x15::zeros();
    for i in 0..3 {
        h[(i, i + 3)] = 1.0;
    }
    (y, h)
}

pub fn update_zupt(state: &mut EskfState, r_zupt: &Matrix3<f64>) -> Result<(), EngineError> {
    let (y, h) = build_zupt_system(state);
    let s = h * state.cov * h.transpose() + r_zupt;
    let s_inv = s.try_inverse().ok_or(EngineError::InversionError)?;

    let k = state.cov * h.transpose() * s_inv;
    let mut dx = k * y;
    dx.fixed_rows_mut::<6>(9).fill(0.0);

    apply_error_injection(state, &dx);
    state.cov = joseph_form_update(&state.cov, &h, &k, r_zupt);
    Ok(())
}

fn build_body_vel_system(
    state: &EskfState,
    v_body_meas: &Vector3<f64>,
) -> (Vector3<f64>, Matrix3x15<f64>) {
    let r_e2b = state.attitude.to_rotation_matrix().into_inner().transpose();
    let v_b = r_e2b * state.vel_ecef;
    let y = v_body_meas - v_b;

    let mut h = Matrix3x15::zeros();
    h.fixed_view_mut::<3, 3>(0, 3).copy_from(&r_e2b);

    if state.vel_ecef.norm() > 0.5 {
        let v_skew = skew_symmetric(&state.vel_ecef);
        let h_att = r_e2b * v_skew;
        h.fixed_view_mut::<3, 3>(0, 6).copy_from(&h_att);
    }
    (y, h)
}

pub fn update_body_velocity(
    state: &mut EskfState,
    v_body_meas: &Vector3<f64>,
    r_v: &Matrix3<f64>,
) -> Result<(), EngineError> {
    let (y, h) = build_body_vel_system(state, v_body_meas);
    let s = h * state.cov * h.transpose() + r_v;
    let s_inv = s.try_inverse().ok_or(EngineError::InversionError)?;

    let k = state.cov * h.transpose() * s_inv;
    let mut dx = k * y;
    dx.fixed_rows_mut::<6>(9).fill(0.0);

    apply_error_injection(state, &dx);
    state.cov = joseph_form_update(&state.cov, &h, &k, r_v);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::UnitQuaternion;

    #[test]
    fn test_nhc_attitude_coupling_reduces_heading_misalignment() {
        // Vehicle moving forward along body X, but attitude is rotated 5 deg around Z
        let forward_vel_b = Vector3::new(15.0, 0.0, 0.0);
        let yaw_err = 5.0_f64.to_radians();
        let att_misaligned = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), yaw_err);
        // True velocity in ECEF is along body X of the true attitude (identity)
        let vel_e = forward_vel_b;

        let mut state = EskfState::new(Vector3::zeros(), vel_e, att_misaligned);
        let r_nhc = Matrix2::from_diagonal(&Vector2::new(0.01, 0.01));

        update_nhc(&mut state, &Vector3::zeros(), &r_nhc).expect("NHC update failed");

        let angle_after = state.attitude.angle();
        assert!(angle_after < yaw_err, "NHC should reduce attitude error, was {yaw_err} rad, now {angle_after} rad");
    }

    #[test]
    fn test_zupt_locks_velocity_to_zero() {
        let vel = Vector3::new(0.8, -0.5, 0.2);
        let mut state = EskfState::new(Vector3::zeros(), vel, UnitQuaternion::identity());
        let r_zupt = Matrix3::from_diagonal(&Vector3::new(1e-4, 1e-4, 1e-4));

        update_zupt(&mut state, &r_zupt).expect("ZUPT failed");
        assert!(state.vel_ecef.norm() < 0.05);
        assert!(state.cov[(3, 3)] < 1e-3);
    }

    #[test]
    fn test_body_velocity_update() {
        let vel_b = Vector3::new(10.0, 0.0, 0.0);
        let mut state = EskfState::new(
            Vector3::zeros(),
            Vector3::new(5.0, 1.0, -1.0),
            UnitQuaternion::identity(),
        );
        let r_v = Matrix3::from_diagonal(&Vector3::new(0.01, 0.01, 0.01));

        update_body_velocity(&mut state, &vel_b, &r_v).expect("Body vel update failed");
        assert!(state.vel_ecef.x > 8.5, "velocity should be pulled toward 10.0");
        assert!(state.cov[(3, 3)] < 0.05, "velocity covariance should decrease");
    }
}
