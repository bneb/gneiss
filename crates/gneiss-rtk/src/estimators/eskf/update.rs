use nalgebra::{Matrix3, SMatrix, SVector, UnitQuaternion, Vector3};

use super::types::{clamp_vector, skew_symmetric, EngineError, EskfState, Matrix15, Vector15};

pub type Vector6<T = f64> = SVector<T, 6>;
pub type Matrix6<T = f64> = SMatrix<T, 6, 6>;
pub type Matrix6x15<T = f64> = SMatrix<T, 6, 15>;

pub fn apply_error_injection(state: &mut EskfState, dx: &Vector15<f64>) {
    state.pos_ecef += dx.fixed_rows::<3>(0);
    state.vel_ecef += dx.fixed_rows::<3>(3);

    let d_theta = dx.fixed_rows::<3>(6).into_owned();
    let dq = UnitQuaternion::from_scaled_axis(d_theta);
    // Left-multiplied global frame attitude reset: q <- dq * q
    state.attitude = UnitQuaternion::new_normalize((dq * state.attitude).into_inner());

    let mut d_ba = dx.fixed_rows::<3>(9).into_owned();
    let mut d_bg = dx.fixed_rows::<3>(12).into_owned();
    clamp_vector(&mut d_ba, 0.05);
    clamp_vector(&mut d_bg, 0.002);

    state.accel_bias += d_ba;
    state.gyro_bias += d_bg;
    clamp_vector(&mut state.accel_bias, 0.5);
    clamp_vector(&mut state.gyro_bias, 0.05);
}

pub fn joseph_form_update<const M: usize>(
    cov: &Matrix15<f64>,
    h: &SMatrix<f64, M, 15>,
    k: &SMatrix<f64, 15, M>,
    r: &SMatrix<f64, M, M>,
) -> Matrix15<f64> {
    let i_kh = Matrix15::identity() - k * h;
    let p_new = i_kh * cov * i_kh.transpose() + k * r * k.transpose();
    0.5 * (p_new + p_new.transpose())
}

fn build_gnss_jacobian(l_e: &Vector3<f64>) -> Matrix6x15<f64> {
    let mut h = Matrix6x15::zeros();
    for i in 0..3 {
        h[(i, i)] = 1.0;
        h[(i + 3, i + 3)] = 1.0;
    }
    let l_e_skew = skew_symmetric(l_e);
    for r in 0..3 {
        for c in 0..3 {
            h[(r, c + 6)] = -l_e_skew[(r, c)];
        }
    }
    h
}

fn build_gnss_system(
    state: &EskfState,
    pos_meas: &Vector3<f64>,
    vel_meas: &Vector3<f64>,
    lever_arm: &Vector3<f64>,
    r_pos: &Matrix3<f64>,
    r_vel: &Matrix3<f64>,
) -> (Vector6<f64>, Matrix6x15<f64>, Matrix6<f64>) {
    let r_b2e = state.attitude.to_rotation_matrix().into_inner();
    let l_e = r_b2e * lever_arm;
    let mut y = Vector6::zeros();
    y.fixed_rows_mut::<3>(0).copy_from(&(pos_meas - (state.pos_ecef + l_e)));
    y.fixed_rows_mut::<3>(3).copy_from(&(vel_meas - state.vel_ecef));

    let h = build_gnss_jacobian(&l_e);
    let mut r_mat = Matrix6::zeros();
    r_mat.fixed_view_mut::<3, 3>(0, 0).copy_from(r_pos);
    r_mat.fixed_view_mut::<3, 3>(3, 3).copy_from(r_vel);
    (y, h, r_mat)
}

pub fn update_gnss_pos_vel(
    state: &mut EskfState,
    pos_meas: &Vector3<f64>,
    vel_meas: &Vector3<f64>,
    lever_arm: &Vector3<f64>,
    r_pos: &Matrix3<f64>,
    r_vel: &Matrix3<f64>,
) -> Result<(), EngineError> {
    let (y, h, r_mat) = build_gnss_system(state, pos_meas, vel_meas, lever_arm, r_pos, r_vel);
    let s = h * state.cov * h.transpose() + r_mat;
    let s_inv = s.try_inverse().ok_or(EngineError::InversionError)?;

    let k = state.cov * h.transpose() * s_inv;
    let dx = k * y;

    apply_error_injection(state, &dx);
    state.cov = joseph_form_update(&state.cov, &h, &k, &r_mat);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gnss_position_update_pulls_state_towards_measurement() {
        let pos_init = Vector3::new(100.0, 200.0, 300.0);
        let vel_init = Vector3::zeros();
        let att_init = UnitQuaternion::identity();
        let mut state = EskfState::new(pos_init, vel_init, att_init);

        let pos_meas = Vector3::new(101.0, 200.0, 300.0);
        let vel_meas = Vector3::zeros();
        let lever_arm = Vector3::zeros();
        let r_pos = Matrix3::from_diagonal(&Vector3::new(0.01, 0.01, 0.01));
        let r_vel = Matrix3::from_diagonal(&Vector3::new(0.01, 0.01, 0.01));

        update_gnss_pos_vel(&mut state, &pos_meas, &vel_meas, &lever_arm, &r_pos, &r_vel)
            .expect("update failed");

        assert!(state.pos_ecef.x > 100.5);
        assert!(state.cov[(0, 0)] < 1.0);
    }

    #[test]
    fn test_lever_arm_measurement_geometry() {
        let state = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
        let lever_arm = Vector3::new(0.0, 1.0, 0.0); // 1m in Y
        let r_pos = Matrix3::identity();
        let r_vel = Matrix3::identity();
        let (y, h, _) = build_gnss_system(&state, &Vector3::new(0.0, 1.0, 0.0), &Vector3::zeros(), &lever_arm, &r_pos, &r_vel);
        assert!(y.norm() < 1e-12);
        assert_eq!(h[(0, 6 + 2)], -1.0); // - [l_e x] term
    }
}
