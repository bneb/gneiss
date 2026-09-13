use gneiss_core::imu::ImuMeasurement;
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

use super::types::{
    earth_rotation_rate_ecef, normal_gravity_ecef, skew_symmetric, EngineError, EskfState,
    Matrix15, Vector15,
};

fn fill_position_transition(
    phi: &mut Matrix15<f64>,
    r_b2e: &Matrix3<f64>,
    f_e_skew: &Matrix3<f64>,
    dt: f64,
) {
    let half_dt2 = 0.5 * dt * dt;
    for r in 0..3 {
        phi[(r, r + 3)] = dt;
        for c in 0..3 {
            phi[(r, c + 6)] = 0.5 * f_e_skew[(r, c)] * dt * dt;
            phi[(r, c + 9)] = -r_b2e[(r, c)] * half_dt2;
        }
    }
}

fn fill_velocity_transition(
    phi: &mut Matrix15<f64>,
    r_b2e: &Matrix3<f64>,
    f_e_skew: &Matrix3<f64>,
    omega_skew: &Matrix3<f64>,
    dt: f64,
) {
    let eye = Matrix3::<f64>::identity();
    for r in 0..3 {
        for c in 0..3 {
            let coriolis = eye[(r, c)] - 2.0 * omega_skew[(r, c)] * dt;
            phi[(r + 3, c + 3)] = coriolis;
            // MANDATORY POSITIVE SIGN: vel_att = +f_e_skew * dt
            phi[(r + 3, c + 6)] = f_e_skew[(r, c)] * dt;
            phi[(r + 3, c + 9)] = -r_b2e[(r, c)] * dt;
        }
    }
}

fn fill_attitude_transition(
    phi: &mut Matrix15<f64>,
    r_b2e: &Matrix3<f64>,
    omega_skew: &Matrix3<f64>,
    dt: f64,
) {
    let eye = Matrix3::<f64>::identity();
    for r in 0..3 {
        for c in 0..3 {
            let att_dyn = eye[(r, c)] - omega_skew[(r, c)] * dt;
            phi[(r + 6, c + 6)] = att_dyn;
            phi[(r + 6, c + 12)] = -r_b2e[(r, c)] * dt;
        }
    }
}

pub fn compute_transition_matrix(
    r_b2e: &Matrix3<f64>,
    f_e: &Vector3<f64>,
    dt: f64,
) -> Matrix15<f64> {
    let mut phi = Matrix15::identity();
    let f_e_skew = skew_symmetric(f_e);
    let omega_skew = skew_symmetric(&earth_rotation_rate_ecef());

    fill_position_transition(&mut phi, r_b2e, &f_e_skew, dt);
    fill_velocity_transition(&mut phi, r_b2e, &f_e_skew, &omega_skew, dt);
    fill_attitude_transition(&mut phi, r_b2e, &omega_skew, dt);
    phi
}

pub fn compute_process_noise(q_diag: &Vector15<f64>, dt: f64) -> Matrix15<f64> {
    let mut q = Matrix15::zeros();
    let dt2 = dt * dt;
    let dt3 = dt2 * dt;

    for i in 0..3 {
        let q_a = q_diag[i + 3];
        q[(i, i)] = (1.0 / 3.0) * q_a * dt3 + q_diag[i] * dt;
        let cross = 0.5 * q_a * dt2;
        q[(i, i + 3)] = cross;
        q[(i + 3, i)] = cross;
        q[(i + 3, i + 3)] = q_a * dt;
        q[(i + 6, i + 6)] = q_diag[i + 6] * dt;
        q[(i + 9, i + 9)] = q_diag[i + 9] * dt;
        q[(i + 12, i + 12)] = q_diag[i + 12] * dt;
    }
    q
}

pub fn propagate_nominal_state(
    state: &mut EskfState,
    accel_meas: &Vector3<f64>,
    gyro_meas: &Vector3<f64>,
    dt: f64,
) -> (Matrix3<f64>, Vector3<f64>) {
    let a_b = accel_meas - state.accel_bias;
    let omega_b = gyro_meas - state.gyro_bias;
    let r_b2e = state.attitude.to_rotation_matrix().into_inner();
    let f_e = r_b2e * a_b;
    let g_e = normal_gravity_ecef(&state.pos_ecef);
    let omega_ie = earth_rotation_rate_ecef();
    let a_coriolis = -2.0 * omega_ie.cross(&state.vel_ecef);
    let a_total = f_e + g_e + a_coriolis;

    state.pos_ecef += state.vel_ecef * dt + 0.5 * a_total * dt * dt;
    state.vel_ecef += a_total * dt;
    let dq = UnitQuaternion::from_scaled_axis(omega_b * dt);
    state.attitude = UnitQuaternion::new_normalize((state.attitude * dq).into_inner());
    (r_b2e, f_e)
}

pub fn predict(
    state: &mut EskfState,
    imu: &ImuMeasurement,
    dt: f64,
    q_diag: &Vector15<f64>,
) -> Result<(), EngineError> {
    predict_with_phi(state, &imu.accel, &imu.gyro, dt, q_diag).map(|_| ())
}

pub fn predict_with_phi(
    state: &mut EskfState,
    accel: &Vector3<f64>,
    gyro: &Vector3<f64>,
    dt: f64,
    q_diag: &Vector15<f64>,
) -> Result<Matrix15<f64>, EngineError> {
    let dt_eff = dt.max(1e-5);
    let (r_b2e, f_e) = propagate_nominal_state(state, accel, gyro, dt_eff);
    let phi = compute_transition_matrix(&r_b2e, &f_e, dt_eff);
    let q = compute_process_noise(q_diag, dt_eff);

    state.cov = phi * state.cov * phi.transpose() + q;
    state.cov = 0.5 * (state.cov + state.cov.transpose());
    Ok(phi)
}

pub fn predict_preintegrated(
    state: &mut EskfState,
    dp: &Vector3<f64>,
    dv: &Vector3<f64>,
    dq: &UnitQuaternion<f64>,
    dt: f64,
    q_diag: &Vector15<f64>,
) -> Result<Matrix15<f64>, EngineError> {
    let dt_eff = dt.max(1e-4);
    let r_b2e = state.attitude.to_rotation_matrix().into_inner();
    let f_e = r_b2e * (dv / dt_eff);
    let g_e = normal_gravity_ecef(&state.pos_ecef);
    let omega_ie = earth_rotation_rate_ecef();
    let a_coriolis = -2.0 * omega_ie.cross(&state.vel_ecef);

    state.pos_ecef += state.vel_ecef * dt_eff + 0.5 * (g_e + a_coriolis) * dt_eff * dt_eff + r_b2e * dp;
    state.vel_ecef += (g_e + a_coriolis) * dt_eff + r_b2e * dv;
    state.attitude = UnitQuaternion::new_normalize((state.attitude * dq).into_inner());

    let phi = compute_transition_matrix(&r_b2e, &f_e, dt_eff);
    let q = compute_process_noise(q_diag, dt_eff);
    state.cov = phi * state.cov * phi.transpose() + q;
    state.cov = 0.5 * (state.cov + state.cov.transpose());
    Ok(phi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transition_matrix_velocity_attitude_coupling_is_strictly_positive() {
        let r = Matrix3::identity();
        let f_e = Vector3::new(1.0, 2.0, 3.0);
        let dt = 0.1;
        let phi = compute_transition_matrix(&r, &f_e, dt);

        // f_e_skew has (0, 1) = -f_z = -3, (1, 0) = +f_z = +3, (0, 2) = +f_y = +2
        let f_e_skew = skew_symmetric(&f_e);
        for row in 0..3 {
            for col in 0..3 {
                let expected = f_e_skew[(row, col)] * dt;
                let actual = phi[(3 + row, 6 + col)];
                assert!((actual - expected).abs() < 1e-12);
            }
        }
        assert_eq!(phi[(3 + 1, 6)], 3.0 * dt);
    }

    #[test]
    fn test_process_noise_positive_semidefinite() {
        let mut q_diag = Vector15::zeros();
        for i in 0..15 {
            q_diag[i] = 0.1;
        }
        let q = compute_process_noise(&q_diag, 0.02);
        assert!((q - q.transpose()).norm() < 1e-12);
        for i in 0..15 {
            assert!(q[(i, i)] > 0.0);
        }
    }

    #[test]
    fn test_predict_stationary_preserves_position_under_reaction_accel() {
        let pos = Vector3::new(6378137.0, 0.0, 0.0);
        let vel = Vector3::zeros();
        let att = UnitQuaternion::identity();
        let mut state = EskfState::new(pos, vel, att);

        // Body accelerometer measures specific force (upward reaction against gravity)
        let g = normal_gravity_ecef(&pos);
        let accel_meas = -g; // specific force cancels gravity
        let gyro_meas = earth_rotation_rate_ecef(); // co-rotating with Earth
        let imu = ImuMeasurement::new(0, accel_meas, gyro_meas);

        let q_diag = Vector15::from_element(1e-6);
        let dt = 0.01;
        predict(&mut state, &imu, dt, &q_diag).expect("predict failed");

        assert!((state.pos_ecef - pos).norm() < 1e-4);
        assert!(state.vel_ecef.norm() < 1e-3);
    }

    #[test]
    fn test_predict_preintegrated_linear_motion() {
        let pos = Vector3::new(6378137.0, 0.0, 0.0);
        let vel = Vector3::new(10.0, 0.0, 0.0);
        let mut state = EskfState::new(pos, vel, UnitQuaternion::identity());

        let dt = 0.1;
        let dp = Vector3::zeros();
        let dv = Vector3::zeros();
        let dq = UnitQuaternion::identity();
        let q_diag = Vector15::from_element(1e-6);

        predict_preintegrated(&mut state, &dp, &dv, &dq, dt, &q_diag).expect("preintegrated failed");
        assert!((state.pos_ecef.x - (pos.x + 1.0)).abs() < 0.1);
    }
}
