use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion, Matrix3};

use crate::engine::{DynamicsModel, EngineConfig};

pub fn integrate_imu_mechanization(state: &mut RtkState, dt: f64, imu_buffer: &[gneiss_core::imu::ImuMeasurement]) {
    let imu_dt = dt / (imu_buffer.len() as f64);
    let omega_ie = Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
    
    for meas in imu_buffer {
        let f_b = meas.accel - state.accel_bias;
        let omega_b = meas.gyro - state.gyro_bias;
        
        let zeta = (omega_b - state.attitude.inverse() * omega_ie) * imu_dt;
        let angle = zeta.norm();
        let dq = if angle > 1e-12 {
            UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(zeta / angle), angle)
        } else {
            UnitQuaternion::identity()
        };

        let dq_mid = if angle > 1e-12 {
            UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(zeta / angle), angle * 0.5)
        } else {
            UnitQuaternion::identity()
        };
        let r_mid = state.attitude * dq_mid;
        
        state.attitude *= dq;
        state.attitude.renormalize();

        let f_e = r_mid * f_b;
        let gravity = gravity_wgs84(state.position.vector);
        let coriolis = 2.0 * omega_ie.cross(&state.velocity);
        let centrifugal = omega_ie.cross(&(omega_ie.cross(&state.position.vector)));
        
        let v_dot = f_e + gravity - coriolis - centrifugal;
        
        let v_mid = state.velocity + v_dot * (imu_dt * 0.5);
        state.velocity += v_dot * imu_dt;
        state.position.vector += v_mid * imu_dt;
    }
}

pub fn compute_transition_matrix(state: &RtkState, dt: f64, imu_buffer: &[gneiss_core::imu::ImuMeasurement]) -> DMatrix<f64> {
    let n = state.covariance.nrows();
    let mut phi = DMatrix::<f64>::identity(n, n);
    let omega_ie = Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
    
    if imu_buffer.is_empty() {
        phi[(0, 3)] = dt; phi[(1, 4)] = dt; phi[(2, 5)] = dt;
    } else {
        let r_b_e = state.attitude.to_rotation_matrix();
        let f_e = state.attitude * (imu_buffer.last().unwrap().accel - state.accel_bias);
        let f_e_skew = skew_symmetric(&f_e);
        let omega_ie_skew = skew_symmetric(&omega_ie);
        
        for i in 0..3 { phi[(i, 3 + i)] = dt; }
        
        let vel_att = -f_e_skew * dt;
        for r in 0..3 { for c in 0..3 { phi[(3 + r, 6 + c)] = vel_att[(r, c)]; } }
        
        let vel_vel = Matrix3::identity() - 2.0 * omega_ie_skew * dt;
        for r in 0..3 { for c in 0..3 { phi[(3 + r, 3 + c)] = vel_vel[(r, c)]; } }
        
        let vel_abias = -r_b_e.matrix() * dt;
        for r in 0..3 { for c in 0..3 { phi[(3 + r, 9 + c)] = vel_abias[(r, c)]; } }
        
        let att_att = Matrix3::identity() - omega_ie_skew * dt;
        for r in 0..3 { for c in 0..3 { phi[(6 + r, 6 + c)] = att_att[(r, c)]; } }
        
        let att_gbias = -r_b_e.matrix() * dt;
        for r in 0..3 { for c in 0..3 { phi[(6 + r, 12 + c)] = att_gbias[(r, c)]; } }
    }
    
    if crate::filter::CORE_STATE_SIZE > 15 {
        phi[(15, 16)] = dt;
    }
    
    phi
}

pub fn compute_process_noise(dt: f64, config: &EngineConfig, is_imu_active: bool, is_fixed: bool, num_amb: usize) -> DMatrix<f64> {
    let n = crate::filter::CORE_STATE_SIZE + num_amb;
    let mut q = DMatrix::<f64>::zeros(n, n);
    let dt_abs = dt.abs();
    
    if !is_imu_active {
        let q_acc = match config.dynamics_model {
            DynamicsModel::Static => 0.001,
            DynamicsModel::Pedestrian => 1.0,
            DynamicsModel::Marine => 2.0,
            DynamicsModel::Automotive => 10.0,
            DynamicsModel::Airborne => 50.0,
        };
        let q_pos = q_acc * dt_abs.powi(3) / 3.0; 
        let q_vel = q_acc * dt_abs;
        let q_pos_vel = q_acc * dt_abs.powi(2) / 2.0;
        for i in 0..3 { 
            q[(i, i)] = q_pos; 
            q[(i+3, i+3)] = q_vel; 
            q[(i, i+3)] = q_pos_vel;
            q[(i+3, i)] = q_pos_vel;
        }
        for i in 6..9 { q[(i, i)] = 1e-7 * dt_abs; } 
    } else {
        let q_vel = config.tuning.sigma_v * config.tuning.sigma_v * dt_abs;
        let q_att = config.tuning.sigma_phi * config.tuning.sigma_phi * dt_abs;
        let q_ab = config.tuning.sigma_ab * config.tuning.sigma_ab * dt_abs;
        let q_gb = config.tuning.sigma_gb * config.tuning.sigma_gb * dt_abs;
        for i in 0..3 {
            q[(3+i, 3+i)] = q_vel;
            q[(6+i, 6+i)] = q_att;
            q[(9+i, 9+i)] = q_ab;
            q[(12+i, 12+i)] = q_gb;
        }
    }
    
    if crate::filter::CORE_STATE_SIZE > 15 {
        q[(15, 15)] = config.process_noise_cb * dt_abs;
        q[(16, 16)] = config.process_noise_cd * dt_abs;
        q[(17, 17)] = config.process_noise_zwd * dt_abs;
    }

    for i in crate::filter::CORE_STATE_SIZE..n {
        q[(i, i)] = if is_fixed { config.process_noise_amb_fixed } else { config.process_noise_amb_float * dt_abs };
    }
    
    q
}

pub fn predict(state: &mut RtkState, dt: f64, config: &EngineConfig, imu_buffer: &[gneiss_core::imu::ImuMeasurement]) {
    if imu_buffer.is_empty() {
        state.position.vector += state.velocity * dt;
    } else {
        integrate_imu_mechanization(state, dt, imu_buffer);
    }
    
    if crate::filter::CORE_STATE_SIZE > 15 {
        state.rcv_clk_bias += state.rcv_clk_drift * dt;
    }
    
    let phi = compute_transition_matrix(state, dt, imu_buffer);
    let q = compute_process_noise(dt, config, !imu_buffer.is_empty(), state.is_fixed, state.ambiguities.len());
    
    state.core_phi = Some(phi.view((0, 0), (crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE)).into_owned());
    
    let mut phi_full = DMatrix::identity(state.covariance.nrows(), state.covariance.ncols());
    phi_full.view_mut((0, 0), (crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE)).copy_from(&state.core_phi.as_ref().unwrap());
    
    state.covariance = &phi_full * &state.covariance * phi_full.transpose() + q;
    state.full_p_predict = Some(state.covariance.clone());
    
    let mut x_pred = DVector::zeros(state.covariance.nrows());
    x_pred.rows_mut(0, 3).copy_from(&state.position.vector);
    x_pred.rows_mut(3, 3).copy_from(&state.velocity);
    if state.covariance.nrows() > 6 {
        x_pred.rows_mut(9, 3).copy_from(&state.accel_bias);
        x_pred.rows_mut(12, 3).copy_from(&state.gyro_bias);
    }
    if crate::filter::CORE_STATE_SIZE > 15 {
        x_pred[15] = state.rcv_clk_bias;
        x_pred[16] = state.rcv_clk_drift;
        x_pred[17] = state.zwd;
    }
    for i in 0..state.ambiguities.len() {
        x_pred[crate::filter::CORE_STATE_SIZE + i] = state.ambiguities[i];
    }
    state.full_x_predict = Some(x_pred);
}

pub fn gravity_wgs84(pos_ecef: Vector3<f64>) -> Vector3<f64> {
    let x = pos_ecef.x;
    let y = pos_ecef.y;
    let z = pos_ecef.z;
    let r = pos_ecef.norm();
    if r < 1.0 { return Vector3::zeros(); }
    
    let r2 = r * r;
    let r3 = r2 * r;
    let a = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
    let mu = 3.986005e14;
    let j2 = 1.082627e-3;
    
    let a_r_2 = (a / r) * (a / r);
    let z_r_2 = (z / r) * (z / r);
    
    let g_base = -mu / r3;
    let g_j2_common = 1.5 * j2 * a_r_2;
    
    let gx = g_base * x * (1.0 - g_j2_common * (5.0 * z_r_2 - 1.0));
    let gy = g_base * y * (1.0 - g_j2_common * (5.0 * z_r_2 - 1.0));
    let gz = g_base * z * (1.0 - g_j2_common * (5.0 * z_r_2 - 3.0));
    
    Vector3::new(gx, gy, gz)
}

fn skew_symmetric(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(
        0.0, -v.z,  v.y,
        v.z,  0.0, -v.x,
       -v.y,  v.x,  0.0
    )
}
