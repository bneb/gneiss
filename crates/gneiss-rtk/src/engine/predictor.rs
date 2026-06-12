use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion, Matrix3};

use crate::engine::{DynamicsModel, EngineConfig};

pub fn predict(state: &mut RtkState, dt: f64, config: &EngineConfig, imu_buffer: &[gneiss_core::imu::ImuMeasurement]) {
    let n = state.covariance.nrows();
    let mut phi = DMatrix::<f64>::identity(n, n);
    let omega_ie = Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S); // Earth rotation rate
    
    if imu_buffer.is_empty() {
        // Standard GNSS-only kinematic model
        state.position.vector += state.velocity * dt;
        phi[(0, 3)] = dt; phi[(1, 4)] = dt; phi[(2, 5)] = dt;
    } else {
        // High-Rate INS Mechanization
        let imu_dt = dt / (imu_buffer.len() as f64);
        
        for meas in imu_buffer {
            // 1. Correct measurements with current bias estimates
            let f_b = meas.accel - state.accel_bias;
            let omega_b = meas.gyro - state.gyro_bias;
            
            // 2. Attitude Update (Strapdown)
            // zeta is the rotation vector over imu_dt
            let zeta = (omega_b - state.attitude.inverse() * omega_ie) * imu_dt;
            let angle = zeta.norm();
            let dq = if angle > 1e-12 {
                UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(zeta / angle), angle)
            } else {
                UnitQuaternion::identity()
            };

            // Midpoint attitude for velocity update (Sculling correction)
            let dq_mid = if angle > 1e-12 {
                UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(zeta / angle), angle * 0.5)
            } else {
                UnitQuaternion::identity()
            };
            let r_mid = state.attitude * dq_mid;
            
            // Final attitude update
            state.attitude *= dq;
            state.attitude.renormalize();

            // 3. Velocity Update
            let f_e = r_mid * f_b;
            let gravity = gravity_wgs84(state.position.vector);
            let coriolis = 2.0 * omega_ie.cross(&state.velocity);
            let centrifugal = omega_ie.cross(&(omega_ie.cross(&state.position.vector)));
            
            let v_dot = f_e + gravity - coriolis - centrifugal;
            
            // Trapezoidal velocity integration
            let v_mid = state.velocity + v_dot * (imu_dt * 0.5);
            state.velocity += v_dot * imu_dt;
            
            // 4. Position Update
            state.position.vector += v_mid * imu_dt;
        }

        // --- Error-State Propagation (Phi Matrix) ---
        // Derived from linearized INS error equations in ECEF frame
        
        let r_b_e = state.attitude.to_rotation_matrix();
        let f_e = state.attitude * (imu_buffer.last().unwrap().accel - state.accel_bias);
        let f_e_skew = skew_symmetric(&f_e);
        let omega_ie_skew = skew_symmetric(&omega_ie);
        
        // Sub-matrices of the 15x15 F matrix (F is continuous-time transition)
        // States: [Pos(0:3), Vel(3:6), Att(6:9), AccelBias(9:12), GyroBias(12:15)]
        
        // dPos/dVel = I
        for i in 0..3 { phi[(i, 3 + i)] = dt; }
        
        // dVel/dAtt = -[f_e x]
        let vel_att = -f_e_skew * dt;
        for r in 0..3 { for c in 0..3 { phi[(3 + r, 6 + c)] = vel_att[(r, c)]; } }
        
        // dVel/dVel = I - 2 [omega_ie x] dt
        let vel_vel = Matrix3::identity() - 2.0 * omega_ie_skew * dt;
        for r in 0..3 { for c in 0..3 { phi[(3 + r, 3 + c)] = vel_vel[(r, c)]; } }
        
        // dVel/dAccelBias = -R_b_e
        let vel_abias = -r_b_e.matrix() * dt;
        for r in 0..3 { for c in 0..3 { phi[(3 + r, 9 + c)] = vel_abias[(r, c)]; } }
        
        // dAtt/dAtt = I - [omega_ie x] dt
        let att_att = Matrix3::identity() - omega_ie_skew * dt;
        for r in 0..3 { for c in 0..3 { phi[(6 + r, 6 + c)] = att_att[(r, c)]; } }
        
        // dAtt/dGyroBias = -R_b_e
        let att_gbias = -r_b_e.matrix() * dt;
        for r in 0..3 { for c in 0..3 { phi[(6 + r, 12 + c)] = att_gbias[(r, c)]; } }
        
        // Biases are modeled as Random Walks (Phi_bias = I, already identity)
    }
    
    // Process Noise Q
    let mut q = DMatrix::<f64>::zeros(n, n);
    
    let dt_abs = dt.abs();
    if imu_buffer.is_empty() {
        // Standard GNSS-only dynamic model (predicts position using velocity)
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
        // IMU-specific process noise (based on VRW, ARW)
        let sigma_v = config.tuning.sigma_v; // Velocity Random Walk (m/s/sqrt(s))
        let sigma_phi = config.tuning.sigma_phi; // Angular Random Walk (rad/sqrt(s))
        let sigma_ab = config.tuning.sigma_ab; // Accel bias instability
        let sigma_gb = config.tuning.sigma_gb; // Gyro bias instability
        
        let q_vel = sigma_v * sigma_v * dt_abs;
        let q_att = sigma_phi * sigma_phi * dt_abs;
        let q_ab = sigma_ab * sigma_ab * dt_abs;
        let q_gb = sigma_gb * sigma_gb * dt_abs;
        
        for i in 0..3 {
            q[(3+i, 3+i)] = q_vel;
            q[(6+i, 6+i)] = q_att;
            q[(9+i, 9+i)] = q_ab;
            q[(12+i, 12+i)] = q_gb;
        }
    }
    
    if crate::filter::CORE_STATE_SIZE > 15 {
        state.rcv_clk_bias += state.rcv_clk_drift * dt;
        phi[(15, 16)] = dt;
        
        // Clock models: random walk + integrated random walk
        // Unsteered TCXO can drift significantly between epochs
        let q_cb = config.process_noise_cb * dt_abs;
        let q_cd = config.process_noise_cd * dt_abs;
        let q_zwd = config.process_noise_zwd * dt_abs;
        
        q[(15, 15)] = q_cb;
        q[(16, 16)] = q_cd;
        q[(17, 17)] = q_zwd;
    }

    // Ambiguity noise
    for i in crate::filter::CORE_STATE_SIZE..n {
        if state.is_fixed { 
            q[(i, i)] = config.process_noise_amb_fixed; 
        } else { 
            q[(i, i)] = config.process_noise_amb_float * dt_abs; 
        }
    }
    
    state.core_phi = Some(phi.view((0, 0), (crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE)).into_owned());
    
    let mut phi_full = DMatrix::identity(n, n);
    let core_size = crate::filter::CORE_STATE_SIZE;
    phi_full.view_mut((0, 0), (core_size, core_size)).copy_from(&phi.view((0, 0), (core_size, core_size)));
    
    let p_new = &phi_full * &state.covariance * phi_full.transpose();
    state.covariance = p_new + q;
    
    state.full_p_predict = Some(state.covariance.clone());
    
    // Save the predicted state vector for the RTS smoother
    let mut x_pred = DVector::zeros(n);
    x_pred.rows_mut(0, 3).copy_from(&state.position.vector);
    x_pred.rows_mut(3, 3).copy_from(&state.velocity);
    // Attitude error is 0 since the reference attitude was just updated
    if n > 6 {
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
