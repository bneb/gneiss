use crate::filter::RtkState;
use nalgebra::{DMatrix, Vector3, UnitQuaternion, Matrix3};

pub fn predict(state: &mut RtkState, dt: f64, _q_var: f64, imu_buffer: &[gneiss_core::imu::ImuMeasurement]) {
    let n = state.covariance.nrows();
    let mut phi = DMatrix::<f64>::identity(n, n);
    let omega_ie = Vector3::new(0.0, 0.0, 7.2921151467e-5); // Earth rotation rate
    
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
        let q_acc = _q_var;
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
        let sigma_v = 0.01; // Velocity Random Walk (m/s/sqrt(s))
        let sigma_phi = 0.001; // Angular Random Walk (rad/sqrt(s))
        let sigma_ab = 1e-4; // Accel bias instability
        let sigma_gb = 1e-5; // Gyro bias instability
        
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
        let q_cb = 1e6 * dt_abs;
        let q_cd = 1e4 * dt_abs;
        let q_zwd = 1e-8 * dt_abs;
        
        q[(15, 15)] = q_cb;
        q[(16, 16)] = q_cd;
        q[(17, 17)] = q_zwd;
    }

    // Ambiguity noise
    for i in crate::filter::CORE_STATE_SIZE..n {
        if state.is_fixed { q[(i, i)] = 1e-12; } else { q[(i, i)] = 1e-8 * dt_abs; }
    }
    
    state.core_phi = Some(phi.view((0, 0), (crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE)).into_owned());
    
    let core_size = crate::filter::CORE_STATE_SIZE;
    let phi_core = phi.view((0, 0), (core_size, core_size));
    
    let p_core = state.covariance.view((0, 0), (core_size, core_size));
    let p_core_new = phi_core * p_core * phi_core.transpose();
    
    if n > core_size {
        let p_cross = state.covariance.view((0, core_size), (core_size, n - core_size));
        let p_cross_new = phi_core * p_cross;
        
        state.covariance.view_mut((0, 0), (core_size, core_size)).copy_from(&p_core_new);
        state.covariance.view_mut((0, core_size), (core_size, n - core_size)).copy_from(&p_cross_new);
        state.covariance.view_mut((core_size, 0), (n - core_size, core_size)).copy_from(&p_cross_new.transpose());
    } else {
        state.covariance.view_mut((0, 0), (core_size, core_size)).copy_from(&p_core_new);
    }
    state.covariance += q;
    
    state.core_p_predict = Some(state.covariance.view((0, 0), (crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE)).into_owned());
}

pub fn gravity_wgs84(pos_ecef: Vector3<f64>) -> Vector3<f64> {
    let x = pos_ecef.x;
    let y = pos_ecef.y;
    let z = pos_ecef.z;
    let r = pos_ecef.norm();
    if r < 1.0 { return Vector3::zeros(); }
    
    let r2 = r * r;
    let r3 = r2 * r;
    let a = 6378137.0;
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
