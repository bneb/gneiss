use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

/// Estimates the Lever Arm (IMU to Antenna offset in the Body frame) 
/// using a Least Squares optimization over dynamic maneuvers.
///
/// `omega_ib_b`: Angular rate of the body frame relative to inertial frame, expressed in body frame (Gyros)
/// `omega_dot_b`: Angular acceleration of the body frame (Derivative of Gyros)
/// `v_gnss_e`: Velocity measured by GNSS in ECEF frame
/// `v_imu_e`: Velocity integrated by IMU in ECEF frame (at the IMU center of navigation)
/// `r_b_e`: Rotation matrix from Body to ECEF frame
pub fn estimate_lever_arm(
    omega_ib_b: &[Vector3<f64>],
    _omega_dot_b: &[Vector3<f64>],
    v_gnss_e: &[Vector3<f64>],
    v_imu_e: &[Vector3<f64>],
    r_b_e: &[Matrix3<f64>],
) -> Result<Vector3<f64>, &'static str> {
    let n = v_gnss_e.len();
    if n < 3 || omega_ib_b.len() != n || v_imu_e.len() != n || r_b_e.len() != n {
        return Err("Insufficient or mismatched data for lever arm estimation");
    }

    // We are solving: Z = H * x
    // where x is the 3x1 lever arm vector.
    // Z = V_gnss - V_imu
    // V_gnss = V_imu + R_b_e * (omega_ib_b x lever_arm)
    // omega x lever_arm = [omega x] * lever_arm
    // So H = R_b_e * [omega x]
    
    let mut h_matrix = DMatrix::<f64>::zeros(n * 3, 3);
    let mut z_vector = DVector::<f64>::zeros(n * 3);

    for i in 0..n {
        let w = omega_ib_b[i];
        let w_skew = Matrix3::new(
            0.0, -w.z,  w.y,
             w.z,  0.0, -w.x,
            -w.y,  w.x,  0.0
        );
        
        let h_i = r_b_e[i] * w_skew;
        let z_i = v_gnss_e[i] - v_imu_e[i];

        h_matrix.fixed_view_mut::<3, 3>(i * 3, 0).copy_from(&h_i);
        z_vector.fixed_rows_mut::<3>(i * 3).copy_from(&z_i);
    }

    // Solve using normal equations: x = (H^T * H)^-1 * H^T * Z
    let h_t = h_matrix.transpose();
    let h_t_h = &h_t * &h_matrix;
    
    // Check if the matrix is invertible (requires dynamic excitation/turning)
    let h_t_h_inv = h_t_h.try_inverse().ok_or("Matrix is singular; insufficient dynamic excitation to observe lever arm")?;
    
    let lever_arm = h_t_h_inv * &h_t * z_vector;

    Ok(Vector3::new(lever_arm[0], lever_arm[1], lever_arm[2]))
}

/// Evaluates a specific GNSS lever arm configuration using a simplified processing pass
/// over a subset of the data (the first 1000 epochs).
pub fn calibrate_lever_arms_grid_search<F>(
    base_config: &crate::engine::EngineConfig,
    mut evaluate_fn: F,
) -> Result<([f64; 3], [f64; 3]), &'static str> 
where
    F: FnMut(&crate::engine::EngineConfig) -> f64,
{
    tracing::info!("Starting Grid Search for GNSS Lever Arm...");
    let x_range = [0.0, 0.5, 1.0, 1.5, 2.0];
    let z_range = [-1.0, -0.5, 0.0, 0.5, 1.0];
    
    let mut combinations = Vec::new();
    for &x in &x_range {
        for &z in &z_range {
            combinations.push((x, z));
        }
    }
    
    let mut results = Vec::new();
    for (x, z) in combinations {
        let mut cfg = base_config.clone();
        cfg.imu_to_antenna_lever_arm = [x, 0.0, z];
        cfg.enable_nhc = false; // Disable NHC while tuning GNSS lever arm
        
        let error = evaluate_fn(&cfg);
        results.push((x, z, error));
    }
    
    let best_gnss = results.into_iter().min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal)).unwrap();
    tracing::info!("Best GNSS Lever Arm: [{:.2}, 0.0, {:.2}] with error metric {:.4}", best_gnss.0, best_gnss.1, best_gnss.2);
    
    let best_gnss_arm = [best_gnss.0, 0.0, best_gnss.1];
    
    tracing::info!("Starting Grid Search for NHC Lever Arm...");
    let mut nhc_combinations = Vec::new();
    for &x in &[0.0, 0.5, 1.0, 1.5, 2.0] {
        for &z in &[0.0, 0.5, 1.0, 1.5, 2.0] {
            nhc_combinations.push((x, z));
        }
    }
    
    let mut nhc_results = Vec::new();
    for (x, z) in nhc_combinations {
        let mut cfg = base_config.clone();
        cfg.imu_to_antenna_lever_arm = best_gnss_arm;
        cfg.enable_nhc = true;
        cfg.imu_to_nhc_lever_arm = [x, 0.0, z];
        
        let error = evaluate_fn(&cfg);
        nhc_results.push((x, z, error));
    }
    
    let best_nhc = nhc_results.into_iter().min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal)).unwrap();
    tracing::info!("Best NHC Lever Arm: [{:.2}, 0.0, {:.2}] with error metric {:.4}", best_nhc.0, best_nhc.1, best_nhc.2);
    
    let best_nhc_arm = [best_nhc.0, 0.0, best_nhc.1];
    
    Ok((best_gnss_arm, best_nhc_arm))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Vector3, Rotation3};

    #[test]
    fn test_estimate_lever_arm() {
        let true_lever_arm = Vector3::new(1.5, -0.5, 2.0);
        
        let mut omega_ib_b = Vec::new();
        let mut omega_dot_b = Vec::new();
        let mut v_gnss_e = Vec::new();
        let mut v_imu_e = Vec::new();
        let mut r_b_e = Vec::new();

        // Simulate 100 epochs of dynamic turning across multiple axes
        for i in 0..100 {
            let t = (i as f64) * 0.1;
            
            // Spinning around multiple axes to excite all lever arm states
            let w = Vector3::new(t.sin(), t.cos(), 1.0);
            let w_dot = Vector3::new(t.cos(), -t.sin(), 0.0);
            
            let r = *Rotation3::from_euler_angles(t.sin(), t.cos(), t).matrix();
            
            // Base velocity at IMU
            let v_i = Vector3::new(10.0, 0.0, 0.0);
            
            // The GNSS velocity is V_imu + R_b_e * (omega x lever_arm)
            let v_g = v_i + r * w.cross(&true_lever_arm);

            omega_ib_b.push(w);
            omega_dot_b.push(w_dot);
            v_imu_e.push(v_i);
            v_gnss_e.push(v_g);
            r_b_e.push(r);
        }

        let estimated_arm = estimate_lever_arm(&omega_ib_b, &omega_dot_b, &v_gnss_e, &v_imu_e, &r_b_e).unwrap();
        
        assert!((estimated_arm - true_lever_arm).norm() < 1e-3, "Expected {:?}, got {:?}", true_lever_arm, estimated_arm);
    }
}
