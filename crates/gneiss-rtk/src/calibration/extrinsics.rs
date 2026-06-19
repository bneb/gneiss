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
        let w_skew = Matrix3::new(0.0, -w.z, w.y, w.z, 0.0, -w.x, -w.y, w.x, 0.0);

        let h_i = r_b_e[i] * w_skew;
        let z_i = v_gnss_e[i] - v_imu_e[i];

        h_matrix.fixed_view_mut::<3, 3>(i * 3, 0).copy_from(&h_i);
        z_vector.fixed_rows_mut::<3>(i * 3).copy_from(&z_i);
    }

    // Solve using normal equations: x = (H^T * H)^-1 * H^T * Z
    let h_t = h_matrix.transpose();
    let h_t_h = &h_t * &h_matrix;

    // Check if the matrix is invertible (requires dynamic excitation/turning)
    let h_t_h_inv = h_t_h
        .try_inverse()
        .ok_or("Matrix is singular; insufficient dynamic excitation to observe lever arm")?;

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

    let best_gnss = results
        .into_iter()
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();
    tracing::info!(
        "Best GNSS Lever Arm: [{:.2}, 0.0, {:.2}] with error metric {:.4}",
        best_gnss.0,
        best_gnss.1,
        best_gnss.2
    );

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

    let best_nhc = nhc_results
        .into_iter()
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();
    tracing::info!(
        "Best NHC Lever Arm: [{:.2}, 0.0, {:.2}] with error metric {:.4}",
        best_nhc.0,
        best_nhc.1,
        best_nhc.2
    );

    let best_nhc_arm = [best_nhc.0, 0.0, best_nhc.1];

    Ok((best_gnss_arm, best_nhc_arm))
}

/// Optimizes the 6-DOF extrinsics (lever arm X, Y, Z and mounting angles Roll, Pitch, Yaw)
/// using the Nelder-Mead method to minimize the provided evaluation function.
pub fn calibrate_extrinsics_6dof<F>(
    base_config: &crate::engine::EngineConfig,
    mut evaluate_fn: F,
) -> Result<([f64; 3], [f64; 3]), &'static str>
where
    F: FnMut(&crate::engine::EngineConfig) -> f64,
{
    tracing::info!("Starting 6-DOF Optimization for Extrinsics (Nelder-Mead)...");

    // Starting point: if base_config has values, use them, otherwise 0
    let start_arm = base_config.imu_to_antenna_lever_arm;
    let start_angles = base_config.imu_mounting_angles.unwrap_or([0.0, 0.0, 0.0]);

    // Initial step sizes: e.g., 0.5m for lever arm, 0.1 rad (~5.7 deg) for angles

    let mut objective = |p: &[f64; 6]| -> f64 {
        let mut cfg = base_config.clone();
        cfg.imu_to_antenna_lever_arm = [p[0], p[1], p[2]];
        cfg.imu_mounting_angles = Some([p[3], p[4], p[5]]);
        evaluate_fn(&cfg)
    };

    // Coarse 1D Grid Search for the X lever arm to find the global minimum basin
    tracing::info!("Running coarse 1D grid search for X lever arm to seed optimization...");
    let mut best_x = start_arm[0];
    let mut best_x_score = f64::MAX;
    for x_candidate in [0.0, 1.0, 2.0, 3.0] {
        let p_candidate = [
            x_candidate,
            start_arm[1],
            start_arm[2],
            start_angles[0],
            start_angles[1],
            start_angles[2],
        ];
        let score = objective(&p_candidate);
        if score < best_x_score {
            best_x_score = score;
            best_x = x_candidate;
        }
    }
    tracing::info!("Coarse search selected X lever arm: {:.1}", best_x);

    let start = [
        best_x,
        start_arm[1],
        start_arm[2],
        start_angles[0],
        start_angles[1],
        start_angles[2],
    ];

    let max_iter = 500;
    let tol = 1e-4;

    let mut simplex = vec![start; 7];
    let step = [0.5, 0.5, 0.5, 0.1, 0.1, 0.1]; // appropriate step for each dim
    for i in 0..6 {
        simplex[i + 1][i] += step[i];
    }

    let mut evals: Vec<([f64; 6], f64)> = simplex
        .into_iter()
        .map(|p| {
            let val = objective(&p);
            (p, val)
        })
        .collect();

    let alpha = 1.0;
    let gamma = 2.0;
    let rho = 0.5;
    let sigma = 0.5;

    for _iter in 0..max_iter {
        evals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let diff = evals.last().unwrap().1 - evals.first().unwrap().1;
        if diff < tol {
            break;
        }

        let mut centroid = [0.0; 6];
        for i in 0..6 {
            for j in 0..6 {
                centroid[j] += evals[i].0[j];
            }
        }
        for j in 0..6 {
            centroid[j] /= 6.0;
        }

        let worst = evals[6].0;
        let worst_score = evals[6].1;

        let mut reflected = [0.0; 6];
        for j in 0..6 {
            reflected[j] = centroid[j] + alpha * (centroid[j] - worst[j]);
        }
        let reflected_score = objective(&reflected);

        if reflected_score >= evals[0].1 && reflected_score < evals[5].1 {
            evals[6] = (reflected, reflected_score);
            continue;
        }

        if reflected_score < evals[0].1 {
            let mut expanded = [0.0; 6];
            for j in 0..6 {
                expanded[j] = centroid[j] + gamma * (reflected[j] - centroid[j]);
            }
            let expanded_score = objective(&expanded);
            if expanded_score < reflected_score {
                evals[6] = (expanded, expanded_score);
            } else {
                evals[6] = (reflected, reflected_score);
            }
            continue;
        }

        let mut contracted = [0.0; 6];
        let contract_outside = reflected_score < worst_score;
        if contract_outside {
            for j in 0..6 {
                contracted[j] = centroid[j] + rho * (reflected[j] - centroid[j]);
            }
        } else {
            for j in 0..6 {
                contracted[j] = centroid[j] + rho * (worst[j] - centroid[j]);
            }
        }

        let contracted_score = objective(&contracted);

        if contract_outside {
            if contracted_score <= reflected_score {
                evals[6] = (contracted, contracted_score);
                continue;
            }
        } else {
            if contracted_score < worst_score {
                evals[6] = (contracted, contracted_score);
                continue;
            }
        }

        let best = evals[0].0;
        for i in 1..7 {
            for j in 0..6 {
                evals[i].0[j] = best[j] + sigma * (evals[i].0[j] - best[j]);
            }
            let p = evals[i].0;
            evals[i].1 = objective(&p);
        }
    }

    evals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let best_params = evals[0].0;

    let best_lever_arm = [best_params[0], best_params[1], best_params[2]];
    let best_angles = [best_params[3], best_params[4], best_params[5]];

    tracing::info!(
        "Best 6-DOF Extrinsics: Lever Arm [{:.3}, {:.3}, {:.3}], Angles [{:.3}, {:.3}, {:.3}] with error {:.4}",
        best_lever_arm[0], best_lever_arm[1], best_lever_arm[2],
        best_angles[0], best_angles[1], best_angles[2], evals[0].1
    );

    Ok((best_lever_arm, best_angles))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Rotation3, Vector3};

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

        let estimated_arm =
            estimate_lever_arm(&omega_ib_b, &omega_dot_b, &v_gnss_e, &v_imu_e, &r_b_e).unwrap();

        assert!(
            (estimated_arm - true_lever_arm).norm() < 1e-3,
            "Expected {:?}, got {:?}",
            true_lever_arm,
            estimated_arm
        );
    }

    #[test]
    fn test_calibrate_extrinsics_6dof() {
        let base_config = crate::engine::EngineConfig::default();

        let true_lever_arm = [1.2, -0.3, 0.8];
        let true_angles = [0.05, -0.02, 0.1];

        let evaluate_fn = |cfg: &crate::engine::EngineConfig| -> f64 {
            let dx = cfg.imu_to_antenna_lever_arm[0] - true_lever_arm[0];
            let dy = cfg.imu_to_antenna_lever_arm[1] - true_lever_arm[1];
            let dz = cfg.imu_to_antenna_lever_arm[2] - true_lever_arm[2];

            let angles = cfg.imu_mounting_angles.unwrap_or([0.0, 0.0, 0.0]);
            let droll = angles[0] - true_angles[0];
            let dpitch = angles[1] - true_angles[1];
            let dyaw = angles[2] - true_angles[2];

            dx * dx + dy * dy + dz * dz + droll * droll + dpitch * dpitch + dyaw * dyaw
        };

        let result = calibrate_extrinsics_6dof(&base_config, evaluate_fn).unwrap();

        let lever_arm = result.0;
        let angles = result.1;

        assert!((lever_arm[0] - true_lever_arm[0]).abs() < 1e-2);
        assert!((lever_arm[1] - true_lever_arm[1]).abs() < 1e-2);
        assert!((lever_arm[2] - true_lever_arm[2]).abs() < 1e-2);

        assert!((angles[0] - true_angles[0]).abs() < 1e-2);
        assert!((angles[1] - true_angles[1]).abs() < 1e-2);
        assert!((angles[2] - true_angles[2]).abs() < 1e-2);
    }
}
