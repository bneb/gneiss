//! Photogrammetric Camera & LiDAR to IMU Boresight and Lever-Arm Calibration.
//!
//! Solves for the 3D lever-arm vector $\mathbf{l}_c = [l_x, l_y, l_z]^T$ and
//! boresight rotation misalignment $\Delta \mathbf{R} = R_z(\Delta \psi) R_y(\Delta \theta) R_x(\Delta \phi)$
//! relating the optical/LiDAR sensor frame to the navigation IMU body frame:
//!
//!   \mathbf{p}_\text{GCP}^\text{ECEF} = \mathbf{p}_\text{IMU}^\text{ECEF} + \mathbf{R}_b^e \left( \mathbf{l}_c + \Delta \mathbf{R} \cdot \mathbf{p}_\text{target}^\text{sensor} \right)

use nalgebra::{DMatrix, DVector, Matrix3, UnitQuaternion, Vector3};

/// A sensor-to-ground target observation pairing during flight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoresightObservation {
    /// Navigation IMU position in ECEF at shutter / scan epoch (meters).
    pub pos_imu_ecef: Vector3<f64>,
    /// Navigation IMU attitude (body-to-ECEF rotation).
    pub q_b2e: UnitQuaternion<f64>,
    /// Target coordinates measured in the sensor coordinate frame (meters).
    pub pos_sensor: Vector3<f64>,
    /// Known Ground Control Point (GCP) or surveyed tie-point in ECEF (meters).
    pub pos_gcp_ecef: Vector3<f64>,
}

/// Result of boresight and lever-arm calibration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoresightCalibrationResult {
    /// Estimated sensor lever arm [Forward, Right, Down] in body frame (meters).
    pub lever_arm_body: Vector3<f64>,
    /// Estimated boresight misalignment Euler angles [roll, pitch, yaw] (radians).
    pub boresight_rpy_rad: Vector3<f64>,
    /// Residual 3D RMS error (meters).
    pub residual_rms_m: f64,
    /// Number of valid observations used in solve.
    pub observations_used: usize,
}

/// Batch least-squares solver for photogrammetric boresight & lever arm.
pub struct BoresightEstimator;

impl BoresightEstimator {
    /// Solve for boresight misalignment angles and lever-arm vector.
    pub fn estimate(
        obs: &[BoresightObservation],
        initial_lever_arm: Vector3<f64>,
    ) -> Option<BoresightCalibrationResult> {
        if obs.len() < 6 {
            return None;
        }

        let mut lever_arm = initial_lever_arm;
        let mut q_sensor = UnitQuaternion::identity();

        for _iter in 0..20 {
            let mut jtj = DMatrix::<f64>::zeros(6, 6);
            let mut jtr = DVector::<f64>::zeros(6);

            for o in obs {
                let r_b2e = o.q_b2e.to_rotation_matrix().into_inner();
                let p_cam_body = lever_arm + q_sensor * o.pos_sensor;
                let pred_gcp = o.pos_imu_ecef + r_b2e * p_cam_body;
                let residual = o.pos_gcp_ecef - pred_gcp; // 3x1

                // Jacobian w.r.t lever arm: J_l = R_b^e (3x3)
                // Jacobian w.r.t small rotation d_theta:
                // R(d_theta) * p_rot = p_rot - [p_rot x] * d_theta
                // so d(pred)/d(d_theta) = -R_b^e * [p_rot x]
                let p_rot = q_sensor * o.pos_sensor;
                let p_skew = Matrix3::new(
                    0.0, -p_rot.z, p_rot.y,
                    p_rot.z, 0.0, -p_rot.x,
                    -p_rot.y, p_rot.x, 0.0,
                );
                let j_ang = -r_b2e * p_skew;

                let mut j_block = DMatrix::<f64>::zeros(3, 6);
                for r in 0..3 {
                    for c in 0..3 {
                        j_block[(r, c)] = r_b2e[(r, c)];
                        j_block[(r, c + 3)] = j_ang[(r, c)];
                    }
                }

                let r_vec = DVector::from_column_slice(residual.as_slice());
                jtj += &j_block.transpose() * &j_block;
                jtr += &j_block.transpose() * &r_vec;
            }

            // Damped normal equations (Levenberg-Marquardt)
            for i in 0..6 {
                jtj[(i, i)] += 1e-6;
            }

            let delta = jtj.try_inverse()? * jtr;

            lever_arm += Vector3::new(delta[0], delta[1], delta[2]);
            let d_theta = Vector3::new(delta[3], delta[4], delta[5]);
            let dq = UnitQuaternion::from_scaled_axis(d_theta);
            q_sensor = dq * q_sensor;

            if delta.norm() < 1e-9 {
                break;
            }
        }

        // Compute final RMS
        let mut sum_sq = 0.0;
        for o in obs {
            let r_b2e = o.q_b2e.to_rotation_matrix().into_inner();
            let p_cam_body = lever_arm + q_sensor * o.pos_sensor;
            let pred_gcp = o.pos_imu_ecef + r_b2e * p_cam_body;
            sum_sq += (o.pos_gcp_ecef - pred_gcp).norm_squared();
        }
        let rms = libm::sqrt(sum_sq / (obs.len() as f64));

        let (roll, pitch, yaw) = q_sensor.euler_angles();

        Some(BoresightCalibrationResult {
            lever_arm_body: lever_arm,
            boresight_rpy_rad: Vector3::new(roll, pitch, yaw),
            residual_rms_m: rms,
            observations_used: obs.len(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boresight_and_lever_arm_estimation_exact_recovery() {
        let true_lever_arm = Vector3::new(0.15, -0.05, 0.25);
        let true_boresight = Vector3::new(0.01, -0.02, 0.03); // rad
        let delta_r = UnitQuaternion::from_euler_angles(true_boresight.x, true_boresight.y, true_boresight.z);

        let mut obs = Vec::new();
        for i in 0..24 {
            let yaw = (i as f64) * 0.4;
            let roll = 0.08 * ((i as f64) * 0.7).sin();
            let pitch = 0.05 * ((i as f64) * 0.5).cos();
            let q_b2e = UnitQuaternion::from_euler_angles(roll, pitch, yaw);
            let r_b2e = q_b2e.to_rotation_matrix().into_inner();

            let pos_imu = Vector3::new(100.0 + 50.0 * (i as f64).cos(), 200.0 + 50.0 * (i as f64).sin(), 150.0);
            let p_sensor = Vector3::new(
                15.0 * ((i as f64) * 1.3).sin(),
                10.0 * ((i as f64) * 1.7).cos(),
                50.0 + 5.0 * ((i as f64) * 0.9).sin(),
            );
            let p_cam_body = true_lever_arm + delta_r * p_sensor;
            let pos_gcp = pos_imu + r_b2e * p_cam_body;

            obs.push(BoresightObservation {
                pos_imu_ecef: pos_imu,
                q_b2e,
                pos_sensor: p_sensor,
                pos_gcp_ecef: pos_gcp,
            });
        }

        let res = BoresightEstimator::estimate(&obs, Vector3::zeros()).expect("estimate");
        assert!((res.lever_arm_body.x - true_lever_arm.x).abs() < 1e-4);
        assert!((res.lever_arm_body.y - true_lever_arm.y).abs() < 1e-4);
        assert!((res.lever_arm_body.z - true_lever_arm.z).abs() < 1e-4);

        assert!((res.boresight_rpy_rad.x - true_boresight.x).abs() < 1e-4);
        assert!((res.boresight_rpy_rad.y - true_boresight.y).abs() < 1e-4);
        assert!((res.boresight_rpy_rad.z - true_boresight.z).abs() < 1e-4);
        assert!(res.residual_rms_m < 1e-6);
    }
}
