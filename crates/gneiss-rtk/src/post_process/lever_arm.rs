//! Automatic Antenna-to-IMU Lever-Arm Estimation.
//!
//! Solves for the 3D lever-arm vector $\mathbf{l}_b = [l_x, l_y, l_z]^T$ relating the GNSS
//! antenna phase center to the IMU center of navigation during vehicle dynamic maneuvers:
//!
//!   \mathbf{a}_\text{GNSS} - \mathbf{a}_\text{IMU} = \left( [\boldsymbol{\omega} \times]^2 + [\dot{\boldsymbol{\omega}} \times] \right) \mathbf{l}_b
//!
//! Implements least-squares batch estimation across maneuvering intervals.

use nalgebra::{Matrix3, Vector3};

/// Coupled GNSS and IMU kinematic observation at epoch $k$.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LeverArmObservation {
    /// Measured angular velocity in body frame (rad/s).
    pub omega_body: Vector3<f64>,
    /// Angular acceleration in body frame (rad/s^2).
    pub alpha_body: Vector3<f64>,
    /// Measured IMU specific force / linear acceleration in body frame (m/s^2).
    pub accel_imu_body: Vector3<f64>,
    /// Differentiated GNSS kinematic acceleration rotated into body frame (m/s^2).
    pub accel_gnss_body: Vector3<f64>,
}

/// Result of lever-arm estimation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LeverArmEstimate {
    /// Estimated lever arm vector [Forward, Right, Down] in meters.
    pub lever_arm_body: Vector3<f64>,
    /// Standard deviations in meters.
    pub std_body: Vector3<f64>,
    /// Number of valid dynamic epochs used in solve.
    pub epochs_used: usize,
    /// Residual root-mean-square error (m/s^2).
    pub residual_rms: f64,
}

/// Batch estimator for antenna-to-IMU lever arm.
pub struct LeverArmEstimator;

impl LeverArmEstimator {
    /// Solves for the 3D lever-arm vector from a collection of kinematic epochs.
    ///
    /// Requires dynamic rotation (turning/banking maneuvers) to resolve the cross-product terms.
    pub fn estimate(obs: &[LeverArmObservation]) -> Option<LeverArmEstimate> {
        if obs.len() < 10 {
            return None;
        }

        let mut ata = Matrix3::<f64>::zeros();
        let mut atb = Vector3::<f64>::zeros();
        let mut count = 0;

        for o in obs {
            // Require minimum angular rate or acceleration to avoid degenerate geometry
            let omega_norm = o.omega_body.norm();
            let alpha_norm = o.alpha_body.norm();
            if omega_norm < 0.05 && alpha_norm < 0.05 {
                continue;
            }

            // Cross-coupling matrix H = [omega x][omega x] + [alpha x]
            let w_skew = skew_symmetric(o.omega_body);
            let a_skew = skew_symmetric(o.alpha_body);
            let h = w_skew * w_skew + a_skew;

            let z = o.accel_gnss_body - o.accel_imu_body;

            ata += h.transpose() * h;
            atb += h.transpose() * z;
            count += 1;
        }

        if count < 10 {
            return None;
        }

        let cov = ata.try_inverse()?;
        let lever_arm = cov * atb;

        // Compute residuals
        let mut sum_res2 = 0.0;
        for o in obs {
            let w_skew = skew_symmetric(o.omega_body);
            let a_skew = skew_symmetric(o.alpha_body);
            let h = w_skew * w_skew + a_skew;
            let pred = h * lever_arm;
            let diff = (o.accel_gnss_body - o.accel_imu_body) - pred;
            sum_res2 += diff.norm_squared();
        }
        let rms = libm::sqrt(sum_res2 / (count as f64 * 3.0));

        let std_x = libm::sqrt(cov[(0, 0)].max(0.0) * rms * rms);
        let std_y = libm::sqrt(cov[(1, 1)].max(0.0) * rms * rms);
        let std_z = libm::sqrt(cov[(2, 2)].max(0.0) * rms * rms);

        Some(LeverArmEstimate {
            lever_arm_body: lever_arm,
            std_body: Vector3::new(std_x, std_y, std_z),
            epochs_used: count,
            residual_rms: rms,
        })
    }
}

fn skew_symmetric(v: Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(
        0.0, -v.z, v.y,
        v.z, 0.0, -v.x,
        -v.y, v.x, 0.0,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lever_arm_estimation_synthetic_rotation() {
        let true_arm = Vector3::new(0.5, -0.2, -1.2);
        let mut obs = Vec::new();

        // Simulate multi-axis turning and pitching maneuvers
        for i in 0..100 {
            let t = i as f64 * 0.1;
            let omega = Vector3::new(
                0.3 * libm::cos(t),
                0.4 * libm::sin(0.7 * t),
                0.5 * libm::sin(t),
            );
            let alpha = Vector3::new(
                -0.3 * libm::sin(t),
                0.28 * libm::cos(0.7 * t),
                0.5 * libm::cos(t),
            );
            let w_skew = skew_symmetric(omega);
            let a_skew = skew_symmetric(alpha);
            let h = w_skew * w_skew + a_skew;

            let a_imu = Vector3::new(1.0, 0.0, 9.81);
            let a_gnss = a_imu + h * true_arm;

            obs.push(LeverArmObservation {
                omega_body: omega,
                alpha_body: alpha,
                accel_imu_body: a_imu,
                accel_gnss_body: a_gnss,
            });
        }

        let est = LeverArmEstimator::estimate(&obs).expect("estimation converges");
        assert!((est.lever_arm_body[0] - true_arm[0]).abs() < 1e-4);
        assert!((est.lever_arm_body[1] - true_arm[1]).abs() < 1e-4);
        assert!((est.lever_arm_body[2] - true_arm[2]).abs() < 1e-4);
        assert!(est.epochs_used > 50);
    }
}
