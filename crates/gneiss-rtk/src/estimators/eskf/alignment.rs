use nalgebra::{Rotation3, UnitQuaternion, Vector3};

use crate::estimators::eskf::types::EskfState;
use crate::swfg::imu_preintegration::ImuSample;

/// Estimate stationary gyroscope bias by averaging static samples.
pub fn compute_gyro_bias(imu_samples: &[ImuSample], max_samples: usize) -> Vector3<f64> {
    let n = imu_samples.len().clamp(1, max_samples);
    let mut sum_g = Vector3::zeros();
    for s in &imu_samples[..n] {
        sum_g += s.gyro;
    }
    sum_g / n as f64
}

/// Compute leveling roll and pitch angles from specific force under gravity.
pub fn compute_leveling_angles(imu_samples: &[ImuSample], max_samples: usize) -> (f64, f64) {
    let n = imu_samples.len().clamp(1, max_samples);
    let mut sum_a = Vector3::zeros();
    for s in &imu_samples[..n] {
        sum_a += s.accel;
    }
    let mean_a = sum_a / n as f64;
    let pitch = (mean_a.x / 9.7803).clamp(-0.5, 0.5);
    let roll = (-mean_a.y / 9.7803).clamp(-0.5, 0.5);
    (roll, pitch)
}

/// Compute initial body-to-ECEF attitude quaternion from leveling angles, position, and heading.
pub fn compute_initial_attitude(
    imu_samples: &[ImuSample],
    init_pos: Vector3<f64>,
    heading_rad: f64,
    max_samples: usize,
) -> UnitQuaternion<f64> {
    let (roll, pitch) = compute_leveling_angles(imu_samples, max_samples);
    let llh = gneiss_core::coords::ecef_to_llh(init_pos);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(llh).transpose();
    let r_body = Rotation3::from_euler_angles(roll, pitch, heading_rad);
    let rot = Rotation3::from_matrix_unchecked(ned_to_ecef * r_body.matrix());
    UnitQuaternion::from_rotation_matrix(&rot)
}

/// Initialize an ESKF state with specified position, velocity, attitude, and gyro bias.
pub fn init_eskf_filter(
    init_pos: Vector3<f64>,
    init_vel: Vector3<f64>,
    init_att: UnitQuaternion<f64>,
    gyro_bias: Vector3<f64>,
) -> EskfState {
    let mut state = EskfState::new(init_pos, init_vel, init_att);
    state.gyro_bias = gyro_bias;
    state
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_sample(accel: Vector3<f64>, gyro: Vector3<f64>) -> ImuSample {
        ImuSample {
            accel,
            gyro,
            time_us: 0,
        }
    }

    #[test]
    fn test_compute_gyro_bias_averaging() {
        let samples = vec![
            make_sample(Vector3::zeros(), Vector3::new(0.01, -0.02, 0.03)),
            make_sample(Vector3::zeros(), Vector3::new(0.03, -0.04, 0.05)),
        ];
        let bias = compute_gyro_bias(&samples, 10);
        assert!((bias.x - 0.02).abs() < 1e-12);
        assert!((bias.y - (-0.03)).abs() < 1e-12);
        assert!((bias.z - 0.04).abs() < 1e-12);
    }

    #[test]
    fn test_compute_leveling_angles_level() {
        // Specific force when upright and stationary: +9.7803 upwards in NED or z-axis
        let samples = vec![make_sample(Vector3::new(0.0, 0.0, 9.7803), Vector3::zeros())];
        let (roll, pitch) = compute_leveling_angles(&samples, 10);
        assert!(roll.abs() < 1e-6);
        assert!(pitch.abs() < 1e-6);
    }

    #[test]
    fn test_init_eskf_filter() {
        let p = Vector3::new(100.0, 200.0, 300.0);
        let v = Vector3::new(1.0, 2.0, 3.0);
        let q = UnitQuaternion::identity();
        let bg = Vector3::new(0.001, -0.002, 0.003);
        let state = init_eskf_filter(p, v, q, bg);
        assert_eq!(state.pos_ecef, p);
        assert_eq!(state.vel_ecef, v);
        assert_eq!(state.gyro_bias, bg);
    }
}
