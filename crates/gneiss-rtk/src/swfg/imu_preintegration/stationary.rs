//! IMU specific force variance and gyro rate detector for Zero Velocity Updates (ZUPT).
//!
//! Detects stationary periods (traffic stops, red lights, parking) by evaluating:
//! 1. Sample variance of specific force magnitude: s_f^2 < 0.05 (m/s^2)^2
//! 2. Mean angular velocity norm: ||omega|| < 0.05 rad/s
//! 3. Mean specific force norm consistency with 1g: |||f|| - 9.80665| < 1.0 m/s^2

use super::ImuSample;

/// Maximum specific force magnitude variance for stationary condition ((m/s^2)^2).
pub const MAX_ACCEL_VAR: f64 = 0.05;

/// Maximum mean angular velocity magnitude for stationary condition (rad/s).
pub const MAX_GYRO_NORM: f64 = 0.05;

/// Tolerance around nominal gravity magnitude (m/s^2).
pub const GRAVITY_TOLERANCE: f64 = 1.0;

/// Standard acceleration due to gravity (m/s^2).
pub const NOMINAL_GRAVITY: f64 = 9.80665;

/// Summary metrics computed over an IMU sample window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImuStationaryMetrics {
    pub accel_mean: f64,
    pub accel_var: f64,
    pub gyro_mean: f64,
}

impl ImuStationaryMetrics {
    /// Check whether these metrics meet the stationary criteria.
    pub fn is_stationary(&self) -> bool {
        self.accel_var < MAX_ACCEL_VAR
            && self.gyro_mean < MAX_GYRO_NORM
            && (self.accel_mean - NOMINAL_GRAVITY).abs() < GRAVITY_TOLERANCE
    }
}

/// Compute specific force and angular velocity metrics over a slice of IMU samples.
pub fn compute_stationary_metrics(samples: &[ImuSample]) -> Option<ImuStationaryMetrics> {
    let n = samples.len();
    if n < 2 {
        return None;
    }

    let n_f = n as f64;
    let accel_mean = samples.iter().map(|s| s.accel.norm()).sum::<f64>() / n_f;
    let accel_var = samples
        .iter()
        .map(|s| {
            let diff = s.accel.norm() - accel_mean;
            diff * diff
        })
        .sum::<f64>()
        / (n_f - 1.0);
    let gyro_mean = samples.iter().map(|s| s.gyro.norm()).sum::<f64>() / n_f;

    Some(ImuStationaryMetrics {
        accel_mean,
        accel_var,
        gyro_mean,
    })
}

/// Returns true if the IMU sample window indicates the vehicle is stationary.
pub fn detect_stationary(samples: &[ImuSample]) -> bool {
    compute_stationary_metrics(samples).is_some_and(|m| m.is_stationary())
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn mock_sample(ax: f64, ay: f64, az: f64, gx: f64, gy: f64, gz: f64, t_us: u64) -> ImuSample {
        ImuSample {
            accel: Vector3::new(ax, ay, az),
            gyro: Vector3::new(gx, gy, gz),
            time_us: t_us,
        }
    }

    #[test]
    fn test_stationary_detected_on_clean_static_imu() {
        let samples: Vec<ImuSample> = (0..20)
            .map(|i| {
                let noise = ((i % 3) as f64 - 1.0) * 0.005;
                mock_sample(0.01 + noise, -0.01 + noise, 9.805 + noise, 0.001, 0.002, -0.001, i * 50_000)
            })
            .collect();
        assert!(detect_stationary(&samples));
    }

    #[test]
    fn test_rejects_moving_due_to_accel_variance() {
        let samples: Vec<ImuSample> = (0..20)
            .map(|i| {
                let a = (i as f64) * 0.5;
                mock_sample(a, 0.0, 9.8, 0.001, 0.001, 0.001, i * 50_000)
            })
            .collect();
        assert!(!detect_stationary(&samples));
    }

    #[test]
    fn test_rejects_turning_due_to_high_gyro() {
        let samples: Vec<ImuSample> = (0..20)
            .map(|i| mock_sample(0.0, 0.0, 9.81, 0.0, 0.0, 0.25, i * 50_000))
            .collect();
        assert!(!detect_stationary(&samples));
    }

    #[test]
    fn test_rejects_insufficient_samples() {
        let sample = mock_sample(0.0, 0.0, 9.81, 0.0, 0.0, 0.0, 0);
        assert!(!detect_stationary(&[sample]));
        assert!(!detect_stationary(&[]));
    }
}
