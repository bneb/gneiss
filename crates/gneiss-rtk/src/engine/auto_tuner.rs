use crate::filter::RtkState;
use crate::engine::config::EkfTuningConfig;
use nalgebra::Vector3;

/// Evaluates a preliminary processing pass and returns a dynamically tuned EkfTuningConfig.
pub fn tune_ekf_parameters(
    state_history: &[RtkState],
    imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
    mut base_config: EkfTuningConfig,
) -> EkfTuningConfig {
    if !base_config.auto_tune.enabled {
        return base_config;
    }

    // 1. Extract IMU Bias Instabilities
    let (sigma_ab, sigma_gb) = extract_imu_statistics(state_history, imu_history);
    
    // 2. Extract GNSS Multipath/Noise Scaling
    // let pr_scale = extract_gnss_statistics(state_history);
    // We will just do a simple multipath estimation based on residual norms if we had them.
    // For now, we focus on IMU since it dictates the INS drift during outages.
    
    let bounds = &base_config.auto_tune;
    
    if let Some(ab) = sigma_ab {
        base_config.sigma_ab = ab.clamp(bounds.min_sigma_ab, bounds.max_sigma_ab);
        tracing::info!("AutoTuner: Configured sigma_ab to {:.3e}", base_config.sigma_ab);
    }
    
    if let Some(gb) = sigma_gb {
        base_config.sigma_gb = gb.clamp(bounds.min_sigma_gb, bounds.max_sigma_gb);
        tracing::info!("AutoTuner: Configured sigma_gb to {:.3e}", base_config.sigma_gb);
    }

    base_config
}

/// Identifies stationary periods in the dataset to estimate the actual accelerometer
/// and gyroscope bias instabilities (random walk).
pub fn extract_imu_statistics(
    state_history: &[RtkState],
    imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>]
) -> (Option<f64>, Option<f64>) {
    if state_history.len() < 100 || imu_history.len() < 100 {
        return (None, None);
    }

    let mut stationary_accel_variances = Vec::new();
    let mut stationary_gyro_variances = Vec::new();

    for (k, state) in state_history.iter().enumerate() {
        // Find epochs where the vehicle is clearly stationary
        // Velocity < 0.05 m/s, and SPP/RTK is locked (position var is relatively low)
        let is_stationary = state.velocity.norm() < 0.05 && state.covariance[(0,0)] < 5.0;

        if is_stationary && k < imu_history.len() {
            let imu_buf = &imu_history[k];
            if imu_buf.len() > 10 {
                let mut sum_a = Vector3::zeros();
                let mut sum_g = Vector3::zeros();
                for m in imu_buf {
                    sum_a += m.accel;
                    sum_g += m.gyro;
                }
                let mean_a = sum_a / (imu_buf.len() as f64);
                let mean_g = sum_g / (imu_buf.len() as f64);

                let mut var_a = 0.0;
                let mut var_g = 0.0;
                for m in imu_buf {
                    var_a += (m.accel - mean_a).norm_squared();
                    var_g += (m.gyro - mean_g).norm_squared();
                }
                var_a /= imu_buf.len() as f64;
                var_g /= imu_buf.len() as f64;

                stationary_accel_variances.push(var_a);
                stationary_gyro_variances.push(var_g);
            }
        }
    }

    if stationary_accel_variances.is_empty() {
        return (None, None);
    }

    // Use median variance during stationary periods to estimate bias instability
    stationary_accel_variances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    stationary_gyro_variances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let median_var_a = stationary_accel_variances[stationary_accel_variances.len() / 2];
    let median_var_g = stationary_gyro_variances[stationary_gyro_variances.len() / 2];

    // Convert observed variance to random walk sigma
    let sigma_ab = median_var_a.sqrt() * 0.01; // heuristic scaling
    let sigma_gb = median_var_g.sqrt() * 0.01;

    (Some(sigma_ab), Some(sigma_gb))
}

pub fn extract_gnss_statistics(_state_history: &[RtkState]) -> f64 {
    // Placeholder for extracting GNSS multipath variance scaling
    1.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use nalgebra::Vector3;
    use gneiss_core::imu::ImuMeasurement;

    #[test]
    fn test_auto_tuner_extracts_imu_bias() {
        let mut state_history = Vec::new();
        let mut imu_history = Vec::new();

        let t0 = GpsTime::new(2000, 100000.0);
        let pos = Coordinate::new(Vector3::new(0.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, t0);

        for _i in 0..150 {
            let mut state = RtkState::new(t0, pos, 1.0);
            state.velocity = Vector3::new(0.01, 0.01, 0.0); // Stationary
            state.covariance[(0,0)] = 1.0;
            state_history.push(state);

            let mut imu_buf = Vec::new();
            for _ in 0..20 {
                imu_buf.push(ImuMeasurement {
                    time_tag: 0,
                    accel: Vector3::new(0.0, 0.0, 9.81) + Vector3::new(0.05, -0.02, 0.01), // some variance
                    gyro: Vector3::new(0.0, 0.0, 0.0) + Vector3::new(0.001, 0.002, -0.001),
                    temperature: Some(20.0),
                });
                imu_buf.push(ImuMeasurement {
                    time_tag: 0,
                    accel: Vector3::new(0.0, 0.0, 9.81) + Vector3::new(-0.05, 0.02, -0.01),
                    gyro: Vector3::new(0.0, 0.0, 0.0) + Vector3::new(-0.001, -0.002, 0.001),
                    temperature: Some(20.0),
                });
            }
            imu_history.push(imu_buf);
        }

        let (ab, gb) = extract_imu_statistics(&state_history, &imu_history);
        assert!(ab.is_some());
        assert!(gb.is_some());

        let ab_val = ab.unwrap();
        let gb_val = gb.unwrap();
        
        assert!(ab_val > 1e-6);
        assert!(gb_val > 1e-6);
    }
}
