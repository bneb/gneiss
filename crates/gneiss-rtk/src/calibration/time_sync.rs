use nalgebra::Vector3;

/// Computes the time offset between GNSS and IMU sensors using cross-correlation of velocity profiles.
/// Returns the offset in seconds. A positive offset means GNSS time is lagging behind IMU time.
pub fn cross_correlate_time_offset(
    gnss_time: &[f64],
    gnss_vel: &[Vector3<f64>],
    imu_time: &[f64],
    imu_vel: &[Vector3<f64>],
) -> f64 {
    if gnss_time.is_empty() || imu_time.len() < 2 {
        return 0.0;
    }

    let mut best_tau = 0.0;
    let mut max_corr = f64::NEG_INFINITY;

    // Sweep from -0.5 to +0.5 seconds with 1 ms resolution
    let steps = 500;
    for step in -steps..=steps {
        let tau = (step as f64) * 0.001;
        let mut sum_xy = 0.0;
        let mut sum_xx = 0.0;
        let mut sum_yy = 0.0;
        let mut count = 0;

        for i in 0..gnss_time.len() {
            let t_eval = gnss_time[i] - tau;
            
            let idx = match imu_time.binary_search_by(|t| t.partial_cmp(&t_eval).unwrap_or(core::cmp::Ordering::Equal)) {
                Ok(j) => j,
                Err(j) => j,
            };

            if idx > 0 && idx < imu_time.len() {
                let t0 = imu_time[idx - 1];
                let t1 = imu_time[idx];
                let v0 = imu_vel[idx - 1];
                let v1 = imu_vel[idx];

                let dt = t1 - t0;
                if dt > 0.0 {
                    let alpha = (t_eval - t0) / dt;
                    let v_interp = v0 + (v1 - v0) * alpha;
                    
                    let x = gnss_vel[i];
                    let y = v_interp;
                    
                    sum_xy += x.dot(&y);
                    sum_xx += x.norm_squared();
                    sum_yy += y.norm_squared();
                    count += 1;
                }
            }
        }

        if count > 10 && sum_xx > 1e-6 && sum_yy > 1e-6 { 
            let corr = sum_xy / (sum_xx.sqrt() * sum_yy.sqrt());
            if corr > max_corr {
                max_corr = corr;
                best_tau = tau;
            }
        }
    }

    best_tau
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_cross_correlate_time_offset() {
        let mut gnss_time = Vec::new();
        let mut gnss_vel = Vec::new();
        let mut imu_time = Vec::new();
        let mut imu_vel = Vec::new();

        // Simulate a sine wave velocity profile
        let freq = 0.5; // 0.5 Hz
        let true_offset = 0.098; // 98 ms offset
        
        // IMU runs at 100 Hz
        for i in 0..1000 {
            let t = i as f64 * 0.01;
            imu_time.push(t);
            let v = (t * core::f64::consts::TAU * freq).sin();
            imu_vel.push(Vector3::new(v, v * 0.5, -v * 0.2));
        }

        // GNSS runs at 5 Hz, delayed by true_offset
        for i in 0..50 {
            let t = i as f64 * 0.2;
            gnss_time.push(t);
            // GNSS measures the velocity at the delayed time
            let v = ((t - true_offset) * core::f64::consts::TAU * freq).sin();
            gnss_vel.push(Vector3::new(v, v * 0.5, -v * 0.2));
        }

        let estimated_offset = cross_correlate_time_offset(&gnss_time, &gnss_vel, &imu_time, &imu_vel);
        
        assert!((estimated_offset - true_offset).abs() < 0.005, "Expected offset {}, got {}", true_offset, estimated_offset);
    }
}
