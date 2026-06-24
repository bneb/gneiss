use gneiss_core::imu::ImuMeasurement;
use nalgebra::{Rotation3, Vector3};

#[derive(Debug, Clone, Copy)]
pub struct MountingAngles {
    pub roll: f64,
    pub pitch: f64,
    pub yaw: f64,
}

impl MountingAngles {
    pub fn to_rotation(&self) -> Rotation3<f64> {
        Rotation3::from_euler_angles(self.roll, self.pitch, self.yaw)
    }

    pub fn apply(&self, vector: Vector3<f64>) -> Vector3<f64> {
        self.to_rotation() * vector
    }
}

/// Solves for the mounting Roll and Pitch by analyzing the gravity vector
/// during a period of relative stillness.
pub fn estimate_gravity_alignment(
    measurements: &[ImuMeasurement],
) -> Result<(f64, f64), &'static str> {
    if measurements.is_empty() {
        return Err("No IMU measurements for gravity alignment");
    }

    // 1. Calculate average acceleration vector
    let mut avg_accel = Vector3::zeros();
    for m in measurements {
        avg_accel += m.accel;
    }
    avg_accel /= measurements.len() as f64;

    // 2. We assume the vehicle is on level ground (or average level).
    // In FRD, Gravity is [0, 0, g].
    // Our sensor measures [ax, ay, az] = R_m_v * [0, 0, -g] (ignoring non-gravitational accel)
    // where R_m_v is Mounting-to-Vehicle rotation.

    let ax = avg_accel.x;
    let ay = avg_accel.y;
    let az = avg_accel.z;

    // Pitch: theta = atan2(ax, sqrt(ay^2 + az^2))
    let pitch = f64::atan2(ax, f64::sqrt(ay * ay + az * az));

    // Roll: phi = atan2(-ay, -az)
    let roll = f64::atan2(-ay, -az);

    Ok((roll, pitch))
}

/// Solves for the mounting Yaw (heading offset) by correlating the longitudinal
/// acceleration spikes with GNSS-derived velocity changes.
pub fn estimate_heading_alignment(
    imu_measurements: &[ImuMeasurement],
    gnss_velocities_ned: &[(f64, Vector3<f64>)], // (TOW, Vel_NED)
    mounting_roll: f64,
    mounting_pitch: f64,
) -> Result<f64, &'static str> {
    if imu_measurements.len() < 100 || gnss_velocities_ned.len() < 10 {
        return Err("Insufficient dynamic data for heading alignment");
    }

    // 1. Rotate all IMU measurements into the intermediate Level-Frame (Roll/Pitch corrected)
    // but with unknown Yaw.
    let _r_lev = Rotation3::from_euler_angles(mounting_roll, mounting_pitch, 0.0);

    let mut best_yaw = 0.0;
    let mut max_corr = -1.0;

    // 2. Search for the Yaw offset that maximizes correlation between
    // Horizontal IMU Accel and GNSS Acceleration.
    for y_deg in 0..360 {
        let yaw = (y_deg as f64).to_radians();
        let r_m_v = Rotation3::from_euler_angles(mounting_roll, mounting_pitch, yaw);

        let mut correlation = 0.0;

        // Pick a few high-dynamic segments
        for i in 1..gnss_velocities_ned.len() {
            let (t0, v0) = gnss_velocities_ned[i - 1];
            let (t1, v1) = gnss_velocities_ned[i];
            let dt = t1 - t0;
            if dt <= 0.0 || dt > 1.0 {
                continue;
            }

            let gnss_accel_ned = (v1 - v0) / dt;
            let gnss_accel_mag = f64::sqrt(
                gnss_accel_ned.x * gnss_accel_ned.x + gnss_accel_ned.y * gnss_accel_ned.y,
            );

            if gnss_accel_mag < 0.5 {
                continue;
            } // Need some dynamic force

            // Find matching IMU window
            let imu_seg: Vec<_> = imu_measurements
                .iter()
                .filter(|m| (m.time_tag as f64) >= t0 && (m.time_tag as f64) <= t1)
                .collect();

            if imu_seg.is_empty() {
                continue;
            }

            let mut avg_imu_accel_m = Vector3::zeros();
            for m in &imu_seg {
                avg_imu_accel_m += m.accel;
            }
            avg_imu_accel_m /= imu_seg.len() as f64;

            // Transform IMU to Vehicle Frame at this test Yaw
            let accel_v = r_m_v * avg_imu_accel_m;

            // In a forward-driving vehicle, the forward accel (accel_v.x)
            // should correlate with the magnitude of the horizontal GNSS accel
            // (assuming driving mostly forward).
            // This is a simplified "detection" pass logic.
            correlation += accel_v.x * gnss_accel_mag;
        }

        if correlation > max_corr {
            max_corr = correlation;
            best_yaw = yaw;
        }
    }

    Ok(best_yaw)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gravity_alignment() {
        // Sensor mounted with 10 deg Pitch, 5 deg Roll
        let roll = 5.0f64.to_radians();
        let pitch = 10.0f64.to_radians();
        let r_m_v = Rotation3::from_euler_angles(roll, pitch, 0.0);

        // Gravity in vehicle frame is [0, 0, -9.81]
        let gravity_v = Vector3::new(0.0, 0.0, -9.81);

        // Sensor measures gravity_m = R_v_m * gravity_v
        let gravity_m = r_m_v.inverse() * gravity_v;

        let mut measurements = Vec::new();
        for _ in 0..10 {
            measurements.push(ImuMeasurement {
                time_tag: 0,
                accel: gravity_m,
                gyro: Vector3::zeros(),
                temperature: Some(20.0),
            });
        }

        let (est_roll, est_pitch) = estimate_gravity_alignment(&measurements).unwrap();

        assert!((est_roll - roll).abs() < 1e-6);
        assert!((est_pitch - pitch).abs() < 1e-6);
    }

    #[test]
    fn test_gravity_alignment_empty_measurements() {
        let measurements = Vec::new();
        let result = estimate_gravity_alignment(&measurements);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "No IMU measurements for gravity alignment");
    }

    #[test]
    fn test_gravity_alignment_single_measurement() {
        let measurements = vec![ImuMeasurement {
            time_tag: 0,
            accel: Vector3::new(0.0, 0.0, -9.81),
            gyro: Vector3::zeros(),
            temperature: Some(25.0),
        }];
        let result = estimate_gravity_alignment(&measurements);
        assert!(result.is_ok());
    }

    #[test]
    fn test_gravity_alignment_non_standard_gravity() {
        // Sensor with gravity measured at [5.0, 2.0, -8.0] (not purely vertical)
        // This avoids the atan2(0,0) degenerate case.
        let measurements = vec![ImuMeasurement {
            time_tag: 0,
            accel: Vector3::new(5.0, 2.0, -8.0),
            gyro: Vector3::zeros(),
            temperature: Some(25.0),
        }; 5];
        let (roll, pitch) = estimate_gravity_alignment(&measurements).unwrap();
        // ax = 5.0, ay = 2.0, az = -8.0
        // pitch = atan2(5.0, sqrt(4 + 64)) = atan2(5.0, sqrt(68)) = atan2(5.0, 8.246) ≈ 0.546
        let expected_pitch = f64::atan2(5.0, f64::sqrt(2.0 * 2.0 + (-8.0) * (-8.0)));
        assert!((pitch - expected_pitch).abs() < 1e-10);
        // roll = atan2(-2.0, 8.0) = atan2(-2.0, 8.0) ≈ -0.245
        let expected_roll = f64::atan2(-2.0, 8.0);
        assert!((roll - expected_roll).abs() < 1e-10);
    }

    #[test]
    fn test_heading_alignment_insufficient_data() {
        let imu_measurements = vec![];
        let gnss_velocities_ned = vec![];
        let result = estimate_heading_alignment(&imu_measurements, &gnss_velocities_ned, 0.0, 0.0);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Insufficient dynamic data for heading alignment"
        );
    }

    #[test]
    fn test_mounting_angles_to_rotation() {
        let angles = MountingAngles { roll: 0.1, pitch: 0.2, yaw: 0.3 };
        let r = angles.to_rotation();
        let v = Vector3::new(1.0, 0.0, 0.0);
        let rotated = angles.apply(v);
        let expected = r * v;
        assert!((rotated - expected).norm() < 1e-10);
    }

    #[test]
    fn test_mounting_angles_zero() {
        let angles = MountingAngles { roll: 0.0, pitch: 0.0, yaw: 0.0 };
        let v = Vector3::new(1.0, 2.0, 3.0);
        let rotated = angles.apply(v);
        assert!((rotated - v).norm() < 1e-10);
    }

    #[test]
    fn test_heading_alignment_success() {
        let true_yaw = 0.5; // radians

        // Generate GNSS velocities at 1 Hz, accelerating along North for 5 s then constant
        let mut gnss_vels: Vec<(f64, Vector3<f64>)> = Vec::new();
        for i in 0..11 {
            let t = i as f64;
            let v_north = if i < 6 { (i as f64) * 2.0 } else { 10.0 };
            gnss_vels.push((t, Vector3::new(v_north, 0.0, 0.0)));
        }

        // Generate IMU measurements at 10 Hz (100 points over 10 seconds)
        let mut imu_meas = Vec::new();
        let r_true = Rotation3::from_euler_angles(0.0, 0.0, true_yaw);
        for i in 0..100 {
            let t = i as f64 * 0.1;
            let time_tag = t as u32; // truncate to seconds to match GNSS time units
            // Vehicle forward acceleration: 2 m/s^2 while accelerating, 0 at constant speed
            let accel_forward = if t < 5.0 { 2.0 } else { 0.0 };
            let accel_vehicle = Vector3::new(accel_forward, 0.0, 0.0);
            // IMU measures in mounting frame: rotated by inverse of true yaw
            let accel_mount = r_true.inverse() * accel_vehicle;
            imu_meas.push(ImuMeasurement {
                time_tag,
                accel: accel_mount,
                gyro: Vector3::zeros(),
                temperature: Some(25.0),
            });
        }

        let est_yaw =
            estimate_heading_alignment(&imu_meas, &gnss_vels, 0.0, 0.0).unwrap();

        // The search is in 1-degree increments, so accuracy is ~0.017 rad
        assert!(
            (est_yaw - true_yaw).abs() < 0.02,
            "Expected yaw ~{}, got {}",
            true_yaw,
            est_yaw
        );
    }

    #[test]
    fn test_heading_alignment_large_dt_skipped() {
        // dt > 1.0 between GNSS epochs should be skipped (no correlation contribution)
        let true_yaw = 0.3;

        let mut gnss_vels: Vec<(f64, Vector3<f64>)> = Vec::new();
        // dt=2.0 between these epochs (skipped because dt > 1.0)
        for i in 0..11 {
            gnss_vels.push((i as f64 * 2.0, Vector3::new(10.0, 0.0, 0.0)));
        }

        let mut imu_meas = Vec::new();
        let r_true = Rotation3::from_euler_angles(0.0, 0.0, true_yaw);
        for i in 0..100 {
            let t = i as f64 * 0.1;
            let time_tag = t as u32;
            let accel_vehicle = Vector3::new(2.0, 0.0, 0.0);
            let accel_mount = r_true.inverse() * accel_vehicle;
            imu_meas.push(ImuMeasurement {
                time_tag,
                accel: accel_mount,
                gyro: Vector3::zeros(),
                temperature: Some(25.0),
            });
        }

        // dt > 1.0 for all pairs -> correlation stays at 0 for all yaw candidates
        // best_yaw will be 0 (default, since max_corr stays at -1 and never updates)
        let est_yaw = estimate_heading_alignment(&imu_meas, &gnss_vels, 0.0, 0.0).unwrap();
        assert!((est_yaw - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_heading_alignment_low_dynamics_skipped() {
        // gnss_accel_mag < 0.5 should be skipped
        let true_yaw = 0.3;

        let mut gnss_vels: Vec<(f64, Vector3<f64>)> = Vec::new();
        // Small velocity differences -> accel < 0.5
        for i in 0..11 {
            gnss_vels.push((i as f64, Vector3::new(0.1, 0.0, 0.0)));
        }

        let mut imu_meas = Vec::new();
        let r_true = Rotation3::from_euler_angles(0.0, 0.0, true_yaw);
        for i in 0..100 {
            let t = i as f64 * 0.1;
            let time_tag = t as u32;
            let accel_vehicle = Vector3::new(2.0, 0.0, 0.0);
            let accel_mount = r_true.inverse() * accel_vehicle;
            imu_meas.push(ImuMeasurement {
                time_tag,
                accel: accel_mount,
                gyro: Vector3::zeros(),
                temperature: Some(25.0),
            });
        }

        // All GNSS pairs have low accel -> skipped -> correlation stays 0
        let est_yaw = estimate_heading_alignment(&imu_meas, &gnss_vels, 0.0, 0.0).unwrap();
        assert!((est_yaw - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_heading_alignment_time_jump_skipped() {
        // Negative dt should be skipped; need at least 10 entries
        let mut gnss_vels: Vec<(f64, Vector3<f64>)> = Vec::new();
        // Most pairs have normal forward time
        for i in 0..10 {
            gnss_vels.push((i as f64, Vector3::new(i as f64 * 2.0, 0.0, 0.0)));
        }
        // Add a pair with time going backwards
        gnss_vels.push((100.0, Vector3::new(20.0, 0.0, 0.0)));
        gnss_vels.push((5.0, Vector3::new(22.0, 0.0, 0.0))); // time goes backwards

        let mut imu_meas = Vec::new();
        for i in 0..100 {
            let t = i as f64 * 0.1;
            imu_meas.push(ImuMeasurement {
                time_tag: t as u32,
                accel: Vector3::new(1.0, 0.0, 0.0),
                gyro: Vector3::zeros(),
                temperature: Some(25.0),
            });
        }

        let result = estimate_heading_alignment(&imu_meas, &gnss_vels, 0.0, 0.0);
        assert!(result.is_ok());
    }

    #[test]
    fn test_mounting_angles_apply_negative() {
        let angles = MountingAngles { roll: -0.1, pitch: 0.2, yaw: -0.3 };
        let v = Vector3::new(1.0, -2.0, 3.0);
        let rotated = angles.apply(v);
        let expected = angles.to_rotation() * v;
        assert!((rotated - expected).norm() < 1e-10);
    }

    #[test]
    fn test_gravity_alignment_zero_rotation() {
        // Sensor perfectly aligned with vehicle frame — gravity reads [0, 0, -g]
        let measurements = vec![ImuMeasurement {
            time_tag: 0,
            accel: Vector3::new(0.0, 0.0, -9.81),
            gyro: Vector3::zeros(),
            temperature: Some(25.0),
        }; 10];
        let (roll, pitch) = estimate_gravity_alignment(&measurements).unwrap();
        assert!((roll - 0.0).abs() < 1e-10, "zero roll for perfect alignment");
        assert!((pitch - 0.0).abs() < 1e-10, "zero pitch for perfect alignment");
    }

    #[test]
    fn test_gravity_alignment_all_positive_accel() {
        // All positive accelerations — edge case for atan2
        let measurements = vec![ImuMeasurement {
            time_tag: 0,
            accel: Vector3::new(1.0, 2.0, 3.0),
            gyro: Vector3::zeros(),
            temperature: Some(25.0),
        }; 5];
        let (roll, pitch) = estimate_gravity_alignment(&measurements).unwrap();
        let expected_pitch = f64::atan2(1.0, f64::sqrt(4.0 + 9.0));
        let expected_roll = f64::atan2(-2.0, -3.0);
        assert!((pitch - expected_pitch).abs() < 1e-10);
        assert!((roll - expected_roll).abs() < 1e-10);
    }

    #[test]
    fn test_heading_alignment_single_gnss_epoch_insufficient() {
        // Fewer than 10 GNSS velocities -> error
        let imu_meas = vec![ImuMeasurement {
            time_tag: 0,
            accel: Vector3::new(1.0, 0.0, 0.0),
            gyro: Vector3::zeros(),
            temperature: None,
        }; 100];
        let gnss_vels = vec![(0.0, Vector3::new(10.0, 0.0, 0.0))]; // only 1
        let result = estimate_heading_alignment(&imu_meas, &gnss_vels, 0.0, 0.0);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Insufficient dynamic data for heading alignment"
        );
    }

    #[test]
    fn test_heading_alignment_fewer_than_100_imu_insufficient() {
        let gnss_vels = vec![(0.0, Vector3::new(10.0, 0.0, 0.0)); 11];
        let imu_meas = vec![ImuMeasurement {
            time_tag: 0,
            accel: Vector3::new(1.0, 0.0, 0.0),
            gyro: Vector3::zeros(),
            temperature: None,
        }; 99]; // fewer than 100
        let result = estimate_heading_alignment(&imu_meas, &gnss_vels, 0.0, 0.0);
        assert!(result.is_err());
    }
}
