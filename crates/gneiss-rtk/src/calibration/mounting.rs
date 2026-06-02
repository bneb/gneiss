use nalgebra::{Vector3, Rotation3};
use gneiss_core::imu::ImuMeasurement;

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
pub fn estimate_gravity_alignment(measurements: &[ImuMeasurement]) -> Result<(f64, f64), &'static str> {
    if measurements.is_empty() { return Err("No IMU measurements for gravity alignment"); }
    
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
    let pitch = f64::atan2(ax, f64::sqrt(ay*ay + az*az));
    
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
            let (t0, v0) = gnss_velocities_ned[i-1];
            let (t1, v1) = gnss_velocities_ned[i];
            let dt = t1 - t0;
            if dt <= 0.0 || dt > 1.0 { continue; }
            
            let gnss_accel_ned = (v1 - v0) / dt;
            let gnss_accel_mag = f64::sqrt(gnss_accel_ned.x * gnss_accel_ned.x + gnss_accel_ned.y * gnss_accel_ned.y);
            
            if gnss_accel_mag < 0.5 { continue; } // Need some dynamic force

            // Find matching IMU window
            let imu_seg: Vec<_> = imu_measurements.iter()
                .filter(|m| (m.time_tag as f64) >= t0 && (m.time_tag as f64) <= t1)
                .collect();
            
            if imu_seg.is_empty() { continue; }
            
            let mut avg_imu_accel_m = Vector3::zeros();
            for m in &imu_seg { avg_imu_accel_m += m.accel; }
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
}
