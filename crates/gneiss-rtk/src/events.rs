//! High-precision UAV / Photogrammetry Camera Shutter Event Mark Interpolator.
//!
//! Synchronizes camera exposure hardware time marks (`MARK3`, `TIMEMARK`, UBX `TIM_TM2`)
//! with the post-processed trajectory, applying antenna-to-camera lever-arm offsets
//! and interpolating exact camera position, velocity, and covariance at exposure epoch.

use gneiss_core::coords::ecef_to_llh;
use gneiss_core::gnss_time::GnssTime;
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

/// A single trajectory epoch point from the PPK forward/backward smoother.
#[derive(Debug, Clone)]
pub struct TrajectoryEpoch {
    /// GPS time of the epoch.
    pub time: GnssTime,
    /// ECEF position of the antenna (meters).
    pub pos_ecef: Vector3<f64>,
    /// ECEF velocity of the antenna (m/s).
    pub vel_ecef: Vector3<f64>,
    /// Diagonal standard deviations East, North, Up (meters).
    pub std_enu: Vector3<f64>,
    /// Optional attitude rotation matrix R_body_to_ecef.
    pub r_body_ecef: Option<Matrix3<f64>>,
}

/// Photogrammetry camera event configuration.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CameraEventConfig {
    /// Lever-arm vector from GNSS antenna to camera optical center in body frame [Forward, Right, Down] (meters).
    pub lever_arm_body: [f64; 3],
    /// Shutter electronic delay offset (seconds, added to event timestamp).
    pub shutter_delay_s: f64,
}

impl Default for CameraEventConfig {
    fn default() -> Self {
        Self {
            lever_arm_body: [0.0, 0.0, 0.0],
            shutter_delay_s: 0.0,
        }
    }
}

/// An interpolated camera exposure event center.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CameraEventRecord {
    /// Event sequence index (1-based).
    pub event_id: usize,
    /// GPS time of the camera exposure.
    pub time_gpst_s: f64,
    /// Camera optical center ECEF X, Y, Z (meters).
    pub pos_ecef: [f64; 3],
    /// Geodetic latitude (degrees).
    pub lat_deg: f64,
    /// Geodetic longitude (degrees).
    pub lon_deg: f64,
    /// Ellipsoidal height (meters).
    pub height_m: f64,
    /// Standard deviations East, North, Up (meters).
    pub std_enu: [f64; 3],
}

/// Interpolates camera exposure events across a continuous trajectory.
pub struct CameraEventInterpolator;

impl CameraEventInterpolator {
    /// Interpolates camera position and covariance for each event timestamp.
    pub fn interpolate_events(
        trajectory: &[TrajectoryEpoch],
        event_times: &[GnssTime],
        config: &CameraEventConfig,
    ) -> Vec<CameraEventRecord> {
        if trajectory.len() < 2 || event_times.is_empty() {
            return Vec::new();
        }

        let mut results = Vec::with_capacity(event_times.len());

        for (idx, &event_time) in event_times.iter().enumerate() {
            let t_eff = event_time.to_gpst().tow + config.shutter_delay_s;
            if let Some(record) = Self::interpolate_single(trajectory, t_eff, idx + 1, config) {
                results.push(record);
            }
        }

        results
    }

    fn interpolate_single(
        traj: &[TrajectoryEpoch],
        t_target: f64,
        event_id: usize,
        config: &CameraEventConfig,
    ) -> Option<CameraEventRecord> {
        // Find surrounding trajectory interval [k, k+1]
        let k = traj.windows(2).position(|w| {
            let t0 = w[0].time.to_gpst().tow;
            let t1 = w[1].time.to_gpst().tow;
            t_target >= t0 && t_target <= t1
        })?;

        let e0 = &traj[k];
        let e1 = &traj[k + 1];
        let t0 = e0.time.to_gpst().tow;
        let t1 = e1.time.to_gpst().tow;
        let dt_span = t1 - t0;

        if dt_span <= 1e-6 {
            return None;
        }

        let dt = t_target - t0;
        let frac = dt / dt_span;

        // Cubic Hermite trajectory interpolation
        let h00 = 2.0 * frac * frac * frac - 3.0 * frac * frac + 1.0;
        let h10 = frac * frac * frac - 2.0 * frac * frac + frac;
        let h01 = -2.0 * frac * frac * frac + 3.0 * frac * frac;
        let h11 = frac * frac * frac - frac * frac;

        let interp_ant_ecef = h00 * e0.pos_ecef
            + h10 * dt_span * e0.vel_ecef
            + h01 * e1.pos_ecef
            + h11 * dt_span * e1.vel_ecef;

        let lever_vec = Vector3::new(
            config.lever_arm_body[0],
            config.lever_arm_body[1],
            config.lever_arm_body[2],
        );

        // Apply lever-arm rotation if body attitude available; otherwise direct offset
        let r_ecef = if let Some(r0) = e0.r_body_ecef {
            r0 * lever_vec
        } else {
            lever_vec
        };

        let camera_pos_ecef = interp_ant_ecef + r_ecef;
        let llh = ecef_to_llh(camera_pos_ecef);

        // Linear interpolation of uncertainties
        let std_enu = (1.0 - frac) * e0.std_enu + frac * e1.std_enu;

        Some(CameraEventRecord {
            event_id,
            time_gpst_s: t_target,
            pos_ecef: [camera_pos_ecef.x, camera_pos_ecef.y, camera_pos_ecef.z],
            lat_deg: llh[0].to_degrees(),
            lon_deg: llh[1].to_degrees(),
            height_m: llh[2],
            std_enu: [std_enu.x, std_enu.y, std_enu.z],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::gnss_time::TimeSystem;
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_camera_event_interpolation_exact_midpoint() {
        let t0 = GnssTime::from_gpst(TimeSystem::Gps, GpsTime::new(2000, 100.0));
        let t1 = GnssTime::from_gpst(TimeSystem::Gps, GpsTime::new(2000, 102.0));

        let p0 = Vector3::new(1000.0, 2000.0, 3000.0);
        let p1 = Vector3::new(1020.0, 2040.0, 3000.0);
        let v = (p1 - p0) / 2.0;

        let traj = vec![
            TrajectoryEpoch {
                time: t0,
                pos_ecef: p0,
                vel_ecef: v,
                std_enu: Vector3::new(0.01, 0.01, 0.02),
                r_body_ecef: None,
            },
            TrajectoryEpoch {
                time: t1,
                pos_ecef: p1,
                vel_ecef: v,
                std_enu: Vector3::new(0.01, 0.01, 0.02),
                r_body_ecef: None,
            },
        ];

        let config = CameraEventConfig {
            lever_arm_body: [0.1, 0.0, -0.2],
            shutter_delay_s: 0.0,
        };

        let event_time = GnssTime::from_gpst(TimeSystem::Gps, GpsTime::new(2000, 101.0)); // Exact midpoint
        let records = CameraEventInterpolator::interpolate_events(&traj, &[event_time], &config);

        assert_eq!(records.len(), 1);
        let rec = &records[0];
        assert_eq!(rec.event_id, 1);
        // Midpoint pos = (1010.0 + 0.1, 2020.0 + 0.0, 3000.0 - 0.2)
        assert!((rec.pos_ecef[0] - 1010.1).abs() < 1e-4);
        assert!((rec.pos_ecef[1] - 2020.0).abs() < 1e-4);
        assert!((rec.pos_ecef[2] - 2999.8).abs() < 1e-4);
    }
}
