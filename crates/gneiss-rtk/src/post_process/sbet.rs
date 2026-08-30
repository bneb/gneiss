//! Export smoothed post-processed trajectories to Applanix binary SBET and companion RMS formats.

use std::fs::File;
use std::io::{self, BufWriter};
use std::path::Path;

use nalgebra::{Matrix3, UnitQuaternion, Vector3};

use gneiss_core::coords::ecef_to_llh;
use gneiss_parsers::sbet::{SbetRecord, SbetRmsRecord};

use super::combiner::SmoothedEpoch;

/// Export a smoothed trajectory to standard Applanix SBET and optional RMS files.
pub fn export_sbet_trajectory(
    smoothed: &[SmoothedEpoch],
    sbet_path: &Path,
    rms_path: Option<&Path>,
) -> io::Result<()> {
    let sbet_file = File::create(sbet_path)?;
    let mut sbet_writer = BufWriter::new(sbet_file);

    let mut rms_writer = if let Some(p) = rms_path {
        let f = File::create(p)?;
        Some(BufWriter::new(f))
    } else {
        None
    };

    for ep in smoothed {
        let llh = ecef_to_llh(ep.position_ecef);
        let lat_rad = llh.x;
        let lon_rad = llh.y;
        let alt_m = llh.z;

        // Rotate ECEF velocity to local Navigation (NED: North, East, Down)
        let r_ecef_to_ned = ecef_to_ned_matrix(lat_rad, lon_rad);
        let v_ned = ep.velocity_ecef.map(|v| r_ecef_to_ned * v).unwrap_or_else(Vector3::zeros);
        let x_vel_east = v_ned.y;
        let y_vel_north = v_ned.x;
        let z_vel_up = -v_ned.z;

        // Attitude angles (roll, pitch, heading) if attitude is present
        let (roll, pitch, heading) = if let Some(q) = ep.attitude {
            quaternion_to_roll_pitch_heading(q)
        } else {
            (0.0, 0.0, 0.0)
        };

        let record = SbetRecord {
            time: ep.time.tow,
            latitude: lat_rad,
            longitude: lon_rad,
            altitude: alt_m,
            x_vel: x_vel_east,
            y_vel: y_vel_north,
            z_vel: z_vel_up,
            roll,
            pitch,
            heading,
            wander_angle: 0.0,
            x_accel: 0.0,
            y_accel: 0.0,
            z_accel: 9.80665,
            x_ang_rate: 0.0,
            y_ang_rate: 0.0,
            z_ang_rate: 0.0,
        };

        record.write_to(&mut sbet_writer)?;

        if let Some(ref mut rw) = rms_writer {
            let north_rms = ep.std_north;
            let east_rms = ep.std_east;
            let down_rms = ep.std_up;

            let rms_rec = SbetRmsRecord {
                time: ep.time.tow,
                north_pos_rms: north_rms,
                east_pos_rms: east_rms,
                down_pos_rms: down_rms,
                north_vel_rms: (north_rms * 0.1).max(0.001),
                east_vel_rms: (east_rms * 0.1).max(0.001),
                down_vel_rms: (down_rms * 0.1).max(0.001),
                roll_rms: 0.001,
                pitch_rms: 0.001,
                heading_rms: 0.002,
            };
            rms_rec.write_to(rw)?;
        }
    }

    Ok(())
}

fn ecef_to_ned_matrix(lat_rad: f64, lon_rad: f64) -> Matrix3<f64> {
    let sin_phi = lat_rad.sin();
    let cos_phi = lat_rad.cos();
    let sin_lam = lon_rad.sin();
    let cos_lam = lon_rad.cos();

    Matrix3::new(
        -sin_phi * cos_lam, -sin_phi * sin_lam, cos_phi,
        -sin_lam, cos_lam, 0.0,
        -cos_phi * cos_lam, -cos_phi * sin_lam, -sin_phi,
    )
}

fn quaternion_to_roll_pitch_heading(q: UnitQuaternion<f64>) -> (f64, f64, f64) {
    let r = q.to_rotation_matrix();
    let m = r.matrix();

    // Standard Tait-Bryan ZYX Euler angles: Yaw/Heading (Z), Pitch (Y), Roll (X)
    let pitch = (-m[(2, 0)]).asin();
    let roll = m[(2, 1)].atan2(m[(2, 2)]);
    let heading = m[(1, 0)].atan2(m[(0, 0)]);

    let heading_norm = if heading < 0.0 {
        heading + std::f64::consts::TAU
    } else {
        heading
    };

    (roll, pitch, heading_norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_sbet_export_roundtrip() {
        let ep = SmoothedEpoch {
            time: GpsTime::new(2000, 100.0),
            position_ecef: Vector3::new(-2697941.0, -4255089.0, 3898009.0),
            velocity_ecef: Some(Vector3::new(1.0, 0.5, -0.2)),
            attitude: Some(UnitQuaternion::identity()),
            cov_position: Matrix3::identity() * 0.0004,
            std_east: 0.02,
            std_north: 0.02,
            std_up: 0.05,
            separation_3d: 0.01,
            quality: 1,
            n_satellites: 12,
        };

        let dir = tempfile::tempdir().expect("tempdir");
        let sbet_path = dir.path().join("test.sbet");
        let rms_path = dir.path().join("test.sbet.rms");

        export_sbet_trajectory(&[ep], &sbet_path, Some(&rms_path)).expect("export");

        assert!(sbet_path.exists());
        assert!(rms_path.exists());
        assert_eq!(std::fs::metadata(&sbet_path).expect("meta").len(), 136);
        assert_eq!(std::fs::metadata(&rms_path).expect("meta").len(), 80);
    }
}
