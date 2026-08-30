//! Applanix POSPac Smoothed Best Estimate of Trajectory (SBET) binary reader & writer.
//!
//! SBET (.sbet / .out): 17 IEEE-754 64-bit float fields = 136 bytes per record (little-endian).
//! SBET RMS (.smrms / .sbet.rms): 10 IEEE-754 64-bit float fields = 80 bytes per record (little-endian).

#[cfg(test)]
mod tests;

use std::io::{self, Read, Write};

/// Standard 17-field Applanix SBET trajectory record (136 bytes).
#[derive(Debug, Clone, PartialEq)]
pub struct SbetRecord {
    /// GPS time of week (seconds).
    pub time: f64,
    /// Latitude (radians).
    pub latitude: f64,
    /// Longitude (radians).
    pub longitude: f64,
    /// Ellipsoidal altitude (meters).
    pub altitude: f64,
    /// East velocity (m/s).
    pub x_vel: f64,
    /// North velocity (m/s).
    pub y_vel: f64,
    /// Up velocity (m/s).
    pub z_vel: f64,
    /// Roll angle (radians).
    pub roll: f64,
    /// Pitch angle (radians).
    pub pitch: f64,
    /// Heading / Azimuth angle (radians clockwise from True North).
    pub heading: f64,
    /// Wander angle (radians).
    pub wander_angle: f64,
    /// Body X acceleration (m/s^2).
    pub x_accel: f64,
    /// Body Y acceleration (m/s^2).
    pub y_accel: f64,
    /// Body Z acceleration (m/s^2).
    pub z_accel: f64,
    /// Body X angular rate (rad/s).
    pub x_ang_rate: f64,
    /// Body Y angular rate (rad/s).
    pub y_ang_rate: f64,
    /// Body Z angular rate (rad/s).
    pub z_ang_rate: f64,
}

impl SbetRecord {
    pub const RECORD_SIZE_BYTES: usize = 136;

    /// Write this record as 136 bytes (17 little-endian f64 values).
    pub fn write_to<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let fields = [
            self.time,
            self.latitude,
            self.longitude,
            self.altitude,
            self.x_vel,
            self.y_vel,
            self.z_vel,
            self.roll,
            self.pitch,
            self.heading,
            self.wander_angle,
            self.x_accel,
            self.y_accel,
            self.z_accel,
            self.x_ang_rate,
            self.y_ang_rate,
            self.z_ang_rate,
        ];
        for f in fields {
            w.write_all(&f.to_le_bytes())?;
        }
        Ok(())
    }

    /// Read a single record from 136 bytes.
    pub fn read_from<R: Read>(r: &mut R) -> io::Result<Self> {
        let mut buf = [0u8; Self::RECORD_SIZE_BYTES];
        r.read_exact(&mut buf)?;

        let mut f = [0.0f64; 17];
        for i in 0..17 {
            let chunk: [u8; 8] = buf[i * 8..(i + 1) * 8].try_into().expect("valid chunk slice");
            f[i] = f64::from_le_bytes(chunk);
        }

        Ok(Self {
            time: f[0],
            latitude: f[1],
            longitude: f[2],
            altitude: f[3],
            x_vel: f[4],
            y_vel: f[5],
            z_vel: f[6],
            roll: f[7],
            pitch: f[8],
            heading: f[9],
            wander_angle: f[10],
            x_accel: f[11],
            y_accel: f[12],
            z_accel: f[13],
            x_ang_rate: f[14],
            y_ang_rate: f[15],
            z_ang_rate: f[16],
        })
    }
}

/// Standard 10-field Applanix SBET companion RMS accuracy record (80 bytes).
#[derive(Debug, Clone, PartialEq)]
pub struct SbetRmsRecord {
    /// GPS time of week (seconds).
    pub time: f64,
    /// North position standard deviation / RMS (meters).
    pub north_pos_rms: f64,
    /// East position standard deviation / RMS (meters).
    pub east_pos_rms: f64,
    /// Down / Vertical position standard deviation / RMS (meters).
    pub down_pos_rms: f64,
    /// North velocity standard deviation / RMS (m/s).
    pub north_vel_rms: f64,
    /// East velocity standard deviation / RMS (m/s).
    pub east_vel_rms: f64,
    /// Down velocity standard deviation / RMS (m/s).
    pub down_vel_rms: f64,
    /// Roll attitude standard deviation / RMS (radians).
    pub roll_rms: f64,
    /// Pitch attitude standard deviation / RMS (radians).
    pub pitch_rms: f64,
    /// Heading attitude standard deviation / RMS (radians).
    pub heading_rms: f64,
}

impl SbetRmsRecord {
    pub const RECORD_SIZE_BYTES: usize = 80;

    /// Write this record as 80 bytes (10 little-endian f64 values).
    pub fn write_to<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let fields = [
            self.time,
            self.north_pos_rms,
            self.east_pos_rms,
            self.down_pos_rms,
            self.north_vel_rms,
            self.east_vel_rms,
            self.down_vel_rms,
            self.roll_rms,
            self.pitch_rms,
            self.heading_rms,
        ];
        for f in fields {
            w.write_all(&f.to_le_bytes())?;
        }
        Ok(())
    }

    /// Read a single record from 80 bytes.
    pub fn read_from<R: Read>(r: &mut R) -> io::Result<Self> {
        let mut buf = [0u8; Self::RECORD_SIZE_BYTES];
        r.read_exact(&mut buf)?;

        let mut f = [0.0f64; 10];
        for i in 0..10 {
            let chunk: [u8; 8] = buf[i * 8..(i + 1) * 8].try_into().expect("valid chunk slice");
            f[i] = f64::from_le_bytes(chunk);
        }

        Ok(Self {
            time: f[0],
            north_pos_rms: f[1],
            east_pos_rms: f[2],
            down_pos_rms: f[3],
            north_vel_rms: f[4],
            east_vel_rms: f[5],
            down_vel_rms: f[6],
            roll_rms: f[7],
            pitch_rms: f[8],
            heading_rms: f[9],
        })
    }
}
