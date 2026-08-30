#![allow(clippy::unwrap_used)]

use super::*;
use std::io::Cursor;

#[test]
fn test_sbet_record_binary_roundtrip() {
    let rec = SbetRecord {
        time: 345600.0,
        latitude: 0.654321,
        longitude: -2.123456,
        altitude: 125.75,
        x_vel: 1.5,
        y_vel: -0.8,
        z_vel: 0.05,
        roll: 0.012,
        pitch: -0.005,
        heading: std::f64::consts::FRAC_PI_2,
        wander_angle: 0.0,
        x_accel: 0.02,
        y_accel: 0.01,
        z_accel: 9.806,
        x_ang_rate: 0.001,
        y_ang_rate: -0.002,
        z_ang_rate: 0.0005,
    };

    let mut buf = Vec::new();
    rec.write_to(&mut buf).expect("write record");
    assert_eq!(buf.len(), 136);

    let mut cursor = Cursor::new(buf);
    let decoded = SbetRecord::read_from(&mut cursor).expect("read record");
    assert_eq!(rec, decoded);
}

#[test]
fn test_sbet_rms_record_binary_roundtrip() {
    let rms = SbetRmsRecord {
        time: 345600.0,
        north_pos_rms: 0.008,
        east_pos_rms: 0.007,
        down_pos_rms: 0.015,
        north_vel_rms: 0.002,
        east_vel_rms: 0.002,
        down_vel_rms: 0.005,
        roll_rms: 0.0001,
        pitch_rms: 0.0001,
        heading_rms: 0.0005,
    };

    let mut buf = Vec::new();
    rms.write_to(&mut buf).expect("write rms");
    assert_eq!(buf.len(), 80);

    let mut cursor = Cursor::new(buf);
    let decoded = SbetRmsRecord::read_from(&mut cursor).expect("read rms");
    assert_eq!(rms, decoded);
}
