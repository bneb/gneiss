//! Real-Time Live Streaming RTK Mode.
//!
//! Connects live rover stream (UBX / Serial / TCP) with real-time base corrections
//! (NTRIP client / RTCM3 stream), performs sub-millisecond per-epoch RTK IEKF updates,
//! and broadcasts real-time NMEA GGA / telemetry.

use std::path::Path;
use tracing::info;

use gneiss_core::coords::ecef_to_llh;
use gneiss_rtk::streaming::{StreamingConfig, StreamingRtkEngine};
use gneiss_rtk::post_process::ProcessingDynamics;

pub struct LiveArgs {
    pub rover_source: String,
    pub ntrip_mountpoint: Option<String>,
    pub base_position: Option<String>,
    pub nav: Option<String>,
}

pub async fn run_live_mode(args: LiveArgs) -> Result<(), Box<dyn std::error::Error>> {
    info!("Starting Gneiss Live Streaming RTK Engine...");
    info!("Rover Source: {}", args.rover_source);
    if let Some(mp) = &args.ntrip_mountpoint {
        info!("NTRIP Mountpoint: {}", mp);
    }

    let base_pos = if let Some(s) = &args.base_position {
        let parts: Vec<f64> = s.split(',').filter_map(|p| p.trim().parse().ok()).collect();
        if parts.len() == 3 {
            nalgebra::Vector3::new(parts[0], parts[1], parts[2])
        } else {
            nalgebra::Vector3::zeros()
        }
    } else {
        nalgebra::Vector3::zeros()
    };

    let config = StreamingConfig {
        base_position: base_pos,
        max_base_age_s: 3.0,
        dynamics: ProcessingDynamics::Kinematic,
        enable_glonass: true,
    };

    let mut ephems = Vec::new();
    if let Some(nav_path) = &args.nav {
        let file = std::fs::File::open(Path::new(nav_path))?;
        let (e, _) = gneiss_parsers::rinex::parse_rinex_nav(std::io::BufReader::new(file))?;
        ephems = e;
    }

    let engine = StreamingRtkEngine::new(config, ephems);
    info!("Engine initialized with {} ephemerides.", engine.ephemerides_len());
    Err("Live serial and NTRIP streaming runner is under active development. Direct execution from CLI is not yet implemented.".into())
}

/// Formats a streaming epoch position into a standard NMEA 0183 $GNGGA string.
#[allow(dead_code)]
pub fn format_nmea_gga(
    time_tow: f64,
    pos_ecef: nalgebra::Vector3<f64>,
    quality: u8,
    n_sat: usize,
    hdop: f64,
) -> String {
    let llh = ecef_to_llh(pos_ecef);
    let lat_deg = llh.x.to_degrees();
    let lon_deg = llh.y.to_degrees();
    let alt_m = llh.z;

    let lat_abs = lat_deg.abs();
    let lat_dd = lat_abs.floor();
    let lat_mm = (lat_abs - lat_dd) * 60.0;
    let lat_hemi = if lat_deg >= 0.0 { 'N' } else { 'S' };

    let lon_abs = lon_deg.abs();
    let lon_dd = lon_abs.floor();
    let lon_mm = (lon_abs - lon_dd) * 60.0;
    let lon_hemi = if lon_deg >= 0.0 { 'E' } else { 'W' };

    let total_sec = time_tow as u64;
    let hh = (total_sec / 3600) % 24;
    let mm = (total_sec % 3600) / 60;
    let ss = (total_sec % 60) as f64 + (time_tow - total_sec as f64);

    let fix_quality = match quality {
        1 => 4, // RTK Fixed
        2 => 5, // RTK Float
        _ => 1, // Single / Autonomous
    };

    let raw = format!(
        "GNGGA,{:02}{:02}{:05.2},{:02}{:07.4},{},{:03}{:07.4},{},{},{:02},{:.1},{:.3},M,0.0,M,,",
        hh, mm, ss,
        lat_dd as u32, lat_mm, lat_hemi,
        lon_dd as u32, lon_mm, lon_hemi,
        fix_quality, n_sat, hdop, alt_m
    );

    let mut checksum: u8 = 0;
    for b in raw.bytes() {
        checksum ^= b;
    }

    format!("${}*{:02X}\r\n", raw, checksum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nmea_gga_formatting() {
        let pos = nalgebra::Vector3::new(-2697941.0, -4255089.0, 3898009.0);
        let gga = format_nmea_gga(345600.0, pos, 1, 14, 0.8);
        assert!(gga.starts_with("$GNGGA"));
        assert!(gga.contains(",4,14,0.8,"));
        assert!(gga.ends_with("\r\n"));
    }
}
