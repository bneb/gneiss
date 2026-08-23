//! SWFG accuracy evaluation on Odaiba dataset.
//!
//! Usage: cargo run --release --bin eval_odaiba

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::engine::SwfgEngine;

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("warn")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dataset = Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba");

    let nav_file = File::open(dataset.join("base.nav")).expect("open base.nav");
    let (ephemerides, klobuchar) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_file))
        .expect("parse nav");

    let ref_csv = std::fs::read_to_string(dataset.join("reference.csv")).expect("read reference.csv");
    let mut truth: BTreeMap<u32, (f64, f64, f64)> = BTreeMap::new();
    for line in ref_csv.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 8 { continue; }
        if let (Ok(tow), Ok(x), Ok(y), Ok(z)) = (
            parts[0].trim().parse::<f64>(),
            parts[5].trim().parse::<f64>(),
            parts[6].trim().parse::<f64>(),
            parts[7].trim().parse::<f64>(),
        ) {
            truth.insert(tow.floor() as u32, (x, y, z));
        }
    }

    let rover_file = File::open(dataset.join("rover_trimble.obs")).expect("open rover_trimble.obs");
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rover_file))
        .expect("parse rover RINEX");

    let base_file = File::open(dataset.join("base_trimble.obs")).expect("open base_trimble.obs");
    let (base_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_file))
        .expect("parse base RINEX");

    let base_index: BTreeMap<u32, &gneiss_core::obs::EpochObs> = base_epochs.iter()
        .map(|e| ((e.time.tow * 1000.0) as u32, e))
        .collect();

    let base_arp = nalgebra::Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let base_llh = gneiss_core::coords::ecef_to_llh(base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(base_llh).transpose();
    // Antenna delta H/E/N = -0.0855m H (Down = -H = +0.0855m)
    let base_pos = base_arp + ned_to_ecef * nalgebra::Vector3::new(0.0, 0.0, 0.0855);

    let config = EngineConfig::Rtk(gneiss_rtk::swfg::config::RtkConfig::default());
    let mut engine = SwfgEngine::new(&config, ephemerides);
    if let Some(ref k) = klobuchar {
        engine.set_klobuchar(k.alpha, k.beta);
    }

    let mut h_errors = Vec::new();
    let mut err_count = 0usize;
    let mut _ok_count = 0usize;
    let mut processed = 0usize;
    for epoch in &rover_epochs {
        if processed >= 500 { break; }
        let tow_sec = epoch.time.tow.floor() as u32;
        let exact_ms = (epoch.time.tow * 1000.0).round() as u32;
        let base_epoch = match base_index.get(&exact_ms).copied() {
            Some(b) => b,
            None => continue,
        };

        let sol_res = engine.process_rtk_epoch(epoch, base_epoch, base_pos);

        match sol_res {
            Ok(sol) => {
                _ok_count += 1;
                processed += 1;
                if let Some((tx, ty, tz)) = truth.get(&tow_sec).copied() {
                    let truth_ecef = nalgebra::Vector3::new(tx, ty, tz);
                    let truth_llh = gneiss_core::coords::ecef_to_llh(truth_ecef);
                    let ned = gneiss_core::coords::ecef_to_ned_matrix(truth_llh) * (sol.position_ecef - truth_ecef);
                    let h = (ned.x * ned.x + ned.y * ned.y).sqrt();

                    println!("MATCH tow={} swfg=({:.0},{:.0},{:.0}) truth=({:.0},{:.0},{:.0}) err={:.1}",
                        tow_sec, sol.position_ecef.x, sol.position_ecef.y, sol.position_ecef.z, tx, ty, tz, h);

                    if h < 500.0 {
                        h_errors.push(h);
                    }
                }
            }
            Err(e) => {
                err_count += 1;
                if err_count <= 3 { println!("ERR: {}", e); }
            }
        }
    }

    let n = h_errors.len();
    if n == 0 { eprintln!("No results!"); return; }
    h_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    println!("=== Odaiba SWFG RTK Accuracy ===");
    println!("Epochs: {n}");
    println!("Horizontal error:");
    println!("  p50:  {:.3}m", h_errors[((n as f64) * 0.50) as usize]);
    println!("  p68:  {:.3}m", h_errors[((n as f64) * 0.68) as usize]);
    println!("  p95:  {:.3}m", h_errors[((n as f64) * 0.95) as usize]);
    println!("  RMS:  {:.3}m", (h_errors.iter().map(|e| e*e).sum::<f64>() / n as f64).sqrt());
}
