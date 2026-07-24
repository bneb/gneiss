//! Zero-baseline RTK test: TMG2 as both base and rover.
//!
//! Uses TMG2 CORS 1-second data from July 2, 2026.
//! Base: 08:00-09:00 UTC. Rover: 09:00-10:00 UTC.
//! Both at same fixed location — expected error = 0m.
//! Ground truth: NGS published TMG2 coordinate.
//!
//! Usage: cargo run --release --bin eval_tmg2_zb

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use gneiss_core::obs::EpochObs;
use gneiss_rtk::engine::{EngineConfig, EngineMode, ProcessingEngine};

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("info")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dir = Path::new("datasets/tmg2_test");

    // --- Read full daily RINEX, filter to two hours ---
    let rinex_path = dir.join("tmg2_day183.26o");
    let rinex_file = File::open(&rinex_path).expect("open RINEX");
    let rinex_reader = BufReader::new(rinex_file);
    let (all_epochs, header) = gneiss_parsers::rinex::parse_rinex_obs(rinex_reader)
        .expect("parse RINEX");
    eprintln!("Parsed {} total epochs", all_epochs.len());

    // Hour 08:00-09:00 = TOW 72000-75600 on day 3 (July 2, 2026 = Thursday = day 4)
    // Actually: July 2 2026 = Thursday. GPS day 4. TOW = 4*86400 + UTC_sod
    // 08:00 = 4*86400 + 28800 = 345600 + 28800 = 374400
    // 09:00 = 4*86400 + 32400 = 345600 + 32400 = 378000
    // Use overlapping windows: base 08:00-09:00, rover 08:30-09:00
    // The last base epoch at 09:00 covers rover epochs through 09:00
    let base_start_tow = 374400.0; // 08:00
    let base_end_tow = 378000.0;   // 09:00
    let rover_start_tow = 376200.0; // 08:30 (overlaps with base)
    let rover_end_tow = 378000.0;   // 09:00

    let base_epochs: Vec<&EpochObs> = all_epochs.iter()
        .filter(|e| e.time.tow >= base_start_tow && e.time.tow < base_end_tow)
        .collect();
    let rover_epochs: Vec<&EpochObs> = all_epochs.iter()
        .filter(|e| e.time.tow >= rover_start_tow && e.time.tow < rover_end_tow)
        .collect();
    eprintln!("Base epochs: {} ({}:00-{}:00)", base_epochs.len(), 8, 9);
    eprintln!("Rover epochs: {} ({}:00-{}:00)", rover_epochs.len(), 9, 10);

    if base_epochs.is_empty() || rover_epochs.is_empty() {
        eprintln!("No epochs in selected time range!");
        return;
    }

    // --- Read nav file ---
    let nav_path = dir.join("brdc1830.26n");
    let ephemerides: Vec<gneiss_core::ephemeris::Ephemeris>;
    let klobuchar;
    if nav_path.exists() {
        let nav_file = File::open(&nav_path).expect("open nav");
        let nav_reader = BufReader::new(nav_file);
        let (ephs, klob) = gneiss_parsers::rinex::parse_rinex_nav(nav_reader).expect("parse nav");
        ephemerides = ephs;
        klobuchar = klob;
        eprintln!("Loaded {} ephemerides, klobuchar={}", ephemerides.len(), klobuchar.is_some());
    } else {
        eprintln!("BRDC nav not found at {} — download from CDDIS FTP", nav_path.display());
        eprintln!("  ftp://gdc.cddis.eosdis.nasa.gov/gnss/data/daily/2026/183/25n/brdc1830.26n.gz");
        return;
    }

    // --- TMG2 published coordinate (NGS, IGS14 frame) ---
    let tmg2_xyz = [-1283433.934, -4713073.304, 4090105.087];

    // --- Index base epochs by TOW ---
    let base_index: std::collections::BTreeMap<u32, &EpochObs> = base_epochs
        .iter()
        .map(|e| { let k = (e.time.tow * 1000.0) as u32; (k, *e) })
        .collect();
    // Debug keys
    let mut keys: Vec<u32> = base_index.keys().copied().collect();
    keys.sort();
    eprintln!("Base keys first 3: {:?}, last 3: {:?}", &keys[..3.min(keys.len())], &keys[keys.len().saturating_sub(3)..]);
    eprintln!("First rover TOW: {:.6}", rover_epochs[0].time.tow);
    eprintln!("Rover key range: {}..={}",
        ((rover_epochs[0].time.tow - 2.0) * 1000.0) as u32,
        ((rover_epochs[0].time.tow + 2.0) * 1000.0) as u32);

    let config = EngineConfig {
        mode: EngineMode::Rtk,
        base_position: Some(tmg2_xyz),
        initial_position: Some(tmg2_xyz),
        elevation_mask_deg: 15.0,
        min_snr_dbhz: 25.0,
        chi_square_pr_threshold: 3.0,
        chi_square_cp_threshold: 3.0,
        dynamics_model: gneiss_rtk::engine::DynamicsModel::Static,
        enable_ar: true,
        enable_tdcp: false,
        lambda_min_ratio: 1.6,
        lambda_min_subset: 4,
        ar_min_epoch_count: 30,
        ar_min_lock: 5,
        ar_ffrt_prob: 0.001,
        pr_window_size: 0,
        initial_ambiguity_variance: 4.0,
        max_base_age_s: 2.0,
        enable_pr_validation: false,
        enable_ins_validation: false,
        ..Default::default()
    };

    let mut engine = ProcessingEngine::new(config);
    engine.ephemerides = ephemerides;
    engine.klobuchar_params = klobuchar;

    // --- Process rover epochs ---
    let mut results: Vec<(f64, f64, f64, bool)> = Vec::new(); // (tow, h_err, v_err, fixed)
    let mut num_fixed = 0usize;
    let mut num_processed = 0usize;

    for rover_epoch in &rover_epochs {
        let tow_ms = (rover_epoch.time.tow * 1000.0) as u32;
        let base_epoch = base_index
            .range(tow_ms.saturating_sub(2000)..=tow_ms.saturating_add(2000))
            .next()
            .map(|(_, e)| *e);
        if base_epoch.is_none() { continue; }

        match engine.process_rtk(rover_epoch, base_epoch) {
            Ok(state) => {
                num_processed += 1;
                if state.is_fixed { num_fixed += 1; }
                let dx = state.position.vector.x - tmg2_xyz[0];
                let dy = state.position.vector.y - tmg2_xyz[1];
                let dz = state.position.vector.z - tmg2_xyz[2];
                let h_err = (dx*dx + dy*dy).sqrt();
                let v_err = dz.abs();
                results.push((rover_epoch.time.tow, h_err, v_err, state.is_fixed));
                if num_processed <= 10 || state.is_fixed {
                    eprintln!("  tow={:.1} h_err={:.3}m v_err={:.3}m fixed={}",
                        rover_epoch.time.tow, h_err, v_err, state.is_fixed);
                }
            }
            Err(_) => {}
        }
    }

    let n = results.len();
    if n < 5 { eprintln!("Too few results ({n})"); return; }

    let mut h_errs: Vec<f64> = results.iter().map(|r| r.1).collect();
    h_errs.sort_by(|a,b| a.partial_cmp(b).unwrap());
    let mut v_errs: Vec<f64> = results.iter().map(|r| r.2).collect();
    v_errs.sort_by(|a,b| a.partial_cmp(b).unwrap());

    let fix_rate = num_fixed as f64 / n as f64 * 100.0;
    println!("=== TMG2 Zero-Baseline RTK (1-second data, July 2026) ===");
    println!("Epochs: {n} processed");
    println!("Fix rate: {fix_rate:.1}% ({num_fixed}/{n})");
    println!("Horizontal error vs NGS published coordinate:");
    println!("  p50: {:.4}m  p95: {:.4}m  RMS: {:.4}m",
        h_errs[((n as f64)*0.50) as usize],
        h_errs[((n as f64)*0.95) as usize],
        (h_errs.iter().map(|e| e*e).sum::<f64>()/n as f64).sqrt());
    println!("Vertical error:");
    println!("  p50: {:.4}m  p95: {:.4}m",
        v_errs[((n as f64)*0.50) as usize],
        v_errs[((n as f64)*0.95) as usize]);
}
