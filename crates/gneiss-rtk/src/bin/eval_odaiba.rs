//! RTK accuracy evaluation on Odaiba 4km baseline.
//!
//! Reads UBX rover + RINEX base + nav, processes through RTK engine,
//! compares against reference truth, reports p50/p95/fix rate.
//!
//! Usage:
//!   cargo run --release --bin eval_odaiba

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use gneiss_rtk::engine::{EngineConfig, EngineMode, ProcessingEngine};

fn main() {
    // Enable tracing for diagnostics
    tracing_subscriber::fmt()
        .with_env_filter("info")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dataset = Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba");

    // --- Read navigation ---
    let nav_file = File::open(dataset.join("base.nav")).expect("open base.nav");
    let nav_reader = BufReader::new(nav_file);
    let (ephemerides, klobuchar_params) = gneiss_parsers::rinex::parse_rinex_nav(nav_reader)
        .expect("parse nav");
    eprintln!("Loaded {} ephemerides", ephemerides.len());

    // --- Read reference truth ---
    let ref_csv =
        std::fs::read_to_string(dataset.join("reference.csv")).expect("read reference.csv");
    let mut truth: BTreeMap<u32, (f64, f64, f64)> = BTreeMap::new();
    for line in ref_csv.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 8 {
            continue;
        }
        let tow: f64 = parts[0].trim().parse().unwrap();
        let x: f64 = parts[5].trim().parse().unwrap();
        let y: f64 = parts[6].trim().parse().unwrap();
        let z: f64 = parts[7].trim().parse().unwrap();
        let tow_millis = (tow * 1000.0) as u32;
        truth.insert(tow_millis, (x, y, z));
    }
    eprintln!("Loaded {} reference epochs", truth.len());

    // --- Read rover RINEX ---
    let rover_file =
        File::open(dataset.join("rover_ublox.obs")).expect("open rover_ublox.obs");
    let rover_reader = BufReader::new(rover_file);
    let (rover_epochs, _) =
        gneiss_parsers::rinex::parse_rinex_obs(rover_reader).expect("parse rover RINEX");
    eprintln!("Parsed {} rover epochs", rover_epochs.len());

    // --- Read base RINEX ---
    let base_file =
        File::open(dataset.join("base_trimble.obs")).expect("open base_trimble.obs");
    let base_reader = BufReader::new(base_file);
    let (base_epochs, _) =
        gneiss_parsers::rinex::parse_rinex_obs(base_reader).expect("parse base RINEX");
    eprintln!("Parsed {} base epochs", base_epochs.len());

    // --- Index base epochs by TOW ---
    let base_index: BTreeMap<u32, &gneiss_core::obs::EpochObs> = base_epochs
        .iter()
        .map(|e| {
            let tow_ms = (e.time.tow * 1000.0) as u32;
            (tow_ms, e)
        })
        .collect();

    // --- Configure engine ---
    // Get initial position from reference (first epoch truth)
    let first_tow = rover_epochs[0].time.tow;
    let initial_pos = truth
        .get(&((first_tow * 1000.0) as u32))
        .copied()
        .unwrap_or((-3963426.8, 3350882.2, 3694865.5));

    let config = EngineConfig {
        mode: EngineMode::Rtk,
        base_position: Some([-3961904.4341, 3348994.266, 3698211.7067]),
        initial_position: Some([initial_pos.0, initial_pos.1, initial_pos.2]),
        elevation_mask_deg: 15.0,
        min_snr_dbhz: 25.0,
        chi_square_pr_threshold: 3.0,
        chi_square_cp_threshold: 3.0,
        dynamics_model: gneiss_rtk::engine::DynamicsModel::Automotive,
        enable_ar: true,
        enable_tdcp: false,
        lambda_min_ratio: 1.6,
        lambda_min_subset: 4,
        ar_min_epoch_count: 20,
        ar_min_lock: 3,
        ar_ffrt_prob: 0.001,
        pr_window_size: 200,
        process_noise_cb: 1.0,
        process_noise_cd: 10.0,
        initial_ambiguity_variance: 4.0,
        max_base_age_s: 5.0,
        enable_pr_validation: true,
        enable_ins_validation: false,
        ..Default::default()
    };
    eprintln!("Initial position from truth: {initial_pos:?}");

    let mut engine = ProcessingEngine::new(config);
    engine.ephemerides = ephemerides;
    engine.klobuchar_params = klobuchar_params;

    // --- Process ---
    let mut results: Vec<(f64, f64, f64, f64, f64, f64, f64, bool)> = Vec::new();
    let mut err_count = 0u64;
    let mut ok_count = 0u64;

    for rover_epoch in &rover_epochs {
        let tow_ms = (rover_epoch.time.tow * 1000.0) as u32;

        // Find nearest base epoch (within 5s to handle 1Hz base)
        let base_epoch = base_index
            .range(tow_ms.saturating_sub(5000)..=tow_ms.saturating_add(5000))
            .next()
            .map(|(_, e)| *e);

        if base_epoch.is_none() {
            continue;
        }

        match engine.process_rtk(rover_epoch, base_epoch) {
            Ok(state) => {
                ok_count += 1;
                if let Some((tx, ty, tz)) = truth.get(&tow_ms).copied() {
                    let is_fixed = state.is_fixed;
                    let h_err = ((state.position.vector.x - tx).powi(2)
                        + (state.position.vector.y - ty).powi(2))
                    .sqrt();
                    results.push((
                        rover_epoch.time.tow,
                        state.position.vector.x,
                        state.position.vector.y,
                        state.position.vector.z,
                        tx,
                        ty,
                        tz,
                        is_fixed,
                    ));
                    // Print first 10 results for debugging
                    if results.len() <= 10 {
                        eprintln!(
                            "tow={:.1} pos=({:.1},{:.1},{:.1}) true=({:.1},{:.1},{:.1}) h_err={:.2}m fixed={}",
                            rover_epoch.time.tow,
                            state.position.vector.x,
                            state.position.vector.y,
                            state.position.vector.z,
                            tx,
                            ty,
                            tz,
                            h_err,
                            is_fixed
                        );
                    }
                }
            }
            Err(e) => {
                err_count += 1;
                if err_count <= 5 {
                    eprintln!(
                        "Error at tow={:.1} (error #{err_count})",
                        rover_epoch.time.tow
                    );
                }
            }
        }
    }
    eprintln!("OK: {ok_count}, Errors: {err_count}");

    // --- Compute statistics ---
    let n = results.len();
    if n == 0 {
        eprintln!("No results!");
        return;
    }

    let mut h_errors: Vec<f64> = results
        .iter()
        .map(|(_tow, ex, ey, _ez, tx, ty, _tz, _fixed)| {
            let dx = ex - tx;
            let dy = ey - ty;
            (dx * dx + dy * dy).sqrt()
        })
        .collect();
    h_errors.sort_by(|a: &f64, b: &f64| a.partial_cmp(b).unwrap());

    let fixed_count = results.iter().filter(|r| r.7).count();
    let fix_rate = fixed_count as f64 / n as f64 * 100.0;

    let p50 = h_errors[((n as f64) * 0.50) as usize];
    let p68 = h_errors[((n as f64) * 0.68) as usize];
    let p95 = h_errors[((n as f64) * 0.95) as usize];
    let rms = (h_errors.iter().map(|e: &f64| e * e).sum::<f64>() / n as f64).sqrt();

    println!("=== Odaiba RTK Accuracy ===");
    println!("Epochs: {n}");
    println!("Fix rate: {fix_rate:.1}% ({fixed_count}/{n})");
    println!("Horizontal error:");
    println!("  p50:  {:.3}m", p50);
    println!("  p68:  {:.3}m", p68);
    println!("  p95:  {:.3}m", p95);
    println!("  RMS:  {:.3}m", rms);

    if fix_rate > 0.0 {
        let mut fixed_errors: Vec<f64> = results
            .iter()
            .filter(|r| r.7)
            .map(|(_tow, ex, ey, _ez, tx, ty, _tz, _fixed)| {
                let dx = ex - tx;
                let dy = ey - ty;
                (dx * dx + dy * dy).sqrt()
            })
            .collect();
        fixed_errors.sort_by(|a: &f64, b: &f64| a.partial_cmp(b).unwrap());
        let fn_ = fixed_errors.len();
        if fn_ > 0 {
            println!("Fixed-only (n={fn_}):");
            println!(
                "  p50:  {:.3}m",
                fixed_errors[((fn_ as f64) * 0.50) as usize]
            );
            println!(
                "  p95:  {:.3}m",
                fixed_errors[((fn_ as f64) * 0.95) as usize]
            );
        }
    }
}

/// Parse UBX binary file into EpochObs vector.
fn parse_ubx_epochs(data: &[u8]) -> Vec<gneiss_core::obs::EpochObs> {
    let mut epochs: Vec<gneiss_core::obs::EpochObs> = Vec::new();
    let mut pos = 0;
    while pos < data.len() {
        match gneiss_parsers::ubx::parse_ubx_frame(&data[pos..]) {
            Ok((rem, frame)) => {
                if frame.class == 0x02 && frame.id == 0x15 {
                    // UBX-RXM-RAWX
                    if let Ok(rawx) = gneiss_parsers::ubx::parse_rxm_rawx(frame.payload) {
                        epochs.push(rawx.into_epoch_obs());
                    }
                }
                pos = data.len() - rem.len();
            }
            Err(_) => {
                pos += 1;
            }
        }
    }
    epochs
}
