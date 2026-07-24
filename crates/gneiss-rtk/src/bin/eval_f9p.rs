//! RTK accuracy evaluation on F9P driving dataset (~8.2km baseline to TMG2).
//!
//! Rover: u-blox F9P, Base: TMG2 CORS, Ground truth: RTKPOST PPK fixed solution.
//! Usage: cargo run --release --bin eval_f9p

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use gneiss_rtk::engine::{EngineConfig, EngineMode, ProcessingEngine};

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("info")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dir = Path::new("datasets/rtkexplorer/sample_1/f9p_ppp_1224");

    // --- Read navigation ---
    let nav_file = File::open(dir.join("rover.nav")).expect("open rover.nav");
    let nav_reader = BufReader::new(nav_file);
    let (ephemerides, klobuchar_params) =
        gneiss_parsers::rinex::parse_rinex_nav(nav_reader).expect("parse nav");
    eprintln!("Loaded {} ephemerides", ephemerides.len());

    // --- Read ground truth (RTKPOST PPK .pos) ---
    let pos_content = std::fs::read_to_string(dir.join("rover_ppk.pos")).expect("read pos");
    let mut truth: BTreeMap<u32, (f64, f64, f64)> = BTreeMap::new();
    for line in pos_content.lines() {
        if line.starts_with('%') || line.is_empty() { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 8 { continue; }
        // Format: YYYY/MM/DD HH:MM:SS.SSS X Y Z Q NSAT SDX SDY ...
        let q: i32 = parts[5].parse().unwrap_or(0);
        if q != 1 { continue; } // fixed only
        let time_str = parts[1];
        let time_parts: Vec<&str> = time_str.split(':').collect();
        if time_parts.len() != 3 { continue; }
        let h: f64 = time_parts[0].parse().unwrap_or(0.0);
        let m: f64 = time_parts[1].parse().unwrap_or(0.0);
        let s: f64 = time_parts[2].parse().unwrap_or(0.0);
        // Dec 24, 2020 = day 4 of GPS week (Thu). TOW = day4*86400 + sod
        let tow = 345600.0 + h * 3600.0 + m * 60.0 + s;
        let x: f64 = parts[2].parse().unwrap();
        let y: f64 = parts[3].parse().unwrap();
        let z: f64 = parts[4].parse().unwrap();
        let tow_ms = (tow * 1000.0) as u32;
        truth.insert(tow_ms, (x, y, z));
    }
    eprintln!("Loaded {} ground truth epochs (fixed only)", truth.len());

    // --- Read rover RINEX ---
    let rover_file = File::open(dir.join("rover.obs")).expect("open rover.obs");
    let rover_reader = BufReader::new(rover_file);
    let (rover_epochs, _) =
        gneiss_parsers::rinex::parse_rinex_obs(rover_reader).expect("parse rover RINEX");
    eprintln!("Parsed {} rover epochs", rover_epochs.len());

    // --- Read base RINEX ---
    let base_file = File::open(dir.join("tmg23590.20o")).expect("open base obs");
    let base_reader = BufReader::new(base_file);
    let (base_epochs, _) =
        gneiss_parsers::rinex::parse_rinex_obs(base_reader).expect("parse base RINEX");
    eprintln!("Truth time range: {:?} - {:?}",
        truth.keys().min(), truth.keys().max());
    eprintln!("Rover time range: tow={:.3} - {:.3}",
        rover_epochs.first().map(|e| e.time.tow).unwrap_or(0.0),
        rover_epochs.last().map(|e| e.time.tow).unwrap_or(0.0));
    eprintln!("Rover week: {:?}", rover_epochs.first().map(|e| e.time.week));

    // --- Index base epochs ---
    let base_index: BTreeMap<u32, &gneiss_core::obs::EpochObs> = base_epochs
        .iter()
        .map(|e| { let tow_ms = (e.time.tow * 1000.0) as u32; (tow_ms, e) })
        .collect();

    // --- Load SP3/CLK (optional, for precise orbits) ---
    let sp3_path = dir.join("COD0MGXFIN_20203590000_01D_05M_ORB.SP3");
    let sp3_epochs = if sp3_path.exists() {
        let sp3_file = File::open(&sp3_path).expect("open SP3");
        let sp3 = gneiss_parsers::sp3::parse_sp3(BufReader::new(sp3_file)).expect("parse SP3");
        eprintln!("Loaded {} SP3 epochs", sp3.len());
        sp3
    } else { Vec::new() };
    let clk_data = None; // skip CLK for speed (SP3 clock is sufficient)

    // --- Get initial position from first truth ---
    // Note: first truth is at TOW 422939 (21:28:59), rover starts at 422922 (21:28:42).
    // We start processing at the first truth epoch so initial position is correct.
    let first_truth_tow: f64 = truth.keys().min().map(|k| *k as f64 / 1000.0).unwrap_or(422939.0);
    let initial_pos = truth.values().next().copied().unwrap_or((-1276972.33, -4717195.55, 4087248.78));
    eprintln!("First truth TOW: {first_truth_tow:.1}, initial pos: {initial_pos:?}");

    let config = EngineConfig {
        mode: EngineMode::Rtk,
        base_position: Some([-1283434.625, -4713071.983, 4090105.048]),
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
        initial_ambiguity_variance: 4.0,
        max_base_age_s: 5.0,
        enable_pr_validation: true,
        enable_ins_validation: false,
        ..Default::default()
    };

    let mut engine = ProcessingEngine::new(config);
    engine.ephemerides = ephemerides;
    engine.klobuchar_params = klobuchar_params;
    engine.sp3_epochs = sp3_epochs;
    engine.clk_data = clk_data;

    // --- Process ---
    let mut results: Vec<(f64, f64, f64, f64, f64, f64, f64, bool)> = Vec::new();
    let mut err_count = 0u64;
    let mut ok_count = 0u64;

    for rover_epoch in &rover_epochs {
        // Skip epochs before first truth (EKF needs correct initial position)
        if rover_epoch.time.tow < first_truth_tow - 0.1 { continue; }
        // Stop after 300 epochs for open-sky segment analysis
        // if results.len() >= 300 { break; }
        let tow_ms = (rover_epoch.time.tow * 1000.0) as u32;
        let base_epoch = base_index
            .range(tow_ms.saturating_sub(5000)..=tow_ms.saturating_add(5000))
            .next()
            .map(|(_, e)| *e);
        if base_epoch.is_none() { continue; }

        match engine.process_rtk(rover_epoch, base_epoch) {
            Ok(state) => {
                ok_count += 1;
                // Truth TOW is integer seconds; rover TOW has fractional.
                // Round to nearest second for matching.
                let truth_tow_ms = ((rover_epoch.time.tow + 0.5) as u32) * 1000;
                if let Some((tx, ty, tz)) = truth.get(&truth_tow_ms).copied() {
                    let h_err = ((state.position.vector.x - tx).powi(2)
                        + (state.position.vector.y - ty).powi(2))
                    .sqrt();
                    results.push((
                        rover_epoch.time.tow,
                        state.position.vector.x, state.position.vector.y, state.position.vector.z,
                        tx, ty, tz,
                        state.is_fixed,
                    ));
                    if results.len() <= 10 {
                        eprintln!("tow={:.1} h_err={:.2}m fixed={}",
                            rover_epoch.time.tow, h_err, state.is_fixed);
                    }
                }
            }
            Err(_e) => { err_count += 1; }
        }
    }
    eprintln!("OK: {ok_count}, Errors: {err_count}");

    // --- Statistics ---
    let n = results.len();
    if n == 0 { eprintln!("No results!"); return; }

    let mut h_errors: Vec<f64> = results.iter().map(|(_, ex, ey, _, tx, ty, _, _)| {
        ((ex - tx).powi(2) + (ey - ty).powi(2)).sqrt()
    }).collect();
    h_errors.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let fixed_count = results.iter().filter(|r| r.7).count();
    let fix_rate = fixed_count as f64 / n as f64 * 100.0;
    let p50 = h_errors[((n as f64) * 0.50) as usize];
    let p68 = h_errors[((n as f64) * 0.68) as usize];
    let p95 = h_errors[((n as f64) * 0.95) as usize];
    let rms = (h_errors.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();

    println!("=== F9P RTK Accuracy (8.2km baseline to TMG2, California US) ===");
    println!("Epochs: {n}");
    println!("Fix rate: {fix_rate:.1}% ({fixed_count}/{n})");
    println!("Horizontal error (all epochs):");
    println!("  p50:  {:.3}m", p50);
    println!("  p68:  {:.3}m", p68);
    println!("  p95:  {:.3}m", p95);
    println!("  RMS:  {:.3}m", rms);

    // Compute speed-filtered statistics using truth positions
    let mut speed_errors: Vec<(f64, f64)> = Vec::new(); // (speed_m_s, h_err)
    for i in 1..results.len() {
        let (_, _, _, _, tx1, ty1, tz1, _) = results[i-1];
        let (_, _, _, _, tx2, ty2, tz2, _) = results[i];
        let dt = 1.0; // 1Hz
        let dist = ((tx2-tx1).powi(2) + (ty2-ty1).powi(2) + (tz2-tz1).powi(2)).sqrt();
        let speed = dist / dt;
        let (_, ex, ey, _, tx, ty, _, _) = results[i];
        let h_err = ((ex - tx).powi(2) + (ey - ty).powi(2)).sqrt();
        speed_errors.push((speed, h_err));
    }

    // Highway: >20 m/s (72 km/h), Open road: >5 m/s (18 km/h), Stopped: <1 m/s
    for (label, min_speed, max_speed) in &[
        ("Highway (>20 m/s)", 20.0, 999.0),
        ("Open road (>10 m/s)", 10.0, 999.0),
        ("Moving (>5 m/s)", 5.0, 999.0),
    ] {
        let mut errs: Vec<f64> = speed_errors.iter()
            .filter(|(s, _)| *s >= *min_speed && *s < *max_speed)
            .map(|(_, e)| *e)
            .collect();
        if errs.len() >= 10 {
            errs.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let m = errs.len();
            println!("{} (n={}):", label, m);
            println!("  p50: {:.3}m  p95: {:.3}m  RMS: {:.3}m",
                errs[((m as f64)*0.50) as usize],
                errs[((m as f64)*0.95) as usize],
                (errs.iter().map(|e| e*e).sum::<f64>() / m as f64).sqrt()
            );
        }
    }

    if fix_rate > 0.0 {
        let mut fixed_errors: Vec<f64> = results.iter()
            .filter(|r| r.7)
            .map(|(_, ex, ey, _, tx, ty, _, _)| {
                ((ex - tx).powi(2) + (ey - ty).powi(2)).sqrt()
            }).collect();
        fixed_errors.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let fn_ = fixed_errors.len();
        if fn_ > 0 {
            println!("Fixed-only (n={fn_}):");
            println!("  p50:  {:.3}m", fixed_errors[((fn_ as f64) * 0.50) as usize]);
            println!("  p95:  {:.3}m", fixed_errors[((fn_ as f64) * 0.95) as usize]);
        }
    }
}
