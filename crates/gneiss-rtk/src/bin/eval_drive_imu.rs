//! RTK+INS accuracy evaluation on drive_imu dataset (Boulder, CO).
//!
//! Rover: u-blox GNSS (UBX) + external IMU (100Hz, 6-axis)
//! Base: TMG2 CORS (Septentrio PolaRx5, 30s interval, ~8.2km baseline)
//! Ground truth: RTKPOST SPAN solution
//!
//! Usage: cargo run --release --bin eval_drive_imu

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use gneiss_core::imu::ImuMeasurement;
use gneiss_rtk::engine::{EngineConfig, EngineMode, ProcessingEngine};

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("info")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dir = Path::new("datasets/drive_imu");

    // --- Read ground truth (RTKPOST SPAN .pos) ---
    let pos_content = std::fs::read_to_string(dir.join("gnss_1934_sf.pos")).expect("read pos");
    let mut truth: BTreeMap<u32, (f64, f64, f64)> = BTreeMap::new();
    for line in pos_content.lines() {
        if line.starts_with('%') || line.is_empty() { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 8 { continue; }
        let q: f64 = parts[5].parse().unwrap_or(0.0);
        if q < 1.0 { continue; }
        let time_str = parts[1];
        let time_parts: Vec<&str> = time_str.split(':').collect();
        if time_parts.len() != 3 { continue; }
        let h: f64 = time_parts[0].parse().unwrap_or(0.0);
        let m: f64 = time_parts[1].parse().unwrap_or(0.0);
        let s: f64 = time_parts[2].parse().unwrap_or(0.0);
        // July 8, 2025 = Tuesday = day 2 of GPS week (Sunday=0)
        let tow = 2.0 * 86400.0 + h * 3600.0 + m * 60.0 + s;
        // This .pos format has LLH, not ECEF
        let lat: f64 = parts[2].parse().unwrap();
        let lon: f64 = parts[3].parse().unwrap();
        let hgt: f64 = parts[4].parse().unwrap();
        // Convert LLH to ECEF
        let (x, y, z) = llh_to_ecef(lat, lon, hgt);
        let tow_ms = (tow * 1000.0) as u32;
        truth.insert(tow_ms, (x, y, z));
    }
    eprintln!("Loaded {} ground truth epochs (fixed only)", truth.len());

    // --- Parse UBX rover data (RAWX + SFRBX in single pass) ---
    let ubx_data = std::fs::read(dir.join("gnss_1934.ubx")).expect("read UBX");
    let (rover_epochs, _sfrbx) = parse_ubx_combined(&ubx_data);
    eprintln!("Parsed {} rover GNSS epochs", rover_epochs.len());

    // --- Read navigation (BRDC broadcast ephemeris from CDDIS FTP) ---
    let nav_path = dir.join("brdc1890.25n");
    let ephemerides: Vec<gneiss_core::ephemeris::Ephemeris> = if nav_path.exists() {
        let nav_file = File::open(&nav_path).expect("open nav");
        let nav_reader = BufReader::new(nav_file);
        let (ephs, _) = gneiss_parsers::rinex::parse_rinex_nav(nav_reader).expect("parse nav");
        eprintln!("Loaded {} ephemerides from {}", ephs.len(), nav_path.display());
        ephs
    } else {
        eprintln!("BRDC nav not found — positioning will fail");
        Vec::new()
    };

    // --- Read base station RINEX (TMG2 CORS, 30s interval) ---
    let base_path = dir.join("tmg21890.25o");
    let base_epochs: Vec<gneiss_core::obs::EpochObs> = if base_path.exists() {
        let base_file = File::open(&base_path).expect("open base RINEX");
        let base_reader = BufReader::new(base_file);
        let (epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(base_reader)
            .expect("parse base RINEX");
        eprintln!("Parsed {} base epochs (TMG2 CORS, 30s interval)", epochs.len());
        epochs
    } else {
        eprintln!("Base RINEX not found at {}", base_path.display());
        return;
    };

    // --- Read IMU data ---
    // IMU time tags are in microseconds, aligned to GNSS time.
    // The first GNSS epoch gives the base TOW; IMU tags are offsets from that.
    let mut imu_measurements: Vec<ImuMeasurement> = Vec::new();
    let imu_csv = std::fs::read_to_string(dir.join("imu_1934.csv")).expect("read IMU CSV");
    for line in imu_csv.lines() {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 7 { continue; }
        let ax: f64 = parts[0].trim().parse().unwrap_or(0.0);
        let ay: f64 = parts[1].trim().parse().unwrap_or(0.0);
        let az: f64 = parts[2].trim().parse().unwrap_or(0.0);
        let gx: f64 = parts[3].trim().parse().unwrap_or(0.0);
        let gy: f64 = parts[4].trim().parse().unwrap_or(0.0);
        let gz: f64 = parts[5].trim().parse().unwrap_or(0.0);
        let time_tag: u32 = parts[6].trim().parse().unwrap_or(0);
        imu_measurements.push(ImuMeasurement {
            accel: nalgebra::Vector3::new(ax, ay, az),
            gyro: nalgebra::Vector3::new(gx, gy, gz),
            time_tag,
            temperature: None,
        });
    }
    eprintln!("Loaded {} IMU measurements", imu_measurements.len());
    // Align IMU time tags to GPS TOW: first GNSS TOW corresponds to first IMU time tag
    let imu_base_tag = imu_measurements.first().map(|m| m.time_tag).unwrap_or(0);
    let gnss_base_tow = rover_epochs.first().map(|e| e.time.tow).unwrap_or(0.0);
    eprintln!("IMU base tag={imu_base_tag}, GNSS base TOW={gnss_base_tow:.3}");

    // --- Index base epochs ---
    let base_index: BTreeMap<u32, &gneiss_core::obs::EpochObs> = base_epochs
        .iter()
        .map(|e| { let tow_ms = (e.time.tow * 1000.0) as u32; (tow_ms, e) })
        .collect();

    // --- Use first truth position as initial position ---
    let initial_pos = truth.values().next().copied()
        .unwrap_or((-1277000.1, -4717237.1, 4087230.1));
    eprintln!("Initial position: {initial_pos:?}");

    let config = EngineConfig {
        mode: EngineMode::Rtk,  // RTK-only for comparison
        base_position: Some([-1283433.936, -4713073.299, 4090105.086]),
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
        max_base_age_s: 15.0, // TMG2 is 30s interval, allow matching
        enable_pr_validation: true,
        enable_ins_validation: false,
        enable_nhc: true,
        ..Default::default()
    };

    let mut engine = ProcessingEngine::new(config);
    engine.ephemerides = ephemerides;

    // --- Process ---
    let mut results: Vec<(f64, f64, f64, f64, f64, f64, f64, bool)> = Vec::new();
    let mut skipped_no_base = 0usize;
    let mut skipped_no_imu = 0usize;
    let mut first = true;
    let mut first_tow_debug = true;

    // Sort IMU measurements by time tag
    imu_measurements.sort_by_key(|m| m.time_tag);

    for rover_epoch in &rover_epochs {
        let tow_ms = (rover_epoch.time.tow * 1000.0) as u32;
        // Push all IMU measurements up to this GNSS epoch's time tag.
        // This is called before process_rtk, which clears the buffer via predict_state.
        // So we repopulate the buffer for each epoch. Slightly redundant but simple.
        let gnss_tag = ((rover_epoch.time.tow - gnss_base_tow) * 1_000_000.0) as u32 + imu_base_tag;
        let end_idx = imu_measurements.partition_point(|m| m.time_tag <= gnss_tag);
        for m in imu_measurements[..end_idx].iter().rev().take(50) {
            engine.add_imu_measurement(m.clone());
        }
        // Skip if still no IMU data
        if engine.imu_buffer.is_empty() { skipped_no_imu += 1; continue; }

        // Find nearest base epoch (within 15s for 30s CORS data)
        let base_epoch = base_index
            .range(tow_ms.saturating_sub(15000)..=tow_ms.saturating_add(15000))
            .next()
            .map(|(_, e)| *e);
        if base_epoch.is_none() {
            skipped_no_base += 1; continue;
        }
        if first {
            let base_age = (rover_epoch.time.tow - base_epoch.unwrap().time.tow).abs();
            eprintln!("  First base match: rover_tow={:.3} base_tow={:.3} age={:.1}s",
                rover_epoch.time.tow, base_epoch.unwrap().time.tow, base_age);
            first = false;
        }

        match engine.process_rtk(rover_epoch, base_epoch) {
            Ok(state) => {
                // Match by truncation — rover TOW has fractional part, truth is integer seconds
                // Match truth by nearest TOW (±1 second tolerance) — rover and
                // truth may have sub-second offsets that cross rounding boundaries
                let rover_tow = rover_epoch.time.tow;
                let mut best_tow_ms: Option<u32> = None;
                let mut best_dist = 9999u32;
                // Search ±1 second range
                for dt in 0..=1000u32 {
                    for &offset in &[dt, dt.wrapping_neg()] {
                        let key = (rover_tow as u32).wrapping_mul(1000).wrapping_add(offset);
                        if truth.contains_key(&key) { best_tow_ms = Some(key); best_dist = dt; break; }
                    }
                    if best_tow_ms.is_some() { break; }
                }
                if let Some(key) = best_tow_ms {
                    if let Some((tx, ty, tz)) = truth.get(&key).copied() {
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
                }  // best_tow_ms
            }
            Err(_e) => {}
        }
    }
    eprintln!("Processed {} epochs with results (skipped: {} no_base, {} no_imu)",
        results.len(), skipped_no_base, skipped_no_imu);

    let n = results.len();
    if n < 10 { eprintln!("Too few results ({n})"); return; }

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

    println!("=== Drive IMU RTK+INS Accuracy (8.2km to TMG2 CORS, Boulder CO) ===");
    println!("Epochs: {n}");
    println!("Fix rate: {fix_rate:.1}% ({fixed_count}/{n})");
    println!("Horizontal error:");
    println!("  p50:  {:.3}m", p50);
    println!("  p68:  {:.3}m", p68);
    println!("  p95:  {:.3}m", p95);
    println!("  RMS:  {:.3}m", rms);

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

fn llh_to_ecef(lat_deg: f64, lon_deg: f64, h_m: f64) -> (f64, f64, f64) {
    let a = 6378137.0;
    let f = 1.0 / 298.257223563;
    let e2 = 2.0 * f - f * f;
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let n = a / (1.0 - e2 * lat.sin().powi(2)).sqrt();
    let x = (n + h_m) * lat.cos() * lon.cos();
    let y = (n + h_m) * lat.cos() * lon.sin();
    let z = (n * (1.0 - e2) + h_m) * lat.sin();
    (x, y, z)
}

fn parse_ubx_combined(data: &[u8]) -> (Vec<gneiss_core::obs::EpochObs>, Vec<gneiss_parsers::ubx::UbxRxmSfrbx>) {
    let mut epochs: Vec<gneiss_core::obs::EpochObs> = Vec::new();
    let mut sfrbx_msgs: Vec<gneiss_parsers::ubx::UbxRxmSfrbx> = Vec::new();
    let mut pos = 0;
    while pos < data.len() {
        match gneiss_parsers::ubx::parse_ubx_frame(&data[pos..]) {
            Ok((rem, frame)) => {
                if frame.class == 0x02 {
                    if frame.id == 0x15 {
                        if let Ok(rawx) = gneiss_parsers::ubx::parse_rxm_rawx(frame.payload) {
                            epochs.push(rawx.into_epoch_obs());
                        }
                    } else if frame.id == 0x13 {
                        if let Ok(sfrbx) = gneiss_parsers::ubx::parse_rxm_sfrbx(frame.payload) {
                            sfrbx_msgs.push(sfrbx);
                        }
                    }
                }
                pos = data.len() - rem.len();
            }
            Err(_) => { pos += 1; }
        }
    }
    (epochs, sfrbx_msgs)
}
