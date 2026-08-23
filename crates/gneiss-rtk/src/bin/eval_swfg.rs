//! SWFG engine evaluation on dataset trajectories with ground-truth accuracy metrics.

#![allow(clippy::unwrap_used)]

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

use nalgebra::Vector3;

use gneiss_core::coords::{ecef_to_llh, llh_to_ecef};
use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::engine::SwfgEngine;
use gneiss_rtk::swfg::smoothing::{BatchFactorGraph, BatchSmoothingConfig};

/// Helper to compute ENU offset from target ECEF relative to reference ECEF/LLH.
fn ecef_to_enu(target_ecef: Vector3<f64>, ref_ecef: Vector3<f64>, ref_llh: Vector3<f64>) -> Vector3<f64> {
    let lat = ref_llh.x;
    let lon = ref_llh.y;
    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_lon = lon.sin();
    let cos_lon = lon.cos();

    let dx = target_ecef.x - ref_ecef.x;
    let dy = target_ecef.y - ref_ecef.y;
    let dz = target_ecef.z - ref_ecef.z;

    let e = -sin_lon * dx + cos_lon * dy;
    let n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
    let u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;
    Vector3::new(e, n, u)
}

/// Simple parser for RTKLIB style .pos ground truth files.
fn load_ground_truth_pos(path: &Path) -> std::collections::BTreeMap<u32, Vector3<f64>> {
    let mut truth = std::collections::BTreeMap::new();
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => return truth,
    };
    let reader = BufReader::new(file);

    for line in reader.lines().map_while(Result::ok) {
        if line.starts_with('%') || line.trim().is_empty() || line.starts_with("GPS") {
            continue;
        }
        if line.contains(',') {
            let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
            if parts.len() >= 5 {
                if let (Ok(tow_f), Ok(lat_deg), Ok(lon_deg), Ok(h_m)) = (
                    parts[0].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                    parts[4].parse::<f64>(),
                ) {
                    let llh = Vector3::new(lat_deg.to_radians(), lon_deg.to_radians(), h_m);
                    let ecef = llh_to_ecef(llh);
                    truth.insert(tow_f.floor() as u32, ecef);
                }
            }
        } else {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 5 {
                if let (Ok(lat_deg), Ok(lon_deg), Ok(h_m)) = (
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                    parts[4].parse::<f64>(),
                ) {
                    let llh = Vector3::new(lat_deg.to_radians(), lon_deg.to_radians(), h_m);
                    let ecef = llh_to_ecef(llh);
                    if let Ok(sec) = parts[1].parse::<f64>() {
                        let tow = sec.floor() as u32;
                        truth.insert(tow, ecef);
                    }
                }
            }
        }
    }
    truth
}

fn percentile(mut data: Vec<f64>, p: f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = ((data.len() - 1) as f64 * p).round() as usize;
    data[idx]
}

fn evaluate_dataset(name: &str, nav_path: &Path, obs_path: &Path, base_obs_path: Option<&Path>, gt_path: Option<&Path>, _initial_ecef: [f64; 3]) {
    eprintln!("\n==========================================");
    eprintln!(" Evaluating Dataset: {}", name);
    eprintln!("==========================================");

    let nav_file = match File::open(nav_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("  Skipping {}: could not open nav file: {}", name, e);
            return;
        }
    };
    let (ephemerides, klobuchar) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_file))
        .expect("parse nav file");
    eprintln!("Loaded {} ephemerides, klobuchar={:?}", ephemerides.len(), klobuchar.is_some());

    let obs_file = match File::open(obs_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("  Skipping {}: could not open obs file: {}", name, e);
            return;
        }
    };
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(obs_file))
        .expect("parse obs file");
    eprintln!("Loaded {} observation epochs", rover_epochs.len());

    let mut base_epochs_map = std::collections::HashMap::new();
    let mut base_position_ecef = Vector3::zeros();
    
    if let Some(base_path) = base_obs_path {
        let base_file = match File::open(base_path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("  Skipping {}: could not open base obs file: {}", name, e);
                return;
            }
        };
        let (b_epochs, b_header) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_file))
            .expect("parse base obs file");
        if let Some(approx) = b_header.approx_position {
            base_position_ecef = Vector3::new(approx[0], approx[1], approx[2]);
            eprintln!("Parsed base position from header: {:?}", base_position_ecef);
        } else {
            eprintln!("WARNING: No APPROX POSITION XYZ found in base header!");
        }
        for be in b_epochs {
            let tow = be.time.tow.floor() as u32;
            base_epochs_map.insert(tow, be);
        }
        eprintln!("Loaded {} base observation epochs", base_epochs_map.len());
    }

    let gt = gt_path.map(load_ground_truth_pos).unwrap_or_default();
    if !gt.is_empty() {
        eprintln!("Loaded {} ground truth positions", gt.len());
    }

    let config = EngineConfig::Spp(gneiss_rtk::swfg::config::SppConfig::default());
    let mut engine = SwfgEngine::new(&config, ephemerides.clone());
    if let Some(ref k) = klobuchar {
        engine.set_klobuchar(k.alpha, k.beta);
    }

    let start = Instant::now();
    let mut n_processed = 0usize;
    let mut horizontal_errors = Vec::new();
    let mut errors_3d = Vec::new();

    for (i, epoch) in rover_epochs.iter().enumerate() {
        let tow = epoch.time.tow.floor() as u32;
        if i == 0 {
            let mut test_map = std::collections::HashMap::new();
            for e in &ephemerides { test_map.insert(e.sat(), e.clone()); }
            
            let ephs: Vec<_> = test_map.values().cloned().collect();
            let config = gneiss_rtk::estimators::spp::SppConfig::default();
            let meas = gneiss_rtk::estimators::spp::build_measurements(epoch, &ephs, &config);
            eprintln!("eval_swfg DIRECT SPP: epoch has {} measurements built", meas.len());

            let spp_res = gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(epoch, &ephs);
            eprintln!("eval_swfg DIRECT SPP TEST for epoch 0: {:?}", spp_res);
        }
        let res = engine.process_epoch(epoch);

        match res {
            Ok(sol) => {
                n_processed += 1;
                
                if let Some(&ref_ecef) = gt.get(&tow) {
                    let ref_llh = ecef_to_llh(ref_ecef);
                    let enu = ecef_to_enu(sol.position_ecef, ref_ecef, ref_llh);
                    let horiz_err = (enu.x * enu.x + enu.y * enu.y).sqrt();
                    let err_3d = enu.norm();
                    horizontal_errors.push(horiz_err);
                    errors_3d.push(err_3d);
                }

                if i < 2 || i % 50 == 0 {
                    let llh = ecef_to_llh(sol.position_ecef);
                    eprintln!(
                        "epoch {:4} | {:.6}°,{:.6}°,{:.1}m | {} sats",
                        i,
                        llh.x.to_degrees(),
                        llh.y.to_degrees(),
                        llh.z,
                        sol.n_satellites,
                    );
                }
            }
            Err(e) => {
                eprintln!("epoch {} failed: {}", i, e);
            }
        }
        if i >= 200 { break; }
    }

    let elapsed = start.elapsed();
    eprintln!("--- Online SWFG Results ---");
    eprintln!(
        "Processed {} epochs in {:.2}s ({:.2} ms/epoch)",
        n_processed,
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1000.0 / n_processed.max(1) as f64
    );

    if !horizontal_errors.is_empty() {
        let p50_h = percentile(horizontal_errors.clone(), 0.50);
        let p95_h = percentile(horizontal_errors.clone(), 0.95);
        let p50_3d = percentile(errors_3d.clone(), 0.50);
        let p95_3d = percentile(errors_3d.clone(), 0.95);

        eprintln!("Online SWFG Horizontal Error: 50th = {:.3} m | 95th = {:.3} m", p50_h, p95_h);
        eprintln!("Online SWFG 3D Position Error: 50th = {:.3} m | 95th = {:.3} m", p50_3d, p95_3d);
    }

    eprintln!("\nRunning Batch Factor Graph Trajectory Smoothing...");
    let batch_config = BatchSmoothingConfig {
        max_iterations: 25,
        convergence_tol: 1e-4,
        enable_ar: true,
        ar_ratio_threshold: 2.0,
    };
    let mut batch = BatchFactorGraph::new(batch_config, ephemerides);
    let mut epoch_idx = 0u32;
    for epoch in &rover_epochs {
        let tow = epoch.time.tow.floor() as u32;
        let base_opt = base_epochs_map.get(&tow).map(|be| (be, base_position_ecef));
        batch.add_epoch(epoch_idx, epoch, base_opt, None);
        epoch_idx += 1;
        if epoch_idx >= 101 { break; }
    }
    batch.add_smoothness_factors(30.0, 1.0);

    match batch.solve() {
        Ok(smoothed_traj) => {
            let mut batch_h_errors = Vec::new();
            for (time, pos) in smoothed_traj {
                let tow = time.tow.floor() as u32;
                if let Some(&ref_ecef) = gt.get(&tow) {
                    let ref_llh = ecef_to_llh(ref_ecef);
                    let enu = ecef_to_enu(pos, ref_ecef, ref_llh);
                    let horiz_err = (enu.x * enu.x + enu.y * enu.y).sqrt();
                    batch_h_errors.push(horiz_err);
                }
            }
            if !batch_h_errors.is_empty() {
                let p50_batch = percentile(batch_h_errors.clone(), 0.50);
                let p95_batch = percentile(batch_h_errors.clone(), 0.95);
                eprintln!("Batch Smoothed Horizontal Error: 50th = {:.3} m | 95th = {:.3} m", p50_batch, p95_batch);
            }
        }
        Err(e) => {
            eprintln!("Batch Smoothing Failed: {:?}", e);
        }
    }
}

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("warn")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    eprintln!("=== Gneiss Multi-Dataset SWFG Benchmark Suite ===");

    // 1. UrbanNav Tokyo Odaiba
    evaluate_dataset(
        "UrbanNav Tokyo Odaiba",
        Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base.nav"),
        Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba/rover_trimble.obs"),
        Some(Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base_trimble.obs")),
        Some(Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba/reference.csv")),
        [-3955635.0, 3349605.0, 3696805.0],
    );
}
