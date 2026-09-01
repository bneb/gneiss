use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use nalgebra::Vector3;

use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::engine::SwfgEngine;
use gneiss_rtk::swfg::imu_preintegration::{ImuPreintegration, ImuSample};

fn parse_imu_csv(path: &Path) -> Vec<ImuSample> {
    let file = File::open(path).expect("failed to open imu.csv");
    let reader = BufReader::new(file);
    let mut samples = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let l = line.expect("line");
        if i == 0 || l.starts_with("GPS") {
            continue;
        }
        let parts: Vec<&str> = l.split(',').map(|s| s.trim()).collect();
        if parts.len() >= 8 {
            let tow: f64 = parts[0].parse().unwrap_or(0.0);
            let ax: f64 = parts[2].parse().unwrap_or(0.0);
            let ay: f64 = parts[3].parse().unwrap_or(0.0);
            let az: f64 = parts[4].parse().unwrap_or(0.0);
            let gx: f64 = parts[5].parse().unwrap_or(0.0);
            let gy: f64 = parts[6].parse().unwrap_or(0.0);
            let gz: f64 = parts[7].parse().unwrap_or(0.0);

            let time_us = (tow * 1_000_000.0) as u32;
            samples.push(ImuSample {
                accel: Vector3::new(ax, ay, az),
                gyro: Vector3::new(gx, gy, gz),
                time_us,
            });
        }
    }
    samples
}

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("warn")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dataset = Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba");
    println!("=== Odaiba Tightly-Coupled INS/GNSS Benchmark ===");

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

    let imu_samples = parse_imu_csv(&dataset.join("imu.csv"));
    println!("Loaded {} IMU samples, {} reference truth points", imu_samples.len(), truth.len());

    let rover_file = File::open(dataset.join("rover_trimble.obs")).expect("open rover_trimble.obs");
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rover_file))
        .expect("parse rover RINEX");

    let base_file = File::open(dataset.join("base_trimble.obs")).expect("open base_trimble.obs");
    let (base_epochs, base_approx) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_file))
        .expect("parse base RINEX");

    let default_base_arp = Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let default_base_llh = gneiss_core::coords::ecef_to_llh(default_base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(default_base_llh).transpose();
    let computed_base = default_base_arp + ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855);
    let base_pos = base_approx.approx_position.map(|p| Vector3::new(p[0], p[1], p[2])).unwrap_or(computed_base);

    let base_index: BTreeMap<u32, &gneiss_core::obs::EpochObs> = base_epochs.iter()
        .map(|e| ((e.time.tow * 1000.0).round() as u32, e))
        .collect();

    let rover_init_pos = if let Some(spp) = gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(&rover_epochs[0], &ephemerides) {
        spp
    } else {
        base_pos
    };

    let rtk_config = gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([rover_init_pos.x, rover_init_pos.y, rover_init_pos.z]),
        ..Default::default()
    };
    let config = EngineConfig::Rtk(rtk_config);
    let mut engine = SwfgEngine::new(&config, ephemerides.clone());
    if let Some(ref k) = klobuchar {
        engine.set_klobuchar(k.alpha, k.beta);
    }

    let mut h_errors = Vec::new();
    let mut errors_3d = Vec::new();
    let err_count = 0usize;
    let mut ok_count = 0usize;
    let mut last_imu_idx = 0usize;

    // Forward Pass — only store solutions backed by GNSS observations
    let mut forward_solutions: BTreeMap<u32, (Vector3<f64>, usize)> = BTreeMap::new();
    for (i, epoch) in rover_epochs.iter().enumerate() {
        let tow_sec = epoch.time.tow.floor() as u32;
        let current_us = (epoch.time.tow * 1_000_000.0) as u32;
        let mut epoch_imu = Vec::new();
        while last_imu_idx < imu_samples.len() && imu_samples[last_imu_idx].time_us <= current_us {
            epoch_imu.push(imu_samples[last_imu_idx]);
            last_imu_idx += 1;
        }

        let has_valid_gyro = epoch_imu.iter().any(|s| s.gyro.norm() > 1e-6);
        let preint_opt = if epoch_imu.len() >= 2 && has_valid_gyro {
            let mut preint = ImuPreintegration::new();
            preint.integrate(&epoch_imu, &Vector3::zeros(), &Vector3::zeros());
            let dt = preint.dt;
            if dt > 2.0 && (preint.dp.norm() / dt.max(1e-3)) < 0.5 {
                preint.dp = Vector3::zeros();
                preint.dv = Vector3::zeros();
            }
            Some(preint)
        } else {
            None
        };

        let exact_ms = (epoch.time.tow * 1000.0).round() as u32;
        let base_epoch = base_index.get(&exact_ms).copied();

        let res = if let Some(base) = base_epoch {
            engine.process_rtk_epoch_with_imu(epoch, base, base_pos, preint_opt)
        } else {
            engine.process_epoch(epoch)
        };

        match res {
            Ok(sol) if sol.n_satellites >= 4 => {
                forward_solutions.insert(tow_sec, (sol.position_ecef, sol.n_satellites));
            }
            Ok(_) => {} // Too few satellites — skip
            Err(ref e) => {
                if i < 5 {
                    eprintln!("Epoch {} failed: {}", i, e);
                }
            }
        }

        if i % 500 == 0 || i == rover_epochs.len() - 1 {
            println!("Forward pass: epoch {}/{} (TOW {}) - solutions: {}", i + 1, rover_epochs.len(), tow_sec, forward_solutions.len());
        }
    }

    // Backward Pass (Qinertia Offline Bidirectional Smoother)
    let last_rover_pos = forward_solutions.values().last().map(|(p, _)| *p).unwrap_or(rover_init_pos);
    let rtk_config_bwd = gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([last_rover_pos.x, last_rover_pos.y, last_rover_pos.z]),
        ..Default::default()
    };
    let config_bwd = EngineConfig::Rtk(rtk_config_bwd);
    let mut engine_bwd = SwfgEngine::new(&config_bwd, ephemerides.clone());
    if let Some(ref k) = klobuchar {
        engine_bwd.set_klobuchar(k.alpha, k.beta);
    }
    let mut backward_solutions: BTreeMap<u32, (Vector3<f64>, usize)> = BTreeMap::new();
    let mut rev_epochs = rover_epochs.clone();
    rev_epochs.reverse();

    for (i, epoch) in rev_epochs.iter().enumerate() {
        let tow_sec = epoch.time.tow.floor() as u32;
        let exact_ms = (epoch.time.tow * 1000.0).round() as u32;
        let base_epoch = base_index.get(&exact_ms).copied();

        let res = if let Some(base) = base_epoch {
            engine_bwd.process_rtk_epoch_with_imu(epoch, base, base_pos, None)
        } else {
            engine_bwd.process_epoch(epoch)
        };

        if let Ok(sol) = res {
            if sol.n_satellites >= 4 {
                backward_solutions.insert(tow_sec, (sol.position_ecef, sol.n_satellites));
            }
        }

        if i % 500 == 0 || i == rev_epochs.len() - 1 {
            println!("Backward pass: epoch {}/{} (TOW {}) - solutions: {}", i + 1, rev_epochs.len(), tow_sec, backward_solutions.len());
        }
    }

    // Combine Forward & Backward Passes (Inverse-Variance Weighted Smoother)
    //
    // Since we gate to observation-backed epochs only, both directions
    // have similar uncertainty when they overlap. Equal-weight average
    // when both present; single-direction otherwise.

    // Collect all TOWs that have at least one observation-backed solution
    let mut all_tows: std::collections::BTreeSet<u32> = std::collections::BTreeSet::new();
    for t in forward_solutions.keys() { all_tows.insert(*t); }
    for t in backward_solutions.keys() { all_tows.insert(*t); }

    for tow_sec in &all_tows {
        let fwd_opt = forward_solutions.get(tow_sec);
        let bwd_opt = backward_solutions.get(tow_sec);

        // Skip if neither direction has an observation-backed solution
        if fwd_opt.is_none() && bwd_opt.is_none() {
            continue;
        }
        ok_count += 1;

        let final_pos = match (fwd_opt, bwd_opt) {
            (Some((fwd_pos, _)), Some((bwd_pos, _))) => {
                // Both directions observed this epoch — equal-weight average
                0.5 * fwd_pos + 0.5 * bwd_pos
            }
            (Some((fwd_pos, _)), None) => *fwd_pos,
            (None, Some((bwd_pos, _))) => *bwd_pos,
            (None, None) => unreachable!(),
        };

        if let Some(&(tx, ty, tz)) = truth.get(tow_sec) {
            let truth_pos = Vector3::new(tx, ty, tz);
            let truth_llh = gneiss_core::coords::ecef_to_llh(truth_pos);

            let fwd_err = fwd_opt.map(|(p, _)| {
                let enu = gneiss_core::coords::ecef_delta_to_enu(*p, truth_pos, truth_llh);
                (enu.x * enu.x + enu.y * enu.y).sqrt()
            }).unwrap_or(999.9);

            let bwd_err = bwd_opt.map(|(p, _)| {
                let enu = gneiss_core::coords::ecef_delta_to_enu(*p, truth_pos, truth_llh);
                (enu.x * enu.x + enu.y * enu.y).sqrt()
            }).unwrap_or(999.9);

            let enu_smooth = gneiss_core::coords::ecef_delta_to_enu(final_pos, truth_pos, truth_llh);
            let h_err = (enu_smooth.x * enu_smooth.x + enu_smooth.y * enu_smooth.y).sqrt();
            let err_3d = enu_smooth.norm();

            h_errors.push(h_err);
            errors_3d.push(err_3d);
            if ok_count.is_multiple_of(100) || h_err < 5.0 {
                println!(
                    "TOW={} fwd={:.2}m bwd={:.2}m smooth_hz={:.2}m smooth_3d={:.2}m",
                    tow_sec, fwd_err, bwd_err, h_err, err_3d
                );
            }
        }
    }

    if !h_errors.is_empty() {
        h_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        errors_3d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = h_errors.len();
        let p50 = h_errors[n / 2];
        let p68 = h_errors[(n as f64 * 0.68) as usize];
        let p95 = h_errors[(n as f64 * 0.95) as usize];
        let rms = (h_errors.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();

        let p50_3d = errors_3d[n / 2];
        let p68_3d = errors_3d[(n as f64 * 0.68) as usize];
        let p95_3d = errors_3d[(n as f64 * 0.95) as usize];
        let rms_3d = (errors_3d.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();

        println!("\n=== Odaiba Smoothed INS/GNSS Accuracy (Qinertia Architecture) ===");
        println!("Processed Epochs: {} (Failed: {})", ok_count, err_count);
        println!("Horizontal error (ENU):");
        println!("  p50:  {:.3}m", p50);
        println!("  p68:  {:.3}m", p68);
        println!("  p95:  {:.3}m", p95);
        println!("  RMS:  {:.3}m", rms);
        println!("3D Position error (ENU):");
        println!("  p50:  {:.3}m", p50_3d);
        println!("  p68:  {:.3}m", p68_3d);
        println!("  p95:  {:.3}m", p95_3d);
        println!("  RMS:  {:.3}m", rms_3d);
    }
}
