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

fn build_preintegration(accumulated_imu: &[ImuSample]) -> Option<ImuPreintegration> {
    let has_valid_gyro = accumulated_imu.iter().any(|s| s.gyro.norm() > 1e-6);
    if accumulated_imu.len() < 2 || !has_valid_gyro {
        return None;
    }
    let mut preint = ImuPreintegration::new();
    preint.integrate(accumulated_imu, &Vector3::zeros(), &Vector3::zeros());
    let dt = preint.dt;
    if dt > 2.0 && (preint.dp.norm() / dt.max(1e-3)) < 0.5 {
        preint.dp = Vector3::zeros();
        preint.dv = Vector3::zeros();
    }
    Some(preint)
}

fn run_forward_pass(
    rover_epochs: &[gneiss_core::obs::EpochObs],
    base_index: &BTreeMap<u32, &gneiss_core::obs::EpochObs>,
    base_pos: Vector3<f64>,
    imu_samples: &[ImuSample],
    rover_init_pos: Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
) -> BTreeMap<u32, (Vector3<f64>, usize)> {
    let rtk_config = gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([rover_init_pos.x, rover_init_pos.y, rover_init_pos.z]),
        ..Default::default()
    };
    let config = EngineConfig::Rtk(rtk_config);
    let mut engine = SwfgEngine::new(&config, ephemerides.to_vec());
    if let Some(k) = klobuchar {
        engine.set_klobuchar(k.0, k.1);
    }

    let mut solutions = BTreeMap::new();
    let mut last_imu_idx = 0usize;
    let mut accumulated_imu = Vec::new();

    for (i, epoch) in rover_epochs.iter().enumerate() {
        let current_us = (epoch.time.tow * 1_000_000.0) as u32;
        while last_imu_idx < imu_samples.len() && imu_samples[last_imu_idx].time_us <= current_us {
            accumulated_imu.push(imu_samples[last_imu_idx]);
            last_imu_idx += 1;
        }

        let exact_ms = (epoch.time.tow * 1000.0).round() as u32;
        let Some(base) = base_index.get(&exact_ms).copied() else {
            continue;
        };

        let preint_opt = build_preintegration(&accumulated_imu);
        accumulated_imu.clear();

        let tow_sec = epoch.time.tow.floor() as u32;
        if let Ok(sol) = engine.process_rtk_epoch_with_imu(epoch, base, base_pos, preint_opt) {
            if sol.n_satellites >= 4 {
                solutions.insert(tow_sec, (sol.position_ecef, sol.n_satellites));
            }
        }

        if i % 1000 == 0 || i == rover_epochs.len() - 1 {
            println!("Forward pass: epoch {}/{} (TOW {}) - solutions: {}", i + 1, rover_epochs.len(), tow_sec, solutions.len());
        }
    }
    solutions
}

fn run_backward_pass(
    rover_epochs: &[gneiss_core::obs::EpochObs],
    base_index: &BTreeMap<u32, &gneiss_core::obs::EpochObs>,
    base_pos: Vector3<f64>,
    rover_final_pos: Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
) -> BTreeMap<u32, (Vector3<f64>, usize)> {
    let rtk_config = gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([rover_final_pos.x, rover_final_pos.y, rover_final_pos.z]),
        ..Default::default()
    };
    let config = EngineConfig::Rtk(rtk_config);
    let mut engine = SwfgEngine::new(&config, ephemerides.to_vec());
    if let Some(k) = klobuchar {
        engine.set_klobuchar(k.0, k.1);
    }

    let mut solutions = BTreeMap::new();
    let mut rev_epochs = rover_epochs.to_vec();
    rev_epochs.reverse();

    for (i, epoch) in rev_epochs.iter().enumerate() {
        let exact_ms = (epoch.time.tow * 1000.0).round() as u32;
        let Some(base) = base_index.get(&exact_ms).copied() else {
            continue;
        };

        let tow_sec = epoch.time.tow.floor() as u32;
        if let Ok(sol) = engine.process_rtk_epoch_with_imu(epoch, base, base_pos, None) {
            if sol.n_satellites >= 4 {
                solutions.insert(tow_sec, (sol.position_ecef, sol.n_satellites));
            }
        }

        if i % 1000 == 0 || i == rev_epochs.len() - 1 {
            println!("Backward pass: epoch {}/{} (TOW {}) - solutions: {}", i + 1, rev_epochs.len(), tow_sec, solutions.len());
        }
    }
    solutions
}

fn print_pass_stats(name: &str, map: &BTreeMap<u32, (Vector3<f64>, usize)>, truth: &BTreeMap<u32, (f64, f64, f64)>) {
    let mut h_errors = Vec::new();
    let mut errors_3d = Vec::new();
    for (tow, (p, _)) in map {
        if let Some(&(tx, ty, tz)) = truth.get(tow) {
            let truth_pos = Vector3::new(tx, ty, tz);
            let truth_llh = gneiss_core::coords::ecef_to_llh(truth_pos);
            let enu = gneiss_core::coords::ecef_delta_to_enu(*p, truth_pos, truth_llh);
            h_errors.push((enu.x * enu.x + enu.y * enu.y).sqrt());
            errors_3d.push(enu.norm());
        }
    }
    if h_errors.is_empty() { return; }
    h_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    errors_3d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = h_errors.len();
    let p50 = h_errors[n / 2];
    let p68 = h_errors[(n as f64 * 0.68) as usize];
    let p95 = h_errors[(n as f64 * 0.95) as usize];
    let rms = (h_errors.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    println!("=== {} (N={}) ===", name, n);
    println!("Horizontal error: p50={:.3}m, p68={:.3}m, p95={:.3}m, RMS={:.3}m", p50, p68, p95, rms);
}

fn fuse_solutions(
    fwd_map: &BTreeMap<u32, (Vector3<f64>, usize)>,
    bwd_map: &BTreeMap<u32, (Vector3<f64>, usize)>,
) -> BTreeMap<u32, (Vector3<f64>, usize)> {
    let mut fused = BTreeMap::new();
    let mut all_tows = std::collections::BTreeSet::new();
    for t in fwd_map.keys() { all_tows.insert(*t); }
    for t in bwd_map.keys() { all_tows.insert(*t); }

    for tow_sec in &all_tows {
        let pos_opt = match (fwd_map.get(tow_sec), bwd_map.get(tow_sec)) {
            (Some((fp, ns_f)), Some((bp, ns_b))) => {
                let sep = (fp - bp).norm();
                if sep < 10.0 {
                    Some((0.5 * fp + 0.5 * bp, *ns_f.max(ns_b)))
                } else if ns_f >= ns_b {
                    Some((*fp, *ns_f))
                } else {
                    Some((*bp, *ns_b))
                }
            }
            (Some(f), None) => Some(*f),
            (None, Some(b)) => Some(*b),
            (None, None) => None,
        };
        if let Some(res) = pos_opt {
            fused.insert(*tow_sec, res);
        }
    }
    fused
}

fn main() {
    tracing_subscriber::fmt().with_env_filter("warn").with_target(false).without_time().try_init().ok();

    let dataset = Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba");
    let nav_file = File::open(dataset.join("base.nav")).expect("open base.nav");
    let (ephemerides, klobuchar) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_file)).expect("parse nav");
    let klob_pair = klobuchar.map(|k| (k.alpha, k.beta));

    let ref_csv = std::fs::read_to_string(dataset.join("reference.csv")).expect("read reference.csv");
    let mut truth: BTreeMap<u32, (f64, f64, f64)> = BTreeMap::new();
    for line in ref_csv.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 8 {
            if let (Ok(tow), Ok(x), Ok(y), Ok(z)) = (parts[0].trim().parse::<f64>(), parts[5].trim().parse::<f64>(), parts[6].trim().parse::<f64>(), parts[7].trim().parse::<f64>()) {
                truth.insert(tow.floor() as u32, (x, y, z));
            }
        }
    }

    let rover_file = File::open(dataset.join("rover_trimble.obs")).expect("open rover_trimble.obs");
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rover_file)).expect("parse rover RINEX");
    let base_file = File::open(dataset.join("base_trimble.obs")).expect("open base_trimble.obs");
    let (base_epochs, base_approx) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_file)).expect("parse base RINEX");

    let default_base_arp = Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let default_base_llh = gneiss_core::coords::ecef_to_llh(default_base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(default_base_llh).transpose();
    let computed_base = default_base_arp + ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855);
    let base_pos = base_approx.approx_position.map(|p| Vector3::new(p[0], p[1], p[2])).unwrap_or(computed_base);

    let base_index: BTreeMap<u32, &gneiss_core::obs::EpochObs> = base_epochs.iter()
        .map(|e| ((e.time.tow * 1000.0).round() as u32, e))
        .collect();

    let max_epochs = std::env::var("MAX_EPOCHS").ok().and_then(|s| s.parse::<usize>().ok()).unwrap_or(usize::MAX);
    let selected_rover = &rover_epochs[..rover_epochs.len().min(max_epochs)];
    let imu_samples = parse_imu_csv(&dataset.join("imu.csv"));

    let rover_init_pos = selected_rover.first()
        .and_then(|e| gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(e, &ephemerides))
        .unwrap_or(base_pos);
    let rover_final_pos = selected_rover.last()
        .and_then(|e| gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(e, &ephemerides))
        .unwrap_or(base_pos);

    println!("=== Odaiba Parallel INS/GNSS Benchmark ===");
    println!("Rover Epochs: {}, IMU Samples: {}, Truth Points: {}", selected_rover.len(), imu_samples.len(), truth.len());

    let (forward_solutions, backward_solutions) = rayon::join(
        || run_forward_pass(selected_rover, &base_index, base_pos, &imu_samples, rover_init_pos, &ephemerides, klob_pair),
        || run_backward_pass(selected_rover, &base_index, base_pos, rover_final_pos, &ephemerides, klob_pair),
    );

    let fused_solutions = fuse_solutions(&forward_solutions, &backward_solutions);
    print_pass_stats("Forward Pass", &forward_solutions, &truth);
    print_pass_stats("Backward Pass", &backward_solutions, &truth);
    print_pass_stats("Smoothed INS/GNSS", &fused_solutions, &truth);
}
