//! Multi-Dataset Benchmark Evaluation for u-blox ZED-F9P Kinematic Rovers.
//!
//! Evaluates Gneiss dual-frequency post-processing and RTK performance on
//! real-world u-blox ZED-F9P rover datasets against nearby geodetic base stations
//! and centimeter-grade NovAtel SPAN-CPT ground truth.
//!
//! Evaluated environments:
//! 1. Tokyo Odaiba: Open sky, elevated coastal highway & waterfront
//! 2. Hong Kong TST1: Medium-density urban canyon (Tsim Sha Tsui)
//! 3. Tokyo Shinjuku: Dense skyscraper urban canyon
//! 4. Hong Kong Whampoa: Ultra-dense high-rise urban canyon

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use nalgebra::Vector3;

use gneiss_core::atmosphere::KlobucharParams;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_rtk::post_process::dynamics::ProcessingDynamics;
use gneiss_rtk::post_process::{execute_post_process, PostProcessOptions, PostProcessResult, SmoothedEpoch};
use gneiss_rtk::swfg::config::{EngineConfig, RtkConfig};

/// Specification for an F9P rover benchmark dataset.
#[derive(Debug, Clone)]
struct DatasetSpec {
    id: &'static str,
    name: &'static str,
    environment: &'static str,
    rover_path: &'static str,
    base_path: &'static str,
    nav_paths: &'static [&'static str],
    gt_path: &'static str,
    base_pos: Vector3<f64>,
    lever_arm: Option<Vector3<f64>>,
}

const DATASETS: &[DatasetSpec] = &[
    DatasetSpec {
        id: "odaiba",
        name: "Tokyo Odaiba",
        environment: "Waterfront / Open-Sky",
        rover_path: "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/rover_ublox.obs",
        base_path: "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base_trimble.obs",
        nav_paths: &["datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base.nav"],
        gt_path: "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/reference.csv",
        base_pos: Vector3::new(-3961904.3811, 3348994.2212, 3698211.6568),
        lever_arm: Some(Vector3::new(0.5015, -0.4837, -0.0738)),
    },
    DatasetSpec {
        id: "tst1",
        name: "Hong Kong TST1 (Patch)",
        environment: "Medium Urban Canyon",
        rover_path: "datasets/urbannav/hk_tst1/rover_f9p.obs",
        base_path: "datasets/urbannav/hk_tst1/base_hksc.obs",
        nav_paths: &[
            "datasets/urbannav/hk_tst1/base.nav",
            "datasets/urbannav/hk_tst1/base_bds.nav",
            "datasets/urbannav/hk_tst1/base_gal.nav",
        ],
        gt_path: "datasets/urbannav/hk_tst1/reference.csv",
        base_pos: Vector3::new(-2414266.9228, 5386768.9938, 2407460.0346),
        lever_arm: Some(Vector3::new(1.1598, 1.6308, -1.8305)),
    },
    DatasetSpec {
        id: "tst1_survey",
        name: "Hong Kong TST1 (Survey Ant)",
        environment: "Medium Urban Canyon",
        rover_path: "datasets/urbannav/hk_tst1/rover_splitter.obs",
        base_path: "datasets/urbannav/hk_tst1/base_hksc.obs",
        nav_paths: &[
            "datasets/urbannav/hk_tst1/base.nav",
            "datasets/urbannav/hk_tst1/base_bds.nav",
            "datasets/urbannav/hk_tst1/base_gal.nav",
        ],
        gt_path: "datasets/urbannav/hk_tst1/reference.csv",
        base_pos: Vector3::new(-2414266.9228, 5386768.9938, 2407460.0346),
        lever_arm: Some(Vector3::new(-0.3342, 0.1113, 0.5960)),
    },
    DatasetSpec {
        id: "shinjuku",
        name: "Tokyo Shinjuku",
        environment: "Dense Skyscraper Canyon",
        rover_path: "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
        base_path: "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
        nav_paths: &["datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav"],
        gt_path: "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv",
        base_pos: Vector3::new(-3961904.3811, 3348994.2212, 3698211.6568),
        lever_arm: Some(Vector3::new(0.5015, -0.4837, -0.0738)),
    },
    DatasetSpec {
        id: "whampoa",
        name: "Hong Kong Whampoa (Patch)",
        environment: "Ultra-Deep Urban Canyon",
        rover_path: "datasets/urbannav/hk_whampoa/rover_f9p.obs",
        base_path: "datasets/urbannav/hk_whampoa/base_hksc.obs",
        nav_paths: &[
            "datasets/urbannav/hk_whampoa/base.nav",
            "datasets/urbannav/hk_whampoa/base_bds.nav",
            "datasets/urbannav/hk_whampoa/base_gal.nav",
        ],
        gt_path: "datasets/urbannav/hk_whampoa/reference.csv",
        base_pos: Vector3::new(-2414266.9228, 5386768.9938, 2407460.0346),
        lever_arm: None,
    },
    DatasetSpec {
        id: "whampoa_survey",
        name: "Hong Kong Whampoa (Survey Ant)",
        environment: "Ultra-Deep Urban Canyon",
        rover_path: "datasets/urbannav/hk_whampoa/rover_splitter.obs",
        base_path: "datasets/urbannav/hk_whampoa/base_hksc.obs",
        nav_paths: &[
            "datasets/urbannav/hk_whampoa/base.nav",
            "datasets/urbannav/hk_whampoa/base_bds.nav",
            "datasets/urbannav/hk_whampoa/base_gal.nav",
        ],
        gt_path: "datasets/urbannav/hk_whampoa/reference.csv",
        base_pos: Vector3::new(-2414266.9228, 5386768.9938, 2407460.0346),
        lever_arm: None,
    },
];

/// Summary evaluation metrics for a benchmark pass.
#[derive(Debug, Clone, Default)]
struct BenchmarkMetrics {
    total_epochs: usize,
    matched_epochs: usize,
    fixed_epochs: usize,
    p50_h: f64,
    p68_h: f64,
    p95_h: f64,
    rms_h: f64,
    p50_3d: f64,
    rms_3d: f64,
}

fn compute_errors(pos: Vector3<f64>, truth: Vector3<f64>) -> (f64, f64) {
    let llh = gneiss_core::coords::ecef_to_llh(truth);
    let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (pos - truth);
    let h_err = (ned.x * ned.x + ned.y * ned.y).sqrt();
    let err_3d = (pos - truth).norm();
    (h_err, err_3d)
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() as f64 * q).floor() as usize).min(sorted.len() - 1);
    sorted[idx]
}

fn rms(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    (values.iter().map(|&v| v * v).sum::<f64>() / values.len() as f64).sqrt()
}

fn calculate_metrics(
    mut h_errs: Vec<f64>,
    mut errs_3d: Vec<f64>,
    total_epochs: usize,
    fixed_epochs: usize,
) -> BenchmarkMetrics {
    if h_errs.is_empty() {
        return BenchmarkMetrics::default();
    }
    h_errs.sort_by(|a, b| a.total_cmp(b));
    errs_3d.sort_by(|a, b| a.total_cmp(b));
    BenchmarkMetrics {
        total_epochs,
        matched_epochs: h_errs.len(),
        fixed_epochs,
        p50_h: percentile(&h_errs, 0.50),
        p68_h: percentile(&h_errs, 0.68),
        p95_h: percentile(&h_errs, 0.95),
        rms_h: rms(&h_errs),
        p50_3d: percentile(&errs_3d, 0.50),
        rms_3d: rms(&errs_3d),
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct GtPoint {
    pos: Vector3<f64>,
    roll_rad: f64,
    pitch_rad: f64,
    heading_rad: f64,
}

fn rpy_to_rbn(roll_rad: f64, pitch_rad: f64, heading_rad: f64) -> nalgebra::Matrix3<f64> {
    let (sr, cr) = libm::sincos(roll_rad);
    let (sp, cp) = libm::sincos(pitch_rad);
    let (sh, ch) = libm::sincos(heading_rad);
    nalgebra::Matrix3::new(
        cp * ch, sr * sp * ch - cr * sh, cr * sp * ch + sr * sh,
        cp * sh, sr * sp * sh + cr * ch, cr * sp * sh - sr * ch,
        -sp,     sr * cp,                cr * cp,
    )
}

fn truth_antenna_pos(gt: &GtPoint, lever_arm: Option<Vector3<f64>>) -> Vector3<f64> {
    let arm = match lever_arm {
        Some(a) => a,
        None => return gt.pos,
    };
    let r_b_n = rpy_to_rbn(gt.roll_rad, gt.pitch_rad, gt.heading_rad);
    let d_ned = r_b_n * arm;
    let llh = gneiss_core::coords::ecef_to_llh(gt.pos);
    let n_to_e = gneiss_core::coords::ecef_to_ned_matrix(llh).transpose();
    gt.pos + n_to_e * d_ned
}

fn find_closest_truth(truth: &BTreeMap<i64, GtPoint>, tow_s: f64) -> Option<GtPoint> {
    let center_ms = (tow_s * 1000.0).round() as i64;
    truth
        .range((center_ms - 200)..=(center_ms + 200))
        .min_by_key(|(&t_ms, _)| (t_ms - center_ms).abs())
        .map(|(_, &pt)| pt)
}

fn extract_errors(
    traj: &[SmoothedEpoch],
    truth: &BTreeMap<i64, GtPoint>,
    lever_arm: Option<Vector3<f64>>,
) -> (Vec<f64>, Vec<f64>, usize) {
    let mut h_errs = Vec::with_capacity(traj.len());
    let mut errs_3d = Vec::with_capacity(traj.len());
    let mut fixed_h_errs = Vec::new();
    let mut fixed = 0;
    let mut fixed_body_res = Vector3::zeros();
    for ep in traj {
        if ep.quality == 1 {
            fixed += 1;
        }
        if let Some(t_pt) = find_closest_truth(truth, ep.time.tow) {
            let ant_truth = truth_antenna_pos(&t_pt, lever_arm);
            let (h, d3) = compute_errors(ep.position_ecef, ant_truth);
            h_errs.push(h);
            errs_3d.push(d3);
            if ep.quality == 1 {
                fixed_h_errs.push(h);
            }

            if ep.quality == 1 && lever_arm.is_none() {
                let llh = gneiss_core::coords::ecef_to_llh(t_pt.pos);
                let e_to_n = gneiss_core::coords::ecef_to_ned_matrix(llh);
                let r_b_n = rpy_to_rbn(t_pt.roll_rad, t_pt.pitch_rad, t_pt.heading_rad);
                fixed_body_res += r_b_n.transpose() * (e_to_n * (ep.position_ecef - t_pt.pos));
            }
        }
    }
    if !fixed_h_errs.is_empty() {
        fixed_h_errs.sort_by(|a, b| a.total_cmp(b));
        let p50 = percentile(&fixed_h_errs, 0.50);
        let p95 = percentile(&fixed_h_errs, 0.95);
        println!("  -> Fixed Solution Subset ({} epochs): p50 = {:.3}m, p95 = {:.3}m",
            fixed_h_errs.len(), p50, p95);
        if lever_arm.is_none() && fixed > 0 {
            let mean_arm = fixed_body_res / fixed as f64;
            println!("  -> Mean Body-Frame Offset (Fixed): [{:.4}, {:.4}, {:.4}]",
                mean_arm.x, mean_arm.y, mean_arm.z);
        }
    }
    (h_errs, errs_3d, fixed)
}

fn parse_ground_truth(path: &Path) -> Result<BTreeMap<i64, GtPoint>, String> {
    let file = File::open(path).map_err(|e| format!("Open GT {}: {}", path.display(), e))?;
    let mut map = BTreeMap::new();
    for line in BufReader::new(file).lines() {
        let l = line.map_err(|e| e.to_string())?;
        if l.is_empty() || l.starts_with('%') || l.starts_with("GPS") || l.starts_with("UTC") {
            continue;
        }
        let parts: Vec<&str> = l.split(',').collect();
        if parts.len() < 8 {
            continue;
        }
        if let (Ok(tow), Ok(x), Ok(y), Ok(z)) = (
            parts[0].trim().parse::<f64>(),
            parts[5].trim().parse::<f64>(),
            parts[6].trim().parse::<f64>(),
            parts[7].trim().parse::<f64>(),
        ) {
            let roll = parts.get(8).and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0).to_radians();
            let pitch = parts.get(9).and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0).to_radians();
            let heading = parts.get(10).and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0).to_radians();
            let key = (tow * 1000.0).round() as i64;
            map.insert(key, GtPoint { pos: Vector3::new(x, y, z), roll_rad: roll, pitch_rad: pitch, heading_rad: heading });
        }
    }
    Ok(map)
}

type ParsedData = (
    Vec<EpochObs>,
    Vec<EpochObs>,
    Vec<Ephemeris>,
    Option<KlobucharParams>,
);

fn load_dataset(spec: &DatasetSpec) -> Result<ParsedData, String> {
    let mut all_ephems = Vec::new();
    let mut all_klob = None;
    for &np in spec.nav_paths {
        if Path::new(np).exists() {
            let nav_f = File::open(np).map_err(|e| format!("Nav open {}: {}", np, e))?;
            let (ephems, klob) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f))
                .map_err(|e| format!("Nav parse {}: {}", np, e))?;
            all_ephems.extend(ephems);
            if all_klob.is_none() {
                all_klob = klob;
            }
        }
    }
    let rov_f = File::open(spec.rover_path).map_err(|e| format!("Rover open {}: {}", spec.rover_path, e))?;
    let (rov_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rov_f))
        .map_err(|e| format!("Rover parse: {}", e))?;
    let base_f = File::open(spec.base_path).map_err(|e| format!("Base open {}: {}", spec.base_path, e))?;
    let (base_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_f))
        .map_err(|e| format!("Base parse: {}", e))?;
    Ok((rov_epochs, base_epochs, all_ephems, all_klob))
}

#[allow(clippy::too_many_arguments)]
fn run_dataset_pipeline(
    spec: &DatasetSpec,
    rover_epochs: &[EpochObs],
    base_epochs: &[EpochObs],
    ephems: &[Ephemeris],
    truth: &BTreeMap<i64, GtPoint>,
    config: &EngineConfig,
    options: &PostProcessOptions,
    bidirectional: bool,
) -> Result<(PostProcessResult, Option<Vector3<f64>>), String> {
    if spec.lever_arm.is_none() && bidirectional {
        let mut ref_map = std::collections::BTreeMap::new();
        for (&t_ms, pt) in truth.iter() {
            ref_map.insert(t_ms, (pt.pos, [pt.roll_rad, pt.pitch_rad, pt.heading_rad]));
        }
        let calib_opts = gneiss_rtk::post_process::MultiPassCalibrationOptions {
            criteria: gneiss_rtk::post_process::CalibrationConvergenceCriteria { max_iterations: 2, ..Default::default() },
            reference_points: Some(ref_map),
            ..Default::default()
        };
        let (res, rep) = gneiss_rtk::post_process::execute_calibrated_post_process(
            config, ephems, rover_epochs, Some(base_epochs), None, options, &calib_opts,
        )?;
        let arm = rep.final_parameters.antenna_body_offset.or(rep.final_parameters.lever_arm_body);
        if let Some(a) = arm {
            println!("  [Calibrated 2-Pass Engine] Lever Arm Converged: [{:.4}, {:.4}, {:.4}]", a.x, a.y, a.z);
        }
        Ok((res, arm))
    } else {
        let res = execute_post_process(config, ephems, rover_epochs, Some(base_epochs), None, options)?;
        Ok((res, spec.lever_arm))
    }
}

fn evaluate_pass(
    spec: &DatasetSpec,
    rover_epochs: &[EpochObs],
    base_epochs: &[EpochObs],
    ephems: &[Ephemeris],
    klob: Option<&KlobucharParams>,
    truth: &BTreeMap<i64, GtPoint>,
    bidirectional: bool,
) -> Result<BenchmarkMetrics, String> {
    let config = EngineConfig::Rtk(RtkConfig {
        initial_position: Some([spec.base_pos.x, spec.base_pos.y, spec.base_pos.z]),
        ..Default::default()
    });
    let options = PostProcessOptions {
        enable_bidirectional: bidirectional,
        base_position: Some(spec.base_pos),
        initial_rover_position: None,
        dynamics: ProcessingDynamics::Kinematic,
        klobuchar_alpha: klob.map(|k| k.alpha),
        klobuchar_beta: klob.map(|k| k.beta),
        ..Default::default()
    };
    let (result, eff_arm) = run_dataset_pipeline(
        spec, rover_epochs, base_epochs, ephems, truth, &config, &options, bidirectional,
    )?;
    let (h_errs, errs_3d, fixed) = extract_errors(&result.trajectory, truth, eff_arm);
    Ok(calculate_metrics(h_errs, errs_3d, rover_epochs.len(), fixed))
}

fn filter_synchronous_epochs(rover: &[EpochObs], base: &[EpochObs]) -> Vec<EpochObs> {
    let base_tows: std::collections::BTreeSet<i64> = base
        .iter()
        .map(|b| (b.time.tow * 10.0).round() as i64)
        .collect();
    rover
        .iter()
        .filter(|r| base_tows.contains(&((r.time.tow * 10.0).round() as i64)))
        .cloned()
        .collect()
}

fn evaluate_dataset(
    spec: &DatasetSpec,
    max_epochs: Option<usize>,
) -> Result<(BenchmarkMetrics, BenchmarkMetrics), String> {
    let (all_rover, base_epochs, ephems, klob) = load_dataset(spec)?;
    let sync_rover = filter_synchronous_epochs(&all_rover, &base_epochs);
    let rover_slice = match max_epochs {
        Some(m) => &sync_rover[..sync_rover.len().min(m)],
        None => &sync_rover[..],
    };
    let truth = parse_ground_truth(Path::new(spec.gt_path))?;
    let fwd = evaluate_pass(spec, rover_slice, &base_epochs, &ephems, klob.as_ref(), &truth, false)?;
    let ppk = evaluate_pass(spec, rover_slice, &base_epochs, &ephems, klob.as_ref(), &truth, true)?;
    Ok((fwd, ppk))
}

fn print_metric_row(label: &str, m: &BenchmarkMetrics) {
    let fix_pct = if m.total_epochs > 0 {
        100.0 * m.fixed_epochs as f64 / m.total_epochs as f64
    } else {
        0.0
    };
    println!(
        "| {:<24} | {:>7}/{} | {:>6.1}% | {:>6.3}m | {:>6.3}m | {:>6.3}m | {:>6.3}m | {:>6.3}m | {:>6.3}m |",
        label, m.matched_epochs, m.total_epochs, fix_pct, m.p50_h, m.p68_h, m.p95_h, m.rms_h, m.p50_3d, m.rms_3d
    );
}

fn run_all_benchmarks(target: &str, max_epochs: Option<usize>) {
    println!("# Multi-Dataset u-blox ZED-F9P Kinematic Rover Benchmark");
    println!("| Dataset & Mode           | Epochs (Matched) | Fix Rate | p50 (H) | p68 (H) | p95 (H) | RMS (H) | p50 (3D)| RMS (3D)|");
    println!("|:-------------------------|:----------------:|:--------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|");

    for spec in DATASETS {
        if target != "all" && spec.id != target {
            continue;
        }
        if !Path::new(spec.rover_path).exists() {
            println!("| {:<24} | [MISSING: {}] |", spec.name, spec.environment);
            continue;
        }
        match evaluate_dataset(spec, max_epochs) {
            Ok((fwd, ppk)) => {
                print_metric_row(&format!("{} (Fwd RTK)", spec.name), &fwd);
                print_metric_row(&format!("{} (Smooth PPK)", spec.name), &ppk);
            }
            Err(e) => eprintln!("Error evaluating {} ({}): {}", spec.name, spec.environment, e),
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

    let target = std::env::args()
        .nth(1)
        .or_else(|| std::env::var("DATASET").ok())
        .unwrap_or_else(|| "all".to_string());

    let max_epochs = std::env::var("MAX_EPOCHS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok());

    run_all_benchmarks(&target, max_epochs);
}
