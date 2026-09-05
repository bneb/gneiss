//! Benchmark evaluation of Qinertia-grade 4-pass RTK/PPK post-processing engine.
//!
//! Evaluates forward-only RTK vs 4-pass bidirectional smoothed PPK on real datasets.
//!
//! Usage: cargo run --release --bin eval_qinertia_ppk

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use nalgebra::Vector3;

use gneiss_core::atmosphere::KlobucharParams;
use gneiss_core::time::GpsTime;
use gneiss_rtk::post_process::{execute_post_process, PostProcessOptions, SmoothedEpoch};
use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::imu_preintegration::ImuSample;

struct DatasetSpec {
    name: &'static str,
    dir: &'static str,
    rover_file: &'static str,
    base_file: &'static str,
    nav_file: &'static str,
    truth_file: &'static str,
    imu_file: Option<&'static str>,
    base_pos_override: Option<Vector3<f64>>,
    max_epochs: usize,
    dynamics: gneiss_rtk::post_process::ProcessingDynamics,
    widelane_ar: bool,
    enable_glonass: bool,
}

fn parse_truth(path: &Path) -> BTreeMap<u32, Vector3<f64>> {
    let mut truth = BTreeMap::new();
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(_) => return truth,
    };

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('%') {
            continue;
        }
        if trimmed.contains(',') {
            parse_csv_truth_line(trimmed, &mut truth);
        } else {
            parse_pos_truth_line(trimmed, &mut truth);
        }
    }
    truth
}

fn parse_csv_truth_line(line: &str, truth: &mut BTreeMap<u32, Vector3<f64>>) {
    let parts: Vec<&str> = line.split(',').collect();
    if parts.len() >= 8 {
        if let (Ok(tow), Ok(x), Ok(y), Ok(z)) = (
            parts[0].trim().parse::<f64>(),
            parts[5].trim().parse::<f64>(),
            parts[6].trim().parse::<f64>(),
            parts[7].trim().parse::<f64>(),
        ) {
            truth.insert(tow.round() as u32, Vector3::new(x, y, z));
        }
    }
}

fn parse_pos_truth_line(line: &str, truth: &mut BTreeMap<u32, Vector3<f64>>) {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() >= 5 {
        let date_parts: Vec<&str> = parts[0].split('/').collect();
        let time_parts: Vec<&str> = parts[1].split(':').collect();
        if date_parts.len() == 3 && time_parts.len() == 3 {
            if let (Ok(y), Ok(m), Ok(d), Ok(hr), Ok(min), Ok(sec), Ok(px), Ok(py), Ok(pz)) = (
                date_parts[0].parse::<i32>(),
                date_parts[1].parse::<i32>(),
                date_parts[2].parse::<i32>(),
                time_parts[0].parse::<i32>(),
                time_parts[1].parse::<i32>(),
                time_parts[2].parse::<f64>(),
                parts[2].parse::<f64>(),
                parts[3].parse::<f64>(),
                parts[4].parse::<f64>(),
            ) {
                let gps_time = GpsTime::from_calendar(y, m, d, hr, min, sec);
                truth.insert(gps_time.tow.round() as u32, Vector3::new(px, py, pz));
            }
        }
    }
}

fn parse_imu(imu_path: &Path) -> Vec<ImuSample> {
    let mut samples = Vec::new();
    let file = match File::open(imu_path) {
        Ok(f) => f,
        Err(_) => return samples,
    };
    let reader = BufReader::new(file);
    for (i, line) in reader.lines().enumerate() {
        if i == 0 { continue; }
        if let Ok(l) = line {
            let p: Vec<&str> = l.split(',').map(|s| s.trim()).collect();
            if p.len() >= 8 {
                let tow: f64 = p[0].parse().unwrap_or(0.0);
                let ax: f64 = p[2].parse().unwrap_or(0.0);
                let ay: f64 = p[3].parse().unwrap_or(0.0);
                let az: f64 = p[4].parse().unwrap_or(0.0);
                let gx: f64 = p[5].parse().unwrap_or(0.0);
                let gy: f64 = p[6].parse().unwrap_or(0.0);
                let gz: f64 = p[7].parse().unwrap_or(0.0);
                samples.push(ImuSample {
                    accel: Vector3::new(ax, ay, az),
                    gyro: Vector3::new(gx, gy, gz),
                    time_us: (tow * 1_000_000.0) as u32,
                });
            }
        }
    }
    samples
}

fn compute_horizontal_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    let llh = gneiss_core::coords::ecef_to_llh(truth);
    let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (pos - truth);
    (ned.x * ned.x + ned.y * ned.y).sqrt()
}

fn compute_3d_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    (pos - truth).norm()
}

fn format_stats(
    out: &mut String,
    name: &str,
    mut h_errs: Vec<f64>,
    mut d3_errs: Vec<f64>,
    fix_count: usize,
    total_count: usize,
) {
    use std::fmt::Write;
    if h_errs.is_empty() { return; }
    h_errs.sort_by(|a, b| a.total_cmp(b));
    d3_errs.sort_by(|a, b| a.total_cmp(b));
    let n = h_errs.len();
    let p10 = h_errs[(n as f64 * 0.10) as usize];
    let p25 = h_errs[(n as f64 * 0.25) as usize];
    let p50 = h_errs[n / 2];
    let p68 = h_errs[(n as f64 * 0.68) as usize];
    let p75 = h_errs[(n as f64 * 0.75) as usize];
    let p90 = h_errs[(n as f64 * 0.90) as usize];
    let p95 = h_errs[(n as f64 * 0.95) as usize];
    let p99 = h_errs[(n as f64 * 0.99) as usize];
    let max = *h_errs.last().unwrap_or(&0.0);
    let mean = h_errs.iter().sum::<f64>() / n as f64;
    let rms = (h_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();

    let p50_3d = d3_errs[n / 2];
    let p95_3d = d3_errs[(n as f64 * 0.95) as usize];
    let rms_3d = (d3_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    let fix_pct = (fix_count as f64 / total_count.max(1) as f64) * 100.0;

    let _ = writeln!(out, "=== {} (N={}, Fixed={}/{} [{:.1}%]) ===", name, n, fix_count, total_count, fix_pct);
    let _ = writeln!(out, "Horizontal Error:  p50={:.3}m,  p68={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50, p68, p95, rms);
    let _ = writeln!(out, "  Horizontal CDF:  p10={:.3}m, p25={:.3}m, p50={:.3}m, p68={:.3}m, p75={:.3}m, p90={:.3}m, p95={:.3}m, p99={:.3}m, max={:.3}m, mean={:.3}m",
        p10, p25, p50, p68, p75, p90, p95, p99, max, mean);
    let _ = writeln!(out, "3D Position Error: p50={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50_3d, p95_3d, rms_3d);
}

fn collect_trajectory_stats(
    trajectory: &[SmoothedEpoch],
    truth: &BTreeMap<u32, Vector3<f64>>,
) -> (Vec<f64>, Vec<f64>, usize) {
    let mut h_errs = Vec::new();
    let mut d3_errs = Vec::new();
    let mut fix_count = 0;
    for ep in trajectory {
        let tow = ep.time.tow.round() as u32;
        if ep.quality == 1 { fix_count += 1; }
        if let Some(&t) = truth.get(&tow) {
            let h = compute_horizontal_error(ep.position_ecef, t);
            let d3 = compute_3d_error(ep.position_ecef, t);
            if h < 100.0 {
                h_errs.push(h);
                d3_errs.push(d3);
            }
        }
    }
    (h_errs, d3_errs, fix_count)
}

fn make_post_process_options(
    spec: &DatasetSpec,
    base_pos: Vector3<f64>,
    rover_init_pos: Option<Vector3<f64>>,
    klob: Option<&KlobucharParams>,
    bidirectional: bool,
) -> PostProcessOptions {
    PostProcessOptions {
        enable_bidirectional: bidirectional,
        base_position: Some(base_pos),
        initial_rover_position: rover_init_pos,
        klobuchar_alpha: klob.map(|k| k.alpha),
        klobuchar_beta: klob.map(|k| k.beta),
        q_accel: None,
        widelane_ar: spec.widelane_ar,
        tropo_gradients: false,
        network_sat_upd: None,
        receiver_pcv: None,
        dynamics: spec.dynamics,
        enable_glonass: spec.enable_glonass,
        continuity_gate: false,
        precise_orbits: None,
        precise_clocks: None,
        sinex_bias: None,
        antex_database: None,
    }
}

fn evaluate_dataset_spec(spec: &DatasetSpec) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    let _ = writeln!(out, "\n========================================================");
    let _ = writeln!(out, "Evaluating Dataset: {}", spec.name);
    let _ = writeln!(out, "========================================================");

    let dir = Path::new(spec.dir);
    let Ok(nav_f) = File::open(dir.join(spec.nav_file)) else { return out; };
    let Ok((ephemerides, klobuchar)) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f)) else { return out; };
    let truth = parse_truth(&dir.join(spec.truth_file));
    let imu_samples = spec.imu_file.map(|f| parse_imu(&dir.join(f)));

    let Ok(rov_f) = File::open(dir.join(spec.rover_file)) else { return out; };
    let Ok((rover_epochs, rover_header)) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rov_f)) else { return out; };
    let Ok(base_f) = File::open(dir.join(spec.base_file)) else { return out; };
    let Ok((base_epochs, base_header)) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_f)) else { return out; };

    let base_pos = spec.base_pos_override.or_else(|| {
        base_header.approx_position.map(|p| Vector3::new(p[0], p[1], p[2]))
    }).unwrap_or_else(|| Vector3::new(0.0, 0.0, 0.0));

    let rover_init_pos = rover_header.approx_position.map(|p| Vector3::new(p[0], p[1], p[2])).or_else(|| {
        gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(&rover_epochs[0], &ephemerides)
    });

    let config = EngineConfig::Rtk(gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([base_pos.x, base_pos.y, base_pos.z]),
        ..Default::default()
    });
    let selected_rover = &rover_epochs[..rover_epochs.len().min(spec.max_epochs)];

    let fwd_opt = make_post_process_options(spec, base_pos, rover_init_pos, klobuchar.as_ref(), false);
    if let Ok(fwd_res) = execute_post_process(&config, &ephemerides, selected_rover, Some(&base_epochs), imu_samples.as_deref(), &fwd_opt) {
        let (h, d3, fix) = collect_trajectory_stats(&fwd_res.trajectory, &truth);
        format_stats(&mut out, "Forward RTK Solution", h, d3, fix, fwd_res.trajectory.len());
    }

    let smooth_opt = make_post_process_options(spec, base_pos, rover_init_pos, klobuchar.as_ref(), true);
    if let Ok(smooth_res) = execute_post_process(&config, &ephemerides, selected_rover, Some(&base_epochs), imu_samples.as_deref(), &smooth_opt) {
        let (h, d3, fix) = collect_trajectory_stats(&smooth_res.trajectory, &truth);
        format_stats(&mut out, "Qinertia-Grade Smoothed PPK Solution", h, d3, fix, smooth_res.trajectory.len());
    }

    out
}

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("warn")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    println!("========================================================");
    println!("Gneiss vs Qinertia-Grade RTK/PPK Post-Processing Harness");
    println!("========================================================");

    let default_base_arp = Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let default_base_llh = gneiss_core::coords::ecef_to_llh(default_base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(default_base_llh).transpose();
    let odaiba_base = default_base_arp + ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855);

    let odaiba_spec = DatasetSpec {
        name: "Odaiba (UrbanNav Trimble 10Hz)",
        dir: "datasets/urbannav/tokyo/Tokyo_Data/Odaiba",
        rover_file: "rover_trimble.obs",
        base_file: "base_trimble.obs",
        nav_file: "base.nav",
        truth_file: "reference.csv",
        imu_file: Some("imu.csv"),
        base_pos_override: Some(odaiba_base),
        max_epochs: 600,
        dynamics: gneiss_rtk::post_process::ProcessingDynamics::Kinematic,
        widelane_ar: false,
        enable_glonass: false,
    };

    let max_f9p = std::env::var("MAX_F9P_EPOCHS")
        .or_else(|_| std::env::var("MAX_EPOCHS"))
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(usize::MAX);

    let f9p_spec = DatasetSpec {
        name: "RTK Explorer F9P (u-blox ZED-F9P Kinematic 1Hz)",
        dir: "datasets/rtkexplorer/sample_1/f9p_ppp_1224",
        rover_file: "rover.obs",
        base_file: "tmg23590.20o",
        nav_file: "BRDC00IGS_R_20203590000_01D_MN.rnx",
        truth_file: "rover_ppk.pos",
        imu_file: None,
        base_pos_override: Some(Vector3::new(-1283434.6250, -4713071.9830, 4090105.0479)),
        max_epochs: max_f9p,
        dynamics: gneiss_rtk::post_process::ProcessingDynamics::Kinematic,
        widelane_ar: false,
        enable_glonass: true,
    };

    let cors_spec = DatasetSpec {
        name: "NGS Geodetic Baseline (TMG2 Base, TMGO Rover, 112.5m, 30s)",
        dir: "datasets/cors_short_baseline",
        rover_file: "tmgo1350.20o",
        base_file: "tmg21350.20o",
        nav_file: "brdc1350.20n",
        truth_file: "tmgo_truth.pos",
        imu_file: None,
        base_pos_override: Some(Vector3::new(-1283433.9360, -4713073.2930, 4090105.0870)),
        max_epochs: 300,
        dynamics: gneiss_rtk::post_process::ProcessingDynamics::Static,
        widelane_ar: false,
        enable_glonass: false,
    };

    let p181_spec = DatasetSpec {
        name: "NOAA CORS Regional Baseline (P181 Base, P224 Rover, 15.0km, 30s)",
        dir: "datasets/cors_short_baseline",
        rover_file: "p2241350.20o",
        base_file: "p1811350.20o",
        nav_file: "brdc1350.20n",
        truth_file: "p224_truth.pos",
        imu_file: None,
        base_pos_override: Some(Vector3::new(-2697941.2851, -4255089.1805, 3898009.7146)),
        max_epochs: 600,
        dynamics: gneiss_rtk::post_process::ProcessingDynamics::Static,
        widelane_ar: false,
        enable_glonass: false,
    };

    use rayon::prelude::*;
    let specs = vec![odaiba_spec, f9p_spec, cors_spec, p181_spec];
    let reports: Vec<String> = specs
        .into_par_iter()
        .filter(|s| Path::new(s.dir).exists())
        .map(|s| evaluate_dataset_spec(&s))
        .collect();

    for report in reports {
        print!("{}", report);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_stats_empty_is_noop() {
        let mut out = String::new();
        format_stats(&mut out, "Test", Vec::new(), Vec::new(), 0, 0);
        assert!(out.is_empty());
    }

    #[test]
    fn test_collect_trajectory_stats_counts_fixes() {
        let mut truth = BTreeMap::new();
        truth.insert(100, Vector3::new(100.0, 200.0, 300.0));
        let ep1 = SmoothedEpoch {
            time: gneiss_core::time::GpsTime::new(2000, 100.0),
            position_ecef: Vector3::new(100.01, 200.0, 300.0),
            velocity_ecef: Some(Vector3::zeros()),
            attitude: None,
            cov_position: nalgebra::Matrix3::identity(),
            std_east: 0.01,
            std_north: 0.01,
            std_up: 0.01,
            separation_3d: 0.005,
            quality: 1,
            n_satellites: 8,
        };
        let (h, d3, fixes) = collect_trajectory_stats(&[ep1], &truth);
        assert_eq!(fixes, 1);
        assert_eq!(h.len(), 1);
        assert_eq!(d3.len(), 1);
        assert!(h[0] < 0.02);
    }
}
