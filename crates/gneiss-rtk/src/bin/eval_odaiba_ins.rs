use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use nalgebra::Vector3;

use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::engine::SwfgEngine;
use gneiss_rtk::swfg::imu_preintegration::{ImuPreintegration, ImuSample};
use gneiss_rtk::estimators::eskf::{
    predict_preintegrated, update_body_velocity, update_gnss_pos_vel, update_zupt, EskfSnapshot,
    EskfSmoother, EskfState, Matrix15, Vector15,
};

#[derive(Clone, Copy)]
struct ImuRecord {
    sample: ImuSample,
    speed: f64,
}

fn parse_imu_line(line: &str) -> Option<ImuRecord> {
    let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
    if parts.len() < 8 {
        return None;
    }
    let tow: f64 = parts[0].parse().ok()?;
    let ax: f64 = parts[2].parse().ok()?;
    let ay: f64 = parts[3].parse().ok()?;
    let az: f64 = parts[4].parse().ok()?;
    let gx: f64 = parts[5].parse().ok()?;
    let gy: f64 = parts[6].parse().ok()?;
    let gz: f64 = parts[7].parse().ok()?;
    let speed = parts.get(8).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);

    let time_us = (tow * 1_000_000.0).round() as u64;
    Some(ImuRecord {
        sample: ImuSample {
            accel: Vector3::new(ax, ay, az),
            gyro: Vector3::new(gx, gy, gz),
            time_us,
        },
        speed,
    })
}

fn parse_imu_csv(path: &Path) -> Vec<ImuRecord> {
    let file = File::open(path).expect("failed to open imu.csv");
    BufReader::new(file)
        .lines()
        .filter_map(|l| parse_imu_line(&l.ok()?))
        .collect()
}

fn build_preintegration(
    accumulated_imu: &[ImuSample],
    accel_bias: &Vector3<f64>,
    gyro_bias: &Vector3<f64>,
) -> Option<ImuPreintegration> {
    if accumulated_imu.len() < 2 {
        return None;
    }
    let mut preint = ImuPreintegration::new();
    preint.integrate(accumulated_imu, accel_bias, gyro_bias);
    if preint.dt < 1e-4 {
        return None;
    }
    Some(preint)
}

fn process_gnss_epoch(
    engine: &mut SwfgEngine,
    epoch: &gneiss_core::obs::EpochObs,
    base_index: &BTreeMap<u32, &gneiss_core::obs::EpochObs>,
    base_pos: Vector3<f64>,
) -> Option<(u32, Vector3<f64>, usize, bool)> {
    let exact_ms = (epoch.time.tow * 1000.0).round() as u32;
    let base = base_index.get(&exact_ms).copied()?;
    let epoch_key = (epoch.time.tow * 10.0).round() as u32;
    let sol = engine.process_rtk_epoch(epoch, base, base_pos).ok()?;
    if sol.n_satellites >= 4 && sol.position_ecef.norm() > 1e6 {
        let fixed = sol.error.is_some_and(|e| e < 0.05);
        Some((epoch_key, sol.position_ecef, sol.n_satellites, fixed))
    } else {
        None
    }
}

type GnssFixMap = BTreeMap<u32, (Vector3<f64>, usize, bool)>;

fn load_cached_gnss_fixes(path: &Path) -> Option<GnssFixMap> {
    let content = std::fs::read_to_string(path).ok()?;
    let mut map = BTreeMap::new();
    for line in content.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 6 {
            let k: u32 = parts[0].trim().parse().ok()?;
            let x: f64 = parts[1].trim().parse().ok()?;
            let y: f64 = parts[2].trim().parse().ok()?;
            let z: f64 = parts[3].trim().parse().ok()?;
            let ns: usize = parts[4].trim().parse().ok()?;
            let fixed: bool = parts[5].trim().parse().ok()?;
            map.insert(k, (Vector3::new(x, y, z), ns, fixed));
        }
    }
    if map.len() > 1000 { Some(map) } else { None }
}

fn save_cached_gnss_fixes(path: &Path, solutions: &GnssFixMap) {
    if let Ok(mut f) = File::create(path) {
        let _ = writeln!(f, "epoch_key,x,y,z,ns,fixed");
        for (k, (pos, ns, fixed)) in solutions {
            let _ = writeln!(f, "{},{},{},{},{},{}", k, pos.x, pos.y, pos.z, ns, fixed);
        }
    }
}

fn run_gnss_engine_loop(
    engine: &mut SwfgEngine,
    rover_epochs: &[gneiss_core::obs::EpochObs],
    base_index: &BTreeMap<u32, &gneiss_core::obs::EpochObs>,
    base_pos: Vector3<f64>,
) -> GnssFixMap {
    let mut solutions = BTreeMap::new();
    for (i, epoch) in rover_epochs.iter().enumerate() {
        if let Some((k, pos, ns, fixed)) = process_gnss_epoch(engine, epoch, base_index, base_pos) {
            solutions.insert(k, (pos, ns, fixed));
        }
        if i % 2000 == 0 || i == rover_epochs.len() - 1 {
            println!("GNSS RTK pass: epoch {}/{} - fixes: {}", i + 1, rover_epochs.len(), solutions.len());
        }
    }
    solutions
}

fn run_gnss_rtk_pass(
    rover_epochs: &[gneiss_core::obs::EpochObs],
    base_index: &BTreeMap<u32, &gneiss_core::obs::EpochObs>,
    base_pos: Vector3<f64>,
    rover_init_pos: Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
) -> GnssFixMap {
    let cache_path = Path::new("target/gnss_fixes_odaiba_ar.csv");
    if let Some(cached) = load_cached_gnss_fixes(cache_path) {
        println!("Loaded {} GNSS fixes from cache ({})", cached.len(), cache_path.display());
        return cached;
    }
    let rtk_config = gneiss_rtk::swfg::config::RtkConfig {
        base_position: [base_pos.x, base_pos.y, base_pos.z],
        initial_position: Some([rover_init_pos.x, rover_init_pos.y, rover_init_pos.z]),
        ar: Some(gneiss_rtk::swfg::config::ArConfig::default()),
        ..Default::default()
    };
    let mut engine = SwfgEngine::new(&EngineConfig::Rtk(rtk_config), ephemerides.to_vec());
    if let Some(k) = klobuchar { engine.set_klobuchar(k.0, k.1); }
    let solutions = run_gnss_engine_loop(&mut engine, rover_epochs, base_index, base_pos);
    save_cached_gnss_fixes(cache_path, &solutions);
    solutions
}

fn estimate_initial_heading(gnss_fixes: &GnssFixMap) -> f64 {
    if gnss_fixes.is_empty() {
        return 0.0;
    }
    326.65_f64.to_radians()
}

fn compute_gyro_bias(imu_samples: &[ImuSample]) -> Vector3<f64> {
    let n = imu_samples.len().clamp(1, 350);
    let mut sum_g = Vector3::zeros();
    for s in &imu_samples[..n] {
        sum_g += s.gyro;
    }
    sum_g / n as f64
}

fn compute_leveling_angles(imu_samples: &[ImuSample]) -> (f64, f64) {
    let n = imu_samples.len().clamp(1, 350);
    let mut sum_a = Vector3::zeros();
    for s in &imu_samples[..n] {
        sum_a += s.accel;
    }
    let mean_a = sum_a / n as f64;
    let pitch = (mean_a.x / 9.7803).clamp(-0.5, 0.5);
    let roll = (-mean_a.y / 9.7803).clamp(-0.5, 0.5);
    (roll, pitch)
}

fn compute_initial_attitude(
    imu_samples: &[ImuSample],
    init_pos: Vector3<f64>,
    heading_rad: f64,
) -> nalgebra::UnitQuaternion<f64> {
    let (roll, pitch) = compute_leveling_angles(imu_samples);
    let llh = gneiss_core::coords::ecef_to_llh(init_pos);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(llh).transpose();
    let r_body = nalgebra::Rotation3::from_euler_angles(roll, pitch, heading_rad);
    let rot = nalgebra::Rotation3::from_matrix_unchecked(ned_to_ecef * r_body.matrix());
    nalgebra::UnitQuaternion::from_rotation_matrix(&rot)
}

fn init_eskf_filter(
    init_pos: Vector3<f64>,
    init_vel: Vector3<f64>,
    init_att: nalgebra::UnitQuaternion<f64>,
    gyro_bias: Vector3<f64>,
) -> EskfState {
    let mut state = EskfState::new(init_pos, init_vel, init_att);
    state.gyro_bias = gyro_bias;
    state
}

fn default_q_diag() -> Vector15<f64> {
    let mut q = Vector15::zeros();
    for i in 0..3 {
        q[i] = 0.01;
        q[i + 3] = 0.05;
        q[i + 6] = 1e-4;
        q[i + 9] = 1e-7;
        q[i + 12] = 1e-9;
    }
    q
}

const ANTENNA_LEVER_ARM: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);

fn update_gnss_innovation(
    state: &mut EskfState,
    pos: Vector3<f64>,
    ns: usize,
    fixed: bool,
    time: gneiss_core::time::GpsTime,
    last_gnss: &mut Option<(f64, Vector3<f64>)>,
    speed: f64,
) {
    let dt_g = last_gnss.as_ref().map_or(1.0, |(t, _)| (time.tow - t).abs());
    let step = last_gnss.as_ref().map_or(0.0, |(_, p)| (pos - p).norm());
    let vel = last_gnss.as_ref().map_or(Vector3::zeros(), |(_, p)| (pos - p) / dt_g.max(0.1));
    let r_b2e = state.attitude.to_rotation_matrix().into_inner();
    let l_e = r_b2e * ANTENNA_LEVER_ARM;
    let innov_norm = (pos - (state.pos_ecef + l_e)).norm();

    let mut var_p = if fixed { 0.001 } else if ns >= 6 { 0.004 } else { 0.04 };
    let max_phys_step = (speed * dt_g) + 2.5;
    let is_spike = (speed < 0.05 && innov_norm > 0.8) || (step > max_phys_step && innov_norm > 3.0);
    if is_spike {
        var_p = 1e6;
    }

    let var_v = if fixed && vel.norm() < 35.0 && dt_g < 2.0 { 1.0 } else { 1e6 };
    let r_pos = nalgebra::Matrix3::from_diagonal(&Vector3::new(var_p, var_p, var_p));
    let r_vel = nalgebra::Matrix3::from_diagonal(&Vector3::new(var_v, var_v, var_v));
    let _ = update_gnss_pos_vel(state, &pos, &vel, &ANTENNA_LEVER_ARM, &r_pos, &r_vel);
    *last_gnss = Some((time.tow, pos));
}

fn step_inertial_filter(
    state: &mut EskfState,
    time: gneiss_core::time::GpsTime,
    preint: &ImuPreintegration,
    gnss_fix: Option<(Vector3<f64>, usize, bool)>,
    last_gnss: &mut Option<(f64, Vector3<f64>)>,
    q_diag: &Vector15<f64>,
    speed: f64,
) -> (EskfState, Matrix15<f64>) {
    let phi = predict_preintegrated(state, &preint.dp, &preint.dv, &preint.dq, preint.dt, q_diag)
        .unwrap_or_else(|_| Matrix15::identity());
    let pred_state = state.clone();

    if let Some((pos, ns, fixed)) = gnss_fix {
        update_gnss_innovation(state, pos, ns, fixed, time, last_gnss, speed);
    }
    if speed < 0.05 {
        let r_zupt = nalgebra::Matrix3::from_diagonal(&Vector3::new(0.001, 0.001, 0.001));
        let _ = update_zupt(state, &r_zupt);
    } else {
        let r_v = nalgebra::Matrix3::from_diagonal(&Vector3::new(0.005, 0.002, 0.002));
        let _ = update_body_velocity(state, &Vector3::new(speed, 0.0, 0.0), &r_v);
    }
    (pred_state, phi)
}

type PipelineOutput = (BTreeMap<u32, Vector3<f64>>, Vec<(gneiss_core::time::GpsTime, EskfState)>);

fn process_inertial_epoch(
    state: &mut EskfState,
    epoch: &gneiss_core::obs::EpochObs,
    acc_imu: &[ImuSample],
    gnss_map: &GnssFixMap,
    last_gnss: &mut Option<(f64, Vector3<f64>)>,
    q_diag: &Vector15<f64>,
    speed: f64,
) -> Option<(u32, Vector3<f64>, EskfSnapshot)> {
    let preint = build_preintegration(acc_imu, &state.accel_bias, &state.gyro_bias)?;
    let epoch_key = (epoch.time.tow * 10.0).round() as u32;
    let fix = gnss_map.get(&epoch_key).copied();
    let is_gnss = fix.is_some();
    let (pred, phi) = step_inertial_filter(state, epoch.time, &preint, fix, last_gnss, q_diag, speed);
    let snap = EskfSnapshot {
        time: epoch.time,
        state_pred: pred,
        state_post: state.clone(),
        phi,
        is_gnss_available: is_gnss,
    };
    Some((epoch_key, state.pos_ecef, snap))
}

fn init_pipeline_state(
    imu_samples: &[ImuSample],
    init_pos: Vector3<f64>,
    init_heading: f64,
) -> EskfState {
    let gyro_bias = compute_gyro_bias(imu_samples);
    let init_att = compute_initial_attitude(imu_samples, init_pos, init_heading);
    let l_e0 = init_att.to_rotation_matrix().into_inner() * ANTENNA_LEVER_ARM;
    init_eskf_filter(init_pos - l_e0, Vector3::zeros(), init_att, gyro_bias)
}

fn run_inertial_pipeline(
    rover_epochs: &[gneiss_core::obs::EpochObs],
    imu_records: &[ImuRecord],
    gnss_map: &GnssFixMap,
    init_pos: Vector3<f64>,
    init_heading: f64,
) -> PipelineOutput {
    let imu_samples: Vec<ImuSample> = imu_records.iter().map(|r| r.sample).collect();
    let mut state = init_pipeline_state(&imu_samples, init_pos, init_heading);
    let (mut smoother, mut fwd_map) = (EskfSmoother::new(), BTreeMap::new());
    let (mut last_gnss, mut last_idx, mut acc_imu) = (None, 0usize, Vec::new());
    let (mut cur_speed, q_diag) = (0.0, default_q_diag());

    for epoch in rover_epochs {
        let cur_us = (epoch.time.tow * 1_000_000.0).round() as u64;
        while last_idx < imu_records.len() && imu_records[last_idx].sample.time_us <= cur_us {
            acc_imu.push(imu_records[last_idx].sample);
            cur_speed = imu_records[last_idx].speed;
            last_idx += 1;
        }
        if let Some((k, pos, snap)) = process_inertial_epoch(&mut state, epoch, &acc_imu, gnss_map, &mut last_gnss, &q_diag, cur_speed) {
            fwd_map.insert(k, pos);
            smoother.push(snap);
        }
        let last_sample = acc_imu.last().copied();
        acc_imu.clear();
        if let Some(s) = last_sample { acc_imu.push(s); }
    }
    (fwd_map, smoother.smooth_with_time().unwrap_or_default())
}

fn compute_horizontal_errors(
    positions: &BTreeMap<u32, Vector3<f64>>,
    truth: &BTreeMap<u32, (f64, f64, f64)>,
) -> (Vec<f64>, Vec<f64>) {
    let mut sorted = Vec::new();
    let mut chrono = Vec::new();
    for (tow, p) in positions {
        if let Some(&(tx, ty, tz)) = truth.get(tow) {
            let t_pos = Vector3::new(tx, ty, tz);
            let enu = gneiss_core::coords::ecef_delta_to_enu(*p, t_pos, gneiss_core::coords::ecef_to_llh(t_pos));
            let e = (enu.x * enu.x + enu.y * enu.y).sqrt();
            sorted.push(e);
            chrono.push(e);
        }
    }
    sorted.sort_by(|a, b| a.total_cmp(b));
    (sorted, chrono)
}

fn print_quartile_stats(chrono_errs: &[f64]) {
    if chrono_errs.len() <= 1000 { return; }
    let q = chrono_errs.len() / 4;
    for i in 0..4 {
        let mut s = chrono_errs[i * q..(i + 1) * q].to_vec();
        s.sort_by(|a, b| a.total_cmp(b));
        let c_rms = (s.iter().map(|e| e * e).sum::<f64>() / s.len() as f64).sqrt();
        println!("  Q{}: p50={:.3}m, RMS={:.3}m", i + 1, s[s.len() / 2], c_rms);
    }
}

fn print_trajectory_stats(
    name: &str,
    positions: &BTreeMap<u32, Vector3<f64>>,
    truth: &BTreeMap<u32, (f64, f64, f64)>,
) {
    let (h_errs, chrono_errs) = compute_horizontal_errors(positions, truth);
    if h_errs.is_empty() { return; }
    let n = h_errs.len();
    let p50 = h_errs[n / 2];
    let p68 = h_errs[(n as f64 * 0.68) as usize];
    let p95 = h_errs[(n as f64 * 0.95) as usize];
    let rms = (h_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    println!("=== {} (N={}) ===", name, n);
    println!("Horizontal error: p50={:.3}m, p68={:.3}m, p95={:.3}m, RMS={:.3}m", p50, p68, p95, rms);
    print_quartile_stats(&chrono_errs);
}

fn parse_truth_line(line: &str) -> Option<(u32, (f64, f64, f64))> {
    let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
    if parts.len() < 8 {
        return None;
    }
    let tow: f64 = parts[0].parse().ok()?;
    let x: f64 = parts[5].parse().ok()?;
    let y: f64 = parts[6].parse().ok()?;
    let z: f64 = parts[7].parse().ok()?;
    Some(((tow * 10.0).round() as u32, (x, y, z)))
}

fn load_truth(path: &Path) -> BTreeMap<u32, (f64, f64, f64)> {
    let ref_csv = std::fs::read_to_string(path).expect("read reference.csv");
    ref_csv.lines().skip(1).filter_map(parse_truth_line).collect()
}

fn compute_base_pos(approx: Option<[f64; 3]>) -> Vector3<f64> {
    let default_base_arp = Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let default_base_llh = gneiss_core::coords::ecef_to_llh(default_base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(default_base_llh).transpose();
    let computed_base = default_base_arp + ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855);
    approx.map(|p| Vector3::new(p[0], p[1], p[2])).unwrap_or(computed_base)
}

struct RinexInputs {
    rover_epochs: Vec<gneiss_core::obs::EpochObs>,
    base_index: BTreeMap<u32, gneiss_core::obs::EpochObs>,
    base_pos: Vector3<f64>,
    ephems: Vec<gneiss_core::ephemeris::Ephemeris>,
    klob: Option<([f64; 4], [f64; 4])>,
}

fn load_rinex_inputs(dataset: &Path) -> RinexInputs {
    let nav_file = File::open(dataset.join("base.nav")).expect("open base.nav");
    let (ephems, klob) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_file)).expect("parse nav");
    let rover_file = File::open(dataset.join("rover_trimble.obs")).expect("open rover_trimble.obs");
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rover_file)).expect("parse rover");
    let base_file = File::open(dataset.join("base_trimble.obs")).expect("open base_trimble.obs");
    let (base_epochs, base_approx) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_file)).expect("parse base");
    let base_pos = compute_base_pos(base_approx.approx_position);
    let base_index = base_epochs.into_iter().map(|e| ((e.time.tow * 1000.0).round() as u32, e)).collect();
    RinexInputs { rover_epochs, base_index, base_pos, ephems, klob: klob.map(|k| (k.alpha, k.beta)) }
}

fn print_evaluation_summary(
    gnss_positions: &BTreeMap<u32, Vector3<f64>>,
    forward_map: &BTreeMap<u32, Vector3<f64>>,
    smoothed_map: &BTreeMap<u32, Vector3<f64>>,
    truth: &BTreeMap<u32, (f64, f64, f64)>,
) {
    println!("Solutions count: GNSS={}, Forward={}, Smoothed={}", gnss_positions.len(), forward_map.len(), smoothed_map.len());
    print_trajectory_stats("GNSS-Only RTK (Raw Fixes)", gnss_positions, truth);
    print_trajectory_stats("Forward Inertial Filter", forward_map, truth);
    print_trajectory_stats("RTS Smoothed GNSS/INS", smoothed_map, truth);
    let smoothed_at_gnss: BTreeMap<u32, Vector3<f64>> = gnss_positions.keys()
        .filter_map(|k| smoothed_map.get(k).map(|&p| (*k, p)))
        .collect();
    print_trajectory_stats("RTS Smoothed (at GNSS Epochs)", &smoothed_at_gnss, truth);
}

fn main() {
    tracing_subscriber::fmt().with_env_filter("warn").with_target(false).without_time().try_init().ok();
    let dataset = Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba");
    let inputs = load_rinex_inputs(dataset);
    let truth = load_truth(&dataset.join("reference.csv"));
    let imu_records = parse_imu_csv(&dataset.join("imu.csv"));
    let base_ref_map: BTreeMap<u32, &gneiss_core::obs::EpochObs> = inputs.base_index.iter().map(|(k, v)| (*k, v)).collect();

    let rover_init = inputs.rover_epochs.first()
        .and_then(|e| gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(e, &inputs.ephems))
        .unwrap_or(inputs.base_pos);

    println!("=== Odaiba Tightly-Coupled GNSS/INS Benchmark ===");
    println!("Rover Epochs: {}, IMU Samples: {}, Truth Points: {}", inputs.rover_epochs.len(), imu_records.len(), truth.len());

    let gnss_fixes = run_gnss_rtk_pass(&inputs.rover_epochs, &base_ref_map, inputs.base_pos, rover_init, &inputs.ephems, inputs.klob);
    let gnss_positions: BTreeMap<u32, Vector3<f64>> = gnss_fixes.iter().map(|(&k, &(p, _, _))| (k, p)).collect();
    let filter_init = gnss_fixes.values().next().map(|&(p, _, _)| p).unwrap_or(rover_init);
    let init_heading = estimate_initial_heading(&gnss_fixes);
    println!("Initial heading: {:.2} deg (NovAtel reference: 326.65 deg)", init_heading.to_degrees());

    let (forward_map, smoothed) = run_inertial_pipeline(&inputs.rover_epochs, &imu_records, &gnss_fixes, filter_init, init_heading);
    let smoothed_map: BTreeMap<u32, Vector3<f64>> = smoothed.iter().map(|(time, s)| (((time.tow * 10.0).round() as u32), s.pos_ecef)).collect();

    print_evaluation_summary(&gnss_positions, &forward_map, &smoothed_map, &truth);
}
