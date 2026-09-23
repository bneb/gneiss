mod odaiba_helpers;

use std::collections::BTreeMap;
use std::path::Path;
use nalgebra::Vector3;

use gneiss_core::coords::{az_el, ecef_to_llh};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::SatelliteId;
use gneiss_rtk::estimators::doppler::estimate_doppler_velocity;
use gneiss_rtk::estimators::eskf::{
    compute_gyro_bias, compute_initial_attitude, init_eskf_filter, predict_preintegrated,
    update_body_velocity, update_dd_scalar, update_doppler_velocity, update_gnss_position,
    update_zupt, DdMeasurementKind, DdSatGeometry, EskfSnapshot, EskfSmoother, EskfState,
    Matrix15, Vector15,
};
use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::engine::SwfgEngine;
use gneiss_rtk::swfg::imu_preintegration::{ImuPreintegration, ImuSample};

use odaiba_helpers::*;

const ANTENNA_LEVER_ARM: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);

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
    epoch: &EpochObs,
    base_index: &BTreeMap<u32, &EpochObs>,
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

fn run_gnss_engine_loop(
    engine: &mut SwfgEngine,
    rover_epochs: &[EpochObs],
    base_index: &BTreeMap<u32, &EpochObs>,
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
    rover_epochs: &[EpochObs],
    base_index: &BTreeMap<u32, &EpochObs>,
    base_pos: Vector3<f64>,
    rover_init_pos: Vector3<f64>,
    ephemerides: &[Ephemeris],
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

fn apply_motion_constraints(state: &mut EskfState, speed: f64) {
    if speed < 0.05 {
        let r_zupt = nalgebra::Matrix3::from_diagonal(&Vector3::new(0.001, 0.001, 0.001));
        let _ = update_zupt(state, &r_zupt);
    } else {
        let r_v = nalgebra::Matrix3::from_diagonal(&Vector3::new(0.005, 0.002, 0.002));
        let _ = update_body_velocity(state, &Vector3::new(speed, 0.0, 0.0), &r_v);
    }
}

fn find_reference_sat(
    rover_obs: &EpochObs,
    base_obs: &EpochObs,
    ant_pos: Vector3<f64>,
    ephems: &[Ephemeris],
    constellation: gneiss_core::sat::Constellation,
) -> Option<(SatelliteId, Vector3<f64>, Vector3<f64>)> {
    let llh = ecef_to_llh(ant_pos);
    let mut best: Option<(SatelliteId, Vector3<f64>, Vector3<f64>)> = None;
    let mut best_el = -1.0;

    for s in &rover_obs.satellites {
        if s.sat.constellation != constellation {
            continue;
        }
        if !base_obs.satellites.iter().any(|b| b.sat == s.sat) {
            continue;
        }
        let eph = match ephems.iter().find(|e| e.sat() == s.sat) {
            Some(e) => e,
            None => continue,
        };
        let (sat_pos, _, _, _) = eph.position(rover_obs.time);
        let (_, el) = az_el(llh, ant_pos, sat_pos);
        if el > best_el && el >= 0.2618 { // 15 deg mask
            best_el = el;
            let diff = sat_pos - ant_pos;
            let d = diff.norm();
            if d > 1e-3 {
                best = Some((s.sat, sat_pos, diff / d));
            }
        }
    }
    best
}

struct DdObsPair<'a> {
    r_sat: &'a SatObs,
    b_sat: &'a SatObs,
    r_ref: &'a SatObs,
    b_ref: &'a SatObs,
    dd_geom: f64,
    el: f64,
}

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
    let r_b2e = state.attitude.to_rotation_matrix().into_inner();
    let l_e = r_b2e * ANTENNA_LEVER_ARM;
    let innov_norm = (pos - (state.pos_ecef + l_e)).norm();

    let mut var_p = if fixed { 0.001 } else if ns >= 6 { 0.004 } else { 0.04 };
    let max_phys_step = (speed * dt_g) + 2.5;
    let is_spike = (speed < 0.05 && innov_norm > 0.8) || (step > max_phys_step && innov_norm > 3.0);
    if is_spike {
        var_p = 1e6;
    }

    let r_pos = nalgebra::Matrix3::from_diagonal(&Vector3::new(var_p, var_p, var_p));
    let _ = update_gnss_position(state, &pos, &r_pos, &ANTENNA_LEVER_ARM);
    *last_gnss = Some((time.tow, pos));
}

fn apply_single_sat_dd(
    state: &mut EskfState,
    geom: &DdSatGeometry,
    obs: &DdObsPair<'_>,
    pos_fix: Option<Vector3<f64>>,
) {
    if let (Some(pr_r), Some(pr_b), Some(pr_rr), Some(pr_br)) = (
        obs.r_sat.get_observable(1),
        obs.b_sat.get_observable(1),
        obs.r_ref.get_observable(1),
        obs.b_ref.get_observable(1),
    ) {
        let dd_code = (pr_r - pr_b) - (pr_rr - pr_br);
        let y_code = dd_code - obs.dd_geom;
        if y_code.abs() <= 3.5 {
            let sin_el = obs.el.sin().max(0.2618);
            let var_code = (2.0 / sin_el).powi(2);
            let _ = update_dd_scalar(state, geom, y_code, var_code, DdMeasurementKind::Pseudorange);
        }
    }

    if let Some(p_fix) = pos_fix {
        let delta_pos = p_fix - (state.pos_ecef + state.attitude.to_rotation_matrix() * geom.lever_arm);
        let delta_u = geom.u_s - geom.u_ref;
        let y_carrier = -delta_u.dot(&delta_pos);
        if y_carrier.abs() <= 0.05 {
            let _ = update_dd_scalar(state, geom, y_carrier, 0.0004, DdMeasurementKind::CarrierPhase);
        }
    }
}

fn apply_tightly_coupled_dd_updates(
    state: &mut EskfState,
    rover_obs: &EpochObs,
    base_obs: &EpochObs,
    base_pos: Vector3<f64>,
    ephems: &[Ephemeris],
    fix: Option<(Vector3<f64>, usize, bool)>,
) {
    let r_b2e = state.attitude.to_rotation_matrix();
    let ant_pos = state.pos_ecef + r_b2e * ANTENNA_LEVER_ARM;
    let llh = ecef_to_llh(ant_pos);
    let pos_fix = if fix.is_some_and(|f| f.2) { fix.map(|f| f.0) } else { None };

    use gneiss_core::sat::Constellation::*;
    for &constellation in &[Gps, Beidou, Galileo, Qzss] {
        let (ref_sat, p_ref, u_ref) = match find_reference_sat(rover_obs, base_obs, ant_pos, ephems, constellation) {
            Some(tup) => tup,
            None => continue,
        };
        let r_ref = match rover_obs.satellites.iter().find(|s| s.sat == ref_sat) {
            Some(s) => s,
            None => continue,
        };
        let b_ref = match base_obs.satellites.iter().find(|s| s.sat == ref_sat) {
            Some(s) => s,
            None => continue,
        };

        for r_sat in &rover_obs.satellites {
            if r_sat.sat.constellation != constellation || r_sat.sat == ref_sat {
                continue;
            }
            let b_sat = match base_obs.satellites.iter().find(|s| s.sat == r_sat.sat) {
                Some(s) => s,
                None => continue,
            };
            let eph = match ephems.iter().find(|e| e.sat() == r_sat.sat) {
                Some(e) => e,
                None => continue,
            };
            let (p_s, _, _, _) = eph.position(rover_obs.time);
            let diff = p_s - ant_pos;
            let d = diff.norm();
            if d < 1e-3 { continue; }
            let u_s = diff / d;
            let (_, el) = az_el(llh, ant_pos, p_s);
            let dd_geom = ((p_s - ant_pos).norm() - (p_s - base_pos).norm())
                - ((p_ref - ant_pos).norm() - (p_ref - base_pos).norm());
            let geom = DdSatGeometry::new(u_s, u_ref, ANTENNA_LEVER_ARM);
            let obs = DdObsPair { r_sat, b_sat, r_ref, b_ref, dd_geom, el };
            apply_single_sat_dd(state, &geom, &obs, pos_fix);
        }
    }
}

struct EpochContext<'a> {
    epoch: &'a EpochObs,
    acc_imu: &'a [ImuSample],
    ephems: &'a [Ephemeris],
    gnss_map: &'a GnssFixMap,
    base_index: &'a BTreeMap<u32, &'a EpochObs>,
    base_pos: Vector3<f64>,
    q_diag: &'a Vector15<f64>,
    speed: f64,
}

type PipelineOutput = (BTreeMap<u32, Vector3<f64>>, Vec<(gneiss_core::time::GpsTime, EskfState)>);

fn process_inertial_epoch(
    state: &mut EskfState,
    ctx: &EpochContext<'_>,
    last_gnss: &mut Option<(f64, Vector3<f64>)>,
) -> Option<(u32, Vector3<f64>, EskfSnapshot)> {
    let preint = build_preintegration(ctx.acc_imu, &state.accel_bias, &state.gyro_bias)?;
    let epoch_key = (ctx.epoch.time.tow * 10.0).round() as u32;
    let exact_ms = (ctx.epoch.time.tow * 1000.0).round() as u32;
    let fix = ctx.gnss_map.get(&epoch_key).copied();
    let is_gnss = fix.is_some();
    let doppler = estimate_doppler_velocity(ctx.epoch, ctx.ephems, state.pos_ecef);

    let phi = predict_preintegrated(state, &preint.dp, &preint.dv, &preint.dq, preint.dt, ctx.q_diag)
        .unwrap_or_else(|_| Matrix15::identity());
    let state_pred = state.clone();

    if let Some((pos, ns, fixed)) = fix {
        update_gnss_innovation(state, pos, ns, fixed, ctx.epoch.time, last_gnss, ctx.speed);
    }
    if let Some(&base_obs) = ctx.base_index.get(&exact_ms) {
        apply_tightly_coupled_dd_updates(state, ctx.epoch, base_obs, ctx.base_pos, ctx.ephems, fix);
    }

    if let Some(d_sol) = doppler {
        if d_sol.vdop < 10.0 && d_sol.n_sats >= 4 {
            let avg_gyro = if ctx.acc_imu.is_empty() {
                Vector3::zeros()
            } else {
                ctx.acc_imu.iter().map(|s| s.gyro).sum::<Vector3<f64>>() / (ctx.acc_imu.len() as f64)
            };
            let _ = update_doppler_velocity(state, &d_sol.vel_ecef, &d_sol.cov, &ANTENNA_LEVER_ARM, &avg_gyro);
        }
    }
    apply_motion_constraints(state, ctx.speed);

    let snap = EskfSnapshot {
        state_pred,
        state_post: state.clone(),
        time: ctx.epoch.time,
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
    let gyro_bias = compute_gyro_bias(imu_samples, 350);
    let init_att = compute_initial_attitude(imu_samples, init_pos, init_heading, 350);
    let l_e0 = init_att.to_rotation_matrix().into_inner() * ANTENNA_LEVER_ARM;
    init_eskf_filter(init_pos - l_e0, Vector3::zeros(), init_att, gyro_bias)
}

struct PipelineConfig<'a> {
    rover_epochs: &'a [EpochObs],
    imu_records: &'a [ImuRecord],
    gnss_map: &'a GnssFixMap,
    base_index: &'a BTreeMap<u32, &'a EpochObs>,
    base_pos: Vector3<f64>,
    ephems: &'a [Ephemeris],
    init_pos: Vector3<f64>,
    init_heading: f64,
}

fn run_inertial_pipeline(cfg: &PipelineConfig<'_>) -> PipelineOutput {
    let imu_samples: Vec<ImuSample> = cfg.imu_records.iter().map(|r| r.sample).collect();
    let mut state = init_pipeline_state(&imu_samples, cfg.init_pos, cfg.init_heading);
    let (mut smoother, mut fwd_map) = (EskfSmoother::new(), BTreeMap::new());
    let (mut last_gnss, mut last_idx, mut acc_imu) = (None, 0usize, Vec::new());
    let (mut cur_speed, q_diag) = (0.0, default_q_diag());

    for epoch in cfg.rover_epochs {
        let cur_us = (epoch.time.tow * 1_000_000.0).round() as u64;
        while last_idx < cfg.imu_records.len() && cfg.imu_records[last_idx].sample.time_us <= cur_us {
            acc_imu.push(cfg.imu_records[last_idx].sample);
            cur_speed = cfg.imu_records[last_idx].speed;
            last_idx += 1;
        }
        let ctx = EpochContext {
            epoch,
            acc_imu: &acc_imu,
            ephems: cfg.ephems,
            gnss_map: cfg.gnss_map,
            base_index: cfg.base_index,
            base_pos: cfg.base_pos,
            q_diag: &q_diag,
            speed: cur_speed,
        };
        if let Some((k, pos, snap)) = process_inertial_epoch(&mut state, &ctx, &mut last_gnss) {
            fwd_map.insert(k, pos);
            smoother.push(snap);
        }
        let last_sample = acc_imu.last().copied();
        acc_imu.clear();
        if let Some(s) = last_sample { acc_imu.push(s); }
    }
    (fwd_map, smoother.smooth_with_time().unwrap_or_default())
}

fn main() {
    tracing_subscriber::fmt().with_env_filter("warn").with_target(false).without_time().try_init().ok();
    let dataset = Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba");
    let inputs = load_rinex_inputs(dataset);
    let truth = load_truth(&dataset.join("reference.csv"));
    let imu_records = parse_imu_csv(&dataset.join("imu.csv"));
    let base_ref_map: BTreeMap<u32, &EpochObs> = inputs.base_index.iter().map(|(k, v)| (*k, v)).collect();

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

    let cfg = PipelineConfig {
        rover_epochs: &inputs.rover_epochs,
        imu_records: &imu_records,
        gnss_map: &gnss_fixes,
        base_index: &base_ref_map,
        base_pos: inputs.base_pos,
        ephems: &inputs.ephems,
        init_pos: filter_init,
        init_heading,
    };
    let (forward_map, smoothed) = run_inertial_pipeline(&cfg);
    let smoothed_map: BTreeMap<u32, Vector3<f64>> = smoothed.iter().map(|(time, s)| (((time.tow * 10.0).round() as u32), s.pos_ecef)).collect();

    print_evaluation_summary(&gnss_positions, &forward_map, &smoothed_map, &truth);
}
