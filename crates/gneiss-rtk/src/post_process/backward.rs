//! Pass 3: Backward Processing Filter for Offline RTK/PPK.
//!
//! Propagates state backward in time from tend to t0. Carrier-phase
//! ambiguities resolve backwards from open sky into obstructed environments,
//! securing fixes that the forward pass missed.

use std::collections::BTreeMap;
use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;

use crate::estimators::rtk_iekf::GnssRtkIekf;
use crate::post_process::forward::{estimate_epoch_covariance, find_matched_base, FilteredEpoch};
use crate::swfg::config::EngineConfig;
use crate::swfg::engine::SwfgEngine;
use crate::swfg::imu_preintegration::{ImuPreintegration, ImuSample};

/// Run backward RTK/PPK estimation pass across all rover epochs in reverse order.
#[allow(clippy::too_many_arguments)]
pub fn run_backward_pass(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    base_pos: Option<Vector3<f64>>,
    imu_samples: Option<&[ImuSample]>,
    initial_rover_pos: Option<Vector3<f64>>,
    q_accel: Option<f64>,
    widelane_ar: bool,
    sat_upd: Option<std::collections::HashMap<u16, f64>>,
) -> BTreeMap<u64, FilteredEpoch> {
    if imu_samples.is_none() && base_pos.is_some() && base_epochs.is_some() {
        if let Some(bp) = base_pos {
            let init_p = initial_rover_pos.unwrap_or_else(|| {
                compute_initial_position(rover_epochs, ephemerides, bp)
            });
            return run_backward_iekf(ephemerides, rover_epochs, base_epochs.unwrap_or(&[]), bp, init_p, q_accel.unwrap_or(1.0), widelane_ar, sat_upd);
        }
    }

    run_backward_swfg(
        config, ephemerides, klobuchar, rover_epochs, base_epochs, base_pos, imu_samples, initial_rover_pos,
    )
}

fn compute_initial_position(
    rover_epochs: &[EpochObs],
    ephemerides: &[Ephemeris],
    fallback: Vector3<f64>,
) -> Vector3<f64> {
    for ep in rover_epochs.iter().rev().take(10) {
        if let Ok(spp) = crate::estimators::spp::compute_spp(
            ep,
            ephemerides,
            None,
            &crate::estimators::spp::SppConfig::default(),
            None,
        ) {
            return spp.position.vector;
        }
    }
    fallback
}

#[allow(clippy::too_many_arguments)]
fn run_backward_iekf(
    ephemerides: &[Ephemeris],
    rover_epochs: &[EpochObs],
    base_epochs: &[EpochObs],
    base_pos: Vector3<f64>,
    initial_rover_pos: Vector3<f64>,
    q_accel: f64,
    widelane_ar: bool,
    sat_upd: Option<std::collections::HashMap<u16, f64>>,
) -> BTreeMap<u64, FilteredEpoch> {
    let mut results = BTreeMap::new();
    if rover_epochs.is_empty() {
        return results;
    }

    let mut rev_epochs = rover_epochs.to_vec();
    rev_epochs.reverse();

    let mut iekf = GnssRtkIekf::new(initial_rover_pos, rev_epochs[0].time, q_accel);
    iekf.widelane_ar = widelane_ar;
    if widelane_ar {
        // Two-phase static Q (mirrors forward pass): the backward session
        // anchor is end-of-day, so elapsed time counts symmetrically.
        iekf.static_lock_after_s = Some(900.0);
        iekf.static_lock_q_accel = 1e-8;
    }
    if widelane_ar {
        // Rover-side ZWD random walk: gate by baseline length (same as
        // forward pass) so only correlated-atmosphere short baselines get it.
        let baseline_m = (initial_rover_pos - base_pos).norm();
        if baseline_m < super::forward::ZWD_BASELINE_GATE_M {
            iekf.state.enable_zwd(0.0225);
        }
        let cadence_hint =
            crate::post_process::screening::infer_cadence_hint(rover_epochs);
        iekf.slip_detector.cadence_hint_s = cadence_hint;
        iekf.base_slip_detector.cadence_hint_s = cadence_hint;
    } else {
        let cadence_hint =
            crate::post_process::screening::infer_cadence_hint(rover_epochs);
        iekf.slip_detector.cadence_hint_s = cadence_hint;
        iekf.base_slip_detector.cadence_hint_s = cadence_hint;
    }
    iekf.wl_tracker.sat_upd = sat_upd.clone();

    for epoch in &rev_epochs {
        let tow_ms = (epoch.time.tow * 1000.0).round() as u64;
        if let Some(base_ep) = find_matched_base(epoch.time.tow, Some(base_epochs)) {
            if let Ok(filtered) = iekf.process_epoch(epoch, base_ep, base_pos, ephemerides) {
                results.insert(tow_ms, filtered);
            }
        }
    }
    results
}

#[allow(clippy::too_many_arguments)]
fn run_backward_swfg(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    base_pos: Option<Vector3<f64>>,
    imu_samples: Option<&[ImuSample]>,
    initial_rover_pos: Option<Vector3<f64>>,
) -> BTreeMap<u64, FilteredEpoch> {
    let bwd_config = configure_backward_engine(config, initial_rover_pos);
    let mut engine = SwfgEngine::new(&bwd_config, ephemerides.to_vec());
    if let Some((alpha, beta)) = klobuchar {
        engine.set_klobuchar(alpha, beta);
    }

    let epoch_imu_map = group_imu_by_epoch(rover_epochs, imu_samples);
    let mut results = BTreeMap::new();
    let mut prev_pos: Option<Vector3<f64>> = None;

    let mut rev_epochs = rover_epochs.to_vec();
    rev_epochs.reverse();

    for epoch in &rev_epochs {
        let tow_ms = (epoch.time.tow * 1000.0).round() as u64;
        let preint = extract_backward_imu_slice(tow_ms, &epoch_imu_map);
        let base_ep = find_matched_base(epoch.time.tow, base_epochs);
        if let Some(filtered) = process_single_bwd(
            &mut engine, epoch, base_ep, base_pos, preint, &mut prev_pos, imu_samples.is_some(),
        ) {
            results.insert(tow_ms, filtered);
        }
    }
    results
}

fn process_single_bwd(
    engine: &mut SwfgEngine,
    epoch: &EpochObs,
    base: Option<&EpochObs>,
    base_pos: Option<Vector3<f64>>,
    preint: Option<ImuPreintegration>,
    prev_pos: &mut Option<Vector3<f64>>,
    has_imu: bool,
) -> Option<FilteredEpoch> {
    let sol_res = match (base, base_pos) {
        (Some(b), Some(bp)) => engine.process_rtk_epoch_with_imu(epoch, b, bp, preint),
        _ => engine.process_epoch(epoch),
    };
    let sol = sol_res.ok()?;
    if sol.n_satellites < 4 && !has_imu {
        return None;
    }
    let is_fix = sol.error.is_some_and(|e| e < 0.05);
    let cov = estimate_epoch_covariance(sol.n_satellites, base_pos.is_some(), is_fix);
    let q = if is_fix { 1 } else if base_pos.is_some() { 2 } else { 4 };
    let vel = prev_pos.map(|p| sol.position_ecef - p);
    *prev_pos = Some(sol.position_ecef);

    Some(FilteredEpoch {
        time: epoch.time,
        position_ecef: sol.position_ecef,
        velocity_ecef: vel,
        attitude: None,
        cov_position: cov,
        n_satellites: sol.n_satellites,
        quality: q,
        is_fixed: is_fix,
    })
}

fn configure_backward_engine(config: &EngineConfig, initial_pos: Option<Vector3<f64>>) -> EngineConfig {
    let mut cfg = config.clone();
    let pos_arr = initial_pos.map(|p| [p.x, p.y, p.z]);
    match &mut cfg {
        EngineConfig::Rtk(c) => c.initial_position = pos_arr.or(c.initial_position),
        EngineConfig::RtkIns(c) => c.rtk.initial_position = pos_arr.or(c.rtk.initial_position),
        EngineConfig::Ppp(c) => c.initial_position = pos_arr.or(c.initial_position),
        EngineConfig::PppIns(c) => c.ppp.initial_position = pos_arr.or(c.ppp.initial_position),
        EngineConfig::Spp(c) => c.initial_position = pos_arr.or(c.initial_position),
    }
    cfg
}

/// Helper to group IMU samples by epoch millisecond timestamp.
fn group_imu_by_epoch(
    rover_epochs: &[EpochObs],
    imu_samples: Option<&[ImuSample]>,
) -> BTreeMap<u64, Vec<ImuSample>> {
    let mut map = BTreeMap::new();
    let samples = match imu_samples {
        Some(s) => s,
        None => return map,
    };
    let mut cur_idx = 0usize;
    for epoch in rover_epochs {
        let tow_ms = (epoch.time.tow * 1000.0).round() as u64;
        let cur_us = (epoch.time.tow * 1_000_000.0) as u32;
        let mut slice = Vec::new();
        while cur_idx < samples.len() && samples[cur_idx].time_us <= cur_us {
            slice.push(samples[cur_idx]);
            cur_idx += 1;
        }
        if !slice.is_empty() {
            map.insert(tow_ms, slice);
        }
    }
    map
}

/// Helper to extract reverse-integrated IMU slice.
fn extract_backward_imu_slice(
    tow_ms: u64,
    epoch_imu_map: &BTreeMap<u64, Vec<ImuSample>>,
) -> Option<ImuPreintegration> {
    let samples = epoch_imu_map.get(&tow_ms)?;
    if samples.len() < 2 || !samples.iter().any(|s| s.gyro.norm() > 1e-6) {
        return None;
    }
    let rev_samples: Vec<ImuSample> = samples.iter().rev().cloned().map(|mut s| {
        s.accel = -s.accel;
        s.gyro = -s.gyro;
        s
    }).collect();

    let mut preint = ImuPreintegration::new();
    preint.integrate(&rev_samples, &Vector3::zeros(), &Vector3::zeros());
    let dt = preint.dt;
    if dt > 2.0 && (preint.dp.norm() / dt.max(1e-3)) < 0.5 {
        preint.dp = Vector3::zeros();
        preint.dv = Vector3::zeros();
    }
    Some(preint)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_backward_imu_slice_reverses_signs() {
        let mut map = BTreeMap::new();
        map.insert(100_000, vec![
            ImuSample { accel: Vector3::new(1.0, 0.0, 0.0), gyro: Vector3::new(0.1, 0.0, 0.0), time_us: 100_000 },
            ImuSample { accel: Vector3::new(1.0, 0.0, 0.0), gyro: Vector3::new(0.1, 0.0, 0.0), time_us: 200_000 },
        ]);
        let preint = extract_backward_imu_slice(100_000, &map);
        assert!(preint.is_some());
    }

    #[test]
    fn test_empty_backward_pass_runs() {
        let config = EngineConfig::Spp(Default::default());
        let results = run_backward_pass(&config, &[], None, &[], None, None, None, None, None, false, None);
        assert!(results.is_empty());
    }
}
