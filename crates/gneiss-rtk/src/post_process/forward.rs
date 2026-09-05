//! Pass 2: Forward Processing Filter for Offline RTK/PPK.
//!
//! Propagates state forward from t0 to tend, estimating positions,
//! covariances, and carrier-phase integer ambiguities.

use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;

use crate::post_process::dynamics::{ProcessingDynamics, Q_ACCEL_UNSET_FALLBACK};
use crate::post_process::iekf_pass::{
    configure_iekf, dump_amb_history, estimate_epoch_covariance, find_matched_base, FilteredEpoch,
};
use crate::swfg::config::EngineConfig;
use crate::swfg::engine::SwfgEngine;
use crate::swfg::imu_preintegration::{ImuPreintegration, ImuSample};

/// Run forward RTK/PPK estimation pass across all rover epochs.
#[allow(clippy::too_many_arguments)]
pub fn run_forward_pass(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    base_pos: Option<Vector3<f64>>,
    imu_samples: Option<&[ImuSample]>,
    initial_rover_pos: Option<Vector3<f64>>,
    q_accel: Option<f64>,
    dynamics: ProcessingDynamics,
    widelane_ar: bool,
    tropo_grad: bool,
    sat_upd: Option<std::collections::HashMap<u16, f64>>,
    receiver_pcv: Option<std::sync::Arc<super::ReceiverPcvPair>>,
    enable_glonass: bool,
    precise_orbits: Option<std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>>,
    precise_clocks: Option<std::sync::Arc<gneiss_parsers::rinex_clk::RinexClock>>,
    sinex_bias: Option<std::sync::Arc<gneiss_parsers::sinex_bia::SinexBias>>,
    antex_database: Option<std::sync::Arc<gneiss_parsers::antex::AntexDatabase>>,
) -> Vec<FilteredEpoch> {
    if imu_samples.is_none() && base_pos.is_some() && base_epochs.is_some() {
        if let Some(bp) = base_pos {
            // An explicit q_accel always wins; `None` keeps the legacy
            // 1.0 fallback in BOTH profiles so unset behaviour stays
            // byte-identical. Profiles choose their Q at the options
            // level (eval binaries pass STATIC/KINEMATIC_Q_ACCEL).
            let q_eff = q_accel.unwrap_or(Q_ACCEL_UNSET_FALLBACK);
            let (epochs, _wl, _pw) = run_forward_iekf(ephemerides, rover_epochs, base_epochs.unwrap_or(&[]), bp, initial_rover_pos, q_eff, dynamics, widelane_ar, tropo_grad, sat_upd.clone(), receiver_pcv, enable_glonass);
            return epochs;
        }
    }

    let mut cfg = config.clone();
    if let EngineConfig::Ppp(ref mut c) = cfg {
        c.is_kinematic = dynamics.is_kinematic();
    }

    run_forward_swfg(&cfg, ephemerides, klobuchar, rover_epochs, base_epochs, base_pos, imu_samples, precise_orbits, precise_clocks, sinex_bias, antex_database)
}

#[allow(clippy::too_many_arguments)]
fn run_forward_iekf(
    ephemerides: &[Ephemeris],
    rover_epochs: &[EpochObs],
    base_epochs: &[EpochObs],
    base_pos: Vector3<f64>,
    initial_rover_pos: Option<Vector3<f64>>,
    q_accel: f64,
    dynamics: ProcessingDynamics,
    widelane_ar: bool,
    tropo_grad: bool,
    sat_upd: Option<std::collections::HashMap<u16, f64>>,
    receiver_pcv: Option<std::sync::Arc<super::ReceiverPcvPair>>,
    enable_glonass: bool,
) -> (Vec<FilteredEpoch>, crate::estimators::rtk_iekf::mw::WidelaneTracker, crate::estimators::rtk_iekf::mw::WidelaneTracker) {
    if rover_epochs.is_empty() {
        return (
            Vec::new(),
            Default::default(),
            Default::default(),
        );
    }
    let init_pos = initial_rover_pos.unwrap_or_else(|| {
        compute_initial_position(rover_epochs, ephemerides, base_pos)
    });
    let mut iekf = configure_iekf(
        init_pos, rover_epochs[0].time, q_accel, base_pos, dynamics,
        widelane_ar, tropo_grad, enable_glonass, receiver_pcv,
        rover_epochs, sat_upd,
    );

    let mut results = Vec::with_capacity(rover_epochs.len());

    for epoch in rover_epochs {
        if let Some(base_ep) = find_matched_base(epoch.time.tow, Some(base_epochs)) {
            if let Ok(filtered) = iekf.process_epoch(epoch, base_ep, base_pos, ephemerides) {
                results.push(filtered);
            }
        }
    }
    if iekf.track_ambiguity_keys && !iekf.history.is_empty() {
        dump_amb_history(&iekf, "Forward");
    }
    (results, iekf.wl_tracker.clone(), iekf.pw_tracker.clone())
}

fn compute_initial_position(
    rover_epochs: &[EpochObs],
    ephemerides: &[Ephemeris],
    fallback: Vector3<f64>,
) -> Vector3<f64> {
    for ep in rover_epochs.iter().take(10) {
        let seed = crate::estimators::spp::SppState::new(
            gneiss_core::coords::Coordinate::new(
                fallback,
                gneiss_core::coords::Datum::WGS84,
                gneiss_core::coords::Frame::ECEF,
                ep.time,
            ),
            0.0,
            0.0,
            0.0,
            0.0,
        );
        if let Ok(spp) = crate::estimators::spp::compute_spp(
            ep,
            ephemerides,
            None,
            &crate::estimators::spp::SppConfig::default(),
            Some(&seed),
        ) {
            println!("SPP initialization successful: [{:.2}, {:.2}, {:.2}]", spp.position.vector.x, spp.position.vector.y, spp.position.vector.z);
            return spp.position.vector;
        }
    }
    fallback
}

/// Forward-only pass that also returns the converged wide-lane (MW) and
/// phase-only wide-lane trackers — used by the network UPD pre-pass.
#[allow(clippy::too_many_arguments)]
pub fn run_forward_pass_collecting(
    ephemerides: &[Ephemeris],
    rover_epochs: &[EpochObs],
    base_epochs: &[EpochObs],
    base_pos: Vector3<f64>,
    q_accel: f64,
) -> (
    Vec<FilteredEpoch>,
    crate::estimators::rtk_iekf::mw::WidelaneTracker,
    crate::estimators::rtk_iekf::mw::WidelaneTracker,
) {
    let (epochs, wl, pw) = run_forward_iekf(
        ephemerides,
        rover_epochs,
        base_epochs,
        base_pos,
        None,
        q_accel,
        // UPD pre-pass stays STATIC regardless of the run profile:
        // satellite wide-lane biases are receiver-motion independent,
        // and static Q maximises arc-mean convergence.
        ProcessingDynamics::Static,
        true,
        false,
        None,
        None,
        // GLONASS pairs are excluded from MW wide-lane arcs by design
        // (code inter-channel biases don't cancel between receivers,
        // docs/NETWORK_RTK_NEXT_STEPS.md "GLONASS: FDMA plumbing landed")
        // -- this pre-pass computes exactly those arcs, independent of
        // the caller's own enable_glonass setting for DD formation.
        false,
    );
    (epochs, wl, pw)
}

#[allow(clippy::too_many_arguments)]
fn configure_swfg_engine(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
    precise_orbits: Option<std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>>,
    precise_clocks: Option<std::sync::Arc<gneiss_parsers::rinex_clk::RinexClock>>,
    sinex_bias: Option<std::sync::Arc<gneiss_parsers::sinex_bia::SinexBias>>,
    antex_database: Option<std::sync::Arc<gneiss_parsers::antex::AntexDatabase>>,
) -> SwfgEngine {
    let mut engine = SwfgEngine::new(config, ephemerides.to_vec());
    if let Some((alpha, beta)) = klobuchar { engine.set_klobuchar(alpha, beta); }
    if let Some(orbits) = precise_orbits { engine.set_precise_products(orbits, precise_clocks); }
    if let Some(bias) = sinex_bias { engine.set_sinex_bias(bias); }
    if let Some(antex) = antex_database { engine.set_antex_database(antex); }
    engine
}

#[allow(clippy::too_many_arguments)]
fn run_forward_swfg(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    base_pos: Option<Vector3<f64>>,
    imu_samples: Option<&[ImuSample]>,
    precise_orbits: Option<std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>>,
    precise_clocks: Option<std::sync::Arc<gneiss_parsers::rinex_clk::RinexClock>>,
    sinex_bias: Option<std::sync::Arc<gneiss_parsers::sinex_bia::SinexBias>>,
    antex_database: Option<std::sync::Arc<gneiss_parsers::antex::AntexDatabase>>,
) -> Vec<FilteredEpoch> {
    let mut engine = configure_swfg_engine(
        config, ephemerides, klobuchar, precise_orbits, precise_clocks, sinex_bias, antex_database,
    );
    let mut imu_idx = 0usize;
    let mut results = Vec::with_capacity(rover_epochs.len());
    let mut prev_pos: Option<Vector3<f64>> = None;

    for epoch in rover_epochs {
        let preint = extract_imu_slice(epoch, imu_samples, &mut imu_idx);
        let base_ep = find_matched_base(epoch.time.tow, base_epochs);
        if let Some(filtered) = process_single_fwd(
            &mut engine, epoch, base_ep, base_pos, preint, &mut prev_pos, imu_samples.is_some(),
        ) {
            results.push(filtered);
        }
    }
    results
}

fn process_single_fwd(
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
    // Per-epoch SWFG outcome trace (docs/SOLVER_MODE_MATRIX.md): the
    // rover-only/PPP path silently drops every failed epoch with no log
    // line anywhere, which is how the IfbGlonass orphan-variable crash
    // went unnoticed as "0 epochs processed" instead of a visible error.
    if std::env::var("GNEISS_SWFG_DEBUG").is_ok() {
        match &sol_res {
            Ok(s) => eprintln!("SWFG tow={:.0} n_sat={} err={:?}", epoch.time.tow, s.n_satellites, s.error),
            Err(e) => eprintln!("SWFG tow={:.0} ERR={}", epoch.time.tow, e),
        }
    }
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

/// Helper to extract IMU samples between epochs for preintegration.
fn extract_imu_slice(
    epoch: &EpochObs,
    imu_samples: Option<&[ImuSample]>,
    imu_idx: &mut usize,
) -> Option<ImuPreintegration> {
    let samples = imu_samples?;
    let cur_us = (epoch.time.tow * 1_000_000.0) as u32;
    let mut epoch_samples = Vec::new();

    while *imu_idx < samples.len() && samples[*imu_idx].time_us <= cur_us {
        epoch_samples.push(samples[*imu_idx]);
        *imu_idx += 1;
    }

    if epoch_samples.len() >= 2 && epoch_samples.iter().any(|s| s.gyro.norm() > 1e-6) {
        let mut preint = ImuPreintegration::new();
        preint.integrate(&epoch_samples, &Vector3::zeros(), &Vector3::zeros());
        let dt = preint.dt;
        if dt > 2.0 && (preint.dp.norm() / dt.max(1e-3)) < 0.5 {
            preint.dp = Vector3::zeros();
            preint.dv = Vector3::zeros();
        }
        Some(preint)
    } else {
        None
    }
}
