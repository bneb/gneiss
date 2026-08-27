//! Pass 2: Forward Processing Filter for Offline RTK/PPK.
//!
//! Propagates state forward from t0 to tend, estimating positions,
//! covariances, and carrier-phase integer ambiguities.

use nalgebra::{Matrix3, Vector3};

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;

use crate::estimators::rtk_iekf::GnssRtkIekf;
use crate::post_process::dynamics::{ProcessingDynamics, Q_ACCEL_UNSET_FALLBACK};
use crate::swfg::config::EngineConfig;
use crate::swfg::engine::SwfgEngine;
use crate::swfg::imu_preintegration::{ImuPreintegration, ImuSample};

/// Baselines shorter than this have correlated wet delay between stations,
/// so a single rover-side ZWD state captures the residual. Above it, the
/// rover/base wet delay decouples and the extra state degrades fix rates
/// (measured: SLAC 49.7 km dropped 68% -> 47.5% when ungated).
pub(crate) const ZWD_BASELINE_GATE_M: f64 = 25_000.0;

/// Filtered epoch output from a single directional pass.
#[derive(Debug, Clone)]
pub struct FilteredEpoch {
    pub time: GpsTime,
    pub position_ecef: Vector3<f64>,
    pub velocity_ecef: Option<Vector3<f64>>,
    pub attitude: Option<nalgebra::UnitQuaternion<f64>>,
    pub cov_position: Matrix3<f64>,
    pub n_satellites: usize,
    pub quality: u8,
    pub is_fixed: bool,
}

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

    run_forward_swfg(config, ephemerides, klobuchar, rover_epochs, base_epochs, base_pos, imu_samples)
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
    );
    if widelane_ar {
        let cadence_hint =
            crate::post_process::screening::infer_cadence_hint(rover_epochs);
        iekf.slip_detector.cadence_hint_s = cadence_hint;
        iekf.base_slip_detector.cadence_hint_s = cadence_hint;
        iekf.wl_tracker.sat_upd = sat_upd.clone();
    } else {
        iekf.wl_tracker.sat_upd = None;
    }

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

/// Shared `GnssRtkIekf` construction and configuration for both directional
/// passes. Extracted after finding five settings -- elevation mask, iono
/// states, per-satellite iono, precise orbits, precise clocks -- that were
/// wired into the forward pass only, via bare env-var reads inside what
/// used to be forward-only setup code. They silently never reached the
/// backward pass at all: the same independently-drifting-duplicate-setup
/// shape as the `enable_glonass` and `GNEISS_PCV`/`GNEISS_RECV_PCV` bugs
/// fixed earlier (docs/NETWORK_RTK_NEXT_STEPS.md). All five are read from
/// unset-by-default env vars, so this refactor is behavior-preserving
/// whenever none of them are set.
///
/// Deliberately NOT unified here: the slip-detector cadence hint and
/// `wl_tracker.sat_upd`'s gating on `widelane_ar`. Both differ between
/// the two passes today (backward sets them unconditionally; forward
/// only when `widelane_ar` is on), but unifying them would change
/// `eval_qinertia_ppk.rs`'s bit-identical walkthrough output (it runs
/// both passes with `widelane_ar: false`), which needs its own
/// deliberate round with reference regeneration, not a side effect of
/// this one. See docs/NETWORK_RTK_NEXT_STEPS.md.
#[allow(clippy::too_many_arguments)]
pub(crate) fn configure_iekf(
    init_pos: Vector3<f64>,
    start_time: GpsTime,
    q_accel: f64,
    base_pos: Vector3<f64>,
    dynamics: ProcessingDynamics,
    widelane_ar: bool,
    tropo_grad: bool,
    enable_glonass: bool,
    receiver_pcv: Option<std::sync::Arc<super::ReceiverPcvPair>>,
) -> GnssRtkIekf {
    let mut iekf = GnssRtkIekf::new(init_pos, start_time, q_accel);
    // Profile-gated robust weighting: 1.0 (Static) is bit-identical legacy.
    iekf.robust_innov_scale = dynamics.innovation_gate_scale();
    if let Ok(deg) = std::env::var("GNEISS_ELEV_DEG") {
        if let Ok(rad) = deg.parse::<f64>() {
            iekf.min_elevation_rad = rad.to_radians();
        }
    }
    iekf.state.iono_enabled = std::env::var("GNEISS_IONO_STATES").as_deref() == Ok("1");
    iekf.state.sat_iono_enabled =
        std::env::var("GNEISS_SAT_IONO").as_deref() == Ok("1");
    iekf.widelane_ar = widelane_ar;
    // Independent of widelane_ar and baseline length -- previously nested
    // inside the ZWD baseline gate below (env-var read), which meant a
    // >25km baseline never got GLONASS regardless of the caller's
    // request, and diverged from the backward pass's own (different)
    // gating. See PostProcessOptions::enable_glonass.
    iekf.enable_glonass = enable_glonass;
    if let Some(pair) = receiver_pcv {
        iekf.receiver_pcv = Some((pair.rover.clone(), pair.base.clone()));
    }
    if widelane_ar && !dynamics.is_kinematic() {
        // Two-phase static Q: converge loosely, then lock the monument.
        // Kinematic rovers have no monument to lock; the loose phase Q
        // stays for the whole session.
        iekf.static_lock_after_s = Some(900.0);
        iekf.static_lock_q_accel = 1e-8;
    }
    // Rover-side ZWD random walk helps short baselines (atmosphere correlated)
    // but hurts long baselines (>25 km) where rover/base wet delay decouples.
    // Gate by baseline length so only correlated-atmosphere cases get the state.
    let baseline_m = (init_pos - base_pos).norm();
    if widelane_ar && baseline_m < ZWD_BASELINE_GATE_M {
        iekf.state.enable_zwd(0.0225); // ~15 cm zenith wet init uncertainty
        // Experimental tropo gradients: opt-in via env while the
        // OHLN interaction is unresolved (v_p95 -6mm pooled, but
        // OHLN h_p95 degrades when unconditional).
        if tropo_grad {
            iekf.state.enable_gradients(crate::estimators::rtk_iekf::update::GRAD_INIT_VAR_M2);
        }
        load_precise_orbits(&mut iekf);
        load_precise_clocks(&mut iekf);
    }
    // Debug-dump toggle: previously also required baseline < 25km in the
    // forward pass only (accidental -- it was just physically nested
    // inside the ZWD-gate block above, not a principled restriction) and
    // required nothing from the backward pass beyond widelane_ar. Now
    // matches backward's (correct) gating in both passes.
    if widelane_ar {
        iekf.track_ambiguity_keys = std::env::var("GNEISS_AMB_DUMP").is_ok();
    }
    iekf
}

/// Precise orbits: loaded when `GNEISS_SP3` points to a valid SP3 file.
fn load_precise_orbits(iekf: &mut GnssRtkIekf) {
    let Ok(sp3_path) = std::env::var("GNEISS_SP3") else { return };
    let result = std::fs::File::open(&sp3_path)
        .map_err(|e| format!("open failed: {}", e))
        .and_then(|f| gneiss_parsers::sp3::parse_sp3(std::io::BufReader::new(f)));
    match result {
        Ok(epochs) => {
            let store = gneiss_parsers::precise_orbit::PreciseOrbit::new(epochs);
            println!("PRECISE: {} sats from {}", store.len(), sp3_path);
            iekf.precise_orbits = Some(std::sync::Arc::new(store));
        }
        Err(e) => eprintln!("SP3: {}", e),
    }
}

/// Precise clock products (RINEX CLK), used with precise orbits.
fn load_precise_clocks(iekf: &mut GnssRtkIekf) {
    let Ok(clk_path) = std::env::var("GNEISS_CLK") else { return };
    match std::fs::read_to_string(&clk_path) {
        Ok(content) => {
            let rc = gneiss_parsers::rinex_clk::RinexClock::parse(&content);
            println!("PRECISE-CLK: {} satellites tracked", rc.satellites.len());
            iekf.precise_clocks = Some(std::sync::Arc::new(rc));
        }
        Err(e) => eprintln!("CLK: {}", e),
    }
}

/// Write per-key float DD ambiguity trajectories from engine history.
pub(crate) fn dump_amb_history(iekf: &GnssRtkIekf, label: &str) {
    let snaps: Vec<_> = iekf.history.iter().filter(|s| !s.amb_keys.is_empty()).collect();
    if snaps.is_empty() { return; }
    let Ok(dir) = std::env::var("GNEISS_AMB_DUMP_DIR") else { return };
    let _ = std::fs::create_dir_all(&dir);
    // stable key union across all snapshots (keys enter/exit as sats rise/set)
    let mut all_keys: Vec<crate::estimators::rtk_iekf::state::DoubleDiffKey> = Vec::new();
    for s in &iekf.history {
        for k in &s.amb_keys {
            if !all_keys.contains(k) { all_keys.push(*k); }
        }
    }
    let path = format!("{dir}/amb_{label}.csv");
    let Ok(mut f) = std::fs::File::create(&path) else { return };
    use std::io::Write;
    let _ = write!(f, "tow");
    for k in &all_keys {
        let _ = write!(f, ",{}_{}_{}_b{}", k.constellation_id, k.sat, k.ref_sat, k.freq_band);
    }
    let _ = writeln!(f);
    for s in &iekf.history {
        let _ = write!(f, "{:.0}", s.time.tow);
        let extra = if iekf.state.zwd_enabled { 1 } else { 0 }
            + if iekf.state.grad_enabled { 2 } else { 0 };
        let offset = 6 + extra;
        for k in &all_keys {
            match s.amb_keys.iter().position(|kk| kk == k) {
                Some(local_idx) => {
                    let idx = offset + local_idx;
                    if idx < s.x_post.len() {
                        let _ = write!(f, ",{:.6}", s.x_post[idx]);
                    } else {
                        let _ = write!(f, ",");
                    }
                }
                None => { let _ = write!(f, ","); }
            }
        }
        let _ = writeln!(f);
    }
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

fn run_forward_swfg(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    klobuchar: Option<([f64; 4], [f64; 4])>,
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    base_pos: Option<Vector3<f64>>,
    imu_samples: Option<&[ImuSample]>,
) -> Vec<FilteredEpoch> {
    let mut engine = SwfgEngine::new(config, ephemerides.to_vec());
    if let Some((alpha, beta)) = klobuchar {
        engine.set_klobuchar(alpha, beta);
    }

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

/// Helper to find closest base epoch within synchronous window (0.1s).
pub(crate) fn find_matched_base(tow: f64, base_epochs: Option<&[EpochObs]>) -> Option<&EpochObs> {
    let epochs = base_epochs?;
    epochs.iter()
        .filter(|b| (b.time.tow - tow).abs() < 0.1)
        .min_by(|a, b| {
            (a.time.tow - tow).abs().total_cmp(&(b.time.tow - tow).abs())
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

/// Helper to compute representative covariance matrix for position.
pub(crate) fn estimate_epoch_covariance(n_sats: usize, is_rtk: bool, is_fixed: bool) -> Matrix3<f64> {
    let geom_factor = (8.0 / (n_sats.max(4) as f64)).max(0.5);
    let sigma = if is_fixed {
        0.01 * geom_factor // 1cm fixed RTK
    } else if is_rtk {
        0.50 * geom_factor // 50cm float RTK
    } else {
        2.5 * geom_factor // 2.5m SPP
    };
    let var = sigma * sigma;
    Matrix3::new(
        var, 0.0, 0.0,
        0.0, var, 0.0,
        0.0, 0.0, var * 2.25,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_epoch_covariance() {
        let cov_fixed = estimate_epoch_covariance(8, true, true);
        let cov_float = estimate_epoch_covariance(8, true, false);
        assert!(cov_fixed[(0, 0)] < cov_float[(0, 0)]);
    }

    #[test]
    fn test_find_matched_base() {
        let time = GpsTime::new(2000, 100.05);
        let ep = EpochObs { time, satellites: Vec::new() };
        let base_list = vec![ep];

        let matched = find_matched_base(100.06, Some(&base_list));
        assert!(matched.is_some());
        assert_eq!(matched.unwrap().time.tow, 100.05);
    }
}
