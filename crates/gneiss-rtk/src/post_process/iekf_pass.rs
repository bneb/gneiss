//! Shared `GnssRtkIekf` pass infrastructure for both directional passes.
//!
//! Extracted out of `forward.rs` because `backward.rs` depended on it via
//! `crate::post_process::forward::{configure_iekf, find_matched_base, ...}`
//! -- a layering smell (one directional pass reaching into the other for
//! logic that belongs to neither) as much as a file-size one. Everything
//! here is genuinely direction-agnostic: constructing and configuring a
//! `GnssRtkIekf`, matching base epochs, estimating output covariance, and
//! dumping ambiguity history.

use nalgebra::{Matrix3, Vector3};

use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;

use crate::estimators::rtk_iekf::GnssRtkIekf;
use crate::post_process::dynamics::ProcessingDynamics;

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
/// The slip-detector cadence hint and `wl_tracker.sat_upd` were also found
/// drifted (backward set them unconditionally; forward only when
/// `widelane_ar` was on) and are now folded in here too, unconditionally:
/// - `cadence_hint_s` reflects the *data stream's* sampling rate
///   (screening.rs's gap-detection threshold) and has no principled
///   dependence on `widelane_ar` at all -- forward's gating was the actual
///   bug, not a deliberate restriction.
/// - `wl_tracker.sat_upd` is proven inert whenever `widelane_ar` is off:
///   its only reader, `WidelaneTracker::fixed_widelane`, is reachable
///   solely through `far_matches_widelanes` and `resolve_cascade`, both
///   called only from `if self.widelane_ar` in `process_epoch`. Setting it
///   when `widelane_ar` is false cannot change any output.
///
/// This does change `eval_qinertia_ppk.rs`'s walkthrough output for the
/// forward pass specifically (it runs `widelane_ar: false`, so forward
/// previously never received a cadence hint even on slow-cadence
/// streams) -- see docs/NETWORK_RTK_NEXT_STEPS.md for the measured delta.
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
    rover_epochs: &[EpochObs],
    sat_upd: Option<std::collections::HashMap<u16, f64>>,
) -> GnssRtkIekf {
    let mut iekf = GnssRtkIekf::new(init_pos, start_time, q_accel);
    // Profile-gated robust weighting: 1.0 (Static) is bit-identical legacy.
    iekf.robust_innov_scale = dynamics.innovation_gate_scale();
    apply_env_overrides(&mut iekf);
    iekf.widelane_ar = widelane_ar;
    apply_cadence_and_upd(&mut iekf, rover_epochs, sat_upd);
    // Independent of widelane_ar and baseline length -- previously nested
    // inside the ZWD baseline gate below (env-var read), which meant a
    // >25km baseline never got GLONASS regardless of the caller's
    // request, and diverged from the backward pass's own (different)
    // gating. See PostProcessOptions::enable_glonass.
    iekf.enable_glonass = enable_glonass;
    if let Some(pair) = receiver_pcv {
        iekf.receiver_pcv = Some((pair.rover.clone(), pair.base.clone()));
    }
    configure_static_lock(&mut iekf, dynamics, widelane_ar);
    configure_atmosphere_and_precise_products(&mut iekf, init_pos, base_pos, widelane_ar, tropo_grad);
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

/// Debug/experimental env-var overrides: elevation mask and per-satellite
/// / global ionosphere state estimation. All default off.
fn apply_env_overrides(iekf: &mut GnssRtkIekf) {
    if let Ok(deg) = std::env::var("GNEISS_ELEV_DEG") {
        if let Ok(rad) = deg.parse::<f64>() {
            iekf.min_elevation_rad = rad.to_radians();
        }
    }
    iekf.state.iono_enabled = std::env::var("GNEISS_IONO_STATES").as_deref() == Ok("1");
    iekf.state.sat_iono_enabled = std::env::var("GNEISS_SAT_IONO").as_deref() == Ok("1");
}

/// Unconditional: `cadence_hint_s` reflects the data stream's own sampling
/// rate, not an AR setting, and `wl_tracker.sat_upd` is proven inert
/// whenever `widelane_ar` is off (see `configure_iekf`'s doc comment) --
/// so setting both the same way regardless of `widelane_ar` is safe.
fn apply_cadence_and_upd(
    iekf: &mut GnssRtkIekf,
    rover_epochs: &[EpochObs],
    sat_upd: Option<std::collections::HashMap<u16, f64>>,
) {
    let cadence_hint = crate::post_process::screening::infer_cadence_hint(rover_epochs);
    iekf.slip_detector.cadence_hint_s = cadence_hint;
    iekf.base_slip_detector.cadence_hint_s = cadence_hint;
    iekf.wl_tracker.sat_upd = sat_upd;
}

/// Two-phase static Q: converge loosely, then lock the monument. Kinematic
/// rovers have no monument to lock; the loose phase Q stays for the whole
/// session.
fn configure_static_lock(iekf: &mut GnssRtkIekf, dynamics: ProcessingDynamics, widelane_ar: bool) {
    if widelane_ar && !dynamics.is_kinematic() {
        iekf.static_lock_after_s = Some(900.0);
        iekf.static_lock_q_accel = 1e-8;
    }
}

/// Rover-side ZWD random walk helps short baselines (atmosphere correlated)
/// but hurts long baselines (>25 km) where rover/base wet delay decouples,
/// so gate it -- and the tropo gradients and precise orbit/clock products
/// that only make sense alongside it -- by baseline length.
fn configure_atmosphere_and_precise_products(
    iekf: &mut GnssRtkIekf,
    init_pos: Vector3<f64>,
    base_pos: Vector3<f64>,
    widelane_ar: bool,
    tropo_grad: bool,
) {
    let baseline_m = (init_pos - base_pos).norm();
    if widelane_ar && baseline_m < ZWD_BASELINE_GATE_M {
        iekf.state.enable_zwd(0.0225); // ~15 cm zenith wet init uncertainty
        // Experimental tropo gradients: opt-in via env while the
        // OHLN interaction is unresolved (v_p95 -6mm pooled, but
        // OHLN h_p95 degrades when unconditional).
        if tropo_grad {
            iekf.state.enable_gradients(crate::estimators::rtk_iekf::update::GRAD_INIT_VAR_M2);
        }
        load_precise_orbits(iekf);
        load_precise_clocks(iekf);
    }
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

/// Helper to find closest base epoch within synchronous window (0.1s).
pub(crate) fn find_matched_base(tow: f64, base_epochs: Option<&[EpochObs]>) -> Option<&EpochObs> {
    let epochs = base_epochs?;
    epochs.iter()
        .filter(|b| (b.time.tow - tow).abs() < 0.1)
        .min_by(|a, b| {
            (a.time.tow - tow).abs().total_cmp(&(b.time.tow - tow).abs())
        })
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
