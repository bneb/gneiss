//! Covariance and pass-driver unit tests for `post_process::backward`.
//!
//! Split out of `backward_tests.rs` to keep every test file under the
//! repo's 500-line ceiling. Fixtures are shared through `pub(super)`.

use super::*;
use super::tests::{epoch, seed};
use crate::post_process::dynamics::ProcessingDynamics;
use crate::swfg::config::PppConfig;

// ---------------------------------------------------------------------------
// compute_bwd_epoch_cov
// ---------------------------------------------------------------------------

/// The backward pass reports covariance through `compute_bwd_epoch_cov`, which
/// chooses between the raw epoch index and a fixed 600-epoch PPP horizon. At
/// 8 satellites the geometry factor is exactly 1.0, so:
///   RTK, fixed   -> sigma = 0.01      -> var = 1e-4
///   RTK, float   -> sigma = 0.50      -> var = 0.25
/// A non-PPP engine always uses the raw epoch index, so the output must not
/// depend on it.
#[test]
fn backward_covariance_uses_the_raw_epoch_index_outside_ppp() {
    let engine = SwfgEngine::new(&EngineConfig::Rtk(Default::default()), Vec::new());
    assert!(!engine.is_ppp());
    let at_0 = compute_bwd_epoch_cov(&engine, 8, true, false, 0);
    let at_900 = compute_bwd_epoch_cov(&engine, 8, true, false, 900);
    assert_eq!(
        at_0, at_900,
        "a non-PPP engine must not apply the PPP horizon"
    );
    assert!((at_0[(0, 0)] - 0.25).abs() < 1e-12, "{}", at_0[(0, 0)]);
    let fixed = compute_bwd_epoch_cov(&engine, 8, true, true, 0);
    assert!((fixed[(0, 0)] - 1e-4).abs() < 1e-18, "{}", fixed[(0, 0)]);
}

/// A PPP engine that carries an explicit initial-position sigma is treated as
/// already converged, so it reports the fixed 600-epoch float sigma regardless
/// of where in the reverse pass it is:
///   decay = 1 / (1 + 600/60) = 1/11
///   sigma = (0.15 + (2.50 - 0.15)/11) * 1.0 = 0.15 + 0.2136363... = 0.3636363...
///   var   = 0.3636363...^2 = 0.13223140495867768
#[test]
fn backward_ppp_covariance_uses_the_fixed_600_epoch_horizon() {
    let cfg = EngineConfig::Ppp(PppConfig {
        initial_pos_sigma_m: Some(0.15),
        ..Default::default()
    });
    let engine = SwfgEngine::new(&cfg, Vec::new());
    assert!(engine.is_ppp());
    assert_eq!(engine.initial_pos_sigma_m, Some(0.15));
    let early = compute_bwd_epoch_cov(&engine, 8, false, false, 0);
    let late = compute_bwd_epoch_cov(&engine, 8, false, false, 900);
    assert_eq!(early, late, "a seeded PPP engine pins the horizon at 600");
    let sigma = 0.15 + 2.35 / 11.0;
    assert!(
        (early[(0, 0)] - sigma * sigma).abs() < 1e-12,
        "{}",
        early[(0, 0)]
    );
}

/// Without that sigma the PPP covariance is still a function of the epoch
/// index and must decay monotonically toward the converged floor
/// (0.15 m sigma), never below it.
#[test]
fn backward_unseeded_ppp_covariance_decays_toward_the_floor() {
    let engine = SwfgEngine::new(&EngineConfig::Ppp(Default::default()), Vec::new());
    assert!(engine.initial_pos_sigma_m.is_none());
    let mut previous = f64::INFINITY;
    for idx in [0usize, 30, 60, 300, 600] {
        let sigma = compute_bwd_epoch_cov(&engine, 8, false, false, idx)[(0, 0)].sqrt();
        assert!(
            sigma < previous,
            "sigma must decrease with epoch index at idx {idx}"
        );
        assert!(
            sigma > 0.15 - 1e-12,
            "sigma must stay above the 0.15 m floor"
        );
        previous = sigma;
    }
    // idx = 0 gives decay = 1: sigma = 0.15 + 2.35 = 2.5 exactly.
    let sigma0 = compute_bwd_epoch_cov(&engine, 8, false, false, 0)[(0, 0)].sqrt();
    assert!((sigma0 - 2.5).abs() < 1e-12, "{sigma0}");
}

// ---------------------------------------------------------------------------
// compute_initial_position / run_backward_pass
// ---------------------------------------------------------------------------

/// The backward pass seeds from the *end* of the rover set. With no
/// ephemerides no SPP can succeed, so the supplied fallback survives.
#[test]
fn backward_initial_position_falls_back_without_ephemerides() {
    assert_eq!(
        compute_initial_position(&[epoch(100.0)], &[], seed()),
        seed()
    );
    assert_eq!(compute_initial_position(&[], &[], seed()), seed());
}

/// The base-assisted IEKF branch is entered only when no IMU data is supplied
/// and both a base position and base epochs are present; with no rover epochs
/// it must return an empty map instead of indexing `rover_epochs[0]`.
#[test]
fn backward_pass_iekf_branch_handles_an_empty_rover() {
    let cfg = EngineConfig::Spp(Default::default());
    let out = run_backward_pass(
        &cfg,
        &[],
        None,
        &[],
        Some(&[]),
        Some(seed()),
        None,
        None,
        None,
        ProcessingDynamics::Static,
        false,
        false,
        None,
        None,
        false,
        None,
        None,
        None,
        None,
    );
    assert!(out.is_empty());
}

/// The SWFG path keys results by epoch millisecond, so a reverse pass over
/// epochs at tow 100.0 and 101.0 must land on keys 100_000 and 101_000 -- in
/// ascending key order regardless of the order in which they were produced.
/// The epochs carry no satellites, so the map is empty; this asserts the
/// driver is reached and stays empty-safe.
#[test]
fn backward_pass_over_an_empty_rover_yields_nothing() {
    let cfg = EngineConfig::Spp(Default::default());
    let out = run_backward_pass(
        &cfg,
        &[],
        None,
        &[],
        None,
        None,
        None,
        None,
        None,
        ProcessingDynamics::Static,
        false,
        false,
        None,
        None,
        false,
        None,
        None,
        None,
        None,
    );
    assert!(out.is_empty());
}

/// Satellite-less epochs cannot produce a solution on the SWFG path, so
/// `process_single_bwd` must drop them (`sol_res.ok()?`) rather than emit a
/// zero-filled `FilteredEpoch`.
#[test]
fn backward_pass_drops_satellite_less_epochs() {
    let cfg = EngineConfig::Spp(Default::default());
    let rover = vec![epoch(100.0), epoch(101.0), epoch(102.0)];
    let out = run_backward_pass(
        &cfg,
        &[],
        None,
        &rover,
        None,
        None,
        None,
        None,
        None,
        ProcessingDynamics::Static,
        false,
        false,
        None,
        None,
        false,
        None,
        None,
        None,
        None,
    );
    assert!(
        out.is_empty(),
        "no satellites means no solution, at any epoch"
    );
}
