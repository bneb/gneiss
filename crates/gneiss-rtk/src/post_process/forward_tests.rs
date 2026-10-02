//! Unit tests for `post_process::forward` (pass 2).
//!
//! Companion file: `forward.rs` is over the 300-line threshold, so its tests
//! live here and are pulled in with `#[cfg(test)] #[path = ...] mod tests`.
//!
//! The preintegration expectations below are the closed-form midpoint
//! trapezoid that `ImuPreintegration::integrate` implements for
//! `dq = identity`, `ba = bg = 0`:
//!   dv_{k+1} = dv_k + a * dt
//!   dp_{k+1} = dp_k + dv_{k+1} * dt + 0.5 * a * dt^2

use super::*;
use crate::post_process::dynamics::ProcessingDynamics;

fn epoch(tow: f64) -> EpochObs {
    EpochObs {
        time: gneiss_core::time::GpsTime::new(2000, tow),
        satellites: Vec::new(),
    }
}

fn imu(t_us: u64, accel: Vector3<f64>, gyro: Vector3<f64>) -> ImuSample {
    ImuSample {
        accel,
        gyro,
        time_us: t_us,
    }
}

fn preint_of(samples: &[ImuSample]) -> ImuPreintegration {
    let mut p = ImuPreintegration::new();
    p.integrate(samples, &Vector3::zeros(), &Vector3::zeros());
    p
}

// ---------------------------------------------------------------------------
// extract_imu_slice
// ---------------------------------------------------------------------------

/// With no IMU stream there is nothing to preintegrate.
#[test]
fn imu_slice_is_none_without_imu_data() {
    let mut idx = 0usize;
    assert!(extract_imu_slice(&epoch(100.0), None, &mut idx).is_none());
    assert_eq!(idx, 0, "a missing stream must not advance the cursor");
}

/// `time_us <= cur_us` drains the buffer, so the cursor is a true monotonic
/// drain point: the same epoch asked twice sees an empty slice the second
/// time and yields `None` (fewer than two samples).
#[test]
fn imu_slice_drains_the_cursor_monotonically() {
    let samples = vec![
        imu(
            50_000_000,
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
        ),
        imu(
            60_000_000,
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
        ),
    ];
    let mut idx = 0usize;
    assert!(extract_imu_slice(&epoch(100.0), Some(&samples), &mut idx).is_some());
    assert_eq!(idx, 2, "both samples precede tow 100 s and are consumed");
    assert!(extract_imu_slice(&epoch(100.0), Some(&samples), &mut idx).is_none());
    assert_eq!(idx, 2);
}

/// Fewer than two samples cannot span an interval, and a window whose gyro is
/// everywhere below 1e-6 rad/s carries no rotational information, so both
/// windows are rejected by the same guard.
#[test]
fn imu_slice_rejects_short_and_zero_rate_windows() {
    let one = vec![imu(
        50_000_000,
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.1, 0.0, 0.0),
    )];
    assert!(extract_imu_slice(&epoch(100.0), Some(&one), &mut 0).is_none());
    let still = vec![
        imu(50_000_000, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros()),
        imu(60_000_000, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros()),
    ];
    assert!(extract_imu_slice(&epoch(100.0), Some(&still), &mut 0).is_none());
}

/// Three samples 50 ms apart with a constant body-X specific force of 1 m/s^2
/// (dq stays identity because the gyro is about X, parallel to `a`).
/// `integrate` advances `dp` from the *pre-update* velocity and then updates
/// `dv`, i.e. semi-implicit:
///   step 1: dp += dv(=0)*0.05 + 0.5*1*0.05^2 = 0.00125 ; dv = 0.05
///   step 2: dp += dv(=0.05)*0.05 + 0.5*1*0.05^2 = 0.0025 + 0.00125
///   total: dt = 0.1 s, dv = (0.1, 0, 0), dp = 0.00125 + 0.00375 = (0.005, 0, 0)
/// `dt = 0.1 <= 2.0`, so the long-window zeroing guard is not taken.
#[test]
fn imu_slice_integrates_constant_specific_force() {
    let samples = (0..3)
        .map(|i| {
            imu(
                50_000_000 + i * 50_000,
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect::<Vec<_>>();
    let p = extract_imu_slice(&epoch(100.0), Some(&samples), &mut 0).expect("preintegration runs");
    assert!((p.dt - 0.1).abs() < 1e-12, "dt = {}", p.dt);
    assert!(
        (p.dv - Vector3::new(0.1, 0.0, 0.0)).norm() < 1e-12,
        "dv = {:?}",
        p.dv.as_slice()
    );
    assert!(
        (p.dp - Vector3::new(0.005, 0.0, 0.0)).norm() < 1e-12,
        "dp = {:?}",
        p.dp.as_slice()
    );
}

/// The zeroing guard fires when `dt > 2.0` while the implied mean velocity
/// `|dp| / dt` stays below 0.5 m/s -- i.e. a long window with no real motion.
/// 30 samples at 100 ms give dt = 29 * 0.1 = 2.9 s > 2.0; with `a = 0.01 m/s^2`
/// the integrated |dp| stays around 4.5e-2 m, so |dp|/dt ~ 1.6e-2 < 0.5.
#[test]
fn imu_slice_zeroes_a_long_but_motionless_window() {
    let n = 30usize;
    let samples = (0..n)
        .map(|i| {
            imu(
                50_000_000 + i as u64 * 100_000,
                Vector3::new(0.01, 0.0, 0.0),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect::<Vec<_>>();
    let raw = preint_of(&samples);
    assert!(
        raw.dt > 2.0,
        "dt = {} must exceed the 2 s threshold",
        raw.dt
    );
    assert!(
        raw.dp.norm() / raw.dt < 0.5,
        "implied speed must be sub-0.5 m/s"
    );
    let p = extract_imu_slice(&epoch(100.0), Some(&samples), &mut 0).expect("preintegration runs");
    assert!(p.dt > 2.0, "dt itself is preserved");
    assert_eq!(p.dp, Vector3::zeros());
    assert_eq!(p.dv, Vector3::zeros());
}

/// The other half of the guard: a long window that *does* accumulate velocity
/// keeps its deltas. 30 samples with `a = 1 m/s^2` give dv = 2.9 m/s and
/// |dp| / dt ~ 1.55 m/s, well above the 0.5 m/s trigger.
#[test]
fn imu_slice_keeps_a_long_window_with_real_acceleration() {
    let samples = (0..30)
        .map(|i| {
            imu(
                50_000_000 + i as u64 * 100_000,
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect::<Vec<_>>();
    let raw = preint_of(&samples);
    assert!(raw.dt > 2.0 && raw.dp.norm() / raw.dt > 0.5);
    let p = extract_imu_slice(&epoch(100.0), Some(&samples), &mut 0).expect("preintegration runs");
    assert!(
        (p.dv - raw.dv).norm() < 1e-12,
        "dv must pass through untouched"
    );
    assert!(
        (p.dp - raw.dp).norm() < 1e-12,
        "dp must pass through untouched"
    );
}

// ---------------------------------------------------------------------------
// compute_initial_position
// ---------------------------------------------------------------------------

/// With no ephemerides every SPP attempt fails, so the documented fallback is
/// returned verbatim. This is the branch that keeps a caller-supplied base
/// position from being replaced by garbage when seeding is impossible.
#[test]
fn initial_position_falls_back_when_spp_cannot_run() {
    let fallback = Vector3::new(-3_961_907.0, 3_351_057.0, 3_694_313.0);
    assert_eq!(
        compute_initial_position(&[epoch(100.0)], &[], fallback),
        fallback
    );
    // Even with an empty epoch list the fallback survives the loop.
    assert_eq!(compute_initial_position(&[], &[], fallback), fallback);
}

// ---------------------------------------------------------------------------
// configure_swfg_engine
// ---------------------------------------------------------------------------

/// `configure_swfg_engine` must hand the caller's `EngineConfig` through to
/// the engine unchanged: an RTK profile must not be silently promoted to PPP
/// (which would change the covariance model) and vice versa. The Klobuchar
/// and precise-product setters are exercised on the same call so the whole
/// configuration body runs.
#[test]
fn swfg_engine_configuration_preserves_the_profile() {
    let cfg = EngineConfig::Rtk(Default::default());
    let engine = configure_swfg_engine(
        &cfg,
        &[],
        Some(([1.0; 4], [2.0; 4])),
        None,
        None,
        None,
        None,
    );
    assert!(!engine.is_ppp(), "an RTK profile must not report as PPP");
    assert!(
        !engine.enable_glonass && !engine.enable_galileo,
        "RTK never enables GLONASS/Galileo"
    );

    let ppp = EngineConfig::Ppp(Default::default());
    let engine = configure_swfg_engine(&ppp, &[], None, None, None, None, None);
    assert!(engine.is_ppp(), "a PPP profile must report as PPP");
}

// ---------------------------------------------------------------------------
// run_forward_pass / run_forward_pass_collecting
// ---------------------------------------------------------------------------

/// An empty rover set produces no epochs on the SWFG path, and the
/// bidirectional-independent signature makes the same call shape as the
/// production caller.
#[test]
fn forward_pass_over_an_empty_rover_yields_nothing() {
    let cfg = EngineConfig::Spp(Default::default());
    let out = run_forward_pass(
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

/// The base-assisted IEKF branch is selected only when IMU data is absent and
/// both a base position and base epochs are present. With no rover epochs it
/// must return an empty vector rather than panicking on `rover_epochs[0]`.
#[test]
fn forward_pass_iekf_branch_handles_an_empty_rover() {
    let cfg = EngineConfig::Spp(Default::default());
    let base_pos = Vector3::new(-3_961_907.0, 3_351_057.0, 3_694_313.0);
    let out = run_forward_pass(
        &cfg,
        &[],
        None,
        &[],
        Some(&[]),
        Some(base_pos),
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
        "empty rover epochs short-circuit the IEKF loop"
    );
}

/// Supplying IMU data deliberately *disables* the base-assisted IEKF branch
/// (`imu_samples.is_none()` guards it), so the run is routed to the SWFG path
/// instead. With a satellite-less epoch that path yields nothing, which is how
/// this test observes the routing decision without a real GNSS dataset.
#[test]
fn forward_pass_routes_to_swfg_when_imu_is_present() {
    let cfg = EngineConfig::Spp(Default::default());
    let base_pos = Vector3::new(-3_961_907.0, 3_351_057.0, 3_694_313.0);
    let imu: Vec<ImuSample> = (0..10)
        .map(|i| {
            imu(
                50_000_000 + i * 50_000,
                Vector3::new(0.0, 0.0, 9.8),
                Vector3::new(0.0, 0.0, 0.0),
            )
        })
        .collect();
    let rover = vec![epoch(100.0)];
    let out = run_forward_pass(
        &cfg,
        &[],
        None,
        &rover,
        Some(&[]),
        Some(base_pos),
        Some(&imu),
        None,
        None,
        ProcessingDynamics::Kinematic,
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
        "a satellite-less epoch cannot produce a solution"
    );
}

/// The UPD pre-pass entry point shares `run_forward_iekf` with the ordinary
/// forward pass, so it must be empty-safe too. The wide-lane trackers are
/// returned by value in their default (empty) state; they are not observable
/// beyond their public `sat_upd` field, which the pre-pass leaves unset.
#[test]
fn forward_pass_collecting_is_empty_safe() {
    let base_pos = Vector3::new(-3_961_907.0, 3_351_057.0, 3_694_313.0);
    let (epochs, wl, pw) = run_forward_pass_collecting(&[], &[], &[], base_pos, 1e-6);
    assert!(epochs.is_empty());
    assert!(wl.sat_upd.is_none());
    assert!(pw.sat_upd.is_none());
}

// ---------------------------------------------------------------------------
// Quality mapping shared by both passes
// ---------------------------------------------------------------------------

/// `estimate_swfg_epoch_covariance` is the single covariance model used by
/// this pass. At 8 satellites the geometry factor is exactly
/// `8 / max(8, 4) = 1.0`, so the horizontal variance is `sigma^2` directly:
/// fixed -> (0.01)^2 = 1e-4, RTK float -> (0.50)^2 = 0.25, and the vertical
/// axis is inflated by 2.25 -> 0.25 * 2.25 = 0.5625 for the float case.
#[test]
fn forward_covariance_model_is_exact_at_eight_satellites() {
    let fixed = estimate_swfg_epoch_covariance(8, true, true, false, 0);
    assert!((fixed[(0, 0)] - 1e-4).abs() < 1e-18, "{}", fixed[(0, 0)]);
    assert!((fixed[(2, 2)] - 2.25e-4).abs() < 1e-18);
    let float = estimate_swfg_epoch_covariance(8, true, false, false, 0);
    assert!((float[(0, 0)] - 0.25).abs() < 1e-12, "{}", float[(0, 0)]);
    assert!((float[(2, 2)] - 0.5625).abs() < 1e-12);
    assert_eq!(fixed, fixed.transpose(), "covariance must be symmetric");
    assert!(
        fixed.try_inverse().is_some(),
        "covariance must be invertible"
    );
    // More satellites can only tighten the geometry factor, never loosen it.
    let dense = estimate_swfg_epoch_covariance(32, true, false, false, 0);
    assert!(
        dense[(0, 0)] < float[(0, 0)],
        "sigma must shrink with more satellites"
    );
    // The `max(4)` floor means 1..4 satellites all collapse to the same value.
    assert_eq!(
        estimate_swfg_epoch_covariance(1, true, false, false, 0),
        estimate_swfg_epoch_covariance(4, true, false, false, 0)
    );
}
