//! Pipeline-stage unit tests for `post_process::calibration`.
//!
//! Companion to `calibration_tests.rs`: the calibration module has more
//! surface than one sub-500-line test file can hold, so the trajectory /
//! parameter / report orchestration is split out here. Fixtures are shared
//! with the sibling module through `pub(super)` helpers.

use super::tests::{assert_vec, empty_result, epoch_at, imu, origin};
use super::*;
use nalgebra::UnitQuaternion;

/// One satellite-less rover epoch: enough to get past `execute_post_process`'s
/// empty-dataset guard, and nothing else -- no trajectory ever comes out, so
/// no calibration parameter can move.
fn synthetic_rover() -> Vec<EpochObs> {
    vec![EpochObs {
        time: gneiss_core::time::GpsTime::new(2000, 100.0),
        satellites: Vec::new(),
    }]
}

/// `synthetic_rover` through a default `PostProcessOptions` + calibration
/// budget, for the tests that only care about the report.
fn run_synthetic(rover: &[EpochObs], calib: &MultiPassCalibrationOptions, imu: Option<&[ImuSample]>) -> MultiPassCalibrationReport {
    execute_calibrated_post_process(
        &EngineConfig::Spp(Default::default()),
        &[],
        rover,
        None,
        imu,
        &PostProcessOptions::default(),
        calib,
    )
    .expect("a non-empty rover set runs")
    .1
}

fn fixed_traj_at(tow: f64, n: usize) -> Vec<SmoothedEpoch> {
    let disp = Vector3::new(-0.05, -0.15, 0.30);
    (0..n)
        .map(|_| {
            let mut e = epoch_at(
                tow,
                Some(Vector3::zeros()),
                Some(UnitQuaternion::identity()),
            );
            e.position_ecef = origin() + disp;
            e.quality = 1;
            e
        })
        .collect()
}

// ---------------------------------------------------------------------------
// estimate_body_antenna_offset
// ---------------------------------------------------------------------------

// Fixed epochs whose ECEF displacement from the reference is exactly
// `R_b_e * offset` with `R_b_e` the zero-rpy DCM at the equator, i.e.
//   R_b_e = [[0,0,-1],[0,1,0],[1,0,0]],  offset = (0.30, -0.15, 0.05)
//   R_b_e * offset = (-0.05, -0.15, 0.30).
// The estimator back-rotates with the same `R_b_e`, so
// `R_b_e^T (R_b_e * offset) = offset` for any orthonormal `R_b_e`.
fn reference_map(tow_ms: i64) -> ReferencePointMap {
    let mut refs = ReferencePointMap::new();
    refs.insert(tow_ms, (origin(), [0.0, 0.0, 0.0]));
    refs
}

#[test]
fn body_antenna_offset_recovers_the_injected_offset() {
    let traj = fixed_traj_at(1000.0, 10);
    let got = estimate_body_antenna_offset(&traj, &reference_map(1_000_000))
        .expect("ten fixed epochs clear the count >= 10 gate");
    assert_vec(
        got,
        Vector3::new(0.30, -0.15, 0.05),
        1e-9,
        "antenna body offset",
    );
}

/// The count gate is `count >= 10`, so nine fixed epochs must be rejected.
#[test]
fn body_antenna_offset_needs_ten_fixed_epochs() {
    let traj = fixed_traj_at(1000.0, 9);
    assert!(estimate_body_antenna_offset(&traj, &reference_map(1_000_000)).is_none());
}

/// Only `quality == 1` epochs contribute; float epochs (quality 2/4) are
/// skipped, so ten float epochs still yield nothing.
#[test]
fn body_antenna_offset_ignores_float_epochs() {
    let mut traj = fixed_traj_at(1000.0, 12);
    for e in traj.iter_mut() {
        e.quality = 2;
    }
    assert!(estimate_body_antenna_offset(&traj, &reference_map(1_000_000)).is_none());
}

/// The reference lookup is a +-150 ms window. An epoch 500 ms away from the
/// only reference point is out of range, so it must not be counted.
#[test]
fn body_antenna_offset_respects_the_150ms_window() {
    let traj = fixed_traj_at(1000.5, 10);
    assert!(estimate_body_antenna_offset(&traj, &reference_map(1_000_000)).is_none());
    let in_window = fixed_traj_at(1000.1, 10);
    assert!(estimate_body_antenna_offset(&in_window, &reference_map(1_000_000)).is_some());
}

// ---------------------------------------------------------------------------
// apply_calibration_to_trajectory
// ---------------------------------------------------------------------------

/// With an identity attitude the shift is exactly `R * offset_body`, subtracted
/// from the position: 100 - 1 = 99, 200 - 0.5 = 199.5, 300 - (-0.2) = 300.2.
#[test]
fn apply_calibration_shifts_positions_by_the_body_offset() {
    let mut traj = vec![SmoothedEpoch {
        time: gneiss_core::time::GpsTime::new(2000, 100.0),
        position_ecef: Vector3::new(100.0, 200.0, 300.0),
        velocity_ecef: Some(Vector3::new(10.0, 0.0, 0.0)),
        attitude: Some(UnitQuaternion::identity()),
        cov_position: Matrix3::identity(),
        std_east: 0.01,
        std_north: 0.01,
        std_up: 0.02,
        separation_3d: 0.01,
        quality: 1,
        n_satellites: 8,
    }];
    apply_calibration_to_trajectory(&mut traj, Vector3::new(1.0, 0.5, -0.2));
    assert!((traj[0].position_ecef.x - 99.0).abs() < 1e-4);
    assert!((traj[0].position_ecef.y - 199.5).abs() < 1e-4);
    assert!((traj[0].position_ecef.z - 300.2).abs() < 1e-4);
}

/// The guard is `offset_body.norm() < 1e-5`; a sub-micrometre offset must be
/// a no-op rather than a numerical perturbation.
#[test]
fn apply_calibration_skips_negligible_offsets() {
    let before = origin();
    let mut traj = vec![epoch_at(
        100.0,
        Some(Vector3::new(10.0, 0.0, 0.0)),
        Some(UnitQuaternion::identity()),
    )];
    apply_calibration_to_trajectory(&mut traj, Vector3::new(1e-9, 0.0, 0.0));
    assert_vec(traj[0].position_ecef, before, 0.0, "position untouched");
    // Exactly 1e-5 is not `< 1e-5`, so the shift is applied.
    let mut traj = vec![epoch_at(
        100.0,
        Some(Vector3::new(10.0, 0.0, 0.0)),
        Some(UnitQuaternion::identity()),
    )];
    apply_calibration_to_trajectory(&mut traj, Vector3::new(1e-5, 0.0, 0.0));
    assert!((traj[0].position_ecef.x - (6_378_137.0 - 1e-5)).abs() < 1e-8);
}

/// Without an attitude and without a fast-enough velocity there is no
/// body-to-ECEF rotation, so the epoch must be left alone rather than
/// falling back to some default frame.
#[test]
fn apply_calibration_leaves_epochs_without_a_rotation_alone() {
    let mut traj = vec![epoch_at(100.0, None, None)];
    apply_calibration_to_trajectory(&mut traj, Vector3::new(1.0, 0.0, 0.0));
    assert_vec(traj[0].position_ecef, origin(), 0.0, "position untouched");
}

/// Velocity-only fallback: driving due East at 10 m/s at the equator gives
/// R_b_e = [[0,0,-1],[1,0,0],[0,-1,0]], so `R * (1, 0, 0)` = (0, 1, 0) and
/// the position moves by exactly -1 m in ECEF Y.
#[test]
fn apply_calibration_uses_velocity_heading_without_attitude() {
    let mut traj = vec![epoch_at(100.0, Some(Vector3::new(0.0, 10.0, 0.0)), None)];
    apply_calibration_to_trajectory(&mut traj, Vector3::new(1.0, 0.0, 0.0));
    assert_vec(
        traj[0].position_ecef,
        origin() - Vector3::new(0.0, 1.0, 0.0),
        1e-9,
        "shifted",
    );
}

// ---------------------------------------------------------------------------
// check_convergence
// ---------------------------------------------------------------------------

/// A 3-4-0 right triangle: the accel-bias change has norm exactly
/// sqrt(0.003^2 + 0.004^2) = sqrt(9e-6 + 16e-6) = sqrt(25e-6) = 0.005, and the
/// default tolerance is `accel_bias_tol_mps2 = 0.005` with a `<=` comparison,
/// so this is exactly on the boundary and must converge.
#[test]
fn convergence_uses_inclusive_tolerances() {
    let crit = CalibrationConvergenceCriteria::default();
    let prev = CalibrationParameters {
        imu_accel_bias: Some(Vector3::zeros()),
        ..Default::default()
    };
    let on_boundary = CalibrationParameters {
        imu_accel_bias: Some(Vector3::new(0.003, 0.004, 0.0)),
        ..Default::default()
    };
    let (conv, d_arm, d_accel, d_gyro) = check_convergence(&prev, &on_boundary, &crit);
    assert!((d_accel - 0.005).abs() < 1e-12, "3-4-5 triangle: {d_accel}");
    assert!((d_arm - 0.0).abs() <= f64::EPSILON && (d_gyro - 0.0).abs() <= f64::EPSILON);
    assert!(conv, "0.005 <= 0.005 must converge");

    // Add 1e-3 in z: norm = sqrt(25e-6 + 1e-6) = sqrt(26e-6) = 0.0050990 > 0.005.
    let past = CalibrationParameters {
        imu_accel_bias: Some(Vector3::new(0.003, 0.004, 0.001)),
        ..Default::default()
    };
    let (conv, _, d_accel, _) = check_convergence(&prev, &past, &crit);
    let want = 26.0_f64.sqrt() * 1e-3; // sqrt(0.003^2 + 0.004^2 + 0.001^2)
    assert!((d_accel - want).abs() < 1e-12, "{d_accel}");
    assert!(!conv, "0.005099 > 0.005 must not converge");

    // A genuinely small lever-arm step: (0.002, 0.001, -0.001) has norm
    // sqrt(4 + 1 + 1) e-3 = sqrt(6) e-3 = 0.0024495 < 0.005.
    let arm_a = CalibrationParameters {
        lever_arm_body: Some(Vector3::new(0.50, 0.20, -0.10)),
        ..Default::default()
    };
    let arm_b = CalibrationParameters {
        lever_arm_body: Some(Vector3::new(0.502, 0.201, -0.101)),
        ..Default::default()
    };
    let (conv, d_arm, _, _) = check_convergence(&arm_a, &arm_b, &crit);
    let want = 6.0_f64.sqrt() * 1e-3; // sqrt(0.002^2 + 0.001^2 + 0.001^2)
    assert!((d_arm - want).abs() < 1e-12, "{d_arm}");
    assert!(conv && d_arm < 0.005);
}

/// A parameter appearing on only one side of the comparison contributes a
/// zero delta rather than being treated as an infinite change: a calibration
/// that first *discovers* a lever arm on pass 2 must not be reported as a
/// huge step.
#[test]
fn convergence_treats_a_newly_appeared_parameter_as_zero_delta() {
    let crit = CalibrationConvergenceCriteria::default();
    let prev = CalibrationParameters::default();
    let curr = CalibrationParameters {
        lever_arm_body: Some(Vector3::new(12.0, -34.0, 56.0)),
        imu_gyro_bias: Some(Vector3::new(1.0, 1.0, 1.0)),
        ..Default::default()
    };
    let (conv, d_arm, d_accel, d_gyro) = check_convergence(&prev, &curr, &crit);
    assert_eq!((d_arm, d_accel, d_gyro), (0.0, 0.0, 0.0));
    assert!(conv);
}

/// Gyro deltas use the same rule on the `gyro_bias_tol_radps = 5e-4` gate:
/// 3e-4 / 4e-4 is exactly at it, 3e-4 / 4e-4 / 1e-4 is sqrt(26)*1e-4 past it.
#[test]
fn convergence_checks_the_gyro_bias_tolerance() {
    let crit = CalibrationConvergenceCriteria::default();
    let prev = CalibrationParameters {
        imu_gyro_bias: Some(Vector3::zeros()),
        ..Default::default()
    };
    let at = CalibrationParameters {
        imu_gyro_bias: Some(Vector3::new(0.0003, 0.0004, 0.0)),
        ..Default::default()
    };
    assert!(check_convergence(&prev, &at, &crit).0);
    let past = CalibrationParameters {
        imu_gyro_bias: Some(Vector3::new(0.0003, 0.0004, 0.0001)),
        ..Default::default()
    };
    let (conv, _, _, d_gyro) = check_convergence(&prev, &past, &crit);
    let want = 26.0_f64.sqrt() * 1e-4; // sqrt(0.0003^2 + 0.0004^2 + 0.0001^2)
    assert!((d_gyro - want).abs() < 1e-12, "{d_gyro}");
    assert!(!conv);
}

// ---------------------------------------------------------------------------
// record_iteration / init_calibration_params / analyze_iteration_result
// ---------------------------------------------------------------------------

/// `record_iteration` counts only `quality == 1` epochs and snapshots the
/// parameters by value, so later mutation of the live parameters cannot
/// rewrite history.
#[test]
fn record_iteration_counts_fixed_epochs_and_snapshots_parameters() {
    let mut res = empty_result();
    let mut traj = fixed_traj_at(1000.0, 12);
    for e in traj.iter_mut().skip(6) {
        e.quality = 2;
    }
    res.trajectory = traj;
    let crit = CalibrationConvergenceCriteria::default();
    let prev = CalibrationParameters::default();
    let mut curr = CalibrationParameters {
        lever_arm_body: Some(Vector3::new(1.0, 0.0, 0.0)),
        ..Default::default()
    };
    let (rec, conv) = record_iteration(2, &curr, &prev, &crit, &res);
    assert_eq!(rec.pass_number, 2);
    assert_eq!((rec.fixed_epochs, rec.total_epochs), (6, 12));
    // `prev` has no lever arm at all, so `check_convergence` takes the
    // mismatched-arms arm and reports a zero delta rather than treating the
    // newly discovered 1 m arm as a 1 m step.
    assert_eq!(rec.lever_arm_delta_m, 0.0, "newly appeared lever arm -> zero delta");
    assert_eq!(rec.accel_bias_delta_mps2, 0.0);
    assert_eq!(rec.gyro_bias_delta_radps, 0.0);
    assert!(conv);
    curr.lever_arm_body = Some(Vector3::new(9.0, 0.0, 0.0));
    assert_eq!(curr.lever_arm_body, Some(Vector3::new(9.0, 0.0, 0.0)));
    assert_eq!(
        rec.parameters.lever_arm_body,
        Some(Vector3::new(1.0, 0.0, 0.0))
    );
}

/// `init_calibration_params` only pre-estimates biases when the option is on
/// *and* IMU samples were supplied.
#[test]
fn init_params_respects_the_imu_bias_option_and_sample_availability() {
    let samples: Vec<ImuSample> = (0..20)
        .map(|i| {
            imu(
                i * 50_000,
                Vector3::new(0.0, 0.0, 9.9),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect();
    let seeded_opts = MultiPassCalibrationOptions {
        estimate_imu_biases: true,
        ..Default::default()
    };
    let seeded = init_calibration_params(&seeded_opts, Some(&samples));
    assert_vec(
        seeded.imu_accel_bias.expect("accel bias"),
        Vector3::new(0.0, 0.0, 0.09335),
        1e-12,
        "pre-estimated accel bias",
    );
    assert!(seeded.imu_gyro_bias.is_some());
    assert!(init_calibration_params(&seeded_opts, None)
        .imu_accel_bias
        .is_none());
    let no_bias_opts = MultiPassCalibrationOptions {
        estimate_imu_biases: false,
        ..Default::default()
    };
    assert!(init_calibration_params(&no_bias_opts, Some(&samples))
        .imu_gyro_bias
        .is_none());
}

/// The antenna-offset branch runs `r_b_e^T` (not `r_b_e`) on the ECEF
/// residual, so an offset injected along the zero-rpy DCM comes back with its
/// original body-frame components -- and, with no lever arm yet, it is also
/// adopted as the initial lever arm.
#[test]
fn analyze_iteration_backfills_lever_arm_from_the_antenna_offset() {
    let mut res = empty_result();
    res.trajectory = fixed_traj_at(1000.0, 10);
    let opts = MultiPassCalibrationOptions {
        estimate_lever_arm: false,
        estimate_imu_biases: false,
        reference_points: Some(reference_map(1_000_000)),
        ..Default::default()
    };
    let mut params = CalibrationParameters::default();
    analyze_iteration_result(&res, None, &opts, &mut params);
    let off = params.antenna_body_offset.expect("offset estimated");
    assert_vec(off, Vector3::new(0.30, -0.15, 0.05), 1e-9, "antenna offset");
    assert_eq!(params.lever_arm_body, params.antenna_body_offset);
}

/// With the antenna option off, or with no reference points, nothing is
/// written back: the analysis step must leave the parameters untouched.
#[test]
fn analyze_iteration_writes_nothing_when_options_are_off() {
    let mut res = empty_result();
    res.trajectory = fixed_traj_at(1000.0, 10);
    let mut opts = MultiPassCalibrationOptions {
        estimate_antenna_offset: false,
        ..Default::default()
    };
    let mut params = CalibrationParameters::default();
    analyze_iteration_result(&res, None, &opts, &mut params);
    assert_eq!(params, CalibrationParameters::default());
    opts.estimate_antenna_offset = true;
    analyze_iteration_result(&res, None, &opts, &mut params);
    assert!(
        params.antenna_body_offset.is_none(),
        "no reference points -> no offset"
    );
}

// ---------------------------------------------------------------------------
// execute_calibrated_post_process
// ---------------------------------------------------------------------------

/// An empty rover set is rejected upstream by `execute_post_process`; the
/// calibration driver propagates that `Err` instead of reporting a pass.
#[test]
fn calibrated_post_process_propagates_an_empty_dataset_error() {
    let err = execute_calibrated_post_process(
        &EngineConfig::Spp(Default::default()),
        &[],
        &[],
        None,
        None,
        &PostProcessOptions::default(),
        &MultiPassCalibrationOptions::default(),
    )
    .expect_err("empty rover epochs must fail");
    assert_eq!(err, "No rover epochs provided");
}

/// With a satellite-less dataset nothing is ever calibratable, so the
/// parameters stay at their initial values, every delta is exactly 0, and the
/// loop therefore satisfies the convergence test at pass 2 -- never at pass 1,
/// because the `pass > 1` guard deliberately withholds convergence from the
/// first pass. Expect exactly two recorded passes and `converged == true`.
#[test]
fn calibrated_post_process_converges_after_two_identical_passes() {
    let report = run_synthetic(&synthetic_rover(), &MultiPassCalibrationOptions::default(), None);
    assert_eq!(
        report.passes_executed, 2,
        "converges as soon as pass > 1 agrees"
    );
    assert!(report.converged);
    assert_eq!(report.initial_parameters, report.final_parameters);
    for rec in &report.history {
        assert_eq!(rec.lever_arm_delta_m, 0.0);
        assert_eq!(rec.accel_bias_delta_mps2, 0.0);
        assert_eq!(rec.gyro_bias_delta_radps, 0.0);
        assert_eq!(rec.fixed_epochs, 0);
    }
}

/// `max_iterations` is the hard cap on recorded passes, and the cap is taken
/// with `.max(1)` so a zero budget still executes one pass rather than
/// returning the "No passes executed" error.
#[test]
fn calibrated_post_process_honours_the_iteration_cap() {
    let mut calib = MultiPassCalibrationOptions::default();
    calib.criteria.max_iterations = 7;
    let report = run_synthetic(&synthetic_rover(), &calib, None);
    assert_eq!(report.passes_executed, 2, "converges before the cap of 7");

    calib.criteria.max_iterations = 0;
    let report = run_synthetic(&synthetic_rover(), &calib, None);
    assert_eq!(
        report.passes_executed, 1,
        "max(1) keeps the driver from erroring out"
    );
    assert!(!report.converged, "pass 1 can never set the converged flag");
}

/// Stationary IMU with `estimate_imu_biases` on seeds the initial parameters,
/// and those seeds must reach the report's `initial_parameters` untouched.
#[test]
fn calibrated_post_process_seeds_initial_parameters_from_stationary_imu() {
    let samples: Vec<ImuSample> = (0..40)
        .map(|i| {
            imu(
                i * 50_000,
                Vector3::new(0.0, 0.0, 9.9),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect();
    let mut calib = MultiPassCalibrationOptions::default();
    calib.criteria.max_iterations = 1;
    let report = run_synthetic(&synthetic_rover(), &calib, Some(&samples));
    assert_vec(
        report.initial_parameters.imu_accel_bias.expect("seeded"),
        Vector3::new(0.0, 0.0, 0.09335),
        1e-12,
        "seeded accel bias",
    );
    assert_eq!(
        report.passes_executed, 1,
        "a one-iteration budget records one pass"
    );
    assert!(!report.converged);
}
