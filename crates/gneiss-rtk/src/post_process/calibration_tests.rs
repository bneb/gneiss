//! Unit tests for `post_process::calibration`.
//!
//! Companion file: `calibration.rs` is over the 300-line threshold, so its
//! tests live here and are pulled in with `#[cfg(test)] #[path = ...] mod tests`.
//!
//! Every expected number below is derived by hand in a `//` comment from the
//! WGS84 equator/prime-meridian ENU axes
//!   North = (0, 0, 1),  East = (0, 1, 0),  Up = (1, 0, 0)
//! at `llh = (0, 0, 0)`, for which `R_e_n` (ECEF -> NED) is the exact
//! permutation [[0,0,1],[0,1,0],[-1,0,0]].

use super::*;
use nalgebra::UnitQuaternion;

/// Surface point on the WGS84 equator / prime meridian (h = 0).
pub(super) fn origin() -> Vector3<f64> {
    Vector3::new(6_378_137.0, 0.0, 0.0)
}

pub(super) fn epoch_at(
    tow: f64,
    vel: Option<Vector3<f64>>,
    attitude: Option<UnitQuaternion<f64>>,
) -> SmoothedEpoch {
    SmoothedEpoch {
        time: gneiss_core::time::GpsTime::new(2000, tow),
        position_ecef: origin(),
        velocity_ecef: vel,
        attitude,
        cov_position: Matrix3::identity(),
        std_east: 0.01,
        std_north: 0.01,
        std_up: 0.02,
        separation_3d: 0.01,
        quality: 1,
        n_satellites: 8,
    }
}

pub(super) fn imu(t_us: u64, accel: Vector3<f64>, gyro: Vector3<f64>) -> ImuSample {
    ImuSample {
        accel,
        gyro,
        time_us: t_us,
    }
}

pub(super) fn assert_vec(actual: Vector3<f64>, expected: Vector3<f64>, tol: f64, what: &str) {
    assert!(
        (actual - expected).norm() <= tol,
        "{what}: got {:?}, expected {:?}",
        actual.as_slice(),
        expected.as_slice()
    );
}

fn assert_mat(actual: Matrix3<f64>, expected: Matrix3<f64>, what: &str) {
    assert!(
        (actual - expected).norm() <= 1e-9,
        "{what}: got {:?}, expected {:?}",
        actual.as_slice(),
        expected.as_slice()
    );
}

/// A satellite-less `execute_post_process` run: the engine returns an empty
/// trajectory but `Ok`. Callers overwrite `res.trajectory` with synthetic
/// epochs so the calibration helpers can be driven directly.
pub(super) fn empty_result() -> PostProcessResult {
    let rover = vec![EpochObs {
        time: gneiss_core::time::GpsTime::new(2000, 100.0),
        satellites: Vec::new(),
    }];
    execute_post_process(
        &EngineConfig::Spp(Default::default()),
        &[],
        &rover,
        None,
        None,
        &PostProcessOptions::default(),
    )
    .expect("satellite-less post-process run succeeds")
}

// ---------------------------------------------------------------------------
// heading_pitch_to_dcm
// ---------------------------------------------------------------------------

/// Heading is measured from North toward East: `atan2(E, N)`.
/// Pitch is `atan2(-D, sqrt(N^2+E^2))` -- the denominator is a `sqrt`, so it
/// is always >= 0 and the pitch stays in [-90, +90] deg even when the
/// horizontal speed component is negative (i.e. driving "backwards").
///
/// Heading East: v_ned = (0, v, 0) -> heading = atan2(v, 0) = +90 deg,
///   pitch = atan2(0, v) = 0.
///     R_b_n (cols = forward, right, down) = [[0,-1,0],[1,0,0],[0,0,1]]
///     R_e_n^T                           = [[0,0,-1],[0,1,0],[1,0,0]]
///     product                           = [[0,0,-1],[1,0,0],[0,-1,0]]
#[test]
fn heading_pitch_to_dcm_east_is_exact_permutation() {
    let dcm = heading_pitch_to_dcm(Vector3::new(0.0, 10.0, 0.0), origin())
        .expect("10 m/s is above the 1 m/s gate");
    assert_mat(
        dcm,
        Matrix3::new(0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0),
        "heading east",
    );
    // Column 1 is body-forward and must equal the unit velocity (level driving).
    assert_vec(
        dcm * Vector3::x(),
        Vector3::new(0.0, 1.0, 0.0),
        1e-9,
        "east forward",
    );
    // Column 3 is body-down: ECEF (-1, 0, 0) = -Up.
    assert_vec(
        dcm * Vector3::z(),
        Vector3::new(-1.0, 0.0, 0.0),
        1e-9,
        "east down",
    );
}

/// Heading West is the 180 deg rotation of heading East: cos(180) = -1 and
/// sin(180) = 0, so in `R_b_n` the forward and right rows flip sign while the
/// down column is untouched.
///
///   R_b_n(west) = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
///   R_e_n^T * R_b_n(west) = [[0,0,-1], [-1,0,0], [0,1,0]]
///
/// Comparing with the east case: column 1 (forward) is negated -- the 180 deg
/// flip -- and column 3 (down) is bit-identical, i.e. vertical is unchanged.
#[test]
fn heading_pitch_to_dcm_180_degree_turn_flips_horizontal_only() {
    let east = heading_pitch_to_dcm(Vector3::new(0.0, 10.0, 0.0), origin()).unwrap();
    let west = heading_pitch_to_dcm(Vector3::new(0.0, -10.0, 0.0), origin()).unwrap();
    assert_mat(
        west,
        Matrix3::new(0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        "heading west",
    );
    let fwd_east = east * Vector3::x();
    let fwd_west = west * Vector3::x();
    let down_east = east * Vector3::z();
    let down_west = west * Vector3::z();
    assert_vec(
        fwd_west,
        -fwd_east,
        1e-9,
        "forward flips sign under 180 deg",
    );
    assert_vec(
        down_west,
        down_east,
        1e-9,
        "down is invariant under 180 deg",
    );
}

/// Heading North (velocity = ECEF (1,0,0)*10 = Up at lat=lon=0) gives
/// v_ned = (0, 0, -10): heading = atan2(0, 0) = 0, pitch = atan2(10, 0) = +90.
///   R_b_n = [[0,0,1],[0,1,0],[-1,0,0]]
///   R_e_n^T * R_b_n = identity   (nose straight up at the equator)
#[test]
fn heading_pitch_to_dcm_upward_velocity_is_identity() {
    let dcm = heading_pitch_to_dcm(Vector3::new(10.0, 0.0, 0.0), origin()).unwrap();
    assert_mat(
        dcm,
        Matrix3::identity(),
        "velocity along ECEF +X = local Up",
    );
}

/// A level vehicle's body-forward axis must coincide with its velocity
/// direction, whatever the azimuth and whatever the (non-degenerate) site.
/// Velocity is built in the local ENU frame and rotated back into ECEF, so
/// the expected forward axis is the velocity itself.
#[test]
fn heading_pitch_to_dcm_forward_axis_equals_velocity_direction() {
    let site = gneiss_core::coords::llh_to_ecef(Vector3::new(0.4, -1.1, 250.0));
    let llh = gneiss_core::coords::ecef_to_llh(site);
    let north = Vector3::new(
        llh.x.cos() * llh.y.cos(),
        llh.x.cos() * llh.y.sin(),
        llh.x.sin(),
    );
    let east = Vector3::new(-llh.y.sin(), llh.y.cos(), 0.0);
    for (az, horiz) in [(0.0_f64, north), (1.0, east), (2.7, -north + 0.4 * east)] {
        let v_ecef = horiz.normalize() * 12.0;
        let dcm = heading_pitch_to_dcm(v_ecef, site)
            .unwrap_or_else(|| panic!("azimuth {az} rad: 12 m/s clears the gate"));
        assert_vec(
            dcm * Vector3::x(),
            v_ecef.normalize(),
            1e-9,
            "forward = unit velocity",
        );
        assert!(
            (dcm.determinant() - 1.0).abs() < 1e-12,
            "DCM determinant must be +1"
        );
    }
}

/// The gate is `speed < 1.0`, evaluated on the *norm*, so 0.999 m/s is
/// rejected and exactly 1.0 m/s is accepted (0.999... is strictly below).
#[test]
fn heading_pitch_to_dcm_speed_gate_is_one_metre_per_second() {
    assert!(heading_pitch_to_dcm(Vector3::new(0.0, 0.999, 0.0), origin()).is_none());
    assert!(heading_pitch_to_dcm(Vector3::new(1.0, 0.0, 0.0), origin()).is_some());
    assert!(heading_pitch_to_dcm(Vector3::zeros(), origin()).is_none());
}

// ---------------------------------------------------------------------------
// estimate_static_imu_biases
// ---------------------------------------------------------------------------

/// Exact static bias: 20 identical samples with |accel| = 9.9 m/s^2 and a
/// constant gyro. Constant accel gives sample variance exactly 0 < 0.05, and
/// `|accel| - 9.80665| = 0.09335 < 1.0`, so the chunk is stationary.
///   gyro bias = mean(gyro) = (0.01, -0.02, 0.03)   (norm 0.0374 < 0.05)
///   accel bias z = 9.9 - 9.80665 = 0.09335 exactly; x and y are 0 by design.
#[test]
fn static_imu_biases_are_the_mean_minus_gravity() {
    let samples: Vec<ImuSample> = (0..20)
        .map(|i| {
            imu(
                i * 50_000,
                Vector3::new(0.0, 0.0, 9.9),
                Vector3::new(0.01, -0.02, 0.03),
            )
        })
        .collect();
    let (ba, bg) = estimate_static_imu_biases(&samples);
    assert_vec(
        bg.expect("gyro bias"),
        Vector3::new(0.01, -0.02, 0.03),
        1e-12,
        "gyro bias",
    );
    assert_vec(
        ba.expect("accel bias"),
        Vector3::new(0.0, 0.0, 0.09335),
        1e-12,
        "accel bias",
    );
}

/// Fewer than 20 samples can never fill a 20-wide chunk, so the function
/// short-circuits before any accumulation.
#[test]
fn static_imu_biases_need_at_least_twenty_samples() {
    let samples: Vec<ImuSample> = (0..19)
        .map(|i| {
            imu(
                i * 50_000,
                Vector3::new(0.0, 0.0, 9.9),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect();
    let (ba, bg) = estimate_static_imu_biases(&samples);
    assert!(
        ba.is_none() && bg.is_none(),
        "19 samples must yield no estimate"
    );
}

/// 39 samples: chunk 0 (20 wide) turns hard (gyro 0.5 > 0.05) so it is not
/// stationary; chunk 1 holds only 19 stationary samples, so the accumulated
/// count stays below the 20 minimum even though the total input is >= 20.
/// This exercises the second early return, which the length guard cannot reach.
#[test]
fn static_imu_biases_reject_when_too_few_samples_are_stationary() {
    let mut samples: Vec<ImuSample> = (0..20)
        .map(|i| {
            imu(
                i * 50_000,
                Vector3::new(0.0, 0.0, 9.9),
                Vector3::new(0.5, 0.0, 0.0),
            )
        })
        .collect();
    samples.extend((20..39).map(|i| {
        imu(
            i * 50_000,
            Vector3::new(0.0, 0.0, 9.9),
            Vector3::new(0.01, 0.0, 0.0),
        )
    }));
    let (ba, bg) = estimate_static_imu_biases(&samples);
    assert!(
        ba.is_none() && bg.is_none(),
        "19 stationary samples is below the minimum"
    );
}

// ---------------------------------------------------------------------------
// extract_lever_arm_observations
// ---------------------------------------------------------------------------

/// IMU grid that brackets tow = 100.5 at indices 1, 2, 3 so
/// `summarize_imu_window` finds `idx = 2` (needs idx >= 2 and idx + 2 < len)
/// and forms the centred difference over indices 1..3.
fn bracket_imu(gyro_mid: f64, gyro_step: f64) -> Vec<ImuSample> {
    vec![
        imu(100_000_000, Vector3::zeros(), Vector3::zeros()),
        imu(100_400_000, Vector3::zeros(), Vector3::new(0.0, 0.0, 0.0)),
        imu(
            100_500_000,
            Vector3::new(0.0, 0.0, 9.80665),
            Vector3::new(gyro_mid, 0.0, 0.0),
        ),
        imu(
            100_600_000,
            Vector3::zeros(),
            Vector3::new(gyro_mid + gyro_step, 0.0, 0.0),
        ),
        imu(101_000_000, Vector3::zeros(), Vector3::zeros()),
    ]
}

/// `build_lever_arm_obs` is fully determined by the three positions, the two
/// endpoint velocities, the attitude, and the bracketing IMU samples:
///
///   dt = 101.0 - 100.0 = 1.0 s                       (0.01 < dt <= 2.0)
///   a_gnss_ecef = (v_next - v_prev)/dt = ((0,2,0) - (0,0,0))/1 = (0, 2, 0)
///   g_ecef = -9.80665 * unit(6378137, 0, 0)          = (-9.80665, 0, 0)
///   R_b2e = I (identity attitude), so a_gnss_body = a_gnss_ecef - g_ecef
///                                                   = (9.80665, 2, 0)
///   dt_us = (100600000 - 100400000) * 1e-6 = 2e5 us -> 0.2 s  (> 1e-4)
///   alpha = (gyro_next - gyro_prev)/dt_us = (2 - 0)/0.2 = 10 rad/s^2
#[test]
fn lever_arm_observation_matches_central_difference() {
    let traj = vec![
        epoch_at(
            100.0,
            Some(Vector3::zeros()),
            Some(UnitQuaternion::identity()),
        ),
        epoch_at(100.5, None, Some(UnitQuaternion::identity())),
        epoch_at(
            101.0,
            Some(Vector3::new(0.0, 2.0, 0.0)),
            Some(UnitQuaternion::identity()),
        ),
    ];
    let obs = extract_lever_arm_observations(&traj, &bracket_imu(1.0, 1.0));
    assert_eq!(
        obs.len(),
        1,
        "only the one interior epoch can yield an observation"
    );
    assert_vec(
        obs[0].omega_body,
        Vector3::new(1.0, 0.0, 0.0),
        1e-12,
        "omega",
    );
    assert_vec(
        obs[0].alpha_body,
        Vector3::new(10.0, 0.0, 0.0),
        1e-9,
        "alpha",
    );
    assert_vec(
        obs[0].accel_imu_body,
        Vector3::new(0.0, 0.0, 9.80665),
        1e-12,
        "imu force",
    );
    assert_vec(
        obs[0].accel_gnss_body,
        Vector3::new(9.80665, 2.0, 0.0),
        1e-12,
        "gnss accel",
    );
}

/// With no attitude the DCM falls back to course-over-ground: velocity
/// (0, 2, 0) is due East, giving R_b_e = [[0,0,-1],[1,0,0],[0,-1,0]].
/// Then `a_gnss_body = R^T (a_gnss_ecef - g_ecef)`:
///   R^T (0, 2, 0)           = (2, 0, 0)
///   R^T (-9.80665, 0, 0)    = (0, 0, 9.80665)
///   a_gnss_body              = (2, 0, -9.80665)
#[test]
fn lever_arm_observation_falls_back_to_velocity_heading() {
    let traj = vec![
        epoch_at(
            100.0,
            Some(Vector3::zeros()),
            Some(UnitQuaternion::identity()),
        ),
        epoch_at(100.5, None, None),
        epoch_at(101.0, Some(Vector3::new(0.0, 2.0, 0.0)), None),
    ];
    let obs = extract_lever_arm_observations(&traj, &bracket_imu(1.0, 1.0));
    assert_eq!(obs.len(), 1);
    assert_vec(
        obs[0].accel_gnss_body,
        Vector3::new(2.0, 0.0, -9.80665),
        1e-12,
        "fallback gnss accel",
    );
}

/// Guard table for `extract_lever_arm_observations` and its callees. Each row
/// exercises a distinct early return; every case must yield no observations.
#[test]
fn lever_arm_extraction_rejects_unusable_trajectories() {
    let good_imu = bracket_imu(1.0, 1.0);
    let stop = Some(UnitQuaternion::identity());
    // fewer than 3 epochs
    assert!(extract_lever_arm_observations(&[], &good_imu).is_empty());
    assert!(extract_lever_arm_observations(
        &[epoch_at(100.0, Some(Vector3::zeros()), stop)],
        &good_imu
    )
    .is_empty());
    // no IMU at all
    let three = vec![
        epoch_at(100.0, Some(Vector3::zeros()), stop),
        epoch_at(100.5, None, stop),
        epoch_at(101.0, Some(Vector3::new(0.0, 2.0, 0.0)), stop),
    ];
    assert!(extract_lever_arm_observations(&three, &[]).is_empty());
    // dt above the 2.0 s ceiling
    let slow = vec![
        epoch_at(100.0, Some(Vector3::zeros()), stop),
        epoch_at(100.5, None, stop),
        epoch_at(103.0, Some(Vector3::new(0.0, 2.0, 0.0)), stop),
    ];
    assert!(extract_lever_arm_observations(&slow, &good_imu).is_empty());
    // endpoint velocity missing
    let no_vel = vec![
        epoch_at(100.0, None, stop),
        epoch_at(100.5, None, stop),
        epoch_at(101.0, Some(Vector3::new(0.0, 2.0, 0.0)), stop),
    ];
    assert!(extract_lever_arm_observations(&no_vel, &good_imu).is_empty());
}

/// The dynamic filter keeps an observation when
/// `|omega| > 0.04 || |alpha| > 0.04`. Constant zero gyro makes both zero, so
/// a perfectly well-formed but perfectly static window is discarded.
#[test]
fn lever_arm_extraction_filters_static_windows() {
    let traj = vec![
        epoch_at(
            100.0,
            Some(Vector3::zeros()),
            Some(UnitQuaternion::identity()),
        ),
        epoch_at(100.5, None, Some(UnitQuaternion::identity())),
        epoch_at(
            101.0,
            Some(Vector3::new(0.0, 2.0, 0.0)),
            Some(UnitQuaternion::identity()),
        ),
    ];
    let static_imu = bracket_imu(0.0, 0.0);
    assert!(extract_lever_arm_observations(&traj, &static_imu).is_empty());
}

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

/// The public defaults are part of the calibration contract; a silent change
/// to `max_iterations` or a tolerance would alter every downstream run.
#[test]
fn convergence_criteria_defaults_are_the_documented_values() {
    let crit = CalibrationConvergenceCriteria::default();
    assert_eq!(crit.max_iterations, 3);
    assert_eq!(crit.lever_arm_tol_m, 0.005);
    assert_eq!(crit.accel_bias_tol_mps2, 0.005);
    assert_eq!(crit.gyro_bias_tol_radps, 0.0005);
    let opts = MultiPassCalibrationOptions::default();
    assert!(opts.estimate_lever_arm && opts.estimate_imu_biases && opts.estimate_antenna_offset);
    assert!(opts.reference_points.is_none());
    assert_eq!(
        CalibrationParameters::default(),
        CalibrationParameters::default()
    );
}
