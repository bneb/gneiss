//! Golden-vector tests for the SWFG epoch setup (priors, IMU coupling, motion
//! constraints).
//!
//! The attitude assertions use the body-frame convention fixed by the code:
//! the body frame is NED, so the initial attitude is the NED -> ECEF rotation
//! and `q * (1,0,0)` must be the local *north* direction, `q * (0,1,0)` the local
//! *east* and `q * (0,0,1)` the local *down* (`-up`).  Those three products are
//! exactly what distinguishes `R` from `R^T`, which is the easy mistake to make
//! here.
#![allow(clippy::unwrap_used)]

use super::*;
use crate::swfg::config::{EngineConfig, RtkConfig};
use crate::swfg::variables::VariableValues;
use nalgebra::{DMatrix, UnitQuaternion};

/// WGS84 semi-major axis; the receiver is placed on the equator / prime
/// meridian unless a test says otherwise.
const RX: f64 = 6_378_137.0;
/// Standard gravity used by the preintegration factor.
const G: f64 = 9.80665;

fn rx_pos() -> Vector3<f64> {
    Vector3::new(RX, 0.0, 0.0)
}

fn solver() -> SlidingWindowSolver {
    SlidingWindowSolver::new(&EngineConfig::Rtk(RtkConfig::default()))
}

struct G {
    solver: SlidingWindowSolver,
    pose: VariableId,
    vel: Option<VariableId>,
}

/// A pose at `epoch`; with `with_imu`, also the velocity states at `epoch - 1`
/// and `epoch` and the session IMU bias that the preintegration factor needs.
fn graph(epoch: u32, with_imu: bool) -> G {
    let mut solver = solver();
    let pose = solver.graph.add_variable(VariableKind::Pose { epoch });
    let vel = if with_imu {
        solver.graph.add_variable(VariableKind::Velocity { epoch: epoch.saturating_sub(1) });
        Some(solver.graph.add_variable(VariableKind::Velocity { epoch }))
    } else {
        None
    };
    if with_imu {
        solver.graph.add_variable(VariableKind::ImuBias);
    }
    G { solver, pose, vel }
}

impl G {
    fn add_prev_pose(&mut self, epoch: u32) -> VariableId {
        self.solver.graph.add_variable(VariableKind::Pose { epoch })
    }

    fn values(&self) -> VariableValues {
        VariableValues::build(&self.solver.graph.variables)
    }

    /// Indices of factors with a 6-row residual (the relative-pose factor).
    fn rel_pose(&self) -> Vec<usize> {
        self.index_by_rows(6)
    }

    /// Indices of 3-row factors (the NHC / ZUPT velocity factors).
    fn motion(&self) -> Vec<usize> {
        self.index_by_rows(3)
    }

    /// Indices of 15-row factors (the IMU preintegration factor).
    fn preint(&self) -> Vec<usize> {
        self.index_by_rows(15)
    }

    fn index_by_rows(&self, rows: usize) -> Vec<usize> {
        let v = self.values();
        self.solver
            .graph
            .factors
            .iter()
            .enumerate()
            .filter(|(_, f)| f.residual(&v).len() == rows)
            .map(|(i, _)| i)
            .collect()
    }
}

/// A preintegration describing one second of free fall in the local vertical,
/// i.e. the specific force an accelerometer at rest on the equator reads:
///   dp = +0.5 g dt^2 * up,   dv = +g dt * up,   dt = 1 s
fn free_fall_preint() -> ImuPreintegration {
    let mut p = ImuPreintegration::new();
    p.dt = 1.0;
    p.dp = Vector3::new(0.5 * G * 1.0, 0.0, 0.0);
    p.dv = Vector3::new(G, 0.0, 0.0);
    p
}

// ---- attitude initialisation ------------------------------------------------

#[test]
fn attitude_is_initialised_to_the_ned_body_frame_in_ecef() {
    // On the equator / prime meridian the NED -> ECEF rotation is the
    // identity, so north = +Z, east = +Y, down = -X.
    let mut g = graph(1, true);
    let mut att: Option<UnitQuaternion<f64>> = None;
    setup_imu_and_rel_factors(&mut g.solver, &mut att, 1, g.pose, None, rx_pos(), &Some(free_fall_preint()), false, false);

    let q = att.expect("attitude must be initialised when IMU data arrives");
    assert!((q * Vector3::new(1.0, 0.0, 0.0) - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-12, "body +X is not north");
    assert!((q * Vector3::new(0.0, 1.0, 0.0) - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12, "body +Y is not east");
    assert!((q * Vector3::new(0.0, 0.0, 1.0) - Vector3::new(-1.0, 0.0, 0.0)).norm() < 1e-12, "body +Z is not down");
}

#[test]
fn attitude_rotation_is_a_proper_rotation_of_the_ned_axes() {
    // At 45 deg N, 30 deg E the local frame is
    //   north = (-sin45 cos30, -sin45 sin30,  cos45)
    //   east  = (-sin30,          cos30,        0    )
    //   up    = ( cos45 cos30,  cos45 sin30,  sin45)
    // A proper rotation has orthonormal columns and det = +1.
    let lat = 45.0_f64.to_radians();
    let lon = 30.0_f64.to_radians();
    let pos = gneiss_core::coords::llh_to_ecef(Vector3::new(lat, lon, 0.0));
    let mut g = graph(1, true);
    let mut att: Option<UnitQuaternion<f64>> = None;
    setup_imu_and_rel_factors(&mut g.solver, &mut att, 1, g.pose, None, pos, &Some(free_fall_preint()), false, false);

    let q = att.expect("attitude");
    let (s, c) = lat.sin_cos();
    let (sl, cl) = lon.sin_cos();
    let north = Vector3::new(-s * cl, -s * sl, c);
    let east = Vector3::new(-sl, cl, 0.0);
    let up = Vector3::new(c * cl, c * sl, s);
    assert!((q * Vector3::new(1.0, 0.0, 0.0) - north).norm() < 1e-12, "body +X must be north");
    assert!((q * Vector3::new(0.0, 1.0, 0.0) - east).norm() < 1e-12, "body +Y must be east");
    assert!((q * Vector3::new(0.0, 0.0, 1.0) + up).norm() < 1e-12, "body +Z must be down = -up");

    let r = q.to_rotation_matrix().into_inner();
    for i in 0..3 {
        for j in 0..3 {
            let dot = (0..3).map(|k| r[(i, k)] * r[(j, k)]).sum::<f64>();
            let want = if i == j { 1.0 } else { 0.0 };
            assert!((dot - want).abs() < 1e-12, "R^T R [{i},{j}] = {dot}");
        }
    }
    assert!((r.determinant() - 1.0).abs() < 1e-12, "det(R) = {}", r.determinant());
}

#[test]
fn attitude_is_advanced_by_the_preintegrated_rotation() {
    // A 90 deg yaw must carry the body +X axis from north (+Z) to east (+Y).
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    let mut preint = free_fall_preint();
    preint.dq = UnitQuaternion::from_scaled_axis(Vector3::new(0.0, 0.0, std::f64::consts::FRAC_PI_2));
    let mut att: Option<UnitQuaternion<f64>> = None;

    setup_imu_and_rel_factors(&mut g.solver, &mut att, 1, g.pose, Some(prev), rx_pos(), &Some(preint), false, false);

    let q = att.expect("attitude");
    assert!(
        (q * Vector3::new(1.0, 0.0, 0.0) - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12,
        "after a 90 deg yaw the body +X axis must point east"
    );
}

// ---- relative pose ----------------------------------------------------------

#[test]
fn relative_pose_prior_weights_position_by_mode() {
    // diag(W) = [q_pos, q_pos, q_pos, 1, 1, 1] with
    //   q_pos = 1/1e-4 = 10 000  for static PPP, 1/25 = 0.04 otherwise.
    let cases = [(false, false, 1.0 / 25.0), (true, false, 1.0 / 1e-4), (true, true, 1.0 / 25.0)];
    for (is_ppp, is_kinematic, want_pos) in cases {
        let mut g = graph(1, false);
        let prev = g.add_prev_pose(0);
        setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(prev), rx_pos(), &None, is_ppp, is_kinematic);
        let idx = g.rel_pose();
        assert_eq!(idx.len(), 1, "is_ppp={is_ppp} is_kinematic={is_kinematic}");
        let w = g.solver.graph.factors[idx[0]].information();
        for i in 0..3 {
            assert!((w[(i, i)] - want_pos).abs() < 1e-12, "W[{i},{i}] = {} want {want_pos}", w[(i, i)]);
        }
        for i in 3..6 {
            assert_eq!(w[(i, i)], 1.0, "attitude block must stay at 1");
        }
        for i in 0..6 {
            for j in 0..6 {
                if i != j {
                    assert_eq!(w[(i, j)], 0.0, "W must be diagonal");
                }
            }
        }
    }
}

#[test]
fn no_relative_pose_factor_for_a_self_edge_or_a_missing_pose() {
    let mut g = graph(1, false);
    // Self edge: previous pose is the current pose.
    setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(g.pose), rx_pos(), &None, false, false);
    assert_eq!(g.solver.graph.n_factors(), 0, "a pose must not be tied to itself");

    // Dangling edge: the previous pose was marginalized away.
    let mut g = graph(1, false);
    let ghost = VariableId::new(9_999);
    setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(ghost), rx_pos(), &None, false, false);
    assert_eq!(g.solver.graph.n_factors(), 0, "a factor on a missing pose must not be built");
}

#[test]
fn imu_path_takes_precedence_over_the_relative_pose_factor() {
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    setup_imu_and_rel_factors(
        &mut g.solver, &mut None, 1, g.pose, Some(prev), rx_pos(), &Some(free_fall_preint()), false, false,
    );
    assert!(g.rel_pose().is_empty(), "the IMU path must not also add a relative pose factor");
    assert_eq!(g.preint().len(), 1);
}

// ---- IMU preintegration factor and motion constraints ----------------------

/// A preintegration factor built at `epoch` with an explicitly identity
/// attitude.
///
/// The attitude is supplied by the caller instead of being derived from the
/// position on purpose: `init_attitude_from_pos` installs the NED -> ECEF
/// rotation, which is a real 90 deg rotation at the equator, and the residual
/// then compares a body-frame prediction against the world-frame vector
/// `ImuPreintegration::integrate` accumulates (see the frame-mixing note on
/// `imu_residual_rotates_only_the_prediction`).
fn preint_factor(g: &mut G, pos: Vector3<f64>, preint: ImuPreintegration) -> usize {
    let prev = g.add_prev_pose(0);
    let mut att: Option<UnitQuaternion<f64>> = Some(UnitQuaternion::identity());
    setup_imu_and_rel_factors(&mut g.solver, &mut att, 1, g.pose, Some(prev), pos, &Some(preint), false, false);
    let idx = g.preint();
    assert_eq!(idx.len(), 1, "the preintegration factor must be built");
    idx[0]
}

#[test]
fn free_fall_preintegration_gives_a_zero_residual_at_the_nominal_state() {
    // Nominal: p_i = p_j = v_i = v_j = bias = 0, q_i = q_j = identity, dt = 1.
    //   dp_pred = p_j - p_i - v_i dt - 0.5 g dt^2 = -0.5 (-9.80665 up)(1) = +4.903325 up
    //   dv_pred = v_j - v_i - g dt               = -(-9.80665 up)(1) = +9.80665 up
    // Choosing dp and dv equal to those values makes the whole 15-row residual
    // vanish, which pins the gravity magnitude, its direction (toward the
    // geocentre), the 1/2 dt^2 term and the sign of the gravity contribution.
    let mut g = graph(1, true);
    let idx = preint_factor(&mut g, rx_pos(), free_fall_preint());
    let v = g.values();
    let r = g.solver.graph.factors[idx].residual(&v);
    assert_eq!(r.len(), 15);
    for (i, val) in r.iter().enumerate() {
        assert!(val.abs() < 1e-12, "residual[{i}] = {val}; free fall must cancel exactly");
    }
}

#[test]
fn gravity_points_at_the_geocentre_not_away_from_it() {
    // Same construction, but the receiver sits on the polar axis, where the
    // ECEF position is exactly (0, 0, b) so `normalize` is exactly +Z and the
    // gravity term is exactly along Z.  At the equator it is exactly along X:
    // the two cases together pin the direction, and a sign flip would double
    // the residual to 2 x 19.6133 m.
    let pos = Vector3::new(0.0, 0.0, 6_356_752.314_245);
    let mut preint = ImuPreintegration::new();
    preint.dt = 2.0;
    preint.dp = Vector3::new(0.0, 0.0, 0.5 * G * 4.0);
    preint.dv = Vector3::new(0.0, 0.0, G * 2.0);

    let mut g = graph(1, true);
    let idx = preint_factor(&mut g, pos, preint);
    let v = g.values();
    let r = g.solver.graph.factors[idx].residual(&v);
    let dp = r.rows(0, 3).iter().fold(0.0_f64, |a, b| a + b * b).sqrt();
    assert!(dp < 1e-9, "position block norm = {dp}");
    let dv = r.rows(3, 3).iter().fold(0.0_f64, |a, b| a + b * b).sqrt();
    assert!(dv < 1e-9, "velocity block norm = {dv}");
}

#[test]
fn motion_constraints_use_nhc_variances_of_25_01_01() {
    // NHC: variances (25, 0.01, 0.01) m^2/s^2 -> W = diag(0.04, 100, 100).
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    // Stationarity is judged on |dp| / dt, so move the *position* delta:
    // 0.5 m over 1 s is 0.5 m/s, well above the 0.15 m/s threshold.
    let mut preint = free_fall_preint();
    preint.dp = Vector3::new(0.5, 0.0, 0.0);
    setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(prev), rx_pos(), &Some(preint), false, false);

    let idx = g.motion();
    assert_eq!(idx.len(), 1, "only the NHC factor for a moving platform");
    let w = g.solver.graph.factors[idx[0]].information();
    let want = [1.0 / 25.0, 100.0, 100.0];
    for i in 0..3 {
        assert!((w[(i, i)] - want[i]).abs() < 1e-12, "W[{i},{i}] = {} want {}", w[(i, i)], want[i]);
    }
    let vel = g.vel.unwrap();
    assert_eq!(g.solver.graph.factors[idx[0]].variables(), [g.pose, vel]);
}

#[test]
fn a_stationary_interval_adds_a_zero_velocity_update() {
    // 1) is_stationary flag wins regardless of the displacement.
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    let mut preint = free_fall_preint();
    preint.is_stationary = true;
    preint.dp = Vector3::new(0.5, 0.0, 0.0); // 0.5 m/s, above the threshold
    setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(prev), rx_pos(), &Some(preint), false, false);
    let idx = g.motion();
    assert_eq!(idx.len(), 2, "NHC + ZUPT expected for a stationary interval");
    // NHC variances (25, 0.01, 0.01) -> W = diag(0.04, 100, 100); ZUPT
    // variances 1e-4 -> W = diag(10 000, 10 000, 10 000).  Both connect
    // [pose, vel_j], so they are told apart by their information matrix.
    let mut weights: Vec<[f64; 3]> = idx
        .iter()
        .map(|&i| {
            let w = g.solver.graph.factors[i].information();
            [w[(0, 0)], w[(1, 1)], w[(2, 2)]]
        })
        .collect();
    weights.sort_by(|a, b| a[0].total_cmp(&b[0]));
    assert!((weights[0][0] - 0.04).abs() < 1e-12, "NHC weights = {:?}", weights[0]);
    assert!((weights[0][1] - 100.0).abs() < 1e-12, "NHC weights = {:?}", weights[0]);
    for (i, w) in weights[1].iter().enumerate() {
        assert!((w - 10_000.0).abs() < 1e-9, "ZUPT weight {i} = {w}, must be 1/1e-4");
    }

    // 2) speed test: |dp| / dt = 0.1 m/s < 0.15 m/s with dt > 0.05 s.
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    let mut preint = free_fall_preint();
    preint.dp = Vector3::new(0.1, 0.0, 0.0);
    preint.dv = Vector3::new(0.1, 0.0, 0.0);
    setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(prev), rx_pos(), &Some(preint), false, false);
    assert_eq!(g.motion().len(), 2, "0.1 m/s is below the 0.15 m/s stationary threshold");

    // 3) 0.5 m/s is still not stationary, pinning the threshold itself.
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    let mut preint = free_fall_preint();
    preint.dp = Vector3::new(0.5, 0.0, 0.0);
    setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(prev), rx_pos(), &Some(preint), false, false);
    assert_eq!(g.motion().len(), 1, "0.5 m/s is above the 0.15 m/s stationary threshold");

    // 4) same speed but a too-short interval: dt must exceed 0.05 s.
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    let mut preint = free_fall_preint();
    preint.dt = 0.04;
    preint.dp = Vector3::new(0.1 * 0.04, 0.0, 0.0);
    setup_imu_and_rel_factors(&mut g.solver, &mut None, 1, g.pose, Some(prev), rx_pos(), &Some(preint), false, false);
    assert_eq!(g.motion().len(), 1, "a 40 ms interval is too short to call stationary");
}

#[test]
fn imu_path_without_velocity_states_adds_no_factor() {
    // The preintegration factor needs velocity at both epochs and the bias;
    // with none of them present the solver must not build a dangling factor.
    let mut g = graph(1, false);
    let prev = g.add_prev_pose(0);
    let mut att: Option<UnitQuaternion<f64>> = None;
    setup_imu_and_rel_factors(
        &mut g.solver, &mut att, 1, g.pose, Some(prev), rx_pos(), &Some(free_fall_preint()), false, false,
    );
    assert_eq!(g.solver.graph.n_factors(), 0);
    assert!(att.is_some(), "the attitude is still initialised and advanced");
}

// ---- priors -----------------------------------------------------------------

#[test]
fn first_epoch_seeds_the_pose_with_the_initial_position() {
    let mut g = graph(0, false);
    setup_priors(&mut g.solver, &None, 0, g.pose, None, rx_pos(), false, None);
    let pose = g.solver.graph.variables[&g.pose].value.clone();
    assert!((pose[0] - RX).abs() < 1e-12 && pose[1].abs() < 1e-12 && pose[2].abs() < 1e-12, "pose = {pose:?}");
    assert!(pose.rows(3, 3).iter().all(|v| v.abs() < 1e-15), "no IMU means a zero attitude seed");
}

#[test]
fn pose_seed_carries_the_scaled_axis_of_the_current_attitude() {
    // A 90 deg yaw is a rotation vector of pi/2 about +Z.
    let mut g = graph(1, true);
    let q = UnitQuaternion::from_scaled_axis(Vector3::new(0.0, 0.0, std::f64::consts::FRAC_PI_2));
    setup_priors(&mut g.solver, &Some(q), 0, g.pose, None, rx_pos(), true, None);
    let pose = g.solver.graph.variables[&g.pose].value.clone();
    assert!((pose[5] - std::f64::consts::FRAC_PI_2).abs() < 1e-12, "yaw seed = {}", pose[5]);
    assert!(pose[3].abs() < 1e-15 && pose[4].abs() < 1e-15);
}

#[test]
fn initial_pose_prior_uses_the_configured_position_sigma() {
    // W = diag(1/sigma^2, ..., 1, 1, 1) without IMU; the attitude block is left
    // at zero with IMU, which is why the separate velocity/attitude priors
    // exist.  Default sigma is 100 m -> 1e-5.
    for (sigma, want) in [(Some(0.01), 10_000.0), (Some(2.0), 0.25), (None, 1.0 / 100_000.0)] {
        let mut g = graph(0, false);
        setup_priors(&mut g.solver, &None, 0, g.pose, None, rx_pos(), false, sigma);
        let w = g.prior(g.pose).expect("pose prior");
        for i in 0..3 {
            assert!((w[(i, i)] - want).abs() < 1e-12, "sigma={sigma:?} W[{i},{i}] = {} want {want}", w[(i, i)]);
        }
        for i in 3..6 {
            assert_eq!(w[(i, i)], 1.0, "attitude block must be 1 without IMU");
        }
    }

    // With IMU the attitude block is deliberately left unconstrained here.
    let mut g = graph(1, true);
    setup_priors(&mut g.solver, &None, 1, g.pose, None, rx_pos(), true, Some(0.01));
    let w = g.prior(g.pose).expect("pose prior");
    for i in 3..6 {
        assert_eq!(w[(i, i)], 0.0, "attitude is not constrained by the pose prior when IMU is present");
    }
    // ... and a zero velocity prior (variance 25 m^2/s^2 -> W = 0.04) is added.
    let vel = g.vel.unwrap();
    let vw = g.prior(vel).expect("velocity prior");
    assert!((vw[(0, 0)] - 0.04).abs() < 1e-12, "velocity prior W = {}", vw[(0, 0)]);
}

#[test]
fn later_epochs_get_an_attitude_only_prior() {
    // From epoch 1 on, the pose is free but the attitude stays pinned:
    // W = diag(0, 0, 0, 1, 1, 1) with the same mean.
    let mut g = graph(1, false);
    let prev = g.add_prev_pose(0);
    setup_priors(&mut g.solver, &None, 1, g.pose, Some(prev), rx_pos(), false, None);
    let w = g.prior(g.pose).expect("attitude prior");
    for i in 0..3 {
        assert_eq!(w[(i, i)], 0.0, "position must not be re-pinned");
    }
    for i in 3..6 {
        assert_eq!(w[(i, i)], 1.0);
    }
    // setup_priors also seeds the pose from init_pos, so the prior mean and
    // the state coincide and the prior residual is exactly zero.
    let r = g.prior_residual(g.pose).expect("attitude prior");
    assert!(r.iter().all(|v| v.abs() < 1e-15), "prior residual = {r:?}");
    let value = g.solver.graph.variables[&g.pose].value.clone();
    assert!((value[0] - RX).abs() < 1e-12, "pose seeded at {}", value[0]);

    // With IMU the attitude comes from the preintegration factor instead.
    let mut g = graph(1, true);
    let prev = g.add_prev_pose(0);
    setup_priors(&mut g.solver, &None, 1, g.pose, Some(prev), rx_pos(), true, None);
    assert!(g.prior(g.pose).is_none(), "IMU epochs must not add an attitude prior");
}

#[test]
fn no_priors_are_repeated_within_the_same_epoch() {
    // A repeated call with the same pose must neither re-seed the value nor
    // stack a second prior.
    let mut g = graph(1, false);
    setup_priors(&mut g.solver, &None, 1, g.pose, Some(g.pose), rx_pos(), false, None);
    assert_eq!(g.solver.graph.n_factors(), 0, "no priors for a self edge");
    let value = g.solver.graph.variables[&g.pose].value.clone();
    assert!(value.iter().all(|v| v.abs() < 1e-15), "the pose must not be re-seeded mid-epoch");
}

impl G {
    /// Information matrix of the (single) prior attached to `var`, if any.
    fn prior(&self, var: VariableId) -> Option<DMatrix<f64>> {
        self.solver
            .graph
            .factors
            .iter()
            .find(|f| f.variables() == [var])
            .map(|f| f.information())
    }

    fn prior_residual(&self, var: VariableId) -> Option<nalgebra::DVector<f64>> {
        let v = self.values();
        self.solver
            .graph
            .factors
            .iter()
            .find(|f| f.variables() == [var])
            .map(|f| f.residual(&v))
    }
}
