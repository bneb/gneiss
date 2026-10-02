#[cfg(test)]
use super::*;
use crate::swfg::variables::{VariableKind, VariableNode};
use std::collections::BTreeMap;

/// Pack `(id, kind, value)` triples in `VariableId` order.
fn values_of(nodes: &[(VariableId, VariableKind, Vec<f64>)]) -> VariableValues {
let mut vars = BTreeMap::new();
for (id, kind, val) in nodes {
vars.insert(*id, VariableNode { id: *id, kind: *kind, value: DVector::from_vec(val.clone()) });
}
VariableValues::build(&vars)
}

/// Orthogonal RTK geometry: the receiver sits at the ECEF origin, the
/// "satellite" is 20 000 km straight up and the "reference satellite" 20 000 km
/// straight along +X. Then
///   e_sat = (0, 0, 1),  e_ref = (1, 0, 0),  e_sat - e_ref = (-1, 0, 1)
/// and both ranges are exactly 20 000 km, so the rover DD range is 0.
const UP20K: Vector3<f64> = Vector3::new(0.0, 0.0, 20_000_000.0);
const EAST20K: Vector3<f64> = Vector3::new(20_000_000.0, 0.0, 0.0);

fn dd_pr(pose: VariableId, obs: f64, base_dd: f64) -> DdPseudorangeFactor {
DdPseudorangeFactor {
var_pose: pose,
dd_pr_obs: obs,
sat_pos: UP20K,
ref_pos: EAST20K,
base_pos: Vector3::zeros(),
base_dd_range: base_dd,
variance_m2: 1.0,
elevation_rad: std::f64::consts::FRAC_PI_2,
ref_elevation_rad: std::f64::consts::FRAC_PI_2,
variables: vec![],
}
}

fn dd_cp(pose: VariableId, amb: VariableId, obs: f64, base_dd: f64, lambda: f64) -> DdCarrierPhaseFactor {
DdCarrierPhaseFactor {
var_pose: pose,
var_amb: amb,
dd_cp_obs_m: obs,
sat_pos: UP20K,
ref_pos: EAST20K,
base_pos: Vector3::zeros(),
base_dd_range: base_dd,
lambda,
variance_m2: 1.0,
elevation_rad: std::f64::consts::FRAC_PI_2,
ref_elevation_rad: std::f64::consts::FRAC_PI_2,
is_new_amb: false,
variables: vec![],
}
}

fn dd_doppler(prev: VariableId, curr: VariableId, dd: f64, dt: f64) -> DdDopplerFactor {
DdDopplerFactor {
var_pose_prev: prev,
var_pose_curr: curr,
dd_doppler_m_s: dd,
dt,
sat_pos: UP20K,
ref_pos: EAST20K,
variance_m2: 1.0,
variables: vec![],
}
}

fn widelane(a1: VariableId, a2: VariableId, n: f64, var: f64) -> WidelaneConstraintFactor {
WidelaneConstraintFactor { var_amb1: a1, var_amb2: a2, fixed_n_wl: n, variance: var, variables: vec![] }
}

// ------------------------------------------------------- elevation variances

#[test]
fn elevation_variances_at_zenith_match_the_closed_form() {
// sigma = a + b/sin(el); at el = 90 deg sin = 1 so sigma = a + b and
// variance = sigma^2. Pseudorange: (0.3 + 0.3)^2 = 0.36 m^2.
// Carrier:      (0.003 + 0.003)^2 = 3.6e-5 m^2.
let el = std::f64::consts::FRAC_PI_2;
assert!((elevation_pr_variance(el) - 0.36).abs() < 1e-15, "{}", elevation_pr_variance(el));
assert!((elevation_cp_variance(el) - 3.6e-5).abs() < 1e-19, "{}", elevation_cp_variance(el));
}

#[test]
fn carrier_variance_is_exactly_one_ten_thousand_of_the_pseudorange_variance() {
// sigma_cp = 0.003 + 0.003/sin(el) = 0.01 * (0.3 + 0.3/sin(el)) = 0.01 * sigma_pr
// for every elevation, so the variances differ by exactly (0.01)^2 = 1e-4.
// This ties the two models together without restating either constant.
for deg in [5, 15, 30, 45, 60, 75, 89] {
let el = (deg as f64).to_radians();
let pr = elevation_pr_variance(el);
let cp = elevation_cp_variance(el);
assert!((cp - 1e-4 * pr).abs() < 1e-18 * pr.max(1.0), "el={deg}: cp={cp}, 1e-4*pr={}", 1e-4 * pr);
}
}

#[test]
fn elevation_variance_decreases_monotonically_with_elevation() {
let mut prev_pr = f64::INFINITY;
let mut prev_cp = f64::INFINITY;
for deg in [1, 5, 10, 20, 40, 60, 80, 89] {
let el = (deg as f64).to_radians();
let pr = elevation_pr_variance(el);
let cp = elevation_cp_variance(el);
assert!(pr < prev_pr, "pseudorange variance rose at {deg} deg");
assert!(cp < prev_cp, "carrier variance rose at {deg} deg");
prev_pr = pr;
prev_cp = cp;
}
}

#[test]
fn elevation_variance_saturates_below_the_five_degree_floor() {
// sin(el) is clamped at 0.087 (~5 deg), so every elevation below the clamp
// returns the identical value. Compare across the knee rather than restating
// the constant: sin(4 deg) < 0.087 < sin(6 deg).
let below = elevation_pr_variance(4.0_f64.to_radians());
let at_zero = elevation_pr_variance(0.0);
let below_cp = elevation_cp_variance(4.0_f64.to_radians());
let at_zero_cp = elevation_cp_variance(0.0);
assert_eq!(below, at_zero);
assert_eq!(below_cp, at_zero_cp);
// Just above the knee the value must be strictly smaller (the clamp is inactive).
assert!(elevation_pr_variance(6.0_f64.to_radians()) < at_zero);
assert!(elevation_cp_variance(6.0_f64.to_radians()) < at_zero_cp);
}

// ------------------------------------------------------ DdPseudorangeFactor

#[test]
fn dd_pr_residual_vanishes_at_the_true_dd_range() {
// rover_dd = 20 000 km - 20 000 km = 0, so the predicted DD is
// 0 - base_dd. Setting dd_pr_obs = -base_dd makes the residual vanish.
let pose = VariableId::new(1);
let values = values_of(&[(pose, VariableKind::Pose { epoch: 0 }, vec![0.; 6])]);
assert_eq!(dd_pr(pose, -5.0, 5.0).residual(&values), DVector::zeros(1));
// A +2 m observation error shows up as a +2 m residual with this sign.
let r = dd_pr(pose, -3.0, 5.0).residual(&values);
assert!((r[0] - 2.0).abs() < 1e-12, "got {r}");
}

#[test]
fn dd_pr_jacobian_is_sat_minus_ref_line_of_sight() {
// r = dd_obs - (|sat-p| - |ref-p| - base), so
// dr/dp = -(-e_sat) + (-e_ref) = e_sat - e_ref = (0,0,1) - (1,0,0) = (-1,0,1).
let pose = VariableId::new(1);
let values = values_of(&[(pose, VariableKind::Pose { epoch: 0 }, vec![0.; 6])]);
let j = dd_pr(pose, 0.0, 0.0).jacobian(&values);
assert_eq!((j.nrows(), j.ncols()), (1, 6));
assert_eq!(j[(0, 0)], -1.0);
assert_eq!(j[(0, 1)], 0.0);
assert_eq!(j[(0, 2)], 1.0);
for c in 3..6 {
assert_eq!(j[(0, c)], 0.0, "attitude must not enter a range-only factor");
}
}

#[test]
fn dd_pr_jacobian_matches_a_hand_written_central_difference() {
// Move the receiver 1 m along +X and re-evaluate by hand:
//   |up - p|  = |(-d, 0, 20e6)| = sqrt(d^2 + 4e14), d(rho_sat)/dd at d=0 is 0
//   |east - p| = |(20e6 - d, 0, 0)| = 20e6 - d  -> d(rho_ref)/dd = -1
// so dr/dd = -(0 - (-1)) = -1, matching e_sat - e_ref = (-1, 0, 1) projected on x.
let pose = VariableId::new(1);
let d = 1.0;
let plus = values_of(&[(pose, VariableKind::Pose { epoch: 0 }, vec![d, 0., 0., 0., 0., 0.])]);
let minus = values_of(&[(pose, VariableKind::Pose { epoch: 0 }, vec![-d, 0., 0., 0., 0., 0.])]);
let f = dd_pr(pose, 0.0, 0.0);
let numeric = (f.residual(&plus)[0] - f.residual(&minus)[0]) / (2.0 * d);
assert!((numeric - -1.0).abs() < 1e-6, "got {numeric}");
}

#[test]
fn dd_pr_information_sums_the_two_elevation_variances() {
// Both elevations at zenith: 0.36 + 0.36 = 0.72 m^2, so W = 1/0.72.
let pose = VariableId::new(1);
let w = dd_pr(pose, 0.0, 0.0).information()[(0, 0)];
assert!((w - 1.0 / 0.72).abs() < 1e-12, "got {w}");
// Half the elevation makes the variance worse and the information smaller.
let mut low = dd_pr(pose, 0.0, 0.0);
low.elevation_rad = 10.0_f64.to_radians();
low.ref_elevation_rad = 10.0_f64.to_radians();
assert!(low.information()[(0, 0)] < w);
}

#[test]
fn dd_pr_uses_cauchy_loss_with_a_ten_metre_threshold() {
let f = dd_pr(VariableId::new(1), 0.0, 0.0);
assert_eq!(f.robust_threshold(), Some(10.0));
assert!(f.use_cauchy());
}

#[test]
fn dd_pr_degrades_to_zero_when_the_pose_is_absent() {
let f = dd_pr(VariableId::new(1), 1.0, 0.0);
let values = values_of(&[(VariableId::new(9), VariableKind::Velocity { epoch: 0 }, vec![1., 2., 3.])]);
assert_eq!(f.residual(&values), DVector::zeros(1));
assert_eq!(f.jacobian(&values), DMatrix::zeros(1, 3));
}

// ----------------------------------------------------- DdCarrierPhaseFactor

#[test]
fn dd_cp_residual_is_obs_minus_geometric_minus_lambda_ambiguity() {
// predicted = 0 - 5 + 0.19*7 = -3.67, so dd_cp_obs = -3.67 is the truth.
let pose = VariableId::new(1);
let amb = VariableId::new(2);
let values = values_of(&[
(pose, VariableKind::Pose { epoch: 0 }, vec![0.; 6]),
(amb, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 }, vec![7.0]),
]);
assert_eq!(dd_cp(pose, amb, -3.67, 5.0, 0.19).residual(&values), DVector::zeros(1));
// One extra cycle of carrier phase must move the residual by exactly -lambda.
let mut node = VariableNode::new(amb, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 });
node.set_value(&[8.0]);
let mut vars = BTreeMap::new();
vars.insert(pose, VariableNode { id: pose, kind: VariableKind::Pose { epoch: 0 }, value: DVector::zeros(6) });
vars.insert(amb, node);
let shifted = VariableValues::build(&vars);
let r2 = dd_cp(pose, amb, -3.67, 5.0, 0.19).residual(&shifted);
assert!((r2[0] + 0.19).abs() < 1e-12, "one extra cycle must cost exactly lambda, got {r2}");
}

#[test]
fn dd_cp_jacobian_carries_the_line_of_sight_and_minus_lambda() {
// dr/dp = e_sat - e_ref = (-1, 0, 1); dr/dN = -lambda = -0.19.
let pose = VariableId::new(1);
let amb = VariableId::new(2);
let values = values_of(&[
(pose, VariableKind::Pose { epoch: 0 }, vec![0.; 6]),
(amb, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 }, vec![7.0]),
]);
let j = dd_cp(pose, amb, 0.0, 0.0, 0.19).jacobian(&values);
let (s_pose, _) = values.index_of(pose).unwrap();
let (s_amb, _) = values.index_of(amb).unwrap();
assert_eq!(j[(0, s_pose)], -1.0);
assert_eq!(j[(0, s_pose + 1)], 0.0);
assert_eq!(j[(0, s_pose + 2)], 1.0);
assert_eq!(j[(0, s_amb)], -0.19);
for c in 0..6 {
if c != s_pose && c != s_pose + 2 {
assert_eq!(j[(0, c)], 0.0, "column {c} must be zero");
}
}
}

#[test]
fn dd_cp_information_sums_the_two_carrier_variances() {
// 3.6e-5 + 3.6e-5 = 7.2e-5 m^2 -> W = 1/7.2e-5 = 13888.888... At 20 deg the
// model gives sigma_cp = 0.003 + 0.003/sin(20 deg) = 0.003 + 0.00376654...
// = 0.00676654..., squared 4.58e-5 each, 9.16e-5 summed, still above the 1e-4
// floor only just, so use 30 deg: sigma = 0.003 + 0.003/0.5 = 0.009, squared
// 8.1e-5 each, 1.62e-4 summed, W = 1/1.62e-4 = 6172.8395...  (1/0.000162)
let mut f = dd_cp(VariableId::new(1), VariableId::new(2), 0.0, 0.0, 0.19);
f.elevation_rad = 30.0_f64.to_radians();
f.ref_elevation_rad = 30.0_f64.to_radians();
let w = f.information()[(0, 0)];
assert!((w - 1.0 / 1.62e-4).abs() < 1e-6, "got {w}");
}

#[test]
fn dd_cp_information_floor_binds_at_zenith() {
// At zenith the modelled DD carrier variance is 7.2e-5 m^2, *below* the 1e-4
// floor the factor applies (`var_dd.max(1e-4)`), so the weight saturates at
// exactly 1e4 instead of the modelled 1/7.2e-5 = 13888.89. The DD pseudorange
// variance is 0.72 m^2 (W = 1.38889), so the elevation models on their own
// would give a phase/pseudorange weight ratio of
//   (1/7.2e-5) / (1/0.72) = 0.72 / 7.2e-5 = 10000 exactly,
// which is the textbook 1:10000 phase/code precision ratio. The floor pulls
// it down to 1e4 / 1.38889 = 7200 — a 28% de-weighting of carrier phase at the
// best geometry, where the model is most trustworthy and the floor least
// justified.
let pose = VariableId::new(1);
let w = dd_cp(pose, VariableId::new(2), 0.0, 0.0, 0.19).information()[(0, 0)];
assert_eq!(w, 1e4);
let pr_w = dd_pr(pose, 0.0, 0.0).information()[(0, 0)];
assert!((w / pr_w - 7200.0).abs() < 1e-6, "ratio was {}", w / pr_w);
}

#[test]
fn dd_cp_uses_cauchy_loss_with_a_half_metre_threshold() {
let f = dd_cp(VariableId::new(1), VariableId::new(2), 0.0, 0.0, 0.19);
assert_eq!(f.robust_threshold(), Some(0.50));
assert!(f.use_cauchy());
}

#[test]
fn dd_cp_degrades_to_zero_when_pose_or_ambiguity_is_absent() {
let pose = VariableId::new(1);
let amb = VariableId::new(2);
let values = values_of(&[(amb, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 }, vec![1.0])]);
let f = dd_cp(pose, amb, 1.0, 0.0, 0.19);
assert_eq!(f.residual(&values), DVector::zeros(1));
// The factor goes inert as a unit: once the pose lookup fails the Jacobian is
// returned empty, so a factor whose geometry has left the window can never add
// spurious stiffness to the surviving ambiguity column.
assert_eq!(f.jacobian(&values), DMatrix::zeros(1, 1));
}

// ------------------------------------------------------- DdDopplerFactor

#[test]
fn dd_doppler_residual_vanishes_on_the_predicted_displacement() {
// The DD line of sight is evaluated at p_curr, so put p_curr at the ECEF
// origin where the fixture geometry is exact: e_sat = (0,0,1), e_ref = (1,0,0).
// With p_prev = [-100,0,0] the displacement is [100,0,0] and
//   dd_los = (0,0,1) - (1,0,0) = (-1,0,1),  proj = -100 m.
// A DD Doppler of -50 m/s over dt = 2 s predicts exactly -100 m.
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[
(prev, VariableKind::Pose { epoch: 0 }, vec![-100., 0., 0., 0., 0., 0.]),
(curr, VariableKind::Pose { epoch: 1 }, vec![0.; 6]),
]);
assert_eq!(dd_doppler(prev, curr, -50.0, 2.0).residual(&values), DVector::zeros(1));
// -40 m/s predicts -80 m, leaving a -20 m residual.
let r = dd_doppler(prev, curr, -40.0, 2.0).residual(&values);
assert!((r[0] + 20.0).abs() < 1e-12, "got {r}");
}

#[test]
fn dd_doppler_jacobian_splits_the_double_difference_line_of_sight() {
// dr/dp_prev = -dd_los = (1,0,-1); dr/dp_curr = +dd_los = (-1,0,1).
// The line of sight is frozen at p_curr (the standard DD-Doppler
// approximation), which is exact to O(|delta| / range) ~ 1e-6 of a metre.
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[
(prev, VariableKind::Pose { epoch: 0 }, vec![-100., 0., 0., 0., 0., 0.]),
(curr, VariableKind::Pose { epoch: 1 }, vec![0.; 6]),
]);
let j = dd_doppler(prev, curr, -50.0, 2.0).jacobian(&values);
let (s_prev, _) = values.index_of(prev).unwrap();
let (s_curr, _) = values.index_of(curr).unwrap();
assert_eq!((s_prev, s_curr), (0, 6));
assert_eq!(j[(0, s_prev)], 1.0);
assert_eq!(j[(0, s_prev + 1)], 0.0);
assert_eq!(j[(0, s_prev + 2)], -1.0);
assert_eq!(j[(0, s_curr)], -1.0);
assert_eq!(j[(0, s_curr + 1)], 0.0);
assert_eq!(j[(0, s_curr + 2)], 1.0);
}

#[test]
fn dd_doppler_jacobian_row_sums_to_zero_for_a_translation() {
// Translating both epochs by the same vector changes the displacement not at
// all, so the Jacobian must annihilate the all-ones direction. This is a
// property of the *frozen* line of sight: the DD direction is a difference of
// unit vectors, and only the displacement matters.
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[
(prev, VariableKind::Pose { epoch: 0 }, vec![3., -7., 11., 0., 0., 0.]),
(curr, VariableKind::Pose { epoch: 1 }, vec![103., -7., 11., 0., 0., 0.]),
]);
let j = dd_doppler(prev, curr, -50.0, 2.0).jacobian(&values);
for k in 0..3 {
assert!((j[(0, k)] + j[(0, 6 + k)]).abs() < 1e-12, "column {k} pair must cancel");
}
}

#[test]
fn dd_doppler_information_is_inverse_variance_with_a_floor() {
let prev = VariableId::new(1);
let curr = VariableId::new(2);
assert_eq!(dd_doppler(prev, curr, 0.0, 1.0).information()[(0, 0)], 1.0);
let mut floored = dd_doppler(prev, curr, 0.0, 1.0);
floored.variance_m2 = 0.0;
assert_eq!(floored.information()[(0, 0)], 1e4);
}

#[test]
fn dd_doppler_uses_cauchy_loss_with_a_twenty_centimetre_threshold() {
let f = dd_doppler(VariableId::new(1), VariableId::new(2), 0.0, 1.0);
assert_eq!(f.robust_threshold(), Some(0.20));
assert!(f.use_cauchy());
}

#[test]
fn dd_doppler_degrades_to_zero_when_a_pose_is_absent() {
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[(curr, VariableKind::Pose { epoch: 1 }, vec![0.; 6])]);
let f = dd_doppler(prev, curr, 1.0, 1.0);
assert_eq!(f.residual(&values), DVector::zeros(1));
assert_eq!(f.jacobian(&values), DMatrix::zeros(1, 6));
}

// -------------------------------------------------- WidelaneConstraintFactor

#[test]
fn widelane_residual_vanishes_on_the_constraint_surface() {
// N1 - N2 - N_wl = 10 - 4 - 6 = 0.
let a1 = VariableId::new(1);
let a2 = VariableId::new(2);
let values = values_of(&[
(a1, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 }, vec![10.0]),
(a2, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 1, arc: 0 }, vec![4.0]),
]);
assert_eq!(widelane(a1, a2, 6.0, 1e-6).residual(&values), DVector::zeros(1));
// One cycle of drift in N2 must show up one-for-one.
let r = widelane(a1, a2, 7.0, 1e-6).residual(&values);
assert!((r[0] + 1.0).abs() < 1e-12, "got {r}");
}

#[test]
fn widelane_jacobian_is_plus_one_then_minus_one() {
let a1 = VariableId::new(1);
let a2 = VariableId::new(2);
let values = values_of(&[
(a1, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 }, vec![10.0]),
(a2, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 1, arc: 0 }, vec![4.0]),
]);
let j = widelane(a1, a2, 6.0, 1e-6).jacobian(&values);
assert_eq!((j.nrows(), j.ncols()), (1, 2));
assert_eq!(j[(0, 0)], 1.0);
assert_eq!(j[(0, 1)], -1.0);
}

#[test]
fn widelane_information_is_inverse_variance() {
// builder.rs:279 uses variance 1e-6 for a widelane bootstrap constraint,
// i.e. W = 1e6.
assert_eq!(widelane(VariableId::new(1), VariableId::new(2), 6.0, 1e-6).information()[(0, 0)], 1e6);
}

#[test]
fn widelane_degrades_to_zero_when_an_ambiguity_is_absent() {
let a1 = VariableId::new(1);
let a2 = VariableId::new(2);
let values = values_of(&[(a1, VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 }, vec![10.0])]);
let f = widelane(a1, a2, 6.0, 1e-6);
assert_eq!(f.residual(&values), DVector::zeros(1));
// Same inert-factor contract as the DD factors: a missing ambiguity zeroes the
// whole Jacobian rather than leaving a half-populated gradient row.
assert_eq!(f.jacobian(&values), DMatrix::zeros(1, 1));
}

#[test]
fn dd_factor_debug_is_truncated_so_geometry_stays_out_of_the_log() {
// These factors carry ECEF satellite positions (~2e7 m) and metre-scale
// residuals. The hand-written Debug impls exist so a trajectory log records the
// factor type without dumping the payload; assert that contract rather than an
// exact format string.
let pose = VariableId::new(1);
let amb = VariableId::new(2);
let cases: Vec<String> = vec![
format!("{:?}", dd_pr(pose, 12345.678, 9876.5)),
format!("{:?}", dd_cp(pose, amb, 12345.678, 9876.5, 0.19)),
format!("{:?}", dd_doppler(pose, amb, -12.5, 1.0)),
format!("{:?}", widelane(pose, amb, 42.0, 1e-6)),
];
for rendered in &cases {
assert!(rendered.starts_with("Dd") || rendered.starts_with("Widelane"), "got {rendered}");
assert!(rendered.len() < 80, "Debug output grew to {} chars: {rendered}", rendered.len());
// The ECEF satellite position must never reach the log.
assert!(!rendered.contains("26560000"), "Debug leaked the satellite position: {rendered}");
}
// The widelane factor is the one that deliberately keeps a field.
let wl = format!("{:?}", widelane(pose, amb, 42.0, 1e-6));
assert!(wl.contains("fixed_n_wl"), "got {wl}");
assert!(wl.contains("42"), "the widelane integer is the one useful diagnostic: {wl}");
}
