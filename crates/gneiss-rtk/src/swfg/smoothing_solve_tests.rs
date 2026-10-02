use super::tests::{batch_cost, constellation, constellation_n, make_obs_n, make_obs_offset, make_test_obs};
use super::*;

// ------------------------------------------------------------------- solve

/// Build a two-epoch batch with a static receiver at a known truth position.
fn static_batch() -> (BatchFactorGraph, Vector3<f64>) {
let t0 = GpsTime::new(2200, 100.0);
let truth = Vector3::new(-3963427.0, 3350882.0, 3694866.0);
let cfg = BatchSmoothingConfig { enable_ar: false, ..Default::default() };
let mut batch = BatchFactorGraph::new(cfg, constellation(t0));
batch.add_epoch(0, &make_test_obs(t0), None, Some(truth));
batch.add_epoch(1, &make_test_obs(t0), None, Some(truth));
(batch, truth)
}

#[test]
fn solve_returns_one_entry_per_epoch_with_the_right_times() {
let (mut batch, _) = static_batch();
let traj = batch.solve().expect("batch solve");
assert_eq!(traj.len(), batch.pose_ids.len());
assert_eq!(traj.len(), 2);
for (i, (t, _)) in traj.iter().enumerate() {
assert_eq!(*t, batch.epochs[i], "trajectory must follow the insertion order");
}
}

#[test]
fn solve_leaves_the_graph_populated_for_a_second_solve() {
// `solve` moves the graph into the solver and back; losing the graph would
// make the batch unusable, and the variable set must not drift across runs.
let (mut batch, _) = static_batch();
batch.solve().expect("first solve");
let after_first = batch.graph.variables.len();
assert!(after_first >= 2 + 2 * 2, "poses, clocks, zeniths plus ambiguities");
batch.solve().expect("second solve");
assert_eq!(batch.graph.variables.len(), after_first);
}

#[test]
fn solve_reduces_the_batch_cost_and_lands_on_a_local_minimum() {
let (mut batch, _) = static_batch();
let cost_before = batch_cost(&batch.graph);
let traj = batch.solve().expect("batch solve");
assert_eq!(traj.len(), 2);
let cost_after = batch_cost(&batch.graph);
assert!(cost_after < cost_before, "Levenberg-Marquardt must descend: {cost_before} -> {cost_after}");

// Perturbing any pose by a decimetre must not lower the cost.
for (i, pose) in batch.pose_ids.clone().into_iter().enumerate() {
for axis in 0..3 {
let original = batch.graph.variables.get(&pose).unwrap().value.clone();
let mut nudged = original.clone();
nudged[axis] += 1.0;
batch.graph.set_value(pose, nudged.as_slice());
let perturbed = batch_cost(&batch.graph);
assert!(
perturbed >= cost_after - 1e-9,
"moving pose epoch {i} axis {axis} by 1 m lowered the cost from {cost_after} to {perturbed}"
);
batch.graph.set_value(pose, original.as_slice());
}
}
assert!((batch_cost(&batch.graph) - cost_after).abs() < 1e-9, "cost must be restored");
}

/// A constellation of `n` satellites tracked over two epochs shares `n`
/// ambiguity variables across the epochs (`ensure_ambiguity` keys on the arc).
fn two_epoch_batch(n: u8, enable_ar: bool) -> BatchFactorGraph {
let t0 = GpsTime::new(2200, 100.0);
let truth = Vector3::new(-3963427.0, 3350882.0, 3694866.0);
let cfg = BatchSmoothingConfig { enable_ar, ..Default::default() };
let mut batch = BatchFactorGraph::new(cfg, constellation_n(t0, n));
batch.add_epoch(0, &make_obs_n(t0, n), None, Some(truth));
batch.add_epoch(1, &make_obs_n(t0, n), None, Some(truth));
batch
}

#[test]
fn ambiguities_are_shared_across_epochs() {
// The same satellite tracked in two epochs is ONE ambiguity variable, so the
// graph carries per-epoch poses/clocks/zeniths but a single ambiguity state
// per satellite. This is what lets the batch smoothing hold the integer fixed.
let batch = two_epoch_batch(8, false);
let count = |k: fn(&VariableKind) -> bool| batch.graph.variables.values().filter(|n| k(&n.kind)).count();
assert_eq!(count(|k| matches!(k, VariableKind::Pose { .. })), 2);
assert_eq!(count(|k| matches!(k, VariableKind::ClockBias { .. })), 2);
assert_eq!(count(|k| matches!(k, VariableKind::TropoZwd { .. })), 2);

// Sharing is the point: two epochs of the same satellites must produce FEWER
// ambiguity states than carrier-phase factors, because `ensure_ambiguity`
// keys on (constellation, satellite, frequency, arc) and not on the epoch.
// Elevation masking means the exact satellite count is geometry-dependent, so
// assert the sharing relation rather than a hard-coded number.
let n_amb = count(|k| matches!(k, VariableKind::Ambiguity { .. }));
let n_cp_factors = batch.graph.factors.iter()
.filter(|f| f.variables().iter().any(|v| matches!(batch.graph.variables[v].kind, VariableKind::Ambiguity { .. })))
.count();
assert!(n_amb > 0, "the fixture must produce carrier phase at all");
assert!(n_cp_factors > n_amb,
"two epochs share one ambiguity per satellite: {n_cp_factors} carrier-phase factors vs {n_amb} ambiguity states");
}

#[test]
fn solve_with_ar_enabled_skips_the_fix_below_four_ambiguities() {
// Three satellites produce three ambiguity states, which cannot support an
// integer validation, so the AR branch must be skipped and the float solution
// returned unchanged.
let mut with_ar = two_epoch_batch(3, true);
let mut without_ar = two_epoch_batch(3, false);
assert!(with_ar.config.enable_ar);
assert!(!without_ar.config.enable_ar);
let a = with_ar.solve().expect("AR-enabled solve");
let b = without_ar.solve().expect("float solve");
assert_eq!(a.len(), b.len());
for (x, y) in a.iter().zip(b.iter()) {
assert_eq!(x.0, y.0);
assert!((x.1 - y.1).norm() < 1e-9, "AR path must be a no-op below four ambiguities");
}
}

#[test]
fn solve_with_ar_enabled_runs_the_ambiguity_branch() {
// Eight ambiguities clear the `>= 4` gate: the float ambiguity state is
// extracted, validated, and — if it passes — the graph is re-solved with the
// integers pinned. Whatever the validation decides, the batch must still
// return one finite position per epoch.
let mut batch = two_epoch_batch(8, true);
let traj = batch.solve().expect("AR-enabled solve");
assert_eq!(traj.len(), 2);
for (_, p) in traj.iter() {
assert!(p.iter().all(|c| c.is_finite()), "non-finite position {p:?}");
assert!(p.norm() > 6.0e6, "a solution off the Earth is not a solution: {p:?}");
}
assert_eq!(batch.config.ar_ratio_threshold, 2.0);
}

#[test]
fn solve_propagates_an_unsolvable_graph_as_an_error_string() {
// A graph with an orphan variable cannot be linearised; the error must come
// back as a message, not a panic.
let cfg = BatchSmoothingConfig { enable_ar: false, ..Default::default() };
let mut batch = BatchFactorGraph::new(cfg, Vec::new());
let orphan = batch.graph.add_variable(VariableKind::Velocity { epoch: 0 });
let _ = orphan;
let err = batch.solve().expect_err("an orphan variable must be rejected");
assert!(err.contains("Batch solve failed"), "unexpected message: {err}");
}

/// The double-differenced branch builds each pseudorange as
/// `(P_rover - P_base) + |sat - base_pos|`.  Because `add_epoch` derives its
/// clock seed from `P - (|sat - init| - sat_clock + tropo + iono)`, adding a
/// constant `delta` to *every* base pseudorange must shift the seeded receiver
/// clock by exactly `-delta`, whatever the per-satellite geometry is: the
/// satellite-specific terms cancel inside the median.
///
/// This pins the sign and the side of the subtraction — a `P_base - P_rover`
/// flip or a doubled `range_b` term would both move the seed differently.
#[test]
fn double_differencing_subtracts_the_base_pseudorange() {
let t0 = GpsTime::new(2200, 100.0);
let truth = Vector3::new(-3963427.0, 3350882.0, 3694866.0);
let delta = 1000.0;

let seed_of = |base_pr_offset: f64| -> f64 {
let obs = make_obs_offset(t0, 8, 0.0);
let base = make_obs_offset(t0, 8, base_pr_offset);
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), constellation_n(t0, 8));
batch.add_epoch(0, &obs, Some((&base, truth)), Some(truth));
batch.graph.variables.values()
.find(|n| matches!(n.kind, VariableKind::ClockBias { .. }))
.expect("a clock state exists").value[0]
};

let plain = seed_of(0.0);
let shifted = seed_of(delta);
let expected = -delta;
// Tolerance 1e-3 m: `extract_raw_observations` iterates the satellite
// transmit time from the pseudorange, so a base-side range change of 1000 m
// shifts each computed satellite position by a few tenths of a millimetre and
// the two seeds no longer cancel exactly. That is physics, not round-off, and
// it is still ~6 orders of magnitude tighter than the ±1000 m a sign flip or a
// doubled `range_b` term would produce.
let moved = shifted - plain;
assert!(
(moved - expected).abs() < 1e-3,
"adding {delta} m to every base pseudorange moved the clock seed by {moved} instead of {expected}"
);
}

#[test]
fn a_double_differenced_epoch_still_produces_measurements() {
let t0 = GpsTime::new(2200, 100.0);
let truth = Vector3::new(-3963427.0, 3350882.0, 3694866.0);
let obs = make_obs_n(t0, 8);
let mut batch = BatchFactorGraph::new(BatchSmoothingConfig::default(), constellation_n(t0, 8));
batch.add_epoch(0, &obs, Some((&obs, truth)), Some(truth));
assert!(!batch.graph.variables.values().any(|n| matches!(n.kind, VariableKind::TropoZwd { .. })));
// Every satellite seen by the rover must also be seen by the base, so the
// `find`-and-keep filter must not silently drop the whole constellation.
assert!(batch.graph.factors.len() > 2, "expected the 6-DOF anchor plus measurement factors, got {}", batch.graph.factors.len());
assert!(batch.graph.validate_graph_structure().is_ok(), "a double-differenced epoch must leave no orphan state");
}
