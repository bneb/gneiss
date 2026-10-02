use super::*;
use nalgebra::Vector3;
use crate::swfg::graph::EstimationGraph;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};

fn dd_cp_factor(pose: VariableId, amb: VariableId) -> crate::swfg::pipeline::DdCarrierPhaseFactor {
crate::swfg::pipeline::DdCarrierPhaseFactor {
var_pose: pose,
var_amb: amb,
dd_cp_obs_m: 10.0,
sat_pos: Vector3::new(100.0, 0.0, 0.0),
ref_pos: Vector3::new(0.0, 100.0, 0.0),
base_pos: Vector3::zeros(),
base_dd_range: 100.0,
lambda: 0.19,
variance_m2: 1e-4,
elevation_rad: 1.0,
ref_elevation_rad: 1.0,
is_new_amb: false,
variables: vec![pose, amb],
}
}

// ------------------------------------------------------ partition_variables

#[test]
fn test_partition_variables() {
let mut graph = EstimationGraph::new();

let pose0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
let vel0 = graph.add_variable(VariableKind::Velocity { epoch: 0 });
let pose1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
let amb = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 1, satellite: 1, ref_satellite: 0, frequency: 1, arc: 0 });

// Add a factor connecting pose1 and amb to keep them in unconsumed factors
let factor = dd_cp_factor(pose1, amb);
graph.add_factor(Box::new(factor));

let (marg_ids, kept_ids) = partition_variables(&graph, 0);

assert_eq!(marg_ids.len(), 2);
assert!(marg_ids.contains(&pose0));
assert!(marg_ids.contains(&vel0));

assert_eq!(kept_ids.len(), 2);
assert!(kept_ids.contains(&pose1));
assert!(kept_ids.contains(&amb));
}

#[test]
fn partition_leaves_later_epochs_and_cross_epoch_variables_alone() {
let mut graph = EstimationGraph::new();
let pose0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
let pose1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
let zwd1 = graph.add_variable(VariableKind::TropoZwd { epoch: 1 });
let amb = graph.add_variable(VariableKind::Ambiguity { constellation_id: 0, satellite: 3, frequency: 1, arc: 0 });
let bias = graph.add_variable(VariableKind::ImuBias);
graph.add_factor(Box::new(dd_cp_factor(pose1, amb)));

let (marg, kept) = partition_variables(&graph, 0);
// Epoch 0's pose is the only variable explicitly matched by epoch. The other
// two are swept up by the documented "dead variable" rule (marginalization.rs:
// "no longer observed by any active factor ... safely marginalize out ... to
// keep the state vector small"): `zwd1` and the session-scoped `ImuBias` are
// not referenced by any factor in this graph, so they carry no information and
// are dropped even though `zwd1` still belongs to a live epoch.
assert_eq!(marg, vec![pose0, zwd1, bias]);
// Only variables still coupled to a surviving multi-variable factor are kept.
assert_eq!(kept, vec![pose1, amb]);
}

#[test]
fn partition_marginalizes_variables_left_dead_by_the_window() {
// A variable touched only by *single*-variable factors loses all observation
// when its epoch slides out, so it is swept up with the epoch. A variable in a
// multi-variable factor is still coupled to the survivors and is kept.
let mut graph = EstimationGraph::new();
let pose0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
let dead = graph.add_variable(VariableKind::Attitude { epoch: 0 });
let pose1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
let amb = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 });
graph.add_factor(Box::new(crate::swfg::factor::AttitudePriorFactor::new(dead, 1.0)));
graph.add_factor(Box::new(dd_cp_factor(pose1, amb)));

let (marg, kept) = partition_variables(&graph, 0);
assert!(marg.contains(&pose0));
assert!(marg.contains(&dead), "an unobserved variable must not survive the window");
assert_eq!(kept, vec![pose1, amb]);
}

#[test]
fn partition_of_an_empty_graph_splits_into_nothing() {
let graph = EstimationGraph::new();
let (marg, kept) = partition_variables(&graph, 0);
assert!(marg.is_empty());
assert!(kept.is_empty());
}

// ------------------------------------------------------------- marginalize

#[test]
fn test_marginalize_completely_disconnected() {
let mut graph = EstimationGraph::new();
let _p0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
let p1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
let amb = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 1, satellite: 1, ref_satellite: 0, frequency: 1, arc: 0 });

let factor = dd_cp_factor(p1, amb);
graph.add_factor(Box::new(factor));

let (marg_ids, kept_ids) = partition_variables(&graph, 0);
let values = VariableValues::build(&graph.variables);
let total_dim = values.total_dim();
let hessian = DMatrix::identity(total_dim, total_dim);
let gradient = DVector::zeros(total_dim);

let prior = marginalize(&graph, &hessian, &gradient, &marg_ids, &kept_ids);
assert!(prior.is_some());
let p = prior.unwrap();
assert_eq!(p.variables.len(), kept_ids.len());
}

/// Golden Schur vector.
///
/// State: `pose0` (6-DOF, epoch 0, id 1) is marginalized; `amb1` and `amb2`
/// (1-DOF each, ids 2 and 3) are kept.  Packed columns are
/// `pose0 = 0..6`, `amb1 = 6`, `amb2 = 7`, so the Hessian is partitioned as
///   `H = [[A, B], [B^T, C]]` with A 6x6, B 6x2, C 2x2.
///
/// `A = I + 1 1^T` (identity plus an all-ones matrix), `b1 = [1..6]`,
/// `b2 = [0,1,0,1,0,1]`, `C = [[100, 2], [2, 200]]`.  The retained block is
/// sized so that `S = C - B^T A^-1 B` stays comfortably positive definite:
/// the quadratic forms come out near `b1'A^-1b1 = 28`, `b1'A^-1b2 = 3`,
/// `b2'A^-1b2 = 12/7`.
///
/// The implementation solves with `A_reg = A + 1e-8 I`.  Writing `m = 1 + 1e-8`,
/// Sherman-Morrison on `(m I + 1 1^T)` gives the closed form
///   `A_reg^-1 v = v / m - (sum(v) / (m + 6)) * 1`
/// so every entry below is hand arithmetic on that formula, not a restatement
/// of the solver's block solve.
struct Golden {
h: DMatrix<f64>,
g: DVector<f64>,
}

fn golden_system() -> Golden {
let mut h = DMatrix::zeros(8, 8);
for r in 0..6 {
for c in 0..6 {
h[(r, c)] = 1.0 + if r == c { 1.0 } else { 0.0 };
}
}
let b1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
let b2 = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
for i in 0..6 {
h[(i, 6)] = b1[i];
h[(6, i)] = b1[i];
h[(i, 7)] = b2[i];
h[(7, i)] = b2[i];
}
h[(6, 6)] = 100.0;
h[(6, 7)] = 2.0;
h[(7, 6)] = 2.0;
h[(7, 7)] = 200.0;
let g = DVector::from_vec(vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 5.0, -1.0]);
Golden { h, g }
}

/// `A_reg^-1 v = v/m - (sum v)/(m(m+6)) * 1` for the golden `A = I + 1 1^T`
/// and `m = 1 + 1e-8` (Sherman-Morrison on `(m I + 1 1^T)`).
fn a_reg_inv(v: &[f64; 6]) -> [f64; 6] {
let m = 1.0 + 1e-8;
let s: f64 = v.iter().sum();
let mut out = [0.0; 6];
for i in 0..6 {
out[i] = v[i] / m - s / (m * (m + 6.0));
}
out
}

fn golden_graph() -> (EstimationGraph, VariableId, VariableId, VariableId) {
let mut graph = EstimationGraph::new();
let pose0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
let amb1 = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 });
let amb2 = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 1, arc: 0 });
graph.set_value(pose0, &[1., 2., 3., 0., 0., 0.]);
graph.set_value(amb1, &[11.0]);
graph.set_value(amb2, &[7.0]);
(graph, pose0, amb1, amb2)
}

#[test]
fn marginalize_reproduces_the_hand_computed_schur_complement() {
let (graph, pose0, amb1, amb2) = golden_graph();
let g = golden_system();
let prior = marginalize(&graph, &g.h, &g.g, &[pose0], &[amb1, amb2])
.expect("a rank-2 marginalization is always solvable");

let m = 1.0 + 1e-8;
let b1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
let b2 = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
let ga = [1.0, 0.0, 2.0, 0.0, 3.0, 0.0];
let ainv_b1 = a_reg_inv(&b1);
let ainv_b2 = a_reg_inv(&b2);
let ainv_ga = a_reg_inv(&ga);

// Quadratic forms b^T A^-1 b, expanded from the per-component products.
let quad = |x: [f64; 6], y: [f64; 6]| -> f64 { (0..6).map(|i| x[i] * y[i]).sum::<f64>() };
let bb11 = quad(b1, ainv_b1);
let bb12 = quad(b1, ainv_b2);
let bb22 = quad(b2, ainv_b2);
let bg1 = quad(b1, ainv_ga);
let bg2 = quad(b2, ainv_ga);

// S = C - B^T A^-1 B.
let want_s = [
[100.0 - bb11, 2.0 - bb12],
[2.0 - bb12, 200.0 - bb22],
];
// g_rem = g_b - B^T A^-1 g_a.
let want_g = [5.0 - bg1, -1.0 - bg2];
// sum(b1_i^2) = 1+4+9+16+25+36 = 91, sum(b2_i^2) = 3, sum(b1_i b2_i) = 0+2+0+4+0+6 = 12,
// sum(b1) = 21, sum(b2) = 3.
assert!((bb11 - (91.0 / m - 441.0 / (m * (m + 6.0)))).abs() < 1e-12, "hand check b1'A^-1b1 = {bb11}");
assert!((bb22 - (3.0 / m - 9.0 / (m * (m + 6.0)))).abs() < 1e-12, "hand check b2'A^-1b2 = {bb22}");
assert!((bb12 - (12.0 / m - 63.0 / (m * (m + 6.0)))).abs() < 1e-12, "hand check b1'A^-1b2 = {bb12}");

assert_eq!((prior.hessian.nrows(), prior.hessian.ncols()), (2, 2));
for (r, row) in want_s.iter().enumerate() {
for (c, want) in row.iter().enumerate() {
assert!(
(prior.hessian[(r, c)] - want).abs() < 1e-12,
"S[{r},{c}] = {} want {want}", prior.hessian[(r, c)]
);
}
assert!((prior.gradient[r] - want_g[r]).abs() < 1e-12, "g_rem[{r}] = {}", prior.gradient[r]);
}
// S is the Schur complement of a symmetric block, so it stays symmetric and
// positive definite (S = C - B^T A^-1 B is the inverse-precision of the
// retained block).
assert!((prior.hessian[(0, 1)] - prior.hessian[(1, 0)]).abs() < 1e-15);
assert!(prior.hessian.clone().cholesky().is_some(), "S must be positive definite");
}

#[test]
fn marginalize_packs_x0_from_the_kept_variables_in_order() {
let (graph, pose0, amb1, amb2) = golden_graph();
let g = golden_system();
let prior = marginalize(&graph, &g.h, &g.g, &[pose0], &[amb1, amb2]).expect("prior");
assert_eq!(prior.variables, vec![(amb1, 1), (amb2, 1)]);
// set_value wrote 11.0 and 7.0 into the two ambiguity nodes.
assert_eq!(prior.x0, DVector::from_vec(vec![11.0, 7.0]));
// Reversing the kept order must reverse the packing, proving x0 follows
// `kept_ids` and not the graph's BTreeMap order.
let flipped = marginalize(&graph, &g.h, &g.g, &[pose0], &[amb2, amb1]).expect("prior");
assert_eq!(flipped.variables, vec![(amb2, 1), (amb1, 1)]);
assert_eq!(flipped.x0, DVector::from_vec(vec![7.0, 11.0]));
}

#[test]
fn schur_complement_reproduces_the_full_matrix_inverse() {
// Canonical identity: for `H = [[A,B],[B^T,C]]` the retained block of the full
// inverse is `(H^-1)_bb = S^-1`. Computed here from a full 8x8 inverse, i.e. a
// completely different numerical route from the implementation's 6x6 block
// solve, so an index transposition or a sign slip cannot survive.
let (graph, pose0, amb1, amb2) = golden_graph();
let g = golden_system();
let prior = marginalize(&graph, &g.h, &g.g, &[pose0], &[amb1, amb2]).expect("prior");

let full_inv = g.h.clone().try_inverse().expect("golden H is invertible");
let mut bb = DMatrix::zeros(2, 2);
for r in 0..2 {
for c in 0..2 {
bb[(r, c)] = full_inv[(6 + r, 6 + c)];
}
}
let product = prior.hessian.clone() * bb;
for r in 0..2 {
for c in 0..2 {
let want = if r == c { 1.0 } else { 0.0 };
// 1e-8 of relative slack: the Schur complement is built from `A + 1e-8 I`,
// so `S` differs from the exact `C - B^T A^-1 B` by that regularisation.
assert!(
(product[(r, c)] - want).abs() < 1e-6,
"S * (H^-1)_bb [{r},{c}] = {} want {want}", product[(r, c)]
);
}
}
}

#[test]
fn reduced_newton_step_reproduces_the_full_system_step() {
// The prior's (S, g_rem) encodes a reduced normal system. Solving
// `S dx_b = -g_rem` must give exactly the retained block of `-H^-1 g`,
// which is the whole point of the Schur complement: the marginalized
// variables' influence is carried by g_rem, not discarded.
let (graph, pose0, amb1, amb2) = golden_graph();
let g = golden_system();
let prior = marginalize(&graph, &g.h, &g.g, &[pose0], &[amb1, amb2]).expect("prior");

let reduced = prior.hessian.clone().lu().solve(&(-&prior.gradient));
let reduced = reduced.expect("S is nonsingular");
let full_step = -(g.h.clone().try_inverse().expect("invertible") * g.g.clone());

assert!((reduced[0] - full_step[6]).abs() < 1e-7, "dx_amb1: reduced {} vs full {}", reduced[0], full_step[6]);
assert!((reduced[1] - full_step[7]).abs() < 1e-7, "dx_amb2: reduced {} vs full {}", reduced[1], full_step[7]);
// Sanity: the full step is genuinely non-trivial, so the check has teeth.
assert!(full_step.norm() > 1e-3);
}

#[test]
fn marginalization_of_an_indefinite_block_falls_back_to_qr() {
// After several Schur complements the accumulated Hessian can be indefinite.
// The Cholesky branch then fails and the QR branch must still return the same
// S. With `A = -I(6)` the regularized block is `-(1 - 1e-8) I`, whose inverse
// is exactly `-1/(1-1e-8) I`, so `S = C + (1/(1-1e-8)) * sum b_i^2` by hand.
let mut graph = EstimationGraph::new();
let pose0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
let amb = graph.add_variable(VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 2, frequency: 0, arc: 0 });

let mut h = DMatrix::zeros(7, 7);
for i in 0..6 {
h[(i, i)] = -1.0;
}
let b = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
for i in 0..6 {
h[(i, 6)] = b[i];
h[(6, i)] = b[i];
}
h[(6, 6)] = 10.0;
let grad = DVector::from_element(7, 1.0);

let prior = marginalize(&graph, &h, &grad, &[pose0], &[amb]).expect("QR fallback must still solve");
let sum_sq: f64 = b.iter().map(|x| x * x).sum();
let want = 10.0 + sum_sq / (1.0 - 1e-8);
assert!((prior.hessian[(0, 0)] - want).abs() < 1e-9, "S = {} want {want}", prior.hessian[(0, 0)]);
// g_rem = g_b - B^T A_reg^-1 g_a = 1 - (-1/(1-1e-8)) * sum b_i * 1 = 1 + 21/(1-1e-8)
let want_g = 1.0 + 21.0 / (1.0 - 1e-8);
assert!((prior.gradient[0] - want_g).abs() < 1e-9, "g_rem = {} want {want_g}", prior.gradient[0]);
// The precondition: this Hessian really is indefinite.
assert!(h.clone().cholesky().is_none(), "the fixture must defeat the Cholesky branch");
}

#[test]
fn marginalize_is_a_no_op_for_a_diagonal_system() {
// A marginalized variable uncoupled from the kept block must leave the
// retained information exactly untouched (S = C, g_rem = g_b).
let (graph, pose0, amb1, amb2) = golden_graph();
let mut h = DMatrix::zeros(8, 8);
for i in 0..8 {
h[(i, i)] = 2.0 + i as f64;
}
let grad = DVector::from_vec(vec![1., 2., 3., 4., 5., 6., 7., 8.]);
let prior = marginalize(&graph, &h, &grad, &[pose0], &[amb1, amb2]).expect("prior");
assert_eq!(prior.hessian, DMatrix::from_diagonal(&DVector::from_vec(vec![8.0, 9.0])));
assert_eq!(prior.gradient, DVector::from_vec(vec![7.0, 8.0]));
}

#[test]
fn marginalize_returns_none_for_empty_or_unresolvable_index_sets() {
let graph = EstimationGraph::new();
let hessian = DMatrix::identity(1, 1);
let gradient = DVector::zeros(1);
// Empty marginalized set.
assert!(marginalize(&graph, &hessian, &gradient, &[], &[VariableId::new(0)]).is_none());
// Empty kept set.
assert!(marginalize(&graph, &hessian, &gradient, &[VariableId::new(0)], &[]).is_none());
// Both empty.
assert!(marginalize(&graph, &hessian, &gradient, &[], &[]).is_none());

let (graph, pose0, amb1, amb2) = golden_graph();
let g = golden_system();
// Ids that exist as `VariableId`s but were never added to the graph.
let ghost_a = VariableId::new(900);
let ghost_b = VariableId::new(901);
assert!(marginalize(&graph, &g.h, &g.g, &[ghost_a], &[amb1]).is_none());
assert!(marginalize(&graph, &g.h, &g.g, &[pose0], &[ghost_b]).is_none());
assert!(marginalize(&graph, &g.h, &g.g, &[pose0], &[amb1, amb2]).is_some());
}

// ------------------------------------------------- apply_marginalization

#[test]
fn apply_marginalization_removes_variables() {
let mut graph = EstimationGraph::new();
let _p0 = graph.add_variable(VariableKind::Pose { epoch: 0 });
let p1 = graph.add_variable(VariableKind::Pose { epoch: 1 });
let (marg_ids, _kept_ids) = partition_variables(&graph, 0);

let prior = MarginalPriorFactor {
variables: vec![(p1, 6)],
x0: DVector::zeros(6),
hessian: DMatrix::identity(6, 6),
gradient: DVector::zeros(6),
};
let n_vars_before = graph.n_variables();
apply_marginalization(&mut graph, prior, &marg_ids);
assert!(graph.n_variables() < n_vars_before);
assert!(graph.marginal_prior.is_some());
}

#[test]
fn apply_marginalization_keeps_the_survivors_and_their_prior() {
// The exact count matters: marginalizing epoch 0 must leave precisely the
// epoch-1 pose, and the prior must still address it by id and dimension.
let (mut graph, pose0, amb1, amb2) = golden_graph();
let g = golden_system();
let prior = marginalize(&graph, &g.h, &g.g, &[pose0], &[amb1, amb2]).expect("prior");
apply_marginalization(&mut graph, prior, &[pose0]);

assert!(!graph.variables.contains_key(&pose0));
assert!(graph.variables.contains_key(&amb1));
assert!(graph.variables.contains_key(&amb2));
assert_eq!(graph.n_variables(), 2);
assert_eq!(graph.total_dim(), 2);

let stored = graph.marginal_prior.clone().expect("prior stored");
assert_eq!(stored.variables, vec![(amb1, 1), (amb2, 1)]);
assert_eq!(stored.hessian.nrows(), 2);
assert_eq!(stored.x0.len(), 2);
}

#[test]
fn apply_marginalization_overwrites_any_previous_prior() {
// Window management marginalizes once per epoch; the newest condensed Hessian
// is the only correct one because it has absorbed the previous prior.
let (mut graph, pose0, amb1, amb2) = golden_graph();
let g = golden_system();
let first = marginalize(&graph, &g.h, &g.g, &[pose0], &[amb1, amb2]).expect("prior");
let mut stale = first.clone();
stale.hessian = DMatrix::zeros(2, 2);
apply_marginalization(&mut graph, stale, &[pose0]);
assert!(graph.marginal_prior.as_ref().unwrap().hessian.iter().all(|v| *v == 0.0));

let second = marginalize(&graph, &g.h, &g.g, &[amb1], &[amb2]).expect("prior");
apply_marginalization(&mut graph, second, &[amb1]);
assert_eq!(graph.marginal_prior.as_ref().unwrap().variables, vec![(amb2, 1)]);
assert!(graph.marginal_prior.as_ref().unwrap().hessian.iter().all(|v| *v != 0.0));
}
