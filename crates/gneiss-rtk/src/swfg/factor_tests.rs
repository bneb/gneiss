
use super::*;
use crate::swfg::variables::{VariableKind, VariableNode};
use std::collections::BTreeMap;

/// Build `VariableValues` from `(id, kind, value)` triples. Variables are
/// packed in `VariableId` order, so the caller controls the offsets. The
/// value length is taken verbatim, which lets a fixture declare a
/// deliberately short "pose" and exercise the `dim.min(...)` guards.
fn values_of(nodes: &[(VariableId, VariableKind, Vec<f64>)]) -> VariableValues {
let mut vars = BTreeMap::new();
for (id, kind, val) in nodes {
let node = VariableNode { id: *id, kind: *kind, value: DVector::from_vec(val.clone()) };
vars.insert(*id, node);
}
VariableValues::build(&vars)
}

/// M-estimator theory ties the loss and the IRLS weight together: for the loss
/// `rho` actually minimised, `rho'(r) = r * w(r)` (the influence function
/// `psi(r) = rho'(r)` equals `r * w(r)`). The two helpers here are written
/// independently, so this cross-check ties them together without restating
/// either closed form.
///
/// `compute_robust_error` carries a harmless global factor of two: the solver
/// multiplies it by `r^T W r / r^2` to form a *cost*, and a constant factor on
/// a cost never moves its minimiser. Both branches carry the same 2:
///   Huber   2kr - k^2       = 2 * (k r - k^2/2)  [standard Huber cost]
///   Cauchy  k^2 ln(1+r^2/k^2) = 2 * (k^2/2) ln(1+r^2/k^2)  [standard Cauchy]
fn assert_loss_and_weight_are_consistent(k: f64, cauchy: bool, scale: f64) {
for &r in &[0.0, 0.5, k, 2.0 * k, 7.0, 100.0] {
let h = 1e-6 * k.max(r.max(1.0));
let d_rho = (compute_robust_error(r + h, k, cauchy)
- compute_robust_error(r - h, k, cauchy))
/ (2.0 * h);
let expected = scale * r * compute_robust_weight(r, k, cauchy);
assert!(
(d_rho - expected).abs() <= 1e-6 * expected.abs().max(1.0),
"rho'(r) != scale*r*w(r) at r={r}, k={k}: d_rho={d_rho}, expected={expected}"
);
}
}

#[test]
fn robust_weight_huber_within_threshold_is_one() {
    assert_eq!(compute_robust_weight(2.0, 3.0, false), 1.0);
}

#[test]
fn robust_weight_huber_above_threshold_downweights() {
    let w = compute_robust_weight(6.0, 3.0, false);
    assert!((w - 0.5).abs() < 1e-12);
}

#[test]
fn robust_weight_cauchy_downweights_squared() {
    // Cauchy with k=3 at r=3 => 1 / (1 + 1) = 0.5
    let w = compute_robust_weight(3.0, 3.0, true);
    assert!((w - 0.5).abs() < 1e-12);
}

#[test]
fn prior_factor_residual_is_zero_at_mu() {
    let kind = VariableKind::Pose { epoch: 0 };
    let id = VariableId::new(0);
    let mut vars = BTreeMap::new();
    let mut node = VariableNode::new(id, kind);
    node.value.copy_from(&DVector::from_vec(vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0]));
    vars.insert(id, node);

    let values = VariableValues::build(&vars);
    let factor = PriorFactor {
        variable: id,
        mu: DVector::from_vec(vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0]),
        information: DMatrix::identity(6, 6),
    };

    let r = factor.residual(&values);
    assert!(r.norm() < 1e-12, "residual at prior mean should be zero");
}

#[test]
fn prior_factor_residual_is_nonzero_away_from_mu() {
    let kind = VariableKind::Velocity { epoch: 0 };
    let id = VariableId::new(0);
    let mut vars = BTreeMap::new();
    let mut node = VariableNode::new(id, kind);
    node.value.copy_from(&DVector::from_vec(vec![0.0, 0.0, 0.0]));
    vars.insert(id, node);

    let values = VariableValues::build(&vars);
    let factor = PriorFactor {
        variable: id,
        mu: DVector::from_vec(vec![1.0, 0.0, 0.0]),
        information: DMatrix::identity(3, 3),
    };

    let r = factor.residual(&values);
    assert!((r.norm() - 1.0).abs() < 1e-12);
}

#[test]
fn relative_pose_factor_residual_and_jacobian() {
    let p1 = VariableId::new(0);
    let p2 = VariableId::new(1);
    let mut vars = BTreeMap::new();

    let mut node1 = VariableNode::new(p1, VariableKind::Pose { epoch: 0 });
    node1.value.copy_from(&DVector::from_vec(vec![10.0, 20.0, 30.0, 0.0, 0.0, 0.0]));
    vars.insert(p1, node1);

    let mut node2 = VariableNode::new(p2, VariableKind::Pose { epoch: 1 });
    node2.value.copy_from(&DVector::from_vec(vec![12.0, 20.0, 30.0, 0.0, 0.0, 0.0]));
    vars.insert(p2, node2);

    let values = VariableValues::build(&vars);
    let rel_factor = RelativePoseFactor::new(p1, p2, 1.0);

    let res = rel_factor.residual(&values);
    assert_eq!(res.len(), 6);
    assert!((res[0] - 2.0).abs() < 1e-12);

    let j = rel_factor.jacobian(&values);
    assert_eq!(j.nrows(), 6);
    assert_eq!(j.ncols(), 12);
    assert_eq!(j[(0, 0)], -1.0);
    assert_eq!(j[(0, 6)], 1.0);
}

// ---------------------------------------------------------------- robust M

#[test]
fn robust_error_huber_is_quadratic_below_threshold() {
    // Huber: ρ(r) = r² for r ≤ k. At r=2, k=3: ρ = 4.
    assert!((compute_robust_error(2.0, 3.0, false) - 4.0).abs() < 1e-12);
}

#[test]
fn robust_error_huber_is_linear_above_threshold() {
    // Huber: ρ(r) = 2kr - k² for r > k. At r=6, k=3: 2·3·6 - 9 = 27.
    assert!((compute_robust_error(6.0, 3.0, false) - 27.0).abs() < 1e-12);
}

#[test]
fn robust_error_huber_is_continuous_at_the_threshold() {
    // Both branches meet at r=k: quadratic branch gives k², linear gives
    // 2k² - k² = k². A sign flip or a k↔k² slip in either branch breaks it.
    let k = 3.0;
    assert!((compute_robust_error(k, k, false) - k * k).abs() < 1e-12);
    // Just below and just above the knee the two branches differ by exactly
    // (k-h)² - (2k(k+h) - k²) = -4kh + h², so the jump is bounded by 4kh + h².
    let h = 1e-10;
    let jump = (compute_robust_error(k - h, k, false) - compute_robust_error(k + h, k, false)).abs();
    // + a few ulps of k^2 for the rounding of each of the three operations.
let ulp = 16.0 * f64::EPSILON * k * k;
assert!(jump <= 4.0 * k * h + h * h + ulp, "knee jump {jump} exceeds the analytic 4kh + h^2");
}

#[test]
fn robust_error_cauchy_matches_hand_logarithm() {
    // Cauchy: ρ(r) = k²·ln(1 + r²/k²). At r=k=3 the argument is 2, so
    // ρ = 9·ln 2 with ln 2 = 0.693147180559945309...
    let expected = 9.0 * std::f64::consts::LN_2;
    assert!((compute_robust_error(3.0, 3.0, true) - expected).abs() < 1e-14);
    // At r=0 the Cauchy loss vanishes, unlike the Huber linear branch's
    // constant offset k² that only applies above the knee.
    assert!(compute_robust_error(0.0, 3.0, true).abs() < 1e-15);
}

#[test]
fn robust_loss_and_weight_obey_the_irls_relation() {
    assert_loss_and_weight_are_consistent(3.0, false, 2.0);
    assert_loss_and_weight_are_consistent(0.5, false, 2.0);
    assert_loss_and_weight_are_consistent(3.0, true, 2.0);
    assert_loss_and_weight_are_consistent(0.2, true, 2.0);
}

#[test]
fn factor_trait_defaults_are_quadratic_l2() {
    // Documented defaults: `robust_threshold() == None` ⇒ L2, and
    // `use_cauchy() == false` ⇒ Huber when a threshold *is* supplied.
    let f = PriorFactor::new(VariableId::new(0), DVector::zeros(3), 1.0);
    assert_eq!(f.robust_threshold(), None);
    assert!(!f.use_cauchy());
}

// ------------------------------------------------------------ PriorFactor

#[test]
fn prior_information_is_inverse_variance() {
    // W = I / variance.  1/0.25 = 4 exactly.
    let f = PriorFactor::new(VariableId::new(0), DVector::zeros(3), 0.25);
    assert_eq!(f.information(), DMatrix::identity(3, 3) * 4.0);
    assert_eq!(f.information(), f.information().transpose());
}

#[test]
fn prior_information_floor_prevents_infinite_weight() {
    // variance = 0 is clamped to 1e-12, so W = 1e12·I — finite but the
    // strongest admissible prior. Without the clamp this would be inf.
    let f = PriorFactor::new(VariableId::new(0), DVector::zeros(2), 0.0);
    assert_eq!(f.information(), DMatrix::identity(2, 2) * 1e12);
    assert!(f.information().iter().all(|v| v.is_finite()));
}

#[test]
fn prior_jacobian_selects_the_right_columns() {
    // Pack a 1-DOF clock (id 1) ahead of a 6-DOF pose (id 2): the pose
    // starts at offset 1, so d r_i / d x_i is the identity at columns 1..7.
    let clock = VariableId::new(1);
    let pose = VariableId::new(2);
    let values = values_of(&[
        (clock, VariableKind::ClockBias { epoch: 0, constellation_id: 0 }, vec![7.0]),
        (pose, VariableKind::Pose { epoch: 0 }, vec![1., 2., 3., 0.1, 0.2, 0.3]),
    ]);
    let f = PriorFactor::new(pose, DVector::zeros(6), 1.0);
    let j = f.jacobian(&values);
    assert_eq!((j.nrows(), j.ncols()), (6, 7));
    for i in 0..6 {
        for c in 0..7 {
            let expect = if c == 1 + i { 1.0 } else { 0.0 };
            assert_eq!(j[(i, c)], expect, "J[{i},{c}]");
        }
    }
    // r = x - mu is exactly matched by J (the prior is linear in x).
    let x = values.get(pose).unwrap().into_owned();
    let mu = DVector::from_vec(vec![0.5, 0.5, 0.5, 0.5, 0.5, 0.5]);
    let f2 = PriorFactor::new(pose, mu.clone(), 1.0);
    assert_eq!(f2.residual(&values), &x - &mu);
}

#[test]
fn prior_jacobian_is_empty_when_the_variable_is_absent() {
    // A factor may outlive its variable during window management. The
    // documented fallback is a 0-row Jacobian, not a panic.
    let other = VariableId::new(9);
    let values = values_of(&[(other, VariableKind::Velocity { epoch: 0 }, vec![1., 2., 3.])]);
    let f = PriorFactor::new(VariableId::new(0), DVector::zeros(6), 1.0);
    assert_eq!(f.jacobian(&values), DMatrix::zeros(0, 3));
}

// ------------------------------------------------------ RelativePoseFactor

#[test]
fn relative_pose_information_is_identity_over_variance() {
    // W = I(6) / 4 → every diagonal entry is exactly 0.25.
    let f = RelativePoseFactor::new(VariableId::new(0), VariableId::new(1), 4.0);
    assert_eq!(f.information(), DMatrix::<f64>::identity(6, 6) * 0.25);
    assert_eq!(f.variables().len(), 2);
}

#[test]
fn relative_pose_jacobian_is_minus_identity_then_identity() {
    // r = x_{p2} - x_{p1} ⇒ J = [-I | +I] over the packed 12-vector.
    let p1 = VariableId::new(1);
    let p2 = VariableId::new(2);
    let values = values_of(&[
        (p1, VariableKind::Pose { epoch: 0 }, vec![0.; 6]),
        (p2, VariableKind::Pose { epoch: 1 }, vec![0.; 6]),
    ]);
    let j = RelativePoseFactor::new(p1, p2, 1.0).jacobian(&values);
    for r in 0..6 {
        for c in 0..12 {
            let expect = match c {
                x if x == r => -1.0,
                x if x == 6 + r => 1.0,
                _ => 0.0,
            };
            assert_eq!(j[(r, c)], expect, "J[{r},{c}]");
        }
    }
}

#[test]
fn relative_pose_jacobian_clamps_rows_beyond_the_variable_dimension() {
    // Pose is nominal 6-DOF, but `dim.min(6)` guards a short variable: a
    // 3-DOF node leaves rows 3..5 of J untouched rather than panicking.
    let p1 = VariableId::new(1);
    let p2 = VariableId::new(2);
    let values = values_of(&[
        (p1, VariableKind::Pose { epoch: 0 }, vec![1., 2., 3.]),
        (p2, VariableKind::Pose { epoch: 1 }, vec![0.; 6]),
    ]);
    let j = RelativePoseFactor::new(p1, p2, 1.0).jacobian(&values);
    assert_eq!((j.nrows(), j.ncols()), (6, 9));
    for r in 0..3 {
        assert_eq!(j[(r, r)], -1.0);
        assert_eq!(j[(r, 3 + r)], 1.0);
    }
    for r in 3..6 {
        assert_eq!(j[(r, r.min(5))], 0.0, "row {r} must stay zero for a 3-DOF p1");
    }
}

// ------------------------------------------------------ AttitudePriorFactor

#[test]
fn attitude_prior_pins_only_the_rotation_block() {
    let pose = VariableId::new(1);
    let values = values_of(&[
        (pose, VariableKind::Pose { epoch: 0 }, vec![10., 20., 30., 0.4, -0.5, 0.6]),
    ]);
    let f = AttitudePriorFactor::new(pose, 0.5);
    // W = I(3)/0.5 = 2·I.
    assert_eq!(f.information(), DMatrix::<f64>::identity(3, 3) * 2.0);
    // r = the last three components, in order.
    let r = f.residual(&values);
    assert_eq!(r, DVector::from_vec(vec![0.4, -0.5, 0.6]));
    // J picks columns 3, 4, 5 — the translation columns stay zero.
    let j = f.jacobian(&values);
    assert_eq!((j.nrows(), j.ncols()), (3, 6));
    assert_eq!(j, DMatrix::from_element(3, 6, 0.0)
        + DMatrix::from_fn(3, 6, |r, c| if r + 3 == c { 1.0 } else { 0.0 }));
    assert_eq!(f.variables(), &[pose]);
}

#[test]
fn attitude_prior_jacobian_offsets_past_leading_variables() {
    // Two 1-DOF scalars packed ahead of the pose push it to offset 2; the
    // rotation columns must move with it (start+3..start+5), not stay at 3..5.
    let zwd = VariableId::new(1);
    let clk = VariableId::new(2);
    let pose = VariableId::new(3);
    let values = values_of(&[
        (zwd, VariableKind::TropoZwd { epoch: 0 }, vec![0.1]),
        (clk, VariableKind::ClockBias { epoch: 0, constellation_id: 0 }, vec![12.0]),
        (pose, VariableKind::Pose { epoch: 0 }, vec![0.; 6]),
    ]);
    let j = AttitudePriorFactor::new(pose, 1.0).jacobian(&values);
    assert_eq!((j.nrows(), j.ncols()), (3, 8));
    assert_eq!(j[(0, 5)], 1.0);
    assert_eq!(j[(1, 6)], 1.0);
    assert_eq!(j[(2, 7)], 1.0);
    assert_eq!(j[(0, 3)], 0.0, "must not use the un-offset columns");
}
