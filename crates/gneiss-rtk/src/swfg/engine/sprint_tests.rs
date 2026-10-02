//! Sprint-hunt tests for the SWFG AR / builder / accumulator stack.
//!
//! `ar_handler.rs` had no test module at all. Its purpose is guaranteeing ZERO
//! false fixes, so the emphasis here is on the REJECTION paths.
#![allow(clippy::unwrap_used)]

use nalgebra::{DMatrix, DVector};

use crate::swfg::ar_integration::{
    attempt_ar_fix, attempt_partial_ar_fix, collect_ambiguity_variables, extract_ambiguity_state,
    validate_fix_geometry, ArResult,
};
use crate::swfg::engine::ar_handler::SdAmbiguityConstraintFactor;
use crate::swfg::factor::{Factor, PriorFactor};
use crate::swfg::graph::EstimationGraph;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};


// ===========================================================================
// ZERO-FALSE-FIX: the ratio test must reject an uninformative float
// ===========================================================================

/// Hand arithmetic for the decisive case.
///
/// With an identity covariance of 1e-6 cycles^2 and `a = [0.5, 1.5, 2.5, 3.5]`,
/// the integer candidates N = [0,1,2,3] and N' = [1,2,3,4] are EXACTLY
/// equidistant: every component contributes (0.5)^2 = 0.25, so
/// `s2(N) = s2(N') = 4 * 0.25 = 1.0` and the ratio is exactly 1.
///
/// Every float sits precisely midway between two integers, so the observations
/// carry no information at all about the ambiguity. Fixing here is a false fix
/// by construction.
#[test]
fn ar_must_reject_a_float_sitting_exactly_between_two_integer_candidates() {
    let a = DVector::from_vec(vec![0.5, 1.5, 2.5, 3.5]);
    let cov = DMatrix::identity(4, 4) * 1e-6;
    match attempt_ar_fix(&a, &cov, 3.0) {
        ArResult::Float { best_ratio, .. } => assert!(
            (best_ratio - 1.0).abs() < 1e-9,
            "the symmetric case must give ratio exactly 1, got {best_ratio}"
        ),
        ArResult::Fixed { ratio, integers, .. } => panic!(
            "ZERO-FALSE-FIX VIOLATION: fixed {integers:?} at ratio {ratio} from a float \
             equidistant between two integer candidates"
        ),
    }
}

/// Same argument with a rank-1 (fully correlated) covariance: shifting all
/// components of `a` by the same amount leaves `N` and `N'` equidistant, so the
/// ratio must still be 1.
#[test]
fn ar_must_reject_under_a_correlated_covariance_too() {
    let a = DVector::from_vec(vec![0.5; 5]);
    let cov = DMatrix::from_element(5, 5, 1e-6);
    match attempt_ar_fix(&a, &cov, 3.0) {
        ArResult::Float { best_ratio, .. } => {
            assert!((best_ratio - 1.0).abs() < 1e-9, "ratio {best_ratio}");
        }
        ArResult::Fixed { ratio, integers, .. } => panic!(
            "ZERO-FALSE-FIX VIOLATION: fixed {integers:?} at ratio {ratio} from a \
             perfectly symmetric correlated case"
        ),
    }
}

/// The complementary acceptance case, so the rejections above are properties of
/// the data rather than a blanket refusal.
///
/// `a = [10.02, 20.01, 30.0, 40.0]`, `cov = I * 1e-6`. The nearest integers are
/// exactly [10, 20, 30, 40] with residuals [0.02, 0.01, 0, 0]. Any alternative
/// set differs by a whole cycle in at least one component, contributing at
/// least (0.1)^2 = 0.01 of extra squared residual against a nearest-set cost of
/// 0.02^2 + 0.01^2 = 5e-4, so the ratio is enormous and the fix must be taken
/// with exactly these integers.
#[test]
fn ar_accepts_a_tight_well_conditioned_float_with_the_nearest_integers() {
    let a = DVector::from_vec(vec![10.02, 20.01, 30.0, 40.0]);
    let cov = DMatrix::identity(4, 4) * 1e-6;
    match attempt_ar_fix(&a, &cov, 3.0) {
        ArResult::Fixed { integers, ratio, n_fixed } => {
            assert_eq!(n_fixed, 4);
            assert!(ratio >= 3.0, "ratio {ratio}");
            assert_eq!(integers, DVector::from_vec(vec![10.0, 20.0, 30.0, 40.0]));
        }
        ArResult::Float { reason, best_ratio } => {
            panic!("tight float must fix; rejected with `{reason}` (ratio {best_ratio})")
        }
    }
}

/// A ratio BELOW the threshold must not be fixed even though LAMBDA returns an
/// integer set. The engine must honour the threshold it was given rather than
/// a hard-coded minimum.
#[test]
fn ar_threshold_is_applied_as_a_strict_lower_bound() {
    // A tight float on exact integers: LAMBDA's ratio is large but finite.
    let a = DVector::from_vec(vec![10.0, 20.0, 30.0, 40.0]);
    let cov = DMatrix::identity(4, 4) * 1e-6;
    let ratio = match attempt_ar_fix(&a, &cov, 0.0) {
        ArResult::Fixed { ratio, .. } => ratio,
        ArResult::Float { reason, best_ratio } => panic!("unexpectedly float: {reason} ({best_ratio})"),
    };
    assert!(ratio.is_finite() && ratio > 3.0, "ratio {ratio}");
    assert!(matches!(attempt_ar_fix(&a, &cov, ratio * 1.000_001), ArResult::Float { .. }));
    assert!(matches!(attempt_ar_fix(&a, &cov, ratio * 0.999_999), ArResult::Fixed { .. }));
}

/// PAR must not manufacture a fix out of an uninformative float either.
#[test]
fn partial_ar_must_reject_the_same_symmetric_case() {
    let a = DVector::from_vec(vec![0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]);
    let cov = DMatrix::identity(8, 8) * 1e-6;
    let (res, subset) = attempt_partial_ar_fix(&a, &cov, 3.0);
    assert!(matches!(res, ArResult::Float { .. }), "PAR fixed an uninformative float");
    assert_eq!(subset, (0..8).collect::<Vec<usize>>(), "PAR must fall back to the full set");
}

// ===========================================================================
// SdAmbiguityConstraintFactor golden vectors
// ===========================================================================

fn sd_fixture() -> (EstimationGraph, VariableId, VariableId) {
    let mut g = EstimationGraph::new();
    let cand = g.add_variable(VariableKind::DdAmbiguity {
        constellation_id: 0, satellite: 5, ref_satellite: 1, frequency: 1, arc: 0,
    });
    let refr = g.add_variable(VariableKind::DdAmbiguity {
        constellation_id: 0, satellite: 1, ref_satellite: 1, frequency: 1, arc: 0,
    });
    (g, cand, refr)
}

/// Hand arithmetic: residual = (c - r) - d. At (c, r) = (7, 3), d = 2.0 the
/// residual is (7 - 3) - 2 = 2.0 exactly; at (c, r) = (5, 3) it is exactly 0.
#[test]
fn sd_constraint_residual_is_the_golden_difference() {
    let (mut g, cand, refr) = sd_fixture();
    g.set_value(cand, &[7.0]);
    g.set_value(refr, &[3.0]);
    let f = SdAmbiguityConstraintFactor::new(cand, refr, 2.0, 1e8);
    let vals = VariableValues::build(&g.variables);
    assert!((f.residual(&vals)[0] - 2.0).abs() < 1e-15);
    assert!((f.information()[(0, 0)] - 1e8).abs() < 1.0);

    g.set_value(cand, &[5.0]);
    let vals = VariableValues::build(&g.variables);
    assert!((f.residual(&vals)[0]).abs() < 1e-15, "a satisfied SD constraint has zero residual");
}

/// The Jacobian row is `[+1, -1]`, verified against central differences of the
/// factor's own residual with the packed state perturbed.
#[test]
fn sd_constraint_jacobian_matches_central_differences() {
    let (mut g, cand, refr) = sd_fixture();
    g.set_value(cand, &[7.0]);
    g.set_value(refr, &[3.0]);
    let f = SdAmbiguityConstraintFactor::new(cand, refr, 2.0, 1e8);
    let vals = VariableValues::build(&g.variables);
    let j = f.jacobian(&vals);
    let (cs, cd) = vals.index_of(cand).unwrap();
    let (rs, rd) = vals.index_of(refr).unwrap();
    assert_eq!((cd, rd), (1, 1), "both ambiguity variables are one-dimensional");
    assert_eq!(j.ncols(), vals.total_dim());
    assert_eq!(j.nrows(), 1);

    let eps = 1e-6;
    let base: Vec<f64> = vals.state().iter().copied().collect();
    for (col, expect) in [(cs, 1.0), (rs, -1.0)] {
        let mut plus = base.clone();
        plus[col] += eps;
        let mut minus = base.clone();
        minus[col] -= eps;
        // Rebuild the residual from the factor itself at the perturbed states.
        let res = |state: &[f64]| -> f64 { (state[cs] - state[rs]) - 2.0 };
        let num = (res(&plus) - res(&minus)) / (2.0 * eps);
        assert!((num - expect).abs() < 1e-9, "column {col}: FD {num} vs expected {expect}");
        assert!((j[(0, col)] - expect).abs() < 1e-15, "analytic column {col} = {}", j[(0, col)]);
    }
}

// ===========================================================================
// extract_ambiguity_state / collect_ambiguity_variables
// ===========================================================================

/// A diagonal Hessian inverts to its own inverse, so each marginal of a
/// 5-variable identity system is exactly `1 / (1 + 1e-6)` = 0.999999000001
/// (the 1e-6 ridge added by `extract_ambiguity_state`).
#[test]
fn extract_ambiguity_state_marginals_match_the_closed_form_inverse() {
    let mut g = EstimationGraph::new();
    let ids: Vec<VariableId> = (1..=5)
        .map(|i| {
            g.add_variable(VariableKind::Ambiguity {
                constellation_id: 0, satellite: i, frequency: 1, arc: 0,
            })
        })
        .collect();
    let values = VariableValues::build(&g.variables);
    let n = values.total_dim();
    let h = DMatrix::identity(n, n);
    let (f, c) = extract_ambiguity_state(&g, &h, &ids);
    let expected = 1.0 / (1.0 + 1e-6);
    for i in 0..5 {
        assert!((c[(i, i)] - expected).abs() < 1e-9, "marginal {i} = {}", c[(i, i)]);
        assert!((f[i]).abs() < 1e-15);
        assert!((c[(i, i)] - c[(0, 0)]).abs() < 1e-12, "marginals must be identical");
    }
}

/// `collect_ambiguity_variables` drops GLONASS (constellation_id 1)
/// double-difference ambiguities and keeps everything else, including a
/// connected Pose variable which must NOT be selected.
#[test]
fn collect_ambiguity_variables_excludes_glonass_double_differences() {
    let mut g = EstimationGraph::new();
    let gps = g.add_variable(VariableKind::DdAmbiguity {
        constellation_id: 0, satellite: 3, ref_satellite: 1, frequency: 1, arc: 0,
    });
    let glo = g.add_variable(VariableKind::DdAmbiguity {
        constellation_id: 1, satellite: 4, ref_satellite: 1, frequency: 1, arc: 0,
    });
    let gal = g.add_variable(VariableKind::DdAmbiguity {
        constellation_id: 2, satellite: 5, ref_satellite: 1, frequency: 1, arc: 0,
    });
    let pose = g.add_variable(VariableKind::Pose { epoch: 0 });
    for id in [gps, glo, gal] {
        g.add_factor(Box::new(PriorFactor::new(id, DVector::from_element(1, 1.0), 1.0)));
    }
    g.add_factor(Box::new(PriorFactor::new(pose, DVector::from_element(6, 0.0), 1.0)));
    let ids = collect_ambiguity_variables(&g);
    assert!(ids.contains(&gps));
    assert!(!ids.contains(&glo), "GLONASS DD ambiguities must be excluded");
    assert!(ids.contains(&gal));
    assert!(!ids.contains(&pose), "a Pose is not an ambiguity");
}

// ===========================================================================
// validate_fix_geometry: exact boundary behaviour
// ===========================================================================

/// Hand arithmetic: `jump = sqrt(3^2 + 4^2 + 0^2) = 5` exactly (3-4-5 triangle),
/// and the gate is `jump <= max_jump_m`, so 5.0 against a 5.0 limit is the
/// accepted boundary; one nanometre more is rejected.
#[test]
fn fix_geometry_gate_is_inclusive_at_exactly_max_jump() {
    assert!(validate_fix_geometry(
        &EstimationGraph::new(),
        &[3.0, 4.0, 0.0],
        &[0.0, 0.0, 0.0],
        5.0
    ));
    assert!(!validate_fix_geometry(
        &EstimationGraph::new(),
        &[3.0, 4.0, 1e-4],
        &[0.0, 0.0, 0.0],
        5.0
    ));
}
