//! Tests for the PPP-AR solver, kept out of `ppp_ar.rs` for the <500 LOC rule.

use nalgebra::{DMatrix, DVector};

use super::ppp_ar::*;

/// Well-conditioned float ambiguities sitting `eps` cycles from known integers.
fn floats_with_eps(integers: &[f64], eps: f64, var: f64) -> (DVector<f64>, DMatrix<f64>) {
    let a = DVector::from_vec(integers.iter().map(|v| v + eps).collect::<Vec<_>>());
    let q = DMatrix::from_diagonal(&DVector::from_element(integers.len(), var));
    (a, q)
}

fn wl_fixed(sat_idx: &[usize]) -> Vec<FixedWideLane> {
    sat_idx.iter().map(|&s| FixedWideLane { sat_idx: s, n_wl: 1, fractional_residual: 0.01 }).collect()
}

// ---------------------------------------------------------------------------
// fix_wide_lane / fix_sd_wide_lane
// ---------------------------------------------------------------------------

#[test]
fn wide_lane_fix_applies_the_bias_before_rounding() {
    // Corrected = 105.42 - 0.40 = 105.02 -> rounds to 105, residual 0.02.
    // A 0.10 gate accepts it; the bias must be subtracted, not added:
    // rounding the raw 105.42 would give 105 with residual 0.42 (rejected).
    let cands = vec![WideLaneCandidate { sat_idx: 10, mw_cycles: 105.42, mw_std_cycles: 0.05, bias_wl_cycles: 0.40 }];
    let fixed = PppArSolver::fix_wide_lane(&cands, 0.10, 0.20);
    assert_eq!(fixed.len(), 1);
    assert_eq!(fixed[0].n_wl, 105);
    assert!((fixed[0].fractional_residual - 0.02).abs() < 1e-12);
}

#[test]
fn wide_lane_fix_rejects_a_noisy_arc_even_when_the_mean_is_near_integer() {
    // Mean is fine (12.23 - 0.25 = 11.98) but sigma 0.30 > max_sigma 0.20:
    // an unconverged arc must not be fixed regardless of its mean.
    let cands = vec![WideLaneCandidate { sat_idx: 1, mw_cycles: 12.23, mw_std_cycles: 0.30, bias_wl_cycles: 0.25 }];
    assert!(PppArSolver::fix_wide_lane(&cands, 0.15, 0.20).is_empty());
}

#[test]
fn sd_wide_lane_differences_out_the_reference_bias() {
    // ref corrected = 50.35 - 0.35 = 50.00
    // sat 10 corrected = 105.02  -> sd = 55.02 -> 55, residual 0.02
    // sat 14 corrected =  82.02  -> sd = 32.02 -> 32, residual 0.02
    let cands = vec![
        WideLaneCandidate { sat_idx: 10, mw_cycles: 105.42, mw_std_cycles: 0.05, bias_wl_cycles: 0.40 },
        WideLaneCandidate { sat_idx: 14, mw_cycles: 82.15, mw_std_cycles: 0.06, bias_wl_cycles: 0.13 },
    ];
    let got = PppArSolver::fix_sd_wide_lane(1, 50.35, 0.35, &cands, 0.10);
    assert_eq!(got.len(), 2);
    assert_eq!(got[0].n_wl, 55);
    assert_eq!(got[1].n_wl, 32);
}

#[test]
fn sd_wide_lane_never_fixes_the_reference_against_itself() {
    let cands = vec![WideLaneCandidate { sat_idx: 1, mw_cycles: 50.35, mw_std_cycles: 0.05, bias_wl_cycles: 0.35 }];
    assert!(PppArSolver::fix_sd_wide_lane(1, 50.35, 0.35, &cands, 0.10).is_empty());
}

// ---------------------------------------------------------------------------
// fix_narrow_lane
// ---------------------------------------------------------------------------

#[test]
fn narrow_lane_fix_requires_four_ambiguities_and_four_wide_lanes() {
    let (a3, q3) = floats_with_eps(&[1.0, 2.0, 3.0], 0.001, 1e-8);
    assert!(PppArSolver::fix_narrow_lane(&a3, &q3, &wl_fixed(&[1, 2, 3, 4]), 3.0).is_none());
    let (a4, q4) = floats_with_eps(&[1.0, 2.0, 3.0, 4.0], 0.001, 1e-8);
    assert!(PppArSolver::fix_narrow_lane(&a4, &q4, &wl_fixed(&[1, 2, 3]), 3.0).is_none());
    assert!(PppArSolver::fix_narrow_lane(&a4, &q4, &wl_fixed(&[1, 2, 3, 4]), 3.0).is_some());
}

#[test]
fn narrow_lane_fix_accepts_a_crisp_set_and_reports_the_wide_lane_ids() {
    let (a, q) = floats_with_eps(&[1.0, 2.0, -3.0, 5.0], 0.002, 1e-8);
    let res = PppArSolver::fix_narrow_lane(&a, &q, &wl_fixed(&[11, 12, 13, 14]), 3.0)
        .expect("crisp set must fix");
    assert!(res.is_fixed);
    assert_eq!(res.fixed_nl, vec![(11, 1), (12, 2), (13, -3), (14, 5)]);
    assert!(res.lambda_ratio >= 3.0);
    assert_eq!(res.fixed_wl.len(), 4);
}

#[test]
fn narrow_lane_fix_refuses_a_set_sitting_exactly_halfway_between_integers() {
    // Floats at exactly n + 0.5: the nearest integer is a perfect tie, so the
    // best and second-best lattice points are equidistant and the ratio is 1.
    // Any ratio gate above 1 must therefore REFUSE. This is the zero-false-fix
    // property in its purest form.
    let (a, q) = floats_with_eps(&[1.0, 2.0, -3.0, 5.0], 0.5, 1e-6);
    let res = PppArSolver::fix_narrow_lane(&a, &q, &wl_fixed(&[1, 2, 3, 4]), 3.0)
        .expect("result is reported, not dropped");
    assert!(!res.is_fixed, "a 50/50 lattice tie must never be fixed");
    assert!(res.fixed_nl.is_empty());
    assert!(res.lambda_ratio < 3.0, "ratio {} should sit near 1", res.lambda_ratio);
}

// ---------------------------------------------------------------------------
// fix_single_diff_ambiguities
// ---------------------------------------------------------------------------

#[test]
fn single_diff_rejects_shape_mismatches_with_a_descriptive_error() {
    let cov4 = DMatrix::identity(4, 4) * 0.01;
    let wl4: Vec<f64> = vec![1.0; 4];
    let amb4: Vec<f64> = vec![1.01, 2.02, 3.03, 4.04];
    let cov3 = DMatrix::identity(3, 3) * 0.01;
    let cov5 = DMatrix::identity(5, 5) * 0.01;
    let wl3: Vec<f64> = vec![1.0; 3];
    let amb3: Vec<f64> = vec![1.01, 2.02, 3.03];
    assert!(PppArSolver.fix_single_diff_ambiguities(&amb3, &cov3, &wl3).is_err());
    assert!(PppArSolver.fix_single_diff_ambiguities(&amb4, &cov5, &wl4).is_err());
    assert!(PppArSolver.fix_single_diff_ambiguities(&amb4, &cov4, &wl3).is_err());
    assert!(PppArSolver.fix_single_diff_ambiguities(&amb4, &cov4, &wl4).is_ok());
}

#[test]
fn single_diff_refuses_an_ambiguous_lattice_instead_of_guessing() {
    // Exactly halfway between integers -> ratio 1 -> must be an error, never a
    // silent fix on a coin flip.
    let amb = vec![1.5, 2.5, 3.5, 4.5];
    let cov = DMatrix::identity(4, 4) * 1e-6;
    let err = PppArSolver
        .fix_single_diff_ambiguities(&amb, &cov, &[1.0_f64; 4])
        .expect_err("ratio test must reject an ambiguous set");
    assert!(matches!(err, crate::estimators::eskf::types::EngineError::Internal(m) if m.contains("Ratio")));
}

#[test]
fn single_diff_scales_metres_into_cycles_by_the_wavelength() {
    // Wavelength 0.19 m: a 1.0 m float is 1/0.19 = 5.2632 cycles. Feeding
    // 1.0 m must resolve to cycle 5, not to metre-rounded 1.
    let amb = vec![5.0 * 0.19 + 0.001, 6.0 * 0.19, 7.0 * 0.19, 8.0 * 0.19];
    let cov = DMatrix::identity(4, 4) * (0.19 * 0.19 * 1e-8);
    let got = PppArSolver
        .fix_single_diff_ambiguities(&amb, &cov, &[0.19_f64; 4])
        .expect("must fix");
    assert_eq!(got, vec![5, 6, 7, 8]);
}

#[test]
fn a_non_positive_wavelength_falls_back_to_one_metre() {
    // The guard substitutes 1.0 rather than dividing by zero, so a metre-scale
    // float set with wavelength 0 must behave exactly as if lambda = 1.
    let amb = vec![1.01, 2.02, 3.03, 4.04];
    let cov = DMatrix::identity(4, 4) * 1e-8;
    let with_zero = PppArSolver
        .fix_single_diff_ambiguities(&amb, &cov, &[0.0_f64; 4])
        .expect("zero wavelength must not divide by zero");
    let with_one = PppArSolver
        .fix_single_diff_ambiguities(&amb, &cov, &[1.0_f64; 4])
        .expect("reference");
    assert_eq!(with_zero, with_one);
    assert_eq!(with_zero, vec![1, 2, 3, 4]);
}

// ---------------------------------------------------------------------------
// form_single_differences / backsubstitute_sd_fix
// ---------------------------------------------------------------------------

#[test]
fn single_differencing_matrix_is_the_expected_pivot_operator() {
    // Reference is index 0; row k must be -1 at the reference and +1 at k+1.
    let (n, m) = (4, 5);
    let amb = vec![0.0; m];
    let cov = DMatrix::identity(m, m) * 0.01;
    let (_sd, _sd_cov, d) = PppArSolver::form_single_differences(&amb, &cov);
    assert_eq!(d.nrows(), n);
    assert_eq!(d.ncols(), m);
    for k in 0..n {
        assert_eq!(d[(k, 0)], -1.0);
        assert_eq!(d[(k, k + 1)], 1.0);
        let row_nz = (0..m).filter(|c| d[(k, *c)] != 0.0).count();
        assert_eq!(row_nz, 2, "row {k} must have exactly two non-zeros");
    }
}

#[test]
fn backsubstitution_reports_failure_on_a_singular_covariance() {
    let m = 4;
    let amb = vec![0.0; m];
    let cov = DMatrix::identity(m, m) * 0.01;
    let (sd_float, sd_cov, d) = PppArSolver::form_single_differences(&amb, &cov);
    let singular = DMatrix::zeros(m - 1, m - 1);
    let out = PppArSolver::backsubstitute_sd_fix(
        &amb,
        &cov,
        &d,
        &singular,
        &sd_float,
        &DVector::zeros(m - 1),
    );
    assert!(out.is_none(), "Cholesky of a zero matrix must fail, not panic");
    assert_eq!(sd_cov.nrows(), m - 1);
}

// ---------------------------------------------------------------------------
// resolve_sd_cascade
// ---------------------------------------------------------------------------

#[test]
fn sd_cascade_rounds_the_wide_lane_biases_before_resolving() {
    let f1 = 1575.42e6_f64;
    let f2 = 1227.60e6_f64;
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let lambda_nl = c / (f1 + f2);
    let alpha = c * f2 / (f1 * f1 - f2 * f2);
    let wl = [3i32, -2, 5, 1];
    let nl = [12i64, -8, 15, 7];
    let m = wl.len() + 1;
    let mut amb = vec![10.0; m];
    for k in 0..wl.len() {
        amb[k + 1] = amb[0] + lambda_nl * nl[k] as f64 + alpha * wl[k] as f64 + 0.001;
    }
    let cov = DMatrix::identity(m, m) * 1e-4;
    // Biases carry a 0.02 cycle residue that must be rounded away first.
    let biases: Vec<f64> = wl.iter().map(|w| *w as f64 + 0.02).collect();
    let (fixed, ratio) = PppArSolver::resolve_sd_cascade(&amb, &cov, f1, f2, &biases, 2.0)
        .expect("cascade must resolve a clean set");
    assert!(ratio >= 2.0);
    for k in 0..wl.len() {
        let got = fixed[k + 1] - fixed[0];
        let want = lambda_nl * nl[k] as f64 + alpha * wl[k] as f64;
        assert!((got - want).abs() < 1e-3, "pair {k}: {got} vs {want}");
    }
}

#[test]
fn sd_cascade_rejects_an_empty_or_short_ambiguity_vector() {
    let f1 = 1575.42e6_f64;
    let f2 = 1227.60e6_f64;
    assert!(PppArSolver::resolve_sd_cascade(&[], &DMatrix::identity(0, 0), f1, f2, &[], 2.0).is_none());
    let amb = vec![1.0, 2.0, 3.0];
    let cov = DMatrix::identity(3, 3) * 0.01;
    assert!(PppArSolver::resolve_sd_cascade(&amb, &cov, f1, f2, &[1.0, 1.0], 2.0).is_none());
}

// ---------------------------------------------------------------------------
// PppMwTracker
// ---------------------------------------------------------------------------

#[test]
fn mw_tracker_averages_consecutive_epochs_and_resets_on_a_slip() {
    let mut t = PppMwTracker::new();
    t.update((0, 1), 50.0, 0, false);
    assert_eq!(t.get_smoothed((0, 1)), Some((50.0, 1)));
    // Running mean: 50 + (52 - 50)/2 = 51 after the second epoch.
    t.update((0, 1), 52.0, 1, false);
    assert_eq!(t.get_smoothed((0, 1)), Some((51.0, 2)));
    // A declared slip discards the arc: count restarts at 1.
    t.update((0, 1), 999.0, 2, true);
    assert_eq!(t.get_smoothed((0, 1)), Some((999.0, 1)));
    assert!(t.get_smoothed((0, 9)).is_none(), "unknown satellite has no arc");
}

#[test]
fn mw_tracker_resets_when_the_epoch_gap_exceeds_two() {
    let mut t = PppMwTracker::new();
    t.update((0, 2), 10.0, 0, false);
    t.update((0, 2), 11.0, 1, false);
    // Epoch 4 is more than two after epoch 1: the arc is not continuous.
    t.update((0, 2), 77.0, 4, false);
    assert_eq!(t.get_smoothed((0, 2)), Some((77.0, 1)));
    // Epoch 3 would still be inside the window and must NOT reset.
    let mut u = PppMwTracker::new();
    u.update((0, 3), 10.0, 1, false);
    u.update((0, 3), 11.0, 3, false);
    assert_eq!(u.get_smoothed((0, 3)), Some((10.5, 2)));
}

#[test]
fn sd_wide_lane_subset_needs_enough_epochs_on_both_ends() {
    let mut t = PppMwTracker::new();
    // Reference with too short an arc -> nothing can be differenced.
    t.update((0, 1), 50.0, 0, false);
    for _ in 0..10 {
        t.update((0, 2), 70.0, 1, false);
    }
    assert!(t.fix_sd_wide_lane_subset((0, 1), &[(0, 2)], 5, 0.35).is_empty());

    // Reference long enough, candidate too short -> candidate skipped.
    let mut u = PppMwTracker::new();
    for _ in 0..10 {
        u.update((0, 1), 50.0, 1, false);
        u.update((0, 2), 70.0, 1, false);
    }
    u.update((0, 2), 70.0, 2, true); // candidate arc restarts at count 1
    let got = u.fix_sd_wide_lane_subset((0, 1), &[(0, 2)], 5, 0.35);
    assert!(got.is_empty(), "one-sided evidence must not fix: {got:?}");
}

#[test]
fn sd_wide_lane_subset_rejects_a_candidate_that_does_not_land_near_an_integer() {
    // 50 -> 70.4 differs by 20.4 cycles: 0.4 from the nearest integer, past
    // a 0.35 gate, so the difference must be refused.
    let mut t = PppMwTracker::new();
    for _ in 0..10 {
        t.update((0, 1), 50.0, 1, false);
        t.update((0, 2), 70.4, 1, false);
    }
    assert!(t.fix_sd_wide_lane_subset((0, 1), &[(0, 2)], 5, 0.35).is_empty());
    // The same arc at 70.35 (0.35 exactly) is on the gate and must fix to 20.
    let mut u = PppMwTracker::new();
    for _ in 0..10 {
        u.update((0, 1), 50.0, 1, false);
        u.update((0, 2), 70.35, 1, false);
    }
    assert_eq!(u.fix_sd_wide_lane_subset((0, 1), &[(0, 2)], 5, 0.35), vec![((0, 2), 20)]);
}