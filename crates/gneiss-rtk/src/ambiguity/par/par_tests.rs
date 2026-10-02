//! Tests for the PAR candidate-quality metric and subset selector.
//!
//! Declared from `par.rs` so the private LDL^T helpers are reachable.

use nalgebra::{DMatrix, DVector};

use super::*;

/// Standard normal CDF at 1 sigma: Phi(1) = 0.8413447460685429.
const PHI_1: f64 = 0.8413447460685429;

// ---------------------------------------------------------------------------
// compute_cqm
// ---------------------------------------------------------------------------

#[test]
fn cqm_matches_its_weighted_definition_by_hand() {
    // el = 30 deg -> sin = 0.5
    //   0.25 * 0.5                       = +0.125
    // snr = 50 -> (50-20)/30 = 1.0 (clamped)
    //   0.25 * 1.0                       = +0.250
    // lock = 60 s -> min(1, 60/30) = 1.0
    //   0.20 * 1.0                       = +0.200
    // cmc = 0.5 m -> 0.5/2 = 0.25
    //  -0.15 * 0.25                      = -0.0375
    // q_ii = 0.04 -> sigma = 0.2; float 10.25 -> frac = 0.25
    //   var_term = min(0.2 + 0.125, 2.0) = 0.325
    //  -0.15 * 0.325                     = -0.04875
    // total = 0.125 + 0.250 + 0.200 - 0.0375 - 0.04875 = 0.48875
    let m = AmbiguityMetadata::new(std::f64::consts::PI / 6.0, 50.0, 60.0, 0.5);
    let got = compute_cqm(Some(&m), 0.04, 10.25);
    assert!((got - 0.48875).abs() < 1e-12, "got {got}");
}

#[test]
fn cqm_default_metadata_matches_hand_arithmetic() {
    // defaults: el = 45 deg -> sin = 0.7071067811865476
    //   0.25 * 0.7071067811865476        = +0.1767766952966369
    // snr = 42 -> (42-20)/30 = 0.7333333333
    //   0.25 * 0.7333333333              = +0.1833333333333333
    // lock = 30 -> min(1, 1.0) = 1.0
    //   0.20                            = +0.2
    // cmc = 0.20 -> 0.1 ; -0.15 * 0.1    = -0.015
    // var_term = 0.325 -> -0.15 * 0.325  = -0.04875
    // total = 0.4963600286299702
    let got = compute_cqm(Some(&AmbiguityMetadata::default()), 0.04, 10.25);
    assert!((got - 0.4963600286299702).abs() < 1e-12, "got {got}");
}

#[test]
fn cqm_without_metadata_falls_back_to_neighbourhood_penalty() {
    // fallback = -(frac + 0.5 * sigma)
    // frac = |10.25 - 10| = 0.25 ; sigma = sqrt(0.04) = 0.2 -> 0.5*0.2 = 0.1
    // => -(0.25 + 0.1) = -0.35
    assert!((compute_cqm(None, 0.04, 10.25) + 0.35).abs() < 1e-12);
    // A float that is already an integer has frac = 0 -> -(0 + 0.1) = -0.1.
    assert!((compute_cqm(None, 0.04, 10.0) + 0.1).abs() < 1e-12);
}

#[test]
fn cqm_clamps_out_of_range_quality_inputs() {
    // snr 10 -> (10-20)/30 = -1/3, clamped to 0 -> no snr credit.
    // lock 0 -> 0. cmc 10 -> clamp(10/2, 0, 2) = 2.
    let lo = AmbiguityMetadata::new(0.0, 10.0, 0.0, 10.0);
    // 0.25*0 + 0.25*0 + 0.20*0 - 0.15*2 - 0.15*var_term
    // var_term = min(0.2 + 0.5*0.05, 2.0) = 0.225
    // total = -0.30 - 0.03375 = -0.33375
    let expect = -0.30 - 0.15 * (0.2 + 0.5 * 0.05);
    assert!((compute_cqm(Some(&lo), 0.04, 10.05) - expect).abs() < 1e-12);
    // snr 60 -> (60-20)/30 > 1, clamped to 1; lock 300 -> clamped to 1.
    let hi = AmbiguityMetadata::new(std::f64::consts::FRAC_PI_2, 60.0, 300.0, 0.0);
    assert!((compute_cqm(Some(&hi), 0.04, 10.0) - (0.25 + 0.25 + 0.20 - 0.15 * 0.2)).abs() < 1e-12);
}

#[test]
fn cqm_var_term_saturates_so_one_bad_pair_cannot_dominate() {
    // sigma huge -> var_term capped at 2.0, so the penalty is at most 0.30.
    let m = AmbiguityMetadata::new(0.0, 10.0, 0.0, 0.0);
    let g1 = compute_cqm(Some(&m), 1.0e12, 10.5);
    let g2 = compute_cqm(Some(&m), 1.0e18, 10.5);
    assert!((g1 - g2).abs() < 1e-12, "cap makes arbitrarily bad variance equal");
    assert!((g1 - (-0.30)).abs() < 1e-12, "g1 = -0.15*2.0 = -0.30");
}

// ---------------------------------------------------------------------------
// cumulative_success_rate
// ---------------------------------------------------------------------------

#[test]
fn success_rate_at_one_sigma_is_the_classical_value() {
    // d = 0.25 -> sigma = 0.5 -> z = 1/(2*0.5) = 1
    // phi = 0.5*(1 + erf(1/sqrt2)) = Phi(1) = 0.8413447460685429
    // rate = 2*phi - 1 = 0.6826894921370859
    let d = DVector::from_vec(vec![0.25]);
    assert!((cumulative_success_rate(&d) - 0.6826894921370859).abs() < 1e-12);
}

#[test]
fn success_rate_multiplies_across_ambiguities() {
    // two independent 1-sigma ambiguities: (2*Phi(1)-1)^2 = 0.68268949214^2
    let d = DVector::from_vec(vec![0.25, 0.25]);
    let expect = 0.6826894921370859_f64 * 0.6826894921370859;
    assert!((cumulative_success_rate(&d) - expect).abs() < 1e-12);
    assert!((cumulative_success_rate(&d) - PHI_1 * 2.0 - 1.0).abs() > 0.4, "two are worse than one");
}

#[test]
fn success_rate_ignores_non_positive_pivot_and_is_clamped() {
    // A zero pivot contributes sigma = max(sqrt(0), 1e-12) -> z huge -> rate 1.
    let d = DVector::from_vec(vec![0.0, 0.25]);
    assert!((cumulative_success_rate(&d) - 0.6826894921370859).abs() < 1e-6);
    // Negative pivots are skipped entirely (rate unchanged), and the result
    // never leaves [0, 1].
    let neg = DVector::from_vec(vec![-1.0]);
    assert_eq!(cumulative_success_rate(&neg), 1.0);
    assert!(cumulative_success_rate(&DVector::from_vec(vec![1.0e12])).abs() >= 0.0);
}

// ---------------------------------------------------------------------------
// select_ils_subset ordering
// ---------------------------------------------------------------------------

#[test]
fn candidates_are_ranked_tightest_covariance_first() {
    // variances 0.09, 0.04, 0.01 -> stdevs 0.3, 0.2, 0.1.
    // Ranking is by descending (-stdev), i.e. ascending sigma:
    //   [2, 1, 0]
    let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.09, 0.04, 0.01]));
    let (idx, sub_a, _) = select_ils_subset(&a, &q, 0.0);
    assert_eq!(idx, vec![2, 1, 0]);
    assert_eq!(sub_a[0], 3.0);
    assert_eq!(sub_a[1], 2.0);
    assert_eq!(sub_a[2], 1.0);
}

#[test]
fn subset_stops_before_the_success_rate_falls_below_target() {
    // sigma = 10 cycles gives a per-ambiguity rate of
    // 2*Phi(1/20)-1 ~ 0.04; three of them multiply to ~6e-5, so with a 0.99
    // target only the first (sigma = 0.1) survives.
    let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 100.0, 100.0]));
    let (idx, sub_a, sub_q) = select_ils_subset(&a, &q, 0.99);
    assert_eq!(idx, vec![0]);
    assert_eq!(sub_a.len(), 1);
    assert_eq!(sub_q.nrows(), 1);
    assert!((sub_q[(0, 0)] - 0.01).abs() < 1e-15);
}

#[test]
fn metadata_overrides_covariance_ranking() {
    // Index 0 has the tighter covariance but poor geometry; index 1 has the
    // looser covariance but a high-elevation, high-SNR, long-lock arc. The CQM
    // ranking must prefer index 1.
    let a = DVector::from_vec(vec![1.001, 2.001]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.0004, 0.04]));
    let meta = vec![
        AmbiguityMetadata::new(0.10, 22.0, 2.0, 3.0),
        AmbiguityMetadata::new(1.30, 48.0, 60.0, 0.1),
    ];
    let (idx, sub_a, _) = select_ils_subset_with_metadata(&a, &q, 0.0, Some(&meta));
    assert_eq!(idx, vec![1, 0], "CQM ranking must win over raw sigma");
    assert_eq!(sub_a[0], 2.001);
}

#[test]
fn missing_metadata_falls_back_per_index_not_per_call() {
    // Metadata shorter than the vector: index 1 has no entry, so it must use
    // the fallback penalty -(frac + 0.5*sigma) = -(0.001 + 0.05) = -0.051,
    // while index 0 gets its CQM credit and must therefore rank first.
    let a = DVector::from_vec(vec![1.001, 2.001]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 0.01]));
    let meta = vec![AmbiguityMetadata::new(1.30, 48.0, 60.0, 0.1)];
    let (idx, _, _) = select_ils_subset_with_metadata(&a, &q, 0.0, Some(&meta));
    assert_eq!(idx, vec![0, 1]);
}

#[test]
fn submatrix_is_exact_for_a_non_contiguous_selection() {
    let q = DMatrix::from_row_slice(4, 4, &[
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
    ]);
    let sub = submatrix(&q, &[3, 1]);
    assert_eq!(sub.nrows(), 2);
    // rows/cols [3, 1] of the 4x4 matrix
    //   [0][0] = Q[3][3] = 16 ; [0][1] = Q[3][1] = 14
    //   [1][0] = Q[1][3] =  8 ; [1][1] = Q[1][1] =  6
    assert_eq!(sub[(0, 0)], 16.0);
    assert_eq!(sub[(0, 1)], 14.0);
    assert_eq!(sub[(1, 0)], 8.0);
    assert_eq!(sub[(1, 1)], 6.0);
}

#[test]
fn ldlt_rejects_a_non_positive_definite_input_without_panicking() {
    // A negative diagonal must be clamped to the 1e-12 floor so the
    // downstream success rate stays finite.
    let q = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, 4.0]);
    let d = ldlt_diagonal(&q);
    assert!(d[0] >= 1e-12);
    assert!(d.iter().all(|v| v.is_finite()));
}