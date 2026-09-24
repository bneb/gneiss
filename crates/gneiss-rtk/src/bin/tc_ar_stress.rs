//! Adversarial empirical stress testing harness for TC-AR.
//!
//! Evaluates:
//! 1. Ill-conditioned Q_aa (cond > 10^10).
//! 2. High correlation between ambiguity states (rho -> 1.0).
//! 3. Injected cycle slips, false candidate vectors, and tracker Q_aa post-fix definiteness.
//! 4. Injected 10-meter position jump gating and float state preservation.
//! 5. 10,000 Monte Carlo perturbations verifying post-condition P_check positive definiteness.

use nalgebra::{DMatrix, DVector, SymmetricEigen, UnitQuaternion, Vector3};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use gneiss_rtk::composite::tc_ambiguity::TcAmbiguityTracker;
use gneiss_rtk::estimators::eskf::condition::apply_integer_conditioning;
use gneiss_rtk::estimators::eskf::{EskfState, Matrix15};
use gneiss_rtk::estimators::rtk_iekf::DoubleDiffKey;

fn make_base_eskf(pos: Vector3<f64>) -> EskfState {
    let mut cov = Matrix15::zeros();
    for i in 0..3 { cov[(i, i)] = 1.0; }
    for i in 3..6 { cov[(i, i)] = 0.1; }
    for i in 6..9 { cov[(i, i)] = 0.01; }
    for i in 9..12 { cov[(i, i)] = 0.04; }
    for i in 12..15 { cov[(i, i)] = 1e-4; }
    EskfState::with_cov(pos, Vector3::zeros(), UnitQuaternion::identity(), Vector3::zeros(), Vector3::zeros(), cov)
}

fn split_joint_blocks(joint: &DMatrix<f64>, n_amb: usize) -> (Matrix15<f64>, DMatrix<f64>, DMatrix<f64>) {
    let mut p_xx = Matrix15::zeros();
    let mut p_xa = DMatrix::zeros(15, n_amb);
    let mut q_aa = DMatrix::zeros(n_amb, n_amb);
    for i in 0..15 {
        for j in 0..15 { p_xx[(i, j)] = joint[(i, j)]; }
        for j in 0..n_amb { p_xa[(i, j)] = joint[(i, 15 + j)]; }
    }
    for i in 0..n_amb {
        for j in 0..n_amb { q_aa[(i, j)] = joint[(15 + i, 15 + j)]; }
    }
    (p_xx, p_xa, q_aa)
}

fn generate_joint_spd(n_amb: usize, rng: &mut StdRng) -> (Matrix15<f64>, DMatrix<f64>, DMatrix<f64>) {
    let tot = 15 + n_amb;
    let mut rand_mat = DMatrix::zeros(tot, tot);
    for i in 0..tot {
        for j in 0..tot { rand_mat[(i, j)] = rng.random_range(-1.0..1.0); }
    }
    let mut joint = &rand_mat * rand_mat.transpose();
    for i in 0..tot { joint[(i, i)] += rng.random_range(0.5..5.0); }
    split_joint_blocks(&joint, n_amb)
}

fn test_ill_conditioned_q_aa() -> bool {
    println!("--- Test 1: Ill-conditioned Q_aa (cond > 10^10) ---");
    let cond_targets = [1e10, 1e11, 1e12, 1e14, 1e16];
    for &target in &cond_targets {
        let mut st = make_base_eskf(Vector3::zeros());
        let prior_st = st.clone();
        let eps = 1.0 / target;
        let q_aa = DMatrix::from_row_slice(2, 2, &[1.0, 1.0 - eps, 1.0 - eps, 1.0]);
        let mut p_xa = DMatrix::zeros(15, 2);
        p_xa[(0, 0)] = 0.01;
        p_xa[(1, 1)] = 0.01;
        let a_flt = DVector::from_column_slice(&[1.05, 2.05]);
        let a_fix = DVector::from_column_slice(&[1.0, 2.0]);
        match apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix) {
            Ok(summary) => {
                println!("  target cond={:.1e} => Inversion succeeded, accepted={}, dx_norm={:.3e}",
                    target, summary.accepted, summary.dx.norm());
                if !summary.accepted {
                    assert_eq!(st.pos_ecef, prior_st.pos_ecef);
                    assert_eq!(st.cov, prior_st.cov);
                }
            }
            Err(e) => {
                println!("  target cond={:.1e} => Graceful inversion error: {:?}", target, e);
                assert_eq!(st.pos_ecef, prior_st.pos_ecef);
                assert_eq!(st.cov, prior_st.cov);
            }
        }
    }
    println!("Test 1 PASSED: Graceful handling under extreme condition numbers.\n");
    true
}

fn test_high_ambiguity_correlation() -> bool {
    println!("--- Test 2: High correlation between ambiguity states ---");
    let rhos = [0.90, 0.99, 0.999, 0.9999, 0.99999, 0.999999];
    for &rho in &rhos {
        let mut st = make_base_eskf(Vector3::zeros());
        let q_aa = DMatrix::from_row_slice(2, 2, &[0.04, 0.04 * rho, 0.04 * rho, 0.04]);
        let mut p_xa = DMatrix::zeros(15, 2);
        p_xa[(0, 0)] = 0.005;
        p_xa[(0, 1)] = 0.005 * rho;
        p_xa[(1, 0)] = 0.005 * rho;
        p_xa[(1, 1)] = 0.005;
        let a_flt = DVector::from_column_slice(&[1.02, 2.02]);
        let a_fix = DVector::from_column_slice(&[1.0, 2.0]);
        let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix)
            .expect("conditioning should complete");
        let min_eig = SymmetricEigen::new(st.cov).eigenvalues.min();
        println!("  rho={:.6} => applied={}, accepted={}, min_eigenvalue={:.6e}",
            rho, res.applied, res.accepted, min_eig);
        assert!(min_eig >= 1e-10, "Covariance lost positive definiteness!");
    }
    println!("Test 2 PASSED: Stable conditioning under extreme cross-correlation.\n");
    true
}

fn test_cycle_slip_and_false_integers() -> bool {
    println!("--- Test 3: Injected cycle slips and false integer candidate vectors ---");
    let mut tracker = TcAmbiguityTracker::new();
    let keys: Vec<DoubleDiffKey> = (2..=6).map(|prn| DoubleDiffKey {
        constellation_id: 0, sat: prn, ref_sat: 1, freq_band: 1,
    }).collect();
    tracker.sync_keys(&keys, &[1.0, 2.0, 3.0, 4.0, 5.0]);
    let mut state = make_base_eskf(Vector3::new(100.0, 200.0, 300.0));
    let prior_state = state.clone();

    let (fixed, _) = tracker.attempt_ar_and_condition(&mut state, |_, key, _| {
        if key.sat == 2 { Some(0.1903) } else { Some(0.005) }
    });
    println!("  Injected 1-cycle slip on Sat 2: fixed={}", fixed);
    assert!(!fixed, "False fix accepted despite 1-cycle carrier residual blunder!");
    assert_eq!(state.pos_ecef, prior_state.pos_ecef, "State modified despite blunder!");
    assert_eq!(state.cov, prior_state.cov, "Covariance modified despite blunder!");

    let (fixed_5, _) = tracker.attempt_ar_and_condition(&mut state, |_, key, _| {
        if key.sat == 3 { Some(0.9515) } else { Some(0.002) }
    });
    println!("  Injected 5-cycle slip on Sat 3: fixed={}", fixed_5);
    assert!(!fixed_5, "False fix accepted on 5-cycle slip!");
    assert_eq!(state.pos_ecef, prior_state.pos_ecef);
    assert_eq!(state.cov, prior_state.cov);

    println!("Test 3 PASSED: Carrier residual gate strictly repels slipped ambiguities.\n");
    true
}

fn test_tracker_post_fix_q_aa_definiteness() -> bool {
    println!("--- Test 3b: Tracker Q_aa positive definiteness after partial/full fix ---");
    let mut tracker = TcAmbiguityTracker::new();
    let keys: Vec<DoubleDiffKey> = (2..=5).map(|prn| DoubleDiffKey {
        constellation_id: 0, sat: prn, ref_sat: 1, freq_band: 1,
    }).collect();
    tracker.sync_keys(&keys, &[0.0, 0.0, 0.0, 0.0]);
    for i in 0..4 {
        for j in 0..4 { tracker.q_aa[(i, j)] = if i == j { 0.08 } else { 0.04 }; }
    }
    let eig_before = SymmetricEigen::new(tracker.q_aa.clone()).eigenvalues.min();
    println!("  Q_aa min eigenvalue before fix: {:.6e}", eig_before);
    assert!(eig_before > 0.0);

    let mut state = make_base_eskf(Vector3::zeros());
    let (fixed, _) = tracker.attempt_ar_and_condition(&mut state, |_, _, _| Some(0.001));
    println!("  Fix attempted: fixed={}", fixed);
    if fixed {
        let eig_after = SymmetricEigen::new(tracker.q_aa.clone()).eigenvalues.min();
        println!("  Q_aa min eigenvalue after fix: {:.6e}", eig_after);
        let (fixed2, _) = tracker.attempt_ar_and_condition(&mut state, |_, _, _| Some(0.001));
        println!("  Epoch 2 fix attempt with corrupted Q_aa: fixed={}", fixed2);
        if eig_after < 0.0 {
            println!("  [EMPIRICAL FINDING]: Q_aa has negative eigenvalue {:.6e} after fix!", eig_after);
            return false;
        }
    }
    true
}

fn test_unscreened_negative_eigenvalue_acceptance() -> bool {
    println!("--- Test 3c: Condition.rs acceptance of negative eigenvalue P_check ---");
    let mut st = make_base_eskf(Vector3::zeros());
    let mut p_xa = DMatrix::zeros(15, 1);
    p_xa[(0, 0)] = 0.20;
    p_xa[(1, 0)] = 0.20;
    let q_aa = DMatrix::from_element(1, 1, 0.01);
    let a_flt = DVector::from_element(1, 1.001);
    let a_fix = DVector::from_element(1, 1.0);

    let res = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix)
        .expect("conditioning returns summary");
    let min_eig = SymmetricEigen::new(st.cov).eigenvalues.min();
    println!("  Conditioning applied={}, accepted={}, min_eigenvalue={:.6e}",
        res.applied, res.accepted, min_eig);
    if res.accepted && min_eig < 0.0 {
        println!("  [EMPIRICAL FINDING]: Condition.rs ACCEPTED invalid P_check with negative eigenvalue {:.6e}!", min_eig);
        return false;
    }
    true
}

fn test_10m_position_jump_gating() -> bool {
    println!("--- Test 4: Injected 10-meter position jump gating ---");
    let mut st = make_base_eskf(Vector3::new(123.456, -789.012, 345.678));
    let prior_st = st.clone();

    let mut p_xa = DMatrix::zeros(15, 1);
    p_xa[(0, 0)] = 0.50; // dx_pos = -P_xa/Q * da = -0.5 / 0.05 * 1.0 = -10.0m
    let q_aa = DMatrix::from_element(1, 1, 0.05);
    let a_flt = DVector::from_element(1, 1.0);
    let a_fix = DVector::from_element(1, 0.0);

    let sigma_3d = (st.cov[(0, 0)] + st.cov[(1, 1)] + st.cov[(2, 2)]).sqrt();
    let gate = (3.0 * sigma_3d).max(0.50);
    println!("  sigma_3d={:.3}m, gate={:.3}m, injected ||dx_p||=10.0m", sigma_3d, gate);

    let summary = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix)
        .expect("gating should return summary");
    println!("  Jump gate result: applied={}, accepted={}, dx[0]={:.2}m",
        summary.applied, summary.accepted, summary.dx[0]);

    assert!(!summary.applied, "10-meter jump was applied!");
    assert!(!summary.accepted, "10-meter jump was accepted!");
    assert_eq!(st.pos_ecef, prior_st.pos_ecef, "Position modified despite 10m jump!");
    assert_eq!(st.vel_ecef, prior_st.vel_ecef, "Velocity modified!");
    assert_eq!(st.attitude, prior_st.attitude, "Attitude modified!");
    assert_eq!(st.accel_bias, prior_st.accel_bias, "Accel bias modified!");
    assert_eq!(st.gyro_bias, prior_st.gyro_bias, "Gyro bias modified!");
    assert_eq!(st.cov, prior_st.cov, "Covariance modified!");

    println!("Test 4 PASSED: 10-meter jump strictly gated, float state perfectly preserved.\n");
    true
}

fn execute_single_mc_trial(
    rng: &mut StdRng,
    min_eig_overall: &mut f64,
    max_eig_overall: &mut f64,
    non_spd_count: &mut usize,
) -> bool {
    let n_amb = rng.random_range(1..=6);
    let (p_xx, p_xa, q_aa) = generate_joint_spd(n_amb, rng);
    let mut st = EskfState::with_cov(
        Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity(),
        Vector3::zeros(), Vector3::zeros(), p_xx,
    );
    let mut a_flt = DVector::zeros(n_amb);
    let mut a_fix = DVector::zeros(n_amb);
    for i in 0..n_amb {
        let int_val = rng.random_range(-10..10) as f64;
        a_fix[i] = int_val;
        a_flt[i] = int_val + rng.random_range(-0.15..0.15);
    }
    if let Ok(summary) = apply_integer_conditioning(&mut st, &p_xa, &q_aa, &a_flt, &a_fix) {
        if summary.accepted {
            let eig = SymmetricEigen::new(st.cov).eigenvalues;
            let m_eig = eig.min();
            if m_eig < *min_eig_overall { *min_eig_overall = m_eig; }
            if eig.max() > *max_eig_overall { *max_eig_overall = eig.max(); }
            if m_eig <= 0.0 { *non_spd_count += 1; }
            return true;
        }
    }
    false
}

fn test_10k_monte_carlo_positive_definiteness() -> bool {
    println!("--- Test 5: 10,000 Randomized Monte Carlo Perturbations ---");
    let mut rng = StdRng::seed_from_u64(0x474e45495353); // "GNEISS"
    let mut min_eig_overall = f64::INFINITY;
    let mut max_eig_overall = f64::NEG_INFINITY;
    let mut non_spd_count = 0usize;
    let mut accepted_count = 0usize;

    for trial in 0..10_000 {
        if execute_single_mc_trial(&mut rng, &mut min_eig_overall, &mut max_eig_overall, &mut non_spd_count) {
            accepted_count += 1;
        }
        if (trial + 1) % 2500 == 0 {
            println!("  Completed {}/10,000 trials (accepted={}, non_spd={})...",
                trial + 1, accepted_count, non_spd_count);
        }
    }

    println!("  Summary: 10,000 trials, accepted={}, non_spd_count={}", accepted_count, non_spd_count);
    println!("  Eigenvalue range of P_check: [{:.6e}, {:.6e}]", min_eig_overall, max_eig_overall);
    assert_eq!(non_spd_count, 0, "Found non-positive definite P_check in Monte Carlo!");
    assert!(min_eig_overall >= 1e-10, "Minimum eigenvalue violates floor!");
    println!("Test 5 PASSED: Post-condition covariance is strictly positive definite across all trials.\n");
    true
}

fn print_verdict(results: &[(&str, bool)]) -> bool {
    println!("============================================================");
    println!("FINAL VERDICT SUMMARY:");
    let mut all_pass = true;
    for (name, passed) in results {
        println!("  {:<38} {}", name, if *passed { "PASS" } else { "FAIL" });
        if !*passed { all_pass = false; }
    }
    println!("============================================================");
    all_pass
}

fn main() {
    println!("============================================================");
    println!("    ADVERSARIAL STRESS TESTING: TC-AR IMPLEMENTATION        ");
    println!("============================================================\n");

    let results = [
        ("Test 1 (Ill-conditioned Q_aa > 10^10):", test_ill_conditioned_q_aa()),
        ("Test 2 (High correlation states):", test_high_ambiguity_correlation()),
        ("Test 3 (Cycle slips & false integers):", test_cycle_slip_and_false_integers()),
        ("Test 3b (Tracker Q_aa definiteness):", test_tracker_post_fix_q_aa_definiteness()),
        ("Test 3c (Condition.rs P_check SPD gate):", test_unscreened_negative_eigenvalue_acceptance()),
        ("Test 4 (10m position jump gate):", test_10m_position_jump_gating()),
        ("Test 5 (10,000 Monte Carlo SPD):", test_10k_monte_carlo_positive_definiteness()),
    ];

    if !print_verdict(&results) {
        std::process::exit(1);
    }
}
