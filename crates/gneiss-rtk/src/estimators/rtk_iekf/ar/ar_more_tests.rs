//! Additional tests for the AR entry points and the integer-conditioning
//! guard rails. Split out because `ar.rs` is already at the file-size limit.

use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

use super::*;
use crate::estimators::rtk_iekf::update::DoubleDiffMeasurement;

fn key(sat: u16) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: 1 }
}

/// State with `n` near-integer ambiguities and the given position variance.
fn state_with(n: usize, amb_var: f64, pos_var: f64) -> RtkState {
    let mut st = RtkState::new(Vector3::new(1.0e6, 2.0e6, 3.0e6), GpsTime::new(2200, 100.0));
    for i in 0..3 {
        st.cov[(i, i)] = pos_var;
    }
    for i in 0..n {
        st.ensure_ambiguity(key(2 + i as u16), 10.0 + i as f64 + 0.001, amb_var);
    }
    st
}

/// DD measurement on the same geometry, carrying `cp` carrier cycles.
fn dd_meas_for(k: DoubleDiffKey, cp: f64, base_pos: Vector3<f64>, pos: Vector3<f64>) -> DoubleDiffMeasurement {
    let sat = pos + Vector3::new(1.0e7, 0.4e7, 1.8e7);
    let refs = pos + Vector3::new(-0.6e7, 1.9e7, 1.1e7);
    DoubleDiffMeasurement {
        key: k,
        dd_pr_m: 0.0,
        dd_cp_cycles: Some(cp),
        sat_pos: sat,
        ref_pos: refs,
        base_pos,
        lambda: 0.19029367,
        pr_var_m2: 0.04,
        cp_var_cycles2: 1.0e-4,
        pr_ref_var_m2: 0.02,
        cp_ref_var_cycles2: 5.0e-5,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m: 0.0,
    }
}

// ---------------------------------------------------------------------------
// Entry-point gating
// ---------------------------------------------------------------------------

#[test]
fn a_loose_float_position_blocks_resolution_entirely() {
    // max_float_trace = 25 m^2 (static) / 12 m^2 (kinematic). A trace of
    // 300 m^2 is a 10 m position uncertainty, which no honest ratio test can
    // rescue, so AR must not even try.
    let st = state_with(6, 0.01, 100.0);
    let res = resolve_ambiguities(&st, 3, 0.001, false);
    assert!(!res.is_fixed);
    assert_eq!(res.num_ambiguities, 6, "the float report still carries the count");
    assert_eq!(res.ratio, 0.0);
    assert!(res.fixed_ambiguities.is_empty());
}

#[test]
fn the_kinematic_profile_uses_the_tighter_float_trace_budget() {
    // trace = 7 + 7 + 6 = 20 m^2: inside the static budget of 25, outside the
    // kinematic budget of 12.
    let mut st = state_with(8, 0.0004, 7.0);
    st.cov[(2, 2)] = 6.0;
    assert!(resolve_ambiguities(&st, 3, 0.001, false).is_fixed);
    assert!(!resolve_ambiguities(&st, 3, 0.001, true).is_fixed);
}

#[test]
fn three_clean_ambiguities_are_enough_to_claim_a_fix() {
    let st = state_with(3, 0.0004, 0.04);
    let res = resolve_ambiguities(&st, 3, 0.001, false);
    assert!(res.is_fixed, "n_amb == 3 clears the floor");
    assert_eq!(res.num_ambiguities, 3);
    // Every fixed integer must be the nearest integer to its float.
    for (k, n) in &res.fixed_ambiguities {
        let i = st.get_amb_idx(k).unwrap() - st.amb_offset();
        assert!((n - st.ambiguities[i].1.round()).abs() < 1e-9);
    }
}

#[test]
fn floats_sitting_on_the_lattice_midpoint_are_never_fixed() {
    // a_hat = n + 0.5 exactly: the nearest integer is a perfect tie, so the
    // best and second-best lattice points are equidistant and the LAMBDA
    // ratio is 1 whatever the covariance. No threshold above 1 can clear it,
    // so a 50/50 candidate must never become a fix.
    let mut st = state_with(8, 1.0e6, 0.04);
    for i in 0..8 {
        st.ambiguities[i].1 = 10.0 + i as f64 + 0.5;
    }
    let res = resolve_ambiguities(&st, 3, 0.001, false);
    assert!(!res.is_fixed, "an equidistant candidate must not be fixed");
    assert!(res.fixed_ambiguities.is_empty());
}

#[test]
fn a_loose_arc_is_excluded_from_the_kinematic_fix_set() {
    // try_full_ar refuses outright when any kinematic ambiguity still carries
    // more than 1 cycle^2, so the seven tight arcs must NOT all be claimed.
    // PAR is allowed to salvage them, and the fixed count proves the loose
    // arc was dropped rather than silently conditioned.
    let mut st = state_with(8, 0.0004, 0.04);
    let loose = key(5);
    let li = st.get_amb_idx(&loose).unwrap();
    st.cov[(li, li)] = 9.0;
    let res = resolve_ambiguities(&st, 3, 0.001, true);
    assert!(res.is_fixed, "seven converged arcs still support a PAR fix");
    assert_eq!(res.num_ambiguities, 7);
    assert!(
        res.fixed_ambiguities.iter().all(|(k, _)| *k != loose),
        "the 9-cycle^2 arc must never be part of a kinematic fix"
    );
}

// ---------------------------------------------------------------------------
// Carrier-residual veto (the zero-false-fix guard)
// ---------------------------------------------------------------------------

#[test]
fn a_one_cycle_carrier_error_vetoes_the_fix() {
    // Build the DD carrier observations so that a *consistent* fix exists,
    // then shift one pair by exactly one L1 cycle (0.19029367 m). The
    // post-fix carrier residual gate is 0.05 m, so a 19 cm discrepancy must
    // veto the whole fix even though LAMBDA itself was happy.
    let st = state_with(6, 0.0004, 0.04);
    let pos = st.pos_ecef;
    let base = pos + Vector3::new(300.0, -100.0, 50.0);
    let build = |shift_cycles: f64| -> Vec<DoubleDiffMeasurement> {
        (0..6).map(|i| {
            let k = key(2 + i as u16);
            let m = dd_meas_for(k, 0.0, base, pos);
            // residual = |cp*lambda - geom - n*lambda|, so the carrier count
            // that reproduces integer n exactly is (geom + n*lambda)/lambda.
            let (a, b, c) = (m.sat_pos, m.ref_pos, base);
            let t = crate::estimators::rtk_iekf::update::compute_tropo_dd(a, b, c, pos);
            let geom = (a - pos).norm() - (b - pos).norm() - ((a - c).norm() - (b - c).norm()) + t;
            let n = 10.0 + i as f64;
            let extra = if i == 3 { shift_cycles } else { 0.0 };
            let mut m2 = dd_meas_for(k, 0.0, base, pos);
            m2.dd_cp_cycles = Some((geom + n * m2.lambda) / m2.lambda + extra);
            m2
        }).collect()
    };
    let consistent = build(0.0);
    assert!(
        resolve_ambiguities_screened(&st, 3, 0.001, false, Some(&consistent)).is_fixed,
        "a self-consistent carrier set must fix"
    );
    let slipped = build(1.0);
    let res = resolve_ambiguities_screened(&st, 3, 0.001, false, Some(&slipped));
    assert!(!res.is_fixed, "a one-cycle carrier slip must veto the fix");
}

#[test]
fn a_fix_must_be_corroborated_by_carrier_observations() {
    // A DD set with code only cannot support an integer fix: the post-fix
    // carrier residual gate has no carrier row to check and refuses, so the
    // whole epoch stays float. That is deliberate -- an integer claimed from
    // code alone would be unfalsifiable.
    let st = state_with(6, 0.0004, 0.04);
    let pos = st.pos_ecef;
    let base = pos + Vector3::new(300.0, -100.0, 50.0);
    let code_only: Vec<DoubleDiffMeasurement> = (0..6)
        .map(|i| {
            let mut m = dd_meas_for(key(2 + i as u16), 0.0, base, pos);
            m.dd_cp_cycles = None;
            m
        })
        .collect();
    let res = resolve_ambiguities_screened(&st, 3, 0.001, false, Some(&code_only));
    assert!(!res.is_fixed, "code alone cannot validate an integer fix");
    // With no DD set at all the same floats do fix: the gate is only
    // restrictive when measurements are supplied to check against.
    assert!(resolve_ambiguities_screened(&st, 3, 0.001, false, None).is_fixed);
}

// ---------------------------------------------------------------------------
// condition_state_on_integers
// ---------------------------------------------------------------------------

#[test]
fn conditioning_on_an_unknown_key_is_a_no_op() {
    let mut st = state_with(3, 0.0004, 0.04);
    let before = st.pos_ecef;
    assert!(!condition_state_on_integers(&mut st, &[(key(99), 10.0)]));
    assert!((st.pos_ecef - before).norm() < 1e-15);
}

#[test]
fn conditioning_refuses_a_jump_beyond_the_two_metre_absolute_bound() {
    // Hand-computed: with Q_aa = 4 the correction is
    //   dx = P_xa * (1/4) * (a_hat - N)
    // A float of 10.5 conditioned to N = 14.0 gives da = -3.5 cycles, so
    // dx = [2, 0, 3] * 0.25 * (-3.5) = [-1.75, 0, -2.625] with
    // |dx| = 3.15 m. The 3-sigma test is wide here (sigma_3d = 5.39 m), so it
    // is the 2 m absolute bound that must reject.
    let mut st = RtkState::new(Vector3::new(1.0, 2.0, 3.0), GpsTime::new(2200, 100.0));
    let k = key(2);
    st.ensure_ambiguity(k, 10.5, 4.0);
    let off = st.amb_offset();
    for (r, v) in [4.0, 9.0, 16.0].iter().enumerate() {
        st.cov[(r, r)] = *v;
    }
    st.cov[(0, off)] = 2.0;
    st.cov[(2, off)] = 3.0;
    st.cov[(off, 0)] = 2.0;
    st.cov[(off, 2)] = 3.0;
    st.cov[(off, off)] = 4.0;
    let before = st.pos_ecef;
    assert!(!condition_state_on_integers(&mut st, &[(k, 14.0)]));
    assert!((st.pos_ecef - before).norm() < 1e-15, "a rejected fix must not move the state");
}

#[test]
fn conditioning_on_the_value_already_held_is_a_pure_variance_reduction() {
    let mut st = RtkState::new(Vector3::new(1.0, 2.0, 3.0), GpsTime::new(2200, 100.0));
    let k = key(2);
    st.ensure_ambiguity(k, 10.0, 4.0);
    let off = st.amb_offset();
    st.cov[(0, off)] = 2.0;
    st.cov[(off, 0)] = 2.0;
    st.cov[(0, 0)] = 4.0;
    st.cov[(off, off)] = 4.0;
    let before = st.pos_ecef;
    assert!(condition_state_on_integers(&mut st, &[(k, 10.0)]));
    assert!((st.pos_ecef - before).norm() < 1e-12);
    // P_xx|a = P_xx - P_xa Q_aa^-1 P_ax = 4 - 2^2/4 = 4 - 1 = 3.
    assert!((st.cov[(0, 0)] - 3.0).abs() < 1e-12, "P_xx = {}", st.cov[(0, 0)]);
}

// ---------------------------------------------------------------------------
// Hysteresis
// ---------------------------------------------------------------------------

#[test]
fn hysteresis_counts_a_stable_set_and_resets_on_the_first_change() {
    let mut consecutive = 2u32;
    let mut last: Vec<(DoubleDiffKey, f64)> = vec![(key(2), 10.0), (key(3), -5.0)];
    // A tiny change of 1e-4 cycles is below the 1e-3 tolerance.
    update_fix_hysteresis(&mut consecutive, &mut last, &[(key(2), 10.0001), (key(3), -5.0)]);
    assert_eq!(consecutive, 3);
    // A 0.01-cycle change is unambiguously outside the 1e-3 tolerance.
    // (Note 10.001 - 10.0 evaluates to 0.0009999999999998899 in binary64,
    // so it would NOT trip the test; do not use it as a boundary case.)
    update_fix_hysteresis(&mut consecutive, &mut last, &[(key(2), 10.01), (key(3), -5.0)]);
    assert_eq!(consecutive, 1);
    assert_eq!(last, vec![(key(2), 10.01), (key(3), -5.0)]);
}

#[test]
fn hysteresis_from_an_empty_history_always_starts_at_one() {
    let mut consecutive = 9u32;
    let mut last: Vec<(DoubleDiffKey, f64)> = Vec::new();
    update_fix_hysteresis(&mut consecutive, &mut last, &[(key(2), 10.0)]);
    assert_eq!(consecutive, 1);
}

#[test]
fn an_empty_current_set_is_vacuously_consistent() {
    // `current.iter().all(..)` over an empty slice is true, so an empty fixed
    // set never contradicts the previous one and the run continues. This is
    // unreachable from the engine (the caller only invokes this on a fixed
    // result, which always carries at least one integer), so it is pinned here
    // to document the vacuous-truth behaviour rather than relied upon.
    let mut consecutive = 5u32;
    let mut last: Vec<(DoubleDiffKey, f64)> = vec![(key(2), 10.0)];
    update_fix_hysteresis(&mut consecutive, &mut last, &[]);
    assert_eq!(consecutive, 6);
    assert!(last.is_empty(), "but the history is cleared for the next epoch");
}