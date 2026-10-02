#![allow(clippy::unwrap_used)]
//! Tests for the robust weighting, screening and scalar-residual filters.
//!
//! Split out of `robust.rs` to respect the repo's <500 LOC rule (this is the
//! same pattern the `swfg` modules use, e.g. `factor_tests.rs`).

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

use crate::estimators::rtk_iekf::state::{DoubleDiffKey, RtkState};

use super::robust::{
    phase_innovation_outliers, update_zwd_scalar, validate_fixed_carrier_residuals,
    validate_fixed_pseudorange_residuals, CmcTracker, MAX_ZWD_STEP_M,
};
use super::system::compute_tropo_dd;
use super::DoubleDiffMeasurement;

/// Nominal L1 wavelength: c / f1 with f1 = 1575.42 MHz (GPS L1).
const L1: f64 = SPEED_OF_LIGHT_M_S / 1575.42e6;

/// A real rover position (Tokyo area) so the troposphere model is exercised
/// at a genuine station height rather than at a degenerate origin.
const ROVER: Vector3<f64> = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);

/// Unit-ish line-of-sight directions used to build the fixture geometry.
const DIRS: [Vector3<f64>; 5] = [
    Vector3::new(0.6, 0.4, 0.7),
    Vector3::new(-0.5, 0.8, 0.3),
    Vector3::new(0.2, -0.6, 0.77),
    Vector3::new(-0.8, -0.3, 0.52),
    Vector3::new(0.35, 0.55, -0.76),
];

const RANGE_M: f64 = 2.1e7;

fn key(sat: u16, band: u8) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: band }
}

fn sat_pos(dir: Vector3<f64>) -> Vector3<f64> {
    ROVER + dir.normalize() * RANGE_M
}

/// Reference satellite straight up.
fn ref_pos() -> Vector3<f64> {
    ROVER + Vector3::new(0.0, 0.0, 2.4e7)
}

/// Base station 1 km from the rover along +x (baseline term in the gates).
fn base_pos() -> Vector3<f64> {
    ROVER + Vector3::new(1000.0, 0.0, 0.0)
}

/// Modelled double-difference geometry (m) at `pos`, including the tropo term
/// the residual gates add. Used only to build *self-consistent* observations:
/// the assertions below are about residual MAGNITUDES (one cycle = 19.029 cm),
/// never about this helper's value.
fn geom_dd(pos: Vector3<f64>, sat: Vector3<f64>, refs: Vector3<f64>, base: Vector3<f64>) -> f64 {
    let r_sat = (sat - pos).norm();
    let r_ref = (refs - pos).norm();
    let base_dd = (sat - base).norm() - (refs - base).norm();
    (r_sat - r_ref) - base_dd + compute_tropo_dd(sat, refs, base, pos)
}

fn meas(key: DoubleDiffKey, sat: Vector3<f64>, refs: Vector3<f64>, base: Vector3<f64>) -> DoubleDiffMeasurement {
    DoubleDiffMeasurement {
        key,
        dd_pr_m: 0.0,
        dd_cp_cycles: None,
        sat_pos: sat,
        ref_pos: refs,
        base_pos: base,
        lambda: L1,
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

/// Code measurement whose residual against `ROVER` is exactly `res_m` metres.
fn code_meas(sat_nr: usize, res_m: f64, pr_var_m2: f64) -> DoubleDiffMeasurement {
    let sat = sat_pos(DIRS[sat_nr]);
    let refs = ref_pos();
    let base = base_pos();
    let mut m = meas(key(2 + sat_nr as u16, 1), sat, refs, base);
    m.dd_pr_m = geom_dd(ROVER, sat, refs, base) + res_m;
    m.pr_var_m2 = pr_var_m2;
    m
}

// ---------------------------------------------------------------------------
// update_zwd_scalar
// ---------------------------------------------------------------------------

/// Closed-form single-observation ZWD update.
///   prior = var + rw*max(|dt|,1e-3);  post = 1/(1/prior + h^2/r);
///   delta = post * (h*y/r), clamped to +/- MAX_ZWD_STEP_M.
fn zwd_closed_form(var: f64, rw: f64, dt: f64, h: f64, y: f64, r: f64) -> (f64, f64) {
    let prior = var + rw * dt.abs().max(1e-3);
    let denom = 1.0 / prior + h * h / r;
    let post = 1.0 / denom;
    let delta = (post * (h * y / r)).clamp(-MAX_ZWD_STEP_M, MAX_ZWD_STEP_M);
    (delta, post)
}

#[test]
fn zwd_scalar_matches_its_closed_form_on_one_observation() {
    // var=0.0225, rw=3e-7, dt=1 -> prior = 0.0225003
    //   1/prior = 44.44329...; h=1,r=0.01 -> +100 -> denom=144.44329...
    //   post = 6.92317e-3; num = 1*0.01/0.01 = 1 -> delta = 6.92317e-3 m
    let (zwd, var) = (0.010, 0.0225);
    let (delta, post) = zwd_closed_form(var, 3e-7, 1.0, 1.0, 0.01, 0.01);
    let (got_z, got_v) = update_zwd_scalar(zwd, var, 3e-7, 1.0, &[(1.0, 0.01, 0.01)]);
    assert!((got_v - post).abs() < 1e-15, "post var {} vs {post}", got_v);
    assert!((got_z - (zwd + delta)).abs() < 1e-15, "zwd {got_z} vs {}", zwd + delta);
    // Hand check of the magnitude: 0.01 m observation shrinks the 0.0225 m^2
    // prior by ~3.4x, because the 0.01 m^2 observation dominates.
    assert!((post - 1.0 / (1.0 / 0.0225003 + 100.0)).abs() < 1e-15);
}

#[test]
fn zwd_scalar_step_saturates_at_the_documented_ceiling() {
    // A 1 m innovation with r=0.01 asks for a 1.0 m step; the per-epoch
    // saturation bound MAX_ZWD_STEP_M = 0.05 m must clip it.
    let (_, _) = update_zwd_scalar(0.0, 0.0225, 3e-7, 1.0, &[(1.0, 1.0, 0.01)]);
    let (z, _) = update_zwd_scalar(0.0, 0.0225, 3e-7, 1.0, &[(1.0, 1.0, 0.01)]);
    assert_eq!(z, MAX_ZWD_STEP_M, "step must be clipped to the ceiling exactly");
    let (z_neg, _) = update_zwd_scalar(0.0, 0.0225, 3e-7, 1.0, &[(1.0, -1.0, 0.01)]);
    assert_eq!(z_neg, -MAX_ZWD_STEP_M, "clipping is symmetric");
}

#[test]
fn zwd_scalar_ignores_pairs_with_non_positive_variance() {
    // r <= 0 rows are unusable and must not enter the information sum; with
    // only such rows the scalar degenerates to its prior.
    let pairs = [(1.0, 1.0, 0.0), (1.0, 1.0, -1.0)];
    let (z, v) = update_zwd_scalar(0.0, 0.0225, 3e-7, 1.0, &pairs);
    assert!((v - 0.0225003).abs() < 1e-15, "prior must be unchanged: {v}");
    assert_eq!(z, 0.0);
}

#[test]
fn zwd_scalar_process_noise_uses_absolute_time_step() {
    // The backward pass runs with negative dt; variance growth must mirror
    // the forward step exactly, so +5 s and -5 s give the same posterior.
    let fwd = update_zwd_scalar(0.0, 0.0225, 3e-7, 5.0, &[(1.0, 0.02, 0.01)]);
    let bwd = update_zwd_scalar(0.0, 0.0225, 3e-7, -5.0, &[(1.0, 0.02, 0.01)]);
    assert_eq!(fwd.0, bwd.0);
    assert_eq!(fwd.1, bwd.1);
}

// ---------------------------------------------------------------------------
// phase_innovation_outliers
// ---------------------------------------------------------------------------

/// State with `n` ambiguities of 1 cycle^2 variance and a tight position prior.
fn phase_state(n: usize) -> RtkState {
    let mut st = RtkState::new(ROVER, GpsTime::new(2000, 100.0));
    for i in 0..3 {
        st.cov[(i, i)] = 1.0e-4; // 1 cm position prior -> s_cp ~ 1 cycle^2
    }
    for i in 0..n {
        st.ensure_ambiguity(key(2 + i as u16, 1), 100.0 + i as f64, 1.0);
    }
    st
}

/// Phase measurement whose predicted cycle count matches the state exactly.
fn phase_meas(st: &RtkState, i: usize, slip_cycles: f64) -> DoubleDiffMeasurement {
    let sat = sat_pos(DIRS[i]);
    let refs = ref_pos();
    let base = base_pos();
    let k = key(2 + i as u16, 1);
    let amb = st.to_dvector()[st.get_amb_idx(&k).unwrap()];
    let mut m = meas(k, sat, refs, base);
    m.dd_cp_cycles = Some(amb + geom_dd(ROVER, sat, refs, base) / L1 + slip_cycles);
    m
}

#[test]
fn phase_innovation_outliers_leaves_clean_pairs_alone() {
    let st = phase_state(4);
    let ms: Vec<_> = (0..4).map(|i| phase_meas(&st, i, 0.0)).collect();
    assert!(phase_innovation_outliers(&st, &ms, 500.0).is_empty());
}

#[test]
fn phase_innovation_outliers_flags_a_single_thousand_cycle_slip() {
    let st = phase_state(4);
    let mut ms: Vec<_> = (0..4).map(|i| phase_meas(&st, i, 0.0)).collect();
    ms[2].dd_cp_cycles = Some(ms[2].dd_cp_cycles.unwrap() + 1000.0);
    let out = phase_innovation_outliers(&st, &ms, 500.0);
    assert_eq!(out, vec![key(4, 1)], "only the slipped pair may be flagged");
}

#[test]
fn phase_innovation_outliers_drops_the_verdict_when_a_majority_flips() {
    // 3 of 4 pairs disagree: a common-mode problem (reference or position
    // prior), not a per-satellite slip. Flagging more than half the pairs is
    // evidence the gate's own premise is wrong, so it must clear.
    let st = phase_state(4);
    let mut ms: Vec<_> = (0..4).map(|i| phase_meas(&st, i, 0.0)).collect();
    for m in ms.iter_mut().take(3) {
        m.dd_cp_cycles = Some(m.dd_cp_cycles.unwrap() + 1000.0);
    }
    assert!(phase_innovation_outliers(&st, &ms, 500.0).is_empty());
}

#[test]
fn phase_innovation_outliers_never_gates_an_unconverged_ambiguity() {
    // A freshly reset arc carries var >= 4 cycles^2; its phase "innovation"
    // is meaningless until it reconverges, so it must be skipped entirely.
    let mut st = phase_state(3);
    let k = key(2, 1);
    let idx = st.get_amb_idx(&k).unwrap();
    st.cov[(idx, idx)] = 5.0;
    let mut ms = vec![phase_meas(&st, 0, 0.0)];
    ms[0].dd_cp_cycles = Some(ms[0].dd_cp_cycles.unwrap() + 5000.0);
    assert!(phase_innovation_outliers(&st, &ms, 500.0).is_empty());
}

#[test]
fn phase_innovation_outliers_skips_pairs_without_carrier_phase() {
    let st = phase_state(2);
    let ms = vec![meas(key(2, 1), sat_pos(DIRS[0]), ref_pos(), base_pos())];
    assert!(phase_innovation_outliers(&st, &ms, 500.0).is_empty());
}

// ---------------------------------------------------------------------------
// validate_fixed_carrier_residuals
// ---------------------------------------------------------------------------

#[test]
fn validate_fixed_carrier_rejects_an_empty_fixed_set() {
    // No integers claimed means nothing was validated: that must NOT read as
    // a pass, otherwise a fix with zero fixed pairs would sail the gate.
    assert!(!validate_fixed_carrier_residuals(ROVER, &[], &[], 0.05));
}

#[test]
fn validate_fixed_carrier_rejects_a_key_with_no_measurement() {
    let st = phase_state(1);
    let ms = vec![phase_meas(&st, 0, 0.0)];
    let absent = key(31, 1);
    assert!(!validate_fixed_carrier_residuals(ROVER, &ms, &[(absent, 100.0)], 0.05));
}

#[test]
fn validate_fixed_carrier_rejects_a_pair_with_no_carrier_phase() {
    let ms = vec![meas(key(2, 1), sat_pos(DIRS[0]), ref_pos(), base_pos())];
    assert!(!validate_fixed_carrier_residuals(ROVER, &ms, &[(key(2, 1), 100.0)], 0.05));
}

#[test]
fn validate_fixed_carrier_boundary_is_quarter_cycles_at_zero_baseline() {
    // Gate: |cp*L - geom - N*L| <= max_res + 4e-6 * |pos - base|.
    // Put base exactly on the rover so the baseline term vanishes and
    // max_res = 0.05 m exactly. One L1 cycle = 0.19029367 m, so the
    // boundary is 0.05 / 0.19029367 = 0.262749... cycles.
    let sat = sat_pos(DIRS[0]);
    let refs = ref_pos();
    let base = ROVER;
    let amb = 100.0;
    // Sanity: base and rover coincide, so the baseline term is exactly zero and
// the effective gate is the bare 0.05 m.
assert!((ROVER - base).norm() < 1e-12);
let g = geom_dd(ROVER, sat, refs, base);
let mut m = meas(key(2, 1), sat, refs, base);
m.dd_cp_cycles = Some((g + amb * L1) / L1);

    let ok = |cycles: f64| {
        let mut mm = m.clone();
        mm.dd_cp_cycles = Some(mm.dd_cp_cycles.unwrap() + cycles);
        validate_fixed_carrier_residuals(ROVER, &[mm], &[(key(2, 1), amb)], 0.05)
    };
    assert!(ok(0.26), "0.26 cyc = 0.04948 m <= 0.05 m must pass");
    assert!(!ok(0.27), "0.27 cyc = 0.05138 m > 0.05 m must fail");
}

#[test]
fn validate_fixed_carrier_widens_the_gate_with_baseline_length() {
    // Same 0.28-cycle error, now over a 1 km baseline:
    // eff_max = 0.05 + 4e-6 * 1000 = 0.054 m, and 0.28 * 0.19029367 =
    // 0.053282 m < 0.054 m, so it must PASS here while the zero-baseline
    // gate would have rejected it. Deleting the baseline term fails this.
    let sat = sat_pos(DIRS[0]);
    let refs = ref_pos();
    let base = base_pos();
    let amb = 100.0;
    let mut m = meas(key(2, 1), sat, refs, base);
    let g = geom_dd(ROVER, sat, refs, base);
    m.dd_cp_cycles = Some((g + amb * L1) / L1 + 0.28);
    assert!(validate_fixed_carrier_residuals(ROVER, &[m], &[(key(2, 1), amb)], 0.05));
}

// ---------------------------------------------------------------------------
// validate_fixed_pseudorange_residuals
// ---------------------------------------------------------------------------

#[test]
fn validate_fixed_pseudorange_accepts_self_consistent_geometry() {
    // Residuals are exactly zero by construction, so rms = 0 and no row is
    // flagged: rms <= 4.0 and n_large = 0 <= max(5/8, 3) = 3.
    let ms: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 0.04)).collect();
    assert!(validate_fixed_pseudorange_residuals(ROVER, &ms, 4.0, 12.0));
}

#[test]
fn validate_fixed_pseudorange_rejects_a_position_displaced_by_100_m() {
    // A false integer fix in weak geometry moves the position but not the
    // observations: every residual then reads tens of metres.
    let truth = code_meas(0, 0.0, 0.04);
    let ms: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 0.04)).collect();
    let truth_dd = truth.dd_pr_m;
    let displaced = ROVER + Vector3::new(100.0, 0.0, 0.0);
    // Rebuild the observations as seen from the true position, then validate
    // them at a position 100 m away: geom at `displaced` differs from geom at
    // ROVER by roughly -(d(pos)/dpos . 100 m) ~ 100 m, so the rms gate fires.
    let ms_shift: Vec<_> = ms.iter().map(|m| {
        let mut c = m.clone();
        c.dd_pr_m = geom_dd(displaced, m.sat_pos, m.ref_pos, m.base_pos) + (truth_dd - geom_dd(ROVER, m.sat_pos, m.ref_pos, m.base_pos));
        c
    }).collect();
    assert!(!validate_fixed_pseudorange_residuals(ROVER, &ms_shift, 4.0, 12.0));
}

#[test]
fn validate_fixed_pseudorange_counts_large_residuals_against_a_fixed_budget() {
    // Four rows at 12 m over a 5-row epoch, with the hard per-row limit set
    // above the residual so the 8 m "outlier" rule does not fire:
    //   res(12) > max_pr_res_m(10)  ->  n_large += 2 each -> 8
    //   budget = max(5/8, 3) = 3                        -> 8 > 3 -> reject
    // rms = sqrt((4*144)/5) = 10.73 m, but the rms gate is opened to 100 m
    // so only the counter can be responsible.
    let mut ms: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 0.04)).collect();
    for m in ms.iter_mut().take(4) {
        m.dd_pr_m += 12.0;
    }
    assert!(!validate_fixed_pseudorange_residuals(ROVER, &ms, 100.0, 10.0));
}

#[test]
fn validate_fixed_pseudorange_keeps_sub_threshold_residuals() {
    // Control for the test above: four rows at 7 m with the limit at 10 m, so
    // neither branch of the counter fires:
    //   rms = sqrt((4*49)/5) = 6.261 m <= 100 m  and  n_large = 0 <= 3.
    let mut ms: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 0.04)).collect();
    for m in ms.iter_mut().take(4) {
        m.dd_pr_m += 7.0;
    }
    assert!(validate_fixed_pseudorange_residuals(ROVER, &ms, 100.0, 10.0));
}

#[test]
fn validate_fixed_pseudorange_outlier_test_scales_with_declared_noise() {
    // 8.5 m > 8.0 m, but 2.5*sqrt(100) = 25 m: with a 10 m-sigma code the
    // residual is NOT an outlier and must not consume the rejection budget.
    let ms: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 100.0)).collect();
    assert!(validate_fixed_pseudorange_residuals(ROVER, &ms, 100.0, 10.0));
    // The same 8.5 m against a 0.2 m-sigma code IS an outlier; four such
    // rows exceed the budget of 3.
    let mut tight: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 0.04)).collect();
    for m in tight.iter_mut().take(4) {
        m.dd_pr_m += 8.5;
    }
    assert!(!validate_fixed_pseudorange_residuals(ROVER, &tight, 100.0, 10.0));
    // Pinning the 8 m floor itself: 9 m consumes the budget, 7.9 m does not.
    let mut nine: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 0.04)).collect();
    for m in nine.iter_mut().take(4) {
        m.dd_pr_m += 9.0;
    }
    assert!(!validate_fixed_pseudorange_residuals(ROVER, &nine, 100.0, 10.0));
    let mut under: Vec<_> = (0..5).map(|i| code_meas(i, 0.0, 0.04)).collect();
    for m in under.iter_mut().take(4) {
        m.dd_pr_m += 7.9;
    }
    assert!(validate_fixed_pseudorange_residuals(ROVER, &under, 100.0, 10.0));
}

#[test]
fn validate_fixed_pseudorange_empty_input_is_vacuously_valid() {
    assert!(validate_fixed_pseudorange_residuals(ROVER, &[], 4.0, 12.0));
}

// ---------------------------------------------------------------------------
// CmcTracker
// ---------------------------------------------------------------------------

#[test]
fn cmc_tracker_reports_the_code_minus_carrier_deviation_in_metres() {
    // cmc = dd_pr_m - lambda * dd_cp_cycles. With dd_pr = 10 m and a phase
    // count that would correspond to 9.5 m of range, cmc = 10 - 9.5 = 0.5 m.
    let k = key(2, 1);
    let mut t = CmcTracker::new();
    let cp_for = |range_m: f64| range_m / L1;
    t.update_pair(k, 10.0, cp_for(9.5), L1, false);
    assert_eq!(t.get_multipath_m(&k), 0.0, "warm-up never declares multipath");
    // Five more clean epochs on a 10.0 m code-minus-carrier arc.
    for _ in 0..5 {
        t.update_pair(k, 10.0, cp_for(10.0), L1, false);
    }
    assert!(!t.is_multipath(&k), "flat arc is not multipath");
    // Now a 3.6 m jump in the CMC observable. The baseline is an exponential
    // average seeded at 0.5 and driven to 0 with alphas 1/2, 1/3, 1/4, 1/5
    // (warm-up) then 1/6 and 1/7 (post-warm-up), which telescopes to
    //   0.5 * (2/3)(3/4)(4/5)(5/6)(6/7) = 0.5 * 2/7 = 1/7
    // so the reported deviation is 3.6 - 1/7 = 3.457142857142857 m.
    let mp = t.update_pair(k, 13.6, cp_for(10.0), L1, false);
    assert!((mp - (3.6 - 1.0 / 7.0)).abs() < 1e-12, "reported multipath {mp}");
    assert!(t.is_multipath(&k));
}

#[test]
fn cmc_tracker_retain_active_drops_pairs_that_disappeared() {
    let (k1, k2) = (key(2, 1), key(3, 1));
    let mut t = CmcTracker::new();
    t.update_pair(k1, 10.0, 10.0 / L1, L1, false);
    t.update_pair(k2, 20.0, 20.0 / L1, L1, false);
    t.retain_active(&[k1]);
    assert!(!t.is_multipath(&k2) && t.get_multipath_m(&k2) == 0.0);
    t.update_pair(k2, 30.0, 30.0 / L1, L1, false);
    assert_eq!(t.get_multipath_m(&k2), 0.0, "dropped pair restarts clean");
}

#[test]
fn cmc_tracker_slip_flag_restarts_the_arc() {
    let k = key(2, 1);
    let mut t = CmcTracker::new();
    for _ in 0..10 {
        t.update_pair(k, 10.0, 10.0 / L1, L1, false);
    }
    let mp = t.update_pair(k, 40.0, 10.0 / L1, L1, true);
    assert_eq!(mp, 0.0, "a declared slip is not multipath evidence");
    assert!(!t.is_multipath(&k));
}

