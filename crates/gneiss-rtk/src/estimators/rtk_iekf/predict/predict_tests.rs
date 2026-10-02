//! Tests for state/covariance time propagation, in particular the process
//! noise injected into every optional state channel.
//!
//! A child module of `predict` so the private Q builder is reachable.

use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

use super::*;
use crate::estimators::rtk_iekf::state::{DoubleDiffKey, RtkState};

const DT: f64 = 10.0;

fn key(sat: u16) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: 1 }
}

fn base_state() -> RtkState {
    RtkState::new(Vector3::zeros(), GpsTime::new(2200, 0.0))
}

fn q_of(state: &RtkState, dt: f64, reverse_safe: bool) -> DMatrix<f64> {
    build_process_noise(state, state.dim(), dt, 1.0, reverse_safe)
}

// ---------------------------------------------------------------------------
// Transition matrix
// ---------------------------------------------------------------------------

#[test]
fn the_transition_matrix_is_the_documented_constant_velocity_block() {
    let f = build_transition_matrix(6, DT);
    for i in 0..3 {
        assert_eq!(f[(i, i + 3)], DT, "position is advanced by dt * velocity");
        for c in 0..6 {
            if c != i && c != i + 3 {
                assert_eq!(f[(i, c)], if i == c { 1.0 } else { 0.0 }, "row {i} col {c}");
            }
        }
    }
    for i in 3..6 {
        for c in 0..6 {
            assert_eq!(f[(i, c)], if i == c { 1.0 } else { 0.0 }, "velocity row {i} col {c}");
        }
    }
}

#[test]
fn an_optional_state_column_is_carried_through_unchanged() {
    // With ZWD enabled the ambiguity block starts one column later; F must
    // still leave that column at identity (the state is not propagated).
    let mut st = base_state();
    st.enable_zwd(0.0225);
    st.ensure_ambiguity(key(2), 1.0, 4.0);
    let f = build_transition_matrix(st.dim(), DT);
    let ai = st.get_amb_idx(&key(2)).unwrap();
    assert_eq!(f[(ai, ai)], 1.0, "ambiguity column must be the identity");
    assert_eq!(f.nrows(), st.dim());
}

// ---------------------------------------------------------------------------
// Process noise per channel
// ---------------------------------------------------------------------------

#[test]
fn a_constant_velocity_kinematic_noise_scales_as_dt_powers() {
    let q = q_of(&base_state(), DT, false);
    // q_pos = dt^3/3 * q_accel = 1000/3 = 333.333...
    assert!((q[(0, 0)] - DT * DT * DT / 3.0).abs() < 1e-12);
    // q_pos_vel = dt^2/2 * q_accel = 50
    assert!((q[(0, 3)] - DT * DT / 2.0).abs() < 1e-12);
    // q_vel = dt * q_accel = 10
    assert!((q[(3, 3)] - DT).abs() < 1e-12);
    // Cross term is symmetric (forward step).
    assert_eq!(q[(3, 0)], q[(0, 3)]);
    // No noise leaks into a constant-velocity pair that shares no state.
    assert_eq!(q[(0, 4)], 0.0);
}

#[test]
fn the_zwd_and_gradient_states_get_their_documented_random_walk_rates() {
    let mut st = base_state();
    st.enable_zwd(0.0225);
    st.enable_gradients(4.0e-6);
    let zi = st.zwd_idx().expect("zwd enabled");
    let (gn, ge) = st.grad_idx().expect("gradients enabled");
    let q = q_of(&st, DT, false);
    assert!((q[(zi, zi)] - super::super::update::ZWD_RW_M2_PER_S * DT).abs() < 1e-18);
    let want = super::super::update::GRAD_RW_M2_PER_S * DT;
    assert!((q[(gn, gn)] - want).abs() < 1e-20);
    assert!((q[(ge, ge)] - want).abs() < 1e-20);
}

#[test]
fn ambiguity_variance_grows_at_its_documented_rate() {
    let mut st = base_state();
    st.ensure_ambiguity(key(2), 1.0, 4.0);
    st.ensure_ambiguity(key(3), 2.0, 4.0);
    let a0 = st.get_amb_idx(&key(2)).unwrap();
    let a1 = st.get_amb_idx(&key(3)).unwrap();
    let q = q_of(&st, DT, false);
    // 1e-7 cycles^2 per second of tracking: at 10 s that is 1e-6 cycles^2.
    assert!((q[(a0, a0)] - 1e-7 * DT).abs() < 1e-18);
    assert!((q[(a1, a1)] - 1e-7 * DT).abs() < 1e-18);
    // Ambiguity process noise never crosses into the position block.
    assert_eq!(q[(a0, 0)], 0.0);
}

#[test]
fn per_pair_and_per_satellite_iono_states_random_walk_too() {
    let mut st = base_state();
    st.ensure_ambiguity(key(2), 1.0, 4.0);
    st.iono_enabled = true;
    st.ensure_iono(key(2), 4.0);
    let ii = st.get_iono_idx(&key(2)).expect("iono state present");
    let q = q_of(&st, DT, false);
    // IONO_RW_M2_PER_S = 1e-8 m^2/s -> 1e-7 at a 10 s step, and the slot
    // starts at zero because no other term writes there.
    assert!((q[(ii, ii)] - 1e-8 * DT).abs() < 1e-18, "q = {}", q[(ii, ii)]);

    let mut st2 = base_state();
    st2.ensure_ambiguity(key(2), 1.0, 4.0);
    st2.sat_iono_enabled = true;
    st2.ensure_sat_iono_key(0, 2);
    st2.ensure_sat_iono_key(0, 1);
    let si = st2.get_sat_iono_key_idx(0, 2).unwrap();
    let q2 = q_of(&st2, DT, false);
    assert!((q2[(si, si)] - 1e-8 * DT).abs() < 1e-18, "sat iono q = {}", q2[(si, si)]);
}

// ---------------------------------------------------------------------------
// Reverse-safe handling
// ---------------------------------------------------------------------------

#[test]
fn a_backward_step_mirrors_the_forward_variances_under_reverse_safe() {
    let st = base_state();
    let back = q_of(&st, -DT, true);
    let fwd = q_of(&st, DT, true);
    for i in 0..6 {
        for j in 0..6 {
            assert!((back[(i, j)] - fwd[(i, j)]).abs() < 1e-15, "Q[{i}][{j}] must mirror");
        }
    }
}

#[test]
fn the_legacy_backward_path_still_produces_negative_variances() {
    // Documented legacy behaviour, kept only so the corrected path stays
    // observable: a negative dt poisons the position variance directly.
    let st = base_state();
    let legacy = q_of(&st, -DT, false);
    assert!(legacy[(0, 0)] < 0.0);
    assert!(legacy[(3, 3)] < 0.0);
    assert!(legacy[(0, 3)] < 0.0, "only the cross term keeps the legacy sign");
}

// ---------------------------------------------------------------------------
// Full propagation
// ---------------------------------------------------------------------------

#[test]
fn propagation_keeps_the_covariance_symmetric_and_positive_definite() {
    let mut st = base_state();
    st.enable_zwd(0.0225);
    st.ensure_ambiguity(key(2), 1.0, 4.0);
    st.ensure_ambiguity(key(3), 2.0, 4.0);
    let before = st.cov.clone();
    predict_state_gated(&mut st, GpsTime::new(2200, DT), 1.0, true);
    for i in 0..st.dim() {
        for j in 0..st.dim() {
            assert!(
                (st.cov[(i, j)] - st.cov[(j, i)]).abs() < 1e-12,
                "P[{i}][{j}] asymmetry"
            );
        }
    }
    let e = nalgebra::linalg::SymmetricEigen::new(st.cov.clone()).eigenvalues;
    assert!(e.min() > 0.0, "P must stay positive definite, min eig {}", e.min());
    assert!(e.max() > nalgebra::linalg::SymmetricEigen::new(before).eigenvalues.max());
}

#[test]
fn a_zero_velocity_state_leaves_the_position_untouched() {
    let mut st = base_state();
    st.pos_ecef = Vector3::new(1.0e6, 2.0e6, 3.0e6);
    st.vel_ecef = Vector3::zeros();
    let before = st.pos_ecef;
    predict_state_gated(&mut st, GpsTime::new(2200, 3600.0), 1.0, false);
    assert!((st.pos_ecef - before).norm() < 1e-12, "a static receiver must not drift");
    assert_eq!(st.time.tow, 3600.0);
}

#[test]
fn a_moving_state_advances_by_exactly_velocity_times_dt() {
    let mut st = base_state();
    st.vel_ecef = Vector3::new(1.0, -2.0, 0.5);
    predict_state(&mut st, GpsTime::new(2200, 2.0), 1.0);
    // 2 s at (1, -2, 0.5) m/s.
    assert!((st.pos_ecef - Vector3::new(2.0, -4.0, 1.0)).norm() < 1e-12);
}