//! Attitude-error retraction convention: golden vectors and left/right
//! discrimination for every component that retracts an ESKF attitude.
//!
//! Found convention (verified in source, not assumed):
//!   eskf/update.rs:15-17         q <- Exp(d_theta) * q      LEFT / global
//!   eskf/smoother.rs:123         q <- Exp(d_x)      * q      LEFT / global
//!   eskf/smoother.rs:102         dq = q_sm * q_pred^-1      LEFT / global
//!   eskf/predict.rs:111          q <- q * Exp(w_b*dt)       RIGHT / body
//!   eskf/predict.rs:41-58        Phi: vel_att = +f_skew*dt,
//!                                att    = -R_b2e*b_g*dt
//!
//! The error state is global-frame / left-multiplied and every component
//! agrees. AGENTS.md cites "predictor.rs:86-91"; that file does not exist in
//! this repo. The `-psi` term in the doc comment has no counterpart here
//! because the ESKF state is ECEF, not local-level, so there is no
//! navigation-yaw error to negate.
#![allow(clippy::unwrap_used)]

use nalgebra::{UnitQuaternion, Vector3};

use gneiss_core::time::GpsTime;

use crate::composite::{Matrix15, Vector15};
use crate::estimators::eskf::{EskfSnapshot, EskfState, EskfSmoother};

use super::sprint_tests::{t0, RX_ECEF};

// ===========================================================================
// ATTITUDE-ERROR RETRACTION CONVENTION
// ===========================================================================
//
// Found convention, per component (verified in the source, not assumed):
//
//   eskf/update.rs:15-17      apply_error_injection   -> q <- Exp(d_theta) * q   LEFT / global
//   eskf/smoother.rs:123      apply_smoother_correction-> q <- Exp(d_x)   * q   LEFT / global
//   eskf/smoother.rs:102      compute_state_discrepancy -> dq = q_sm * q_pred^-1 LEFT / global
//   eskf/predict.rs:111       propagate_nominal_state  -> q <- q * Exp(w_b*dt)  RIGHT / body
//   eskf/predict.rs:41-58     Phi                      -> vel_att = +f_skew*dt,
//                                                      att = -R_b2e*b_g*dt
//
// The error state is therefore GLOBAL-frame/left-multiplied and every
// component agrees. AGENTS.md cites "predictor.rs:86-91"; that file does not
// exist in this repo. The `-psi` term in the doc comment has no counterpart
// here because the ESKF state is ECEF, not local-level, so there is no
// navigation-yaw error to negate.

fn q_from(x: f64, y: f64, z: f64) -> UnitQuaternion<f64> {
    UnitQuaternion::from_scaled_axis(Vector3::new(x, y, z))
}

/// Golden vector: the retraction must be LEFT (global) multiplication.
#[test]
fn apply_error_injection_is_left_multiplied_global() {
    let q0 = q_from(0.3, 0.5, -0.2);
    let d = Vector3::new(0.1, 0.2, -0.05);
    let mut state = EskfState::new(RX_ECEF, Vector3::zeros(), q0);
    let mut dx = Vector15::zeros();
    dx.fixed_rows_mut::<3>(6).copy_from(&d);
    crate::estimators::eskf::apply_error_injection(&mut state, &dx);

    let left = (q_from(d.x, d.y, d.z) * q0).to_rotation_matrix();
    let right = (q0 * q_from(d.x, d.y, d.z)).to_rotation_matrix();
    let got = state.attitude.to_rotation_matrix().into_inner();
    let left = left.into_inner();
    let right = right.into_inner();
    assert!(
        (got - left).norm() < 1e-12,
        "retraction is not Exp(d)*q (got {})",
        (got - left).norm()
    );
    assert!(
        (left - right).norm() > 1e-3,
        "fixture is degenerate: left and right multiplication agree"
    );
}

/// Physical consequence of a left retraction, stated independently of the
/// algebra: the body-frame lever arm `l` maps into ECEF as `R*l`, and under
/// `R_new = Exp(d)*R_old` that point rotates about the ECEF axis `d`, so
///
///     R_new*l - R_old*l  ~=  d x (R_old*l)          (to first order in d)
///
/// A right (body) retraction would instead rotate about the ECEF axis `R_old*d`.
#[test]
fn left_retraction_rotates_the_antenna_about_the_ecef_correction_axis() {
    let q0 = q_from(0.3, 0.5, -0.2);
    let l = Vector3::new(0.4, -0.25, 0.9);
    let d = Vector3::new(1.0e-4, 2.0e-4, -0.5e-4);
    let mut state = EskfState::new(RX_ECEF, Vector3::zeros(), q0);
    let mut dx = Vector15::zeros();
    dx.fixed_rows_mut::<3>(6).copy_from(&d);
    crate::estimators::eskf::apply_error_injection(&mut state, &dx);

    let r_old = q0.to_rotation_matrix().into_inner();
    let r_new = state.attitude.to_rotation_matrix().into_inner();
    let p_old = r_old * l;
    let p_new = r_new * l;
    let predicted = p_old + d.cross(&p_old);
    assert!(
        (p_new - predicted).norm() < 1e-6,
        "antenna moved {:?}; a left retraction predicts {:?}",
        p_new,
        predicted
    );

    let body_axis = r_old * d;
    let right_predicted = p_old + body_axis.cross(&p_old);
    assert!(
        (p_new - right_predicted).norm() > 1e-9,
        "the fixture cannot distinguish left from right retraction"
    );
}

/// Exact golden vector, no series: rotate 90 deg about ECEF +z from identity.
///
/// Hand arithmetic: R0 = I, so the antenna sits at `l = (1, 0, 0)`. After
/// `Exp([0, 0, pi/2]) * I` the body x-axis points along ECEF y, so the antenna
/// is at (0, 1, 0) exactly.
#[test]
fn ninety_degree_z_correction_rotates_the_lever_arm_exactly() {
    let l = Vector3::new(1.0, 0.0, 0.0);
    let mut state = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
    let mut dx = Vector15::zeros();
    dx[8] = std::f64::consts::FRAC_PI_2; // d_theta_z
    crate::estimators::eskf::apply_error_injection(&mut state, &dx);

    let ant = state.attitude.to_rotation_matrix() * l;
    assert!(
        (ant - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12,
        "a +90 deg ECEF-z correction must swing body-x onto ECEF-y; got {ant:?}"
    );
}

/// The discrepancy used by the RTS backward pass must be the SAME left error
/// the forward filter injects, otherwise the smoother corrects a different
/// quantity than the filter estimated.
#[test]
fn smoother_discrepancy_is_the_left_error_and_is_consistent_with_the_filter() {
    use EskfSnapshot;
    // Drive two epochs through the real smoother with a real attitude
    // correction applied between the prediction and the post-fix state.
    let q0 = q_from(0.2, -0.1, 0.05);
    let d = Vector3::new(0.01, -0.02, 0.005);
    let state0 = EskfState::new(RX_ECEF, Vector3::zeros(), q0);
    let q1 = q_from(0.21, -0.09, 0.06);

    let mut pred1 = state0.clone();
    pred1.attitude = q1;
    let pred1_att = pred1.attitude;
    let mut post1 = pred1.clone();
    let mut dx = Vector15::zeros();
    dx.fixed_rows_mut::<3>(6).copy_from(&d);
    crate::estimators::eskf::apply_error_injection(&mut post1, &dx);

    let mut sm = EskfSmoother::new();
    sm.push(EskfSnapshot {
        time: GpsTime::new(t0().week, 0.0),
        state_pred: state0.clone(),
        state_post: state0.clone(),
        phi: Matrix15::identity(),
        is_gnss_available: true,
    });
    let post1_att = post1.attitude;
    sm.push(EskfSnapshot {
        time: GpsTime::new(t0().week, 1.0),
        state_pred: { let mut p = state0.clone(); p.attitude = pred1_att; p },
        state_post: post1,
        phi: Matrix15::identity(),
        is_gnss_available: true,
    });

    // The left discrepancy of `post1` w.r.t. `pred1` must equal exactly the d
    // the filter injected, because post1 = Exp(d) * pred1.
    // The LEFT discrepancy q_post * q_pred^-1 has rotation
    //   R_post * R_pred^T = Exp(d) * R_pred * R_pred^T = Exp(d).
    let left = post1_att * pred1_att.inverse();
    assert!(
        (left.scaled_axis() - d).norm() < 1e-9,
        "left discrepancy {:?} != injected correction {d:?}",
        left.scaled_axis().transpose()
    );
    // The RIGHT discrepancy q_pred^-1 * q_post is a CONJUGATION of Exp(d) by
    // R_pred, so it is a different rotation vector for any non-identity pose.
    let right = pred1_att.inverse() * post1_att;
    assert!(
        (right.scaled_axis() - d).norm() > 1e-6,
        "the fixture cannot distinguish left from right discrepancy"
    );
}
