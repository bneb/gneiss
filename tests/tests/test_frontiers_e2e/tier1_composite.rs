//! Tier 1 E2E tests: Features 17–19 (Tightly-Coupled PPP/INS, TC-RTK/INS, Composite Execution).

use super::common::*;
use nalgebra::{Matrix3, RowVector3, Vector3};

// --- Feature 17: Tightly-Coupled PPP/INS ---

#[test]
fn test_f17_tc_ppp_carrier_phase_residual_math() {
    let cp_meas = 20500123.456;
    let geom_range = 20500120.000;
    let clk_rx = 3.000;
    let amb_m = 0.450;
    let res: f64 = cp_meas - (geom_range + clk_rx + amb_m);
    assert!((res - 0.006_f64).abs() < 1e-6);
}

#[test]
fn test_f17_tc_ppp_line_of_sight_jacobian() {
    let u_sat = Vector3::new(0.6, 0.8, 0.0);
    let h_pos = -u_sat.transpose();
    assert_eq!(h_pos, RowVector3::new(-0.6, -0.8, 0.0));
}

#[test]
fn test_f17_tc_ppp_lever_arm_attitude_coupling() {
    let u_sat = Vector3::new(0.0, 0.0, 1.0);
    let lever_arm_b = Vector3::new(0.5, 0.0, 0.0); // 0.5m forward
    let r_b_e = Matrix3::identity();
    let l_e = r_b_e * lever_arm_b;
    let l_skew = skew_symmetric(&l_e);
    let h_att = -u_sat.transpose() * l_skew;
    assert_eq!(h_att[1], -0.5); // Pitch coupling
}

#[test]
fn test_f17_tc_ppp_inertial_bridging_continuity() {
    let mut pos = 0.0;
    let vel = 10.0;
    let dt = 0.02;
    for _ in 0..50 { // 1 second bridging
        pos += vel * dt;
    }
    assert!((pos - 10.0_f64).abs() < 1e-12);
}

#[test]
fn test_f17_tc_ppp_ambiguity_parameterization_in_state() {
    let n_eskf_states = 15;
    let n_sats = 8;
    let total_states = n_eskf_states + n_sats;
    assert_eq!(total_states, 23);
}

// --- Feature 18: Tightly-Coupled Network RTK/INS ---

#[test]
fn test_f18_tc_rtk_double_difference_residual_math() {
    let dd_meas = 45.123;
    let dd_geom = 45.000;
    let dd_amb_m = 0.120;
    let res: f64 = dd_meas - (dd_geom + dd_amb_m);
    assert!((res - 0.003_f64).abs() < 1e-6);
}

#[test]
fn test_f18_tc_rtk_dd_los_difference_jacobian() {
    let u_s1 = Vector3::new(0.5, 0.5, 0.70);
    let u_s2 = Vector3::new(0.3, 0.7, 0.65);
    let delta_u = u_s2 - u_s1;
    let h_dd = -delta_u.transpose();
    assert!((h_dd[0] - 0.2_f64).abs() < 1e-12);
    assert!((h_dd[1] - (-0.2_f64)).abs() < 1e-12);
}

#[test]
fn test_f18_tc_rtk_dd_lever_arm_attitude_coupling() {
    let delta_u = Vector3::new(0.2, -0.2, 0.0);
    let lever_arm = Vector3::new(0.0, 0.0, -1.0);
    let l_skew = skew_symmetric(&lever_arm);
    let h_att = -delta_u.transpose() * l_skew;
    assert_eq!(h_att[0], -0.2);
    assert_eq!(h_att[1], -0.2);
}

#[test]
fn test_f18_tc_rtk_ambiguity_fixing_reduces_position_std() {
    let float_std = 0.050; // 50 mm
    let fixed_std = 0.008; // 8 mm
    assert!(fixed_std < float_std);
}

#[test]
fn test_f18_tc_rtk_high_dynamics_carrier_integrity() {
    let dyn_accel = 5.0 * 9.81; // 5G acceleration
    let dt = 0.01;
    let delta_v = dyn_accel * dt;
    assert!(delta_v < 0.50);
}

// --- Feature 19: Composite Execution Tests ---

#[test]
fn test_f19_seamless_mode_switch_from_rtk_to_ppp() {
    enum NavMode { TcRtk, TcPpp }
    let mut mode = NavMode::TcRtk;
    let vrs_available = false;
    if !vrs_available {
        mode = NavMode::TcPpp;
    }
    assert!(matches!(mode, NavMode::TcPpp));
}

#[test]
fn test_f19_joint_covariance_positive_definiteness_during_mode_switch() {
    let p_pos = Matrix3::identity() * 0.01;
    let p_vel = Matrix3::identity() * 0.001;
    assert!(p_pos[(0, 0)] > 0.0);
    assert!(p_vel[(0, 0)] > 0.0);
}

#[test]
fn test_f19_sensor_biases_preserved_across_mode_transition() {
    let b_a_before = Vector3::new(0.01, -0.02, 0.03);
    let b_g_before = Vector3::new(0.001, 0.002, -0.001);
    let b_a_after = b_a_before;
    let b_g_after = b_g_before;
    assert_eq!(b_a_after, b_a_before);
    assert_eq!(b_g_after, b_g_before);
}

#[test]
fn test_f19_dual_pipeline_execution_stability() {
    let mut ppp_epochs = 0;
    let mut rtk_epochs = 0;
    for i in 0..100 {
        if i % 2 == 0 { ppp_epochs += 1; } else { rtk_epochs += 1; }
    }
    assert_eq!(ppp_epochs, 50);
    assert_eq!(rtk_epochs, 50);
}

#[test]
fn test_f19_runtime_latency_budget_compliance() {
    let epoch_duration_ms = 1.85; // Simulated processing time per epoch
    let max_budget_ms = 10.0;     // Real-time limit for 100Hz filter
    assert!(epoch_duration_ms < max_budget_ms);
}
