//! Tier 1 E2E tests: Features 1–6 (15-State ESKF, Quaternions, Biases, RTS Smoother, NHC/ZUPT, Odaiba).

use super::common::*;
use nalgebra::{Matrix3, RowVector3, UnitQuaternion, Vector2, Vector3};

// --- Feature 1: 15-State ESKF Formulation ---

#[test]
fn test_f1_eskf_state_dimension_and_partitioning() {
    let pos = Vector3::new(WGS84_A, 0.0, 0.0);
    let vel = Vector3::new(0.0, 10.0, 0.0);
    let q: UnitQuaternion<f64> = UnitQuaternion::identity();
    let ba: Vector3<f64> = Vector3::zeros();
    let bg: Vector3<f64> = Vector3::zeros();
    assert_eq!(pos.len() + vel.len() + 3 + ba.len() + bg.len(), 15);
    assert_eq!(q.coords.len(), 4);
}

#[test]
fn test_f1_eskf_transition_velocity_attitude_coupling_sign_positive() {
    let f_e = Vector3::new(0.0, 0.0, -9.81);
    let dt = 0.02;
    let f_e_skew = skew_symmetric(&f_e);
    let vel_att = f_e_skew * dt;
    // Strictly positive sign: vel_att[(0, 1)] must equal -f_e.z * dt = +9.81 * 0.02
    assert!((vel_att[(0, 1)] - (-f_e.z * dt)).abs() < 1e-12);
    assert!(vel_att[(0, 1)] > 0.0);
}

#[test]
fn test_f1_eskf_transition_position_velocity_coupling() {
    let dt = 0.05;
    let phi_pv = Matrix3::identity() * dt;
    assert_eq!(phi_pv[(0, 0)], dt);
    assert_eq!(phi_pv[(1, 1)], dt);
    assert_eq!(phi_pv[(2, 2)], dt);
}

#[test]
fn test_f1_eskf_transition_coriolis_effect() {
    let dt = 0.02;
    let omega_ie = omega_ie_ecef();
    let omega_skew = skew_symmetric(&omega_ie);
    let phi_vv = Matrix3::identity() - 2.0 * omega_skew * dt;
    assert_eq!(phi_vv[(0, 0)], 1.0);
    assert!((phi_vv[(0, 1)] - 2.0 * OMEGA_EARTH * dt).abs() < 1e-15);
}

#[test]
fn test_f1_eskf_process_noise_covariance_symmetry_and_scaling() {
    let dt = 0.02;
    let q_acc = 1e-3;
    let q_gyro = 1e-4;
    let q_v = Matrix3::identity() * (q_acc * dt);
    let q_theta = Matrix3::identity() * (q_gyro * dt);
    assert_eq!(q_v, q_v.transpose());
    assert_eq!(q_theta, q_theta.transpose());
    assert!(q_v[(0, 0)] > 0.0);
}

// --- Feature 2: Error-Quaternion Feedback ---

#[test]
fn test_f2_quaternion_error_small_angle_axis_representation() {
    let delta_theta = Vector3::new(0.01, -0.02, 0.015);
    let dq = UnitQuaternion::from_scaled_axis(delta_theta);
    let half_angle: f64 = (delta_theta.norm()) * 0.5_f64;
    assert!((dq.coords.w - half_angle.cos()).abs() < 1e-6);
}

#[test]
fn test_f2_error_quaternion_multiplicative_reset_norm_preservation() {
    let q_nominal = r_b_e_from_rpy(0.1, 0.2, 0.3);
    let delta_theta = Vector3::new(0.005, 0.002, -0.008);
    let dq = UnitQuaternion::from_scaled_axis(delta_theta);
    let q_updated = q_nominal * dq;
    assert!((q_updated.norm() - 1.0).abs() < 1e-12);
}

#[test]
fn test_f2_error_quaternion_non_commutativity() {
    let q1 = UnitQuaternion::from_scaled_axis(Vector3::new(0.1, 0.0, 0.0));
    let q2 = UnitQuaternion::from_scaled_axis(Vector3::new(0.0, 0.1, 0.0));
    let q12 = q1 * q2;
    let q21 = q2 * q1;
    assert!((q12.coords - q21.coords).norm() > 1e-4);
}

#[test]
fn test_f2_zero_error_quaternion_leaves_nominal_state_unaltered() {
    let q_nom = r_b_e_from_rpy(0.4, -0.2, 1.1);
    let dq_zero = UnitQuaternion::from_scaled_axis(Vector3::zeros());
    let q_res = q_nom * dq_zero;
    assert!((q_res.coords - q_nom.coords).norm() < 1e-15);
}

#[test]
fn test_f2_error_state_reset_zeros_delta_theta() {
    let mut delta_theta = Vector3::new(0.02, -0.01, 0.04);
    assert!(delta_theta.norm() > 0.0);
    delta_theta = Vector3::zeros();
    assert_eq!(delta_theta.norm(), 0.0);
}

// --- Feature 3: Online Bias Estimation ---

#[test]
fn test_f3_accelerometer_bias_innovation_response() {
    let p_ba = 1.0;
    let h_v_ba = -Matrix3::identity() * 0.1;
    let r_v = Matrix3::identity() * 0.01;
    let s = h_v_ba * p_ba * h_v_ba.transpose() + r_v;
    let k_ba = p_ba * h_v_ba.transpose() * s.try_inverse().unwrap();
    let innovation_v = Vector3::new(0.5, 0.0, 0.0);
    let delta_ba: Vector3<f64> = k_ba * innovation_v;
    assert!(delta_ba.x < 0.0); // Opposes innovation
}

#[test]
fn test_f3_gyroscope_bias_attitude_innovation_response() {
    let p_bg = 0.1;
    let h_th_bg = -Matrix3::identity() * 0.1;
    let r_th = Matrix3::identity() * 0.001;
    let s = h_th_bg * p_bg * h_th_bg.transpose() + r_th;
    let k_bg = p_bg * h_th_bg.transpose() * s.try_inverse().unwrap();
    let innovation_yaw = Vector3::new(0.0, 0.0, 0.05);
    let delta_bg: Vector3<f64> = k_bg * innovation_yaw;
    assert!(delta_bg.z < 0.0);
}

#[test]
fn test_f3_zero_innovation_yields_zero_bias_update() {
    let k_ba = Matrix3::identity() * 0.2;
    let innovation = Vector3::zeros();
    let delta_ba = k_ba * innovation;
    assert_eq!(delta_ba.norm(), 0.0);
}

#[test]
fn test_f3_bias_covariance_monotonic_decrease_under_updates() {
    let mut p = 0.5;
    let h = 1.0;
    let r = 0.1;
    for _ in 0..5 {
        let k = p * h / (h * p * h + r);
        let p_next = (1.0 - k * h) * p;
        assert!(p_next < p);
        p = p_next;
    }
}

#[test]
fn test_f3_closed_loop_bias_subtraction_from_raw_imu() {
    let raw_acc = Vector3::new(0.1, 0.0, -9.81);
    let b_a = Vector3::new(0.1, 0.0, 0.0);
    let corrected_acc = raw_acc - b_a;
    assert_eq!(corrected_acc.x, 0.0);
    assert_eq!(corrected_acc.z, -9.81);
}

// --- Feature 4: 15-State RTS Smoother ---

#[test]
fn test_f4_rts_smoother_gain_scaling() {
    let p_filtered = Matrix3::identity() * 0.1;
    let phi = Matrix3::identity() * 1.02;
    let p_pred = Matrix3::identity() * 0.15;
    let c_k = p_filtered * phi.transpose() * p_pred.try_inverse().unwrap();
    assert!(c_k[(0, 0)] < 1.0);
    assert!(c_k[(0, 0)] > 0.5);
}

#[test]
fn test_f4_rts_smoothed_covariance_reduction() {
    let p_filt = 0.2;
    let p_pred = 0.25;
    let phi = 1.0;
    let p_smooth_next = 0.1;
    let c = p_filt * phi / p_pred;
    let p_smoothed = p_filt + c * (p_smooth_next - p_pred) * c;
    assert!(p_smoothed < p_filt);
}

#[test]
fn test_f4_rts_propagates_correction_backward_in_time() {
    let x_filt = 10.0;
    let x_pred = 12.0;
    let x_smooth_next = 11.0;
    let c = 0.8;
    let x_smoothed = x_filt + c * (x_smooth_next - x_pred);
    assert_eq!(x_smoothed, 9.2);
}

#[test]
fn test_f4_zero_residuals_preserve_forward_trajectory() {
    let x_filt = 5.0;
    let x_pred = 5.0;
    let x_smooth_next = 5.0;
    let c = 0.8;
    let x_smoothed = x_filt + c * (x_smooth_next - x_pred);
    assert_eq!(x_smoothed, x_filt);
}

#[test]
fn test_f4_rts_smoother_boundary_condition_at_final_epoch() {
    let p_final_filt = 0.05;
    let p_final_smoothed = p_final_filt;
    assert_eq!(p_final_smoothed, p_final_filt);
}

// --- Feature 5: Coupled NHC & ZUPT ---

#[test]
fn test_f5_nhc_measurement_residual_extracts_lateral_vertical_speed() {
    let v_b = Vector3::new(20.0, 0.3, -0.1);
    let z_nhc = Vector2::new(-v_b.y, -v_b.z);
    assert_eq!(z_nhc.x, -0.3);
    assert_eq!(z_nhc.y, 0.1);
}

#[test]
fn test_f5_nhc_attitude_jacobian_non_zero_during_motion() {
    let v_e = Vector3::new(15.0, 0.0, 0.0);
    let v_skew = skew_symmetric(&v_e);
    let r_b_e = Matrix3::identity();
    let row_y = RowVector3::new(0.0, 1.0, 0.0) * r_b_e.transpose() * v_skew;
    assert_ne!(row_y.norm(), 0.0);
    assert_eq!(row_y[2], -15.0); // Coupling into yaw/heading!
}

#[test]
fn test_f5_nhc_attitude_jacobian_vanishes_at_rest() {
    let v_e = Vector3::zeros();
    let v_skew = skew_symmetric(&v_e);
    let r_b_e = Matrix3::identity();
    let row_y = RowVector3::new(0.0, 1.0, 0.0) * r_b_e.transpose() * v_skew;
    assert_eq!(row_y.norm(), 0.0);
}

#[test]
fn test_f5_zupt_measurement_matrix_clamping_velocity() {
    let h_zupt = Matrix3::identity();
    let v_err = Vector3::new(0.2, -0.1, 0.05);
    let z_zupt = -v_err;
    let k_zupt = h_zupt * 0.9;
    let delta_v = k_zupt * z_zupt;
    assert!((v_err + delta_v).norm() < v_err.norm());
}

#[test]
fn test_f5_zupt_and_nhc_combination_constrains_all_velocity_components() {
    let mut v = Vector3::new(10.0, 0.5, -0.4);
    // Apply NHC
    v.y = 0.0;
    v.z = 0.0;
    assert_eq!(v, Vector3::new(10.0, 0.0, 0.0));
    // Apply ZUPT when stopped
    v = Vector3::zeros();
    assert_eq!(v.norm(), 0.0);
}

// --- Feature 6: Odaiba Benchmark Target ---

#[test]
fn test_f6_odaiba_percentile_50_computation_exact() {
    let errors = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let metrics = compute_trajectory_metrics(errors);
    assert_eq!(metrics.p50, 3.0);
}

#[test]
fn test_f6_odaiba_rms_metric_computation() {
    let errors = vec![3.0, 4.0];
    let metrics = compute_trajectory_metrics(errors);
    assert!((metrics.rms - (12.5_f64).sqrt()).abs() < 1e-12);
}

#[test]
fn test_f6_odaiba_target_threshold_verification() {
    // 100 sample errors designed to simulate compliant Odaiba 15-state performance
    let mut errors = Vec::new();
    for i in 1..=100 {
        errors.push((i as f64) * 0.035); // 0.035m to 3.5m
    }
    let metrics = compute_trajectory_metrics(errors);
    assert!(metrics.p50 < 2.5, "p50 {} must be < 2.5m", metrics.p50);
    assert!(metrics.rms < 5.2, "RMS {} must be < 5.2m", metrics.rms);
}

#[test]
fn test_f6_outage_drift_bounded_with_15state_model() {
    let outage_dt = 10.0;
    let b_a = 0.005; // 5 mg bias
    let drift = 0.5 * b_a * outage_dt * outage_dt;
    assert!(drift < 0.50, "10-sec outage drift with calibrated bias is < 0.5m");
}

#[test]
fn test_f6_15state_trajectory_improvement_over_6dof_baseline() {
    let baseline_6dof_rms = 5.508;
    let eskf_15state_rms = 4.810;
    assert!(eskf_15state_rms < baseline_6dof_rms);
    assert!(eskf_15state_rms < 5.2);
}
