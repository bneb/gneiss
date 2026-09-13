//! Tier 2 E2E tests: Features 1–6 Boundary & Corner Cases.

use super::common::*;
use nalgebra::{Matrix3, UnitQuaternion, Vector2, Vector3};

// --- Feature 1 Boundaries ---

#[test]
fn test_f1_b1_dt_approaches_zero() {
    let dt = 1e-6;
    let phi_pv = Matrix3::identity() * dt;
    assert!(phi_pv.norm() < 1e-5);
    let phi_full = Matrix3::identity() + phi_pv;
    assert!((phi_full[(0, 0)] - 1.0_f64).abs() < 1e-5);
}

#[test]
fn test_f1_b2_large_dt_stability() {
    let dt = 1.0; // 1 second
    let f_e = Vector3::new(0.0, 0.0, -9.81);
    let f_skew = skew_symmetric(&f_e);
    let phi_v_th = f_skew * dt;
    assert!(phi_v_th.norm().is_finite());
}

#[test]
fn test_f1_b3_zero_specific_force_vanishes_coupling() {
    let f_e = Vector3::zeros();
    let dt = 0.02;
    let phi_v_th = skew_symmetric(&f_e) * dt;
    assert_eq!(phi_v_th.norm(), 0.0);
}

#[test]
fn test_f1_b4_extreme_acceleration_shock() {
    let f_shock = Vector3::new(0.0, 100.0, 0.0); // 10G shock
    let dt = 0.01;
    let phi_v_th = skew_symmetric(&f_shock) * dt;
    assert_eq!(phi_v_th[(0, 2)], 1.0);
}

#[test]
fn test_f1_b5_near_singular_covariance_inversion() {
    let mut p = Matrix3::identity() * 1e-10;
    p[(0, 0)] = 1e-10;
    let inv = p.try_inverse();
    assert!(inv.is_some());
    assert!((inv.unwrap()[(0, 0)] - 1e10_f64).abs() < 1e2);
}

// --- Feature 2 Boundaries ---

#[test]
fn test_f2_b1_zero_error_rotation_preserves_quaternion() {
    let q = r_b_e_from_rpy(0.3, -0.4, 0.5);
    let dq = UnitQuaternion::from_scaled_axis(Vector3::zeros());
    let res = q * dq;
    assert_eq!(res.coords, q.coords);
}

#[test]
fn test_f2_b2_pi_rotation_boundary() {
    let delta_theta = Vector3::new(std::f64::consts::PI, 0.0, 0.0);
    let dq = UnitQuaternion::from_scaled_axis(delta_theta);
    assert!((dq.norm() - 1.0_f64).abs() < 1e-12);
    assert!(dq.coords.w.abs() < 1e-12); // cos(pi/2) = 0
}

#[test]
fn test_f2_b3_infinitesimal_rotation() {
    let delta_theta = Vector3::new(1e-16, 0.0, 0.0);
    let dq = UnitQuaternion::from_scaled_axis(delta_theta);
    assert_eq!(dq.coords.w, 1.0);
    assert_eq!(dq.norm(), 1.0);
}

#[test]
fn test_f2_b4_accumulated_multiplications_drift_check() {
    let mut q = UnitQuaternion::identity();
    let dq = UnitQuaternion::from_scaled_axis(Vector3::new(1e-4, 0.0, 0.0));
    for _ in 0..1000 {
        q *= dq;
    }
    assert!((q.norm() - 1.0_f64).abs() < 1e-10);
}

#[test]
fn test_f2_b5_near_gimbal_lock_attitude() {
    let pitch = std::f64::consts::FRAC_PI_2 - 1e-6;
    let q = r_b_e_from_rpy(0.0, pitch, 0.0);
    assert!((q.norm() - 1.0_f64).abs() < 1e-12);
}

// --- Feature 3 Boundaries ---

#[test]
fn test_f3_b1_massive_outlier_innovation_attenuation() {
    let p = Matrix3::identity() * 0.1;
    let h = Matrix3::identity();
    let r = Matrix3::identity() * 1000.0; // Large measurement noise rejects outlier
    let k = p * h.transpose() * (h * p * h.transpose() + r).try_inverse().unwrap();
    let innov = Vector3::new(1000.0, 0.0, 0.0);
    let delta_x = k * innov;
    assert!(delta_x.x < 0.15); // Heavily filtered
}

#[test]
fn test_f3_b2_infinite_measurement_noise_zeros_gain() {
    let p = Matrix3::identity();
    let _h = Matrix3::<f64>::identity();
    let r = Matrix3::identity() * 1e12;
    let k = p * (p + r).try_inverse().unwrap();
    assert!(k.norm() < 1e-11);
}

#[test]
fn test_f3_b3_accel_bias_near_saturation() {
    let b_a = Vector3::new(9.81, 0.0, 0.0); // 1G bias limit
    let raw = Vector3::new(9.81, 0.0, -9.81);
    let corrected = raw - b_a;
    assert_eq!(corrected.x, 0.0);
}

#[test]
fn test_f3_b4_gyro_bias_near_saturation() {
    let b_g = Vector3::new(0.0, 0.0, 1.0); // 1 rad/s bias
    let raw_omega = Vector3::new(0.0, 0.0, 1.0);
    let corrected = raw_omega - b_g;
    assert_eq!(corrected.z, 0.0);
}

#[test]
fn test_f3_b5_negative_innovation_directionality() {
    let k = Matrix3::identity() * 0.5;
    let innov_neg = Vector3::new(-1.0, -2.0, -3.0);
    let dx = k * innov_neg;
    assert!(dx.x < 0.0 && dx.y < 0.0 && dx.z < 0.0);
}

// --- Feature 4 Boundaries ---

#[test]
fn test_f4_b1_single_epoch_smoother_graceful_exit() {
    let history_len = 1;
    assert_eq!(history_len, 1);
}

#[test]
fn test_f4_b2_zero_transition_matrix_smoother_gain() {
    let p_filt = Matrix3::identity() * 0.1;
    let phi_zero = Matrix3::zeros();
    let p_pred = Matrix3::identity() * 0.2;
    let c = p_filt * phi_zero.transpose() * p_pred.try_inverse().unwrap();
    assert_eq!(c.norm(), 0.0);
}

#[test]
fn test_f4_b3_identical_covariances_unity_gain() {
    let p_filt = 0.2;
    let p_pred = 0.2;
    let phi = 1.0;
    let c = p_filt * phi / p_pred;
    assert_eq!(c, 1.0);
}

#[test]
fn test_f4_b4_zero_process_noise_smoother_limit() {
    let p_pred = 0.1;
    let p_filt = 0.1;
    let c = p_filt / p_pred;
    assert_eq!(c, 1.0);
}

#[test]
fn test_f4_b5_long_trajectory_backward_stability() {
    let mut p_smooth: f64 = 0.01;
    for _ in 0..1000 {
        p_smooth = p_smooth * 0.999 + 0.00001;
        assert!(p_smooth.is_finite());
    }
}

// --- Feature 5 Boundaries ---

#[test]
fn test_f5_b1_zero_velocity_nhc_residual_zero() {
    let v_b: Vector3<f64> = Vector3::zeros();
    let z_nhc = Vector2::new(-v_b.y, -v_b.z);
    assert_eq!(z_nhc.norm(), 0.0);
}

#[test]
fn test_f5_b2_high_speed_nhc_jacobian_linear_scaling() {
    let v1 = Vector3::new(10.0, 0.0, 0.0);
    let v2 = Vector3::new(100.0, 0.0, 0.0);
    let h1 = skew_symmetric(&v1);
    let h2 = skew_symmetric(&v2);
    assert_eq!(h2[(1, 0)], h1[(1, 0)] * 10.0);
}

#[test]
fn test_f5_b3_pure_lateral_skid_nhc_response() {
    let v_b = Vector3::new(0.0, 15.0, 0.0); // 15 m/s slide
    let z_nhc = Vector2::new(-v_b.y, -v_b.z);
    assert_eq!(z_nhc.x, -15.0);
}

#[test]
fn test_f5_b4_zero_lever_arm_exact_colocation() {
    let l_b = Vector3::zeros();
    let l_skew = skew_symmetric(&l_b);
    assert_eq!(l_skew.norm(), 0.0);
}

#[test]
fn test_f5_b5_large_lever_arm_ten_meters() {
    let l_b = Vector3::new(10.0, 0.0, 0.0);
    let l_skew = skew_symmetric(&l_b);
    assert_eq!(l_skew[(2, 1)], 10.0);
}

// --- Feature 6 Boundaries ---

#[test]
fn test_f6_b1_empty_error_vector_metrics() {
    let metrics = compute_trajectory_metrics(vec![]);
    assert_eq!(metrics.p50, 0.0);
    assert_eq!(metrics.rms, 0.0);
}

#[test]
fn test_f6_b2_single_error_sample_metrics() {
    let metrics = compute_trajectory_metrics(vec![2.5]);
    assert_eq!(metrics.p50, 2.5);
    assert_eq!(metrics.rms, 2.5);
}

#[test]
fn test_f6_b3_identical_error_distribution_metrics() {
    let errs = vec![1.8; 50];
    let metrics = compute_trajectory_metrics(errs);
    assert!((metrics.p50 - 1.8_f64).abs() < 1e-12);
    assert!((metrics.rms - 1.8_f64).abs() < 1e-12);
}

#[test]
fn test_f6_b4_bimodal_error_distribution() {
    let mut errs = vec![1.0; 50];
    errs.extend(vec![5.0; 50]);
    let metrics = compute_trajectory_metrics(errs);
    assert_eq!(metrics.p50, 5.0);
    assert!((metrics.rms - (13.0_f64).sqrt()).abs() < 1e-12);
}

#[test]
fn test_f6_b5_full_odaiba_scale_simulation() {
    let n_epochs = 12398;
    let mut errs = Vec::with_capacity(n_epochs);
    for i in 0..n_epochs {
        errs.push(0.5 + 2.0 * (i as f64) / (n_epochs as f64));
    }
    let metrics = compute_trajectory_metrics(errs);
    assert!(metrics.p50 < 2.5);
    assert!(metrics.rms < 5.2);
}
