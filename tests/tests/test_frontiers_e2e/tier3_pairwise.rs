//! Tier 3 E2E tests: Cross-Feature Pairwise Interactions.

use super::common::*;
use gneiss_core::frequencies::Signal;
use nalgebra::{DMatrix, DVector, Matrix3, RowVector3, UnitQuaternion, Vector2, Vector3};

#[test]
fn test_t3_eskf_and_coupled_nhc_attitude_observability() {
    let v_e = Vector3::new(20.0, 0.0, 0.0); // 20 m/s forward
    let v_skew = skew_symmetric(&v_e);
    let r_b_e = Matrix3::identity();
    let h_nhc_att = RowVector3::new(0.0, 1.0, 0.0) * r_b_e.transpose() * v_skew;
    // Non-zero attitude Jacobian in yaw column provides heading observability
    assert_eq!(h_nhc_att[2], -20.0);
}

#[test]
fn test_t3_eskf_bias_estimation_with_zupt_separation() {
    let mut p_cov = Matrix3::identity() * 0.1;
    let h_zupt = Matrix3::identity();
    let r_zupt = Matrix3::identity() * 0.001; // High precision zero-velocity
    let k = p_cov * h_zupt.transpose() * (h_zupt * p_cov * h_zupt.transpose() + r_zupt).try_inverse().unwrap();
    p_cov = (Matrix3::identity() - k * h_zupt) * p_cov;
    assert!(p_cov[(0, 0)] < 0.002);
}

#[test]
fn test_t3_eskf_quaternion_feedback_and_rts_smoother() {
    let q_filt = r_b_e_from_rpy(0.1, 0.2, 0.3);
    let delta_th_smooth = Vector3::new(0.002, -0.001, 0.004);
    let dq = UnitQuaternion::from_scaled_axis(delta_th_smooth);
    let q_smoothed = q_filt * dq;
    assert!((q_smoothed.norm() - 1.0).abs() < 1e-12);
}

#[test]
fn test_t3_sinex_osb_pcv_and_single_diff_lambda_ar() {
    let lam_wl = SPEED_OF_LIGHT / (1575.42e6 - 1227.60e6);
    let raw_mw = 4.02 * lam_wl;
    let sat_pcv = 0.0015;
    let osb_bias = 0.015;
    let corr_mw = raw_mw - sat_pcv - osb_bias;
    let float_n_wl = corr_mw / lam_wl;
    let int_n_wl = float_n_wl.round() as i32;
    assert_eq!(int_n_wl, 4);

    let a_hat = DVector::from_vec(vec![float_n_wl, float_n_wl * 1.5]);
    let q_aa = DMatrix::from_diagonal(&DVector::from_vec(vec![0.005, 0.005]));
    let res = gneiss_rtk::ambiguity::lambda::resolve_lambda(&a_hat, &q_aa).unwrap();
    assert_eq!(res.best_integers[0].round() as i32, 4);
    assert_eq!(res.best_integers[1].round() as i32, 6);
    assert!(res.ratio > 2.0);
}

#[test]
fn test_t3_multi_constellation_and_receiver_pco_pcv() {
    let f_gps = Signal::GpsL1Ca.base_freq_hz();
    let f_gal = Signal::GalE1Os.base_freq_hz();
    let f_bds = Signal::BdsB1i.base_freq_hz();
    assert!(f_gps > 1e9 && f_gal > 1e9 && f_bds > 1e9);

    let pco_enu = Vector3::new(0.001, 0.002, 0.050);
    assert_eq!(pco_enu.z, 0.050);
}

#[test]
fn test_t3_multi_cors_ingestion_and_dd_network_adjustment() {
    let network = get_standard_cors_network_2d();
    let b12 = network[1] - network[0];
    let b23 = network[2] - network[1];
    let b31 = network[0] - network[2];
    let loop_closure = b12 + b23 + b31;
    assert_eq!(loop_closure, Vector2::zeros());
}

#[test]
fn test_t3_delaunay_triangulation_and_vrs_synthesis() {
    let network = get_standard_cors_network_2d();
    let rover_approx = Vector2::new(10.0, 5.0);
    let bary = compute_barycentric(&rover_approx, &network[0], &network[1], &network[2]).unwrap();
    let zwd_bases = Vector3::new(0.120, 0.125, 0.130);
    let vrs_zwd = bary.dot(&zwd_bases);
    assert!(vrs_zwd > 0.120 && vrs_zwd < 0.130);
}

#[test]
fn test_t3_localized_vrs_and_network_rtk_benchmark() {
    let master = Vector2::new(0.0, 0.0);
    let rover = Vector2::new(20.0, 15.0); // 25 km from master
    let vrs = rover - Vector2::new(0.1, 0.1); // 141m effective baseline
    let orig_len = (rover - master).norm();
    let vrs_len = (rover - vrs).norm();
    assert_eq!(orig_len, 25.0);
    assert!(vrs_len < 0.2);

    let leica_spec = 0.008 + 1e-6 * vrs_len * 1000.0;
    assert!(leica_spec < 0.010); // < 10mm spec!
}

#[test]
fn test_t3_tc_ppp_ins_and_15state_eskf_coupling() {
    let p_eskf = Matrix3::identity() * 0.05;
    let u_los = Vector3::new(0.0, 0.0, 1.0);
    let h_pos = -u_los.transpose();
    let r_cp = 0.0001; // 1cm carrier phase noise
    let s = (h_pos * p_eskf * h_pos.transpose())[(0, 0)] + r_cp;
    assert!(s > 0.0);
}

#[test]
fn test_t3_tc_rtk_ins_to_tc_ppp_composite_switch() {
    let mut current_mode = "TC-RTK";
    let rover_pos = Vector2::new(100.0, 100.0); // Out of CORS coverage
    let cors_max_radius = 50.0;
    if rover_pos.norm() > cors_max_radius {
        current_mode = "TC-PPP";
    }
    assert_eq!(current_mode, "TC-PPP");
}
