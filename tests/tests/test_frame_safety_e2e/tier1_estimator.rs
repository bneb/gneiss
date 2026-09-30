//! Tier 1 E2E tests: Estimator state safety, lever arm rotation, and pipeline typestates.

use super::common::*;
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

// --- Feature 14: Typed ESKF Estimator State ---

#[test]
fn test_f14_eskf_15state_partitioning_and_dimension() {
    let pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let vel: SpatialVelocity<Ecef<Wgs84>> = SpatialVelocity::from_vector(Vector3::new(0.0, 15.0, 0.0));
    let q: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(UnitQuaternion::identity());
    let ba: SpatialVector<BodyFrd> = SpatialVector::zero();
    let bg: SpatialVector<BodyFrd> = SpatialVector::zero();
    let dim = pos.coords().len() + vel.vector().len() + 3 + ba.vector().len() + bg.vector().len();
    assert_eq!(dim, 15);
    assert!((q.quaternion().norm() - 1.0).abs() < 1e-12);
}

#[test]
fn test_f14_attitude_quaternion_normalization() {
    let q_raw = r_b_e_from_rpy(0.2, -0.3, 1.4);
    assert!((q_raw.norm() - 1.0).abs() < 1e-12);
}

#[test]
fn test_f14_body_biases_typed_in_body_frd() {
    let ba: SpatialVector<BodyFrd> = SpatialVector::new(0.05, -0.02, 0.01);
    let bg: SpatialVector<BodyFrd> = SpatialVector::new(1e-4, -2e-4, 5e-5);
    assert_eq!(ba.vector().x, 0.05);
    assert_eq!(bg.vector().y, -2e-4);
}

#[test]
fn test_f14_initial_covariance_symmetry_and_positive_diagonal() {
    let mut cov = nalgebra::OMatrix::<f64, nalgebra::U15, nalgebra::U15>::zeros();
    for i in 0..15 {
        cov[(i, i)] = (i as f64 + 1.0) * 0.1;
    }
    assert_eq!(cov, cov.transpose());
    for i in 0..15 {
        assert!(cov[(i, i)] > 0.0);
    }
}

#[test]
fn test_f14_state_propagation_position_step() {
    let pos0: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let vel: SpatialVelocity<Ecef<Wgs84>> = SpatialVelocity::from_vector(Vector3::new(10.0, 20.0, 0.0));
    let dt = 0.02; // 50 Hz
    let displacement = SpatialVector::from_vector(vel.vector() * dt);
    let pos1 = pos0 + displacement;
    assert_eq!(pos1.coords().x, WGS84_A + 0.2);
    assert_eq!(pos1.coords().y, 0.4);
}

// --- Feature 15: Lever Arm Rotation Enforcement ---

#[test]
fn test_f15_gnss_pos_innovation_with_rotated_lever_arm() {
    let pos_imu: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let q = r_b_e_from_rpy(0.0, 0.0, core::f64::consts::FRAC_PI_2); // 90 deg yaw
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let lever_arm = AntennaLeverArm::new(1.0, 0.0, 0.0); // 1m forward in Body
    let l_e = att.rotate_vector(lever_arm.as_body_vector());
    let predicted_ant_pos = pos_imu + l_e;
    let meas_ant_pos = pos_imu + l_e + SpatialVector::new(0.01, -0.01, 0.0);
    let innovation = meas_ant_pos - predicted_ant_pos;
    assert!((innovation.vector().x - 0.01).abs() < 1e-8);
    assert!((innovation.vector().y - (-0.01)).abs() < 1e-8);
}

#[test]
fn test_f15_attitude_error_jacobian_skew_symmetric() {
    let l_e = Vector3::new(0.0, 1.0, 0.0);
    let l_e_skew = skew_symmetric(&l_e);
    let h_att = -l_e_skew; // H_att = -[l^e x]
    // l_e x [dtheta] = - [dtheta] x l_e
    let dtheta = Vector3::new(0.01, 0.02, 0.03);
    let pert1 = h_att * dtheta;
    let pert2 = -l_e.cross(&dtheta);
    assert!((pert1 - pert2).norm() < 1e-12);
}

#[test]
fn test_f15_zero_lever_arm_attitude_jacobian_vanishes() {
    let zero_l_e = Vector3::zeros();
    let h_att = -skew_symmetric(&zero_l_e);
    assert_eq!(h_att, Matrix3::zeros());
}

#[test]
fn test_f15_doppler_velocity_lever_arm_angular_rate_coupling() {
    let q: UnitQuaternion<f64> = UnitQuaternion::identity();
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let lever_arm = AntennaLeverArm::new(1.0, 0.0, 0.0); // 1m forward
    let omega_b = Vector3::new(0.0, 0.0, 1.0); // 1 rad/s yaw rate
    // v_rot_b = omega_b x r_ant_b = [0, 0, 1] x [1, 0, 0] = [0, 1, 0] (lateral velocity)
    let v_rot_b = omega_b.cross(lever_arm.as_body_vector().vector());
    assert_eq!(v_rot_b, Vector3::new(0.0, 1.0, 0.0));
    let v_rot_e = att.rotate_velocity(&SpatialVelocity::from_vector(v_rot_b));
    assert_eq!(v_rot_e.vector(), &Vector3::new(0.0, 1.0, 0.0));
}

#[test]
fn test_f15_doppler_gyro_bias_jacobian() {
    let lever_arm_b = Vector3::new(1.0, 0.0, 0.0);
    let arm_skew = skew_symmetric(&lever_arm_b);
    let bg_err = Vector3::new(0.0, 0.0, 0.05);
    // delta_v_b = - [r_arm x] delta_bg
    let dv = -arm_skew * bg_err;
    assert_eq!(dv, Vector3::new(0.0, 0.05, 0.0));
}

// --- Feature 16: SwfgEngine & PostProcessOptions Refactoring ---

#[test]
fn test_f16_postprocess_options_typed_base_position() {
    let base_pos: Option<EcefPos<Itrf2014>> = Some(Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0)));
    assert!(base_pos.is_some());
    assert!(base_pos.as_ref().unwrap().coords().norm() > 6_000_000.0);
}

#[test]
fn test_f16_postprocess_options_typed_lever_arm() {
    let lever_arm: Option<AntennaLeverArm> = Some(AntennaLeverArm::new(0.1, -0.2, 0.5));
    assert!(lever_arm.is_some());
    assert_eq!(lever_arm.unwrap().as_body_vector().vector().x, 0.1);
}

#[test]
fn test_f16_swfg_pose_datum_consistency() {
    let base_itrf: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let rover_itrf: EcefPos<Itrf2014> = base_itrf + SpatialVector::new(10.0, 20.0, 30.0);
    let baseline = rover_itrf - base_itrf;
    assert!((baseline.norm() - 10.0_f64.hypot(20.0_f64.hypot(30.0))).abs() < 1e-12);
}

#[test]
fn test_f16_swfg_dd_factor_lever_arm_projection() {
    let q = r_b_e_from_rpy(0.0, 0.0, 0.0);
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let arm = AntennaLeverArm::new(0.5, 0.0, -1.0);
    let l_e = att.rotate_vector(arm.as_body_vector());
    assert_eq!(l_e.vector(), arm.as_body_vector().vector());
}

#[test]
fn test_f16_swfg_epoch_delta_nanosecond_precision() {
    let t_prev: Epoch<GpsScale> = Epoch::from_week_tow(2200, 100.0);
    let t_curr: Epoch<GpsScale> = Epoch::from_week_tow(2200, 100.1); // 10 Hz
    let dt = t_curr - t_prev;
    assert_eq!(dt.as_nanos(), 100_000_000);
}

// --- Feature 17: eval_odaiba_ins Refactoring ---

#[test]
fn test_f17_odaiba_base_arp_reference_position() {
    let default_base_arp = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    let base_pos: EcefPos<Wgs84> = Point3::from_coords(default_base_arp);
    assert!((base_pos.coords().norm() - 6371000.0).abs() < 5000.0);
}

#[test]
fn test_f17_odaiba_antenna_height_offset_projection() {
    let base_arp = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    let base_llh = ecef_to_llh_analytical(&base_arp);
    let r_ned = ecef_to_ned_matrix_analytical(&base_llh);
    let ned_offset = Vector3::new(0.0, 0.0, -0.0855); // 85.5 mm Up (-Down)
    let ecef_offset = r_ned.transpose() * ned_offset;
    assert!((ecef_offset.norm() - 0.0855).abs() < 1e-12);
}

#[test]
fn test_f17_odaiba_ned_covariance_to_ecef_projection() {
    let base_arp = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    let base_llh = ecef_to_llh_analytical(&base_arp);
    let r_ned = ecef_to_ned_matrix_analytical(&base_llh);
    let var_h = 0.04; // (0.2m)^2
    let var_v = 0.25; // (0.5m)^2
    let mut cov_ned = Matrix3::zeros();
    cov_ned[(0, 0)] = var_h;
    cov_ned[(1, 1)] = var_h;
    cov_ned[(2, 2)] = var_v;
    let cov_ecef = r_ned.transpose() * cov_ned * r_ned;
    assert!((cov_ecef.trace() - (2.0 * var_h + var_v)).abs() < 1e-12);
}

#[test]
fn test_f17_odaiba_zero_lever_arm_invariant() {
    let lever_arm = AntennaLeverArm::zero();
    assert_eq!(lever_arm.as_body_vector().norm(), 0.0);
}

#[test]
fn test_f17_odaiba_innovation_chi_square_consistency() {
    let residual = Vector3::new(0.02, -0.01, 0.03);
    let cov = Matrix3::identity() * 0.01;
    let chi2 = residual.transpose() * cov.try_inverse().unwrap() * residual;
    assert!(chi2[(0, 0)] > 0.0);
    assert!(chi2[(0, 0)] < 20.0); // Well within acceptance gate
}
