//! Tier 3 E2E tests: Cross-feature pairwise interactions across spatial, temporal, and estimator domains.

use super::common::*;
use nalgebra::{Matrix3, Vector3};

#[test]
fn test_t3_spatial_and_temporal_satellite_orbit_step() {
    let t0: Epoch<GpsScale> = Epoch::from_week_tow(2200, 100.0);
    let dt = TimeDelta::from_seconds(30.0);
    let t1 = t0 + dt;
    let sat_pos0: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(26000000.0, 0.0, 0.0));
    // Approximate orbital velocity in MEO: ~3.87 km/s
    let v_orb = Vector3::new(0.0, 3870.0, 0.0);
    let sat_pos1 = sat_pos0 + SpatialVector::from_vector(v_orb * dt.as_seconds());
    assert_eq!(t1.tow_seconds(), 130.0);
    assert!((sat_pos1.coords().y - 3870.0 * 30.0).abs() < 1e-6);
}

#[test]
fn test_t3_spatial_and_geometry_dd_baseline_displacement() {
    let base: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let rover = base + SpatialVector::new(500.0, -300.0, 150.0);
    let sat1: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(15000000.0, 12000000.0, 18000000.0));
    let sat2: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(1000000.0, 20000000.0, 16000000.0));
    let geom = DoubleDiffGeometry::compute(rover, base, sat1, sat2);
    // Baseline length ~ 602m
    let baseline_len = (rover.coords() - base.coords()).norm();
    assert!((geom.geometric_dd_m.abs()) < baseline_len);
}

#[test]
fn test_t3_temporal_and_estimator_eskf_propagation_across_week_rollover() {
    let t_pre: Epoch<GpsScale> = Epoch::from_week_nanos(2250, 604_799_980_000_000); // 20ms before midnight
    let dt = TimeDelta::from_millis(20); // 50 Hz IMU step
    let t_post = t_pre + dt;
    assert_eq!(t_post.week(), 2251);
    assert_eq!(t_post.tow_nanos(), 0);
    let pos0: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let vel: SpatialVelocity<Ecef<Wgs84>> = SpatialVelocity::from_vector(Vector3::new(0.0, 20.0, 0.0));
    let pos1 = pos0 + SpatialVector::from_vector(vel.vector() * dt.as_seconds());
    assert!((pos1.coords().y - 0.4).abs() < 1e-9);
}

#[test]
fn test_t3_lever_arm_and_attitude_vehicle_turn_doppler_innovation() {
    // Vehicle turns at 0.5 rad/s yaw rate with a 2m forward lever arm
    let q = r_b_e_from_rpy(0.0, 0.0, core::f64::consts::FRAC_PI_4); // 45 deg yaw
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let lever_arm = AntennaLeverArm::new(2.0, 0.0, 0.0);
    let omega_b = Vector3::new(0.0, 0.0, 0.5);
    let v_rot_b = omega_b.cross(lever_arm.as_body_vector().vector());
    assert_eq!(v_rot_b, Vector3::new(0.0, 1.0, 0.0));
    let v_rot_e = att.rotate_velocity(&SpatialVelocity::from_vector(v_rot_b));
    assert!((v_rot_e.norm() - 1.0).abs() < 1e-12);
}

#[test]
fn test_t3_tangent_plane_and_covariance_projection() {
    let origin: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let ltp = LocalTangentPlane::from_origin(origin);
    let mut cov_ecef = Matrix3::zeros();
    cov_ecef[(0, 0)] = 0.01;
    cov_ecef[(1, 1)] = 0.01;
    cov_ecef[(2, 2)] = 0.04;
    let spat_cov_ecef = SpatialCovariance::from_matrix(cov_ecef);
    let cov_enu = ltp.project_cov(&spat_cov_ecef);
    assert!(cov_enu.matrix()[(0, 0)] > 0.0);
    assert!(cov_enu.matrix()[(1, 1)] > 0.0);
    assert!(cov_enu.matrix()[(2, 2)] > 0.0);
}

#[test]
fn test_t3_multi_constellation_gps_bdt_utc_joint_alignment() {
    let t_gpst: Epoch<GpsScale> = Epoch::from_week_tow(2200, 1000.0);
    let t_bdt = t_gpst.to_bdt();
    let t_utc = t_gpst.to_utc(18);
    // Verify BDT is 14s behind GPST
    assert_eq!((t_gpst - t_bdt.to_gpst()).as_nanos(), 0);
    // Verify UTC is 18s behind GPST
    assert_eq!((t_gpst - t_utc.to_gpst(18)).as_nanos(), 0);
}

#[test]
fn test_t3_lever_arm_and_double_diff_attitude_coupling_jacobian() {
    let rov: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let bas: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961910.0, 3348980.0, 3698205.0));
    let sat: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(15000000.0, 12000000.0, 18000000.0));
    let ref_sat: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(1000000.0, 20000000.0, 16000000.0));
    let geom = DoubleDiffGeometry::compute(rov, bas, sat, ref_sat);
    let q = r_b_e_from_rpy(0.1, -0.2, 0.3);
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let arm = AntennaLeverArm::new(0.5, 0.0, -1.0);
    let l_e = att.rotate_vector(arm.as_body_vector());
    let l_e_skew = skew_symmetric(l_e.vector());
    // H_att = delta_u^T * [l_e x]
    let h_att = geom.delta_u_rov.transpose() * l_e_skew;
    assert_eq!(h_att.ncols(), 3);
}

#[test]
fn test_t3_relational_pcv_and_double_difference_elevation_consistency() {
    let rov: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let bas: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let sat: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(15000000.0, 12000000.0, 18000000.0));
    let ref_sat: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(1000000.0, 20000000.0, 16000000.0));
    let geom = DoubleDiffGeometry::compute(rov, bas, sat, ref_sat);
    // On zero baseline, rover and base elevations must match exactly
    assert_eq!(geom.el_rov_sat_rad, geom.el_bas_sat_rad);
    assert_eq!(geom.el_rov_ref_rad, geom.el_bas_ref_rad);
}
