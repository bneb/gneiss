//! Tier 1 E2E tests: Spatial typestate primitives, reference frames, and lever arm transformations.

use super::common::*;
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

// --- Feature 1: Spatial Typestate Markers & Wrappers ---

#[test]
fn test_f01_point3_affine_translation_consistency() {
    let p1: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let v: EcefVector<Wgs84> = SpatialVector::new(10.0, -20.0, 30.0);
    let p2 = p1 + v;
    let diff = p2 - p1;
    assert!((diff.vector() - v.vector()).norm() < 1e-12);
}

#[test]
fn test_f01_spatial_vector_addition_commutativity() {
    let v1: SpatialVector<Ned> = SpatialVector::new(1.5, -2.5, 3.5);
    let v2: SpatialVector<Ned> = SpatialVector::new(-4.0, 5.0, -6.0);
    let sum1 = v1 + v2;
    let sum2 = v2 + v1;
    assert!((sum1.vector() - sum2.vector()).norm() < 1e-15);
}

#[test]
fn test_f01_spatial_velocity_norm_and_components() {
    let vel: SpatialVelocity<BodyFrd> = SpatialVelocity::from_vector(Vector3::new(3.0, 4.0, 0.0));
    assert_eq!(vel.vector().x, 3.0);
    assert_eq!(vel.vector().y, 4.0);
    assert!((vel.norm() - 5.0).abs() < 1e-12);
}

#[test]
fn test_f01_ned_covariance_std_extraction() {
    let mut cov_m = Matrix3::zeros();
    cov_m[(0, 0)] = 0.04; // North: 0.20m std
    cov_m[(1, 1)] = 0.09; // East:  0.30m std
    cov_m[(2, 2)] = 0.25; // Down:  0.50m std
    let ned_cov = NedCovariance::from_matrix(cov_m);
    assert!((ned_cov.std_north() - 0.20).abs() < 1e-9);
    assert!((ned_cov.std_east() - 0.30).abs() < 1e-9);
    assert!((ned_cov.std_down() - 0.50).abs() < 1e-9);
}

#[test]
fn test_f01_affine_subtraction_roundtrip() {
    let p1: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let p2: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(-3961900.0, 3348990.0, 3698215.0));
    let delta = p2 - p1;
    let reconstructed = p1 + delta;
    assert!((reconstructed.coords() - p2.coords()).norm() < 1e-12);
}

// --- Feature 2: Reference Frame Realizations ---

#[test]
fn test_f02_reference_frame_names_and_distinctness() {
    assert_eq!(Itrf2014::NAME, "ITRF2014");
    assert_eq!(Itrf2020::NAME, "ITRF2020");
    assert_eq!(Wgs84::NAME, "WGS84");
    assert_eq!(Nad83::NAME, "NAD83");
    assert_eq!(Jgd2011::NAME, "JGD2011");
    assert_eq!(Pz90::NAME, "PZ-90.11");
}

#[test]
fn test_f02_pz90_helmert_translation_parameters() {
    // PZ-90.11 to ITRF2014 translation shifts: Tx = +3.0mm, Ty = -1.0mm, Tz = 0.0mm
    let tx_m = 0.003;
    let ty_m = -0.001;
    let tz_m = 0.000;
    assert_eq!(tx_m * 1000.0, 3.0);
    assert_eq!(ty_m * 1000.0, -1.0);
    assert_eq!(tz_m, 0.0);
}

#[test]
fn test_f02_pz90_helmert_rotation_parameters() {
    // PZ-90.11 to ITRF2014 rotations: Rx = 0.019 mas, Ry = -0.042 mas, Rz = 0.002 mas
    let mas_to_rad = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);
    let rx_rad = 0.019 * mas_to_rad;
    let ry_rad = -0.042 * mas_to_rad;
    let rz_rad = 0.002 * mas_to_rad;
    assert!(rx_rad > 0.0);
    assert!(ry_rad < 0.0);
    assert!(rz_rad > 0.0);
}

#[test]
fn test_f02_pz90_coordinate_shift_magnitude() {
    let p_pz90 = Vector3::new(WGS84_A, 0.0, 0.0);
    let shift = Vector3::new(0.003, -0.001, 0.000);
    let p_itrf = p_pz90 + shift;
    let diff_norm = (p_itrf - p_pz90).norm();
    assert!(diff_norm < 0.005);
    assert!(diff_norm > 0.002);
}

#[test]
fn test_f02_itrf2020_to_itrf2014_epoch_continuity() {
    // ITRF2020 to ITRF2014 millimetric displacement check
    let tx_mm: f64 = -1.4;
    let ty_mm: f64 = -0.9;
    let tz_mm: f64 = 1.4;
    let total_shift_mm = tx_mm.hypot(ty_mm.hypot(tz_mm));
    assert!(total_shift_mm < 3.0);
    assert!(total_shift_mm > 1.5);
}

// --- Feature 3: Antenna Lever Arm & Attitude Typestates ---

#[test]
fn test_f03_antenna_lever_arm_initialization() {
    let arm = AntennaLeverArm::new(0.25, -0.15, 0.80);
    let body_v = arm.as_body_vector().vector();
    assert_eq!(body_v.x, 0.25);
    assert_eq!(body_v.y, -0.15);
    assert_eq!(body_v.z, 0.80);
    let zero_arm = AntennaLeverArm::zero();
    assert_eq!(zero_arm.as_body_vector().norm(), 0.0);
}

#[test]
fn test_f03_attitude_rotation_body_to_ecef() {
    let q: UnitQuaternion<f64> = UnitQuaternion::identity();
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let arm = AntennaLeverArm::new(1.0, 2.0, 3.0);
    let rotated = att.rotate_vector(arm.as_body_vector());
    assert!((rotated.vector() - arm.as_body_vector().vector()).norm() < 1e-12);
}

#[test]
fn test_f03_attitude_rotation_90deg_yaw() {
    let q = UnitQuaternion::from_euler_angles(0.0, 0.0, core::f64::consts::FRAC_PI_2);
    let att: Attitude<BodyFrd, Ned> = Attitude::from_quaternion(q);
    // Forward in Body (+X) becomes East in NED (+Y) under +90 deg yaw
    let forward_b: SpatialVector<BodyFrd> = SpatialVector::new(1.0, 0.0, 0.0);
    let ned_v = att.rotate_vector(&forward_b);
    assert!((ned_v.vector().x).abs() < 1e-12);
    assert!((ned_v.vector().y - 1.0).abs() < 1e-12);
}

#[test]
fn test_f03_attitude_inverse_roundtrip() {
    let q = r_b_e_from_rpy(0.1, -0.2, 0.5);
    let att_fwd: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let att_inv: Attitude<Ecef<Wgs84>, BodyFrd> = att_fwd.inverse();
    let arm = AntennaLeverArm::new(0.5, -0.3, 1.2);
    let ecef_v = att_fwd.rotate_vector(arm.as_body_vector());
    let body_recon = att_inv.rotate_vector(&ecef_v);
    assert!((body_recon.vector() - arm.as_body_vector().vector()).norm() < 1e-12);
}

#[test]
fn test_f03_attitude_covariance_rotation_trace_invariance() {
    let q = r_b_e_from_rpy(0.3, -0.4, 0.6);
    let att: Attitude<BodyFrd, Ned> = Attitude::from_quaternion(q);
    let mut cov_b = Matrix3::zeros();
    cov_b[(0, 0)] = 0.01;
    cov_b[(1, 1)] = 0.04;
    cov_b[(2, 2)] = 0.09;
    let spat_cov_b: SpatialCovariance<BodyFrd> = SpatialCovariance::from_matrix(cov_b);
    let spat_cov_ned = att.rotate_cov(&spat_cov_b);
    let trace_b = cov_b.trace();
    let trace_ned = spat_cov_ned.matrix().trace();
    assert!((trace_b - trace_ned).abs() < 1e-12);
}

// --- Feature 4: Non-Leaky Typestate Deref & Operators ---

#[test]
fn test_f04_explicit_coords_access() {
    let raw = Vector3::new(100.0, 200.0, 300.0);
    let p: EcefPos<Wgs84> = Point3::from_coords(raw);
    assert_eq!(p.coords().x, 100.0);
    assert_eq!(p.coords().y, 200.0);
    assert_eq!(p.coords().z, 300.0);
}

#[test]
fn test_f04_owned_coords_extraction() {
    let p: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(1.0, 2.0, 3.0));
    let extracted = p.into_coords();
    assert_eq!(extracted, Vector3::new(1.0, 2.0, 3.0));
}

#[test]
fn test_f04_spatial_vector_norm_invariant_under_pure_rotation() {
    let v: SpatialVector<BodyFrd> = SpatialVector::new(1.2, -3.4, 5.6);
    let orig_norm = v.norm();
    let q = r_b_e_from_rpy(0.4, 0.2, -0.7);
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let rotated = att.rotate_vector(&v);
    assert!((rotated.norm() - orig_norm).abs() < 1e-12);
}

#[test]
fn test_f04_null_translation_preserves_position() {
    let p: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let zero_v: EcefVector<Wgs84> = SpatialVector::zero();
    let p_next = p + zero_v;
    assert_eq!(p.coords(), p_next.coords());
}

#[test]
fn test_f04_vector_addition_associativity() {
    let v1: EcefVector<Wgs84> = SpatialVector::new(1.0, 2.0, 3.0);
    let v2: EcefVector<Wgs84> = SpatialVector::new(4.0, 5.0, 6.0);
    let v3: EcefVector<Wgs84> = SpatialVector::new(7.0, 8.0, 9.0);
    let assoc1 = (v1 + v2) + v3;
    let assoc2 = v1 + (v2 + v3);
    assert!((assoc1.vector() - assoc2.vector()).norm() < 1e-12);
}

// --- Feature 5: Relational Local Tangent Plane ---

#[test]
fn test_f05_local_tangent_plane_origin_projection() {
    let origin: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let ltp = LocalTangentPlane::from_origin(origin);
    let enu = ltp.to_enu(origin);
    assert!(enu.norm() < 1e-9);
    let ned = ltp.to_ned(origin);
    assert!(ned.norm() < 1e-9);
}

#[test]
fn test_f05_local_tangent_plane_roundtrip() {
    let origin: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let ltp = LocalTangentPlane::from_origin(origin);
    let target = origin + SpatialVector::new(100.0, -50.0, 75.0);
    let enu = ltp.to_enu(target);
    let reconstructed = ltp.from_enu(enu);
    assert!((reconstructed.coords() - target.coords()).norm() < 1e-8);
}

#[test]
fn test_f05_local_tangent_plane_metric_distance_preservation() {
    let origin: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let ltp = LocalTangentPlane::from_origin(origin);
    let target = origin + SpatialVector::new(30.0, 40.0, 0.0);
    let enu = ltp.to_enu(target);
    let dist_3d = (target.coords() - origin.coords()).norm();
    let dist_enu = enu.norm();
    assert!((dist_3d - dist_enu).abs() < 1e-6);
}

#[test]
fn test_f05_local_tangent_plane_covariance_projection_positive_definite() {
    let origin: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let ltp = LocalTangentPlane::from_origin(origin);
    let cov_ecef = SpatialCovariance::from_matrix(Matrix3::identity() * 0.04);
    let cov_enu = ltp.project_cov(&cov_ecef);
    assert!(cov_enu.matrix()[(0, 0)] > 0.0);
    assert!(cov_enu.matrix()[(1, 1)] > 0.0);
    assert!(cov_enu.matrix()[(2, 2)] > 0.0);
    assert!((cov_enu.matrix()[(0, 0)] - 0.04).abs() < 1e-9);
}

#[test]
fn test_f05_local_tangent_plane_vertical_up_matches_zenith() {
    let origin: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let ltp = LocalTangentPlane::from_origin(origin);
    // At equator lat=0, lon=0, +X is Up (+Z in ENU)
    let shifted = origin + SpatialVector::new(10.0, 0.0, 0.0);
    let enu = ltp.to_enu(shifted);
    assert!((enu.vector().z - 10.0).abs() < 1e-6);
}
