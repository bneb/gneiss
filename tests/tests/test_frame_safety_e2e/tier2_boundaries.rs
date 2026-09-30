//! Tier 2 E2E tests: Boundary and corner cases across spatial, temporal, geometry, and estimator domains.

use super::common::*;
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

// --- Spatial Boundary & Corner Cases ---

#[test]
fn test_t2_antimeridian_crossing_continuity() {
    let p_west = Vector3::new(-WGS84_A, -1.0, 0.0);
    let p_east = Vector3::new(-WGS84_A, 1.0, 0.0);
    let llh_west = ecef_to_llh_analytical(&p_west);
    let llh_east = ecef_to_llh_analytical(&p_east);
    // West lon is close to -pi, East lon is close to +pi
    assert!((llh_west.y.abs() - core::f64::consts::PI).abs() < 1e-4);
    assert!((llh_east.y.abs() - core::f64::consts::PI).abs() < 1e-4);
    // ECEF distance between points across antimeridian is continuous and small
    assert!((p_east - p_west).norm() < 3.0);
}

#[test]
fn test_t2_north_polar_singularity() {
    let p_north = Vector3::new(0.0, 0.0, WGS84_B);
    let llh = ecef_to_llh_analytical(&p_north);
    // Latitude exactly +90 deg (+pi/2)
    assert!((llh.x - core::f64::consts::FRAC_PI_2).abs() < 1e-8);
    // Ellipsoidal height matches semi-minor axis b
    assert!(llh.z.abs() < 1e-6);
}

#[test]
fn test_t2_south_polar_singularity() {
    let p_south = Vector3::new(0.0, 0.0, -WGS84_B);
    let llh = ecef_to_llh_analytical(&p_south);
    // Latitude exactly -90 deg (-pi/2)
    assert!((llh.x - (-core::f64::consts::FRAC_PI_2)).abs() < 1e-8);
    assert!(llh.z.abs() < 1e-6);
}

#[test]
fn test_t2_zero_lever_arm_invariant() {
    let arm = AntennaLeverArm::zero();
    let q = r_b_e_from_rpy(0.5, -0.8, 1.2);
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let l_e = att.rotate_vector(arm.as_body_vector());
    assert_eq!(l_e.norm(), 0.0);
}

#[test]
fn test_t2_extreme_10_000km_intercontinental_baseline() {
    let p1: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let p2: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(0.0, WGS84_A, 0.0));
    let baseline = p2 - p1;
    let dist = baseline.norm();
    // 90 deg chord through Earth: A * sqrt(2) approx 9,020 km
    assert!((dist - WGS84_A * 2.0_f64.sqrt()).abs() < 1.0);
}

// --- Temporal Boundary & Corner Cases ---

#[test]
fn test_t2_midnight_week_rollover_exact_nanosecond_boundary() {
    // Exact last nanosecond of week
    let t_last: Epoch<GpsScale> = Epoch::from_week_nanos(2200, WEEK_NANOS - 1);
    let dt_1ns = TimeDelta::from_nanos(1);
    let t_next = t_last + dt_1ns;
    assert_eq!(t_next.week(), 2201);
    assert_eq!(t_next.tow_nanos(), 0);
}

#[test]
fn test_t2_sub_millisecond_tolerance_exact_edge() {
    let t1: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 100_000_000_000);
    let tol = TimeDelta::from_millis(5);
    // Exactly at tolerance boundary (5ms)
    let t_edge = t1 + tol;
    assert!(t1.is_within(t_edge, tol));
    // 1 nanosecond beyond tolerance
    let t_beyond = t_edge + TimeDelta::from_nanos(1);
    assert!(!t1.is_within(t_beyond, tol));
}

#[test]
fn test_t2_leap_second_transition_instant() {
    // Instant before leap second and instant after
    let t_utc0: Epoch<UtcScale> = Epoch::from_week_tow(2200, 100.0);
    let t_gpst_pre = t_utc0.to_gpst(18);
    let t_gpst_post = t_utc0.to_gpst(19);
    let leap_step = t_gpst_post - t_gpst_pre;
    assert_eq!(leap_step.as_nanos(), 1_000_000_000);
}

#[test]
fn test_t2_negative_subtraction_order_inversion() {
    let t_early: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 100);
    let t_late: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 200);
    let dt_pos = t_late - t_early;
    let dt_neg = t_early - t_late;
    assert_eq!(dt_pos.as_nanos(), -dt_neg.as_nanos());
}

#[test]
fn test_t2_zero_duration_delta() {
    let t: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 500_000);
    let dt = t - t;
    assert_eq!(dt.as_nanos(), 0);
    let t_same = t + TimeDelta::from_nanos(0);
    assert_eq!(t_same, t);
}

// --- Geometry Boundary & Corner Cases ---

#[test]
fn test_t2_zenith_satellite_elevation_limit() {
    let p_station = Vector3::new(-3961904.0, 3348994.0, 3698211.0);
    let p_llh = ecef_to_llh_analytical(&p_station);
    let r_ned = ecef_to_ned_matrix_analytical(&p_llh);
    // Point directly up in NED: [0, 0, -20000km]
    let up_ecef = r_ned.transpose() * Vector3::new(0.0, 0.0, -20_000_000.0);
    let sat_zenith = p_station + up_ecef;
    let (_, _, el_rad, _) = los_and_az_el(&p_station, &sat_zenith);
    assert!((el_rad - core::f64::consts::FRAC_PI_2).abs() < 1e-4);
}

#[test]
fn test_t2_horizon_satellite_elevation_limit() {
    let p_station = Vector3::new(-3961904.0, 3348994.0, 3698211.0);
    let p_llh = ecef_to_llh_analytical(&p_station);
    let r_ned = ecef_to_ned_matrix_analytical(&p_llh);
    // Point directly North in NED on horizon: [20000km, 0, 0]
    let north_ecef = r_ned.transpose() * Vector3::new(20_000_000.0, 0.0, 0.0);
    let sat_horizon = p_station + north_ecef;
    let (_, _, el_rad, _) = los_and_az_el(&p_station, &sat_horizon);
    assert!(el_rad.abs() < 1e-4);
}

#[test]
fn test_t2_co_located_receivers_identity() {
    let pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let sat1: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(20000000.0, 0.0, 10000000.0));
    let sat2: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(10000000.0, 20000000.0, 10000000.0));
    let geom = DoubleDiffGeometry::compute(pos, pos, sat1, sat2);
    assert_eq!(geom.geometric_dd_m, 0.0);
    assert_eq!(geom.u_rov_sat, geom.u_bas_sat);
    assert_eq!(geom.u_rov_ref, geom.u_bas_ref);
}

#[test]
fn test_t2_collinear_satellite_difference_vector() {
    let pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    // Two satellites along the exact same line of sight at different distances
    let sat1: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(20000000.0, 0.0, 0.0));
    let sat2: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(26000000.0, 0.0, 0.0));
    let geom = DoubleDiffGeometry::compute(pos, pos, sat1, sat2);
    // Unit LOS vectors are identical, so delta_u_rov is zero
    assert!(geom.delta_u_rov.norm() < 1e-12);
}

#[test]
fn test_t2_azimuth_0_to_2pi_normalization() {
    let p_station = Vector3::new(WGS84_A, 0.0, 0.0);
    let sat = Vector3::new(WGS84_A, -10000.0, 10000.0); // South-West
    let (_, _, _, az_rad) = los_and_az_el(&p_station, &sat);
    assert!(az_rad >= 0.0);
    assert!(az_rad < 2.0 * core::f64::consts::PI);
}

// --- Estimator Boundary & Corner Cases ---

#[test]
fn test_t2_zero_error_quaternion_identity() {
    let q: UnitQuaternion<f64> = UnitQuaternion::identity();
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let v: SpatialVector<BodyFrd> = SpatialVector::new(1.0, 2.0, 3.0);
    let rotated = att.rotate_vector(&v);
    assert_eq!(rotated.vector(), v.vector());
}

#[test]
fn test_t2_pi_yaw_rotation_boundary() {
    let q = UnitQuaternion::from_euler_angles(0.0, 0.0, core::f64::consts::PI);
    let att: Attitude<BodyFrd, Ned> = Attitude::from_quaternion(q);
    let forward_b: SpatialVector<BodyFrd> = SpatialVector::new(1.0, 0.0, 0.0);
    let ned_v = att.rotate_vector(&forward_b);
    // +180 deg yaw reverses forward (+X) to South (-X in NED)
    assert!((ned_v.vector().x - (-1.0)).abs() < 1e-12);
    assert!(ned_v.vector().y.abs() < 1e-12);
}

#[test]
fn test_t2_high_angular_rate_centrifugal_coupling() {
    let lever_arm = AntennaLeverArm::new(1.0, 0.0, 0.0);
    let omega_extreme = Vector3::new(0.0, 0.0, 10.0); // 10 rad/s (aggressive spin)
    let v_rot = omega_extreme.cross(lever_arm.as_body_vector().vector());
    assert_eq!(v_rot.norm(), 10.0); // 10 m/s tangential velocity
}

#[test]
fn test_t2_near_singular_covariance_positive_eigenvalues() {
    let mut cov = Matrix3::zeros();
    cov[(0, 0)] = 1e-8;
    cov[(1, 1)] = 1e-8;
    cov[(2, 2)] = 1e-8;
    let spat_cov: SpatialCovariance<Ned> = SpatialCovariance::from_matrix(cov);
    assert!(spat_cov.std_north() > 0.0);
    assert!((spat_cov.std_north() - 1e-4).abs() < 1e-8);
}

#[test]
fn test_t2_infinitesimal_dt_position_propagation() {
    let pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(WGS84_A, 0.0, 0.0));
    let vel: SpatialVelocity<Ecef<Wgs84>> = SpatialVelocity::from_vector(Vector3::new(10.0, 0.0, 0.0));
    let dt = 1e-6; // 1 microsecond
    let delta = SpatialVector::from_vector(vel.vector() * dt);
    let pos_next = pos + delta;
    assert!((pos_next.coords().x - (WGS84_A + 1e-5)).abs() < 1e-10);
}
