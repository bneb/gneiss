//! Tier 4 E2E tests: Realistic operational mission scenarios.

use super::common::*;
use nalgebra::Vector3;

#[test]
fn test_t4_scenario1_tokyo_odaiba_kinematic_ins_with_lever_arm() {
    // Scenario 1: Tokyo Odaiba coastal highway vehicle run with lever arm [0.25, 0.10, -0.85]m
    let pos_imu: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let arm = AntennaLeverArm::new(0.25, 0.10, -0.85);
    // Vehicle pitches up 3 deg, rolls 2 deg, headings 45 deg
    let q = r_b_e_from_rpy(0.035, 0.052, core::f64::consts::FRAC_PI_4);
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let l_e = att.rotate_vector(arm.as_body_vector());
    assert!((l_e.norm() - arm.as_body_vector().norm()).abs() < 1e-12);
    let pos_ant = pos_imu + l_e;
    let baseline_ant_imu = pos_ant - pos_imu;
    assert!((baseline_ant_imu.norm() - arm.as_body_vector().norm()).abs() < 1e-8);
}

#[test]
fn test_t4_scenario2_multi_station_cors_baseline_across_datums() {
    // Scenario 2: Baseline processing between Tsukuba CORS (JGD2011) and ITRF2014 network
    let tsukuba_itrf: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(-3957408.0, 3310229.0, 3737494.0));
    let tokyo_cors: EcefPos<Itrf2014> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let baseline = tokyo_cors - tsukuba_itrf;
    let dist_m = baseline.norm();
    // Tsukuba to Tokyo is ~54 km baseline
    assert!(dist_m > 50_000.0);
    assert!(dist_m < 60_000.0);
}

#[test]
fn test_t4_scenario3_saturday_midnight_week_rollover_continuous_kinematic() {
    // Scenario 3: Continuous 50 Hz kinematic tracking across GPS week rollover
    let mut epoch = Epoch::<GpsScale>::from_week_nanos(2250, 604_799_900_000_000); // 100ms before rollover
    let dt_step = TimeDelta::from_millis(20); // 50 Hz
    let mut keys = Vec::new();
    for _ in 0..10 {
        keys.push(epoch.to_key());
        epoch = epoch + dt_step;
    }
    // Verify epoch transitioned into week 2251 cleanly
    assert_eq!(epoch.week(), 2251);
    // Verify keys are strictly monotonic across rollover
    for i in 1..keys.len() {
        assert!(keys[i] > keys[i - 1]);
    }
}

#[test]
fn test_t4_scenario4_beidou_clock_bias_free_baseline_solution() {
    // Scenario 4: Joint GPS + BDS multi-constellation RTK
    let t_gps: Epoch<GpsScale> = Epoch::from_week_tow(2200, 345600.0);
    let t_bdt = t_gps.to_bdt();
    // BDT satellite clock evaluated at emission time:
    let t_emit_bdt = t_bdt + TimeDelta::from_seconds(-0.075); // ~75ms transit
    let t_emit_gpst = t_emit_bdt.to_gpst();
    let transit = t_gps - t_emit_gpst;
    assert!((transit.as_seconds() - 0.075).abs() < 1e-12);
}

#[test]
fn test_t4_scenario5_high_dynamic_uav_aerobatics_with_doppler_updates() {
    // Scenario 5: High-dynamic UAV in 45 deg bank turn at 30 m/s with 1.5 rad/s yaw rate
    let vel_imu: SpatialVelocity<Ecef<Wgs84>> = SpatialVelocity::from_vector(Vector3::new(30.0, 0.0, 0.0));
    let arm = AntennaLeverArm::new(0.5, 0.0, -0.2);
    let omega_b = Vector3::new(0.0, 0.0, 1.5); // 1.5 rad/s yaw rate
    let v_rot_b = omega_b.cross(arm.as_body_vector().vector());
    let q = r_b_e_from_rpy(0.785, 0.1, 1.0); // 45 deg roll
    let att: Attitude<BodyFrd, Ecef<Wgs84>> = Attitude::from_quaternion(q);
    let v_rot_e = att.rotate_velocity(&SpatialVelocity::from_vector(v_rot_b));
    let v_ant = vel_imu.vector() + v_rot_e.vector();
    assert!((v_ant.norm() - 30.0).abs() < 2.0); // Antenna velocity perturbed by < 2 m/s
}
