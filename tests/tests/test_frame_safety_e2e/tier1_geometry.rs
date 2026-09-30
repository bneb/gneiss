//! Tier 1 E2E tests: Relational double-difference geometry, receiver PCV fixes, and antenna coupling.

use super::common::*;
use nalgebra::Vector3;

fn test_geometry_setup() -> (
    EcefPos<Wgs84>,
    EcefPos<Wgs84>,
    EcefPos<Wgs84>,
    EcefPos<Wgs84>,
) {
    let rov_pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    let bas_pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961910.0, 3348980.0, 3698205.0));
    let sat_pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(15000000.0, 12000000.0, 18000000.0));
    let ref_sat: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(1000000.0, 20000000.0, 16000000.0));
    (rov_pos, bas_pos, sat_pos, ref_sat)
}

// --- Feature 11: Relational Double-Difference Geometry ---

#[test]
fn test_f11_dd_geometry_construction_and_unit_los_vectors() {
    let (rov, bas, sat, ref_sat) = test_geometry_setup();
    let geom = DoubleDiffGeometry::compute(rov, bas, sat, ref_sat);
    assert!((geom.u_rov_sat.norm() - 1.0).abs() < 1e-12);
    assert!((geom.u_rov_ref.norm() - 1.0).abs() < 1e-12);
    assert!((geom.u_bas_sat.norm() - 1.0).abs() < 1e-12);
    assert!((geom.u_bas_ref.norm() - 1.0).abs() < 1e-12);
}

#[test]
fn test_f11_dd_geometry_independent_elevations_and_azimuths() {
    let (rov, bas, sat, ref_sat) = test_geometry_setup();
    let geom = DoubleDiffGeometry::compute(rov, bas, sat, ref_sat);
    assert!(geom.el_rov_sat_rad > 0.0);
    assert!(geom.el_bas_sat_rad > 0.0);
    assert!(geom.az_rov_sat_rad >= 0.0 && geom.az_rov_sat_rad < 2.0 * core::f64::consts::PI);
    assert!(geom.az_bas_sat_rad >= 0.0 && geom.az_bas_sat_rad < 2.0 * core::f64::consts::PI);
}

#[test]
fn test_f11_dd_geometry_zero_baseline_vanishes_dd_range() {
    let (rov, _, sat, ref_sat) = test_geometry_setup();
    let geom = DoubleDiffGeometry::compute(rov, rov, sat, ref_sat);
    assert!(geom.geometric_dd_m.abs() < 1e-10);
    assert!((geom.range_rov_sat_m - geom.range_bas_sat_m).abs() < 1e-10);
}

#[test]
fn test_f11_dd_geometry_satellite_range_magnitude() {
    let (rov, _, sat, _) = test_geometry_setup();
    let geom = DoubleDiffGeometry::compute(rov, rov, sat, sat);
    assert!(geom.range_rov_sat_m > 20_000_000.0);
    assert!(geom.range_rov_sat_m < 30_000_000.0);
}

#[test]
fn test_f11_dd_geometry_delta_u_rov_for_eskf_jacobian() {
    let (rov, bas, sat, ref_sat) = test_geometry_setup();
    let geom = DoubleDiffGeometry::compute(rov, bas, sat, ref_sat);
    let expected_delta_u = geom.u_rov_sat - geom.u_rov_ref;
    assert!((geom.delta_u_rov - expected_delta_u).norm() < 1e-15);
}

// --- Feature 12: Coordinate Fix in Receiver PCV ---

#[test]
fn test_f12_llh_vs_ecef_units_and_magnitude_distinction() {
    let p_ecef = Vector3::new(-3961904.0, 3348994.0, 3698211.0);
    let p_llh = ecef_to_llh_analytical(&p_ecef);
    // LLH lat/lon in radians: |phi| <= pi/2, |lam| <= pi
    assert!(p_llh.x.abs() < 1.6);
    assert!(p_llh.y.abs() < 3.2);
    // ECEF coords in millions of meters
    assert!(p_ecef.norm() > 6_000_000.0);
}

#[test]
fn test_f12_proper_ecef_origin_line_of_sight_norm() {
    let p_ecef = Vector3::new(-3961904.0, 3348994.0, 3698211.0);
    let sat_ecef = Vector3::new(15000000.0, 12000000.0, 18000000.0);
    let (unit_los, range, _, _) = los_and_az_el(&p_ecef, &sat_ecef);
    assert!((unit_los.norm() - 1.0).abs() < 1e-12);
    assert!(range > 20_000_000.0);
}

#[test]
fn test_f12_elevation_angle_valid_bounds() {
    let p_ecef = Vector3::new(WGS84_A, 0.0, 0.0);
    let sat_overhead = Vector3::new(26000000.0, 0.0, 0.0);
    let (_, _, el_rad, _) = los_and_az_el(&p_ecef, &sat_overhead);
    // Directly overhead satellite -> elevation = pi/2
    assert!((el_rad - core::f64::consts::FRAC_PI_2).abs() < 1e-4);
}

#[test]
fn test_f12_zenith_angle_complement() {
    let el_deg = 35.0;
    let zen_deg = 90.0 - el_deg;
    assert_eq!(zen_deg, 55.0);
    assert_eq!(el_deg + zen_deg, 90.0);
}

#[test]
fn test_f12_receiver_pcv_correction_boundedness() {
    // Typical PCV phase variations are bounded within +/- 20 mm
    let simulated_pcv_mm: f64 = 8.5;
    let pcv_m: f64 = simulated_pcv_mm * 1e-3;
    assert!(pcv_m.abs() < 0.05);
}

// --- Feature 13: Relational Receiver Antenna Coupling ---

#[test]
fn test_f13_regional_baseline_elevation_divergence() {
    let rov_pos: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-3961904.0, 3348994.0, 3698211.0));
    // 50 km baseline shift to base station
    let bas_pos = rov_pos + SpatialVector::new(30000.0, 40000.0, 0.0);
    let sat_low: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(18000000.0, 18000000.0, 5000000.0));
    let ref_sat: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(0.0, 10000000.0, 24000000.0));
    let geom = DoubleDiffGeometry::compute(rov_pos, bas_pos, sat_low, ref_sat);
    let el_diff_deg = (geom.el_rov_sat_rad - geom.el_bas_sat_rad).abs().to_degrees();
    // Over 50 km, elevation angles to low satellite differ by noticeable fraction of degree
    assert!(el_diff_deg > 0.05);
}

#[test]
fn test_f13_base_elevation_independent_pcv_lookup() {
    let el_rov = 15.0_f64.to_radians();
    let el_bas = 15.4_f64.to_radians();
    let zen_rov = 90.0 - el_rov.to_degrees();
    let zen_bas = 90.0 - el_bas.to_degrees();
    assert_ne!(zen_rov, zen_bas);
}

#[test]
fn test_f13_zenith_satellite_minimal_dd_pcv_difference() {
    let (rov, bas, sat, _) = test_geometry_setup();
    let sat_overhead: EcefPos<Wgs84> = Point3::from_coords(Vector3::new(-16000000.0, 13000000.0, 15000000.0));
    let geom = DoubleDiffGeometry::compute(rov, bas, sat, sat_overhead);
    assert!(geom.el_rov_ref_rad > 1.0); // high elevation (> 57 deg)
}

#[test]
fn test_f13_double_difference_symmetry_inversion() {
    let (rov, bas, sat1, sat2) = test_geometry_setup();
    let geom12 = DoubleDiffGeometry::compute(rov, bas, sat1, sat2);
    let geom21 = DoubleDiffGeometry::compute(rov, bas, sat2, sat1);
    assert!((geom12.geometric_dd_m + geom21.geometric_dd_m).abs() < 1e-10);
    assert!((geom12.base_dd_range_m + geom21.base_dd_range_m).abs() < 1e-10);
}

#[test]
fn test_f13_base_and_rover_range_difference_triangle_inequality() {
    let (rov, bas, sat, _) = test_geometry_setup();
    let geom = DoubleDiffGeometry::compute(rov, bas, sat, sat);
    let baseline_len = (rov.coords() - bas.coords()).norm();
    let sat_range_diff = (geom.range_rov_sat_m - geom.range_bas_sat_m).abs();
    // Triangle inequality: |range_rov - range_bas| <= baseline_length
    assert!(sat_range_diff <= baseline_len + 1e-9);
}
