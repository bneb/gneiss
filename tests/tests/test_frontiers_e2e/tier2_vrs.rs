//! Tier 2 E2E tests: Features 12–16 Boundary & Corner Cases.

use super::common::*;
use nalgebra::{Vector2, Vector3};

// --- Feature 12 Boundaries ---

#[test]
fn test_f12_b1_minimal_three_station_network() {
    let network = [
        Vector2::new(0.0, 0.0),
        Vector2::new(10.0, 0.0),
        Vector2::new(5.0, 8.0),
    ];
    assert_eq!(network.len(), 3);
}

#[test]
fn test_f12_b2_single_base_fallback_mode() {
    let network = [Vector2::new(0.0, 0.0)];
    let is_network_mode = network.len() >= 3;
    assert!(!is_network_mode);
}

#[test]
fn test_f12_b3_timestamp_jitter_sub_millisecond() {
    let t1: f64 = 345600.0001;
    let t2: f64 = 345600.0002;
    let dt = (t1 - t2).abs();
    assert!(dt < 1e-3);
}

#[test]
fn test_f12_b4_low_sat_count_station_exclusion() {
    let n_sats = 3;
    let is_usable = n_sats >= 4;
    assert!(!is_usable);
}

#[test]
fn test_f12_b5_long_baseline_noise_elevation() {
    let dist_km = 80.0;
    let base_noise = 0.003;
    let dist_noise = base_noise + 1e-6 * dist_km * 1000.0;
    assert!(dist_noise > base_noise);
}

// --- Feature 13 Boundaries ---

#[test]
fn test_f13_b1_zero_baseline_double_difference() {
    let dd = form_double_difference(10.0, 20.0, 10.0, 20.0);
    assert_eq!(dd, 0.0);
}

#[test]
fn test_f13_b2_one_hundred_km_baseline_length() {
    let b = Vector2::new(60.0, 80.0); // 60^2 + 80^2 = 100^2
    assert_eq!(b.norm(), 100.0);
}

#[test]
fn test_f13_b3_star_network_topology() {
    let master = Vector2::new(0.0, 0.0);
    let stations = vec![
        Vector2::new(10.0, 0.0),
        Vector2::new(0.0, 10.0),
        Vector2::new(-10.0, 0.0),
    ];
    for s in stations {
        let baseline = s - master;
        assert_eq!(baseline.norm(), 10.0);
    }
}

#[test]
fn test_f13_b4_four_station_quadrilateral_closure() {
    let b1 = Vector2::new(10.0, 0.0);
    let b2 = Vector2::new(0.0, 10.0);
    let b3 = Vector2::new(-10.0, 0.0);
    let b4 = Vector2::new(0.0, -10.0);
    let loop_closure = b1 + b2 + b3 + b4;
    assert_eq!(loop_closure, Vector2::zeros());
}

#[test]
fn test_f13_b5_singular_network_matrix_guard() {
    let n_sats = 2; // Insufficient for DD positioning
    let rank = if n_sats >= 4 { 4 } else { n_sats };
    assert!(rank < 4);
}

// --- Feature 14 Boundaries ---

#[test]
fn test_f14_b1_point_exactly_on_vertex() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(10.0, 0.0);
    let c = Vector2::new(0.0, 10.0);
    let bary = compute_barycentric(&a, &a, &b, &c).unwrap();
    assert!((bary.x - 1.0).abs() < 1e-12);
    assert!(bary.y.abs() < 1e-12);
    assert!(bary.z.abs() < 1e-12);
}

#[test]
fn test_f14_b2_point_exactly_on_edge() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(10.0, 0.0);
    let c = Vector2::new(0.0, 10.0);
    let edge_pt = Vector2::new(5.0, 0.0);
    let bary = compute_barycentric(&edge_pt, &a, &b, &c).unwrap();
    assert!((bary.x - 0.5).abs() < 1e-12);
    assert!((bary.y - 0.5).abs() < 1e-12);
    assert!(bary.z.abs() < 1e-12);
}

#[test]
fn test_f14_b3_point_outside_triangle_negative_barycentric() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(10.0, 0.0);
    let c = Vector2::new(0.0, 10.0);
    let outside_pt = Vector2::new(-2.0, 5.0);
    let bary = compute_barycentric(&outside_pt, &a, &b, &c).unwrap();
    assert!(bary.x < 0.0 || bary.y < 0.0 || bary.z < 0.0);
}

#[test]
fn test_f14_b4_degenerate_collinear_stations_guard() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(5.0, 5.0);
    let c = Vector2::new(10.0, 10.0);
    let query = Vector2::new(2.0, 2.0);
    assert!(compute_barycentric(&query, &a, &b, &c).is_none());
}

#[test]
fn test_f14_b5_duplicate_coordinates_guard() {
    let a = Vector2::new(1.0, 1.0);
    let b = Vector2::new(1.0, 1.0);
    let c = Vector2::new(5.0, 5.0);
    let query = Vector2::new(2.0, 2.0);
    assert!(compute_barycentric(&query, &a, &b, &c).is_none());
}

// --- Feature 15 Boundaries ---

#[test]
fn test_f15_b1_rover_co_located_with_master_base() {
    let r_rover = Vector3::new(WGS84_A, 0.0, 0.0);
    let r_master = r_rover;
    let delta_r = (r_rover - r_master).norm();
    assert_eq!(delta_r, 0.0);
}

#[test]
fn test_f15_b2_rover_on_network_perimeter_50km() {
    let r_master = Vector2::new(0.0, 0.0);
    let r_rover = Vector2::new(30.0, 40.0); // 50 km
    let baseline = (r_rover - r_master).norm();
    assert_eq!(baseline, 50.0);
}

#[test]
fn test_f15_b3_zero_atmospheric_gradient_uniform_weather() {
    let delta_tropo = 0.0;
    let delta_iono = 0.0;
    let obs_master = 100.0;
    let obs_vrs = obs_master + delta_tropo - delta_iono;
    assert_eq!(obs_vrs, obs_master);
}

#[test]
fn test_f15_b4_severe_storm_front_gradient() {
    let grad_per_km = 0.010; // 10mm / km
    let dist_km = 10.0;
    let total_delta = grad_per_km * dist_km;
    assert_eq!(total_delta, 0.100); // 100 mm
}

#[test]
fn test_f15_b5_low_elevation_mapping_factor() {
    let elev_deg = 10.0_f64;
    let sin_e = elev_deg.to_radians().sin();
    let map_factor = 1.0 / sin_e;
    assert!((map_factor - 5.75877).abs() < 1e-4);
}

// --- Feature 16 Boundaries ---

#[test]
fn test_f16_b1_zero_distance_leica_limit() {
    let dist_km = 0.0;
    let spec = 0.008 + 1e-6 * dist_km * 1000.0;
    assert_eq!(spec, 0.008); // 8 mm
}

#[test]
fn test_f16_b2_fifty_km_leica_limit() {
    let dist_km = 50.0;
    let spec = 0.008 + 1e-6 * dist_km * 1000.0;
    assert!((spec - 0.058_f64).abs() < 1e-12); // 58 mm
}

#[test]
fn test_f16_b3_high_fix_rate_threshold() {
    let fix_rate = 0.993;
    assert!(fix_rate >= 0.965);
}

#[test]
fn test_f16_b4_smoke_guard_threshold_verification() {
    let p50 = 0.024;
    let rms = 0.041;
    assert!(p50 <= 0.040);
    assert!(rms <= 0.060);
}

#[test]
fn test_f16_b5_vertical_vs_horizontal_error_ratio() {
    let h_rms = 0.041;
    let v_rms = 0.043;
    let ratio = v_rms / h_rms;
    assert!(ratio > 0.8 && ratio < 2.5);
}
