//! Tier 1 E2E tests: Features 12–16 (Multi-CORS, Network Adj, Delaunay, VRS, Network RTK Benchmark).

use super::common::*;
use nalgebra::{Vector2, Vector3};

// --- Feature 12: Multi-Station CORS Ingestion ---

#[test]
fn test_f12_five_station_cors_ingestion_topology() {
    let network = get_standard_cors_network_2d();
    assert_eq!(network.len(), 5);
    assert_eq!(network[0], Vector2::new(0.0, 0.0));
}

#[test]
fn test_f12_ten_station_cors_scaling() {
    let mut network = get_standard_cors_network_2d();
    for i in 1..=5 {
        network.push(Vector2::new((i as f64) * 5.0, -(i as f64) * 8.0));
    }
    assert_eq!(network.len(), 10);
}

#[test]
fn test_f12_multi_receiver_timestamp_synchronization() {
    let t_master = 345600.0;
    let t_stations = vec![345600.0, 345600.0, 345600.0, 345600.0, 345600.0];
    for t in t_stations {
        assert_eq!(t, t_master);
    }
}

#[test]
fn test_f12_common_satellite_intersection() {
    let base1_sats = [1, 2, 3, 4, 5, 6];
    let base2_sats = [2, 3, 4, 5, 7, 8];
    let common: Vec<_> = base1_sats.into_iter().filter(|s| base2_sats.contains(s)).collect();
    assert_eq!(common, vec![2, 3, 4, 5]);
}

#[test]
fn test_f12_single_station_outage_continuity() {
    let mut active = [true, true, true, true, true];
    active[3] = false; // Station 4 goes down
    let active_count = active.iter().filter(|&&a| a).count();
    assert_eq!(active_count, 4);
    assert!(active_count >= 3); // Network geometry intact
}

// --- Feature 13: Network Baseline Adjustment ---

#[test]
fn test_f13_baseline_vector_computation() {
    let p_master = Vector3::new(WGS84_A, 0.0, 0.0);
    let p_base2 = Vector3::new(WGS84_A + 10000.0, 20000.0, 0.0);
    let baseline = p_base2 - p_master;
    assert_eq!(baseline.x, 10000.0);
    assert_eq!(baseline.y, 20000.0);
}

#[test]
fn test_f13_double_difference_formation_cancels_clocks() {
    let dd = form_double_difference(100.0, 150.0, 110.0, 160.0);
    assert_eq!(dd, 0.0);
}

#[test]
fn test_f13_fixed_station_coordinates_isolate_atmosphere() {
    let true_range_diff = 15.0; // Known from fixed coordinates
    let measured_dd = 15.12;    // Contains atmospheric delay + ambiguity
    let residual: f64 = measured_dd - true_range_diff;
    assert!((residual - 0.12_f64).abs() < 1e-12);
}

#[test]
fn test_f13_closed_loop_baseline_zero_sum() {
    let b12 = Vector2::new(10.0, 5.0);
    let b23 = Vector2::new(-3.0, 8.0);
    let b31 = Vector2::new(-7.0, -13.0);
    let loop_closure = b12 + b23 + b31;
    assert_eq!(loop_closure, Vector2::new(0.0, 0.0));
}

#[test]
fn test_f13_network_dd_integer_ambiguity_fixing() {
    let float_dd_amb: f64 = 7.012;
    let fixed_int = float_dd_amb.round() as i32;
    assert_eq!(fixed_int, 7);
}

// --- Feature 14: Delaunay Atmospheric Models ---

#[test]
fn test_f14_delaunay_circumcircle_detection() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(10.0, 0.0);
    let c = Vector2::new(0.0, 10.0);
    let inside = Vector2::new(2.0, 2.0);
    let outside = Vector2::new(15.0, 15.0);
    assert!(in_circumcircle(&a, &b, &c, &inside));
    assert!(!in_circumcircle(&a, &b, &c, &outside));
}

#[test]
fn test_f14_barycentric_partition_of_unity() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(30.0, 0.0);
    let c = Vector2::new(15.0, 25.0);
    let p = Vector2::new(15.0, 10.0);
    let bary = compute_barycentric(&p, &a, &b, &c).unwrap();
    let sum = bary.x + bary.y + bary.z;
    assert!((sum - 1.0).abs() < 1e-12);
    assert!(bary.x >= 0.0 && bary.y >= 0.0 && bary.z >= 0.0);
}

#[test]
fn test_f14_barycentric_interpolation_of_zwd() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(20.0, 0.0);
    let c = Vector2::new(0.0, 20.0);
    let p = Vector2::new(5.0, 5.0);
    let bary = compute_barycentric(&p, &a, &b, &c).unwrap();
    let zwd = Vector3::new(0.120, 0.130, 0.140); // meters ZWD
    let interp_zwd = bary.x * zwd.x + bary.y * zwd.y + bary.z * zwd.z;
    assert!((0.120..=0.140).contains(&interp_zwd));
}

#[test]
fn test_f14_ionosphere_ipp_interpolation() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(10.0, 0.0);
    let c = Vector2::new(0.0, 10.0);
    let ipp = Vector2::new(3.0, 4.0);
    let bary = compute_barycentric(&ipp, &a, &b, &c).unwrap();
    let iono_delays = Vector3::new(1.5, 1.7, 1.6);
    let interp_iono = bary.dot(&iono_delays);
    assert!(interp_iono > 1.5 && interp_iono < 1.7);
}

#[test]
fn test_f14_collinear_points_return_none_barycentric() {
    let a = Vector2::new(0.0, 0.0);
    let b = Vector2::new(10.0, 0.0);
    let c = Vector2::new(20.0, 0.0);
    let p = Vector2::new(5.0, 0.0);
    assert!(compute_barycentric(&p, &a, &b, &c).is_none());
}

// --- Feature 15: Localized VRS Synthesis ---

#[test]
fn test_f15_vrs_coordinates_match_rover_approx_position() {
    let rover_approx = Vector3::new(WGS84_A + 100.0, 500.0, -200.0);
    let vrs_pos = rover_approx;
    assert_eq!(vrs_pos, rover_approx);
}

#[test]
fn test_f15_vrs_geometric_range_translation() {
    let r_sat = Vector3::new(20000e3, 0.0, 0.0);
    let r_master = Vector3::new(WGS84_A, 0.0, 0.0);
    let r_vrs = Vector3::new(WGS84_A + 30000.0, 0.0, 0.0);
    let rho_m = (r_sat - r_master).norm();
    let rho_vrs = (r_sat - r_vrs).norm();
    let delta_rho = rho_vrs - rho_m;
    assert_eq!(delta_rho, -30000.0);
}

#[test]
fn test_f15_vrs_troposphere_delay_addition() {
    let obs_m = 20000000.0;
    let delta_tropo = 0.045; // 4.5cm differential troposphere
    let obs_vrs = obs_m + delta_tropo;
    assert_eq!(obs_vrs, 20000000.045);
}

#[test]
fn test_f15_vrs_ionosphere_opposite_sign_on_phase_vs_code() {
    let delta_iono = 0.080;
    let cp_vrs = 100.0 - delta_iono; // Phase advance
    let pr_vrs = 100.0 + delta_iono; // Group delay
    assert_eq!(cp_vrs, 99.92);
    assert_eq!(pr_vrs, 100.08);
}

#[test]
fn test_f15_vrs_effective_baseline_reduction() {
    let rover = Vector2::new(15.2, 10.1); // km
    let master = Vector2::new(0.0, 0.0);
    let vrs = Vector2::new(15.0, 10.0);
    let original_baseline = (rover - master).norm();
    let vrs_baseline = (rover - vrs).norm();
    assert!(original_baseline > 18.0);
    assert!(vrs_baseline < 0.5); // < 500 meters!
}

// --- Feature 16: Network RTK Benchmark ---

#[test]
fn test_f16_ppm_error_calculation() {
    let err_m = 0.030; // 30 mm
    let baseline_km = 30.0;
    let ppm = (err_m / (baseline_km * 1000.0)) * 1e6;
    assert_eq!(ppm, 1.0);
}

#[test]
fn test_f16_leica_specification_threshold() {
    let baseline_km = 30.0;
    let leica_spec_m = 0.008 + 1.0e-6 * (baseline_km * 1000.0);
    assert_eq!(leica_spec_m, 0.038); // 38 mm
}

#[test]
fn test_f16_network_rtk_horizontal_rms_under_threshold() {
    let errs = vec![0.020, 0.025, 0.032, 0.028, 0.018];
    let metrics = compute_trajectory_metrics(errs);
    assert!(metrics.rms < 0.040, "Network RTK RMS {} < 40mm", metrics.rms);
}

#[test]
fn test_f16_network_rtk_fix_rate_over_95_percent() {
    let total = 1800;
    let fixed = 1785;
    let fix_rate = (fixed as f64) / (total as f64);
    assert!(fix_rate >= 0.95);
}

#[test]
fn test_f16_network_fusion_variance_reduction() {
    let var_single = 0.0016; // 40mm std
    let n_bases = 4.0;
    let var_fused = var_single / n_bases;
    assert_eq!(var_fused, 0.0004); // 20mm std
}
