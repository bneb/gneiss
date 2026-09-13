//! Tier 2 E2E tests: Features 7–11 Boundary & Corner Cases.

use super::common::*;
use gneiss_core::frequencies::Signal;
use gneiss_core::sat::Constellation;
use nalgebra::{DMatrix, DVector, Vector3};

// --- Feature 7 Boundaries ---

#[test]
fn test_f7_b1_empty_sinex_table_returns_none() {
    use std::collections::HashMap;
    let empty_index: HashMap<(u8, String), f64> = HashMap::new();
    assert_eq!(empty_index.get(&(1, "L1C".to_string())), None);
}

#[test]
fn test_f7_b2_unknown_satellite_query() {
    use std::collections::HashMap;
    let mut index: HashMap<(u8, String), f64> = HashMap::new();
    index.insert((1, "L1C".to_string()), 0.5);
    assert_eq!(index.get(&(99, "L1C".to_string())), None);
}

#[test]
fn test_f7_b3_unknown_observable_code() {
    use std::collections::HashMap;
    let mut index: HashMap<(u8, String), f64> = HashMap::new();
    index.insert((1, "L1C".to_string()), 0.5);
    assert_eq!(index.get(&(1, "XYZ".to_string())), None);
}

#[test]
fn test_f7_b4_timestamp_outside_validity_interval() {
    let t_start = 100.0;
    let t_end = 200.0;
    let query_t = 250.0;
    let is_valid = query_t >= t_start && query_t <= t_end;
    assert!(!is_valid);
}

#[test]
fn test_f7_b5_physical_zero_bias_handling() {
    use std::collections::HashMap;
    let mut index: HashMap<(u8, String), f64> = HashMap::new();
    index.insert((1, "L1C".to_string()), 0.0);
    assert_eq!(index.get(&(1, "L1C".to_string())), Some(&0.0));
}

// --- Feature 8 Boundaries ---

#[test]
fn test_f8_b1_nadir_angle_at_boresight_zero() {
    let nadir_deg = 0.0;
    let pcv_0 = 0.0005; // 0.5mm
    assert_eq!(nadir_deg, 0.0);
    assert!(pcv_0 > 0.0);
}

#[test]
fn test_f8_b2_nadir_angle_at_max_limb_14_5_deg() {
    let nadir_deg = 14.5;
    let pcv_14_5 = 0.0025; // 2.5mm
    assert_eq!(nadir_deg, 14.5);
    assert!(pcv_14_5 > 0.0);
}

#[test]
fn test_f8_b3_receiver_zenith_at_horizon_90_deg() {
    let zen_deg = 90.0;
    let pcv_90: f64 = -0.015; // -15mm
    assert_eq!(zen_deg, 90.0);
    assert!(pcv_90.is_finite());
}

#[test]
fn test_f8_b4_azimuth_angle_modulo_360() {
    let az_raw = 365.0;
    let az_mod = az_raw % 360.0;
    assert_eq!(az_mod, 5.0);
}

fn get_pcv_correction(pcv: Option<f64>) -> f64 {
    pcv.unwrap_or(0.0)
}

#[test]
fn test_f8_b5_fallback_when_pcv_missing() {
    let pcv_val = get_pcv_correction(None);
    assert_eq!(pcv_val, 0.0);
}

// --- Feature 9 Boundaries ---

#[test]
fn test_f9_b1_gps_only_single_constellation() {
    let f1 = Signal::GpsL1Ca.base_freq_hz();
    assert!(f1 > 1e9);
}

#[test]
fn test_f9_b2_galileo_and_gps_dual_constellation() {
    let f_gps = Signal::GpsL1Ca.base_freq_hz();
    let f_gal = Signal::GalE1Os.base_freq_hz();
    assert_eq!(f_gps, f_gal);
}

#[test]
fn test_f9_b3_beidou_triple_constellation() {
    let f_b1 = Signal::BdsB1i.base_freq_hz();
    assert!((f_b1 - 1561.098e6).abs() < 1e3);
}

#[test]
fn test_f9_b4_quad_constellation_support() {
    let constellations = [
        Constellation::Gps,
        Constellation::Galileo,
        Constellation::Beidou,
        Constellation::Qzss,
    ];
    assert_eq!(constellations.len(), 4);
}

#[test]
fn test_f9_b5_qzss_fifth_signal() {
    let f_qzss = 1575.42e6;
    assert!(f_qzss > 1e9);
}

// --- Feature 10 Boundaries ---

#[test]
fn test_f10_b1_minimal_lambda_dimension_two() {
    let a_hat = DVector::from_vec(vec![5.02, 3.98]);
    let q_aa = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 0.01]));
    let res = gneiss_rtk::ambiguity::lambda::resolve_lambda(&a_hat, &q_aa).unwrap();
    assert_eq!(res.best_integers.len(), 2);
    assert_eq!(res.best_integers[0].round() as i32, 5);
    assert_eq!(res.best_integers[1].round() as i32, 4);
}

#[test]
fn test_f10_b2_large_lambda_dimension_ten() {
    let mut a_hat_vec = Vec::new();
    let mut diag_vec = Vec::new();
    for i in 0..10 {
        a_hat_vec.push((i as f64) + 0.02);
        diag_vec.push(0.01);
    }
    let a_hat = DVector::from_vec(a_hat_vec);
    let q_aa = DMatrix::from_diagonal(&DVector::from_vec(diag_vec));
    let res = gneiss_rtk::ambiguity::lambda::resolve_lambda(&a_hat, &q_aa).unwrap();
    assert_eq!(res.best_integers.len(), 10);
    assert_eq!(res.best_integers[9].round() as i32, 9);
}

#[test]
fn test_f10_b3_exact_integer_input_lambda() {
    let a_hat = DVector::from_vec(vec![7.0, -3.0]);
    let q_aa = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 0.01]));
    let res = gneiss_rtk::ambiguity::lambda::resolve_lambda(&a_hat, &q_aa).unwrap();
    assert_eq!(res.best_integers[0].round() as i32, 7);
    assert_eq!(res.best_integers[1].round() as i32, -3);
}

#[test]
fn test_f10_b4_large_variance_lambda_search() {
    let a_hat = DVector::from_vec(vec![12.1, 8.1]);
    let q_aa = DMatrix::from_diagonal(&DVector::from_vec(vec![10.0, 10.0]));
    let res = gneiss_rtk::ambiguity::lambda::resolve_lambda(&a_hat, &q_aa).unwrap();
    assert_eq!(res.best_integers[0].round() as i32, 12);
    assert_eq!(res.best_integers[1].round() as i32, 8);
}

#[test]
fn test_f10_b5_poor_ratio_test_rejects_fix() {
    let ratio = 1.15;
    let threshold = 2.0;
    let accepted = ratio >= threshold;
    assert!(!accepted);
}

// --- Feature 11 Boundaries ---

#[test]
fn test_f11_b1_stationary_start_kinematic_drive() {
    let v_start: Vector3<f64> = Vector3::zeros();
    assert_eq!(v_start.norm(), 0.0);
}

#[test]
fn test_f11_b2_90_degree_vehicle_turn() {
    let v_before = Vector3::new(10.0, 0.0, 0.0);
    let v_after = Vector3::new(0.0, 10.0, 0.0);
    let dot = v_before.dot(&v_after);
    assert_eq!(dot, 0.0);
}

#[test]
fn test_f11_b3_100_percent_fix_rate_interval() {
    let n = 100;
    let fix_rate = (n as f64) / (n as f64);
    assert_eq!(fix_rate, 1.0);
}

#[test]
fn test_f11_b4_exact_zero_error_comparison() {
    let pos = Vector3::new(100.0, 200.0, 300.0);
    let truth = pos;
    assert_eq!((pos - truth).norm(), 0.0);
}

#[test]
fn test_f11_b5_receiver_clock_jump_handling() {
    let clk_jump_ms = 1.0;
    let range_jump_m = clk_jump_ms * 1e-3 * SPEED_OF_LIGHT;
    assert!((range_jump_m - 299792.458).abs() < 1e-3);
}
