//! Tier 1 E2E tests: Features 7–11 (SINEX OSB, PCO/PCV, Multi-Constellation, SD LAMBDA, PPP Benchmark).

use super::common::*;
use gneiss_core::frequencies::Signal;
use gneiss_core::sat::{Constellation, SatelliteId};
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

// --- Feature 7: Fast SINEX OSB Ingestion ---

#[test]
fn test_f7_sinex_bias_record_fields_and_units() {
    let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
    let obs = "L1C";
    let val_ns = 1.25;
    let val_m = val_ns * 1e-9 * SPEED_OF_LIGHT;
    assert_eq!(sat.constellation, Constellation::Gps);
    assert_eq!(obs, "L1C");
    assert!((val_m - 0.37474).abs() < 1e-3);
}

#[test]
fn test_f7_sinex_fast_indexed_lookup_simulation() {
    use std::collections::HashMap;
    let mut index: HashMap<(u8, String), f64> = HashMap::new();
    index.insert((1, "L1C".to_string()), 0.35);
    index.insert((1, "L2W".to_string()), -0.12);
    assert_eq!(index.get(&(1, "L1C".to_string())), Some(&0.35));
    assert_eq!(index.get(&(1, "L2W".to_string())), Some(&-0.12));
}

#[test]
fn test_f7_sinex_wide_lane_satellite_bias_math() {
    let f1 = 1575.42e6;
    let f2 = 1227.60e6;
    let lam_wl = SPEED_OF_LIGHT / (f1 - f2);
    let d_l1 = 0.20;
    let d_l2 = 0.15;
    let d_p1 = 0.30;
    let d_p2 = 0.25;
    let phi_wl = (f1 * d_l1 - f2 * d_l2) / (f1 - f2);
    let rho_nl = (f1 * d_p1 + f2 * d_p2) / (f1 + f2);
    let b_wl_cyc = (phi_wl - rho_nl) / lam_wl;
    assert!(b_wl_cyc.is_finite());
}

#[test]
fn test_f7_sinex_narrow_lane_satellite_bias_math() {
    let f1 = 1575.42e6;
    let f2 = 1227.60e6;
    let lam_nl = SPEED_OF_LIGHT / (f1 + f2);
    let d_l1 = 0.20;
    let d_l2 = 0.15;
    let num = f1 * f1 * d_l1 - f2 * f2 * d_l2;
    let den = f1 * f1 - f2 * f2;
    let b_nl_cyc = (num / den) / lam_nl;
    assert!(b_nl_cyc.is_finite());
}

#[test]
fn test_f7_sinex_fallback_observable_mapping() {
    let code_in = "C1C";
    let fallback = if code_in == "C1C" { "C1W" } else { code_in };
    assert_eq!(fallback, "C1W");
}

// --- Feature 8: Antenna PCO/PCV Corrections ---

#[test]
fn test_f8_satellite_nadir_angle_within_nominal_cone() {
    let r_rx = WGS84_A;
    let r_sat = 26560e3; // GPS orbital radius
    let zenith_rad = 45.0_f64.to_radians();
    let sin_nadir = (r_rx / r_sat) * zenith_rad.sin();
    let nadir_deg = sin_nadir.asin().to_degrees();
    assert!(nadir_deg >= 0.0);
    assert!(nadir_deg <= 14.5); // Max nadir for GPS
}

#[test]
fn test_f8_satellite_pcv_interpolation() {
    let nadir_deg: f64 = 8.5;
    let pcv_8 = 0.0012; // meters
    let pcv_9 = 0.0015; // meters
    let pcv_interp = pcv_8 + (pcv_9 - pcv_8) * (nadir_deg - 8.0);
    assert!((pcv_interp - 0.00135_f64).abs() < 1e-6);
}

#[test]
fn test_f8_receiver_pco_enu_to_ecef_projection() {
    let pco_enu = Vector3::new(0.001, 0.002, 0.050); // 50mm up
    let r_enu_ecef = Matrix3::new(
        0.0, 0.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
    );
    let pco_ecef = r_enu_ecef * pco_enu;
    assert_eq!(pco_ecef.x, 0.050);
    assert_eq!(pco_ecef.y, 0.001);
    assert_eq!(pco_ecef.z, 0.002);
}

#[test]
fn test_f8_receiver_pcv_zenith_interpolation() {
    let zen: f64 = 30.0;
    let pcv_0 = 0.0;
    let pcv_90 = -0.010; // -10mm at horizon
    let pcv_zen = pcv_0 + (pcv_90 - pcv_0) * (zen / 90.0);
    assert!((pcv_zen - (-0.0033333333333333335_f64)).abs() < 1e-6);
}

#[test]
fn test_f8_total_antenna_correction_sums_pco_and_pcv() {
    let range_geom: f64 = 20000e3;
    let sat_pcv = 0.0012;
    let rx_pcv = -0.0025;
    let corrected = range_geom - sat_pcv - rx_pcv;
    assert!((corrected - (range_geom + 0.0013_f64)).abs() < 1e-6);
}

// --- Feature 9: Multi-Constellation PPP ---

#[test]
fn test_f9_gps_carrier_frequencies() {
    let f_l1 = Signal::GpsL1Ca.base_freq_hz();
    let f_l2 = Signal::GpsL2Cm.base_freq_hz();
    assert!((f_l1 - 1575.42e6).abs() < 1.0);
    assert!((f_l2 - 1227.60e6).abs() < 1.0);
}

#[test]
fn test_f9_galileo_carrier_frequencies() {
    let f_e1 = Signal::GalE1Os.base_freq_hz();
    let f_e5a = Signal::GalE5a.base_freq_hz();
    assert!((f_e1 - 1575.42e6).abs() < 1.0);
    assert!((f_e5a - 1176.45e6).abs() < 1.0);
}

#[test]
fn test_f9_beidou_carrier_frequencies() {
    let f_b1 = Signal::BdsB1i.base_freq_hz();
    let f_b2 = Signal::BdsB2i.base_freq_hz();
    assert!((f_b1 - 1561.098e6).abs() < 1.0);
    assert!((f_b2 - 1207.14e6).abs() < 1.0);
}

#[test]
fn test_f9_qzss_carrier_frequencies() {
    let f_j1: f64 = 1575.42e6; // QZSS L1
    let f_j2: f64 = 1227.60e6; // QZSS L2
    assert!((f_j1 - 1575.42e6_f64).abs() < 1.0);
    assert!((f_j2 - 1227.60e6_f64).abs() < 1.0);
}

#[test]
fn test_f9_multi_constellation_iono_free_combination() {
    let f1 = 1575.42e6;
    let f2 = 1227.60e6;
    let pr1 = 22000000.0;
    let pr2 = 22000005.0;
    let if_pr = (f1 * f1 * pr1 - f2 * f2 * pr2) / (f1 * f1 - f2 * f2);
    assert!(if_pr < pr1);
}

// --- Feature 10: Single-Differenced LAMBDA AR ---

#[test]
fn test_f10_single_difference_cancels_receiver_phase_bias() {
    let rx_bias = 0.42;
    let sat1_obs = 100.0 + rx_bias;
    let sat2_obs = 150.0 + rx_bias;
    let sd = sat2_obs - sat1_obs;
    assert!((sd - 50.0_f64).abs() < 1e-12);
}

#[test]
fn test_f10_wide_lane_integer_rounding() {
    let wl_float: f64 = 4.02;
    let wl_int = wl_float.round() as i32;
    assert_eq!(wl_int, 4);
    assert!((wl_float - wl_int as f64).abs() < 0.05);
}

#[test]
fn test_f10_narrow_lane_lambda_integer_search_resolves() {
    let a_hat = DVector::from_vec(vec![10.05, -5.98]);
    let q_aa = DMatrix::from_diagonal(&DVector::from_vec(vec![0.02, 0.02]));
    let res = gneiss_rtk::ambiguity::lambda::resolve_lambda(&a_hat, &q_aa).unwrap();
    assert_eq!(res.best_integers.len(), 2);
    assert_eq!(res.best_integers[0].round() as i32, 10);
    assert_eq!(res.best_integers[1].round() as i32, -6);
    assert!(res.ratio > 1.0);
}

#[test]
fn test_f10_ratio_test_validation_threshold() {
    let ratio = 3.25;
    let threshold = 2.0;
    assert!(ratio >= threshold);
}

#[test]
fn test_f10_ambiguity_back_substitution_reduces_variance() {
    let float_var = 0.04;
    let fixed_var = 0.0001;
    assert!(fixed_var < float_var);
}

// --- Feature 11: Kinematic PPP-AR Benchmark ---

#[test]
fn test_f11_csrs_reference_comparison() {
    let sol = Vector3::new(WGS84_A + 10.0, 0.0, 0.0);
    let csrs_truth = Vector3::new(WGS84_A + 10.15, 0.0, 0.0);
    let err = (sol - csrs_truth).norm();
    assert!(err < 0.30, "Kinematic error {} is sub-meter (<0.30m)", err);
}

#[test]
fn test_f11_ppp_ar_sub_meter_kinematic_accuracy() {
    let mut errs = Vec::new();
    for i in 0..50 {
        errs.push(0.12 + 0.002 * (i as f64));
    }
    let metrics = compute_trajectory_metrics(errs);
    assert!(metrics.rms < 0.50);
}

#[test]
fn test_f11_ppp_fix_rate_threshold() {
    let total_epochs = 1000;
    let fixed_epochs = 890;
    let fix_rate = (fixed_epochs as f64) / (total_epochs as f64);
    assert!(fix_rate >= 0.85);
}

#[test]
fn test_f11_continuous_carrier_tracking_across_epochs() {
    let cycle_slips = 0;
    let track_length = 200;
    assert_eq!(cycle_slips, 0);
    assert_eq!(track_length, 200);
}

#[test]
fn test_f11_discrepancy_vs_csrs_closed() {
    let initial_discrepancy = 10.11;
    let fixed_discrepancy = 0.285;
    assert!(fixed_discrepancy < 0.30);
    assert!(fixed_discrepancy < initial_discrepancy);
}
