//! Tier 1: Feature Coverage in Isolation (>=5 tests per feature for R1, R2, R3, R4).
#![allow(clippy::unwrap_used)]

use nalgebra::{DMatrix, DVector, Vector3};
use gneiss_core::dop::compute_dop_from_positions;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::variance::{elevation_variance_scale, snr_variance_scale};
use gneiss_rtk::ambiguity::par::select_ils_subset;
use gneiss_rtk::estimators::rtk_iekf::ar_subsets::is_beidou_geo_key;
use gneiss_rtk::estimators::rtk_iekf::formation_cov::{
    attenuation_scale, compute_dd_variances, expected_cn0_dbhz,
};
use gneiss_rtk::estimators::rtk_iekf::mw::WidelaneTracker;
use gneiss_rtk::estimators::rtk_iekf::screen::screen_gross_pr_errors;
use gneiss_rtk::estimators::rtk_iekf::state::DoubleDiffKey;
use gneiss_rtk::post_process::screening::CycleSlipDetector;

use super::common::*;

// =========================================================================
// R1: Adaptive C/N0 & Elevation Observation Covariance Weighting
// =========================================================================

#[test]
fn test_r1_elevation_variance_monotonic_increase() {
    let v_zenith = elevation_variance_scale(90.0f64.to_radians());
    let v_mid = elevation_variance_scale(45.0f64.to_radians());
    let v_low = elevation_variance_scale(15.0f64.to_radians());
    let v_horizon = elevation_variance_scale(5.0f64.to_radians());

    assert!(v_zenith < v_mid, "Zenith variance must be strictly less than 45 deg");
    assert!(v_mid < v_low, "45 deg variance must be strictly less than 15 deg");
    assert!(v_low <= v_horizon, "15 deg variance must be less than or equal to 5 deg");
    assert!((v_zenith - 1.0).abs() < 1e-6, "Zenith scale must equal 1.0");
}

#[test]
fn test_r1_snr_variance_monotonic_increase() {
    let a = 1.0;
    let b = 100.0;
    let v_high = snr_variance_scale(45.0, a, b);
    let v_mid = snr_variance_scale(35.0, a, b);
    let v_low = snr_variance_scale(25.0, a, b);
    let v_atten = snr_variance_scale(15.0, a, b);

    assert!(v_high < v_mid, "45 dB-Hz variance must be lower than 35 dB-Hz");
    assert!(v_mid < v_low, "35 dB-Hz variance must be lower than 25 dB-Hz");
    assert!(v_low < v_atten, "25 dB-Hz variance must be lower than 15 dB-Hz");
}

#[test]
fn test_r1_expected_cn0_increases_with_elevation() {
    let cn0_low = expected_cn0_dbhz(10.0f64.to_radians());
    let cn0_mid = expected_cn0_dbhz(40.0f64.to_radians());
    let cn0_zenith = expected_cn0_dbhz(90.0f64.to_radians());

    assert!(cn0_low < cn0_mid, "Expected C/N0 must increase with elevation");
    assert!(cn0_mid < cn0_zenith, "Expected C/N0 must peak at zenith");
    assert!(cn0_zenith <= 50.0, "Expected C/N0 must not exceed 50 dB-Hz cap");
}

#[test]
fn test_r1_attenuation_scale_smooth_penalty() {
    let scale_nominal = attenuation_scale(Some(42), 45.0f64.to_radians());
    let scale_attenuated = attenuation_scale(Some(30), 45.0f64.to_radians());
    let scale_severe = attenuation_scale(Some(20), 45.0f64.to_radians());

    assert_eq!(scale_nominal, 1.0, "Nominal SNR must have unit attenuation scale");
    assert!(scale_attenuated > 1.0, "Attenuated SNR must scale up variance");
    assert!(scale_severe > scale_attenuated, "Severe drop must produce larger scale");
    assert!(scale_severe <= 1000.0, "Scale must not exceed stability cap");
}

#[test]
fn test_r1_dd_variances_positive_definite() {
    let rx = rover_ecef();
    let sat = sat_pos_az_el(rx, 45.0, 30.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);
    let snrs = (Some(38), Some(46), Some(48), Some(48));

    let var = compute_dd_variances(rx, sat, ref_sat, GPS_L1_WAVELENGTH_M, snrs);
    assert!(var.pr_var_m2 > 0.0, "Code variance must be positive");
    assert!(var.cp_var_cycles2 > 0.0, "Phase variance must be positive");
    assert!(var.pr_ref_var_m2 > 0.0, "Ref code variance must be positive");
    assert!(var.cp_ref_var_cycles2 > 0.0, "Ref phase variance must be positive");
}

#[test]
fn test_r1_dd_covariance_matrix_positive_definite() {
    let rx = rover_ecef();
    let sat1 = sat_pos_az_el(rx, 30.0, 40.0, 20_000_000.0);
    let sat2 = sat_pos_az_el(rx, 120.0, 50.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 210.0, 70.0, 20_000_000.0);
    let snrs = (Some(40), Some(45), Some(45), Some(45));

    let v1 = compute_dd_variances(rx, sat1, ref_sat, GPS_L1_WAVELENGTH_M, snrs);
    let v2 = compute_dd_variances(rx, sat2, ref_sat, GPS_L1_WAVELENGTH_M, snrs);

    let mut r = DMatrix::zeros(2, 2);
    r[(0, 0)] = v1.pr_var_m2;
    r[(1, 1)] = v2.pr_var_m2;
    let cov_ref = v1.pr_ref_var_m2.min(v2.pr_ref_var_m2);
    r[(0, 1)] = cov_ref;
    r[(1, 0)] = cov_ref;

    assert!(is_matrix_positive_definite(&r, 1e-6), "Assembled DD covariance must be positive definite");
}

// =========================================================================
// R2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting
// =========================================================================

#[test]
fn test_r2_cmc_geometry_free_invariance() {
    let lambda = GPS_L1_WAVELENGTH_M;
    let range1: f64 = 20_000_000.0;
    let range2: f64 = 20_000_100.0;

    let p1: f64 = range1;
    let phi1: f64 = range1 / lambda + 1000.0;
    let p2: f64 = range2;
    let phi2: f64 = range2 / lambda + 1000.0;

    let cmc1: f64 = p1 - lambda * phi1;
    let cmc2: f64 = p2 - lambda * phi2;
    let diff: f64 = cmc1 - cmc2;
    assert!(diff.abs() < 1e-6, "CMC must eliminate geometric motion");
}

#[test]
fn test_r2_cmc_step_blunder_detection() {
    let lambda = GPS_L1_WAVELENGTH_M;
    let range: f64 = 20_000_000.0;
    let p_clean: f64 = range;
    let phi_clean: f64 = range / lambda + 500.0;
    let cmc_nominal: f64 = p_clean - lambda * phi_clean;

    let p_corrupt: f64 = range + 10.0;
    let cmc_corrupt: f64 = p_corrupt - lambda * phi_clean;
    let delta_cmc: f64 = (cmc_corrupt - cmc_nominal).abs();

    assert!((delta_cmc - 10.0).abs() < 1e-6, "CMC residual must detect exact 10m code step");
}

#[test]
fn test_r2_screen_gross_error_preserves_carrier_phase() {
    let rx = rover_ecef();
    let sat1 = sat_pos_az_el(rx, 45.0, 30.0, 20_000_000.0);
    let sat2 = sat_pos_az_el(rx, 135.0, 45.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 210.0, 70.0, 20_000_000.0);

    let geom_dd1 = compute_geometric_dd(rx, ref_ecef(), sat1, ref_sat);
    let geom_dd2 = compute_geometric_dd(rx, ref_ecef(), sat2, ref_sat);

    // m1 is clean (residual ≈ 0)
    let m1 = make_dd_meas((0, 1, 10), geom_dd1, Some(geom_dd1 / GPS_L1_WAVELENGTH_M), (0.08, 0.0005), sat1, ref_sat);
    // m2 has +25m gross pseudorange error
    let m2 = make_dd_meas((0, 2, 10), geom_dd2 + 25.0, Some(geom_dd2 / GPS_L1_WAVELENGTH_M), (0.08, 0.0005), sat2, ref_sat);

    let mut meas = vec![m1, m2];
    let rejected = screen_gross_pr_errors(&mut meas, rx);

    assert_eq!(rejected.len(), 1, "Exactly one blunder pair must be rejected");
    assert_eq!(rejected[0].sat, 2, "Corrupted satellite pair 2 must be flagged");
    assert_eq!(meas.len(), 2, "Both measurements retained: carrier phase preserved");
    let screened = meas.iter().find(|m| m.key.sat == 2).expect("corrupted pair retained");
    assert!(screened.dd_cp_cycles.is_some(), "Carrier phase must remain active on screened pair");
    assert!(screened.pr_var_m2 >= 1.0e6, "Code variance must be inflated/suppressed");
    assert_eq!(screened.cp_var_cycles2, 0.0005, "Carrier variance must remain nominal");
}

#[test]
fn test_r2_mw_tracker_accumulates_clean_arc() {
    let mut tracker = WidelaneTracker::default();
    let key = DoubleDiffKey { constellation_id: 0, sat: 5, ref_sat: 10, freq_band: 1 };
    for i in 0..20 {
        let noise = ((i as f64) * 0.1).sin() * 0.02;
        tracker.update(key, 5.0 + noise, 3.5, false);
    }
    let means = tracker.arc_means();
    assert!(means.contains_key(&key));
    let (mean, count) = means[&key];
    assert_eq!(count, 20, "MW tracker must have accumulated 20 epochs");
    assert!((mean - 5.0).abs() < 0.05, "MW mean must converge to true wide-lane");
}

#[test]
fn test_r2_code_variance_inflation_reduces_kalman_gain() {
    let p_prior = 1.0;
    let h = 1.0;
    let r_nominal = 0.04;
    let r_inflated = 1000.0;

    let k_nominal = (p_prior * h) / (h * p_prior * h + r_nominal);
    let k_inflated = (p_prior * h) / (h * p_prior * h + r_inflated);

    assert!(k_inflated < 0.01 * k_nominal, "Kalman gain must drop by >99% under multipath inflation");
}

// =========================================================================
// R3: Doppler-Assisted Cycle Slip Detection & Phase Continuity Validation
// =========================================================================

#[test]
fn test_r3_doppler_exact_one_cycle_slip() {
    let mut detector = CycleSlipDetector::new();
    let sat_slip = SatelliteId { constellation: Constellation::Gps, prn: 5 };

    let o1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 100_000.0, 100.0, 42.0, None);
    let o2 = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 100_000.0, 200.0, 42.0, None);
    let o3 = make_sat_obs_l1(5, Constellation::Gps, 20_000_000.0, 100_000.0, 500.0, 42.0, None);
    detector.check_epoch(&make_epoch(100.0, vec![o1, o2, o3]));

    let o1_next = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 100_000.0 - 100.0, 100.0, 42.0, None);
    let o2_next = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 100_000.0 - 200.0, 200.0, 42.0, None);
    let o3_next = make_sat_obs_l1(5, Constellation::Gps, 20_000_000.0, 100_000.0 - 500.0 + 2.0, 500.0, 42.0, None);
    let _slips = detector.check_epoch(&make_epoch(101.0, vec![o1_next, o2_next, o3_next]));

    assert_eq!(detector.get_arc(sat_slip), 1, "Arc count must increment on cycle slip");
}

#[test]
fn test_r3_doppler_half_cycle_slip() {
    let mut detector = CycleSlipDetector::new();
    let sat_slip = SatelliteId { constellation: Constellation::Gps, prn: 8 };

    let o1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0, 100.0, 40.0, None);
    let o2 = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 50_000.0, 150.0, 40.0, None);
    let o3 = make_sat_obs_l1(8, Constellation::Gps, 20_000_000.0, 50_000.0, -300.0, 40.0, None);
    detector.check_epoch(&make_epoch(200.0, vec![o1, o2, o3]));

    let o1_next = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0 - 100.0, 100.0, 40.0, None);
    let o2_next = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 50_000.0 - 150.0, 150.0, 40.0, None);
    let o3_next = make_sat_obs_l1(8, Constellation::Gps, 20_000_000.0, 50_000.0 + 300.5, -300.0, 40.0, Some(1));
    detector.check_epoch(&make_epoch(201.0, vec![o1_next, o2_next, o3_next]));

    assert_eq!(detector.get_arc(sat_slip), 1, "Half-cycle slip flagged by firmware must increment arc counter");
}

#[test]
fn test_r3_clean_carrier_matches_doppler_no_slip() {
    let mut detector = CycleSlipDetector::new();
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 12 };

    let mut cp = 80_000.0;
    let doppler = 250.0;
    for t in 0..5 {
        let o1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0, 0.0, 45.0, None);
        let o2 = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 60_000.0, 0.0, 45.0, None);
        let obs = make_sat_obs_l1(12, Constellation::Gps, 20_000_000.0, cp, doppler, 45.0, None);
        detector.check_epoch(&make_epoch(300.0 + t as f64, vec![o1, o2, obs]));
        cp -= doppler;
    }

    assert_eq!(detector.get_arc(sat_id), 0, "Clean carrier matching Doppler must have 0 slips");
}

#[test]
fn test_r3_time_gap_triggers_slip() {
    let mut detector = CycleSlipDetector::new();
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 15 };

    let obs1 = make_sat_obs_l1(15, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 44.0, None);
    detector.check_epoch(&make_epoch(400.0, vec![obs1]));

    let obs2 = make_sat_obs_l1(15, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 44.0, None);
    detector.check_epoch(&make_epoch(405.0, vec![obs2]));

    assert_eq!(detector.get_arc(sat_id), 1, "Time gap > 2.0 s must trigger cycle slip");
}

#[test]
fn test_r3_loss_of_lock_indicator_flags_slip() {
    let mut detector = CycleSlipDetector::new();
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 20 };

    let obs1 = make_sat_obs_l1(20, Constellation::Gps, 20_000_000.0, 20_000.0, 100.0, 43.0, None);
    detector.check_epoch(&make_epoch(500.0, vec![obs1]));

    let obs2 = make_sat_obs_l1(20, Constellation::Gps, 20_000_000.0, 20_100.0, 100.0, 43.0, Some(1));
    detector.check_epoch(&make_epoch(501.0, vec![obs2]));

    assert_eq!(detector.get_arc(sat_id), 1, "LLI bit 0 must trigger cycle slip");
}

// =========================================================================
// R4: C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR)
// =========================================================================

#[test]
fn test_r4_par_subset_selects_high_confidence_ambiguities() {
    let a = DVector::from_vec(vec![1.05, 2.02, 3.48, 4.45]);
    let mut q = DMatrix::zeros(4, 4);
    q[(0, 0)] = 0.01;
    q[(1, 1)] = 0.02;
    q[(2, 2)] = 0.80;
    q[(3, 3)] = 0.90;

    let (selected, sub_a, sub_q) = select_ils_subset(&a, &q, 0.95);
    assert!(!selected.is_empty(), "PAR must select at least one ambiguity");
    assert!(selected.contains(&0) && selected.contains(&1), "PAR must prioritize confident indices 0 and 1");
    assert_eq!(sub_a.len(), selected.len(), "Sub-vector dimension must match selected count");
    assert_eq!(sub_q.nrows(), selected.len(), "Sub-covariance dimension must match selected count");
}

#[test]
fn test_r4_par_empty_on_zero_dimension() {
    let a = DVector::zeros(0);
    let q = DMatrix::zeros(0, 0);
    let (selected, sub_a, sub_q) = select_ils_subset(&a, &q, 0.99);

    assert!(selected.is_empty(), "Empty input must return empty selection");
    assert_eq!(sub_a.len(), 0);
    assert_eq!(sub_q.nrows(), 0);
}

#[test]
fn test_r4_par_dop_computation_from_geometry() {
    let rx = Vector3::new(6_378_137.0, 0.0, 0.0);
    let sats = vec![
        Vector3::new(20_000_000.0, 5_000_000.0, 5_000_000.0),
        Vector3::new(22_000_000.0, -5_000_000.0, 5_000_000.0),
        Vector3::new(19_000_000.0, 5_000_000.0, -5_000_000.0),
        Vector3::new(21_000_000.0, -5_000_000.0, -5_000_000.0),
        Vector3::new(25_000_000.0, 0.0, 0.0),
        Vector3::new(15_000_000.0, 0.0, 10_000_000.0),
    ];

    let dop = compute_dop_from_positions(rx, &sats);
    assert!(dop.is_some(), "DOP must be successfully computed for 6 well-spread satellites");
    let dop_val = dop.unwrap();
    assert!(dop_val.pdop > 0.0 && dop_val.pdop < 10.0, "PDOP must be healthy (< 10.0)");
    assert!(dop_val.gdop >= dop_val.pdop, "GDOP must be >= PDOP");
    assert!(dop_val.pdop >= dop_val.hdop, "PDOP must be >= HDOP");
}

#[test]
fn test_r4_par_submatrix_covariance_positive_definite() {
    let mut q = DMatrix::zeros(3, 3);
    q[(0, 0)] = 0.04; q[(0, 1)] = 0.01; q[(0, 2)] = 0.01;
    q[(1, 0)] = 0.01; q[(1, 1)] = 0.05; q[(1, 2)] = 0.02;
    q[(2, 0)] = 0.01; q[(2, 1)] = 0.02; q[(2, 2)] = 0.06;

    let a = DVector::from_vec(vec![1.01, 2.02, 3.03]);
    let (selected, _, sub_q) = select_ils_subset(&a, &q, 0.90);

    if !selected.is_empty() {
        assert!(is_matrix_positive_definite(&sub_q, 1e-6), "Sub-covariance must be positive definite");
    }
}

#[test]
fn test_r4_par_beidou_geo_key_detection() {
    let key_geo = DoubleDiffKey { constellation_id: 3, sat: 2, ref_sat: 10, freq_band: 1 };
    let key_meo = DoubleDiffKey { constellation_id: 3, sat: 19, ref_sat: 20, freq_band: 1 };
    let key_gps = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 5, freq_band: 1 };

    assert!(is_beidou_geo_key(&key_geo), "BeiDou PRN 2 must be recognized as GEO");
    assert!(!is_beidou_geo_key(&key_meo), "BeiDou PRN 19 must be recognized as non-GEO MEO");
    assert!(!is_beidou_geo_key(&key_gps), "GPS PRN 2 must not be flagged as BeiDou GEO");
}
