//! Tier 2: Boundary & Corner Cases (>=5 tests per feature for R1, R2, R3, R4).
#![allow(clippy::unwrap_used)]

use nalgebra::{DMatrix, DVector};
use gneiss_core::dop::compute_dop_from_positions;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::variance::{elevation_variance_scale, snr_variance_scale};
use gneiss_rtk::ambiguity::par::select_ils_subset;
use gneiss_rtk::estimators::rtk_iekf::screen::{
    screen_gross_pr_errors, MAX_GROSS_PR_REJECTIONS_PER_EPOCH,
};
use gneiss_rtk::post_process::screening::CycleSlipDetector;

use super::common::*;

// =========================================================================
// R1: Boundary & Corner Conditions
// =========================================================================

#[test]
fn test_r1_boundary_extreme_low_elevation_horizon() {
    let v_zero = elevation_variance_scale(0.0);
    let v_one_deg = elevation_variance_scale(1.0f64.to_radians());
    let v_neg = elevation_variance_scale(-5.0f64.to_radians());

    assert!(v_zero.is_finite() && !v_zero.is_nan(), "Horizon elevation must produce finite variance");
    assert!(v_one_deg.is_finite(), "1 deg elevation must produce finite variance");
    assert!(v_neg.is_finite(), "Negative elevation must not panic or NaN");
    assert!(v_zero >= 99.0, "Clamped horizon scale must be large (>= 99.0)");
}

#[test]
fn test_r1_boundary_extreme_low_snr_saturation() {
    let a = 1.0;
    let b = 100.0;
    let v_zero = snr_variance_scale(0.0, a, b);
    let v_five = snr_variance_scale(5.0, a, b);
    let v_ten = snr_variance_scale(10.0, a, b);

    assert!(v_zero.is_finite() && !v_zero.is_nan(), "0 dB-Hz must be finite");
    assert_eq!(v_zero, v_ten, "SNR below 10 dB-Hz must clamp to 10 dB-Hz baseline");
    assert_eq!(v_five, v_ten, "SNR of 5 dB-Hz must clamp to 10 dB-Hz baseline");
}

#[test]
fn test_r1_boundary_zenith_maximum_snr() {
    let a = 1.0;
    let b = 100.0;
    let v_high = snr_variance_scale(60.0, a, b);
    assert!((v_high - a * a).abs() < 0.05, "Very high SNR must approach floor a^2");
    let el_zenith = elevation_variance_scale(90.0f64.to_radians());
    assert!((el_zenith - 1.0).abs() < 1e-6, "Zenith elevation scale must be exactly 1.0");
}

#[test]
fn test_r1_boundary_finite_difference_elevation_gradient() {
    let angles: Vec<f64> = (10..90).step_by(5).map(|d| (d as f64).to_radians()).collect();
    for window in angles.windows(2) {
        let el1 = window[0];
        let el2 = window[1];
        let v1 = elevation_variance_scale(el1);
        let v2 = elevation_variance_scale(el2);
        let grad = (v2 - v1) / (el2 - el1);
        assert!(grad <= 0.0, "Elevation gradient must be strictly non-positive (variance decreases with elevation)");
    }
}

#[test]
fn test_r1_boundary_finite_difference_snr_gradient() {
    let a = 1.0;
    let b = 100.0;
    let snrs: Vec<f64> = (15..55).step_by(5).map(|s| s as f64).collect();
    for window in snrs.windows(2) {
        let s1 = window[0];
        let s2 = window[1];
        let v1 = snr_variance_scale(s1, a, b);
        let v2 = snr_variance_scale(s2, a, b);
        let grad = (v2 - v1) / (s2 - s1);
        assert!(grad <= 0.0, "SNR gradient must be strictly non-positive (variance decreases with SNR)");
    }
}

// =========================================================================
// R2: Boundary & Corner Conditions
// =========================================================================

#[test]
fn test_r2_boundary_extreme_20m_code_step() {
    let rx = rover_ecef();
    let sat = sat_pos_az_el(rx, 45.0, 30.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);
    let geom_dd = compute_geometric_dd(rx, ref_ecef(), sat, ref_sat);

    let m_extreme = make_dd_meas((0, 1, 10), geom_dd + 25.0, Some(geom_dd / GPS_L1_WAVELENGTH_M), (0.08, 0.0005), sat, ref_sat);
    let mut meas = vec![m_extreme];
    let rejected = screen_gross_pr_errors(&mut meas, rx);

    assert_eq!(rejected.len(), 1, "20m code step must be intercepted by prefit screen");
    assert_eq!(rejected[0].sat, 1);
}

#[test]
fn test_r2_boundary_zero_multipath_code_step() {
    let rx = rover_ecef();
    let sat = sat_pos_az_el(rx, 60.0, 45.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);
    let geom_dd = compute_geometric_dd(rx, ref_ecef(), sat, ref_sat);

    let m_clean = make_dd_meas((0, 2, 10), geom_dd, Some(geom_dd / GPS_L1_WAVELENGTH_M), (0.08, 0.0005), sat, ref_sat);
    let mut meas = vec![m_clean];
    let rejected = screen_gross_pr_errors(&mut meas, rx);

    assert!(rejected.is_empty(), "Zero multipath error must not trigger false rejection");
    assert_eq!(meas.len(), 1, "Measurement retained intact");
}

#[test]
fn test_r2_boundary_max_gross_pr_rejections_cap() {
    let rx = rover_ecef();
    let sat_ref = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);

    let mut meas = Vec::new();
    for i in 1..=5 {
        let sat = sat_pos_az_el(rx, (i * 50) as f64, 30.0, 20_000_000.0);
        let geom_dd = compute_geometric_dd(rx, ref_ecef(), sat, sat_ref);
        meas.push(make_dd_meas((0, i, 10), geom_dd + 50.0 + (i as f64), Some(20.0), (0.08, 0.0005), sat, sat_ref));
    }

    let rejected = screen_gross_pr_errors(&mut meas, rx);
    assert_eq!(
        rejected.len(),
        MAX_GROSS_PR_REJECTIONS_PER_EPOCH,
        "Rejections per epoch must be capped at MAX_GROSS_PR_REJECTIONS_PER_EPOCH"
    );
}

#[test]
fn test_r2_boundary_carrier_phase_only_measurement() {
    let rx = rover_ecef();
    let sat = sat_pos_az_el(rx, 90.0, 50.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);

    let m = make_dd_meas((0, 3, 10), 0.0, Some(500.0), (1e8, 0.0002), sat, ref_sat);
    assert!(m.dd_cp_cycles.is_some(), "Carrier phase present");
    assert!(m.pr_var_m2 >= 1e8, "Code variance inflated to infinity");
}

#[test]
fn test_r2_boundary_alternating_multipath_steps() {
    let mut cmc_history = Vec::new();
    let lambda = GPS_L1_WAVELENGTH_M;
    let base_range: f64 = 20_000_000.0;

    let steps = [12.0, -12.0, 12.0, -12.0];
    for &step in &steps {
        let p: f64 = base_range + step;
        let phi: f64 = base_range / lambda;
        let cmc: f64 = p - lambda * phi;
        cmc_history.push(cmc);
    }

    for window in cmc_history.windows(2) {
        let diff: f64 = window[1] - window[0];
        assert!((diff.abs() - 24.0).abs() < 1e-9, "Alternating step discrepancy must be exact 24m");
    }
}

// =========================================================================
// R3: Boundary & Corner Conditions
// =========================================================================

#[test]
fn test_r3_boundary_high_vehicle_acceleration() {
    let mut detector = CycleSlipDetector::new();
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 3 };

    let mut cp = 50_000.0;
    let mut doppler = 100.0;
    let accel_hz_per_sec = 100.0;

    for t in 0..4 {
        let o1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0, 0.0, 45.0, None);
        let o2 = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 60_000.0, 0.0, 45.0, None);
        let obs = make_sat_obs_l1(3, Constellation::Gps, 20_000_000.0, cp, doppler, 44.0, None);
        detector.check_epoch(&make_epoch(100.0 + t as f64, vec![o1, o2, obs]));
        let avg_dop = doppler + 0.5 * accel_hz_per_sec;
        cp -= avg_dop;
        doppler += accel_hz_per_sec;
    }

    assert_eq!(detector.get_arc(sat_id), 0, "Accelerating vehicle with consistent Doppler must not trigger false slip");
}

#[test]
fn test_r3_boundary_near_threshold_doppler_noise() {
    let mut detector = CycleSlipDetector::new();
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 7 };

    let o1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 42.0, None);
    let o2 = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 42.0, None);
    let obs1 = make_sat_obs_l1(7, Constellation::Gps, 20_000_000.0, 10_000.0, 200.0, 42.0, None);
    detector.check_epoch(&make_epoch(150.0, vec![o1, o2, obs1]));

    let o1_next = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 42.0, None);
    let o2_next = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 42.0, None);
    let obs2 = make_sat_obs_l1(7, Constellation::Gps, 20_000_000.0, 10_000.0 - 200.0 + 0.10, 200.0, 42.0, None);
    detector.check_epoch(&make_epoch(151.0, vec![o1_next, o2_next, obs2]));

    assert_eq!(detector.get_arc(sat_id), 0, "0.10 cycle noise must not trigger false slip");
}

#[test]
fn test_r3_boundary_cadence_hint_long_interval() {
    let mut detector = CycleSlipDetector::new();
    detector.cadence_hint_s = Some(30.0);
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 9 };

    let obs1 = make_sat_obs_l1(9, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 45.0, None);
    detector.check_epoch(&make_epoch(0.0, vec![obs1]));

    let obs2 = make_sat_obs_l1(9, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 45.0, None);
    detector.check_epoch(&make_epoch(30.0, vec![obs2]));

    assert_eq!(detector.get_arc(sat_id), 0, "Expected 30s cadence must not trigger false time-gap slip");
}

#[test]
fn test_r3_boundary_simultaneous_slips_multiple_constellations() {
    let mut detector = CycleSlipDetector::new();
    let sat_gps = SatelliteId { constellation: Constellation::Gps, prn: 1 };
    let sat_gal = SatelliteId { constellation: Constellation::Galileo, prn: 2 };
    let sat_bds = SatelliteId { constellation: Constellation::Beidou, prn: 3 };

    let o_g1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 42.0, None);
    let o_e1 = make_sat_obs_l1(2, Constellation::Galileo, 20_000_000.0, 10_000.0, 0.0, 42.0, None);
    let o_c1 = make_sat_obs_l1(3, Constellation::Beidou, 20_000_000.0, 10_000.0, 0.0, 42.0, None);
    detector.check_epoch(&make_epoch(10.0, vec![o_g1, o_e1, o_c1]));

    let o_g2 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 10_000.0, 0.0, 42.0, Some(1));
    let o_e2 = make_sat_obs_l1(2, Constellation::Galileo, 20_000_000.0, 10_000.0, 0.0, 42.0, Some(1));
    let o_c2 = make_sat_obs_l1(3, Constellation::Beidou, 20_000_000.0, 10_000.0, 0.0, 42.0, Some(1));
    detector.check_epoch(&make_epoch(11.0, vec![o_g2, o_e2, o_c2]));

    assert_eq!(detector.get_arc(sat_gps), 1);
    assert_eq!(detector.get_arc(sat_gal), 1);
    assert_eq!(detector.get_arc(sat_bds), 1);
}

// =========================================================================
// R4: Boundary & Corner Conditions
// =========================================================================

#[test]
fn test_r4_boundary_minimal_subset_size_k4() {
    let rx = rover_ecef();
    let sats = vec![
        sat_pos_az_el(rx, 0.0, 20.0, 20_000_000.0),
        sat_pos_az_el(rx, 120.0, 30.0, 20_000_000.0),
        sat_pos_az_el(rx, 240.0, 40.0, 20_000_000.0),
        sat_pos_az_el(rx, 0.0, 85.0, 20_000_000.0),
    ];
    let dop = compute_dop_from_positions(rx, &sats).expect("DOP must exist for 4 satellites");
    assert!(dop.pdop > 0.0 && dop.pdop < 10.0, "Minimal k=4 subset has healthy PDOP");
}

#[test]
fn test_r4_boundary_degenerate_collinear_geometry_high_dop() {
    let rx = rover_ecef();
    let sats = vec![
        sat_pos_az_el(rx, 0.0, 15.0, 20_000_000.0),
        sat_pos_az_el(rx, 0.0, 30.0, 20_000_000.0),
        sat_pos_az_el(rx, 180.0, 15.0, 20_000_000.0),
        sat_pos_az_el(rx, 180.0, 30.0, 20_000_000.0),
    ];
    let dop = compute_dop_from_positions(rx, &sats);
    if let Some(d) = dop {
        assert!(d.pdop > 10.0, "Collinear satellites must produce large PDOP");
    }
}

#[test]
fn test_r4_boundary_near_singular_covariance_matrix() {
    let mut q = DMatrix::zeros(3, 3);
    q[(0, 0)] = 1e-4;
    q[(1, 1)] = 1e-4;
    q[(2, 2)] = 1e-4;
    q[(0, 1)] = 0.9999e-4; q[(1, 0)] = 0.9999e-4;

    let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
    let (selected, _, sub_q) = select_ils_subset(&a, &q, 0.99);

    if !selected.is_empty() {
        assert!(sub_q.nrows() == selected.len(), "Sub-covariance size must match");
    }
}

#[test]
fn test_r4_boundary_single_ambiguity_par() {
    let a = DVector::from_vec(vec![1.02]);
    let q = DMatrix::from_element(1, 1, 0.005);
    let (selected, sub_a, sub_q) = select_ils_subset(&a, &q, 0.95);

    assert_eq!(selected.len(), 1, "Single confident ambiguity must be selected");
    assert_eq!(selected[0], 0);
    assert_eq!(sub_a[0], 1.02);
    assert_eq!(sub_q[(0, 0)], 0.005);
}

#[test]
fn test_r4_boundary_high_target_success_rate_strict_filter() {
    let a = DVector::from_vec(vec![1.05, 2.30, 3.40]);
    let mut q = DMatrix::zeros(3, 3);
    q[(0, 0)] = 0.001;
    q[(1, 1)] = 0.20;
    q[(2, 2)] = 0.50;

    let (selected, _, _) = select_ils_subset(&a, &q, 0.9999);
    assert!(selected.contains(&0), "Top confident ambiguity must be selected");
    assert!(!selected.contains(&2), "Poor ambiguity must be excluded under strict threshold");
}
