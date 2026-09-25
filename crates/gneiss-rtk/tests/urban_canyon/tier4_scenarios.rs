//! Tier 4: Real-World Urban Canyon Mission Scenarios (>=5 scenarios).
#![allow(clippy::unwrap_used)]

use nalgebra::Vector3;
use gneiss_core::dop::compute_dop_from_positions;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_rtk::estimators::rtk_iekf::formation_cov::compute_dd_variances;
use gneiss_rtk::estimators::rtk_iekf::screen::screen_gross_pr_errors;
use gneiss_rtk::estimators::rtk_iekf::update::DoubleDiffMeasurement;
use gneiss_rtk::post_process::screening::CycleSlipDetector;

use super::common::*;

fn build_shinjuku_meas(rx: Vector3<f64>, sat_ref: Vector3<f64>) -> Vec<DoubleDiffMeasurement> {
    let s_zen1 = sat_pos_az_el(rx, 40.0, 70.0, 20_000_000.0);
    let s_zen2 = sat_pos_az_el(rx, 220.0, 75.0, 20_000_000.0);
    let s_ref1 = sat_pos_az_el(rx, 90.0, 35.0, 20_000_000.0);
    let s_ref2 = sat_pos_az_el(rx, 270.0, 40.0, 20_000_000.0);
    let s_low1 = sat_pos_az_el(rx, 10.0, 12.0, 20_000_000.0);
    let s_low2 = sat_pos_az_el(rx, 190.0, 12.0, 20_000_000.0);

    let g_zen1 = compute_geometric_dd(rx, ref_ecef(), s_zen1, sat_ref);
    let g_zen2 = compute_geometric_dd(rx, ref_ecef(), s_zen2, sat_ref);
    let g_ref1 = compute_geometric_dd(rx, ref_ecef(), s_ref1, sat_ref);
    let g_ref2 = compute_geometric_dd(rx, ref_ecef(), s_ref2, sat_ref);
    let g_low1 = compute_geometric_dd(rx, ref_ecef(), s_low1, sat_ref);
    let g_low2 = compute_geometric_dd(rx, ref_ecef(), s_low2, sat_ref);

    vec![
        make_dd_meas((0, 1, 10), g_zen1, Some(g_zen1 / GPS_L1_WAVELENGTH_M), (0.04, 0.0003), s_zen1, sat_ref),
        make_dd_meas((0, 2, 10), g_zen2, Some(g_zen2 / GPS_L1_WAVELENGTH_M), (0.04, 0.0003), s_zen2, sat_ref),
        make_dd_meas((0, 3, 10), g_ref1 + 25.0, Some(g_ref1 / GPS_L1_WAVELENGTH_M), (0.50, 0.0010), s_ref1, sat_ref),
        make_dd_meas((0, 4, 10), g_ref2 + 25.0, Some(g_ref2 / GPS_L1_WAVELENGTH_M), (0.50, 0.0010), s_ref2, sat_ref),
        make_dd_meas((0, 5, 10), g_low1, Some(g_low1 / GPS_L1_WAVELENGTH_M), (4.00, 0.0050), s_low1, sat_ref),
        make_dd_meas((0, 6, 10), g_low2, Some(g_low2 / GPS_L1_WAVELENGTH_M), (4.00, 0.0050), s_low2, sat_ref),
    ]
}

/// Scenario 1: Tokyo Shinjuku Skyscraper Canyon (Deep Multipath & Specular Reflections).
#[test]
fn test_tier4_scenario1_tokyo_shinjuku_skyscraper() {
    let rx = rover_ecef();
    let sat_ref = sat_pos_az_el(rx, 180.0, 80.0, 20_000_000.0);
    let mut meas = build_shinjuku_meas(rx, sat_ref);
    let rejected = screen_gross_pr_errors(&mut meas, rx);

    assert_eq!(rejected.len(), 2, "Both skyscraper reflected blunders must be caught");
    assert!(rejected.iter().any(|k| k.sat == 3));
    assert!(rejected.iter().any(|k| k.sat == 4));

    assert!(meas.iter().any(|m| m.key.sat == 1 && m.dd_cp_cycles.is_some()));
    assert!(meas.iter().any(|m| m.key.sat == 2 && m.dd_cp_cycles.is_some()));
}

/// Scenario 2: Hong Kong Whampoa High-Rise Urban Canyon (Frequent Half-Cycle and 1-Cycle Slips).
#[test]
fn test_tier4_scenario2_hong_kong_whampoa_frequent_slips() {
    let mut detector = CycleSlipDetector::new();
    let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 14 };
    let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 18 };

    let o_ref_t0 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0, 0.0, 45.0, None);
    let o1_t0 = make_sat_obs_l1(14, Constellation::Gps, 20_000_000.0, 100_000.0, 200.0, 41.0, None);
    let o2_t0 = make_sat_obs_l1(18, Constellation::Gps, 20_000_000.0, 150_000.0, -150.0, 43.0, None);
    detector.check_epoch(&make_epoch(100.0, vec![o_ref_t0, o1_t0, o2_t0]));

    let o_ref_t1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0, 0.0, 45.0, None);
    let o1_t1 = make_sat_obs_l1(14, Constellation::Gps, 20_000_000.0, 100_000.0 - 200.0 + 2.0, 200.0, 38.0, Some(1));
    let o2_t1 = make_sat_obs_l1(18, Constellation::Gps, 20_000_000.0, 150_000.0 + 150.0, -150.0, 43.0, None);
    detector.check_epoch(&make_epoch(101.0, vec![o_ref_t1, o1_t1, o2_t1]));

    assert_eq!(detector.get_arc(sat1), 1, "Sat 1 must have slip recorded");
    assert_eq!(detector.get_arc(sat2), 0, "Sat 2 must maintain continuous arc without slip");

    let o_ref_t2 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0, 0.0, 45.0, None);
    let o1_t2 = make_sat_obs_l1(14, Constellation::Gps, 20_000_000.0, 100_000.0 - 400.0 + 2.0, 200.0, 40.0, None);
    let o2_t2 = make_sat_obs_l1(18, Constellation::Gps, 20_000_000.0, 150_000.0 + 300.0 + 0.5, -150.0, 35.0, Some(1));
    detector.check_epoch(&make_epoch(102.0, vec![o_ref_t2, o1_t2, o2_t2]));

    assert_eq!(detector.get_arc(sat2), 1, "Sat 2 must record slip on half-cycle jump");
}

/// Scenario 3: Highway Overpass Outage & Rapid Reacquisition.
#[test]
fn test_tier4_scenario3_highway_underpass_outage_reacquisition() {
    let mut detector = CycleSlipDetector::new();
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 21 };

    let o_before = make_sat_obs_l1(21, Constellation::Gps, 20_000_000.0, 80_000.0, 50.0, 44.0, None);
    detector.check_epoch(&make_epoch(500.0, vec![o_before]));
    assert_eq!(detector.get_arc(sat_id), 0);

    let o_after = make_sat_obs_l1(21, Constellation::Gps, 20_000_050.0, 80_150.0, 50.0, 28.0, None);
    detector.check_epoch(&make_epoch(503.5, vec![o_after]));

    assert_eq!(detector.get_arc(sat_id), 1, "Underpass outage gap must trigger arc reset");
}

/// Scenario 4: Asymmetric CORS Base vs Urban Canyon Rover Noise.
#[test]
fn test_tier4_scenario4_asymmetric_cors_base_vs_canyon_rover() {
    let rx = rover_ecef();
    let sat = sat_pos_az_el(rx, 45.0, 35.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);

    let snrs = (Some(22), Some(42), Some(48), Some(48));
    let var = compute_dd_variances(rx, sat, ref_sat, GPS_L1_WAVELENGTH_M, snrs);

    assert!(var.pr_var_m2 > var.pr_ref_var_m2 * 2.0, "Attenuated rover noise dominates total variance");
    assert!(var.pr_ref_var_m2 > 0.0, "Base reference variance is non-zero");
    assert!(var.pr_var_m2.is_finite() && !var.pr_var_m2.is_nan());
}

/// Scenario 5: Collinear Street Canyon Geometry & DOP Degeneracy Guard.
#[test]
fn test_tier4_scenario5_collinear_street_canyon_dop_guard() {
    let rx = rover_ecef();
    let collinear_sats = vec![
        sat_pos_az_el(rx, 0.0, 20.0, 20_000_000.0),
        sat_pos_az_el(rx, 0.0, 45.0, 20_000_000.0),
        sat_pos_az_el(rx, 180.0, 30.0, 20_000_000.0),
        sat_pos_az_el(rx, 180.0, 60.0, 20_000_000.0),
    ];

    let dop = compute_dop_from_positions(rx, &collinear_sats);
    let is_degenerate = match dop {
        None => true,
        Some(d) => d.pdop > 10.0 || d.hdop > 10.0,
    };

    assert!(is_degenerate, "Collinear street canyon satellites must be identified as degenerate");
}
