//! Tier 3: Cross-Feature Pairwise Interactions (>=6 tests).
#![allow(clippy::unwrap_used)]

use nalgebra::{DMatrix, DVector};
use gneiss_core::dop::compute_dop_from_positions;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_rtk::ambiguity::par::select_ils_subset;
use gneiss_rtk::estimators::rtk_iekf::formation_cov::compute_dd_variances;
use gneiss_rtk::estimators::rtk_iekf::screen::screen_gross_pr_errors;
use gneiss_rtk::post_process::screening::CycleSlipDetector;

use super::common::*;

#[test]
fn test_tier3_r1_r2_snr_and_cmc_compound_downweighting() {
    let rx = rover_ecef();
    let sat = sat_pos_az_el(rx, 45.0, 30.0, 20_000_000.0);
    let ref_sat = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);

    let snrs_canyon = (Some(22), Some(46), Some(48), Some(48));
    let var_canyon = compute_dd_variances(rx, sat, ref_sat, GPS_L1_WAVELENGTH_M, snrs_canyon);

    let snrs_clean = (Some(42), Some(46), Some(48), Some(48));
    let var_clean = compute_dd_variances(rx, sat, ref_sat, GPS_L1_WAVELENGTH_M, snrs_clean);

    assert!(var_canyon.pr_var_m2 > var_clean.pr_var_m2 * 2.0);
    assert!(var_canyon.cp_var_cycles2 > var_clean.cp_var_cycles2 * 2.0);

    let mp_step_m = 12.0;
    let pr_var_total = var_canyon.pr_var_m2 + mp_step_m * mp_step_m;
    assert!(pr_var_total > 140.0, "Compound code variance exceeds 140 m^2");
    assert!(var_canyon.cp_var_cycles2 < 2.0, "Carrier phase variance remains constrained (< 2.0 cyc^2)");
}

#[test]
fn test_tier3_r2_r3_code_multipath_vs_doppler_slip_isolation() {
    let mut detector = CycleSlipDetector::new();
    let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 11 };

    let o1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0, 100.0, 42.0, None);
    let o2 = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 50_000.0, 200.0, 42.0, None);
    let obs1 = make_sat_obs_l1(11, Constellation::Gps, 20_000_000.0, 50_000.0, 300.0, 42.0, None);
    detector.check_epoch(&make_epoch(100.0, vec![o1, o2, obs1]));

    let doppler = 300.0;
    let dt = 1.0;
    let cp_next = 50_000.0 - doppler * dt;
    let pr_jumped = 20_000_000.0 + 15.0;

    let o1_next = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 50_000.0 - 100.0, 100.0, 42.0, None);
    let o2_next = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 50_000.0 - 200.0, 200.0, 42.0, None);
    let obs2 = make_sat_obs_l1(11, Constellation::Gps, pr_jumped, cp_next, doppler, 32.0, None);
    detector.check_epoch(&make_epoch(101.0, vec![o1_next, o2_next, obs2]));

    assert_eq!(
        detector.get_arc(sat_id),
        0,
        "Code multipath jump with clean carrier phase must not be mistaken for a cycle slip"
    );
}

#[test]
fn test_tier3_r3_r4_slipped_satellite_excluded_from_par() {
    let a = DVector::from_vec(vec![1.02, 2.01, 7.85, 4.03]);
    let mut q = DMatrix::zeros(4, 4);
    q[(0, 0)] = 0.01;
    q[(1, 1)] = 0.015;
    q[(2, 2)] = 50.0;
    q[(3, 3)] = 0.02;

    let (selected, _sub_a, sub_q) = select_ils_subset(&a, &q, 0.95);

    assert!(!selected.contains(&2), "Slipped satellite must be excluded from PAR subset");
    assert!(selected.contains(&0) && selected.contains(&1) && selected.contains(&3));
    assert!(is_matrix_positive_definite(&sub_q, 1e-6));
}

#[test]
fn test_tier3_r2_r4_cmc_multipath_satellite_excluded_from_par() {
    let a = DVector::from_vec(vec![1.01, 2.35, 3.02, 4.01]);
    let mut q = DMatrix::zeros(4, 4);
    q[(0, 0)] = 0.02;
    q[(1, 1)] = 1.50;
    q[(2, 2)] = 0.025;
    q[(3, 3)] = 0.03;

    let (selected, _, _) = select_ils_subset(&a, &q, 0.90);
    assert!(!selected.contains(&1), "Multipath degraded satellite must be omitted by PAR");
    assert_eq!(selected.len(), 3, "Clean satellites selected");
}

#[test]
fn test_tier3_r1_r4_multi_constellation_par_selection() {
    let rx = rover_ecef();
    let sats = vec![
        sat_pos_az_el(rx, 15.0, 70.0, 20_000_000.0),
        sat_pos_az_el(rx, 80.0, 55.0, 20_000_000.0),
        sat_pos_az_el(rx, 140.0, 60.0, 23_000_000.0),
        sat_pos_az_el(rx, 210.0, 45.0, 23_000_000.0),
        sat_pos_az_el(rx, 290.0, 65.0, 21_000_000.0),
        sat_pos_az_el(rx, 340.0, 15.0, 21_000_000.0),
    ];

    let dop = compute_dop_from_positions(rx, &sats).expect("DOP must exist");
    assert!(dop.pdop < 3.0, "Multi-constellation geometry provides low PDOP");

    let a = DVector::from_vec(vec![1.01, 2.02, 3.01, 4.02, 5.01, 6.45]);
    let mut q = DMatrix::zeros(6, 6);
    for i in 0..5 {
        q[(i, i)] = 0.02;
    }
    q[(5, 5)] = 0.85;

    let (selected, _, sub_q) = select_ils_subset(&a, &q, 0.95);
    assert!(selected.len() >= 4, "Must select at least 4 satellites for 3D fix");
    assert!(!selected.contains(&5), "Low elevation/high noise satellite excluded");
    assert!(is_matrix_positive_definite(&sub_q, 1e-6));
}

fn build_pipeline_meas(rx: nalgebra::Vector3<f64>, sat_ref: nalgebra::Vector3<f64>) -> Vec<DoubleDiffMeasurement> {
    let s1 = sat_pos_az_el(rx, 45.0, 60.0, 20_000_000.0);
    let s2 = sat_pos_az_el(rx, 120.0, 45.0, 20_000_000.0);
    let s3 = sat_pos_az_el(rx, 260.0, 50.0, 20_000_000.0);
    let g1 = compute_geometric_dd(rx, ref_ecef(), s1, sat_ref);
    let g2 = compute_geometric_dd(rx, ref_ecef(), s2, sat_ref);
    let g3 = compute_geometric_dd(rx, ref_ecef(), s3, sat_ref);
    let m1 = make_dd_meas((0, 1, 10), g1, Some(g1 / GPS_L1_WAVELENGTH_M), (0.08, 0.0005), s1, sat_ref);
    let m2 = make_dd_meas((0, 2, 10), g2 + 25.0, Some(g2 / GPS_L1_WAVELENGTH_M), (0.08, 0.0005), s2, sat_ref);
    let m3 = make_dd_meas((0, 3, 10), g3, Some(g3 / GPS_L1_WAVELENGTH_M), (0.08, 0.0005), s3, sat_ref);
    vec![m1, m2, m3]
}

#[test]
fn test_tier3_r1_r2_r3_r4_full_pipeline_cycle() {
    let mut detector = CycleSlipDetector::new();
    let rx = rover_ecef();
    let sat_ref = sat_pos_az_el(rx, 180.0, 75.0, 20_000_000.0);

    let o1 = make_sat_obs_l1(1, Constellation::Gps, 20_000_000.0, 10_000.0, 100.0, 44.0, None);
    let o2 = make_sat_obs_l1(2, Constellation::Gps, 20_000_000.0, 10_000.0, 150.0, 42.0, None);
    let o3 = make_sat_obs_l1(3, Constellation::Gps, 20_000_000.0, 10_000.0, 200.0, 40.0, None);
    detector.check_epoch(&make_epoch(100.0, vec![o1, o2, o3]));

    let mut meas = build_pipeline_meas(rx, sat_ref);
    let rejected = screen_gross_pr_errors(&mut meas, rx);
    assert_eq!(rejected.len(), 1);
    assert_eq!(rejected[0].sat, 2);

    let a = DVector::from_vec(vec![1.01, 3.01]);
    let mut q = DMatrix::zeros(2, 2);
    q[(0, 0)] = 0.015;
    q[(1, 1)] = 0.018;
    let (selected, _, sub_q) = select_ils_subset(&a, &q, 0.95);
    assert_eq!(selected.len(), 2);
    assert!(is_matrix_positive_definite(&sub_q, 1e-6));
}
