#![allow(clippy::unwrap_used)]

use std::collections::HashMap;
use nalgebra::Vector3;
use gneiss_core::time::GpsTime;
use crate::estimators::rtk_iekf::state::{DoubleDiffKey, RtkState};
use super::system::build_measurement_system;
use super::robust::{if_residual_outliers, robust_inflate, ZWD_INIT_VAR_M2, ZWD_RW_M2_PER_S};
use super::{iekf_update, iekf_update_gated, pcv_corrected_cp, DoubleDiffMeasurement};

#[test]
fn robust_gate_scale_one_is_bitwise_legacy() {
    let var = 0.04;
    for scale in [1.0_f64, 3.0] {
        assert_eq!(robust_inflate(0.1, var, scale).to_bits(), var.to_bits());
    }
    assert_eq!(robust_inflate(3.0_f64, 1.0_f64, 1.0).to_bits(), 1.0_f64.to_bits());
}

#[test]
fn kinematic_gate_keeps_nominal_weight_between_legacy_and_widened_knee() {
    let var = 0.04_f64;
    let innov = (18.0_f64 * var).sqrt();
    let inflated_static = robust_inflate(innov, var, 1.0);
    assert!(inflated_static > var, "static profile must inflate at NIS=18");
    assert_eq!(
        robust_inflate(innov, var, 3.0).to_bits(),
        var.to_bits(),
        "kinematic profile keeps nominal weight below its widened knee"
    );
    let gross = (1000.0 * var).sqrt();
    assert!(robust_inflate(gross, var, 3.0) < robust_inflate(gross, var, 1.0));
    assert!(robust_inflate(gross, var, 3.0) > var);
}

#[test]
fn gated_update_matches_legacy_update_with_unit_scale() {
    let build = || {
        let mut s = RtkState::new(Vector3::new(101.0, 198.0, 301.0), GpsTime::new(2000, 100.0));
        s.ensure_ambiguity(DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 }, 0.0, 100.0);
        s
    };
    let meas = vec![DoubleDiffMeasurement {
        key: DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 },
        dd_pr_m: 1.5,
        dd_cp_cycles: Some(8.0),
        sat_pos: Vector3::new(10_000.0, 20_000.0, 20_000.0),
        ref_pos: Vector3::new(5_000.0, 25_000.0, 20_000.0),
        base_pos: Vector3::zeros(),
        lambda: 0.190,
        pr_var_m2: 0.04,
        cp_var_cycles2: 0.0001,
        pr_ref_var_m2: 0.02,
        cp_ref_var_cycles2: 0.00005,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m: 0.0,
    }];
    let mut a = build();
    let mut b = build();
    iekf_update(&mut a, &meas).expect("legacy update ok");
    iekf_update_gated(&mut b, &meas, 1.0).expect("gated update ok");
    assert_eq!(a.pos_ecef, b.pos_ecef);
    for i in 0..a.cov.nrows() {
        for j in 0..a.cov.ncols() {
            assert_eq!(a.cov[(i, j)].to_bits(), b.cov[(i, j)].to_bits());
        }
    }
}

#[test]
fn test_iekf_update_reduces_position_error() {
    let true_pos = Vector3::new(100.0, 200.0, 300.0);
    let base_pos = Vector3::new(0.0, 0.0, 0.0);
    let time = GpsTime::new(2000, 100.0);

    let mut state = RtkState::new(true_pos + Vector3::new(2.0, -2.0, 1.0), time);
    let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
    state.ensure_ambiguity(key, 0.0, 100.0);

    let sat_pos = Vector3::new(10_000.0, 20_000.0, 20_000.0);
    let ref_pos = Vector3::new(5_000.0, 25_000.0, 20_000.0);
    let base_dd = (sat_pos - base_pos).norm() - (ref_pos - base_pos).norm();
    let true_dd = (sat_pos - true_pos).norm() - (ref_pos - true_pos).norm() - base_dd;

    let meas = vec![DoubleDiffMeasurement {
        key,
        dd_pr_m: true_dd,
        dd_cp_cycles: Some(true_dd / 0.190),
        sat_pos,
        ref_pos,
        base_pos,
        lambda: 0.190,
        pr_var_m2: 0.04,
        cp_var_cycles2: 0.0001,
        pr_ref_var_m2: 0.02,
        cp_ref_var_cycles2: 0.00005,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dm_wet_rov: 0.0,
        dd_pcv_m: 0.0,
    }];

    let res = iekf_update(&mut state, &meas);
    assert!(res.is_ok());
    assert!((state.pos_ecef.x - true_pos.x).abs() < 2.0);
}

#[test]
fn test_zwd_state_tracks_wet_delay_ramp() {
    use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let f1 = 1575.42e6_f64;
    let l1 = SPEED_OF_LIGHT_M_S / f1;
    let true_pos = Vector3::new(100.0, 200.0, 300.0);
    let base_pos = Vector3::new(0.0, 0.0, 0.0);
    let time = GpsTime::new(2000, 100.0);

    let mut state = RtkState::new(true_pos + Vector3::new(0.3, -0.2, 0.1), time);
    state.enable_zwd(ZWD_INIT_VAR_M2);

    let dirs = [
        (Vector3::new(20_000_000.0, 5_000_000.0, 8_000_000.0), 0.9),
        (Vector3::new(-6_000_000.0, 18_000_000.0, 19_000_000.0), 0.55),
        (Vector3::new(9_000_000.0, -14_000_000.0, 21_000_000.0), 0.35),
        (Vector3::new(2_000_000.0, 7_000_000.0, 22_500_000.0), 0.15),
    ];
    let ref_pos = Vector3::new(5_000.0, 25_000.0, 20_000.0);
    let mut meas = Vec::new();
    for (p, (dir, dm)) in dirs.iter().enumerate() {
        let key = DoubleDiffKey { constellation_id: 0, sat: 2 + p as u16, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(key, 50.0 + p as f64, 100.0);
        meas.push(DoubleDiffMeasurement {
            key,
            dd_pr_m: 0.0,
            dd_cp_cycles: Some(0.0),
            sat_pos: true_pos + *dir * 20.0,
            ref_pos,
            base_pos,
            lambda: l1,
            pr_var_m2: 0.04,
            cp_var_cycles2: 1e-4,
            pr_ref_var_m2: 0.02,
            cp_ref_var_cycles2: 0.5e-4,
            dm_wet_rov: *dm,
            dgrad_n_rov: 0.0,
            dgrad_e_rov: 0.0,
            tide_dd_m: 0.0,
            dd_pcv_m: 0.0,
        });
    }

    for k in 0..40 {
        let zwd_true = 0.0025 * k as f64;
        for (i, m) in meas.iter_mut().enumerate() {
            let pos = true_pos;
            let geom = (m.sat_pos - pos).norm() - (ref_pos - pos).norm()
                - (m.sat_pos - base_pos).norm() + (ref_pos - base_pos).norm();
            m.dd_pr_m = geom + m.dm_wet_rov * zwd_true;
            m.dd_cp_cycles = Some(geom / l1 + 50.0 + i as f64 + m.dm_wet_rov * zwd_true / l1);
        }
        let zi = state.zwd_idx().unwrap();
        state.cov[(zi, zi)] += ZWD_RW_M2_PER_S * 30.0;
        iekf_update(&mut state, &meas).expect("update ok");
    }

    assert!((state.zwd_m - 0.0975).abs() < 0.04,
        "zwd must track the wet ramp, got {:.4} vs {}", state.zwd_m, 0.0975);
    let pos_err = (state.pos_ecef - true_pos).norm();
    assert!(pos_err < 0.2, "position must stay pinned while zwd absorbs drift: {:.3}", pos_err);
}

#[allow(clippy::type_complexity)]
fn if_meas_fixture(
    slip_cycles: Option<f64>,
) -> (Vec<DoubleDiffMeasurement>, Vector3<f64>, HashMap<DoubleDiffKey, f64>, HashMap<DoubleDiffKey, f64>) {
    let true_pos = Vector3::new(0.0, 0.0, 0.0);
    let base_pos = Vector3::new(4_000_000.0, 500_000.0, 1_000_000.0);
    let dirs = [
        Vector3::new(-20_000.0, 15_000.0, 18_000.0),
        Vector3::new(5_000.0, 24_000.0, -14_000.0),
        Vector3::new(16_000.0, -18_000.0, 22_000.0),
        Vector3::new(-9_000.0, -6_000.0, -26_000.0),
        Vector3::new(22_000.0, 9_000.0, -19_000.0),
        Vector3::new(-15_000.0, 21_000.0, 8_000.0),
    ];
    let ref_dir = dirs[0];
    let f1 = 1575.42e6_f64;
    let f2 = 1227.60e6_f64;
    let l1 = 299792458.0 / f1;
    let l2 = 299792458.0 / f2;
    let mut meas = Vec::new();
    for (i, dir) in dirs.iter().enumerate().skip(1) {
        let key = DoubleDiffKey { constellation_id: 0, sat: 1 + i as u16, ref_sat: 1, freq_band: 1 };
        let key2 = DoubleDiffKey { freq_band: 2, ..key };
        let sat_p = *dir * 20_000_000.0;
        let ref_p = ref_dir * 20_000_000.0;
        let rho = |p: Vector3<f64>| (true_pos - p).norm();
        let geom_s = rho(sat_p); let geom_r = rho(ref_p);
        let base_dd = (base_pos - sat_p).norm() - (base_pos - ref_p).norm();
        for (k, lam, slip_here) in [(key, l1, false), (key2, l2, true)] {
            let n = 10.0 + i as f64;
            let mut cp = (geom_s - geom_r) / lam - base_dd / lam + n;
            if slip_here && i == 1 {
                if let Some(sl) = slip_cycles { cp += sl; }
            }
            meas.push(DoubleDiffMeasurement {
                key: k,
                dd_pr_m: (geom_s - geom_r) - base_dd,
                dd_cp_cycles: Some(cp),
                sat_pos: sat_p,
                ref_pos: ref_p,
                base_pos,
                lambda: lam,
                pr_var_m2: 0.04,
                cp_var_cycles2: 1e-4,
                pr_ref_var_m2: 0.02,
                cp_ref_var_cycles2: 0.5e-4,
                dm_wet_rov: 0.0,
                dgrad_n_rov: 0.0,
                dgrad_e_rov: 0.0,
                tide_dd_m: 0.0,
                dd_pcv_m: 0.0,
            });
        }
    }
    let mut n1 = HashMap::new();
    let mut n2 = HashMap::new();
    for (i, _) in dirs.iter().enumerate().skip(1) {
        let key = DoubleDiffKey { constellation_id: 0, sat: 1 + i as u16, ref_sat: 1, freq_band: 1 };
        let key2 = DoubleDiffKey { freq_band: 2, ..key };
        let n = 10.0 + i as f64;
        n1.insert(key, n);
        n2.insert(key2, n);
    }
    (meas, true_pos, n1, n2)
}

#[test]
fn test_if_screen_clean_data_has_no_outliers() {
    let (meas, pos, n1, n2) = if_meas_fixture(None);
    let out = if_residual_outliers(pos, &meas, &n1, &n2);
    assert!(out.is_empty(), "clean data flagged: {:?}", out);
}

#[test]
fn test_if_screen_flags_single_slipped_pair() {
    let (meas, pos, n1, n2) = if_meas_fixture(Some(2.0));
    let out = if_residual_outliers(pos, &meas, &n1, &n2);
    assert_eq!(out.len(), 1, "expected exactly the slipped pair: {:?}", out);
    assert_eq!(out[0].sat, 2, "slipped pair is sat=2 vs ref=1");
}

#[test]
fn test_if_screen_common_mode_position_error_not_flagged() {
    let (meas, _, n1, n2) = if_meas_fixture(None);
    let biased = Vector3::new(0.02, -0.012, 0.008);
    let out = if_residual_outliers(biased, &meas, &n1, &n2);
    assert!(out.is_empty(), "common-mode error flagged pairs: {:?}", out);
}

fn pcv_fixture(dd_pcv_m: f64) -> (RtkState, DoubleDiffMeasurement) {
    let truth = Vector3::new(100.0, 200.0, 300.0);
    let mut state = RtkState::new(truth, GpsTime::new(2000, 100.0));
    let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
    state.ensure_ambiguity(key, 7.25, 100.0);
    let m = DoubleDiffMeasurement {
        key,
        dd_pr_m: 10.0,
        dd_cp_cycles: Some(55.0),
        sat_pos: truth + Vector3::new(1.2e7, 0.4e7, 1.8e7),
        ref_pos: truth + Vector3::new(-0.6e7, 1.9e7, 1.1e7),
        base_pos: truth,
        lambda: 0.190,
        pr_var_m2: 0.04,
        cp_var_cycles2: 1e-4,
        pr_ref_var_m2: 0.02,
        cp_ref_var_cycles2: 0.5e-4,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m,
    };
    (state, m)
}

#[test]
fn pcv_corrected_cp_subtracts_pcv_over_lambda() {
    let (_, m) = pcv_fixture(0.008);
    let corrected = pcv_corrected_cp(&m).expect("phase present");
    assert!(
        (corrected - (55.0 - 0.008 / 0.190)).abs() < 1e-12,
        "corrected={corrected}"
    );
}

#[test]
fn pcv_corrected_cp_zero_pcv_is_identity() {
    let (_, m) = pcv_fixture(0.0);
    assert_eq!(pcv_corrected_cp(&m), Some(55.0));
}

#[test]
fn pcv_corrected_cp_passes_none_through() {
    let (_, mut m) = pcv_fixture(0.008);
    m.dd_cp_cycles = None;
    assert_eq!(pcv_corrected_cp(&m), None);
}

#[test]
fn phase_innovation_shifts_exactly_by_dd_pcv_over_lambda() {
    let pcv_m = 0.008_f64;
    let (state_off, m_off) = pcv_fixture(0.0);
    let (state_on, m_on) = pcv_fixture(pcv_m);
    let x = |s: &RtkState| s.to_dvector();
    let (_, y_off, _) = build_measurement_system(&state_off, &x(&state_off), &[m_off], 1.0);
    let (_, y_on, _) = build_measurement_system(&state_on, &x(&state_on), &[m_on], 1.0);
    assert_eq!(y_off.len(), 2, "code + phase rows expected");
    assert_eq!(y_on.len(), 2);
    let shift = y_off[1] - y_on[1];
    let expected = pcv_m / 0.190;
    assert!(
        (shift - expected).abs() < 1e-12,
        "shift={shift} expected={expected}"
    );
}

#[test]
fn test_correlated_dd_covariance_matrix() {
    let truth = Vector3::new(100.0, 200.0, 300.0);
    let mut state = RtkState::new(truth, GpsTime::new(2000, 100.0));
    let k2 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
    let k3 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };
    state.ensure_ambiguity(k2, 5.0, 100.0);
    state.ensure_ambiguity(k3, 8.0, 100.0);

    let mk_meas = |key: DoubleDiffKey, dx: f64| DoubleDiffMeasurement {
        key,
        dd_pr_m: 10.0,
        dd_cp_cycles: Some(50.0),
        sat_pos: truth + Vector3::new(dx, 2e7, 1e7),
        ref_pos: truth + Vector3::new(0.0, 2.5e7, 1e7),
        base_pos: truth,
        lambda: 0.190,
        pr_var_m2: 0.05,
        cp_var_cycles2: 0.0002,
        pr_ref_var_m2: 0.02,
        cp_ref_var_cycles2: 0.00008,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m: 0.0,
    };
    let meas = vec![mk_meas(k2, 1e6), mk_meas(k3, 2e6)];
    let (_, _, r) = build_measurement_system(&state, &state.to_dvector(), &meas, 1.0);

    assert_eq!(r.nrows(), 4);
    assert!((r[(0, 2)] - 0.02).abs() < 1e-12, "PR cross-covariance");
    assert!((r[(2, 0)] - 0.02).abs() < 1e-12, "PR symmetry");
    assert!((r[(1, 3)] - 0.00008).abs() < 1e-12, "CP cross-covariance");
    assert!((r[(3, 1)] - 0.00008).abs() < 1e-12, "CP symmetry");
    assert_eq!(r[(0, 1)], 0.0, "PR-CP cross is zero");
    assert_eq!(r[(0, 3)], 0.0, "PR-CP cross is zero");
    assert!(r.cholesky().is_some(), "R must be positive definite");
}
