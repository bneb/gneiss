#![allow(clippy::unwrap_used)]

use super::*;
use crate::sim::generator::{generate_simulation_dataset, SimulationConfig};

#[test]
fn test_gate_defaults_off_and_mask_defaults_to_measurement_mask() {
    let eng = test_engine(GpsTime::new(2200, 0.0));
    assert!(!eng.ar_gate, "gate must default off without env");
    assert_eq!(eng.ar_elevation_mask_rad, eng.min_elevation_rad);
    assert!(eng.code_phase_div.is_empty());
}

#[test]
fn test_ar_elevation_mask_unit_view_through_engine() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.ar_gate = true;
    eng.state.pos_ecef = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    let hi = dd_key(2, 1);
    let lo = dd_key(3, 1);
    eng.state.ensure_ambiguity(hi, 10.5, 4.0);
    eng.state.ensure_ambiguity(lo, 20.5, 9.0);
    let up = eng.state.pos_ecef.normalize();
    let east = Vector3::new(-up.z, 0.0, up.x).normalize();
    let mk_meas = |k: DoubleDiffKey, dir: Vector3<f64>| update::DoubleDiffMeasurement {
        key: k,
        dd_pr_m: 0.0,
        dd_cp_cycles: None,
        sat_pos: eng.state.pos_ecef + dir * 2.4e7,
        ref_pos: eng.state.pos_ecef + up * 2.6e7,
        base_pos: eng.state.pos_ecef,
        lambda: 0.19,
        pr_var_m2: 1.0,
        cp_var_cycles2: 1.0,
        pr_ref_var_m2: 0.5,
        cp_ref_var_cycles2: 0.5,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m: 0.0,
    };
    let meas = vec![mk_meas(hi, up), mk_meas(lo, east)];
    eng.ar_elevation_mask_rad = 0.2618;
    let view = ar_gate::elevation_filtered_view(&eng.state, &meas, eng.ar_elevation_mask_rad);
    assert_eq!(view.ambiguities.len(), 1);
    assert_eq!(view.ambiguities[0].0, hi);
    assert_eq!(eng.state.ambiguities.len(), 2);
}

#[test]
fn test_gate_off_extreme_ar_mask_is_ignored_end_to_end() {
    let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    eng.ar_elevation_mask_rad = 1.55;
    let res = run_sim(&mut eng, &sim, cfg.base_ecef);
    let fixed = res.iter().filter(|(f, _)| *f).count();
    assert!(fixed >= 5, "legacy path ignores the AR mask, got {fixed} fixes");
}

#[test]
fn test_gate_on_extreme_ar_mask_suppresses_all_fixes() {
    let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    eng.ar_gate = true;
    eng.ar_elevation_mask_rad = 1.55;
    let res = run_sim(&mut eng, &sim, cfg.base_ecef);
    assert!(res.iter().all(|(f, _)| !f), "AR mask must exclude all pairs; unexpected fix");
}

#[test]
fn test_gate_on_with_default_masks_still_fixes_sim() {
    let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    eng.ar_gate = true;
    let mut fixed = 0;
    for i in 0..sim.rover_epochs.len() {
        let s = eng.process_epoch(
            &sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides,
        ).expect("epoch must process");
        if s.is_fixed {
            fixed += 1;
            let err = (s.position_ecef - sim.truth_positions[i].1).norm();
            assert!(err < 0.05, "gated fixed error {err:.4} m exceeds 5 cm");
        }
    }
    assert!(fixed >= 5, "gate on with defaults should still fix, got {fixed}");
}

fn mean_tail_error(
    engine: &mut GnssRtkIekf,
    sim: &crate::sim::generator::SimulationDataset,
    base: Vector3<f64>,
    skip: usize,
) -> f64 {
    let mut sum = 0.0;
    let mut n = 0;
    for i in 0..sim.rover_epochs.len() {
        let sol = engine
            .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], base, &sim.ephemerides)
            .expect("epoch must process");
        if i >= skip {
            sum += (sol.position_ecef - sim.truth_positions[i].1).norm();
            n += 1;
        }
    }
    sum / n.max(1) as f64
}

#[test]
fn kinematic_profile_tracks_linear_ramp_and_velocity_states_stay_live() {
    use crate::post_process::dynamics::{KINEMATIC_INNOV_GATE_SCALE, KINEMATIC_Q_ACCEL, STATIC_Q_ACCEL};
    use crate::sim::generator::TrajectoryProfile;
    let cfg = SimulationConfig {
        duration_s: 1800.0,
        epoch_rate_hz: 1.0 / 30.0,
        profile: TrajectoryProfile::Linear {
            start_offset_ned: Vector3::new(100.0, 100.0, 0.0),
            velocity_ned: Vector3::new(4.0, 3.0, 0.0),
        },
        ..Default::default()
    };
    let sim = generate_simulation_dataset(&cfg);
    let true_speed = Vector3::new(4.0_f64, 3.0, 0.0).norm();
    let mut static_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, STATIC_Q_ACCEL);
    let mut kin_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, KINEMATIC_Q_ACCEL);
    kin_eng.robust_innov_scale = KINEMATIC_INNOV_GATE_SCALE;

    let err_static = mean_tail_error(&mut static_eng, &sim, cfg.base_ecef, 10);
    let err_kin = mean_tail_error(&mut kin_eng, &sim, cfg.base_ecef, 10);
    let speed_kin = kin_eng.state.vel_ecef.norm();
    assert!(err_kin < 3.0, "kinematic must track the ramp, got {err_kin:.3} m");
    assert!(err_static < 3.0, "static+Huber also tracks steady ramps, got {err_static:.3} m");
    assert!(
        (speed_kin - true_speed).abs() < 0.25 * true_speed,
        "kinematic velocity state must track true speed, got {speed_kin:.2}"
    );
}

#[test]
fn kinematic_profile_outperforms_static_under_acceleration() {
    use crate::post_process::dynamics::{KINEMATIC_INNOV_GATE_SCALE, KINEMATIC_Q_ACCEL, STATIC_Q_ACCEL};
    use crate::sim::generator::TrajectoryProfile;
    let cfg = SimulationConfig {
        duration_s: 3600.0,
        epoch_rate_hz: 1.0 / 30.0,
        profile: TrajectoryProfile::Circular {
            center_offset_ned: Vector3::new(200.0, 200.0, 0.0),
            radius_m: 500.0,
            speed_m_s: 15.0,
        },
        ..Default::default()
    };
    let sim = generate_simulation_dataset(&cfg);
    let mut static_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, STATIC_Q_ACCEL);
    let mut kin_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, KINEMATIC_Q_ACCEL);
    kin_eng.robust_innov_scale = KINEMATIC_INNOV_GATE_SCALE;

    let err_static = mean_tail_error(&mut static_eng, &sim, cfg.base_ecef, 12);
    let err_kin = mean_tail_error(&mut kin_eng, &sim, cfg.base_ecef, 12);
    assert!(err_kin < 3.0, "kinematic must track the turn, got {err_kin:.3} m");
    assert!(
        err_static > 2.0 * err_kin,
        "static-tuned Q must lag under sustained acceleration (static={err_static:.3}, kin={err_kin:.3})"
    );
}

#[test]
fn test_new_bias_init_gate_off_keeps_raw_code_phase_seed() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    let lambda = 0.2_f64;
    let raw = 105.75 - 100.0 / lambda;
    eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false, 1.0);
    assert_eq!(eng.state.get_amb_idx(&dd_key(5, 1)), Some(6));
    assert!((eng.state.ambiguities[0].1 - raw).abs() < 1e-12, "seed must equal raw init");
    assert!(eng.code_phase_div.is_empty(), "no samples recorded when gate off");
}

#[test]
fn test_new_bias_init_applies_median_coherency_offset() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.ar_gate = true;
    eng.state.ensure_ambiguity(dd_key(2, 1), 50.25, 100.0);
    eng.code_phase_div = vec![(1, 2.5), (2, 100.0), (1, 3.5), (1, 4.5)];
    let lambda = 0.2_f64;
    let raw = 105.75 - 100.0 / lambda;
    eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false, 1.0);
    let seeded = eng.state.ambiguities.iter().find(|(k, _)| *k == dd_key(5, 1)).unwrap().1;
    assert!((seeded - (raw - 3.5)).abs() < 1e-12, "seed = raw − median(band-1 divergences)");
    assert_eq!(eng.code_phase_div.len(), 5, "seeded pair joins the epoch's sample set");
    assert_eq!(eng.code_phase_div.last().copied(), Some((1, 3.5)));
}

#[test]
fn test_new_bias_init_without_prior_samples_seeds_raw_under_gate() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.ar_gate = true;
    let lambda = 0.2_f64;
    let raw = 105.75 - 100.0 / lambda;
    eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false, 1.0);
    let seeded = eng.state.ambiguities[0].1;
    assert!((seeded - raw).abs() < 1e-12, "no prior pairs -> no offset");
    assert_eq!(eng.code_phase_div, vec![(1, 0.0)]);
}

#[test]
fn test_coherency_offset_never_crosses_frequency_bands() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.ar_gate = true;
    eng.code_phase_div = vec![(2, 2.5), (2, 3.5)];
    let lambda = 0.2_f64;
    let raw = 105.75 - 100.0 / lambda;
    eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false, 1.0);
    assert!((eng.state.ambiguities[0].1 - raw).abs() < 1e-12);
}

#[test]
fn test_slip_reset_stays_raw_even_under_gate() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.ar_gate = true;
    eng.state.ensure_ambiguity(dd_key(2, 1), 50.25, 100.0);
    eng.code_phase_div = vec![(1, 2.5), (1, 3.5)];
    let lambda = 0.2_f64;
    let raw = 205.75 - 100.0 / lambda;
    eng.update_dd_ambiguity(dd_key(2, 1), Some(205.75), 100.0, lambda, true, 1.0);
    assert!((eng.state.ambiguities[0].1 - raw).abs() < 1e-12);
    assert_eq!(eng.code_phase_div.len(), 3, "established pair still contributes a sample");
}

#[test]
fn test_gate_disabled_run_matches_default_run_bit_for_bit() {
    let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
    let sim = generate_simulation_dataset(&cfg);
    let mut baseline = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    let mut gated = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    gated.ar_gate = false;
    let base_out = run_sim(&mut baseline, &sim, cfg.base_ecef);
    let gate_out = run_sim(&mut gated, &sim, cfg.base_ecef);
    assert_eq!(base_out, gate_out, "gate off must reproduce legacy exactly");
    assert!(base_out.iter().filter(|(f, _)| *f).count() >= 5);
}
