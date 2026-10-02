//! Sprint-hunt tests for the DD-IEKF engine (`rtk_iekf/mod.rs`).
//!
//! Scope: `process_epoch`, `resolve_ar_candidate`, `screen_fixed_residuals`,
//! `finalize_fixed_position`, `update_dd_ambiguity`, `get_or_init_iekf`
//! wiring, and the invariants that must hold across a real run.
#![allow(clippy::unwrap_used)]

use nalgebra::{Vector3};

use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

use super::state::DoubleDiffKey;
use super::update::{iekf_update_gated, DoubleDiffMeasurement};
use super::{GnssRtkIekf, RtkState};
use crate::sim::generator::{generate_simulation_dataset, SimulationConfig, TrajectoryProfile};

const BASE_ECEF: Vector3<f64> = Vector3::new(-3_961_904.434_1, 3_348_994.266, 3_698_211.706_7);

fn dd_key(sat: u16, band: u8) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: band }
}

fn meas(sat: u16, dd_pr_m: f64, dd_cp_cycles: Option<f64>, sat_pos: Vector3<f64>, ref_pos: Vector3<f64>) -> DoubleDiffMeasurement {
    DoubleDiffMeasurement {
        key: dd_key(sat, 1),
        dd_pr_m,
        dd_cp_cycles,
        sat_pos,
        ref_pos,
        base_pos: BASE_ECEF,
        lambda: 0.19029367,
        pr_var_m2: 0.04,
        cp_var_cycles2: 1e-4,
        pr_ref_var_m2: 0.02,
        cp_ref_var_cycles2: 5e-5,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m: 0.0,
    }
}

/// Six satellites spread around the sky, referenced to satellite 1.
fn good_measurements() -> Vec<DoubleDiffMeasurement> {
    let rx = BASE_ECEF;
    let sat_pos = |i: u16| -> Vector3<f64> {
        let az = 2.0 * std::f64::consts::PI * (f64::from(i) - 1.0) / 6.0;
        rx + Vector3::new(az.cos(), az.sin(), 0.6).normalize() * 2.1e7
    };
    let ref_pos = sat_pos(1);
    // dd_pr must equal the geometric double difference at `rx` so the float
    // solution starts with zero innovation:
    //   dd_pr = (|p_i - rx| - |p_ref - rx|) - (|p_i - base| - |p_ref - base|)
    (2..=6u16)
        .map(|i| {
            let p = sat_pos(i);
            let geom = (p - rx).norm() - (ref_pos - rx).norm()
                - ((p - BASE_ECEF).norm() - (ref_pos - BASE_ECEF).norm());
            meas(i, geom, Some(1_000_000.0 + 1000.0 * f64::from(i)), p, ref_pos)
        })
        .collect()
}

// ===========================================================================
// NaN containment at the filter's own update path
// ===========================================================================

/// `iekf_update_gated` must never leave the state carrying a non-finite value.
/// Its inverse comes from `s.cholesky().map(inverse).or_else(try_inverse)`;
/// neither is guaranteed to reject a NaN, and `screen_gross_pr_errors` runs
/// OUTSIDE this function. So the contract has to hold here.
#[test]
#[ignore = "BUG-4 (update/kalman.rs:11-28): neither cholesky() nor try_inverse() rejects a\nNaN, so compute_iekf_step returns Some(NaN) and the filter state is poisoned\npermanently. Re-enable once iekf_update_gated rejects non-finite measurements."]
fn iekf_update_with_a_nan_pseudorange_never_nans_the_state() {
    let mut state = RtkState::new(BASE_ECEF, GpsTime::new(2200, 300_000.0));
    let mut ms = good_measurements();
    ms[2].dd_pr_m = f64::NAN;
    let res = iekf_update_gated(&mut state, &ms, 1.0);
    assert!(
        state.pos_ecef.iter().all(|v| v.is_finite()),
        "a NaN pseudorange produced a non-finite position: {:?}",
        state.pos_ecef
    );
    assert!(
        state.cov.iter().all(|v| v.is_finite()),
        "a NaN pseudorange produced a non-finite covariance"
    );
    // Whether the update is rejected or down-weighted, both outcomes are
    // acceptable; silently returning Ok while the state is NaN is not.
    assert!(
        res.is_err() || state.pos_ecef.iter().all(|v| v.is_finite()),
        "iekf_update_gated returned {res:?} with a NaN position"
    );
}

/// The same non-finite observation, injected through the PUBLIC entry point.
/// `process_epoch` runs `screen_gross_pr_errors` first, whose `worst_pr_residual`
/// ranks candidates with `f64::total_cmp` (NaN sorts last, so NaN is picked as
/// "worst" and then de-weighted) — an undocumented accident, not a guard. This
/// test pins the end-to-end contract.
#[test]
fn process_epoch_recovers_from_a_single_non_finite_observation() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    for i in 0..3 {
        eng.process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
            .unwrap_or_else(|e| panic!("warm-up epoch {i} failed: {e}"));
    }
    for bad in [f64::NAN, f64::INFINITY] {
        let mut eng2 = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        for i in 0..3 {
            eng2
                .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
                .expect("warm-up must process");
        }
        let mut poisoned = sim.rover_epochs[3].clone();
        let s0 = &mut poisoned.satellites[0];
        for o in s0.observations.iter_mut() {
            if o.code.obs_type == ObsType::Pseudorange {
                o.value = bad;
                break;
            }
        }
        let _ = eng2.process_epoch(&poisoned, &sim.base_epochs[3], cfg.base_ecef, &sim.ephemerides);
        for i in 4..8 {
            let sol = eng2
                .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
                .unwrap_or_else(|e| panic!("clean epoch {i} after a {bad} input must still process: {e}"));
            assert!(
                sol.position_ecef.iter().all(|v| v.is_finite()),
                "a {bad} pseudorange permanently poisoned the filter -> {:?}",
                sol.position_ecef
            );
            let err = (sol.position_ecef - sim.truth_positions[i].1).norm();
            assert!(err < 2.0, "epoch {i}: {err:.3} m error after a {bad} input");
        }
    }
}

#[test]
#[ignore = "BUG-4 (update/kalman.rs:11-28): same defect for an infinite pseudorange.\nRe-enable with BUG-4."]
fn iekf_update_with_an_infinite_pseudorange_never_nans_the_state() {
    let mut state = RtkState::new(BASE_ECEF, GpsTime::new(2200, 300_000.0));
    let mut ms = good_measurements();
    ms[1].dd_pr_m = f64::INFINITY;
    let _ = iekf_update_gated(&mut state, &ms, 1.0);
    assert!(state.pos_ecef.iter().all(|v| v.is_finite()), "{:?}", state.pos_ecef);
    assert!(state.cov.iter().all(|v| v.is_finite()));
}

#[test]
fn iekf_update_with_a_nan_carrier_phase_never_nans_the_state() {
    let mut state = RtkState::new(BASE_ECEF, GpsTime::new(2200, 300_000.0));
    let mut ms = good_measurements();
    ms[0].dd_cp_cycles = Some(f64::NAN);
    let _ = iekf_update_gated(&mut state, &ms, 1.0);
    assert!(state.pos_ecef.iter().all(|v| v.is_finite()), "{:?}", state.pos_ecef);
    assert!(state.cov.iter().all(|v| v.is_finite()));
}

/// A clean update must actually MOVE the state. This is the control that stops
/// the NaN tests above from passing merely because the update is a no-op.
///
/// The start offset is deliberately small (0.3 m): `append_dd_code_row` drops
/// a row outright via `is_code_blunder` once the code innovation is large, so
/// a large offset is (correctly) rejected rather than converged from.
#[test]
fn a_clean_dd_update_converges_toward_the_geometry() {
    let start = BASE_ECEF + Vector3::new(0.3, 0.0, 0.0);
    let mut state = RtkState::new(start, GpsTime::new(2200, 300_000.0));
    let ms = good_measurements();
    let before = (state.pos_ecef - BASE_ECEF).norm();
    for _ in 0..5 {
        iekf_update_gated(&mut state, &ms, 1.0).expect("clean update must succeed");
    }
    let after = (state.pos_ecef - BASE_ECEF).norm();
    assert!(before > 0.2, "fixture must start off, got {before}");
    assert!(after < 0.05, "5 consistent DD updates must converge; error still {after} m");
}

/// The covariance must stay symmetric positive-definite through updates.
#[test]
fn dd_update_covariance_stays_symmetric_and_positive_definite() {
    let mut state = RtkState::new(BASE_ECEF + Vector3::new(10.0, 5.0, -3.0), GpsTime::new(2200, 300_000.0));
    let ms = good_measurements();
    for _ in 0..20 {
        iekf_update_gated(&mut state, &ms, 1.0).expect("update must succeed");
        let asym = (&state.cov - state.cov.transpose()).norm();
        assert!(asym < 1e-6, "covariance asymmetry {asym}");
        let eig = state.cov.clone().symmetric_eigenvalues();
        assert!(eig[0] > 0.0, "covariance lost positive definiteness: min eig {}", eig[0]);
    }
}

// ===========================================================================
// Full engine invariants over a simulated run
// ===========================================================================

fn sim_config() -> SimulationConfig {
    SimulationConfig {
        duration_s: 25.0,
        num_satellites: 20,
        profile: TrajectoryProfile::Static { offset_ned: Vector3::new(25.0, 15.0, 0.0) },
        ..Default::default()
    }
}

#[test]
fn engine_run_keeps_state_finite_symmetric_and_dimension_consistent() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    let mut n_ok = 0;
    for i in 0..sim.rover_epochs.len() {
        let sol = eng
            .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
            .unwrap_or_else(|e| panic!("epoch {i} failed: {e}"));
        n_ok += 1;
        assert!(sol.position_ecef.iter().all(|v| v.is_finite()), "epoch {i}");
        assert!(sol.cov_position.iter().all(|v| v.is_finite()), "epoch {i}");
        let cov = &eng.state.cov;
        assert!((cov - cov.transpose()).norm() < 1e-6, "epoch {i} covariance asymmetry");
        assert!(cov.diagonal().iter().all(|v| *v > 0.0), "epoch {i} lost a positive variance");
        // The recorded snapshot must match the live state dimension, or the
        // backward pass silently mixes incompatible spaces.
        let last = eng.history.last().expect("history must record every epoch");
        assert_eq!(last.x_post.len(), eng.state.dim(), "epoch {i} snapshot dimension drift");
        assert_eq!(last.p_post.nrows(), eng.state.dim(), "epoch {i} snapshot dimension drift");
        assert_eq!(last.f_mat.nrows(), eng.state.dim(), "epoch {i} transition dimension drift");
    }
    assert!(n_ok >= 20, "expected every epoch to process, got {n_ok}");
    let smoothed = eng.smooth();
    assert_eq!(smoothed.len(), n_ok);
    for s in &smoothed {
        assert!(s.position_ecef.iter().all(|v| v.is_finite()));
    }
}

#[test]
fn engine_epoch_with_no_common_satellites_is_rejected_not_silently_accepted() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    eng.process_epoch(&sim.rover_epochs[0], &sim.base_epochs[0], cfg.base_ecef, &sim.ephemerides)
        .expect("first epoch must process");
    // A base epoch with nothing in common with the rover forms no DD pair.
    let empty_base = EpochObs {
        time: sim.base_epochs[1].time,
        satellites: vec![SatObs {
            sat: SatelliteId { constellation: Constellation::Galileo, prn: 99 },
            observations: Vec::new(),
        }],
    };
    let err = eng
        .process_epoch(&sim.rover_epochs[1], &empty_base, cfg.base_ecef, &sim.ephemerides)
        .expect_err("no common satellites must be an explicit error");
    assert!(err.contains("double-difference"), "unexpected message: {err}");
}

/// Forcing a total satellite loss (rover sees nothing) must leave the filter
/// usable: finite state, non-negative variances, unchanged dimension, and the
/// NEXT good epoch must still solve.
#[test]
fn engine_survives_a_forced_total_satellite_loss() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    for i in 0..3 {
        eng.process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
            .unwrap_or_else(|e| panic!("epoch {i} failed: {e}"));
    }
    let dim_before = eng.state.dim();
    let blank = EpochObs { time: sim.rover_epochs[3].time, satellites: Vec::new() };
    assert!(
        eng.process_epoch(&blank, &sim.base_epochs[3], cfg.base_ecef, &sim.ephemerides).is_err(),
        "an empty rover epoch must be rejected"
    );
    assert!(eng.state.cov.iter().all(|v| v.is_finite()));
    assert!(eng.state.cov.diagonal().iter().all(|v| *v >= 0.0));
    assert_eq!(eng.state.dim(), dim_before, "a rejected epoch changed the state dimension");

    eng.process_epoch(&sim.rover_epochs[4], &sim.base_epochs[4], cfg.base_ecef, &sim.ephemerides)
        .expect("recovery epoch must process");
    let err = (eng.state.pos_ecef - sim.truth_positions[4].1).norm();
    assert!(err < 2.0, "position error {err:.3} m after a forced outage is unusable");
}

fn strip_sat(epoch: &EpochObs, prn: u8) -> EpochObs {
    let mut e = epoch.clone();
    e.satellites
        .retain(|s| !(s.sat.constellation == Constellation::Gps && s.sat.prn == prn));
    e
}

/// Forcing a reference switch must not corrupt the covariance or teleport the
/// solution. The reference in force after the warm-up is read out of the
/// engine, then that satellite is deleted from every subsequent rover and base
/// epoch, which forces the hysteresis onto a new reference.
#[test]
fn engine_survives_a_forced_reference_switch() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    let warm = 6usize;
    for i in 0..warm {
        eng.process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
            .unwrap_or_else(|e| panic!("warm-up epoch {i} failed: {e}"));
    }
    let old_ref_u16 = eng
        .ref_sats
        .get(&0)
        .copied()
        .expect("the reference must be chosen during warm-up");
    let old_ref = sim.rover_epochs[warm]
        .satellites
        .iter()
        .find(|s| s.sat.constellation == Constellation::Gps && u16::from(s.sat.prn) == old_ref_u16)
        .map(|s| s.sat.prn)
        .expect("the reference satellite must be present in the epoch");
    let mut prev_pos = eng.state.pos_ecef;
    let mut switched_at = None;

    for i in warm..sim.rover_epochs.len() {
        let rover = strip_sat(&sim.rover_epochs[i], old_ref);
        let base = strip_sat(&sim.base_epochs[i], old_ref);
        let sol = eng
            .process_epoch(&rover, &base, cfg.base_ecef, &sim.ephemerides)
            .unwrap_or_else(|e| panic!("epoch {i} failed: {e}"));
        assert!(sol.position_ecef.iter().all(|v| v.is_finite()), "epoch {i}");
        let cov = &eng.state.cov;
        assert!((cov - cov.transpose()).norm() < 1e-6, "epoch {i} asymmetry after ref switch");
        assert!(cov.diagonal().iter().all(|v| *v > 0.0), "epoch {i} lost a positive variance");
        let eig = cov.clone().symmetric_eigenvalues();
        assert!(eig[0] > 0.0, "epoch {i} covariance lost positive definiteness (min eig {})", eig[0]);
        let now = eng.ref_sats.get(&0).copied().unwrap_or(0);
        if switched_at.is_none() && now != old_ref_u16 {
            switched_at = Some(i);
            let jump = (sol.position_ecef - prev_pos).norm();
            assert!(jump < 5.0, "reference switch {old_ref_u16}->{now} at epoch {i} moved the solution {jump:.3} m");
        }
        prev_pos = sol.position_ecef;
        // Once the new reference is locked in, every later epoch must be a
        // normal few-cm solution.
        if let Some(k) = switched_at {
            if i > k + 3 {
                let err = (sol.position_ecef - sim.truth_positions[i].1).norm();
                assert!(err < 1.0, "epoch {i}: {err:.3} m error after the reference switch");
            }
        }
    }
    assert!(
        switched_at.is_some(),
        "dropping PRN {old_ref} never forced a reference switch; the test proved nothing"
    );
}

/// `min_ar_lock_epochs` must actually gate AR: with a threshold above the run
/// length, NO pair ever becomes eligible, so no fix may be declared.
#[test]
fn min_ar_lock_epochs_above_the_run_length_prevents_every_fix() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    eng.min_ar_lock_epochs = 10_000;
    let mut any_fixed = false;
    for i in 0..sim.rover_epochs.len() {
        let sol = eng
            .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
            .expect("epoch must process");
        any_fixed |= sol.is_fixed;
    }
    assert!(
        !any_fixed,
        "no pair can reach 10 000 tracked epochs in a {} s run, yet a fix was declared",
        cfg.duration_s
    );
}

/// With the gate OFF (the default) the same data may fix; that is a
/// configuration question, not an invariant, so this test only records the
/// gate's observable effect on the ambiguity state rather than freezing it.
#[test]
fn min_ar_lock_epochs_gate_drops_young_pairs_from_the_ar_state() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let mut gated = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
    gated.min_ar_lock_epochs = 1_000_000;
    for i in 0..6 {
        gated
            .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides)
            .expect("epoch must process");
    }
    assert_eq!(
        gated.state.ambiguities.len(),
        0,
        "with an unreachable lock threshold every pair is ineligible and the AR \\
         ambiguity state must be empty"
    );
    assert!(gated.state.cov.diagonal().iter().all(|v| *v > 0.0));
}

// ===========================================================================
// Engine construction defaults (wiring, not policy)
// ===========================================================================

#[test]
fn engine_construction_records_its_anchor_time_and_live_states() {
    let t0 = GpsTime::new(2200, 300_000.0);
    let eng = GnssRtkIekf::new(BASE_ECEF, t0, 0.5);
    assert_eq!(eng.start_tow, t0.tow, "the two-phase static Q anchor must be the session start");
    assert_eq!(eng.state.time, t0);
    assert_eq!(eng.q_accel, 0.5);
    assert_eq!(eng.state.dim(), 6, "a fresh state has position + velocity only");
    assert!(eng.ref_sats.is_empty());
    assert!(!eng.track_ambiguity_keys);
    assert_eq!(eng.target_pf, 0.001);
    assert_eq!(eng.min_elevation_rad, 0.1745);
}

/// `process_epoch` must be a pure function of its inputs w.r.t. the reference
/// map: with an identical dataset the same epochs must select the same
/// reference satellites (the hysteresis must not depend on wall-clock or hash
/// iteration order).
#[test]
fn processing_is_reproducible_for_identical_input() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    let run = || -> Vec<(u16, u8, Vector3<f64>)> {
        let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        (0..sim.rover_epochs.len())
            .map(|i| {
                let s = eng
                    .process_epoch(
                        &sim.rover_epochs[i],
                        &sim.base_epochs[i],
                        cfg.base_ecef,
                        &sim.ephemerides,
                    )
                    .expect("epoch must process");
                (eng.ref_sats.get(&0).copied().unwrap_or(0), s.quality, s.position_ecef)
            })
            .collect()
    };
    assert_eq!(run(), run());
}

// ===========================================================================
// Observation plumbing sanity (the fixture itself)
// ===========================================================================

#[test]
fn the_simulated_fixture_produces_real_observations_for_both_receivers() {
    let cfg = sim_config();
    let sim = generate_simulation_dataset(&cfg);
    for e in [&sim.base_epochs[0], &sim.rover_epochs[0]] {
        assert!(e.satellites.len() >= 10);
        for s in &e.satellites {
            let pr = s.observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange);
            let cp = s.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase);
            assert!(pr.is_some() && cp.is_some(), "fixture satellite {:?} lacks observables", s.sat);
            assert!(pr.unwrap().value > 1.0e7 && cp.unwrap().value > 1.0e7);
        }
    }
    let _ = Observation {
        code: ObsCode { obs_type: ObsType::Snr, signal: SignalCode { freq_band: 1, attribute: 'C' } },
        value: 0.0,
        lock_time: None,
        lli: None,
    };
}