//! Sprint-hunt tests for the real-time streaming RTK engine
//! (`streaming.rs`). The three pre-existing tests all drive `dummy_epoch`,
//! whose SatObs carries `observations: vec![]`, so DD formation, screening, the
//! Kalman update and the smoother are never entered. Everything below drives
//! the engine with the real observation simulator.
#![allow(clippy::unwrap_used)]

use super::*;

use crate::sim::generator::{generate_simulation_dataset, SimulationConfig, TrajectoryProfile};

const WEEK: u32 = 2100;

fn t(tow: f64) -> GpsTime {
    GpsTime::new(WEEK, tow)
}


    fn sim_config() -> SimulationConfig {
        SimulationConfig {
            duration_s: 8.0,
            num_satellites: 16,
            profile: TrajectoryProfile::Static { offset_ned: Vector3::new(30.0, 20.0, 0.0) },
            ..Default::default()
        }
    }

    fn engine() -> (StreamingRtkEngine, crate::sim::generator::SimulationDataset, SimulationConfig) {
        let cfg = sim_config();
        let sim = generate_simulation_dataset(&cfg);
        let sc = StreamingConfig {
            base_position: cfg.base_ecef,
            max_base_age_s: 2.0,
            ..Default::default()
        };
        (StreamingRtkEngine::new(sc, sim.ephemerides.clone()), sim, cfg)
    }

    /// Baseline: a clean streaming run must produce solutions and keep the
    /// covariance symmetric positive-definite and the state dimension fixed.
    #[test]
    fn clean_streaming_run_is_finite_symmetric_and_fixed_dimension() {
        let (mut eng, sim, cfg) = engine();
        let mut dim_seen = None;
        let mut n = 0;
        for i in 0..sim.rover_epochs.len() {
            eng.push_base_obs(sim.base_epochs[i].clone());
            let Some(sol) = eng.process_rover_obs(&sim.rover_epochs[i]) else { continue };
            n += 1;
            assert!(sol.position_ecef.iter().all(|v| v.is_finite()), "epoch {i} position non-finite");
            assert!(sol.std_east.is_finite() && sol.std_north.is_finite() && sol.std_up.is_finite());
            let state = &eng.iekf.as_ref().unwrap().state;
            if let Some(d) = dim_seen {
                assert_eq!(state.dim(), d, "state dimension changed at epoch {i}");
            } else {
                dim_seen = Some(state.dim());
            }
            let cov = &state.cov;
            assert!((cov - cov.transpose()).norm() < 1e-6, "epoch {i} covariance asymmetry");
            let e = cov.clone().symmetric_eigenvalues();
            assert!(e[0] > 0.0, "epoch {i} covariance lost positive definiteness (min eig {})", e[0]);
        }
        assert!(n >= 6, "expected solutions for most epochs, got {n}");
        let _ = cfg;
    }

    // ---------------------------------------------------------------------
    // BUG C — the base buffer is keyed on time-of-week only; the GPS week is
    //          discarded, so an epoch one week away is accepted as "current".
    // ---------------------------------------------------------------------

    /// `push_base_obs` keys the buffer on `(tow * 1000).round() as u64` and
    /// drops `GpsTime.week` entirely. A base epoch from week 2100 and one from
    /// week 2101 therefore collide on the same key, and `find_closest_base`
    /// reports a match across a 604 800 s separation.
    #[test]
    #[ignore = "BUG-3 (streaming.rs:86-97, 123-136): push_base_obs keys the buffer on\n(tow*1000).round() and drops GpsTime.week, so find_closest_base matches a base\nepoch 604800 s away. Re-enable once the buffer key includes the GPS week."]
    fn base_buffer_must_not_match_across_a_gps_week_boundary() {
        let mut eng = StreamingRtkEngine::new(StreamingConfig::default(), Vec::new());
        let old = EpochObs { time: t(100.0), satellites: Vec::new() };
        eng.push_base_obs(old);

        // Rover exactly one week later, same time-of-week.
        let rover_next_week = EpochObs { time: GpsTime::new(WEEK + 1, 100.0), satellites: Vec::new() };
        let tow_ms = (rover_next_week.time.tow * 1000.0).round() as u64;
        assert!(
            eng.find_closest_base(tow_ms).is_none(),
            "a base epoch 604800 s stale must not satisfy max_base_age_s = {:?}",
            eng.config.max_base_age_s
        );
    }

    /// The two weeks collide on the same buffer key, so the second base epoch
    /// silently evicts the first. Both are "the last base epoch received", but
    /// they are a week apart.
    #[test]
    #[ignore = "BUG-3 (streaming.rs:86-87): epochs one week apart collide on the same u64 buffer\nkey, so the newer base epoch silently evicts the older one. Re-enable with BUG-3."]
    fn base_buffer_distinguishes_weeks() {
        let mut eng = StreamingRtkEngine::new(StreamingConfig::default(), Vec::new());
        eng.push_base_obs(EpochObs { time: t(100.0), satellites: Vec::new() });
        let key_w0 = eng.base_buffer.keys().next().copied().expect("one entry");
        eng.push_base_obs(EpochObs { time: GpsTime::new(WEEK + 1, 100.0), satellites: Vec::new() });
        assert_eq!(
            eng.base_buffer.len(),
            2,
            "one week apart must be two distinct buffered epochs, not one overwritten key"
        );
        let _ = key_w0;
    }

    // ---------------------------------------------------------------------
    // BUG D — no re-initialisation after a gap: the constant-velocity model is
    //          integrated straight across it and the position jumps by v*dt.
    // ---------------------------------------------------------------------

    /// `process_rover_obs` has no gap detector at all. Between two epochs a
    /// full week apart (still perfectly matchable by a base epoch, because of
    /// BUG C) `predict_state_gated` applies `pos += vel * dt` with
    /// dt = 604 800 s. A filtered RTK stream must not integrate a
    /// constant-velocity model across a gap of that size; it must drop and
    /// re-initialise the filter.
    #[test]
    fn streaming_must_reinitialise_across_a_long_gap() {
        let (mut eng, sim, _cfg) = engine();
        for i in 0..4 {
            eng.push_base_obs(sim.base_epochs[i].clone());
            eng.process_rover_obs(&sim.rover_epochs[i]).expect("warm-up epoch must solve");
        }
        let pos_before = eng.iekf.as_ref().unwrap().state.pos_ecef;
        let vel = eng.iekf.as_ref().unwrap().state.vel_ecef;
        assert!(vel.norm() > 1e-3, "filter should have learned a non-zero velocity for this test to bite");

        // Same simulator data, replayed a week later with matching base epochs.
        let week_later = 7.0 * 86_400.0;
        for i in 0..4 {
            let mut b = sim.base_epochs[i].clone();
            b.time = GpsTime::new(WEEK + 1, b.time.tow + week_later);
            eng.push_base_obs(b);
        }
        let mut rover = sim.rover_epochs[0].clone();
        rover.time = GpsTime::new(WEEK + 1, rover.time.tow + week_later);
        let sol = eng.process_rover_obs(&rover).expect("engine must still answer");

        let straight_line_jump = (vel * week_later).norm();
        assert!(
            (sol.position_ecef - pos_before).norm() < straight_line_jump * 0.5,
            "position jumped {:.1} m across a 7-day gap; dead reckoning at v={:.3} m/s \
             alone would move {:.1} m, so the gap was integrated rather than reset",
            (sol.position_ecef - pos_before).norm(),
            vel.norm(),
            straight_line_jump
        );
    }

    // ---------------------------------------------------------------------
    // BUG E — a NaN observation poisons the filter permanently.
    // ---------------------------------------------------------------------

    /// `screen::screen_gross_pr_errors` breaks out of its loop on the first
    /// non-comparable residual (`if residual.abs() <= 15.0 { break }` is false
    /// for NaN), so a NaN pseudorange is NOT screened. It then flows into
    /// `compute_iekf_step`, whose `S` pivot check is `diag.is_zero()`, also
    /// false for NaN, so no error is raised either — the state and covariance
    /// simply become NaN and stay NaN forever.
    #[test]
    fn a_single_nan_observation_must_not_poison_the_filter() {
        let (mut eng, sim, _cfg) = engine();
        for i in 0..3 {
            eng.push_base_obs(sim.base_epochs[i].clone());
            eng.process_rover_obs(&sim.rover_epochs[i]).expect("warm-up epoch must solve");
        }

        // Inject one NaN pseudorange. There is NO explicit is_finite guard on
        // the DD update path; survival here relies on `worst_pr_residual`'s
        // `total_cmp` ordering putting NaN last in the max_by scan, which then
        // de-weights (or drops) that row. See
        // rtk_iekf::sprint_tests::iekf_update_rejects_non_finite_measurements
        // for the lower-level guarantee this must rest on.
        let mut poisoned = sim.rover_epochs[3].clone();
        for s in poisoned.satellites.iter_mut() {
            if let Some(o) = s.observations.iter_mut().find(|o| {
                matches!(o.code.obs_type, gneiss_core::obs::ObsType::Pseudorange)
            }) {
                o.value = f64::NAN;
                break;
            }
        }

        // Then three perfectly clean epochs. The filter must recover.
        for i in 4..7 {
            eng.push_base_obs(sim.base_epochs[i].clone());
            let sol = eng
                .process_rover_obs(&sim.rover_epochs[i])
                .unwrap_or_else(|| panic!("clean epoch {i} must still produce a solution"));
            assert!(
                sol.position_ecef.iter().all(|v| v.is_finite()),
                "epoch {i}: one NaN pseudorange permanently NaN'd the filter -> {:?}",
                sol.position_ecef
            );
        }
    }

    // ---------------------------------------------------------------------
    // Buffer bounds / satellite-loss behaviour
    // ---------------------------------------------------------------------

    /// `push_base_obs` documents a 30 s retention horizon. The invariant is
    /// therefore "every buffered epoch is within 30 s of the newest", NOT a
    /// fixed entry count: at 20 Hz that is ~601 entries, and the code's prune
    /// condition (`tow_ms - oldest > 30_000`) is what enforces it.
    #[test]
    fn base_buffer_never_retains_epochs_older_than_30_seconds() {
        const DT_MS: u64 = 50; // 20 Hz base feed
        let mut eng = StreamingRtkEngine::new(StreamingConfig::default(), Vec::new());
        for k in 0..4000i64 {
            let tow = 100.0 + k as f64 * (DT_MS as f64) / 1000.0;
            eng.push_base_obs(EpochObs { time: t(tow), satellites: Vec::new() });
            let (newest, oldest) = {
                let newest = *eng.base_buffer.keys().next_back().unwrap();
                let oldest = *eng.base_buffer.keys().next().unwrap();
                (newest, oldest)
            };
            assert!(
                newest - oldest <= 30_000 + DT_MS,
                "after {} epochs the buffer spans {} ms (> 30 s horizon + one epoch)",
                k,
                newest - oldest
            );
        }
    }

    /// Out-of-order arrival: an old base epoch must not be retained forever and
    /// must not evict the current one.
    #[test]
    fn base_buffer_survives_out_of_order_arrival() {
        let mut eng = StreamingRtkEngine::new(StreamingConfig::default(), Vec::new());
        for k in 0..40i64 {
            eng.push_base_obs(EpochObs { time: t(100.0 + k as f64 * 0.1), satellites: Vec::new() });
        }
        let newest = *eng.base_buffer.keys().next_back().unwrap();
        // A straggler from 5 s ago.
        eng.push_base_obs(EpochObs { time: t(100.0 + 39.0 * 0.1 - 5.0), satellites: Vec::new() });
        assert_eq!(
            *eng.base_buffer.keys().next_back().unwrap(),
            newest,
            "a straggler base epoch must not evict the newest one"
        );
    }

    /// Losing every satellite on one epoch must not corrupt the engine: the
    /// next good epoch must still solve and the state dimension must be
    /// unchanged.
    #[test]
    fn satellite_loss_then_recovery_keeps_a_usable_state() {
        let (mut eng, sim, _cfg) = engine();
        for i in 0..3 {
            eng.push_base_obs(sim.base_epochs[i].clone());
            eng.process_rover_obs(&sim.rover_epochs[i]).expect("warm-up epoch must solve");
        }
        let dim_before = eng.iekf.as_ref().unwrap().state.dim();
        let cov_before = eng.iekf.as_ref().unwrap().state.cov.clone();

        // Epoch 3: rover sees nothing at all.
        let blank = EpochObs { time: sim.rover_epochs[3].time, satellites: Vec::new() };
        eng.push_base_obs(sim.base_epochs[3].clone());
        assert!(eng.process_rover_obs(&blank).is_none(), "no DD data must yield no solution");
        let state = &eng.iekf.as_ref().unwrap().state;
        assert!(state.cov.iter().all(|v| v.is_finite()), "total outage made covariance non-finite");
        assert!(
            state.cov.diagonal().iter().all(|v| *v >= 0.0),
            "total outage produced a negative variance: {:?}",
            state.cov.diagonal()
        );
        assert_eq!(state.dim(), dim_before, "outage changed the state dimension");

        // A rejected epoch carries no information, so no state may have moved.
        let st = eng.iekf.as_ref().unwrap();
        assert!(
            st.state.cov[(0, 0)] >= cov_before[(0, 0)],
            "a fully-rejected epoch reduced position variance from {} to {}",
            cov_before[(0, 0)],
            st.state.cov[(0, 0)]
        );

        // Epoch 4: normal again.
        eng.push_base_obs(sim.base_epochs[4].clone());
        let sol = eng.process_rover_obs(&sim.rover_epochs[4]).expect("recovery epoch must solve");
        assert!(sol.position_ecef.iter().all(|v| v.is_finite()));
        let truth = sim.truth_positions[4].1;
        let err = (sol.position_ecef - truth).norm();
        assert!(err < 2.0, "position error {err:.3} m after a total outage is unusable");
    }

    /// `push_ephemeris` replaces by satellite; the count must not grow when the
    /// same satellite is re-pushed, and must grow for a new one.
    #[test]
    fn ephemeris_catalogue_replaces_by_satellite() {
        let sim = generate_simulation_dataset(&sim_config());
        let mut eng = StreamingRtkEngine::new(StreamingConfig::default(), Vec::new());
        let a = sim.ephemerides[0].clone();
        let b = sim.ephemerides[1].clone();
        eng.push_ephemeris(a.clone());
        eng.push_ephemeris(b.clone());
        assert_eq!(eng.ephemerides_len(), 2);
        eng.push_ephemeris(a);
        assert_eq!(eng.ephemerides_len(), 2, "re-pushing a known satellite must replace, not append");
    }
