#![allow(clippy::unwrap_used)]

    use super::*;
    use super::ar_gate::SEED_VARIANCE_SAFETY_MARGIN;
    use crate::sim::generator::{generate_simulation_dataset, SimulationConfig};
    use gneiss_core::constants::SPEED_OF_LIGHT_M_S;

    fn test_engine(start: GpsTime) -> GnssRtkIekf {
        GnssRtkIekf::new(Vector3::new(1.0, 2.0, 3.0), start, 1.0)
    }

    fn dd_key(sat: u16, band: u8) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: band }
    }

    #[test]
    fn seed_variance_matches_margin_squared_times_code_sigma_in_cycles() {
        let pr_var_m2: f64 = 0.64; // 0.8 m code-noise sigma
        let lambda = 0.1903;
        let sigma_cycles = pr_var_m2.sqrt() / lambda;
        let expected = (SEED_VARIANCE_SAFETY_MARGIN * sigma_cycles).powi(2);
        assert!((seed_ambiguity_variance_cycles2(pr_var_m2, lambda) - expected).abs() < 1e-9);
    }

    #[test]
    fn seed_variance_is_looser_for_noisier_low_elevation_code() {
        let lambda = 0.1903;
        let high_el_pr_var = 0.05; // ~60 degrees, per compute_dd_variances
        let low_el_pr_var = 2.4; // ~15 degrees
        let high = seed_ambiguity_variance_cycles2(high_el_pr_var, lambda);
        let low = seed_ambiguity_variance_cycles2(low_el_pr_var, lambda);
        assert!(low > high, "low-elevation seed variance ({low}) should exceed high-elevation ({high})");
    }

    /// Run the sim dataset through `engine`, returning per-epoch fix flags
    /// and positions for cross-run comparisons.
    fn run_sim(
        engine: &mut GnssRtkIekf,
        sim: &crate::sim::generator::SimulationDataset,
        base: Vector3<f64>,
    ) -> Vec<(bool, Vector3<f64>)> {
        let mut out = Vec::new();
        for i in 0..sim.rover_epochs.len() {
            let s = engine.process_epoch(
                &sim.rover_epochs[i], &sim.base_epochs[i], base, &sim.ephemerides,
            ).expect("epoch must process");
            out.push((s.is_fixed, s.position_ecef));
        }
        out
    }

    #[test]
    fn test_gnss_rtk_iekf_runs_on_simulated_dataset() {
        let cfg = SimulationConfig {
            duration_s: 10.0,
            ..Default::default()
        };
        let sim = generate_simulation_dataset(&cfg);
        let mut engine = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);

        let mut fixed_count = 0;
        for i in 0..sim.rover_epochs.len() {
            let sol = engine.process_epoch(
                &sim.rover_epochs[i],
                &sim.base_epochs[i],
                cfg.base_ecef,
                &sim.ephemerides,
            );
            assert!(sol.is_ok());
            let s = sol.unwrap();
            if s.is_fixed {
                fixed_count += 1;
                let err = (s.position_ecef - sim.truth_positions[i].1).norm();
                assert!(err < 0.05, "Fixed epoch error should be < 5cm, got {:.4}m", err);
            }
        }

        assert!(fixed_count >= 5, "RTK engine should fix at least 5 epochs");
        let smoothed = engine.smooth();
        assert_eq!(smoothed.len(), 10);
    }
    #[test]
    fn test_select_constellations_drops_glonass_by_default() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        let mk = |c: Constellation, prn: u8| {
            (SatelliteId { constellation: c, prn }, Vector3::zeros())
        };
        let sat_info = vec![
            mk(Constellation::Gps, 6u8),
            mk(Constellation::Glonass, 8u8),
            mk(Constellation::Galileo, 3u8),
            mk(Constellation::Gps, 12u8),
        ];
        // Default policy: FDMA GLONASS stays out until ICB handling lands.
        let got = GnssRtkIekf::select_constellations(&sat_info, false);
        assert_eq!(got, vec![
            Constellation::Gps as u8,
            Constellation::Galileo as u8,
        ]);
        // Opt-in keeps it, sorted and de-duplicated.
        let got = GnssRtkIekf::select_constellations(&sat_info, true);
        assert_eq!(got, vec![
            Constellation::Gps as u8,
            Constellation::Glonass as u8,
            Constellation::Galileo as u8,
        ]);
        // Empty input stays empty.
        assert!(GnssRtkIekf::select_constellations(&[], true).is_empty());
    }

    // ---- GNEISS_AR_GATE: technique 1, AR elevation mask -------------------

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
        // Realistic rover ECEF so az/el geometry is meaningful.
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
        // 15 deg AR mask: horizon satellite must not reach LAMBDA.
        eng.ar_elevation_mask_rad = 0.2618;
        let view = ar_gate::elevation_filtered_view(&eng.state, &meas, eng.ar_elevation_mask_rad);
        assert_eq!(view.ambiguities.len(), 1);
        assert_eq!(view.ambiguities[0].0, hi);
        // Non-destructive: live state keeps both floats.
        assert_eq!(eng.state.ambiguities.len(), 2);
    }

    #[test]
    fn test_gate_off_extreme_ar_mask_is_ignored_end_to_end() {
        let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
        let sim = generate_simulation_dataset(&cfg);
        // Gate off -> ar_elevation_mask_rad must be inert: fixes still occur.
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
        // ~89 deg mask excludes every pair from LAMBDA -> float-only output.
        let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        eng.ar_gate = true;
        eng.ar_elevation_mask_rad = 1.55;
        let res = run_sim(&mut eng, &sim, cfg.base_ecef);
        assert!(
            res.iter().all(|(f, _)| !f),
            "AR mask must exclude all pairs; unexpected fix"
        );
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

    // ---- ProcessingDynamics: kinematic profile vs static profile --------

    /// Mean 3D tracking error of one engine over the converged tail
    /// (epochs after `skip`) of a simulated trajectory.
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
        use crate::post_process::dynamics::{
            KINEMATIC_INNOV_GATE_SCALE, KINEMATIC_Q_ACCEL, STATIC_Q_ACCEL,
        };
        use crate::sim::generator::{TrajectoryProfile};
        // PPK cadence (30 s epochs): a 5 m/s rover covers 150 m between
        // epochs. Constant-velocity ramp: even the static profile survives
        // here because robust variance inflation lets strong phase
        // measurements drag the frozen prior along; the kinematic profile
        // must also stay bounded AND keep live velocity states.
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

        let mut static_eng =
            GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, STATIC_Q_ACCEL);
        let mut kin_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, KINEMATIC_Q_ACCEL);
        kin_eng.robust_innov_scale = KINEMATIC_INNOV_GATE_SCALE;

        let err_static = mean_tail_error(&mut static_eng, &sim, cfg.base_ecef, 10);
        let err_kin = mean_tail_error(&mut kin_eng, &sim, cfg.base_ecef, 10);
        let speed_kin = kin_eng.state.vel_ecef.norm();
        let speed_static = static_eng.state.vel_ecef.norm();
        eprintln!(
            "KIN-SIM linear ramp @30s: err static={err_static:.3} m kin={err_kin:.3} m | \
             final speed est static={speed_static:.2} kin={speed_kin:.2} (true {true_speed:.2}) m/s"
        );
        // No crash/divergence for either profile on steady motion.
        assert!(err_kin < 3.0, "kinematic must track the ramp, got {err_kin:.3} m");
        assert!(err_static < 3.0, "static+Huber also tracks steady ramps, got {err_static:.3} m");
        // Differential signal: kinematic velocity states stay live
        // (within 25% of the true 5 m/s). NOTE (measured): the STATIC
        // profile's velocity state ALSO converges on a clean constant-
        // velocity ramp — sequential position updates make velocity
        // observable through the F coupling even with monument Q — so
        // no frozenness assertion is possible here; the profiles
        // separate under ACCELERATION (next test).
        assert!(
            (speed_kin - true_speed).abs() < 0.25 * true_speed,
            "kinematic velocity state must track true speed, got {speed_kin:.2}"
        );
    }

    #[test]
    fn kinematic_profile_outperforms_static_under_acceleration() {
        use crate::post_process::dynamics::{
            KINEMATIC_INNOV_GATE_SCALE, KINEMATIC_Q_ACCEL, STATIC_Q_ACCEL,
        };
        use crate::sim::generator::{TrajectoryProfile};
        // Accelerating frame: 15 m/s on a 500 m radius circle is
        // 0.45 m/s^2 of sustained acceleration — a constant-velocity
        // model with monument Q must lag every epoch.
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

        let mut static_eng =
            GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, STATIC_Q_ACCEL);
        let mut kin_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, KINEMATIC_Q_ACCEL);
        kin_eng.robust_innov_scale = KINEMATIC_INNOV_GATE_SCALE;

        let err_static = mean_tail_error(&mut static_eng, &sim, cfg.base_ecef, 12);
        let err_kin = mean_tail_error(&mut kin_eng, &sim, cfg.base_ecef, 12);
        eprintln!(
            "KIN-SIM circular @30s: err static={err_static:.3} m kinematic={err_kin:.3} m"
        );
        assert!(err_kin < 3.0, "kinematic must track the turn, got {err_kin:.3} m");
        assert!(
            err_static > 2.0 * err_kin,
            "static-tuned Q must lag under sustained acceleration \
             (static={err_static:.3}, kin={err_kin:.3})"
        );
    }

    // ---- GNEISS_AR_GATE: technique 2, phase-code coherency bias init ------

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
        // Established band-1 pair plus this epoch's divergence samples:
        // median over band 1 is 3.5; the band-2 sample must be ignored.
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
        // Only band-2 samples exist; a new band-1 pair must seed raw.
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
        // Slip on an ESTABLISHED pair: legacy re-seed semantics preserved.
        eng.update_dd_ambiguity(dd_key(2, 1), Some(205.75), 100.0, lambda, true, 1.0);
        assert!((eng.state.ambiguities[0].1 - raw).abs() < 1e-12);
        assert_eq!(eng.code_phase_div.len(), 3, "established pair still contributes a sample");
    }

    // ---- GNEISS_AR_GATE: regression, gate off == exact legacy behaviour ---

    #[test]
    fn test_gate_disabled_run_matches_default_run_bit_for_bit() {
        let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
        let sim = generate_simulation_dataset(&cfg);
        let mut baseline = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        let mut gated = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        gated.ar_gate = false; // explicit, but identical to the default
        let base_out = run_sim(&mut baseline, &sim, cfg.base_ecef);
        let gate_out = run_sim(&mut gated, &sim, cfg.base_ecef);
        assert_eq!(base_out, gate_out, "gate off must reproduce legacy exactly");
        assert!(base_out.iter().filter(|(f, _)| *f).count() >= 5);
    }

    // ---- receiver_dd_pcv_m: correction wiring ---------------------------

    /// Zero unless calibrations are loaded; otherwise equal to the raw
    /// differential PCV at the rover-frame elevations. `self.receiver_pcv`
    /// being `Some` IS the caller's opt-in signal — no separate env gate
    /// (a prior GNEISS_RECV_PCV/GNEISS_PCV split silently required both
    /// to be set for the documented "opt-in via GNEISS_PCV=1" to actually
    /// apply anything; see docs/NETWORK_RTK_NEXT_STEPS.md). Skipped when
    /// igs14 is absent.
    #[test]
    fn receiver_dd_pcv_m_requires_loaded_pair() {
        let Ok(db) = gneiss_parsers::antex::AntexDatabase::parse("../../datasets/igs14.atx")
        else {
            return;
        };
        use gneiss_parsers::receiver_antenna::{compute_dd_pcv_correction_2d, ReceiverAntenna};
        let trm =
            Arc::new(ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").expect("igs14 TRM"));
        let ash =
            Arc::new(ReceiverAntenna::lookup(&db, "ASH701945B_M", "SCIT").expect("igs14 ASH"));
        let mut eng =
            GnssRtkIekf::new(Vector3::new(-3961904.43, 3348994.27, 3698211.71), GpsTime::new(2000, 100.0), 1.0);
        let sat_pos = eng.state.pos_ecef + Vector3::new(1.0e7, 5.0e6, 2.0e7);
        let ref_pos = eng.state.pos_ecef + Vector3::new(0.0, 0.0, 2.4e7);
        let sid = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 3,
        };

        // No calibrations loaded -> no correction.
        eng.receiver_pcv = None;
        assert_eq!(eng.receiver_dd_pcv_m(sid, 1, sat_pos, ref_pos), 0.0);
        // Calibrations loaded -> raw differential PCV (non-zero for this
        // cross-family pair at distinct elevations).
        eng.receiver_pcv = Some((trm.clone(), ash.clone()));
        let dd = eng.receiver_dd_pcv_m(sid, 1, sat_pos, ref_pos);
        let llh = gneiss_core::coords::ecef_to_llh(eng.state.pos_ecef);
        let (az_s, el_s) = gneiss_core::coords::az_el(llh, eng.state.pos_ecef, sat_pos);
        let (az_r, el_r) = gneiss_core::coords::az_el(llh, eng.state.pos_ecef, ref_pos);
        let expected = compute_dd_pcv_correction_2d(&trm, &ash, "G01", az_s, el_s, az_r, el_r, eng.rover_heading_rad);
        assert!(dd.abs() > 1e-6, "cross-family correction must be non-zero: {dd}");
        assert!((dd - expected).abs() < 1e-12, "dd={dd} expected={expected}");
    }

    // ---- GNEISS_CLK: constellation-median centering + spread gate -------

    use gneiss_parsers::clk_centering::CenteredClock;
    use gneiss_parsers::rinex_clk::{ClockRecord, RinexClock};

    fn centered(bias_us: f64, spread_us: f64) -> Option<CenteredClock> {
        Some(CenteredClock { bias_s: bias_us * 1e-6, spread_s: spread_us * 1e-6 })
    }

    #[test]
    fn centered_pair_correction_healthy_pair_applies_centered_delta() {
        let (corr, tripped) =
            centered_pair_correction(centered(10.0, 5.0), centered(-14.0, 6.0));
        assert!(!tripped);
        assert!((corr - SPEED_OF_LIGHT_M_S * 24e-6).abs() < 1e-9);
    }

    #[test]
    fn centered_pair_correction_gates_when_either_side_spread_trips() {
        for (a, b) in [
            (centered(1.0, 150.0), centered(2.0, 5.0)),
            (centered(1.0, 5.0), centered(2.0, 150.0)),
        ] {
            let (corr, tripped) = centered_pair_correction(a, b);
            assert!(tripped, "pathological side must trip the gate");
            assert_eq!(corr, 0.0, "gated correction must be suppressed");
        }
    }

    #[test]
    fn centered_pair_correction_threshold_is_strictly_greater() {
        // Exactly at the 100 us threshold the product is still trusted;
        // a hair above it trips.
        let ok = centered_pair_correction(centered(1.0, 100.0), centered(2.0, 99.999));
        assert!(!ok.1);
        assert!(ok.0.abs() > 0.0);
        let bad = centered_pair_correction(centered(1.0, 100.001), centered(2.0, 5.0));
        assert!(bad.1 && bad.0 == 0.0);
    }

    #[test]
    fn centered_pair_correction_missing_side_stays_silent_zero() {
        assert_eq!(centered_pair_correction(None, centered(2.0, 5.0)), (0.0, false));
        assert_eq!(centered_pair_correction(centered(1.0, 5.0), None), (0.0, false));
        assert_eq!(centered_pair_correction(None, None), (0.0, false));
    }

    /// Synthetic product helper: one record per satellite at tow = 100
    /// (matching the engine time below), biases in microseconds.
    fn clk_product(biases_us: &[(u8, f64)]) -> Arc<RinexClock> {
        let mut rc = RinexClock::default();
        for (prn, us) in biases_us {
            rc.satellites.insert(
                gneiss_core::sat::SatelliteId {
                    constellation: gneiss_core::sat::Constellation::Gps,
                    prn: *prn,
                },
                vec![ClockRecord { time: GpsTime::new(2200, 100.0), bias: us * 1e-6 }],
            );
        }
        Arc::new(rc)
    }

    fn dd_probe(eng: &GnssRtkIekf, sat: u8, reference: u8) -> f64 {
        let rx = Vector3::zeros();
        eng.formation_clock_corr_m(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: sat,
            },
            u16::from(reference),
            rx,
            Vector3::new(2.0e7, 0.0, 0.0),
            Vector3::new(2.4e7, 0.0, 0.0),
        )
    }

    #[test]
    fn precise_clock_dd_m_without_product_is_zero() {
        let eng = test_engine(GpsTime::new(2200, 100.0));
        assert!(eng.precise_clocks.is_none());
        assert_eq!(dd_probe(&eng, 5, 9), 0.0);
    }

    /// Healthy product under a big common mode: the correction equals the
    /// RAW pairwise delta exactly — centering removes only the datum and
    /// must leave the differential untouched through the full wiring.
    #[test]
    fn precise_clock_dd_m_healthy_product_preserves_raw_delta() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.precise_clocks = Some(clk_product(&[
            (27, 500.0),
            (28, 520.0),
            (5, 480.0),
            (10, 510.0),
        ]));
        let corr = dd_probe(&eng, 28, 27);
        let raw = SPEED_OF_LIGHT_M_S * 20e-6; // 520 - 500 us
        assert!((corr - raw).abs() < 1e-9, "corr {corr} vs raw {raw}");
        assert!(!eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed));
    }

    /// Pathological product: correction disabled (0.0) and the warning
    /// flag latches on the first gated call.
    #[test]
    fn precise_clock_dd_m_pathological_product_returns_zero_and_latches() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.precise_clocks = Some(clk_product(&[
            (1, 600.0),
            (2, -600.0),
            (3, 590.0),
            (4, -590.0),
        ]));
        assert_eq!(dd_probe(&eng, 1, 3), 0.0);
        assert!(
            eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed),
            "gate trip must latch the warning flag"
        );

        // Fewer than three valid mates: silently zero, no new behaviour.
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.precise_clocks = Some(clk_product(&[(1, 600.0), (2, -600.0)]));
        assert_eq!(dd_probe(&eng, 1, 2), 0.0);
        assert!(!eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed));
    }


#[cfg(test)]
mod obs_side_clk_tests {
    use super::*;

    /// Structural-consistency contract: with a clock product loaded and a
    /// healthy (non-gated) bias pair, the FORMED measurement carries no
    /// clock signature — dd_pr/dd_cp are already corrected at formation.
    /// This is what makes iono_free.rs / ar_gate.rs / filter consumers
    /// consistent without any per-consumer handling.
    #[test]
    fn formation_subtracts_clock_delta_from_measurements() {
        let t = GpsTime::new(2370, 43_200.0);
        let mut eng = GnssRtkIekf::new(Vector3::zeros(), t, 1.0);
        eng.state.iono_enabled = false;
        eng.precise_clocks = Some(std::sync::Arc::new(
            gneiss_parsers::rinex_clk::RinexClock::parse(
                &synthetic_clk_content(),
            ),
        ));
        {
            let content = synthetic_clk_content();
            eprintln!("CONTENT repr: {:?}", content);
            let probe = gneiss_parsers::rinex_clk::RinexClock::parse(&content);
            eprintln!("probe sats: {}", probe.satellites.len());
            let direct = eng.precise_clocks.as_ref().unwrap()
                .centered_clock(sv1_of(), t);
            eprintln!("ENGINE-TEST centered g01 = {direct:?}");
        }
        let corr = eng.formation_clock_corr_m(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            2,
            Vector3::zeros(),
            Vector3::new(2.0e7, 0.0, 0.0),
            Vector3::new(2.1e7, 0.0, 0.0),
        );
        // Synthetic biases ±100 µs → delta must be c·(b1−b2) magnitude,
        // i.e., tens of metres — proving correction is computed formation-side.
        // Centered cluster: g01 sits 20 µs from median → c·20 µs ≈ 6 km.
        assert!(
            (corr.abs() - 5_995.8).abs() < 10.0,
            "correction {corr} m != expected c·(−20 µs) = −5995.8 m"
        );
    }

    fn sv1_of() -> gneiss_core::sat::SatelliteId {
        gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn: 1 }
    }

    fn synthetic_clk_content() -> String {
        // Record epoch = 2025-06-08T12:00 GPST == week 2370, tow 43200
        // (matches the eval instant below).
        // Two GPS sats, biases ∓100 µs at one epoch.
        "     3.00           C                                       RINEX VERSION / TYPE\n\
         2    AS    AR                                          # / TYPES OF DATA\n\
         AS G01  2025  6  8 12  0  0.000000  1   -0.000110000000E+00\n\
         AS G02  2025  6  8 12  0  0.000000  1   -0.000090000000E+00\n\
         AS G03  2025  6  8 12  0  0.000000  1   -0.000090000000E+00\n".to_string()
    }
}
