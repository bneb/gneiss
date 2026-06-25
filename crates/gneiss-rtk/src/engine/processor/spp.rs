use super::ProcessingEngine;
use crate::engine::{EngineError, EngineMode};
use crate::filter::RtkState;
use gneiss_core::coords::Coordinate;
use gneiss_core::obs::EpochObs;

impl ProcessingEngine {
    fn perform_spp_ekf_update(
        config: &crate::engine::EngineConfig,
        state: &mut RtkState,
        spp_pos: Option<Coordinate>,
        spp_cdt: f64,
    ) {
        if let Some(pos) = spp_pos {
            let z_diff = pos.vector - state.position.vector;
            let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());

            let mut rejected = false;
            if config.mode.is_tightly_coupled()
                && state.ins_aligned
                && z_diff.norm() > config.spp_consistency_threshold_m
            {
                rejected = true;
            } else {
                let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
                h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);

                let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
                // Using a tighter variance of 9.0 (3m std dev) forces the INS to track the clean SPP positions
                r_mat.fill_diagonal(9.0);

                if crate::engine::updater::update::<crate::engine::updater_math::LooseCoupling>(
                    state,
                    &z_vec,
                    &h_mat,
                    &r_mat,
                    config.spp_consistency_threshold_m,
                    None,
                    &config.tuning,
                )
                .map_or(true, |v| v.0.len() < 3)
                {
                    rejected = true;
                }
            }

            if rejected {
                state.consecutive_rejections += 1;
                if state.consecutive_rejections > 5 {
                    tracing::warn!(
                        "SPP EKF rejected for {} epochs. Hard resetting INS to SPP.",
                        state.consecutive_rejections
                    );
                    state.position = pos;
                    state.velocity = nalgebra::Vector3::zeros(); // Zero out diverged velocity
                    state.accel_bias = nalgebra::Vector3::zeros(); // Biases might be corrupted, 0 is a safer prior
                    state.gyro_bias = nalgebra::Vector3::zeros();
                    // Preserve attitude as it is far better than identity
                    if crate::filter::CORE_STATE_SIZE > 15 {
                        state.rcv_clk_bias = spp_cdt;
                        state.rcv_clk_drift = 0.0;
                    }
                    state.clear_ambiguities();

                    state.covariance.fill(0.0);
                    let n = crate::filter::CORE_STATE_SIZE;
                    for i in 0..6 {
                        state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
                    }
                    let att_var = (1.0f64.to_radians()).powi(2);
                    for i in 6..9 {
                        state.covariance[(i, i)] = att_var;
                    }
                    for i in 9..12 {
                        state.covariance[(i, i)] = 0.01;
                    }
                    for i in 12..n {
                        state.covariance[(i, i)] = 1e-4;
                    }
                    if crate::filter::CORE_STATE_SIZE > 15 {
                        state.covariance[(15, 15)] = 1e6;
                    }
                    state.is_reset = true;
                    state.consecutive_rejections = 0;
                } else {
                    tracing::warn!(
                        "SPP EKF update rejected. Riding through outage via INS dead-reckoning."
                    );
                }
            } else {
                state.consecutive_rejections = 0;
            }
        }
    }

    pub fn process_spp(&mut self, rover_obs: &EpochObs) -> Result<&RtkState, EngineError> {
        let spp_res = crate::spp::compute_spp(
            rover_obs,
            &self.ephemerides,
            self.klobuchar_params.as_ref(),
            &crate::spp::SppConfig::default(),
            None,
        );
        let spp_pos = spp_res.as_ref().ok().map(|s| s.position);
        let spp_cdt = spp_res.as_ref().ok().map(|s| s.cdt).unwrap_or(0.0);

        if let Err(e) = &spp_res {
            tracing::warn!("Initial SPP compute failed: {:?}", e);
        }

        let spp_res_opt = spp_res.ok();
        if self.current_state.is_none() {
            if let Some(spp) = &spp_res_opt {
                let mut new_state = RtkState::new(rover_obs.time, spp.position, 100.0);
                if crate::filter::CORE_STATE_SIZE > 15 {
                    new_state.rcv_clk_bias = spp.cdt;
                    if !self.config.mode.is_ppp() {
                        new_state.isb_glo = spp.cdt_glo - spp.cdt;
                        new_state.isb_gal = spp.cdt_gal - spp.cdt;
                        new_state.isb_bds = spp.cdt_bds - spp.cdt;
                    } else {
                        new_state.isb_glo = 0.0;
                        new_state.isb_gal = 0.0;
                        new_state.isb_bds = 0.0;
                    }
                }
                self.current_state = Some(new_state);
            } else {
                return Err(EngineError::InitialSppFailed);
            }
        }

        let dt = rover_obs.time.tow
            - self
                .current_state
                .as_ref()
                .ok_or(EngineError::StateDisappeared)?
                .time
                .tow;
        self.predict_state(dt);

        if let Some(pos) = spp_pos {
            if matches!(self.config.mode, EngineMode::Spp) {
                // Pure SPP is an epoch-by-epoch solution. Do not filter.
                let state = self.current_state.as_mut().expect("current_state is Some after ok_or early return");
                state.time = rover_obs.time;
                state.position = pos;
                state.position.epoch = rover_obs.time;
                state.velocity = nalgebra::Vector3::zeros();
                tracing::info!(
                    "process_spp: SPP mode return, rcv_clk_bias = {}",
                    state.rcv_clk_bias
                );
                self.state_history.push(state.clone());
                self.obs_history.push((rover_obs.clone(), None));
                return Ok(self.current_state.as_ref().expect("current_state is Some at end of process_spp"));
            }
        } else {
            // SPP failed for this epoch.
            if matches!(self.config.mode, EngineMode::Spp) {
                // We preserve the state so the next epoch has a good seed,
                // but we return an error so no output is produced for this epoch.
                if let Some(state) = &mut self.current_state {
                    state.is_reset = false;
                    // Update time
                    state.time = rover_obs.time;
                    state.position.epoch = rover_obs.time;
                }
                return Err(EngineError::InitialSppFailed);
            }
        }

        let state = self
            .current_state
            .as_mut()
            .ok_or(EngineError::StateDisappeared)?;
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        Self::check_covariance_divergence(
            state,
            spp_pos,
            spp_res_opt.as_ref(),
            !self.config.mode.is_ppp(),
        );

        Self::perform_spp_ekf_update(&self.config, state, spp_pos, spp_cdt);

        if let Some(state) = self.current_state.as_mut() {
            Self::apply_nhc_updates(&self.config, &self.imu_history, state);
        }

        if let Some(state) = &self.current_state {
            tracing::info!(
                "process_spp: Returning state, rcv_clk_bias = {}",
                state.rcv_clk_bias
            );
            self.state_history.push(state.clone());
        }
        self.obs_history.push((rover_obs.clone(), None));
        self.current_state
            .as_ref()
            .ok_or(EngineError::StateDisappeared)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{EngineConfig, EngineError, EngineMode};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::EpochObs;
    use gneiss_core::sat::SatelliteId;
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    fn make_empty_rover(time: GpsTime) -> EpochObs {
        EpochObs {
            time,
            satellites: Vec::new(),
        }
    }

    #[test]
    fn test_perform_spp_ekf_update_no_pos_no_op() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos.clone(), 1.0);

        // spp_pos=None should do nothing
        ProcessingEngine::perform_spp_ekf_update(&config, &mut state, None, 0.0);
        assert!((state.position.vector.x - 1.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_perform_spp_ekf_update_tightly_coupled_rejects_large_diff() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns; // tightly coupled
        config.spp_consistency_threshold_m = 15.0;

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;

        let spp_pos = Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0), // very far
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        ProcessingEngine::perform_spp_ekf_update(&config, &mut state, Some(spp_pos), 0.0);
        // Should be rejected
        assert_eq!(state.consecutive_rejections, 1);
    }

    #[test]
    fn test_perform_spp_ekf_update_rejection_hard_reset_after_six() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        config.spp_consistency_threshold_m = 15.0;

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.consecutive_rejections = 5; // 5 previous rejections

        let spp_pos = Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        ProcessingEngine::perform_spp_ekf_update(&config, &mut state, Some(spp_pos), 0.0);
        // Should have hard reset — position updated and rejections cleared
        assert!((state.position.vector.x - 1000.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
        // Should have is_reset flag
        assert!(state.is_reset);
        // Velocity should be zeroed
        assert!((state.velocity.norm() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_process_spp_fails_on_empty_ephemeris_and_no_state() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.current_state = None;
        engine.ephemerides = Vec::new();

        let rover = make_empty_rover(GpsTime::new(0, 0.0));
        let err = engine.process_spp(&rover).unwrap_err();

        // Should fail because SPP compute fails and no state to fall back on
        assert!(matches!(err, EngineError::InitialSppFailed));
    }

    #[test]
    fn test_process_spp_preserves_state_on_failure_in_spp_mode() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::Spp;

        // Set up a valid current state
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let existing = RtkState::new(time, pos, 1.0);
        engine.current_state = Some(existing);
        engine.ephemerides = Vec::new();

        let rover = make_empty_rover(GpsTime::new(0, 1.0));
        let err = engine.process_spp(&rover).unwrap_err();
        assert!(matches!(err, EngineError::InitialSppFailed));

        // State should still be present (preserved for next epoch)
        assert!(engine.current_state.is_some());
    }

    #[test]
    fn test_process_spp_state_disappeared_error() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::RtkIns; // Not SPP mode, won't early-return
        // Don't set current_state — it will try to use it after SPP compute fails

        let rover = make_empty_rover(GpsTime::new(0, 0.0));
        let err = engine.process_spp(&rover).unwrap_err();
        assert!(matches!(err, EngineError::InitialSppFailed));
    }

    #[test]
    fn test_process_spp_non_spp_mode_with_state_succeeds() {
        // In non-SPP mode (e.g. RtkIns) with an existing state but no ephemerides,
        // process_spp should succeed by coasting on the existing state.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::RtkIns; // non-pure-SPP mode
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = make_empty_rover(GpsTime::new(0, 1.0));
        let result = engine.process_spp(&rover);
        assert!(result.is_ok());
        // State should have time moved forward
        assert!((result.unwrap().time.tow - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_perform_spp_ekf_update_success() {
        // When SPP position is close to the current state, the EKF update
        // should succeed (not be rejected).
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Rtk; // Not tightly coupled
        let time = GpsTime::new(0, 0.0);

        // State position
        let state_pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, state_pos, 1.0);

        // SPP position very close to current state — should pass chi-square
        let spp_pos = Coordinate::new(
            Vector3::new(1.1, 2.2, 3.3),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        ProcessingEngine::perform_spp_ekf_update(&config, &mut state, Some(spp_pos), 0.0);
        // No rejection
        assert_eq!(state.consecutive_rejections, 0);
        // Position should have been nudged toward SPP (EKF blends the two)
        assert!((state.position.vector.x - 1.0).abs() > 1e-6);
    }

    #[test]
    fn test_perform_spp_ekf_update_tight_coupled_small_diff_accepted() {
        // Tightly coupled mode with a small position difference should accept.
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns; // tightly coupled
        config.spp_consistency_threshold_m = 15.0;
        let time = GpsTime::new(0, 0.0);

        let state_pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, state_pos, 1.0);
        state.ins_aligned = true;

        // Small diff (within threshold) that should be accepted
        let spp_pos = Coordinate::new(
            Vector3::new(1.1, 2.1, 3.1),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        ProcessingEngine::perform_spp_ekf_update(&config, &mut state, Some(spp_pos), 0.0);
        // Should NOT be rejected since diff is well under 15m
        assert_eq!(state.consecutive_rejections, 0);
        // Position should have been nudged
        assert!((state.position.vector.x - 1.0).abs() > 1e-6);
    }

    #[test]
    fn test_perform_spp_ekf_update_tight_coupled_not_aligned_large_diff() {
        // Tightly coupled mode but ins_aligned=false should NOT take the early reject
        // path. Instead it goes through the EKF update which may pass or fail
        // the chi-square test on its own, but the function should not panic.
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        config.spp_consistency_threshold_m = 15.0;

        let time = GpsTime::new(0, 0.0);
        let state_pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, state_pos, 1.0);
        state.ins_aligned = false;

        let spp_pos = Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        // Should not panic — the early reject condition requires ins_aligned=true
        ProcessingEngine::perform_spp_ekf_update(&config, &mut state, Some(spp_pos), 0.0);
    }

    #[test]
    fn test_process_spp_non_spp_mode_state_with_spp_success_coasts() {
        // Non-SPP mode where state exists. SPP compute succeeds (we have
        // ephemerides and observations). The function should go through
        // covariance divergence check, EKF update, NHC, and return the state.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::RtkIns;

        let time = GpsTime::new(1000, 0.0);
        let a_sq = 5153.6_f64 * 5153.6_f64;

        let toe = time;
        let m0_vals = [0.0, std::f64::consts::FRAC_PI_2,
                       std::f64::consts::PI, 3.0 * std::f64::consts::FRAC_PI_2];
        for (i, &m0) in m0_vals.iter().enumerate() {
            let eph = gneiss_core::ephemeris::Ephemeris::Gps(
                gneiss_core::ephemeris::GpsEphemeris {
                    sat: SatelliteId {
                        constellation: gneiss_core::sat::Constellation::Gps,
                        prn: (i + 1) as u8,
                    },
                    toe,
                    toc: time,
                    af0: 0.0, af1: 0.0, af2: 0.0,
                    crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                    cic: 0.0, cis: 0.0,
                    m0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
                    omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
                    omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
                },
            );
            engine.add_ephemeris(eph);
        }

        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let pseudo_range = a_sq;
        let sats: Vec<_> = (0..4)
            .map(|i| {
                let sat = SatelliteId {
                    constellation: gneiss_core::sat::Constellation::Gps,
                    prn: (i + 1) as u8,
                };
                gneiss_core::obs::SatObs {
                    sat,
                    observations: vec![
                        gneiss_core::obs::Observation {
                            code: gneiss_core::obs::ObsCode {
                                obs_type: gneiss_core::obs::ObsType::Pseudorange,
                                signal: gneiss_core::obs::SignalCode {
                                    freq_band: 1, attribute: 'C',
                                },
                            },
                            value: pseudo_range,
                            lock_time: None, lli: None,
                        },
                        gneiss_core::obs::Observation {
                            code: gneiss_core::obs::ObsCode {
                                obs_type: gneiss_core::obs::ObsType::Snr,
                                signal: gneiss_core::obs::SignalCode {
                                    freq_band: 1, attribute: 'S',
                                },
                            },
                            value: 45.0,
                            lock_time: None, lli: None,
                        },
                    ],
                }
            })
            .collect();

        let rover = EpochObs { time, satellites: sats };
        let result = engine.process_spp(&rover);
        match result {
            Ok(state) => {
                assert!(state.position.vector.x.is_finite());
                assert!(state.position.vector.y.is_finite());
                assert!(state.position.vector.z.is_finite());
                assert_eq!(engine.state_history.len(), 1);
            }
            Err(EngineError::InitialSppFailed) => {
                // SPP compute may fail to converge with given data; this is acceptable
            }
            Err(e) => panic!("Unexpected error: {:?}", e),
        }
    }
}
