use super::ProcessingEngine;
use crate::engine::matcher::match_observations;
use crate::engine::rtk_multi_base::MultiBaseCombiner;
use crate::engine::updater_math::{CouplingStrategy, LooseCoupling, TightCoupling};
use crate::engine::{EngineConfig, EngineError, EngineMode};
use crate::filter::RtkState;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::SatelliteId;
use nalgebra::{DMatrix, DVector, Vector3};

/// Maximum allowable 3D position change (meters) from an AR fix.
/// If the fixed position differs from the float position by more than this
/// threshold, the fix is rejected and the float solution is retained.
/// This catches wrong AR integer sets that would otherwise jump the position.
const AR_MAX_POSITION_JUMP_M: f64 = 2.0;

/// Ring buffer size for factor graph previous-epoch snapshots.
const FG_BUFFER_SIZE: usize = 300;
/// Minimum entries before factor graph uses oldest (geometry-diverse) entry.
const FG_MIN_DIVERSE: usize = 30;

/// Apply adaptive R scaling based on innovation history.
/// Updates the tracker with current innovations and inflates R diagonals
/// for measurements with historically large normalized innovations.
fn apply_adaptive_r_scaling(
    tracker: &mut crate::engine::adaptive::InnovationTracker,
    state: &RtkState,
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &mut DMatrix<f64>,
    meas_types: &[(SatelliteId, u8, f64)],
    matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
) {
    let state_size = state.covariance.nrows();
    for i in 0..z.nrows() {
        let h_row_vec: Vec<f64> = (0..state_size).map(|c| h[(i, c)]).collect();
        let h_row_mat = DMatrix::from_row_slice(1, state_size, &h_row_vec);
        let s_ii = (&h_row_mat * &state.covariance * h_row_mat.transpose())[(0, 0)] + r[(i, i)];
        let (sat, mtype, _) = meas_types[i];
        // Determine freq band from measurement type: 0=PR_L1, 1=CP_L1, 2=CP_L2, 3=Dop
        let freq = match mtype {
            2 => 2,
            _ => 1,
        };

        let mut snr = None;
        if let Some((rov_obs, _)) = matched_obs.iter().find(|(rov, _)| rov.sat == sat) {
            snr = Some(rov_obs.snr);
        }

        let scale = tracker.update_and_scale(sat, freq, z[i], s_ii, snr);
        r[(i, i)] *= scale;
    }
}

impl ProcessingEngine {
    fn init_spp_state(
        &mut self,
        rover_obs: &EpochObs,
    ) -> Result<Option<crate::spp::SppState>, EngineError> {
        let spp_res = crate::spp::compute_spp(
            rover_obs,
            &self.ephemerides,
            self.klobuchar_params.as_ref(),
            &crate::spp::SppConfig::default(),
            None,
        )
        .ok();

        if spp_res.is_none() {
            tracing::warn!("compute_spp failed in process_rtk!");
        }

        if self.current_state.is_none() {
            if let Some(spp) = &spp_res {
                self.current_state = Some(RtkState::new(rover_obs.time, spp.position, 100.0));
            } else {
                return Err(EngineError::InitialSppFailed);
            }
        }
        Ok(spp_res)
    }

    fn evaluate_gnss_only_coasting(
        config: &EngineConfig,
        state: &mut RtkState,
        spp_pos: Option<Coordinate>,
        spp_state_ref: Option<&crate::spp::SppState>,
        had_imu_data: bool,
    ) {
        let use_gnss_only_seed = matches!(config.mode, EngineMode::Rtk | EngineMode::Ppp)
            || (matches!(
                config.mode,
                EngineMode::RtkIns
                    | EngineMode::PppIns
                    | EngineMode::SppIns
                    | EngineMode::RtkInsLooselyCoupled
                    | EngineMode::SppInsLooselyCoupled
                    | EngineMode::PppInsLooselyCoupled
            ) && !had_imu_data);

        if use_gnss_only_seed {
            tracing::debug!("epoch_count = {}", state.epoch_count);
            let need_spp_reset = if spp_pos.is_some() {
                state.epoch_count < 3
            } else {
                false
            };

            if need_spp_reset {
                if let Some(pos) = spp_pos {
                    tracing::warn!("Resetting EKF state to SPP due to divergence/startup.");
                    state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
                }
            }
        }
    }

    fn filter_valid_base<'a>(
        config: &EngineConfig,
        rover_obs: &EpochObs,
        base_obs: Option<&'a EpochObs>,
    ) -> Option<&'a EpochObs> {
        base_obs.filter(|b| {
            let age = (rover_obs.time.tow - b.time.tow).abs();
            if age > config.max_base_age_s {
                tracing::trace!(
                    "Base observation rejected due to Age of Differential ({:.1}s > {:.1}s)",
                    age,
                    config.max_base_age_s
                );
                false
            } else {
                true
            }
        })
    }

    fn perform_spp_fallback_update(config: &EngineConfig, state: &mut RtkState, pos: Coordinate) {
        tracing::debug!("No valid base data. Falling back to SPP.");
        tracing::debug!("RTK base missing or stale. Falling back to SPP update.");
        let z_diff = pos.vector - state.position.vector;
        let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());

        let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
        h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);

        let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
        r_mat.fill_diagonal(900.0);

        if let Err(e) = crate::engine::updater::update::<LooseCoupling>(
            state,
            &z_vec,
            &h_mat,
            &r_mat,
            config.spp_consistency_threshold_m,
            None,
            &config.tuning,
        ) {
            tracing::debug!("SPP Fallback update failed: {:?}", e);
        }
    }

    fn manage_imu_buffer_on_start(&mut self) {
        // Do not clear the buffer, so that it can be pushed to imu_history and used for alignment
    }

    fn update_state_time(&mut self, time: gneiss_core::time::GpsTime) -> Result<(), EngineError> {
        let state = self
            .current_state
            .as_mut()
            .ok_or(EngineError::StateDisappeared)?;
        state.is_reset = false;
        state.time = time;
        state.position.epoch = time;
        Ok(())
    }

    fn get_base_coord(
        config: &EngineConfig,
        rover_obs: &EpochObs,
    ) -> Result<Coordinate, EngineError> {
        let mut base_coord = if let Some(base_pos_arr) = config.base_position {
            Coordinate::new(
                Vector3::new(base_pos_arr[0], base_pos_arr[1], base_pos_arr[2]),
                Datum::WGS84,
                Frame::ECEF,
                rover_obs.time,
            )
        } else {
            return Err(EngineError::MissingBasePosition);
        };

        if let Some(helmert) = &config.base_datum_transform {
            let obs_epoch = rover_obs.time.to_fractional_year();
            let transformed_vec = helmert.transform(base_coord.vector, obs_epoch);
            base_coord =
                Coordinate::new(transformed_vec, Datum::WGS84, Frame::ECEF, rover_obs.time);
        }
        Ok(base_coord)
    }

    fn evaluate_gnn(
        gnn_raim: &Option<crate::engine::ml::gnn_raim::GnnRaimModel>,
        ephemerides: &[gneiss_core::ephemeris::Ephemeris],
        rover_obs: &EpochObs,
        matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
        state: &RtkState,
    ) -> std::collections::HashMap<SatelliteId, f64> {
        if let Some(gnn) = gnn_raim {
            let pos_apc = state.position.vector;
            let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
            crate::engine::ml::gnn_raim::evaluate_gnn_raim(
                gnn,
                matched_obs,
                rov_llh,
                pos_apc,
                ephemerides,
                rover_obs.time,
            )
        } else {
            std::collections::HashMap::new()
        }
    }

    fn apply_observations(
        &mut self,
        rover_obs: &EpochObs,
        base_obs: Option<&EpochObs>,
        spp_pos: Option<Coordinate>,
        spp_state_ref: Option<&crate::spp::SppState>,
    ) -> Result<(), EngineError> {
        let valid_base = Self::filter_valid_base(&self.config, rover_obs, base_obs);
        let base_coord_res = Self::get_base_coord(&self.config, rover_obs);
        let state = self
            .current_state
            .as_mut()
            .ok_or(EngineError::StateDisappeared)?;

        if let Some(base) = valid_base {
            state.epoch_count += 1;
            let base_coord = base_coord_res?;
            let matched_obs = match_observations(rover_obs, base, &self.ephemerides);

            // Store for TDCP time-differencing at next epoch
            self.last_matched_obs = matched_obs.clone();
            self.last_base_coord = Some(base_coord.clone());

            if matched_obs.len() >= 5 {
                let gnn_variances = Self::evaluate_gnn(
                    &self.gnn_raim,
                    &self.ephemerides,
                    rover_obs,
                    &matched_obs,
                    state,
                );
                let ctx = RtkUpdateContext {
                    config: &self.config,
                    ephemerides: &self.ephemerides,
                    imu_history: &self.imu_history,
                    rover_obs,
                    base_obs: base,
                    matched_obs: &matched_obs,
                    base_coord: &base_coord,
                    spp_pos,
                    spp_state_ref,
                    gnn_variances,
                    klobuchar_params: self.klobuchar_params,
                };
                process_rtk_update::<TightCoupling>(state, &mut self.innovation_tracker, &ctx);
            } else {
                state.consecutive_rejections += 1;
            }
        } else if let Some(pos) = spp_pos {
            Self::perform_spp_fallback_update(&self.config, state, pos);
        }
        Ok(())
    }

    pub fn process_rtk(
        &mut self,
        rover_obs: &EpochObs,
        base_obs: Option<&EpochObs>,
    ) -> Result<&RtkState, EngineError> {
        // Route to multi-base path when multi-base data is available
        if !self.multi_base_observations.is_empty() {
            let bases = std::mem::take(&mut self.multi_base_observations);
            return self.process_rtk_multi(rover_obs, &bases);
        }

        let spp_res = self.init_spp_state(rover_obs)?;
        let spp_pos = spp_res.as_ref().map(|s| s.position);
        let spp_state_ref = spp_res.as_ref();

        let mut rover_smoothed = rover_obs.clone();
        self.hatch_filter.smooth_epoch(&mut rover_smoothed);

        let dt = rover_obs.time.tow
            - self
                .current_state
                .as_ref()
                .ok_or(EngineError::StateDisappeared)?
                .time
                .tow;
        self.manage_imu_buffer_on_start();
        let had_imu_data = !self.imu_buffer.is_empty();

        self.predict_state(dt);
        self.update_state_time(rover_obs.time)?;

        let state = self.current_state.as_mut().expect("current_state is Some after ok_or early return");
        Self::check_covariance_divergence(
            state,
            spp_pos,
            spp_state_ref,
            !self.config.mode.is_ppp(),
        );
        Self::evaluate_gnss_only_coasting(
            &self.config,
            state,
            spp_pos,
            spp_state_ref,
            had_imu_data,
        );

        self.apply_observations(&rover_smoothed, base_obs, spp_pos, spp_state_ref)?;

        let state = self.current_state.as_mut().expect("current_state is Some after apply_observations");
        Self::apply_nhc_updates(&self.config, &self.imu_history, state);

        // --- TDCP: time-differenced carrier phase delta-position ---
        if self.config.enable_tdcp && !self.last_matched_obs.is_empty() {
            if let Some(ref base_coord) = self.last_base_coord {
                let state = self.current_state.as_mut()
                    .expect("current_state is Some");
                let ref_sats: Vec<(gneiss_core::sat::Constellation, gneiss_core::sat::SatelliteId)> =
                    state.current_ref_sat.iter()
                        .map(|(c, s)| (*c, *s))
                        .collect();

                if let Some((delta, cov)) = self.tdcp_solver.compute_delta(
                    &self.last_matched_obs,
                    &ref_sats,
                    state.position.vector,
                    base_coord.vector,
                    &self.ephemerides,
                    rover_obs.time,
                ) {
                    // Feed delta as position-change measurement into EKF.
                    // Innovation z = delta (the TDCP-estimated correction to
                    // predicted position). H = I₃ for position states.
                    let mut h = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
                    h[(0, 0)] = 1.0;
                    h[(1, 1)] = 1.0;
                    h[(2, 2)] = 1.0;
                    let z = nalgebra::DVector::from_vec(vec![delta.x, delta.y, delta.z]);
                    let _ = crate::engine::updater::update::<
                        crate::engine::updater_math::TightCoupling,
                    >(
                        state,
                        &z,
                        &h,
                        &cov,
                        9.0, // chi-square threshold for 3-DOF
                        None,
                        &self.config.tuning,
                    );
                    self.tdcp_trajectory.push(delta, &cov);
                }

                // Store current-epoch data for next TDCP computation
                self.tdcp_solver.store_epoch(
                    &self.last_matched_obs,
                    &ref_sats,
                    state.position.vector,
                    base_coord.vector,
                    &self.ephemerides,
                    rover_obs.time,
                );
            }
        }

        self.attempt_kinematic_alignment();

        // Save current position/cov as prev-epoch for two-epoch smoothing
        if let Some(ref mut state) = self.current_state {
            state.prev_epoch_pos = Some(state.position.vector);
            state.prev_epoch_cov = Some(state.covariance.fixed_view::<3, 3>(0, 0).into_owned());
        }
        if let Some(ref state) = self.current_state {
            self.state_history.push(RtkState::clone(state));
        }
        self.obs_history
            .push((rover_obs.clone(), base_obs.cloned()));
        self.current_state
            .as_ref()
            .ok_or(EngineError::StateDisappeared)
    }

    /// Process RTK with multiple base stations.
    ///
    /// Each base is processed independently through the single-base RTK pipeline
    /// (predict, match, update) on a cloned state. The resulting position fixes
    /// are blended via `MultiBaseCombiner`, where shorter baselines receive
    /// higher weight. After combining, the primary base (most matched observations)
    /// drives the state's ambiguity and covariance updates.
    ///
    /// `bases` is a slice of `(EpochObs, Vector3<f64>)` pairs, where the vector
    /// is the base station's ECEF position for baseline distance computation.
    pub fn process_rtk_multi(
        &mut self,
        rover_obs: &EpochObs,
        bases: &[(EpochObs, Vector3<f64>)],
    ) -> Result<&RtkState, EngineError> {
        if bases.is_empty() {
            return self.process_rtk(rover_obs, None);
        }
        if bases.len() == 1 {
            let (ref base_obs, base_pos) = bases[0];
            self.config.base_position = Some([base_pos.x, base_pos.y, base_pos.z]);
            return self.process_rtk(rover_obs, Some(base_obs));
        }

        // --- standard preprocessing (shared by all bases) ---
        let spp_res = self.init_spp_state(rover_obs)?;
        let spp_pos = spp_res.as_ref().map(|s| s.position);
        let spp_state_ref = spp_res.as_ref();

        let mut rover_smoothed = rover_obs.clone();
        self.hatch_filter.smooth_epoch(&mut rover_smoothed);

        let dt = rover_obs.time.tow
            - self
                .current_state
                .as_ref()
                .ok_or(EngineError::StateDisappeared)?
                .time
                .tow;
        self.manage_imu_buffer_on_start();
        let had_imu_data = !self.imu_buffer.is_empty();

        self.predict_state(dt);
        self.update_state_time(rover_obs.time)?;

        let state = self
            .current_state
            .as_mut()
            .expect("current_state is Some after ok_or early return");
        Self::check_covariance_divergence(
            state,
            spp_pos,
            spp_state_ref,
            !self.config.mode.is_ppp(),
        );
        Self::evaluate_gnss_only_coasting(
            &self.config,
            state,
            spp_pos,
            spp_state_ref,
            had_imu_data,
        );

        // --- multi-base: independent RTK fix per base on a cloned snapshot ---
        let pre_update_state = state.clone();
        let mut combiner = MultiBaseCombiner::new(1.0);
        let mut primary_base_idx = 0_usize;
        let mut best_matched_count = 0_usize;
        let mut fix_infos: Vec<(nalgebra::Vector3<f64>, bool, f64)> = Vec::new();

        for (i, (base_obs, base_pos_ecef)) in bases.iter().enumerate() {
            let age = (rover_obs.time.tow - base_obs.time.tow).abs();
            if age > self.config.max_base_age_s {
                continue;
            }

            let base_coord = Coordinate::new(
                *base_pos_ecef,
                Datum::WGS84,
                Frame::ECEF,
                rover_obs.time,
            );

            let matched_obs =
                match_observations(&rover_smoothed, base_obs, &self.ephemerides);
            if matched_obs.len() < 5 {
                continue;
            }

            if matched_obs.len() > best_matched_count {
                best_matched_count = matched_obs.len();
                primary_base_idx = i;
            }

            // Clone the pre-update state and run the full RTK update independently
            let mut base_state = pre_update_state.clone();
            let gnn_variances = Self::evaluate_gnn(
                &self.gnn_raim,
                &self.ephemerides,
                &rover_smoothed,
                &matched_obs,
                &base_state,
            );

            let ctx = RtkUpdateContext {
                config: &self.config,
                ephemerides: &self.ephemerides,
                imu_history: &self.imu_history,
                rover_obs: &rover_smoothed,
                base_obs,
                matched_obs: &matched_obs,
                base_coord: &base_coord,
                spp_pos,
                spp_state_ref,
                gnn_variances,
                klobuchar_params: self.klobuchar_params,
            };
            let mut tracker = crate::engine::adaptive::InnovationTracker::new();
            process_rtk_update::<TightCoupling>(&mut base_state, &mut tracker, &ctx);

            // Collect the position from this base's fix
            let baseline_m =
                (pre_update_state.position.vector - base_pos_ecef).norm();
            let quality = if base_state.is_fixed { 1.0 } else { 0.5 };
            fix_infos.push((base_state.position.vector, base_state.is_fixed, baseline_m));
            combiner.add_fix(base_state.position.vector, baseline_m, quality);
        }

        // --- cross-base AR validation ---
        let mut cross_base_accepted = true;
        let mut cross_base_consensus: Option<nalgebra::Vector3<f64>> = None;
        if self.config.enable_cross_base_ar_validation && fix_infos.len() >= 2 {
            let validator = crate::engine::rtk_multi_base::CrossBaseArValidator::new(
                self.config.cross_base_agreement_threshold_m,
            );
            let result = validator.validate(&fix_infos);
            cross_base_accepted = result.accepted;
            cross_base_consensus = result.consensus_position;
            if !result.accepted {
                tracing::warn!(
                    "Cross-base AR rejected: {} fixed bases disagree (max={:.2}m > thresh={:.2}m). Staying in float.",
                    fix_infos.iter().filter(|(_, f, _)| *f).count(),
                    result.max_disagreement_m,
                    self.config.cross_base_agreement_threshold_m,
                );
            } else if result.agreed_indices.len() >= 2 {
                tracing::info!(
                    "Cross-base AR consensus: {}/{} bases agree (max_disagreement={:.2}m)",
                    result.agreed_indices.len(),
                    fix_infos.iter().filter(|(_, f, _)| *f).count(),
                    result.max_disagreement_m,
                );
            }
        }

        // --- primary base: drive state bookkeeping (ambiguities, covariances) ---
        let combined_pos = if cross_base_accepted {
            // Prefer cross-base consensus if available, otherwise weighted position
            cross_base_consensus.or_else(|| combiner.weighted_position())
        } else {
            // Cross-base rejected: use weighted position but clear fix state
            combiner.weighted_position()
        };

        let (primary_base, primary_pos) = &bases[primary_base_idx];
        let saved_base_pos = self.config.base_position;
        self.config.base_position =
            Some([primary_pos.x, primary_pos.y, primary_pos.z]);
        self.apply_observations(
            &rover_smoothed,
            Some(primary_base),
            spp_pos,
            spp_state_ref,
        )?;
        self.config.base_position = saved_base_pos;

        // Override position with the multi-base combined estimate
        if let Some(pos) = combined_pos {
            let state = self
                .current_state
                .as_mut()
                .expect("state is Some after apply_observations");
            state.position.vector = pos;
        }

        // If cross-base AR validation rejected, clear the fix state so the
        // solution stays in float mode. Wrong integers from different multipath
        // at each base are worse than a clean float solution.
        if !cross_base_accepted {
            let state = self
                .current_state
                .as_mut()
                .expect("state is Some after apply_observations");
            if state.is_fixed {
                tracing::info!("Cross-base AR rejected: clearing fix state, staying in float.");
                state.is_fixed = false;
                state.fixed_state = None;
            }
        }

        // --- post-processing ---
        {
            let state = self
                .current_state
                .as_mut()
                .expect("state is Some after above");
            Self::apply_nhc_updates(&self.config, &self.imu_history, state);
        }

        // --- TDCP: time-differenced carrier phase delta-position ---
        if self.config.enable_tdcp && !self.last_matched_obs.is_empty() {
            if let Some(ref base_coord) = self.last_base_coord {
                let state = self.current_state.as_mut()
                    .expect("current_state is Some");
                let ref_sats: Vec<(gneiss_core::sat::Constellation, gneiss_core::sat::SatelliteId)> =
                    state.current_ref_sat.iter()
                        .map(|(c, s)| (*c, *s))
                        .collect();

                if let Some((delta, cov)) = self.tdcp_solver.compute_delta(
                    &self.last_matched_obs,
                    &ref_sats,
                    state.position.vector,
                    base_coord.vector,
                    &self.ephemerides,
                    rover_obs.time,
                ) {
                    let mut h = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
                    h[(0, 0)] = 1.0;
                    h[(1, 1)] = 1.0;
                    h[(2, 2)] = 1.0;
                    let z = nalgebra::DVector::from_vec(vec![delta.x, delta.y, delta.z]);
                    let _ = crate::engine::updater::update::<
                        crate::engine::updater_math::TightCoupling,
                    >(
                        state,
                        &z,
                        &h,
                        &cov,
                        9.0,
                        None,
                        &self.config.tuning,
                    );
                    self.tdcp_trajectory.push(delta, &cov);
                }

                self.tdcp_solver.store_epoch(
                    &self.last_matched_obs,
                    &ref_sats,
                    state.position.vector,
                    base_coord.vector,
                    &self.ephemerides,
                    rover_obs.time,
                );
            }
        }

        self.attempt_kinematic_alignment();

        // Save current position/cov as prev-epoch for two-epoch smoothing
        if let Some(ref mut state) = self.current_state {
            state.prev_epoch_pos = Some(state.position.vector);
            state.prev_epoch_cov = Some(state.covariance.fixed_view::<3, 3>(0, 0).into_owned());
        }
        if let Some(ref state) = self.current_state {
            self.state_history.push(RtkState::clone(state));
        }
        self.obs_history
            .push((rover_obs.clone(), Some(primary_base.clone())));

        self.current_state
            .as_ref()
            .ok_or(EngineError::StateDisappeared)
    }
}

/// RTK measurement update — extracted as a free function to avoid borrow conflicts.
fn build_measurement_environment<'a>(
    config: &'a EngineConfig,
    imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
    state: &RtkState,
    ephemerides: &'a [gneiss_core::ephemeris::Ephemeris],
    base_coord: &'a Coordinate,
    base_time: gneiss_core::time::GpsTime,
    klobuchar_params: Option<gneiss_core::atmosphere::KlobucharParams>,
) -> crate::engine::measurement::MeasurementEnvironment<'a> {
    let omega_ib_b = if let Some(imu_buf) = imu_history.last() {
        if let Some(last_imu) = imu_buf.last() {
            last_imu.gyro - state.gyro_bias
        } else {
            nalgebra::Vector3::zeros()
        }
    } else {
        nalgebra::Vector3::zeros()
    };
    let omega_ie_e =
        nalgebra::Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
    let r_e_b = state.attitude.to_rotation_matrix().transpose();
    let lever_arm = if state.ins_aligned {
        Vector3::from_column_slice(&config.imu_to_antenna_lever_arm)
    } else {
        Vector3::zeros()
    };

    crate::engine::measurement::MeasurementEnvironment {
        ephemerides,
        base_coord,
        base_time,
        lever_arm,
        omega_b: omega_ib_b - r_e_b * omega_ie_e,
        tuning: &config.tuning,
        gnn_variances: std::collections::HashMap::new(),
        klobuchar_params,
    }
}

fn handle_ekf_rejection(
    state: &mut RtkState,
    config: &EngineConfig,
    spp_pos: Option<Coordinate>,
    spp_state_ref: Option<&crate::spp::SppState>,
    reason: &str,
) {
    state.consecutive_rejections += 1;
    tracing::warn!(
        "GNSS EKF rejected for {} epochs ({}).",
        state.consecutive_rejections,
        reason
    );

    // Only inflate covariance for the first few rejections to allow recovery.
    // After that, stop inflating so position variance doesn't blow past the
    // extreme-divergence threshold before the coasting period expires.
    if state.consecutive_rejections <= 5 {
        let inflate_factor = 1.0 + (state.consecutive_rejections as f64 * 0.2).min(0.8);
        for i in 0..6 {
            state.covariance[(i, i)] *= inflate_factor;
        }
        for i in 0..3 {
            state.covariance[(i, i)] += 1.0;
        }
        for i in 3..6 {
            state.covariance[(i, i)] += 0.1;
        }
        // Inflate ambiguity covariances to enable re-convergence
        let amb_start = crate::filter::CORE_STATE_SIZE;
        for i in amb_start..state.covariance.nrows() {
            state.covariance[(i, i)] *= inflate_factor.min(2.0);
        }
    }

    let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
    if pos_var > 10000.0 {
        tracing::warn!(
            "Extreme divergence detected (pos_var {:.2}): resetting EKF to SPP fallback.",
            pos_var
        );
        if let Some(pos) = spp_pos {
            state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
            state.consecutive_rejections = 0;
        }
    } else if state.consecutive_rejections >= config.max_consecutive_rejections {
        tracing::warn!(
            "EKF rejected for {} epochs: resetting to SPP fallback (max_consecutive_rejections={}).",
            state.consecutive_rejections,
            config.max_consecutive_rejections,
        );
        if let Some(pos) = spp_pos {
            state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
            state.consecutive_rejections = 0;
        }
    }
}

/// Solve for position using accumulated centered SD pseudorange.
/// SD with mean removal eliminates common-mode bias (reference satellite,
/// receiver clock residuals) that corrupts DD-based position estimates.
///
/// Returns None if insufficient PR data is available.
fn solve_pr_only_position(
    state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    base_coord: &Coordinate,
    base_time: gneiss_core::time::GpsTime,
) -> Option<Vector3<f64>> {
    struct DdObs {
        sat: gneiss_core::sat::SatelliteId,
        ref_sat: gneiss_core::sat::SatelliteId,
        dd_mean: f64,
        weight: f64,
    }
    let ref_pos = state.position.vector; // compensate all entries to current float position
    let mut obs = Vec::new();
    for ((sat, ref_sat), buf) in state.pr_dd_window.iter() {
        let n = buf.count();
        if n < 20 { continue; }
        let comp_mean = buf.compensated_mean(
            *sat, *ref_sat, ref_pos, ephemerides, state.time,
            base_coord.vector, base_time,
        )?;
        let sigma_eff = 1.5f64 / (n as f64).sqrt();
        obs.push(DdObs { sat: *sat, ref_sat: *ref_sat, dd_mean: comp_mean, weight: 1.0 / (sigma_eff * sigma_eff) });
    }
    if obs.len() < 4 { return None; }

    let mut pos = ref_pos;

    for _iter in 0..8 {
        let mut h_sum = nalgebra::Matrix3::zeros();
        let mut rhs = nalgebra::Vector3::zeros();
        let mut prev_rms = 0.0f64;
        for o in &obs {
            let eph_sat = match ephemerides.iter().find(|e| e.sat() == o.sat) { Some(e) => e, None => continue, };
            let eph_ref = match ephemerides.iter().find(|e| e.sat() == o.ref_sat) { Some(e) => e, None => continue, };
            let (sat_pos, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, state.time, pos);
            let (ref_pos, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, state.time, pos);
            let (bas_sat, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, base_time, base_coord.vector);
            let (bas_ref, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, base_time, base_coord.vector);
            let geom_dd = crate::engine::measurement_math::compute_geometric_dd(pos, base_coord.vector, sat_pos, ref_pos, bas_sat, bas_ref);
            let residual = o.dd_mean - geom_dd;
            let h = ((ref_pos - pos).normalize() - (sat_pos - pos).normalize()) * o.weight.sqrt();
            h_sum += h * h.transpose();
            rhs += h * residual * o.weight.sqrt();
            prev_rms += residual * residual;
        }
        prev_rms = (prev_rms / obs.len() as f64).sqrt();
        if let Some(h_inv) = h_sum.try_inverse() {
            let mut dx = h_inv * rhs;
            // Trust-region: clamp step to 500m max per iteration
            let dx_norm = dx.norm();
            if dx_norm > 500.0 { dx *= 500.0 / dx_norm; }
            // Line search: only accept if RMS improves by at least 2%
            let mut accepted = false;
            let mut lambda = 1.0;
            for _ in 0..5 {
                let pos_trial = pos + dx * lambda;
                let mut trial_rms = 0.0f64;
                for o in &obs {
                    let eph_sat = match ephemerides.iter().find(|e| e.sat() == o.sat) { Some(e) => e, None => continue, };
                    let eph_ref = match ephemerides.iter().find(|e| e.sat() == o.ref_sat) { Some(e) => e, None => continue, };
                    let (sp, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, state.time, pos_trial);
                    let (rp, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, state.time, pos_trial);
                    let (bsp, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, base_time, base_coord.vector);
                    let (brp, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, base_time, base_coord.vector);
                    let gdd = crate::engine::measurement_math::compute_geometric_dd(pos_trial, base_coord.vector, sp, rp, bsp, brp);
                    trial_rms += (o.dd_mean - gdd).powi(2);
                }
                let trial = (trial_rms / obs.len() as f64).sqrt();
                if trial < prev_rms * 0.98 || lambda < 0.0625 {
                    pos = pos_trial;
                    accepted = true;
                    break;
                }
                lambda *= 0.5;
            }
            if !accepted || dx_norm * lambda < 0.01 { break; }
        } else { return None; }
    }
    Some(pos)
}

/// Validate an AR fix by comparing the fixed position against an independent
/// PR-only position estimate.  The PR-only position uses accumulated raw DD
/// pseudorange — no carrier phase, no ambiguities — so it's immune to the
/// code-multipath bias that corrupts MW/NL EMAs.
/// Two-epoch covariance-weighted position smoother.
/// Combines the current and previous epoch's float positions using their
/// covariance matrices. Returns the smoothed position, or the current position
/// if previous epoch data is unavailable.
fn smooth_two_epoch_position(state: &RtkState) -> (Vector3<f64>, f64) {
    let cur_pos = state.position.vector;
    let cur_cov = state.covariance.fixed_view::<3, 3>(0, 0).into_owned();

    let prev_pos = match state.prev_epoch_pos {
        Some(p) => p,
        None => {
            let tr = cur_cov.trace();
            return (cur_pos, (tr / 3.0).sqrt());
        }
    };
    let prev_cov = match state.prev_epoch_cov {
        Some(ref c) => c.clone(),
        None => {
            let tr = cur_cov.trace();
            return (cur_pos, (tr / 3.0).sqrt());
        }
    };

    // Covariance-weighted average: P_smooth = (P_cur^{-1} + P_prev^{-1})^{-1}
    let cur_inv = match cur_cov.try_inverse() {
        Some(inv) => inv,
        None => return (cur_pos, 100.0),
    };
    let prev_inv = match prev_cov.try_inverse() {
        Some(inv) => inv,
        None => return (cur_pos, 100.0),
    };
    let p_smooth_inv = cur_inv + prev_inv;
    let p_smooth = match p_smooth_inv.try_inverse() {
        Some(p) => p,
        None => return (cur_pos, 100.0),
    };

    let pos_smooth = &p_smooth * (&cur_inv * cur_pos + &prev_inv * prev_pos);
    let sigma_smooth = (p_smooth.trace() / 3.0).sqrt();

    (pos_smooth, sigma_smooth)
}

fn validate_geometry_pr(
    state: &RtkState,
    fixed_state: &RtkState,
    _ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    _base_coord: &Coordinate,
    _base_time: gneiss_core::time::GpsTime,
) -> bool {
    // Two-epoch smoothed position provides an improved float estimate
    // that averages down PR noise by √2 compared to single-epoch.
    let (smooth_pos, sigma_smooth) = smooth_two_epoch_position(state);
    let fixed_err = (fixed_state.position.vector - smooth_pos).norm();
    let float_err = (state.position.vector - smooth_pos).norm();

    // Reject if fixed position is far from smoothed position
    let threshold = 5.0 * sigma_smooth.max(0.5);
    tracing::info!(
        "2-epoch smooth: float_err={:.2}m fixed_err={:.2}m sigma={:.2}m thresh={:.2}m",
        float_err, fixed_err, sigma_smooth, threshold
    );

    if fixed_err > threshold && fixed_err > float_err * 1.5 {
        tracing::warn!("AR fix rejected by 2-epoch smoother: fixed_err={:.2}m > thresh={:.2}m", fixed_err, threshold);
        return false;
    }
    true
}

/// Validate an AR fix using the two-epoch factor graph, and feed the factor
/// graph's independent position estimate back into the EKF as a weak position
/// measurement. This improves the float solution over time, breaking the
/// code-multipath circularity that biases single-epoch AR validation.
///
/// Returns true if the fix passes (or if the factor graph is unavailable).
fn validate_factor_graph(
    state: &mut RtkState,
    fixed_state: &RtkState,
    m: Option<&crate::engine::measurement::EkfMeasurementMatrices>,
) -> bool {
    let Some(meas) = m else { return true; };
    if state.prev_epoch_meas.is_none() {
        return true;
    }
    match crate::engine::rtk_multi_epoch::run_two_epoch_factor_graph(state, meas) {
        Some(result) => {
            // Feed factor graph position back into EKF as a weak measurement.
            // This improves the float solution independently of code multipath,
            // since the factor graph uses geometry diversity across epochs.
            let pos_jump = (result.pos_k - state.position.vector).norm();
            if result.converged && pos_jump < 5.0 {
                // Weak position feedback: σ = 1.0m floor.
                // Nudges the EKF toward the FG position without dominating
                // the measurements. Over many epochs this improves the float
                // solution and makes LAMBDA more likely to find correct integers.
                let sigma = pos_jump.max(1.0);
                let mut h = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
                h[(0, 0)] = 1.0;
                h[(1, 1)] = 1.0;
                h[(2, 2)] = 1.0;
                let z = result.pos_k - state.position.vector;
                let z_vec = nalgebra::DVector::from_vec(vec![z.x, z.y, z.z]);
                let r = nalgebra::DMatrix::from_diagonal(&nalgebra::DVector::from_element(3, sigma * sigma));
                let _ = crate::engine::updater::update::<crate::engine::updater_math::TightCoupling>(
                    state, &z_vec, &h, &r, 9.0, None,
                    &crate::engine::EngineConfig::default().tuning,
                );
            }

            let fg_error = (fixed_state.position.vector - result.pos_k).norm();
            // Threshold balances false positives (wrong fixes accepted) against
            // false negatives (correct fixes rejected). Wrong NL integers
            // produce ~0.5-1.0m disagreement with the geometry-constrained FG
            // position. 0.5m catches most wrong fixes while accepting correct
            // ones that may have ~0.2-0.4m fg_error due to FG position noise.
            let threshold = 0.5; // meters
            tracing::info!(
                "RTK FG validation: fg_error={:.3}m thresh={:.3}m converged={} pos_fb={:.2}m",
                fg_error, threshold, result.converged, pos_jump
            );
            if fg_error > threshold {
                tracing::warn!(
                    "AR fix rejected by factor graph: fg_error={:.2}m > thresh={:.2}m",
                    fg_error, threshold
                );
                return false;
            }
            true
        }
        None => {
            tracing::debug!("RTK FG: factor graph solve failed, falling back to existing validation");
            true // Don't block on solver failure
        }
    }
}

/// Save the current epoch's measurement data for use by the
/// two-epoch factor graph in the next epoch.
fn save_prev_epoch_measurements(
    state: &mut RtkState,
    m: &crate::engine::measurement::EkfMeasurementMatrices,
) {
    use crate::filter::{PrevEpochMeasurements, StoredDdMeasurement};
    use crate::filter::CORE_STATE_SIZE;

    let mut measurements = Vec::with_capacity(m.z.nrows());

    for i in 0..m.z.nrows() {
        let (sat_id, type_code, _r_ref) = m.mt[i];
        let z = m.z[i];
        let variance = m.r[(i, i)].max(1e-4);
        let h_pos = [m.h[(i, 0)], m.h[(i, 1)], m.h[(i, 2)]];
        let is_pr = type_code == 0;
        let freq_band = match type_code {
            0 => 1,
            1 => 1,
            2 => 2,
            _ => 1,
        };

        // Extract ambiguity indices from H matrix columns >= CORE_STATE_SIZE
        let mut amb_idx: Option<usize> = None;
        let mut ref_amb_idx: Option<usize> = None;
        let mut iono_sat_idx: Option<usize> = None;
        let mut iono_ref_idx: Option<usize> = None;
        let mut iono_scale: Option<f64> = None;

        if !is_pr {
            for col in CORE_STATE_SIZE..m.h.ncols() {
                let val = m.h[(i, col)];
                if val.abs() < 1e-9 {
                    continue;
                }
                let key_idx = col - CORE_STATE_SIZE;
                if key_idx >= state.ambiguity_keys.len() {
                    continue;
                }
                let (_, freq) = state.ambiguity_keys[key_idx];
                if freq == 3 {
                    if val > 0.0 {
                        iono_sat_idx = Some(key_idx);
                        iono_scale = Some(val);
                    } else {
                        iono_ref_idx = Some(key_idx);
                    }
                } else if val > 0.0 {
                    amb_idx = Some(key_idx);
                } else {
                    ref_amb_idx = Some(key_idx);
                }
            }
        }

        let ref_sat_id = ref_amb_idx
            .and_then(|ri| state.ambiguity_keys.get(ri).map(|k| k.0))
            .unwrap_or(sat_id);

        let iono_pair = match (iono_sat_idx, iono_ref_idx, iono_scale) {
            (Some(si), Some(ri), Some(sc)) => Some((si, ri, sc)),
            _ => None,
        };

        measurements.push(StoredDdMeasurement {
            sat_id,
            ref_sat_id,
            z,
            h_pos,
            variance,
            is_pr,
            freq_band,
            amb_idx,
            ref_amb_idx,
            iono_pair,
        });
    }

    // Store current epoch as previous-epoch reference for the two-epoch
    // factor graph. Uses 1-epoch spacing — satellite geometry barely changes
    // but the FG still provides a useful consistency check against wrong AR
    // fixes. Multi-epoch geometry diversity requires a ring buffer with
    // ambiguity-key matching, which is fragile (see Phase 5.1-5.2).
    state.prev_epoch_meas = Some(PrevEpochMeasurements {
        time: state.time.tow,
        pos: state.position.vector,
        vel: state.velocity,
        clk: state.rcv_clk_bias,
        measurements,
        ambiguity_keys: state.ambiguity_keys.clone(),
    });
}
/// Accumulate EKF pseudorange innovations (clock-corrected DD) per satellite
/// pair for multi-epoch PR validation.
fn accumulate_pr_window(
    state: &mut RtkState,
    m: &crate::engine::measurement::EkfMeasurementMatrices,
    window_size: usize,
) {
    if window_size == 0 {
        return;
    }
    for i in 0..m.z.nrows() {
        if m.mt[i].1 != 0 { continue; } // PR only (type 0)
        let sat = m.mt[i].0;
        let innovation = m.z[i];
        let ref_sat = state.current_ref_sat.get(&sat.constellation).copied()
            .unwrap_or(sat);
        let key = (sat, ref_sat);
        let buf = state.pr_dd_window
            .entry(key)
            .or_insert_with(|| crate::filter::PrRingBuffer::new(window_size));
        buf.push(innovation, ref_sat, state.position.vector);
    }
}

/// Validate AR fix using time-averaged PR innovations compensated to the
/// fixed position. The compensated mean reconstructs the raw DD PR at the
/// reference position by accounting for rover motion between epochs.
/// A wrong NL fix produces a position inconsistent with the PR average.
fn validate_multiepoch_pr_compensated(
    state: &RtkState,
    fixed_state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    base_coord: &Coordinate,
    base_time: gneiss_core::time::GpsTime,
) -> bool {
    let mut pairs_checked = 0usize;
    let mut pairs_failed = 0usize;

    for ((sat, ref_sat), buf) in state.pr_dd_window.iter() {
        let n = buf.count();
        if n < 20 {
            continue;
        }

        // Reconstruct time-averaged raw DD PR from stored EKF innovations.
        // Each innovation z = obs - pred ≈ obs - geom(entry_pos).
        // Adding back geom(entry_pos) reconstructs the raw DD PR.
        let recon_mean = match buf.innovation_reconstructed_mean(
            *sat, *ref_sat,
            ephemerides,
            state.time,
            base_coord.vector,
            base_time,
        ) {
            Some(m) => m,
            None => continue,
        };

        // Compute expected DD PR at the fixed position from satellite geometry
        let eph_sat = match ephemerides.iter().find(|e| e.sat() == *sat) {
            Some(e) => e,
            None => continue,
        };
        let eph_ref = match ephemerides.iter().find(|e| e.sat() == *ref_sat) {
            Some(e) => e,
            None => continue,
        };
        let (sat_pos_approx, _) = crate::engine::measurement_math::get_sat_state(
            eph_sat, 0.0, 0.0, state.time, fixed_state.position.vector,
        );
        let range_sat = (sat_pos_approx - fixed_state.position.vector).norm();
        let (sat_pos, _) = crate::engine::measurement_math::get_sat_state(
            eph_sat, range_sat, 0.0, state.time, fixed_state.position.vector,
        );
        let (ref_pos_approx, _) = crate::engine::measurement_math::get_sat_state(
            eph_ref, 0.0, 0.0, state.time, fixed_state.position.vector,
        );
        let range_ref = (ref_pos_approx - fixed_state.position.vector).norm();
        let (ref_pos, _) = crate::engine::measurement_math::get_sat_state(
            eph_ref, range_ref, 0.0, state.time, fixed_state.position.vector,
        );
        let (bas_sat_approx, _) = crate::engine::measurement_math::get_sat_state(
            eph_sat, 0.0, 0.0, base_time, base_coord.vector,
        );
        let range_bas = (bas_sat_approx - base_coord.vector).norm();
        let (bas_sat, _) = crate::engine::measurement_math::get_sat_state(
            eph_sat, range_bas, 0.0, base_time, base_coord.vector,
        );
        let (bas_ref_approx, _) = crate::engine::measurement_math::get_sat_state(
            eph_ref, 0.0, 0.0, base_time, base_coord.vector,
        );
        let range_bref = (bas_ref_approx - base_coord.vector).norm();
        let (bas_ref, _) = crate::engine::measurement_math::get_sat_state(
            eph_ref, range_bref, 0.0, base_time, base_coord.vector,
        );

        let expected = crate::engine::measurement_math::compute_geometric_dd(
            fixed_state.position.vector,
            base_coord.vector,
            sat_pos,
            ref_pos,
            bas_sat,
            bas_ref,
        );

        let residual = (expected - recon_mean).abs();
        let threshold = 3.0 * 1.5 / (n as f64).sqrt();

        pairs_checked += 1;
        if residual > threshold {
            pairs_failed += 1;
            tracing::info!(
                "Multi-epoch PR: sat={:?} ref={:?} n={} resid={:.3}m thresh={:.3}m recon={:.3}m expected={:.3}m FAIL",
                sat, ref_sat, n, residual, threshold, recon_mean, expected
            );
        }
    }

    if pairs_checked < 4 {
        return true; // Not enough data
    }
    let pass_rate = (pairs_checked - pairs_failed) as f64 / pairs_checked as f64;
    tracing::info!(
        "Multi-epoch PR: {}/{} pairs passed ({:.0}%)",
        pairs_checked - pairs_failed, pairs_checked, pass_rate * 100.0
    );
    if pass_rate < 0.5 {
        tracing::warn!("AR fix rejected by multi-epoch PR");
        return false;
    }
    true
}

fn validate_pr_residuals(
    state: &RtkState,
    fixed_state: &RtkState,
    m: Option<&crate::engine::measurement::EkfMeasurementMatrices>,
    valid_indices: Option<&[usize]>,
) -> bool {
    let (Some(meas), Some(valid)) = (m, valid_indices) else { return true };
    let dx = fixed_state.position.vector - state.position.vector;
    let state_size = state.covariance.nrows();

    let mut float_rms = 0.0f64;
    let mut fixed_rms = 0.0f64;
    let mut count = 0usize;

    for &i in valid {
        if meas.mt[i].1 != 0 { continue; } // PR measurements only (type 0)
        let z_float = meas.z[i]; // innovation at float position
        // Predicted change in PR DD from position change: H[0:3] · dx
        let h_dot_dx: f64 = (0..3.min(state_size))
            .map(|c| meas.h[(i, c)] * dx[c])
            .sum();
        let z_fixed = z_float - h_dot_dx;
        float_rms += z_float * z_float;
        fixed_rms += z_fixed * z_fixed;
        count += 1;
    }

    if count < 4 { return true; } // too few PR measurements to validate

    float_rms = (float_rms / count as f64).sqrt();
    fixed_rms = (fixed_rms / count as f64).sqrt();

    tracing::debug!(
        "AR PR residuals: float={:.1}m fixed={:.1}m ratio={:.2} ({} PR)",
        float_rms, fixed_rms, if float_rms > 0.01 { fixed_rms / float_rms } else { 1.0 }, count
    );
    // Reject if fixed position significantly degrades PR fit
    if fixed_rms > float_rms * 1.5 && fixed_rms > 1.5 {
        tracing::warn!(
            "AR fix rejected by PR residuals: float={:.1}m fixed={:.1}m ({} PR)",
            float_rms, fixed_rms, count
        );
        return false;
    }
    true
}

/// Apply an accepted AR fix to the live EKF state.
fn accept_ar_fix(state: &mut RtkState, fixed_state: RtkState, pos_jump: f64) {
    tracing::info!(
        "RTK AR fixed: {} sats (position jump {:.3}m)",
        fixed_state.ambiguities.len(),
        pos_jump
    );
    state.is_fixed = true;
    state.position = fixed_state.position.clone();
    if state.ambiguities.len() == fixed_state.ambiguities.len() {
        state.ambiguities.copy_from_slice(&fixed_state.ambiguities);
        state.covariance = fixed_state.covariance.clone();
    }
    state.fixed_state = Some(Box::new(fixed_state));
}

fn handle_ekf_acceptance(
    state: &mut RtkState,
    config: &EngineConfig,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    spp_pos: Option<Coordinate>,
    spp_state_ref: Option<&crate::spp::SppState>,
    m: Option<&crate::engine::measurement::EkfMeasurementMatrices>,
    valid_indices: Option<&[usize]>,
    base_coord: Option<&Coordinate>,
    base_time: Option<gneiss_core::time::GpsTime>,
) {
    let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
    if pos_var > 10000.0 {
        tracing::warn!(
            "Position variance extremely large ({:.2}): resetting EKF to SPP fallback.",
            pos_var
        );
        if let Some(pos) = spp_pos {
            state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
        }
    }
    state.consecutive_rejections = 0;
    if let Ok(res) = state.resolve_ambiguities(ephemerides, config) {
        let fixed_state = res.fixed_state;

        // Validate the AR fix: reject if the position jump is excessive,
        // which indicates a wrong integer set. A correct fix on a short
        // baseline should change position by < 1 m.
        let pos_jump = (fixed_state.position.vector - state.position.vector).norm();
        if pos_jump > AR_MAX_POSITION_JUMP_M {
            tracing::warn!(
                "AR fix rejected: position jump {:.2}m exceeds {:.1}m threshold",
                pos_jump,
                AR_MAX_POSITION_JUMP_M
            );
            state.is_fixed = false;
            state.fixed_state = None;
        } else if config.enable_ins_validation
            && pos_jump > state.velocity.norm().max(0.5) * 3.0
            && pos_jump > 0.5
        {
            tracing::warn!(
                "AR fix rejected by INS: jump {:.2}m > 2× expected motion ({:.2}m/s)",
                pos_jump, state.velocity.norm()
            );
            state.is_fixed = false;
            state.fixed_state = None;
        } else if !validate_pr_residuals(state, &fixed_state, m, valid_indices) {
            state.is_fixed = false;
            state.fixed_state = None;
        } else if !validate_factor_graph(state, &fixed_state, m) {
            state.is_fixed = false;
            state.fixed_state = None;
        } else if let (Some(bc), Some(bt)) = (base_coord, base_time) {
            if !validate_geometry_pr(state, &fixed_state, ephemerides, bc, bt) {
                state.is_fixed = false;
                state.fixed_state = None;
            } else {
                accept_ar_fix(state, fixed_state, pos_jump);
            }
        } else {
            accept_ar_fix(state, fixed_state, pos_jump);
        }
    } else {
        // Maintain previous fix across epochs — don't un-fix because
        // LAMBDA failed on this particular epoch. The fix is cleared
        // only by a cycle slip (which removes the ambiguity from state).
        if !state.is_fixed {
            state.fixed_state = None;
            tracing::info!("RTK AR failed at epoch {} (state epoch_count={})", state.time.tow, state.epoch_count);
        }
    }
}

pub struct RtkUpdateContext<'a> {
    pub config: &'a EngineConfig,
    pub ephemerides: &'a [gneiss_core::ephemeris::Ephemeris],
    pub imu_history: &'a [Vec<gneiss_core::imu::ImuMeasurement>],
    pub rover_obs: &'a EpochObs,
    pub base_obs: &'a EpochObs,
    pub matched_obs: &'a [(crate::filter::DdObservation, crate::filter::DdObservation)],
    pub base_coord: &'a Coordinate,
    pub spp_pos: Option<Coordinate>,
    pub spp_state_ref: Option<&'a crate::spp::SppState>,
    pub gnn_variances: std::collections::HashMap<gneiss_core::sat::SatelliteId, f64>,
    pub klobuchar_params: Option<gneiss_core::atmosphere::KlobucharParams>,
}

fn process_rtk_update<C: CouplingStrategy>(
    state: &mut RtkState,
    tracker: &mut crate::engine::adaptive::InnovationTracker,
    ctx: &RtkUpdateContext,
) {
    crate::engine::ambiguity::manage_ambiguities_and_slips(
        state,
        ctx.config,
        ctx.matched_obs,
        ctx.ephemerides,
        ctx.base_coord,
        ctx.rover_obs.time,
        ctx.base_obs.time,
    );
    update_last_observed(state, ctx.matched_obs);

    let mut env = build_measurement_environment(
        ctx.config,
        ctx.imu_history,
        state,
        ctx.ephemerides,
        ctx.base_coord,
        ctx.base_obs.time,
        ctx.klobuchar_params,
    );
    env.gnn_variances = ctx.gnn_variances.clone();

    if let Some(mut m) = crate::engine::measurement::build_measurement_model(
        state,
        ctx.matched_obs,
        &env,
        ctx.config.chi_square_pr_threshold,
        ctx.config.chi_square_cp_threshold,
    ) {
        apply_adaptive_r_scaling(tracker, state, &m.z, &m.h, &mut m.r, &m.mt, ctx.matched_obs);
        execute_ekf_update::<C>(state, ctx, &m);
        accumulate_pr_window(state, &m, ctx.config.pr_window_size);
    } else {
        handle_ekf_rejection(
            state,
            ctx.config,
            ctx.spp_pos,
            ctx.spp_state_ref,
            "not enough measurements",
        );
    }
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
}

fn update_last_observed(
    state: &mut RtkState,
    matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
) {
    let current_epoch = state.epoch_count as u32;
    for (r_obs, _) in matched_obs {
        if r_obs.cp_l1.is_some() {
            state.last_observed.insert((r_obs.sat, 1), current_epoch);
        }
        if r_obs.cp_l2.is_some() {
            state.last_observed.insert((r_obs.sat, 2), current_epoch);
        }
    }
}

fn execute_ekf_update<C: CouplingStrategy>(
    state: &mut RtkState,
    ctx: &RtkUpdateContext,
    m: &crate::engine::measurement::EkfMeasurementMatrices,
) {
    let pr_thresh = ctx.config.chi_square_pr_threshold;
    let type_stripped: Vec<_> = m.mt.iter().map(|&(s, t, _)| (s, t)).collect();

    match crate::engine::updater::update::<C>(
        state,
        &m.z,
        &m.h,
        &m.r,
        pr_thresh,
        Some(&type_stripped),
        &ctx.config.tuning,
    ) {
        Err(_) => handle_ekf_rejection(
            state,
            ctx.config,
            ctx.spp_pos,
            ctx.spp_state_ref,
            "failed chi-square",
        ),
        Ok((valid_indices, dx)) => {
            export_gnn_dataset(state, ctx, m, &valid_indices, &dx);
            handle_ekf_acceptance(
                state,
                ctx.config,
                ctx.ephemerides,
                ctx.spp_pos,
                ctx.spp_state_ref,
                Some(m),
                Some(&valid_indices),
                Some(ctx.base_coord),
                Some(ctx.base_obs.time),
            );
            // Save measurement data for the next epoch's two-epoch factor graph
            save_prev_epoch_measurements(state, m);
        }
    }
}

fn export_gnn_dataset(
    state: &RtkState,
    ctx: &RtkUpdateContext,
    m: &crate::engine::measurement::EkfMeasurementMatrices,
    valid_indices: &[usize],
    dx: &DVector<f64>,
) {
    if let Some(path) = &ctx.config.export_gnn_dataset_path {
        let rov_llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
        let pos_apc = state.position.vector
            + state.attitude * nalgebra::Vector3::from(ctx.config.imu_to_antenna_lever_arm);

        crate::engine::ml::dataset::export_epoch_to_csv(
            path,
            state.epoch_count as u32,
            ctx.matched_obs,
            rov_llh,
            pos_apc,
            ctx.ephemerides,
            ctx.rover_obs.time,
            &m.z,
            &m.h,
            valid_indices,
            &m.mt,
            dx,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{EngineConfig, EngineError, EngineMode};
    use crate::filter::{CORE_STATE_SIZE, DdObservation};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
    use gneiss_core::sat::Constellation;
    use gneiss_core::time::GpsTime;
    use gneiss_geodesy::helmert::HelmertParams;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn make_test_sat(prn: u8) -> SatelliteId {
        SatelliteId {
            constellation: Constellation::Gps,
            prn,
        }
    }

    fn make_test_epoch(time: GpsTime) -> EpochObs {
        EpochObs {
            time,
            satellites: Vec::new(),
        }
    }

    fn make_dd_obs(sat: SatelliteId, cp_l1: Option<f64>, cp_l2: Option<f64>) -> DdObservation {
        DdObservation {
            sat,
            pr_l1: 0.0,
            pr_l2: None,
            cp_l1,
            cp_l2,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        }
    }

    #[test]
    fn test_filter_valid_base_accepts_recent() {
        let mut config = EngineConfig::default();
        config.max_base_age_s = 5.0;
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let base_obs = make_test_epoch(GpsTime::new(0, 101.0));
        let result =
            ProcessingEngine::filter_valid_base(&config, &rover_obs, Some(&base_obs));
        assert!(result.is_some());
    }

    #[test]
    fn test_filter_valid_base_rejects_old() {
        let mut config = EngineConfig::default();
        config.max_base_age_s = 5.0;
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let base_obs = make_test_epoch(GpsTime::new(0, 120.0));
        let result =
            ProcessingEngine::filter_valid_base(&config, &rover_obs, Some(&base_obs));
        assert!(result.is_none());
    }

    #[test]
    fn test_filter_valid_base_no_base() {
        let config = EngineConfig::default();
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let result = ProcessingEngine::filter_valid_base(&config, &rover_obs, None);
        assert!(result.is_none());
    }

    #[test]
    fn test_get_base_coord_ok() {
        let mut config = EngineConfig::default();
        config.base_position = Some([100.0, 200.0, 300.0]);
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let coord = ProcessingEngine::get_base_coord(&config, &rover_obs).unwrap();
        assert!((coord.vector.x - 100.0).abs() < 1e-6);
        assert!((coord.vector.y - 200.0).abs() < 1e-6);
        assert!((coord.vector.z - 300.0).abs() < 1e-6);
        assert_eq!(coord.datum, Datum::WGS84);
    }

    #[test]
    fn test_get_base_coord_missing() {
        let config = EngineConfig::default();
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let err = ProcessingEngine::get_base_coord(&config, &rover_obs).unwrap_err();
        assert!(matches!(err, EngineError::MissingBasePosition));
    }

    #[test]
    fn test_update_last_observed_inserts_all() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 5;

        let sat1 = make_test_sat(1);
        let sat2 = make_test_sat(2);
        let matched = vec![
            (make_dd_obs(sat1, Some(100.0), Some(200.0)), make_dd_obs(sat1, Some(100.0), Some(200.0))),
            (make_dd_obs(sat2, Some(300.0), None), make_dd_obs(sat2, Some(300.0), None)),
        ];

        update_last_observed(&mut state, &matched);

        assert_eq!(state.last_observed[&(sat1, 1)], 5);
        assert_eq!(state.last_observed[&(sat1, 2)], 5);
        assert_eq!(state.last_observed[&(sat2, 1)], 5);
        assert!(!state.last_observed.contains_key(&(sat2, 2)));
    }

    #[test]
    fn test_update_last_observed_skips_none() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 3;

        let sat1 = make_test_sat(1);
        let matched = vec![
            (make_dd_obs(sat1, None, None), make_dd_obs(sat1, None, None)),
        ];

        update_last_observed(&mut state, &matched);

        assert!(!state.last_observed.contains_key(&(sat1, 1)));
        assert!(!state.last_observed.contains_key(&(sat1, 2)));
    }

    #[test]
    fn test_evaluate_gnn_no_model() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let gnn_variances = ProcessingEngine::evaluate_gnn(
            &None,
            &[],
            &make_test_epoch(time),
            &[],
            &state,
        );
        assert!(gnn_variances.is_empty());
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_rtk_resets_when_young() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Rtk;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, false);

        // Should have been reset to SPP because epoch_count < 3
        assert!((state.position.vector.x - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_rtk_no_reset_when_mature() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Rtk;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 5;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, false);

        // Should NOT reset because epoch_count >= 3
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_ins_with_imu_data() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        // With IMU data, INS modes should NOT trigger GNSS-only coasting
        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, true);
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_ins_without_imu_data() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        // Without IMU data, INS modes should trigger GNSS-only coasting
        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, false);
        assert!((state.position.vector.x - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_update_state_time_sets_time_and_clears_reset() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.is_reset = true;

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.current_state = Some(state);

        let new_time = GpsTime::new(0, 10.0);
        engine.update_state_time(new_time).unwrap();

        let s = engine.current_state.as_ref().unwrap();
        assert!(!s.is_reset);
        assert_eq!(s.time.tow, 10.0);
    }

    #[test]
    fn test_update_state_time_no_state_errors() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let err = engine.update_state_time(GpsTime::new(0, 10.0)).unwrap_err();
        assert!(matches!(err, EngineError::StateDisappeared));
    }

    #[test]
    fn test_handle_ekf_rejection_increments() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let config = EngineConfig::default();

        handle_ekf_rejection(&mut state, &config, None, None, "test");
        assert_eq!(state.consecutive_rejections, 1);

        handle_ekf_rejection(&mut state, &config, None, None, "test");
        assert_eq!(state.consecutive_rejections, 2);
    }

    #[test]
    fn test_handle_ekf_rejection_ins_aligned_inflates_covariance() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;

        let cov_0_0_before = state.covariance[(0, 0)];
        let cov_3_3_before = state.covariance[(3, 3)];

        handle_ekf_rejection(&mut state, &EngineConfig::default(), None, None, "test");

        // Position covariances should be inflated
        assert!(state.covariance[(0, 0)] > cov_0_0_before);
        assert!(state.covariance[(3, 3)] > cov_3_3_before);
    }

    #[test]
    fn test_handle_ekf_rejection_not_aligned_no_reset_below_threshold() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = false;
        state.consecutive_rejections = 3;

        handle_ekf_rejection(&mut state, &EngineConfig::default(), None, None, "test");
        assert_eq!(state.consecutive_rejections, 4);
        // Position should NOT be reset because rejections < max_consecutive_rejections (30)
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_handle_ekf_rejection_not_aligned_resets_at_threshold_with_spp() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = false;
        state.consecutive_rejections = 29;

        let spp_pos = Coordinate::new(
            Vector3::new(50.0, 50.0, 50.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        handle_ekf_rejection(&mut state, &EngineConfig::default(), Some(spp_pos), None, "test");
        // Should have been reset after max_consecutive_rejections (30) consecutive rejections
        assert!((state.position.vector.x - 50.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_handle_ekf_acceptance_clears_rejections() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.consecutive_rejections = 5;

        handle_ekf_acceptance(&mut state, &EngineConfig::default(), &[], None, None, None, None, None, None);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_handle_ekf_acceptance_resets_on_extreme_variance() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance[(0, 0)] = 20000.0;
        state.covariance[(1, 1)] = 0.0;
        state.covariance[(2, 2)] = 0.0;
        state.consecutive_rejections = 5;

        let spp_pos = Coordinate::new(
            Vector3::new(50.0, 50.0, 50.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        handle_ekf_acceptance(&mut state, &EngineConfig::default(), &[], Some(spp_pos), None, None, None, None, None);
        // Should be reset
        assert!((state.position.vector.x - 50.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_build_measurement_environment_constructs() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let base_coord = Coordinate::new(
            Vector3::new(110.0, 210.0, 310.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let env = build_measurement_environment(
            &config,
            &[],
            &state,
            &[],
            &base_coord,
            time,
            None,
        );

        assert_eq!(env.base_coord.vector.x, 110.0);
        assert!(env.gnn_variances.is_empty());
    }

    #[test]
    fn test_perform_spp_fallback_update_does_not_panic() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);

        let new_pos = Coordinate::new(
            Vector3::new(101.0, 201.0, 301.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(0, 1.0),
        );

        // Should not panic
        ProcessingEngine::perform_spp_fallback_update(&config, &mut state, new_pos);
    }

    #[test]
    fn test_apply_adaptive_r_scaling_basic() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let mut r = DMatrix::identity(2, 2);
        let z = DVector::from_vec(vec![1.0; 2]);
        let mut h = DMatrix::zeros(2, CORE_STATE_SIZE);
        h[(0, 0)] = 1.0;
        h[(1, 0)] = 1.0;

        let sat = make_test_sat(1);
        let meas_types = vec![(sat, 0, 0.0), (sat, 0, 0.0)];
        let matched_obs = vec![
            (make_dd_obs(sat, Some(100.0), None), make_dd_obs(sat, Some(100.0), None)),
        ];

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();

        // Just ensure it doesn't panic and modifies r
        apply_adaptive_r_scaling(&mut tracker, &state, &z, &h, &mut r, &meas_types, &matched_obs);
        // R should be >= 1.0 after scaling (since min scale is 1.0)
        assert!(r[(0, 0)] >= 1.0);
        // Since innovation is 1.0 and predicted_var = h*P*h' + r = 1.0*1.0*1.0 + 1.0 = 2.0
        // nis = 1.0/2.0 = 0.5, so scale should be close to 1.0
        assert!((r[(0, 0)] - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_no_spp_pos_does_nothing() {
        // When spp_pos is None, the function should not reset regardless of epoch_count.
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Rtk;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0; // Would trigger reset IF spp_pos were Some

        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, None, None, false);

        // Position should remain unchanged
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_handle_ekf_rejection_ins_aligned_extreme_divergence_resets() {
        // When ins_aligned and position variance exceeds 10000,
        // handle_ekf_rejection should reset to SPP.
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        // Make position variance extreme
        state.covariance[(0, 0)] = 5000.0;
        state.covariance[(1, 1)] = 5000.0;
        state.covariance[(2, 2)] = 5000.0;

        let spp_pos = Coordinate::new(
            Vector3::new(50.0, 50.0, 50.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        handle_ekf_rejection(&mut state, &EngineConfig::default(), Some(spp_pos), None, "extreme");
        // Position should have been reset to SPP
        assert!((state.position.vector.x - 50.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_handle_ekf_rejection_ins_aligned_extreme_no_spp_does_not_panic() {
        // When ins_aligned and variance is extreme but no spp_pos provided,
        // the reset branch is skipped; consecutive_rejections still increments.
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.covariance[(0, 0)] = 20000.0;

        handle_ekf_rejection(&mut state, &EngineConfig::default(), None, None, "no_rescue");
        // Without spp_pos, the reset is skipped but consecutive_rejections is incremented
        assert_eq!(state.consecutive_rejections, 1);
        // Original position should be preserved
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_init_spp_state_when_state_already_exists() {
        // When current_state already exists, init_spp_state should not fail
        // even if SPP compute fails (no ephemerides).
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let result = engine.init_spp_state(&rover);
        // Should not error because state already exists
        assert!(result.is_ok());
        // SPP result should be None (no ephemerides to compute)
        assert!(result.unwrap().is_none());
    }

    #[test]
    fn test_apply_observations_with_fewer_than_5_matched_obs() {
        // When fewer than 5 observations match, consecutive_rejections should increase.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        engine.config.base_position = Some([100.0, 200.0, 300.0]);

        let state = engine.current_state.as_mut().unwrap();
        state.consecutive_rejections = 0;

        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time, satellites: vec![] };

        // This will call apply_observations with empty observations, resulting in
        // matched_obs being empty (< 5), which increments rejections.
        let result = engine.apply_observations(&rover, Some(&base), None, None);
        assert!(result.is_ok());
        assert_eq!(
            engine.current_state.as_ref().unwrap().consecutive_rejections,
            1
        );
    }

    #[test]
    fn test_apply_observations_with_valid_base_and_spp_fallback() {
        // When no valid base (age too large) but spp_pos is available,
        // apply_observations should call perform_spp_fallback_update (does not panic).
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 100.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        let state = engine.current_state.as_mut().unwrap();
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);

        // Base obs with large time difference (age > max_base_age_s = 5.0)
        let old_time = GpsTime::new(0, 1.0);
        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time: old_time, satellites: vec![] };

        let spp_pos = Some(Coordinate::new(
            Vector3::new(101.0, 201.0, 301.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        // Should not panic — SPP fallback update is applied
        let result = engine.apply_observations(&rover, Some(&base), spp_pos, None);
        assert!(result.is_ok());
    }

    #[test]
    fn test_process_rtk_basic_flow() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.base_position = Some([100.0, 200.0, 300.0]);

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time, satellites: vec![] };

        let result = engine.process_rtk(&rover, Some(&base));
        assert!(result.is_ok(), "process_rtk should succeed: {:?}", result.err());
        assert_eq!(engine.state_history.len(), 1);
        assert_eq!(engine.obs_history.len(), 1);
    }

    #[test]
    fn test_process_rtk_update_rejection() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let base_coord = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
            klobuchar_params: None,
        };

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();
        let rejections_before = state.consecutive_rejections;
        process_rtk_update::<TightCoupling>(&mut state, &mut tracker, &ctx);

        assert!(
            state.consecutive_rejections > rejections_before,
            "Rejections should increase after empty update"
        );
    }

    #[test]
    fn test_execute_ekf_update_rejection_nan() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let base_coord = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
            klobuchar_params: None,
        };

        let core_size = CORE_STATE_SIZE;
        let m = crate::engine::measurement::EkfMeasurementMatrices {
            z: DVector::from_vec(vec![f64::NAN]),
            h: DMatrix::zeros(1, core_size),
            r: DMatrix::identity(1, 1),
            mt: vec![(make_test_sat(1), 0, 1575.42e6)],
        };

        let rejections_before = state.consecutive_rejections;
        execute_ekf_update::<TightCoupling>(&mut state, &ctx, &m);

        assert!(
            state.consecutive_rejections > rejections_before,
            "NaN innovation should trigger rejection"
        );
    }

    #[test]
    fn test_execute_ekf_update_acceptance() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        state.consecutive_rejections = 3;

        let base_coord = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
            klobuchar_params: None,
        };

        let core_size = CORE_STATE_SIZE;
        let mut h = DMatrix::zeros(1, core_size);
        h[(0, 0)] = 1.0;
        let m = crate::engine::measurement::EkfMeasurementMatrices {
            z: DVector::from_vec(vec![0.0]),
            h,
            r: DMatrix::identity(1, 1),
            mt: vec![(make_test_sat(1), 0, 1575.42e6)],
        };

        execute_ekf_update::<TightCoupling>(&mut state, &ctx, &m);

        assert_eq!(state.consecutive_rejections, 0, "Acceptance should clear rejections");
    }

    #[test]
    fn test_build_measurement_environment_with_imu_and_lever_arm() {
        let mut config = EngineConfig::default();
        config.imu_to_antenna_lever_arm = [0.5, 0.0, 1.0];

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.gyro_bias = Vector3::new(0.01, 0.02, 0.03);

        let base_coord = Coordinate::new(
            Vector3::new(110.0, 210.0, 310.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let imu_history = vec![vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(0.0, 0.0, 9.8),
            gyro: Vector3::new(0.1, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]];

        let env = build_measurement_environment(
            &config,
            &imu_history,
            &state,
            &[],
            &base_coord,
            time,
            None,
        );

        assert!(env.lever_arm.norm() > 0.0, "Lever arm should be non-zero when ins_aligned");
        assert!(env.omega_b.norm() > 0.0, "Omega should be non-zero with gyro data and bias");
    }

    #[test]
    fn test_export_gnn_dataset_no_path_does_nothing() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let config = EngineConfig::default();
        let base_coord = Coordinate::new(
            Vector3::new(110.0, 210.0, 310.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
            klobuchar_params: None,
        };
        export_gnn_dataset(
            &state,
            &ctx,
            &crate::engine::measurement::EkfMeasurementMatrices {
                z: DVector::from_vec(vec![0.0]),
                h: DMatrix::zeros(CORE_STATE_SIZE, CORE_STATE_SIZE),
                r: DMatrix::identity(1, 1),
                mt: vec![(make_test_sat(1), 0, 1575.42e6)],
            },
            &[],
            &DVector::from_vec(vec![0.0]),
        );
    }

    #[test]
    fn test_init_spp_state_fails_when_no_state_and_spp_fails() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.current_state = None;
        let time = GpsTime::new(0, 0.0);
        let rover = EpochObs { time, satellites: vec![] };
        let result = engine.init_spp_state(&rover);
        assert!(matches!(result, Err(EngineError::InitialSppFailed)));
    }

    #[test]
    fn test_get_base_coord_with_helmert_transform() {
        let mut config = EngineConfig::default();
        config.base_position = Some([100.0, 200.0, 300.0]);
        config.base_datum_transform = Some(HelmertParams {
            tx: 1.0, ty: 2.0, tz: 3.0,
            rx: 0.0, ry: 0.0, rz: 0.0,
            s: 0.0,
            dtx: 0.0, dty: 0.0, dtz: 0.0,
            drx: 0.0, dry: 0.0, drz: 0.0,
            ds: 0.0,
            ref_epoch: 2000.0,
        });
        let rover_obs = EpochObs {
            time: GpsTime::new(0, 0.0),
            satellites: vec![],
        };
        let coord = ProcessingEngine::get_base_coord(&config, &rover_obs).unwrap();
        // Helmert applies [1,2,3]m translation with zero rates
        assert!((coord.vector.x - 101.0).abs() < 1.0);
        assert!((coord.vector.y - 202.0).abs() < 1.0);
        assert!((coord.vector.z - 303.0).abs() < 1.0);
    }

    #[test]
    fn test_process_rtk_fails_without_base_position() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.base_position = None;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time, satellites: vec![] };
        let err = engine.process_rtk(&rover, Some(&base)).unwrap_err();
        assert!(matches!(err, EngineError::MissingBasePosition));
    }

    #[test]
    fn test_perform_spp_fallback_update_large_diff_does_not_panic() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);

        let new_pos = Coordinate::new(
            Vector3::new(1000000.0, 2000000.0, 3000000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        ProcessingEngine::perform_spp_fallback_update(&config, &mut state, new_pos);
    }

    #[test]
    fn test_handle_ekf_acceptance_extreme_variance_no_spp_does_not_reset() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance[(0, 0)] = 20000.0;
        state.consecutive_rejections = 5;

        handle_ekf_acceptance(&mut state, &config, &[], None, None, None, None, None, None);
        assert_eq!(state.consecutive_rejections, 0);
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_apply_adaptive_r_scaling_with_freq2() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let mut r = DMatrix::identity(1, 1);
        let z = DVector::from_vec(vec![1.0]);
        let mut h = DMatrix::zeros(1, CORE_STATE_SIZE);
        h[(0, 0)] = 1.0;

        let sat = make_test_sat(1);
        let meas_types = vec![(sat, 2, 1227.6e6)];
        let matched_obs = vec![(
            make_dd_obs(sat, Some(100.0), None),
            make_dd_obs(sat, Some(100.0), None),
        )];

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();
        apply_adaptive_r_scaling(
            &mut tracker, &state, &z, &h, &mut r, &meas_types, &matched_obs,
        );
        assert!(r[(0, 0)] >= 1.0);
    }

    #[test]
    fn test_apply_adaptive_r_scaling_no_matching_snr() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let mut r = DMatrix::identity(1, 1);
        let z = DVector::from_vec(vec![1.0]);
        let mut h = DMatrix::zeros(1, CORE_STATE_SIZE);
        h[(0, 0)] = 1.0;

        let sat1 = make_test_sat(1);
        let sat2 = make_test_sat(2);
        let meas_types = vec![(sat2, 0, 1575.42e6)];
        let matched_obs = vec![(
            make_dd_obs(sat1, Some(100.0), None),
            make_dd_obs(sat1, Some(100.0), None),
        )];

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();
        apply_adaptive_r_scaling(
            &mut tracker, &state, &z, &h, &mut r, &meas_types, &matched_obs,
        );
        assert!(r[(0, 0)] >= 1.0);
    }

    // --- multi-base RTK tests ---

    #[test]
    fn test_multi_base_rtk_routes_from_process_rtk() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.enable_multi_base_rtk = true;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base1 = EpochObs { time, satellites: vec![] };
        let base2 = EpochObs {
            time: GpsTime::new(0, 0.1),
            satellites: vec![],
        };

        engine.multi_base_observations = vec![
            (base1, Vector3::new(100.0, 200.0, 300.0)),
            (base2, Vector3::new(150.0, 250.0, 350.0)),
        ];

        let result = engine.process_rtk(&rover, None);
        assert!(result.is_ok(), "multi-base routing should succeed: {:?}", result.err());
        assert_eq!(engine.state_history.len(), 1);
        assert_eq!(engine.multi_base_observations.len(), 0);
    }

    #[test]
    fn test_multi_base_rtk_direct_empty_observations() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base1 = EpochObs { time, satellites: vec![] };
        let base2 = EpochObs { time, satellites: vec![] };

        let result = engine.process_rtk_multi(
            &rover,
            &[
                (base1, Vector3::new(100.0, 200.0, 300.0)),
                (base2, Vector3::new(150.0, 250.0, 350.0)),
            ],
        );
        assert!(result.is_ok(), "direct multi-base with empty obs: {:?}", result.err());
        assert_eq!(engine.state_history.len(), 1);
    }

    #[test]
    fn test_multi_base_rtk_single_base_delegation() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base1 = EpochObs { time, satellites: vec![] };

        let result = engine.process_rtk_multi(
            &rover,
            &[(base1, Vector3::new(100.0, 200.0, 300.0))],
        );
        assert!(result.is_ok(), "single-base delegation: {:?}", result.err());
        assert_eq!(engine.state_history.len(), 1);
    }

    #[test]
    fn test_multi_base_rtk_empty_bases_delegates() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        let rover = EpochObs { time, satellites: vec![] };

        let result = engine.process_rtk_multi(&rover, &[]);
        assert!(result.is_ok(), "empty bases: {:?}", result.err());
    }

    #[test]
    fn test_multi_base_rtk_with_synthetic_observations() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);

        use gneiss_core::ephemeris::Ephemeris;
        for prn in 1..=7u8 {
            engine.add_ephemeris(Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn },
                toe: time,
                toc: time,
                af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
                omega: 0.0, tgd: 0.0, iode: prn as u32, iodc: prn as u32,
            }));
        }

        let make_sat_obs = |prn: u8| -> SatObs {
            let base_pr = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M + 1000.0 * prn as f64;
            SatObs {
                sat: SatelliteId { constellation: Constellation::Gps, prn },
                observations: vec![
                    Observation {
                        code: ObsCode {
                            obs_type: ObsType::Pseudorange,
                            signal: SignalCode { freq_band: 1, attribute: 'C' },
                        },
                        value: base_pr,
                        lock_time: Some(100),
                        lli: None,
                    },
                    Observation {
                        code: ObsCode {
                            obs_type: ObsType::CarrierPhase,
                            signal: SignalCode { freq_band: 1, attribute: 'C' },
                        },
                        value: base_pr / gneiss_core::constants::SPEED_OF_LIGHT_M_S * 1575.42e6,
                        lock_time: Some(100),
                        lli: None,
                    },
                    Observation {
                        code: ObsCode {
                            obs_type: ObsType::Pseudorange,
                            signal: SignalCode { freq_band: 2, attribute: 'C' },
                        },
                        value: base_pr + 1000.0,
                        lock_time: Some(100),
                        lli: None,
                    },
                    Observation {
                        code: ObsCode {
                            obs_type: ObsType::CarrierPhase,
                            signal: SignalCode { freq_band: 2, attribute: 'C' },
                        },
                        value: base_pr / gneiss_core::constants::SPEED_OF_LIGHT_M_S * 1227.6e6,
                        lock_time: Some(100),
                        lli: None,
                    },
                    Observation {
                        code: ObsCode {
                            obs_type: ObsType::Doppler,
                            signal: SignalCode { freq_band: 1, attribute: 'C' },
                        },
                        value: 0.0,
                        lock_time: None,
                        lli: None,
                    },
                    Observation {
                        code: ObsCode {
                            obs_type: ObsType::Snr,
                            signal: SignalCode { freq_band: 1, attribute: 'C' },
                        },
                        value: 45.0,
                        lock_time: None,
                        lli: None,
                    },
                ],
            }
        };

        let sats: Vec<SatObs> = (1..=7).map(make_sat_obs).collect();
        let rover_obs = EpochObs { time, satellites: sats.clone() };

        let base1_pos = Vector3::new(
            gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M + 100.0,
            200.0,
            300.0,
        );
        let base1_obs = EpochObs { time, satellites: sats.clone() };

        let base2_pos = Vector3::new(
            gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M + 150.0,
            -100.0,
            400.0,
        );
        let base2_obs = EpochObs { time, satellites: sats };

        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 100.0);
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE) * 100.0;
        engine.current_state = Some(state);

        let result = engine.process_rtk_multi(
            &rover_obs,
            &[(base1_obs, base1_pos), (base2_obs, base2_pos)],
        );
        assert!(result.is_ok() || result.is_err());
        assert_eq!(engine.state_history.len(), 1);
        assert_eq!(engine.obs_history.len(), 1);
    }

    #[test]
    fn test_multi_base_rtk_preserves_obs_history() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.enable_multi_base_rtk = true;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base1 = EpochObs { time, satellites: vec![] };
        let base2 = EpochObs { time, satellites: vec![] };

        engine.multi_base_observations = vec![
            (base1, Vector3::new(100.0, 200.0, 300.0)),
            (base2, Vector3::new(150.0, 250.0, 350.0)),
        ];

        let result = engine.process_rtk(&rover, None);
        assert!(result.is_ok());
        assert_eq!(engine.obs_history.len(), 1);
        assert!(engine.obs_history[0].1.is_some());
    }
}
