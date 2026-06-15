use crate::engine::{EngineMode, EngineError, EngineConfig};
use crate::filter::RtkState;
use gneiss_core::obs::EpochObs;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::sat::SatelliteId;
use nalgebra::{Vector3, DVector, DMatrix};
use super::ProcessingEngine;
use crate::engine::matcher::match_observations;

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
) {
    let state_size = state.covariance.nrows();
    for i in 0..z.nrows() {
        let h_row_vec: Vec<f64> = (0..state_size).map(|c| h[(i, c)]).collect();
        let h_row_mat = DMatrix::from_row_slice(1, state_size, &h_row_vec);
        let s_ii = (&h_row_mat * &state.covariance * h_row_mat.transpose())[(0, 0)] + r[(i, i)];
        let (sat, mtype, _) = meas_types[i];
        // Determine freq band from measurement type: 0=PR_L1, 1=CP_L1, 2=CP_L2, 3=Dop
        let freq = match mtype { 2 => 2, _ => 1 };
        let scale = tracker.update_and_scale(sat, freq, z[i], s_ii);
        r[(i, i)] *= scale;
    }
}

impl ProcessingEngine {

    fn init_spp_state(&mut self, rover_obs: &EpochObs) -> Result<Option<crate::spp::SppState>, EngineError> {
        let spp_res = crate::spp::compute_spp(rover_obs, &self.ephemerides, self.klobuchar_params.as_ref(), &crate::spp::SppConfig::default(), None).ok();
        
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

    fn evaluate_gnss_only_coasting(config: &EngineConfig, state: &mut RtkState, spp_pos: Option<Coordinate>, spp_state_ref: Option<&crate::spp::SppState>, had_imu_data: bool) {
        let use_gnss_only_seed = matches!(config.mode, EngineMode::Rtk | EngineMode::Ppp)
            || (matches!(config.mode, EngineMode::RtkIns | EngineMode::PppIns | EngineMode::SppIns | EngineMode::RtkInsLooselyCoupled | EngineMode::SppInsLooselyCoupled | EngineMode::PppInsLooselyCoupled) && !had_imu_data);
        
        if use_gnss_only_seed {
            tracing::warn!("epoch_count = {}", state.epoch_count);
            let need_spp_reset = if spp_pos.is_some() {
                state.epoch_count < 3
            } else { false };

            if need_spp_reset {
                if let Some(pos) = spp_pos {
                    tracing::warn!("Resetting EKF state to SPP due to divergence/startup.");
                    state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
                }
            }
        }
    }

    fn filter_valid_base<'a>(config: &EngineConfig, rover_obs: &EpochObs, base_obs: Option<&'a EpochObs>) -> Option<&'a EpochObs> {
        base_obs.filter(|b| {
            let age = (rover_obs.time.tow - b.time.tow).abs();
            if age > config.max_base_age_s {
                tracing::trace!("Base observation rejected due to Age of Differential ({:.1}s > {:.1}s)", age, config.max_base_age_s);
                false
            } else {
                true
            }
        })
    }

    fn perform_spp_fallback_update(config: &EngineConfig, state: &mut RtkState, pos: Coordinate) {
        tracing::warn!("valid_base IS NONE! Falling back to SPP!");
        tracing::debug!("RTK base missing or stale. Falling back to SPP update.");
        let z_diff = pos.vector - state.position.vector;
        let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());
        
        let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
        h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);
        
        let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
        r_mat.fill_diagonal(900.0);

        if let Err(e) = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, config.spp_consistency_threshold_m, None, config.mode.is_tightly_coupled() && state.ins_aligned, &config.tuning) {
            tracing::debug!("SPP Fallback update failed: {:?}", e);
        }
    }


    pub fn process_rtk(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {
        let spp_res = self.init_spp_state(rover_obs)?;
        let spp_pos = spp_res.as_ref().map(|s| s.position);
        let spp_cdt = spp_res.as_ref().map(|s| s.cdt).unwrap_or(0.0);
        let spp_state_ref = spp_res.as_ref();

        // Carrier-smooth rover pseudoranges before any state access
        let mut rover_smoothed = rover_obs.clone();
        self.hatch_filter.smooth_epoch(&mut rover_smoothed);

        let dt = rover_obs.time - self.current_state.as_ref().ok_or(EngineError::StateDisappeared)?.time;
        
        if let Some(state) = &self.current_state {
            if !state.ins_aligned {
                self.imu_buffer.clear();
            }
        }
        
        let had_imu_data = !self.imu_buffer.is_empty();
        
        if let Some(state) = &self.current_state {
            if state.epoch_count == 59 || state.epoch_count == 60 || state.epoch_count == 61 {
                tracing::info!("BEFORE predict_state({}): Pos={:.2?} Vel={:.2?} AccelBias={:.5?} GyroBias={:.5?}", dt, state.position.vector.as_slice(), state.velocity.as_slice(), state.accel_bias.as_slice(), state.gyro_bias.as_slice());
            }
        }
        
        self.predict_state(dt);
        
        if let Some(state) = &self.current_state {
            if state.epoch_count == 59 || state.epoch_count == 60 || state.epoch_count == 61 {
                tracing::info!("AFTER predict_state({}): Pos={:.2?} Vel={:.2?} AccelBias={:.5?} GyroBias={:.5?}", dt, state.position.vector.as_slice(), state.velocity.as_slice(), state.accel_bias.as_slice(), state.gyro_bias.as_slice());
            }
        }
        
        let state = self.current_state.as_mut().ok_or(EngineError::StateDisappeared)?;
        state.is_reset = false;
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;
        
        if state.epoch_count % 100 == 0 {
            tracing::info!("RTK INS State: Pos={:.2?} Vel={:.2?} AccelBias={:.5?} GyroBias={:.5?}", state.position.vector.as_slice(), state.velocity.as_slice(), state.accel_bias.as_slice(), state.gyro_bias.as_slice());
        }

        Self::check_covariance_divergence(state, spp_pos, spp_state_ref, !self.config.mode.is_ppp());
        Self::evaluate_gnss_only_coasting(&self.config, state, spp_pos, spp_state_ref, had_imu_data);

        let valid_base = Self::filter_valid_base(&self.config, &rover_smoothed, base_obs);

        if let Some(base) = valid_base {
            state.epoch_count += 1;
            tracing::warn!("valid_base IS SOME! incrementing epoch_count to {}", state.epoch_count);
            let mut base_coord = if let Some(base_pos_arr) = self.config.base_position {
                Coordinate::new(Vector3::new(base_pos_arr[0], base_pos_arr[1], base_pos_arr[2]), Datum::WGS84, Frame::ECEF, rover_obs.time)
            } else { return Err(EngineError::MissingBasePosition); };
            
            if let Some(helmert) = &self.config.base_datum_transform {
                let obs_epoch = rover_obs.time.to_fractional_year();
                let transformed_vec = helmert.transform(base_coord.vector, obs_epoch);
                base_coord = Coordinate::new(transformed_vec, Datum::WGS84, Frame::ECEF, rover_obs.time);
            }
            
            let matched_obs = match_observations(&rover_smoothed, base, &self.ephemerides);
            let epoch_num = state.epoch_count;
            if epoch_num % 100 == 0 { tracing::info!("Epoch {}: Matched {} satellites, {} ambiguities tracked", epoch_num, matched_obs.len(), state.ambiguity_keys.len()); }

            if matched_obs.len() >= 5 {
                process_rtk_update(&self.config, &self.ephemerides, &self.imu_history, state, &rover_smoothed, base, &matched_obs, &base_coord, &mut self.innovation_tracker, spp_pos, spp_state_ref);
            } else {
                tracing::warn!("Not enough valid measurements for EKF update. Riding through outage.");
                state.consecutive_rejections += 1;
            }
        } else if let Some(pos) = spp_pos {
            Self::perform_spp_fallback_update(&self.config, state, pos);
        }
        
        Self::apply_nhc_updates(&self.config, &self.imu_history, state);
        
        self.attempt_kinematic_alignment();

        if let Some(state) = &self.current_state { self.state_history.push(state.clone()); }
        self.obs_history.push((rover_obs.clone(), base_obs.cloned()));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
    }
}

/// RTK measurement update — extracted as a free function to avoid borrow conflicts.
fn process_rtk_update<'a>(
    config: &EngineConfig, ephemerides: &[gneiss_core::ephemeris::Ephemeris], imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>], 
    state: &mut RtkState, rover_obs: &EpochObs, base_obs: &'a EpochObs, matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
    base_coord: &Coordinate, tracker: &mut crate::engine::adaptive::InnovationTracker, spp_pos: Option<Coordinate>, spp_state_ref: Option<&crate::spp::SppState>,
) {
    crate::engine::ambiguity::manage_ambiguities_and_slips(state, config, matched_obs, ephemerides, base_coord, rover_obs.time, base_obs.time);
    
    let current_epoch = state.epoch_count as u32;
    for (r_obs, _) in matched_obs {
        if r_obs.cp_l1.is_some() { state.last_observed.insert((r_obs.sat, 1), current_epoch); }
        if r_obs.cp_l2.is_some() { state.last_observed.insert((r_obs.sat, 2), current_epoch); }
    }

    let omega_ib_b = if let Some(imu_buf) = imu_history.last() {
        if let Some(last_imu) = imu_buf.last() { last_imu.gyro - state.gyro_bias } else { nalgebra::Vector3::zeros() }
    } else { nalgebra::Vector3::zeros() };
    
    let omega_ie_e = nalgebra::Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
    let r_e_b = state.attitude.to_rotation_matrix().transpose();
    let omega_b = omega_ib_b - r_e_b * omega_ie_e;

    let lever_arm = if state.ins_aligned {
        Vector3::from_column_slice(&config.imu_to_antenna_lever_arm)
    } else {
        Vector3::zeros()
    };

    let env = crate::engine::measurement::MeasurementEnvironment {
        ephemerides,
        base_coord,
        base_time: base_obs.time,
        lever_arm,
        omega_b,
        tuning: &config.tuning,
    };

    let pr_thresh = config.chi_square_pr_threshold;
    let cp_thresh = config.chi_square_cp_threshold;

    if let Some(m) = crate::engine::measurement::build_measurement_model(
        state, matched_obs, &env, pr_thresh, cp_thresh
    ) {
        let crate::engine::measurement::EkfMeasurementMatrices { z: z_safe, h: h_safe, r: mut r_safe, mt: type_safe } = m;

        apply_adaptive_r_scaling(tracker, state, &z_safe, &h_safe, &mut r_safe, &type_safe);

        let type_stripped: Vec<_> = type_safe.iter().map(|&(s, t, _)| (s, t)).collect();
        
        if crate::engine::updater::update(state, &z_safe, &h_safe, &r_safe, pr_thresh, Some(&type_stripped), config.mode.is_tightly_coupled() && state.ins_aligned, &config.tuning).is_err() { 
            state.consecutive_rejections += 1;
            tracing::warn!("GNSS EKF rejected for {} epochs.", state.consecutive_rejections);
            
            let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
            let max_rejections = if state.ins_aligned { 500 } else { 3 };
            if state.consecutive_rejections >= max_rejections || (state.ins_aligned && pos_var > 900.0) {
                tracing::warn!("Divergence detected after {} epochs (pos_var {:.2}): resetting EKF to SPP fallback.", state.consecutive_rejections, pos_var);
                if let Some(pos) = spp_pos {
                    state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
                    state.consecutive_rejections = 0;
                }
            }
        } else {
            let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
            if pos_var > 900.0 {
                tracing::warn!("Position variance too large ({:.2}): resetting EKF to SPP fallback.", pos_var);
                if let Some(pos) = spp_pos {
                    state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
                }
            }
            
            state.consecutive_rejections = 0;
            match state.resolve_ambiguities(ephemerides, config.lambda_min_subset, config.ar_min_epoch_count, config.ar_min_lock, config.lambda_min_ratio, config.ar_ffrt_prob) {
                Ok((fixed_state, _da, _q_fixed, _ratio, _subset_size)) => {
                    tracing::debug!("Integer ambiguities resolved!");
                    state.fixed_state = Some(Box::new(fixed_state));
                }
                Err(e) => {
                    tracing::debug!("AR Failed: {}", e);
                    state.fixed_state = None;
                }
            }
        }
    } else {
        state.consecutive_rejections += 1;
        tracing::warn!("GNSS EKF rejected for {} epochs (measurement model empty).", state.consecutive_rejections);
        
        let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
        let max_rejections = if state.ins_aligned { 500 } else { 3 };
        if state.consecutive_rejections >= max_rejections || (state.ins_aligned && pos_var > 900.0) {
            tracing::warn!("Divergence detected after {} epochs (pos_var {:.2}): resetting EKF to SPP fallback.", state.consecutive_rejections, pos_var);
            if let Some(pos) = spp_pos {
                state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
                state.consecutive_rejections = 0;
            }
        }
    }
    tracing::warn!("End of RTK loop, epoch_count = {}", state.epoch_count);
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
}
