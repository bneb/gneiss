use crate::engine::{EngineMode, EngineError, EngineConfig};
use crate::filter::RtkState;
use gneiss_core::obs::EpochObs;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::sat::SatelliteId;
use nalgebra::{Vector3, DVector, DMatrix};
use super::ProcessingEngine;
use crate::engine::matcher::match_observations;
use crate::engine::updater_math::{CouplingStrategy, TightCoupling, LooseCoupling};

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
        let freq = match mtype { 2 => 2, _ => 1 };
        
        let mut snr = None;
        if let Some((rov_obs, _)) = matched_obs.iter().find(|(rov, _)| rov.sat == sat) {
            snr = Some(rov_obs.snr);
        }
        
        let scale = tracker.update_and_scale(sat, freq, z[i], s_ii, snr);
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
            tracing::debug!("epoch_count = {}", state.epoch_count);
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
        tracing::debug!("No valid base data. Falling back to SPP.");
        tracing::debug!("RTK base missing or stale. Falling back to SPP update.");
        let z_diff = pos.vector - state.position.vector;
        let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());
        
        let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
        h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);
        
        let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
        r_mat.fill_diagonal(900.0);

        if let Err(e) = crate::engine::updater::update::<LooseCoupling>(state, &z_vec, &h_mat, &r_mat, config.spp_consistency_threshold_m, None, &config.tuning) {
            tracing::debug!("SPP Fallback update failed: {:?}", e);
        }
    }


    pub fn process_rtk(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {
        let spp_res = self.init_spp_state(rover_obs)?;
        let spp_pos = spp_res.as_ref().map(|s| s.position);
        let _spp_cdt = spp_res.as_ref().map(|s| s.cdt).unwrap_or(0.0);
        let spp_state_ref = spp_res.as_ref();

        // Carrier-smooth rover pseudoranges before any state access
        let mut rover_smoothed = rover_obs.clone();
        self.hatch_filter.smooth_epoch(&mut rover_smoothed);

        let dt = rover_obs.time.tow - self.current_state.as_ref().ok_or(EngineError::StateDisappeared)?.time.tow ;
        
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
            tracing::debug!("valid_base found, epoch_count = {}", state.epoch_count);
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
                let mut gnn_variances = std::collections::HashMap::new();
                if let Some(gnn) = &self.gnn_raim {
                    let pos_apc = state.position.vector;
                    let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
                    gnn_variances = crate::engine::ml::gnn_raim::evaluate_gnn_raim(
                        gnn, &matched_obs, rov_llh, pos_apc, &self.ephemerides, rover_obs.time
                    );
                }

                let ctx = RtkUpdateContext {
                    config: &self.config, ephemerides: &self.ephemerides, imu_history: &self.imu_history,
                    rover_obs: &rover_smoothed, base_obs: base, matched_obs: &matched_obs,
                    base_coord: &base_coord, spp_pos, spp_state_ref, gnn_variances,
                };
                process_rtk_update::<TightCoupling>(state, &mut self.innovation_tracker, &ctx);
            } else {
                tracing::warn!("Not enough valid measurements for EKF update. Riding through outage.");
                state.consecutive_rejections += 1;
            }
        } else if let Some(pos) = spp_pos {
            Self::perform_spp_fallback_update(&self.config, state, pos);
        }
        
        Self::apply_nhc_updates(&self.config, &self.imu_history, state);
        
        self.attempt_kinematic_alignment();

        if let Some(state) = &self.current_state { self.state_history.push(RtkState::clone(state)); }
        self.obs_history.push((rover_obs.clone(), base_obs.cloned()));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
    }
}

/// RTK measurement update — extracted as a free function to avoid borrow conflicts.
fn build_measurement_environment<'a>(
    config: &'a EngineConfig, imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>], 
    state: &RtkState, ephemerides: &'a [gneiss_core::ephemeris::Ephemeris], 
    base_coord: &'a Coordinate, base_time: gneiss_core::time::GpsTime
) -> crate::engine::measurement::MeasurementEnvironment<'a> {
    let omega_ib_b = if let Some(imu_buf) = imu_history.last() {
        if let Some(last_imu) = imu_buf.last() { last_imu.gyro - state.gyro_bias } else { nalgebra::Vector3::zeros() }
    } else { nalgebra::Vector3::zeros() };
    let omega_ie_e = nalgebra::Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
    let r_e_b = state.attitude.to_rotation_matrix().transpose();
    let lever_arm = if state.ins_aligned { Vector3::from_column_slice(&config.imu_to_antenna_lever_arm) } else { Vector3::zeros() };

    crate::engine::measurement::MeasurementEnvironment {
        ephemerides, base_coord, base_time, lever_arm, omega_b: omega_ib_b - r_e_b * omega_ie_e, tuning: &config.tuning,
        gnn_variances: std::collections::HashMap::new(),
    }
}

fn handle_ekf_rejection(state: &mut RtkState, config: &EngineConfig, spp_pos: Option<Coordinate>, spp_state_ref: Option<&crate::spp::SppState>, reason: &str) {
    state.consecutive_rejections += 1;
    tracing::warn!("GNSS EKF rejected for {} epochs ({}).", state.consecutive_rejections, reason);
    if state.ins_aligned {
        let inflate_factor = 1.0 + (state.consecutive_rejections as f64 * 0.02).min(0.5);
        for i in 0..6 { state.covariance[(i, i)] *= inflate_factor; }
        for i in 0..3 { state.covariance[(i, i)] += 1.0; }
        for i in 3..6 { state.covariance[(i, i)] += 0.1; }
    }
    let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
    if state.ins_aligned && pos_var > 10000.0 {
        tracing::warn!("Extreme divergence detected (pos_var {:.2}): resetting EKF to SPP fallback.", pos_var);
        if let Some(pos) = spp_pos { state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp()); state.consecutive_rejections = 0; }
    } else if !state.ins_aligned && state.consecutive_rejections >= 3 {
        tracing::warn!("Loosely coupled GNSS EKF rejected for 3 epochs: resetting to SPP fallback.");
        if let Some(pos) = spp_pos { state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp()); state.consecutive_rejections = 0; }
    }
}

fn handle_ekf_acceptance(state: &mut RtkState, config: &EngineConfig, ephemerides: &[gneiss_core::ephemeris::Ephemeris], spp_pos: Option<Coordinate>, spp_state_ref: Option<&crate::spp::SppState>) {
    let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
    if pos_var > 10000.0 {
        tracing::warn!("Position variance extremely large ({:.2}): resetting EKF to SPP fallback.", pos_var);
        if let Some(pos) = spp_pos { state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp()); }
    }
    state.consecutive_rejections = 0;
    if let Ok(res) = state.resolve_ambiguities(ephemerides, config) {
        let fixed_state = res.fixed_state;
        tracing::debug!("Integer ambiguities resolved!");
        state.fixed_state = Some(Box::new(fixed_state));
    } else {
        state.fixed_state = None;
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
}

fn process_rtk_update<C: CouplingStrategy>(
    state: &mut RtkState, tracker: &mut crate::engine::adaptive::InnovationTracker, ctx: &RtkUpdateContext
) {
    crate::engine::ambiguity::manage_ambiguities_and_slips(state, ctx.config, ctx.matched_obs, ctx.ephemerides, ctx.base_coord, ctx.rover_obs.time, ctx.base_obs.time);
    let current_epoch = state.epoch_count as u32;
    for (r_obs, _) in ctx.matched_obs {
        if r_obs.cp_l1.is_some() { state.last_observed.insert((r_obs.sat, 1), current_epoch); }
        if r_obs.cp_l2.is_some() { state.last_observed.insert((r_obs.sat, 2), current_epoch); }
    }

    let mut env = build_measurement_environment(ctx.config, ctx.imu_history, state, ctx.ephemerides, ctx.base_coord, ctx.base_obs.time);
    env.gnn_variances = ctx.gnn_variances.clone();
    let pr_thresh = ctx.config.chi_square_pr_threshold;

    if let Some(mut m) = crate::engine::measurement::build_measurement_model(state, ctx.matched_obs, &env, pr_thresh, ctx.config.chi_square_cp_threshold) {
        apply_adaptive_r_scaling(tracker, state, &m.z, &m.h, &mut m.r, &m.mt, ctx.matched_obs);
        
        let tuning = ctx.config.tuning.clone();

        let type_stripped: Vec<_> = m.mt.iter().map(|&(s, t, _)| (s, t)).collect();
        match crate::engine::updater::update::<C>(state, &m.z, &m.h, &m.r, pr_thresh, Some(&type_stripped), &tuning) {
            Err(_) => {
                handle_ekf_rejection(state, ctx.config, ctx.spp_pos, ctx.spp_state_ref, "failed chi-square");
            }
            Ok((valid_indices, dx)) => {
                if let Some(path) = &ctx.config.export_gnn_dataset_path {
                    let ephemerides = ctx.ephemerides;
                    let rov_llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
                    let pos_apc = state.position.vector + state.attitude * nalgebra::Vector3::from(ctx.config.imu_to_antenna_lever_arm);
                    
                    crate::engine::ml::dataset::export_epoch_to_csv(
                        path,
                        state.epoch_count as u32,
                        ctx.matched_obs,
                        rov_llh,
                        pos_apc,
                        ephemerides,
                        ctx.rover_obs.time,
                        &m.z,
                        &m.h,
                        &valid_indices,
                        &m.mt,
                        &dx,
                    );
                }
                handle_ekf_acceptance(state, ctx.config, ctx.ephemerides, ctx.spp_pos, ctx.spp_state_ref);
            }
        }
    } else {
        handle_ekf_rejection(state, ctx.config, ctx.spp_pos, ctx.spp_state_ref, "not enough measurements");
    }
    tracing::debug!("End of RTK loop, epoch_count = {}", state.epoch_count);
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
}
