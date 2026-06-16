use crate::engine::{EngineMode, EngineError};
use crate::filter::RtkState;
use gneiss_core::obs::EpochObs;
use gneiss_core::coords::Coordinate;
use super::ProcessingEngine;

impl ProcessingEngine {

        fn perform_spp_ekf_update(config: &crate::engine::EngineConfig, state: &mut RtkState, spp_pos: Option<Coordinate>, spp_cdt: f64) {
        if let Some(pos) = spp_pos {
            let z_diff = pos.vector - state.position.vector;
            let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());
            
            let mut rejected = false;
            if config.mode.is_tightly_coupled() && state.ins_aligned && z_diff.norm() > config.spp_consistency_threshold_m {
                rejected = true;
            } else {
                let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
                h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);
                
                let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
                // Using a tighter variance of 9.0 (3m std dev) forces the INS to track the clean SPP positions
                r_mat.fill_diagonal(9.0);

                if crate::engine::updater::update::<crate::engine::updater_math::LooseCoupling>(state, &z_vec, &h_mat, &r_mat, config.spp_consistency_threshold_m, None, &config.tuning).map_or(true, |v| v.0.len() < 3) {
                    rejected = true;
                }
            }

            if rejected {
                state.consecutive_rejections += 1;
                if state.consecutive_rejections > 5 {
                    tracing::warn!("SPP EKF rejected for {} epochs. Hard resetting INS to SPP.", state.consecutive_rejections);
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
                    for i in 6..9 { state.covariance[(i, i)] = att_var; }
                    for i in 9..12 { state.covariance[(i, i)] = 0.01; }
                    for i in 12..n {
                        state.covariance[(i, i)] = 1e-4;
                    }
                    if crate::filter::CORE_STATE_SIZE > 15 {
                        state.covariance[(15, 15)] = 1e6;
                    }
                    state.is_reset = true;
                    state.consecutive_rejections = 0;
                } else {
                    tracing::warn!("SPP EKF update rejected. Riding through outage via INS dead-reckoning.");
                }
            } else {
                state.consecutive_rejections = 0;
            }
        }
    }

    pub fn process_spp(&mut self, rover_obs: &EpochObs) -> Result<&RtkState, EngineError> {
        let spp_res_opt = crate::spp::compute_spp(rover_obs, &self.ephemerides, self.klobuchar_params.as_ref(), &crate::spp::SppConfig::default(), None).ok();
        let spp_pos = spp_res_opt.as_ref().map(|s| s.position);
        let spp_cdt = spp_res_opt.as_ref().map(|s| s.cdt).unwrap_or(0.0);

        if spp_res_opt.is_none() {
            tracing::warn!("Initial SPP compute failed");
        }

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

        let dt = rover_obs.time.tow - self.current_state.as_ref().ok_or(EngineError::StateDisappeared)?.time.tow ;
        self.predict_state(dt);
        
        if let Some(pos) = spp_pos {
            if matches!(self.config.mode, EngineMode::Spp) {
                // Pure SPP is an epoch-by-epoch solution. Do not filter.
                let state = self.current_state.as_mut().unwrap();
                state.time = rover_obs.time;
                state.position = pos;
                state.position.epoch = rover_obs.time;
                state.velocity = nalgebra::Vector3::zeros();
                tracing::info!("process_spp: SPP mode return, rcv_clk_bias = {}", state.rcv_clk_bias);
                self.state_history.push(state.clone());
                self.obs_history.push((rover_obs.clone(), None));
                return Ok(self.current_state.as_ref().unwrap());
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

        let state = self.current_state.as_mut().ok_or(EngineError::StateDisappeared)?;
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        Self::check_covariance_divergence(state, spp_pos, spp_res_opt.as_ref(), !self.config.mode.is_ppp());

        Self::perform_spp_ekf_update(&self.config, state, spp_pos, spp_cdt);

        if let Some(state) = self.current_state.as_mut() {
            Self::apply_nhc_updates(&self.config, &self.imu_history, state);
        }

        if let Some(state) = &self.current_state {
            tracing::info!("process_spp: Returning state, rcv_clk_bias = {}", state.rcv_clk_bias);
            self.state_history.push(state.clone());
        }
        self.obs_history.push((rover_obs.clone(), None));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
    }

}
