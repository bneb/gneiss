use super::ProcessingEngine;
use crate::engine::{EngineError, EngineMode};
use crate::filter::RtkState;
use gneiss_core::obs::EpochObs;

impl ProcessingEngine {
    pub fn process_rtk_loosely_coupled(
        &mut self,
        rover_obs: &EpochObs,
        base_obs: Option<&EpochObs>,
    ) -> Result<&RtkState, EngineError> {
        let prev_state = self.current_state.clone(); // preserve INS state
        let prev_config = self.config.mode;

        // 1. Run standard RTK purely for GNSS (swap states)
        self.current_state = self.gnss_only_state.take();
        self.config.mode = EngineMode::Rtk; // temporarily act as pure RTK
        let gnss_res_cloned = self.process_rtk(rover_obs, base_obs).cloned();
        self.gnss_only_state = self.current_state.take(); // save GNSS ambiguity state
        self.config.mode = prev_config; // restore
        self.current_state = prev_state; // restore INS state

        let gnss_state = gnss_res_cloned?;

        if self.current_state.is_none() {
            // Seed INS filter with first RTK fix
            let mut ins_state = RtkState::new(rover_obs.time, gnss_state.position, 0.1);
            ins_state.velocity = gnss_state.velocity;
            ins_state.is_fixed = gnss_state.is_fixed;
            self.current_state = Some(ins_state);
        }

        let dt = rover_obs.time.tow - self.current_state.as_ref().unwrap().time.tow;
        self.predict_state(dt);
        let state = self.current_state.as_mut().unwrap();
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        let lever_arm = nalgebra::Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm);
        let r_b_e = state.attitude.to_rotation_matrix();
        let omega_ie_e =
            nalgebra::Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
        let omega_b = if let Some(imu_buf) = self.imu_history.last() {
            if let Some(last_imu) = imu_buf.last() {
                last_imu.gyro - state.gyro_bias - r_b_e.transpose() * omega_ie_e
            } else {
                nalgebra::Vector3::zeros()
            }
        } else {
            nalgebra::Vector3::zeros()
        };

        if crate::engine::updater::update_loosely_coupled(
            state,
            &gnss_state,
            lever_arm,
            omega_b,
            &self.config.tuning,
        )
        .is_err()
        {
            state.consecutive_rejections += 1;
            if state.consecutive_rejections > 5 {
                tracing::warn!(
                    "Loose coupling rejected for {} epochs. Hard resetting INS to GNSS.",
                    state.consecutive_rejections
                );
                state.position = gnss_state.position;
                state.velocity = gnss_state.velocity;
                state.accel_bias = nalgebra::Vector3::zeros();
                state.gyro_bias = nalgebra::Vector3::zeros();
                // Preserve attitude as gnss_state.attitude is likely identity

                state.covariance.fill(0.0);
                let n = crate::filter::CORE_STATE_SIZE;
                for i in 0..6 {
                    state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
                }
                for i in 6..n {
                    state.covariance[(i, i)] = 1e-4;
                }
                state.is_reset = true;
                state.consecutive_rejections = 0;
            } else {
                tracing::warn!(
                    "Loose coupling update rejected. Riding through outage via INS dead-reckoning."
                );
            }
        } else {
            state.consecutive_rejections = 0;
        }

        state.is_fixed = gnss_state.is_fixed;

        // We already pushed history in the inner process_rtk if we let it, but actually process_rtk
        // pushes to state_history. To prevent duplicates or mixed histories, we should pop the last ones
        // or just rely on this wrapper for the *final* history.
        // Wait, `process_rtk` pushes to `self.state_history` and `self.obs_history`.
        // We will pop them here and push the true INS state.
        self.state_history.pop();
        self.obs_history.pop();

        self.state_history.push(state.clone());
        self.obs_history
            .push((rover_obs.clone(), base_obs.cloned()));

        Ok(self.current_state.as_ref().unwrap())
    }

    pub fn process_spp_loosely_coupled(
        &mut self,
        rover_obs: &EpochObs,
    ) -> Result<&RtkState, EngineError> {
        let prev_config = self.config.mode;
        self.config.mode = EngineMode::Spp;
        let mut gnss_res_cloned = self.process_spp(rover_obs).cloned()?;
        self.config.mode = prev_config;

        for i in 3..6 {
            gnss_res_cloned.covariance[(i, i)] = 1e6; // Ignore the zero velocity from pure SPP
        }

        self.gnss_only_state = Some(gnss_res_cloned);

        if self.current_state.is_none() {
            let state = self.gnss_only_state.as_ref().unwrap();
            let mut ins_state = RtkState::new(rover_obs.time, state.position, 0.1);
            ins_state.velocity = state.velocity;
            self.current_state = Some(ins_state);
        }

        let dt = rover_obs.time.tow - self.current_state.as_ref().unwrap().time.tow;
        self.predict_state(dt);
        let state = self.current_state.as_mut().unwrap();
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        let gnss_state = self.gnss_only_state.as_ref().unwrap();
        let r_b_e = state.attitude.to_rotation_matrix();
        let omega_ie_e =
            nalgebra::Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
        let omega_b = if let Some(imu_buf) = self.imu_history.last() {
            if let Some(last_imu) = imu_buf.last() {
                last_imu.gyro - state.gyro_bias - r_b_e.transpose() * omega_ie_e
            } else {
                nalgebra::Vector3::zeros()
            }
        } else {
            nalgebra::Vector3::zeros()
        };

        if crate::engine::updater::update_loosely_coupled(
            state,
            gnss_state,
            self.config.imu_to_antenna_lever_arm.into(),
            omega_b,
            &self.config.tuning,
        )
        .is_err()
        {
            state.consecutive_rejections += 1;
            if state.consecutive_rejections > 5 {
                tracing::warn!(
                    "SPP-INS EKF rejected for {} epochs. Hard resetting INS to SPP.",
                    state.consecutive_rejections
                );
                state.position = gnss_state.position;
                state.velocity = gnss_state.velocity; // Zero out diverged velocity
                state.accel_bias = nalgebra::Vector3::zeros(); // Biases might be corrupted, 0 is a safer prior
                state.gyro_bias = nalgebra::Vector3::zeros();
                // Preserve attitude as it is far better than identity
                state.covariance.fill(0.0);
                for i in 0..6 {
                    state.covariance[(i, i)] = if i < 3 { 10.0 } else { 1.0 };
                }
                let att_var = (1.0f64.to_radians()).powi(2);
                for i in 6..9 {
                    state.covariance[(i, i)] = att_var;
                }
                for i in 9..12 {
                    state.covariance[(i, i)] = 0.01;
                }
                for i in 12..15 {
                    state.covariance[(i, i)] = (0.1f64.to_radians()).powi(2);
                }
                state.consecutive_rejections = 0;
                state.is_reset = true;
            }
        } else {
            state.consecutive_rejections = 0;
            state.is_reset = false;
        }

        if state.is_reset || self.state_history.len().is_multiple_of(100) {
            tracing::info!(
                "INS State: AccelBias={:.5?} GyroBias={:.5?}",
                state.accel_bias.as_slice(),
                state.gyro_bias.as_slice()
            );
        }
        self.state_history.push(state.clone());
        self.obs_history.push((rover_obs.clone(), None));
        Ok(self.current_state.as_ref().unwrap())
    }

    pub fn predict_state(&mut self, dt: f64) {
        if let Some(state) = self.current_state.as_mut() {
            let enable_imu = matches!(
                self.config.mode,
                EngineMode::SppIns
                    | EngineMode::RtkIns
                    | EngineMode::PppIns
                    | EngineMode::RtkInsLooselyCoupled
                    | EngineMode::SppInsLooselyCoupled
                    | EngineMode::PppInsLooselyCoupled
                    | EngineMode::RtkInsIekf
            );
            let mut imu_data = if enable_imu {
                &self.imu_buffer[..]
            } else {
                &[]
            };
            if !state.ins_aligned {
                imu_data = &[];
            }
            crate::engine::predictor::predict(state, dt, &self.config, imu_data);

            state.predicted_position = Some(state.position);
            state.predicted_velocity = Some(state.velocity);
            state.predicted_attitude = Some(state.attitude);
            state.predicted_accel_bias = Some(state.accel_bias);
            state.predicted_gyro_bias = Some(state.gyro_bias);
        }
        self.imu_history.push(self.imu_buffer.clone());
        self.imu_buffer.clear();
    }
}
