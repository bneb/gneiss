use crate::engine::{EngineConfig, EngineMode, EngineError, DynamicsModel};
use crate::filter::RtkState;
use gneiss_core::obs::{EpochObs, ObsType};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::SatelliteId;
use gneiss_core::coords::Coordinate;

mod rtk;
mod spp;
mod ins;
mod ppk;


pub struct ProcessingEngine {
    pub config: EngineConfig,
    pub klobuchar_params: Option<gneiss_core::atmosphere::KlobucharParams>,
    pub current_state: Option<RtkState>,
    pub gnss_only_state: Option<RtkState>, // For loosely coupled modes
    pub ephemerides: Vec<Ephemeris>,
    pub state_history: Vec<RtkState>,
    pub obs_history: Vec<(EpochObs, Option<EpochObs>)>,
    pub imu_buffer: Vec<gneiss_core::imu::ImuMeasurement>,
    pub imu_history: Vec<Vec<gneiss_core::imu::ImuMeasurement>>,
    pub ref_sat: Option<SatelliteId>,
    pub hatch_filter: crate::hatch::HatchFilter,
    pub innovation_tracker: crate::engine::adaptive::InnovationTracker,
    pub sp3_epochs: Vec<gneiss_parsers::sp3::Sp3Epoch>,
    pub clk_data: Option<gneiss_parsers::rinex_clk::RinexClock>,
    pub antex: Option<gneiss_parsers::antex::AntexDatabase>,
    pub dcbs: std::collections::HashMap<(gneiss_core::sat::SatelliteId, String), f64>,
}

impl ProcessingEngine {
    pub fn new(mut config: EngineConfig) -> Self {
        if matches!(config.dynamics_model, DynamicsModel::Automotive | DynamicsModel::Pedestrian) {
            config.enable_nhc = true;
        }
        Self {
            config,
            klobuchar_params: None,
            current_state: None,
            gnss_only_state: None,
            ephemerides: Vec::new(),
            state_history: Vec::new(),
            obs_history: Vec::new(),
            imu_buffer: Vec::new(),
            imu_history: Vec::new(),
            ref_sat: None,
            hatch_filter: crate::hatch::HatchFilter::default(),
            innovation_tracker: crate::engine::adaptive::InnovationTracker::default(),
            sp3_epochs: Vec::new(),
            clk_data: None,
            antex: None,
            dcbs: std::collections::HashMap::new(),
        }
    }

    pub fn add_imu_measurement(&mut self, mut meas: gneiss_core::imu::ImuMeasurement) {
        // Apply mounting calibration if provided
        if let Some(angles) = self.config.imu_mounting_angles {
            let r_m_v = nalgebra::Rotation3::from_euler_angles(angles[0], angles[1], angles[2]);
            meas.accel = r_m_v * meas.accel;
            meas.gyro = r_m_v * meas.gyro;
        }
        self.imu_buffer.push(meas);
    }

    pub fn add_ephemeris(&mut self, eph: gneiss_core::ephemeris::Ephemeris) {
        self.ephemerides.push(eph);
    }

    pub(crate) fn check_covariance_divergence(state: &mut RtkState, spp_pos: Option<Coordinate>, spp_state_ref: Option<&crate::spp::SppState>, init_isbs: bool) {
        let pos_var_max = state.covariance[(0,0)].max(state.covariance[(1,1)]).max(state.covariance[(2,2)]);
        if pos_var_max > 10000.0 {
            if let Some(pos) = spp_pos {
                tracing::warn!("Position variance {:.0} m² exceeds integrity limit. Resetting to SPP.", pos_var_max);
                state.reset_to_spp(pos, spp_state_ref, init_isbs);
            }
        }
    }

    pub(crate) fn apply_nhc_updates(config: &EngineConfig, imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>], state: &mut RtkState) {
        let is_ins = matches!(config.mode, EngineMode::SppIns | EngineMode::RtkIns | EngineMode::PppIns | EngineMode::RtkInsLooselyCoupled | EngineMode::SppInsLooselyCoupled | EngineMode::PppInsLooselyCoupled);
        tracing::trace!("apply_nhc_updates: enable_nhc={}, is_ins={}, ins_aligned={}", config.enable_nhc, is_ins, state.ins_aligned);
        if config.enable_nhc && is_ins && state.ins_aligned {
            let mut is_stationary = false;
            let mut accel_var = 1.0;
            if let Some(imu_buf) = imu_history.last() {
                if imu_buf.len() > 10 {
                    let mut sum_a = nalgebra::Vector3::zeros();
                    let mut sum_g = nalgebra::Vector3::zeros();
                    for m in imu_buf { sum_a += m.accel; sum_g += m.gyro; }
                    let mean_a = sum_a / (imu_buf.len() as f64);
                    let mean_g = sum_g / (imu_buf.len() as f64);
                    
                    let mut var_a = 0.0f64;
                    let mut var_g = 0.0f64;
                    for m in imu_buf {
                        var_a += (m.accel - mean_a).norm_squared();
                        var_g += (m.gyro - mean_g).norm_squared();
                    }
                    var_a /= imu_buf.len() as f64;
                    var_g /= imu_buf.len() as f64;
                    
                    if var_a < 0.05 && var_g < 0.005 && state.velocity.norm() < 1.0 { is_stationary = true; }
                    accel_var = var_a.max(0.001f64);
                }
            }
            
            if !is_stationary && state.velocity.norm() < 0.05 { is_stationary = true; }

            if is_stationary {
                let zupt_var = (accel_var * 0.1).clamp(0.001f64, 0.1f64).sqrt();
                let _ = crate::nhc::apply_zupt(state, zupt_var, &config.tuning);
            } else {
                let omega_b = if let Some(imu_buf) = imu_history.last() {
                    if let Some(last_imu) = imu_buf.last() { last_imu.gyro - state.gyro_bias } else { nalgebra::Vector3::<f64>::zeros() }
                } else { nalgebra::Vector3::<f64>::zeros() };
                let _ = crate::nhc::apply_nhc(state, 0.1, 0.1, &config.imu_to_nhc_lever_arm, &omega_b, &config.tuning);
            }
        }
    }

    pub fn reset_for_multipass(&mut self) {
        if let Some(first) = self.state_history.first() {
            let mut reset_state = first.clone();
            reset_state.predicted_position = None;
            reset_state.predicted_velocity = None;
            reset_state.predicted_attitude = None;
            reset_state.predicted_accel_bias = None;
            reset_state.predicted_gyro_bias = None;
            self.current_state = Some(reset_state);
        } else {
            self.current_state = None;
        }
        self.gnss_only_state = None;
        self.state_history.clear();
        self.obs_history.clear();
        self.imu_buffer.clear();
        self.imu_history.clear();
        self.ref_sat = None;
    }

    pub fn process_epoch(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {
        let mut filtered_rover = rover_obs.clone();
        if let Some(enabled) = &self.config.enabled_constellations {
            filtered_rover.satellites.retain(|s| enabled.contains(&s.sat.constellation));
        }
        if self.config.min_snr_dbhz > 0.0 {
            filtered_rover.satellites.retain(|s| {
                let snr = s.observations.iter().find(|o| o.code.obs_type == ObsType::Snr && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(25.0);
                snr >= self.config.min_snr_dbhz
            });
        }

        let mut filtered_base_storage = None;
        if let Some(b) = base_obs {
            let mut clone = b.clone();
            if let Some(enabled) = &self.config.enabled_constellations {
                clone.satellites.retain(|s| enabled.contains(&s.sat.constellation));
            }
            if self.config.min_snr_dbhz > 0.0 {
                clone.satellites.retain(|s| {
                    let snr = s.observations.iter().find(|o| o.code.obs_type == ObsType::Snr && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(25.0);
                    snr >= self.config.min_snr_dbhz
                });
            }
            filtered_base_storage = Some(clone);
        }
        let filtered_base = filtered_base_storage.as_ref();

        let err = match self.config.mode {
            EngineMode::Spp => self.process_spp(&filtered_rover).err(),
            EngineMode::SppIns => crate::engine::spp_tight::process_spp_tightly_coupled(self, &filtered_rover).err(),
            EngineMode::SppInsLooselyCoupled => self.process_spp_loosely_coupled(&filtered_rover).err(),
            EngineMode::Rtk | EngineMode::RtkIns => self.process_rtk(&filtered_rover, filtered_base).err(),
            EngineMode::RtkInsLooselyCoupled => self.process_rtk_loosely_coupled(&filtered_rover, filtered_base).err(),
            EngineMode::Ppp | EngineMode::PppIns | EngineMode::PppInsLooselyCoupled => crate::engine::ppp::process_ppp(self, &filtered_rover).err(),
        };
        if let Some(e) = err {
            if let EngineError::StateDisappeared = e {
                tracing::warn!("EKF unrecoverable divergence. Resetting state...");
                self.current_state = None;
            } else {
                tracing::warn!("Epoch processing failed: {:?}. Preserving state for next epoch.", e);
            }
            Err(e)
        } else {
            Ok(self.current_state.as_ref().unwrap())
        }
    }

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
}

pub fn snr_scale(snr: f64) -> f64 { gneiss_core::variance::snr_variance_scale(snr, 45.0, 10.0) }

#[cfg(test)]
mod tests {
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use nalgebra::Vector3;

    use super::*;
    use nalgebra::{DMatrix, DVector};
    use gneiss_core::time::GpsTime;
    #[test]
    fn test_engine_detects_movement() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let initial_pos = Coordinate::new(Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        engine.current_state = Some(RtkState::new(GpsTime::new(0, 0.0), initial_pos, 1.0));
        engine.predict_state(1.0);
        let pos1 = engine.current_state.as_ref().unwrap().position.vector;
        engine.current_state.as_mut().unwrap().velocity = Vector3::new(10.0, 0.0, 0.0);
        engine.predict_state(1.0);
        let pos2 = engine.current_state.as_ref().unwrap().position.vector;
        assert!(pos2.x > pos1.x);
    }

    #[test]
    fn test_rts_smoother_basic() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time0 = GpsTime::new(2000, 0.0);
        let time1 = GpsTime::new(2000, 1.0);
        
        let pos0 = Coordinate::new(Vector3::new(10.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time0);
        let state0 = RtkState::new(time0, pos0, 1.0);
        
        let pos1 = Coordinate::new(Vector3::new(12.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time1);
        let mut state1 = RtkState::new(time1, pos1, 0.5);
        state1.is_fixed = true;
        
        // Mock prediction values from 0 to 1
        let core_size = crate::filter::CORE_STATE_SIZE;
        state1.core_phi = Some(DMatrix::identity(core_size, core_size));
        state1.full_p_predict = Some(DMatrix::identity(core_size, core_size) * 1.5);
        let mut x_pred = DVector::zeros(core_size);
        x_pred[0] = 10.0; // Assume velocity was 0, so predicted pos is 10
        state1.full_x_predict = Some(x_pred);
        
        engine.state_history.push(state0);
        engine.state_history.push(state1);
        
        let smoothed = engine.run_combined_ppk().expect("Should run RTS smoother");
        
        assert_eq!(smoothed.len(), 2);
        assert!(!smoothed[0].is_fixed, "Fix should NOT propagate backwards anymore");
        
        // P_{0|1} = P_0 + C_0 (P_{1|1} - P_{1|0}) C_0^T
        // C_0 = P_0 * Phi^T * P_{1|0}^-1 = 1.0 * 1 * (1.5)^-1 = 0.666
        // P_{0|1} = 1.0 + 0.666 * (0.5 - 1.5) * 0.666 = 1.0 - 0.444 = 0.555
        let p_0_1 = smoothed[0].covariance[(0,0)];
        assert!((p_0_1 - 0.555555).abs() < 1e-4, "Covariance mismatch: {}", p_0_1);
        
        // x_{0|1} = x_0 + C_0 (x_1 - x_pred)
        // x_{0|1} = 10.0 + 0.666 * (12.0 - 10.0) = 11.333
        let x_0_1 = smoothed[0].position.vector.x;
        assert!((x_0_1 - 11.333333).abs() < 1e-4, "Position mismatch: {}", x_0_1);
    }
}


impl ProcessingEngine {
    pub(crate) fn attempt_kinematic_alignment(&mut self) {
        let state = if let Some(s) = &mut self.current_state { s } else { return };
        if state.ins_aligned || !self.config.mode.is_tightly_coupled() { return; }
        
        let speed = state.velocity.norm();
        if speed > 3.0 && self.state_history.len() >= 5 && self.state_history.iter().rev().take(5).all(|s| s.velocity.norm() > 3.0) {
            let llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
            let ecef_to_ned = gneiss_core::coords::ecef_to_ned_matrix(llh);
            let v_ned = ecef_to_ned * state.velocity;
            let yaw = f64::atan2(v_ned.y, v_ned.x);
            let rot_veh_to_ned = nalgebra::Rotation3::from_euler_angles(0.0, 0.0, yaw);
            
            let ned_to_ecef = ecef_to_ned.transpose();
            // Since imu measurements are rotated to vehicle frame in add_imu_measurement,
            // state.attitude should represent the rotation from VEHICLE to ECEF.
            let rot_mat = nalgebra::Rotation3::from_matrix_unchecked(ned_to_ecef * rot_veh_to_ned.matrix());
            state.attitude = nalgebra::UnitQuaternion::from_rotation_matrix(&rot_mat);
            
            // Shift state position from antenna phase center to IMU center
            let r_e_v = state.attitude.to_rotation_matrix();
            let lever_arm = nalgebra::Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm);
            state.position.vector -= r_e_v * lever_arm;
            
            state.ins_aligned = true;
            tracing::info!("Kinematic alignment successful! Speed: {:.2} m/s, Veh Yaw: {:.2} deg", speed, yaw.to_degrees());
            
            // Set attitude covariance higher because we assume 0 roll/pitch and yaw is based on noisy GNSS velocity
            for i in 6..9 { state.covariance[(i, i)] = (15.0f64.to_radians()).powi(2); }
            self.imu_buffer.clear();
        }
    }
}
