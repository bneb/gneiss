pub mod ppp_fg;
pub mod tight_fg;
pub mod tcar;
pub mod processed_sat;
pub mod matcher;
pub mod predictor;
pub mod updater;
pub mod measurement;
pub mod ppp;
pub mod ambiguity;
pub mod config;
pub mod auto_tuner;
pub mod smoother;

use nalgebra::{Vector3, DMatrix, DVector};
use gneiss_core::obs::EpochObs;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::SatelliteId;
use crate::filter::RtkState;
use crate::engine::matcher::match_observations;

/// Defines the core positioning mode
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub enum EngineMode {
    Spp,
    SppIns,
    SppInsLooselyCoupled,
    #[default]
    Rtk,
    RtkIns,
    RtkInsLooselyCoupled,
    Ppp,
    PppIns,
    PppInsLooselyCoupled,
}

impl EngineMode {
    pub fn is_tightly_coupled(&self) -> bool {
        matches!(self, Self::SppIns | Self::RtkIns | Self::PppIns)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum DynamicsModel {
    Static,
    Pedestrian,
    Automotive,
    Marine,
    Airborne,
}

/// Configuration for the RTK processing engine.
#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
#[serde(default)]
pub struct EngineConfig {
    pub mode: EngineMode,
    pub initial_position: Option<[f64; 3]>,
    pub base_position: Option<[f64; 3]>,
    pub base_datum_transform: Option<gneiss_geodesy::helmert::HelmertParams>,
    pub imu_to_antenna_lever_arm: [f64; 3],
    pub imu_mounting_angles: Option<[f64; 3]>, // [Roll, Pitch, Yaw] in radians
    pub imu_to_nhc_lever_arm: [f64; 3], // [x, y, z] from IMU to NHC point in body frame
    pub enable_nhc: bool,
    pub enable_backward_smoothing: bool,
    pub lambda_min_ratio: f64,
    pub lambda_min_subset: usize,
    pub enabled_constellations: Option<Vec<gneiss_core::sat::Constellation>>,
    
    // Tuning Parameters
    pub raim_pseudorange_outlier_m: f64,
    pub chi_square_pr_threshold: f64,
    pub chi_square_cp_threshold: f64,
    pub nominal_snr_dbhz: f64,
    pub dynamics_model: DynamicsModel,
    pub doppler_slip_threshold_cycles: f64,
    pub max_reject_count: usize,
    pub max_base_age_s: f64,
    pub spp_consistency_threshold_m: f64,
    pub initial_ambiguity_variance: f64,
    pub ar_min_epoch_count: u32,
    pub ar_min_lock: u32,
    
    // Process Noise
    pub process_noise_cb: f64,
    pub process_noise_cd: f64,
    pub process_noise_zwd: f64,
    pub process_noise_amb_float: f64,
    pub process_noise_amb_fixed: f64,
    
    // External Tuning configuration
    pub tuning: crate::engine::config::EkfTuningConfig,
}

impl Default for EngineConfig {
    fn default() -> Self {
        Self {
            mode: EngineMode::Rtk,
            initial_position: None,
            base_position: None,
            base_datum_transform: None,
            imu_to_antenna_lever_arm: [0.0, 0.0, 0.0],
            imu_mounting_angles: None,
            imu_to_nhc_lever_arm: [0.0, 0.0, 0.0],
            enable_nhc: false,
            enable_backward_smoothing: false,
            lambda_min_ratio: 3.0,
            lambda_min_subset: 5,
            enabled_constellations: None,
            raim_pseudorange_outlier_m: 25.0,
            chi_square_pr_threshold: 3.0,
            chi_square_cp_threshold: 1000000.0,
            nominal_snr_dbhz: 40.0,
            dynamics_model: DynamicsModel::Automotive,
            doppler_slip_threshold_cycles: 5.0,
            max_reject_count: 3,
            max_base_age_s: 5.0,
            spp_consistency_threshold_m: 15.0,
            initial_ambiguity_variance: 10000.0,
            ar_min_epoch_count: 5,
            ar_min_lock: 3,
            process_noise_cb: 1e6,
            process_noise_cd: 1e4,
            process_noise_zwd: 1e-8,
            process_noise_amb_float: 1e-8,
            process_noise_amb_fixed: 1e-12,
            tuning: Default::default(),
        }
    }
}


#[derive(Debug, Clone, PartialEq)]
pub enum EngineError {
    NoObservations,
    InitialSppFailed,
    StateDisappeared,
    InsufficientSatellites,
    MissingBasePosition,
    GeodeticMismatch(&'static str),
}

impl std::fmt::Display for EngineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EngineError::NoObservations => write!(f, "No observations available"),
            EngineError::InitialSppFailed => write!(f, "Initial SPP failed"),
            EngineError::StateDisappeared => write!(f, "EKF state disappeared mid-execution"),
            EngineError::InsufficientSatellites => write!(f, "Insufficient satellites for EKF update"),
            EngineError::MissingBasePosition => write!(f, "Base station position must be provided"),
            EngineError::GeodeticMismatch(msg) => write!(f, "Geodetic Gatekeeper Failed: {}", msg),
        }
    }
}
impl std::error::Error for EngineError {}

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

    pub fn predict_state(&mut self, dt: f64) {
        if let Some(state) = self.current_state.as_mut() {
            let enable_imu = matches!(self.config.mode, EngineMode::SppIns | EngineMode::RtkIns | EngineMode::PppIns | EngineMode::RtkInsLooselyCoupled | EngineMode::SppInsLooselyCoupled | EngineMode::PppInsLooselyCoupled);
            let imu_data = if enable_imu { &self.imu_buffer[..] } else { &[] };
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

        let mut filtered_base_storage = None;
        if let Some(b) = base_obs {
            let mut clone = b.clone();
            if let Some(enabled) = &self.config.enabled_constellations {
                clone.satellites.retain(|s| enabled.contains(&s.sat.constellation));
            }
            filtered_base_storage = Some(clone);
        }
        let filtered_base = filtered_base_storage.as_ref();

        let err = match self.config.mode {
            EngineMode::Spp | EngineMode::SppIns => self.process_spp(&filtered_rover).err(),
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

    pub fn process_spp_loosely_coupled(&mut self, rover_obs: &EpochObs) -> Result<&RtkState, EngineError> {
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
        
        let dt = rover_obs.time - self.current_state.as_ref().unwrap().time;
        self.predict_state(dt);
        let state = self.current_state.as_mut().unwrap();
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        let gnss_state = self.gnss_only_state.as_ref().unwrap();
        let omega_b = if let Some(imu_buf) = self.imu_history.last() {
            if let Some(last_imu) = imu_buf.last() {
                last_imu.gyro - state.gyro_bias
            } else { nalgebra::Vector3::zeros() }
        } else { nalgebra::Vector3::zeros() };

        if crate::engine::updater::update_loosely_coupled(state, gnss_state, self.config.imu_to_antenna_lever_arm.into(), omega_b, &self.config.tuning).is_err() {
            state.consecutive_rejections += 1;
            if state.consecutive_rejections > 5 {
                tracing::warn!("SPP-INS EKF rejected for {} epochs. Hard resetting INS to SPP.", state.consecutive_rejections);
                state.position = gnss_state.position;
                state.velocity = gnss_state.velocity; // Zero out diverged velocity
                state.accel_bias = nalgebra::Vector3::zeros(); // Biases might be corrupted, 0 is a safer prior
                state.gyro_bias = nalgebra::Vector3::zeros();
                // Preserve attitude as it is far better than identity
                state.covariance.fill(0.0);
                for i in 0..6 { state.covariance[(i, i)] = if i < 3 { 10.0 } else { 1.0 }; }
                let att_var = (1.0f64.to_radians()).powi(2);
                for i in 6..9 { state.covariance[(i, i)] = att_var; }
                for i in 9..12 { state.covariance[(i, i)] = 0.01; }
                for i in 12..15 { state.covariance[(i, i)] = (0.1f64.to_radians()).powi(2); }
                state.consecutive_rejections = 0;
                state.is_reset = true;
            }
        } else {
            state.consecutive_rejections = 0;
            state.is_reset = false;
        }

        self.state_history.push(state.clone());
        self.obs_history.push((rover_obs.clone(), None));
        Ok(self.current_state.as_ref().unwrap())
    }

    pub fn process_rtk_loosely_coupled(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {
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

        let dt = rover_obs.time - self.current_state.as_ref().unwrap().time;
        self.predict_state(dt);
        let state = self.current_state.as_mut().unwrap();
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        let lever_arm = nalgebra::Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm);
        let omega_b = if let Some(imu_buf) = self.imu_history.last() {
            if let Some(last_imu) = imu_buf.last() {
                last_imu.gyro - state.gyro_bias
            } else {
                nalgebra::Vector3::zeros()
            }
        } else {
            nalgebra::Vector3::zeros()
        };

        if crate::engine::updater::update_loosely_coupled(state, &gnss_state, lever_arm, omega_b, &self.config.tuning).is_err() {
            state.consecutive_rejections += 1;
            if state.consecutive_rejections > 5 {
                tracing::warn!("Loose coupling rejected for {} epochs. Hard resetting INS to GNSS.", state.consecutive_rejections);
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
                tracing::warn!("Loose coupling update rejected. Riding through outage via INS dead-reckoning.");
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
        self.obs_history.push((rover_obs.clone(), base_obs.cloned()));
        
        Ok(self.current_state.as_ref().unwrap())
    }

    pub fn process_spp(&mut self, rover_obs: &EpochObs) -> Result<&RtkState, EngineError> {
        let mut spp_pos = None;
        let mut spp_cdt = 0.0;
        match crate::spp::compute_spp(rover_obs, &self.ephemerides, self.klobuchar_params.as_ref(), &crate::spp::SppConfig::default(), None) {
            Ok(spp_res) => {
                spp_pos = Some(spp_res.position);
                spp_cdt = spp_res.cdt;
            },
            Err(e) => {
                tracing::warn!("Initial SPP compute failed: {:?}", e);
            }
        }

        if self.current_state.is_none() {
            if let Some(pos) = spp_pos {
                let mut new_state = RtkState::new(rover_obs.time, pos, 100.0);
                if crate::filter::CORE_STATE_SIZE > 15 {
                    new_state.rcv_clk_bias = spp_cdt;
                }
                self.current_state = Some(new_state);
            } else {
                return Err(EngineError::InitialSppFailed);
            }
        }

        let dt = rover_obs.time - self.current_state.as_ref().ok_or(EngineError::StateDisappeared)?.time;
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

        // Covariance integrity monitor for SPP-INS
        let pos_var_max = state.covariance[(0,0)].max(state.covariance[(1,1)]).max(state.covariance[(2,2)]);
        if pos_var_max > 10000.0 {
            if let Some(pos) = spp_pos {
                tracing::warn!("SPP-INS position variance {:.0} m² exceeds integrity limit. Resetting.", pos_var_max);
                state.position = pos;
                state.velocity = nalgebra::Vector3::zeros();
                state.decouple_position();
                for i in 0..3 { state.covariance[(i, i)] = 100.0; }
                for i in 3..6 { state.covariance[(i, i)] = 10.0; }
                state.is_reset = true;
                state.consecutive_rejections = 0;
            }
        }

        // EKF update for SPP
        if let Some(pos) = spp_pos {
            let z_diff = pos.vector - state.position.vector;
            let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());
            
            let mut rejected = false;
            if self.config.mode.is_tightly_coupled() && z_diff.norm() > self.config.spp_consistency_threshold_m {
                rejected = true;
            } else {
                let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
                h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);
                
                let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
                // Using a tighter variance of 9.0 (3m std dev) forces the INS to track the clean SPP positions
                r_mat.fill_diagonal(9.0);

                if crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, self.config.spp_consistency_threshold_m, None, self.config.mode.is_tightly_coupled(), &self.config.tuning).map_or(true, |v| v.len() < 3) {
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

        if self.config.enable_nhc {
            if let Some(state) = self.current_state.as_mut() {
                let mut is_stationary = false;
                let mut accel_var = 1.0;
                if let Some(imu_buf) = self.imu_history.last() {
                    if imu_buf.len() > 10 {
                        let mut sum_a = nalgebra::Vector3::zeros();
                        let mut sum_g = nalgebra::Vector3::zeros();
                        for m in imu_buf {
                            sum_a += m.accel;
                            sum_g += m.gyro;
                        }
                        let mean_a = sum_a / (imu_buf.len() as f64);
                        let mean_g = sum_g / (imu_buf.len() as f64);
                        
                        let mut var_a = 0.0;
                        let mut var_g = 0.0;
                        for m in imu_buf {
                            var_a += (m.accel - mean_a).norm_squared();
                            var_g += (m.gyro - mean_g).norm_squared();
                        }
                        var_a /= imu_buf.len() as f64;
                        var_g /= imu_buf.len() as f64;
                        
                        if var_a < 0.05 && var_g < 0.005 {
                            is_stationary = true;
                        }
                        accel_var = var_a.max(0.001);
                    }
                }
                
                if !is_stationary && state.velocity.norm() < 0.05 {
                    is_stationary = true;
                }

                if is_stationary {
                    let zupt_var = (accel_var * 0.1).clamp(0.001, 0.1).sqrt();
                    let _ = crate::nhc::apply_zupt(state, zupt_var, &self.config.tuning);
                } else {
                    let omega_b = if let Some(imu_buf) = self.imu_history.last() {
                        if let Some(last_imu) = imu_buf.last() {
                            last_imu.gyro - state.gyro_bias
                        } else {
                            nalgebra::Vector3::zeros()
                        }
                    } else {
                        nalgebra::Vector3::zeros()
                    };
                    let _ = crate::nhc::apply_nhc(state, 0.1, 0.1, &self.config.imu_mounting_angles, &self.config.imu_to_nhc_lever_arm, &omega_b, &self.config.tuning);
                }
            }
        }

        if let Some(state) = &self.current_state {
            tracing::info!("process_spp: Returning state, rcv_clk_bias = {}", state.rcv_clk_bias);
            self.state_history.push(state.clone());
        }
        self.obs_history.push((rover_obs.clone(), None));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
    }

    pub fn process_rtk(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {

        let mut spp_pos = None;
        let mut spp_cdt = 0.0;
        if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, self.klobuchar_params.as_ref(), &crate::spp::SppConfig::default(), None) {
            spp_pos = Some(spp_res.position);
            spp_cdt = spp_res.cdt;
        } else {
            tracing::warn!("compute_spp failed in process_rtk!");
        }

        if self.current_state.is_none() {
            if let Some(pos) = spp_pos {
                self.current_state = Some(RtkState::new(rover_obs.time, pos, 100.0));
            } else {
                return Err(EngineError::InitialSppFailed);
            }
        }

        let dt = rover_obs.time - self.current_state.as_ref().ok_or(EngineError::StateDisappeared)?.time;
        
        // Track whether INS prediction actually used IMU data
        let had_imu_data = !self.imu_buffer.is_empty();
        self.predict_state(dt);
        
        let state = self.current_state.as_mut().ok_or(EngineError::StateDisappeared)?;
        state.is_reset = false;
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        // Covariance integrity monitor: if position variance has blown up,
        // the filter is diverging. Force a hard reset to SPP.
        let pos_var_max = state.covariance[(0,0)].max(state.covariance[(1,1)]).max(state.covariance[(2,2)]);
        if pos_var_max > 10000.0 {
            if let Some(pos) = spp_pos {
                tracing::warn!("Position variance {:.0} m² exceeds integrity limit. Resetting to SPP.", pos_var_max);
                state.position = pos;
                state.velocity = nalgebra::Vector3::zeros();
                state.rcv_clk_bias = spp_cdt;
                state.rcv_clk_drift = 0.0;
                state.clear_ambiguities();
                state.decouple_position();
                for i in 0..3 { state.covariance[(i, i)] = 100.0; }
                for i in 3..6 { state.covariance[(i, i)] = 10.0; }
                state.covariance[(15, 15)] = 1e6;
                state.is_reset = true;
                state.consecutive_rejections = 0;
            }
        }

        // For GNSS-only modes (Rtk/Ppp without INS), OR for INS modes that had no IMU data:
        // Use velocity-propagated position as the primary estimate. This provides temporal
        // filtering like RTKLIB — the position covariance from the prediction model reflects
        // the actual uncertainty, and carrier phase measurements refine it.
        // Only fall back to SPP when the prediction has clearly diverged.
        let use_gnss_only_seed = matches!(self.config.mode, EngineMode::Rtk | EngineMode::Ppp)
            || (matches!(self.config.mode, EngineMode::RtkIns | EngineMode::PppIns | EngineMode::SppIns | EngineMode::RtkInsLooselyCoupled | EngineMode::SppInsLooselyCoupled | EngineMode::PppInsLooselyCoupled) && !had_imu_data);
        
        if use_gnss_only_seed {
            let need_spp_reset = if let Some(_pos) = spp_pos {
                // Also reset for early epochs when the filter hasn't converged yet
                state.epoch_count < 3
            } else {
                false // No SPP available; keep the predicted position
            };

            if need_spp_reset {
                if let Some(pos) = spp_pos {
                    tracing::warn!("Resetting EKF state to SPP due to divergence/startup.");
                    state.position = pos;
                    state.velocity = nalgebra::Vector3::zeros();
                    state.clear_ambiguities();
                    state.is_reset = true;
                    state.consecutive_rejections = 0;
                    // Reset position covariance when falling back to SPP
                    for i in 0..3 {
                        for j in 0..state.covariance.ncols() {
                            if i == j {
                                state.covariance[(i, j)] = 100.0;
                            } else {
                                state.covariance[(i, j)] = 0.0;
                                state.covariance[(j, i)] = 0.0;
                            }
                        }
                    }
                    for i in 3..6 {
                        for j in 0..state.covariance.ncols() {
                            if i == j {
                                state.covariance[(i, j)] = 10.0;
                            } else {
                                state.covariance[(i, j)] = 0.0;
                                state.covariance[(j, i)] = 0.0;
                            }
                        }
                    }
                }
            }
            // Otherwise: keep velocity-propagated position and covariance from predict_state()
            // The process noise in the predictor handles position uncertainty growth
        }

        let valid_base = base_obs.filter(|b| {
            let age = (rover_obs.time.tow - b.time.tow).abs();
            if age > self.config.max_base_age_s {
                tracing::trace!("Base observation rejected due to Age of Differential ({:.1}s > {:.1}s)", age, self.config.max_base_age_s);
                false
            } else {
                true
            }
        });

        if let Some(base) = valid_base {
            state.epoch_count += 1;
            let mut base_coord = if let Some(base_pos_arr) = self.config.base_position {
                // Base station is physically static, so its coordinate is valid at the rover's epoch time
                Coordinate::new(Vector3::new(base_pos_arr[0], base_pos_arr[1], base_pos_arr[2]), Datum::WGS84, Frame::ECEF, rover_obs.time)
            } else { return Err(EngineError::MissingBasePosition); };
            
            // GEODETIC PIPELINE: Route base station through Helmert Transformation if configured
            if let Some(helmert) = &self.config.base_datum_transform {
                let obs_epoch = rover_obs.time.to_fractional_year();
                let transformed_vec = helmert.transform(base_coord.vector, obs_epoch);
                base_coord = Coordinate::new(transformed_vec, Datum::WGS84, Frame::ECEF, rover_obs.time);
            }
            
            let matched_obs = match_observations(rover_obs, base, &self.ephemerides);
            
            let epoch_num = state.epoch_count;
            if epoch_num % 100 == 0 {
                tracing::info!("Epoch {}: Matched {} satellites, {} ambiguities tracked", epoch_num, matched_obs.len(), state.ambiguity_keys.len());
            }

            if matched_obs.len() >= 5 {
                crate::engine::ambiguity::manage_ambiguities_and_slips(state, &self.config, &matched_obs, &self.ephemerides, &base_coord, rover_obs.time, base.time);
                
                // Update last_observed for all active ambiguities that were matched in this epoch
                let current_epoch = state.epoch_count as u32;
                for (r_obs, _) in &matched_obs {
                    if r_obs.cp_l1.is_some() {
                        state.last_observed.insert((r_obs.sat, 1), current_epoch);
                    }
                    if r_obs.cp_l2.is_some() {
                        state.last_observed.insert((r_obs.sat, 2), current_epoch);
                    }
                }

                    let omega_b = if let Some(imu_buf) = self.imu_history.last() {
                        if let Some(last_imu) = imu_buf.last() {
                            last_imu.gyro - state.gyro_bias
                        } else {
                            nalgebra::Vector3::zeros()
                        }
                    } else {
                        nalgebra::Vector3::zeros()
                    };

                    let env = crate::engine::measurement::MeasurementEnvironment {
                        ephemerides: &self.ephemerides,
                        base_coord: &base_coord,
                        base_time: base.time,
                        lever_arm: Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm),
                        omega_b,
                        tuning: &self.config.tuning,
                    };

                    if let Some((z_safe, h_safe, r_safe, type_safe)) = crate::engine::measurement::build_measurement_model(
                        state,
                        &matched_obs,
                        &env,
                        self.config.chi_square_pr_threshold,
                        self.config.chi_square_cp_threshold,
                    ) {
                    
                    if crate::engine::updater::update(state, &z_safe, &h_safe, &r_safe, self.config.chi_square_pr_threshold, {
                            let type_stripped: Vec<_> = type_safe.iter().map(|&(s, t, _)| (s, t)).collect();
                            Some(type_stripped)
                        }.as_deref(), self.config.mode.is_tightly_coupled(), &self.config.tuning).is_err() { 
                        state.consecutive_rejections += 1;
                        tracing::warn!("GNSS EKF rejected for {} epochs.", state.consecutive_rejections);
                    } else {
                        state.consecutive_rejections = 0;
                        match state.resolve_ambiguities(&self.ephemerides, self.config.lambda_min_subset, self.config.ar_min_epoch_count, self.config.ar_min_lock, self.config.lambda_min_ratio) {
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
                    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
                } else {
                    tracing::warn!("Not enough valid measurements for EKF update. Riding through outage.");
                    state.consecutive_rejections += 1;
                }
            }
        } else if let Some(pos) = spp_pos {
            // RTK Coasting Fallback to SPP
            tracing::debug!("RTK base missing or stale. Falling back to SPP update.");
            let z_diff = pos.vector - state.position.vector;
            let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());
            
            let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
            h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);
            
            let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
            r_mat.fill_diagonal(25.0);

            if let Err(e) = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, self.config.spp_consistency_threshold_m, None, self.config.mode.is_tightly_coupled(), &self.config.tuning) {
                tracing::debug!("SPP Fallback update failed: {:?}", e);
            }
        }
        
        let is_ins = matches!(self.config.mode, EngineMode::SppIns | EngineMode::RtkIns | EngineMode::PppIns | EngineMode::RtkInsLooselyCoupled | EngineMode::SppInsLooselyCoupled | EngineMode::PppInsLooselyCoupled);
        if self.config.enable_nhc && is_ins {
            if let Some(state) = self.current_state.as_mut() {
                // Determine if stationary for ZUPT
                let mut is_stationary = false;
                let mut accel_var = 1.0;
                if let Some(imu_buf) = self.imu_history.last() {
                    if imu_buf.len() > 10 {
                        let mut sum_a = nalgebra::Vector3::zeros();
                        let mut sum_g = nalgebra::Vector3::zeros();
                        for m in imu_buf {
                            sum_a += m.accel;
                            sum_g += m.gyro;
                        }
                        let mean_a = sum_a / (imu_buf.len() as f64);
                        let mean_g = sum_g / (imu_buf.len() as f64);
                        
                        let mut var_a = 0.0;
                        let mut var_g = 0.0;
                        for m in imu_buf {
                            var_a += (m.accel - mean_a).norm_squared();
                            var_g += (m.gyro - mean_g).norm_squared();
                        }
                        var_a /= imu_buf.len() as f64;
                        var_g /= imu_buf.len() as f64;
                        
                        if var_a < 0.05 && var_g < 0.005 {
                            is_stationary = true;
                        }
                        accel_var = var_a.max(0.001);
                    }
                }
                
                if !is_stationary && state.velocity.norm() < 0.05 {
                    is_stationary = true;
                }

                if is_stationary {
                    let zupt_var = (accel_var * 0.1).clamp(0.001, 0.1).sqrt();
                    let _ = crate::nhc::apply_zupt(state, zupt_var, &self.config.tuning);
                } else {
                    let omega_b = if let Some(imu_buf) = self.imu_history.last() {
                        if let Some(last_imu) = imu_buf.last() {
                            last_imu.gyro - state.gyro_bias
                        } else {
                            nalgebra::Vector3::zeros()
                        }
                    } else {
                        nalgebra::Vector3::zeros()
                    };
                    let _ = crate::nhc::apply_nhc(state, 0.1, 0.1, &self.config.imu_mounting_angles, &self.config.imu_to_nhc_lever_arm, &omega_b, &self.config.tuning);
                }
            }
        }

        if let Some(state) = &self.current_state { self.state_history.push(state.clone()); }
        self.obs_history.push((rover_obs.clone(), base_obs.cloned()));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
        }

    pub fn run_combined_ppk(&mut self) -> Result<Vec<RtkState>, EngineError> {
        crate::engine::smoother::run_combined_ppk(self)
    }
}

pub fn snr_scale(snr: f64) -> f64 { gneiss_core::variance::snr_variance_scale(snr, 45.0) }

#[cfg(test)]
mod tests {
    use super::*;
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
mod tests_predictor;
pub mod jacobian_verify;
mod tests_updater;
