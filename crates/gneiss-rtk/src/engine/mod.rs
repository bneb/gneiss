pub mod matcher;
pub mod predictor;
pub mod updater;
pub mod measurement;
pub mod ppp;
pub mod ambiguity;

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
            chi_square_pr_threshold: 15.0,
            chi_square_cp_threshold: 1000000.0,
            nominal_snr_dbhz: 40.0,
            dynamics_model: DynamicsModel::Automotive,
            doppler_slip_threshold_cycles: 5.0,
            max_reject_count: 3,
            max_base_age_s: 5.0,
            spp_consistency_threshold_m: 50.0,
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
            let enable_imu = matches!(self.config.mode, EngineMode::SppIns | EngineMode::RtkIns | EngineMode::PppIns | EngineMode::RtkInsLooselyCoupled | EngineMode::SppInsLooselyCoupled);
            let imu_data = if enable_imu { &self.imu_buffer[..] } else { &[] };
            crate::engine::predictor::predict(state, dt, self.config.dynamics_model, imu_data);
            
            state.predicted_position = Some(state.position.clone());
            state.predicted_velocity = Some(state.velocity.clone());
            state.predicted_attitude = Some(state.attitude.clone());
            state.predicted_accel_bias = Some(state.accel_bias.clone());
            state.predicted_gyro_bias = Some(state.gyro_bias.clone());
        }
        self.imu_history.push(self.imu_buffer.clone());
        self.imu_buffer.clear();
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
            EngineMode::Ppp | EngineMode::PppIns => crate::engine::ppp::process_ppp(self, &filtered_rover).err(),
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

        if crate::engine::updater::update_loosely_coupled(state, gnss_state, self.config.imu_to_antenna_lever_arm.into(), omega_b).is_err() {
            state.consecutive_rejections += 1;
            if state.consecutive_rejections > 15 {
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

        if crate::engine::updater::update_loosely_coupled(state, &gnss_state, lever_arm, omega_b).is_err() {
            state.consecutive_rejections += 1;
            if state.consecutive_rejections > 15 {
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
        match crate::spp::compute_spp(rover_obs, &self.ephemerides, Some(&gneiss_core::atmosphere::KlobucharParams::default()), &crate::spp::SppConfig::default(), None) {
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

        // EKF update for SPP
        if let Some(pos) = spp_pos {
            let z_diff = pos.vector - state.position.vector;
            let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());
            
            let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
            h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);
            
            let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
            // Using a tighter variance of 9.0 (3m std dev) forces the INS to track the clean SPP positions
            r_mat.fill_diagonal(9.0);

            if crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, self.config.chi_square_pr_threshold, None).is_err() {
                state.consecutive_rejections += 1;
                if state.consecutive_rejections > 15 {
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
                    let _ = crate::nhc::apply_zupt(state, zupt_var);
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
                    let _ = crate::nhc::apply_nhc(state, 0.1, 0.1, &self.config.imu_mounting_angles, &self.config.imu_to_nhc_lever_arm, &omega_b);
                }
            }
        }

        if let Some(state) = &self.current_state { self.state_history.push(state.clone()); }
        self.obs_history.push((rover_obs.clone(), None));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
    }

    pub fn process_rtk(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {

        let mut spp_pos = None;
        let mut spp_cdt = 0.0;
        if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, Some(&gneiss_core::atmosphere::KlobucharParams::default()), &crate::spp::SppConfig::default(), None) {
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
        self.predict_state(dt);
        
        let state = self.current_state.as_mut().ok_or(EngineError::StateDisappeared)?;
        state.is_reset = false;
        state.time = rover_obs.time;
        state.position.epoch = rover_obs.time;

        if matches!(self.config.mode, EngineMode::Rtk | EngineMode::Ppp) {
            if let Some(pos) = spp_pos {
                state.position = pos;
                // Reset positional covariance to 100.0 m^2 and decouple from ambiguities
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
                
                // Also decouple velocity
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
        
        if matches!(self.config.mode, EngineMode::Rtk | EngineMode::Ppp) {
            // Decouple position to model it as white noise for standard Kinematic mode
            for i in 0..3 {
                for j in 0..state.covariance.ncols() {
                    if i == j {
                        state.covariance[(i, j)] = 10000.0;
                    } else {
                        state.covariance[(i, j)] = 0.0;
                        state.covariance[(j, i)] = 0.0;
                    }
                }
            }
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
            
            if state.epoch_count < 5 || state.epoch_count % 100 == 0 {
                tracing::info!("Epoch {}: Matched {} satellites", state.epoch_count, matched_obs.len());
            }

            if matched_obs.len() >= 5 {
                crate::engine::ambiguity::manage_ambiguities_and_slips(state, &self.config, &matched_obs, &self.ephemerides, &base_coord, rover_obs.time, base.time);

                    let omega_b = if let Some(imu_buf) = self.imu_history.last() {
                        if let Some(last_imu) = imu_buf.last() {
                            last_imu.gyro - state.gyro_bias
                        } else {
                            nalgebra::Vector3::zeros()
                        }
                    } else {
                        nalgebra::Vector3::zeros()
                    };

                    if let Some((z_safe, h_safe, r_safe, type_safe)) = crate::engine::measurement::build_measurement_model(
                        state,
                        &matched_obs,
                        &self.ephemerides,
                        &base_coord,
                        rover_obs.time,
                        base.time,
                        Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm),
                        omega_b,
                        self.config.chi_square_pr_threshold,
                        self.config.chi_square_cp_threshold,
                    ) {
                    
                    if crate::engine::updater::update(state, &z_safe, &h_safe, &r_safe, self.config.chi_square_pr_threshold, Some(&type_safe)).is_err() { 
                        state.consecutive_rejections += 1;
                        if state.consecutive_rejections > 15 {
                            if let Some(pos) = spp_pos {
                                tracing::warn!("GNSS EKF rejected for {} epochs. Hard resetting INS to SPP.", state.consecutive_rejections);
                                state.position = pos;
                                state.velocity = nalgebra::Vector3::zeros();
                                state.accel_bias = nalgebra::Vector3::zeros();
                                state.gyro_bias = nalgebra::Vector3::zeros();
                                // Preserve attitude as it is far better than identity
                                state.rcv_clk_bias = spp_cdt;
                                state.rcv_clk_drift = 0.0;
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
                                state.covariance[(15, 15)] = 1e6;
                                state.is_reset = true;
                                state.consecutive_rejections = 0;
                            } else {
                                tracing::warn!("GNSS EKF rejected, but SPP is unavailable. Riding through outage.");
                            }
                        } else {
                            tracing::warn!("GNSS EKF update rejected. Riding through outage via INS dead-reckoning.");
                        }
                    } else {
                        state.consecutive_rejections = 0;
                        if let Err(e) = state.resolve_ambiguities(&self.ephemerides, self.config.lambda_min_subset) {
                            tracing::debug!("AR Failed: {}", e);
                        }
                    }
                    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
                } else {
                    tracing::warn!("Not enough valid measurements for EKF update.");
                    if let Some(pos) = spp_pos {
                        tracing::warn!("Resetting EKF position and velocity to SPP due to insufficient measurements. Preserving attitude and biases.");
                        state.position = pos;
                        state.velocity = nalgebra::Vector3::zeros();
                        state.clear_ambiguities();
                        for i in 0..6 {
                            state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
                        }
                    }
                    state.is_reset = true;
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

            if let Err(e) = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, 15.0, None) {
                tracing::debug!("SPP Fallback update failed: {:?}", e);
            }
        }
        
        let is_ins = matches!(self.config.mode, EngineMode::SppIns | EngineMode::RtkIns | EngineMode::PppIns);
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
                    if let Err(e) = crate::nhc::apply_zupt(state, zupt_var) {
                        tracing::debug!("ZUPT failed: {}", e);
                    }
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
                    if let Err(e) = crate::nhc::apply_nhc(state, 0.1, 0.1, &self.config.imu_mounting_angles, &self.config.imu_to_nhc_lever_arm, &omega_b) {
                        tracing::debug!("NHC failed: {}", e);
                    }
                }
            }
        }

        if let Some(state) = &self.current_state { self.state_history.push(state.clone()); }
        self.obs_history.push((rover_obs.clone(), base_obs.cloned()));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
        }

    pub fn run_combined_ppk(&mut self) -> Result<Vec<RtkState>, EngineError> {
        let n_epochs = self.state_history.len();
        if n_epochs == 0 { return Err(EngineError::NoObservations); }
        
        let mut smoothed_states = self.state_history.clone();
        
        // Backward pass
        for k in (0..n_epochs - 1).rev() {
            if smoothed_states[k+1].is_reset {
                tracing::debug!("Epoch {} was reset. Breaking smoothing chain at k={}", k+1, k);
                continue;
            }
            // Extract matrices generated during the forward prediction from k to k+1
            // Note: These are stored in state_k1 during predict_state
            let phi_k = match &smoothed_states[k+1].core_phi {
                Some(p) => p.clone(),
                None => continue,
            };
            let p_pred_k1 = match &smoothed_states[k+1].full_p_predict {
                Some(p) => p.clone(),
                None => continue,
            };
            let x_pred_k1 = match &smoothed_states[k+1].full_x_predict {
                Some(x) => x.clone(),
                None => continue,
            };
            
            let state_k = &smoothed_states[k]; // x_{k|k}
            let state_k1 = &smoothed_states[k+1]; // x_{k+1|N}

            // Only smooth the first 15 states (kinematic + IMU biases) + matched ambiguities to avoid clock jump numerical instability
            let core_size = crate::filter::CORE_STATE_SIZE;
            
            // Find matching ambiguities
            let mut matched_k_indices = Vec::new();
            let mut matched_k1_indices = Vec::new();
            for (i, key_k) in state_k.ambiguity_keys.iter().enumerate() {
                if let Some(j) = state_k1.ambiguity_keys.iter().position(|k| k == key_k) {
                    let cov_idx = crate::filter::CORE_STATE_SIZE + j;
                    if state_k1.covariance[(cov_idx, cov_idx)] < 10.0 {
                        matched_k_indices.push(crate::filter::CORE_STATE_SIZE + i);
                        matched_k1_indices.push(cov_idx);
                    }
                }
            }
            
            let smooth_len = core_size + matched_k_indices.len();
            
            // Build index mapping for state k
            let mut idx_k = Vec::with_capacity(smooth_len);
            for i in 0..core_size { idx_k.push(i); }
            idx_k.extend(matched_k_indices.clone());
            
            // Build index mapping for state k+1
            let mut idx_k1 = Vec::with_capacity(smooth_len);
            for i in 0..core_size { idx_k1.push(i); }
            idx_k1.extend(&matched_k1_indices);

            let mut x_k1_n = DVector::zeros(smooth_len);
            x_k1_n.rows_mut(0, 3).copy_from(&state_k1.position.vector);
            x_k1_n.rows_mut(3, 3).copy_from(&state_k1.velocity);
            x_k1_n.rows_mut(9, 3).copy_from(&state_k1.accel_bias);
            x_k1_n.rows_mut(12, 3).copy_from(&state_k1.gyro_bias);
            x_k1_n[15] = state_k1.rcv_clk_bias;
            x_k1_n[16] = state_k1.rcv_clk_drift;
            x_k1_n[17] = state_k1.zwd;
            for i in 0..matched_k1_indices.len() {
                x_k1_n[core_size + i] = state_k1.ambiguities[matched_k1_indices[i] - crate::filter::CORE_STATE_SIZE];
            }
            
            let mut p_k1_n = DMatrix::zeros(smooth_len, smooth_len);
            for i in 0..smooth_len {
                for j in 0..smooth_len {
                    p_k1_n[(i, j)] = state_k1.covariance[(idx_k1[i], idx_k1[j])];
                }
            }
            
            let mut p_k = DMatrix::zeros(smooth_len, smooth_len);
            for i in 0..smooth_len {
                for j in 0..smooth_len {
                    p_k[(i, j)] = state_k.covariance[(idx_k[i], idx_k[j])];
                }
            }
            
            let mut p_pred_k1_sub = DMatrix::zeros(smooth_len, smooth_len);
            for i in 0..smooth_len {
                for j in 0..smooth_len {
                    p_pred_k1_sub[(i, j)] = p_pred_k1[(idx_k[i], idx_k[j])];
                }
            }
            
            let mut phi_k_sub = DMatrix::zeros(smooth_len, smooth_len);
            for i in 0..core_size {
                for j in 0..core_size {
                    phi_k_sub[(i, j)] = phi_k[(i, j)];
                }
            }
            for i in core_size..smooth_len {
                phi_k_sub[(i, i)] = 1.0;
            }
            
            let mut active_indices = Vec::with_capacity(smooth_len);
            if p_pred_k1_sub.iter().any(|&x| x.is_nan() || x.is_infinite() || x.abs() > 1e10) || p_k.iter().any(|&x| x.is_nan() || x.is_infinite() || x.abs() > 1e10) {
                tracing::debug!("RTS Smoothing skipped at k={} due to non-finite covariance", k);
                continue;
            }
            
            for i in 0..smooth_len {
                if p_pred_k1_sub[(i, i)] > 1e-12 {
                    active_indices.push(i);
                }
            }

            let mut x_pred_k1_sub = DVector::zeros(smooth_len);
            for i in 0..smooth_len {
                x_pred_k1_sub[i] = x_pred_k1[idx_k[i]];
            }

            let m = active_indices.len();
            let p_pred_inv_opt = if m == smooth_len {
                let reg = DMatrix::identity(smooth_len, smooth_len) * 1e-9;
                (p_pred_k1_sub.clone() + reg).try_inverse()
            } else if m > 0 {
                let mut p_active = DMatrix::zeros(m, m);
                for (i, &r) in active_indices.iter().enumerate() {
                    for (j, &c) in active_indices.iter().enumerate() {
                        p_active[(i, j)] = p_pred_k1_sub[(r, c)];
                    }
                }
                
                let reg = DMatrix::identity(m, m) * 1e-9;
                let inv_active_opt = (p_active + reg).try_inverse();
                    
                inv_active_opt.map(|inv_active| {
                    let mut inv_full = DMatrix::zeros(smooth_len, smooth_len);
                    for (i, &r) in active_indices.iter().enumerate() {
                        for (j, &c) in active_indices.iter().enumerate() {
                            inv_full[(r, c)] = inv_active[(i, j)];
                        }
                    }
                    inv_full
                })
            } else {
                None
            };
            
            let p_pred_inv = match p_pred_inv_opt {
                Some(inv) => inv,
                None => {
                    tracing::debug!("RTS Smoothing failed to invert P_pred at k={}", k);
                    continue;
                }
            };
            
            // Smoother gain C_k = P_{k|k} * Phi_k^T * P_{k+1|k}^{-1}
            let c_k = &p_k * phi_k_sub.transpose() * p_pred_inv;
            
            // Reconstruct x_{k|k}
            let mut x_k = DVector::zeros(smooth_len);
            x_k.rows_mut(0, 3).copy_from(&state_k.position.vector);
            x_k.rows_mut(3, 3).copy_from(&state_k.velocity);
            if core_size > 6 {
                x_k.rows_mut(9, 3).copy_from(&state_k.accel_bias);
                x_k.rows_mut(12, 3).copy_from(&state_k.gyro_bias);
            }
            for i in 0..matched_k_indices.len() {
                x_k[core_size + i] = state_k.ambiguities[matched_k_indices[i] - crate::filter::CORE_STATE_SIZE];
            }
            
            let correction = &c_k * (x_k1_n.clone() - x_pred_k1_sub.clone());
            let pos_corr_norm = correction.fixed_rows::<3>(0).norm();
            if k % 1000 == 0 {
                tracing::info!("Smoother at k={}. pos_corr: {:.3}m", k, pos_corr_norm);
                tracing::info!("P_k pos var: {:.3e}, vel var: {:.3e}", p_k[(0, 0)], p_k[(3, 3)]);
                tracing::info!("P_pred pos var: {:.3e}, vel var: {:.3e}", p_pred_k1_sub[(0, 0)], p_pred_k1_sub[(3, 3)]);
                tracing::info!("C_k pos,pos: {:.3e}, pos,vel: {:.3e}", c_k[(0, 0)], c_k[(0, 3)]);
            }
            if pos_corr_norm > 10.0 {
                tracing::warn!("Smoother huge pos correction: {:.1}m at k={}. x_k1_n pos: {:.1}, x_pred_k1 pos: {:.1}", pos_corr_norm, k, x_k1_n.fixed_rows::<3>(0).norm(), x_pred_k1_sub.fixed_rows::<3>(0).norm());
                tracing::warn!("P_k: {}", p_k);
                tracing::warn!("P_pred: {}", p_pred_k1_sub);
                tracing::warn!("C_k: {}", c_k);
                tracing::warn!("diff: {}", (x_k1_n.clone() - x_pred_k1_sub.clone()));
            }
            
            let x_k_n = x_k + correction;
            let p_k_n = p_k + &c_k * (p_k1_n - p_pred_k1_sub) * c_k.transpose();
            
            // Update smoothed_states[k]
            let s_k_mut = &mut smoothed_states[k];
            s_k_mut.position.vector = x_k_n.fixed_rows::<3>(0).into_owned();
            s_k_mut.velocity = x_k_n.fixed_rows::<3>(3).into_owned();
            s_k_mut.accel_bias = x_k_n.fixed_rows::<3>(9).into_owned();
            s_k_mut.gyro_bias = x_k_n.fixed_rows::<3>(12).into_owned();
            s_k_mut.rcv_clk_bias = x_k_n[15];
            s_k_mut.rcv_clk_drift = x_k_n[16];
            s_k_mut.zwd = x_k_n[17];
            
            for i in 0..matched_k_indices.len() {
                let amb_idx = matched_k_indices[i] - crate::filter::CORE_STATE_SIZE;
                s_k_mut.ambiguities[amb_idx] = x_k_n[core_size + i];
            }
            for i in 0..smooth_len {
                for j in 0..smooth_len {
                    s_k_mut.covariance[(idx_k[i], idx_k[j])] = p_k_n[(i, j)];
                }
            }
        }
        Ok(smoothed_states)
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
        let initial_pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
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
        let mut state0 = RtkState::new(time0, pos0, 1.0);
        
        let pos1 = Coordinate::new(Vector3::new(12.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time1);
        let mut state1 = RtkState::new(time1, pos1, 0.5);
        state1.is_fixed = true;
        
        // Mock prediction values from 0 to 1
        state1.core_phi = Some(DMatrix::identity(18, 18));
        state1.full_p_predict = Some(DMatrix::identity(18, 18) * 1.5);
        let mut x_pred = DVector::zeros(18);
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
mod tests_measurement;
mod tests_predictor;
pub mod jacobian_verify;
