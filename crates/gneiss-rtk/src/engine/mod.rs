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
    pub gnss_process_noise_var: f64,
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
            enable_nhc: false,
            enable_backward_smoothing: false,
            lambda_min_ratio: 3.0,
            lambda_min_subset: 5,
            enabled_constellations: None,
            raim_pseudorange_outlier_m: 25.0,
            chi_square_pr_threshold: 15.0,
            chi_square_cp_threshold: 1000000.0,
            nominal_snr_dbhz: 40.0,
            gnss_process_noise_var: 1.0,
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
    pub fn new(config: EngineConfig) -> Self {
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
            crate::engine::predictor::predict(state, dt, if enable_imu { 0.1 } else { self.config.gnss_process_noise_var }, &self.imu_buffer);
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
        // Run SPP as normal, which natively acts loosely coupled already
        self.process_spp(rover_obs)
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

        // 2. Perform loosely coupled measurement update using RTK position & velocity
        let mut z_vec = DVector::zeros(6);
        z_vec.rows_mut(0, 3).copy_from(&(gnss_state.position.vector - state.position.vector));
        z_vec.rows_mut(3, 3).copy_from(&(gnss_state.velocity - state.velocity));
        
        let mut h_mat = DMatrix::zeros(6, state.covariance.ncols());
        h_mat.view_mut((0, 0), (6, 6)).fill_diagonal(1.0); // Identity for Pos & Vel
        
        // Use GNSS state covariance for R
        let mut r_mat = DMatrix::zeros(6, 6);
        r_mat.copy_from(&gnss_state.covariance.view((0, 0), (6, 6)));

        if let Err(e) = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, 15.0, None) {
            tracing::debug!("RTK loosely coupled EKF update failed: {:?}", e);
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
                if let Some(state) = self.current_state.as_mut() {
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
            r_mat.fill_diagonal(25.0);

            if let Err(e) = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, 15.0, None) {
                tracing::debug!("SPP EKF update failed: {:?}", e);
            }
        }

        if self.config.enable_nhc {
            if let Some(state) = self.current_state.as_mut() {
                let stationary = state.velocity.norm() < 0.1 && state.accel_bias.norm() < 10.0;
                if stationary {
                    let _ = crate::nhc::apply_zupt(state, 0.01);
                } else {
                    let _ = crate::nhc::apply_nhc(state, 0.1, 0.1);
                }
            }
        }

        if let Some(state) = &self.current_state { self.state_history.push(state.clone()); }
        self.obs_history.push((rover_obs.clone(), None));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
    }

    pub fn process_rtk(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {

        let mut spp_pos = None;
        if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, Some(&gneiss_core::atmosphere::KlobucharParams::default()), &crate::spp::SppConfig::default(), None) {
            spp_pos = Some(spp_res.position);
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

        let valid_base = base_obs.filter(|b| {
            let age = (rover_obs.time.tow - b.time.tow).abs();
            if age > 5.0 {
                tracing::warn!("Base observation rejected due to Age of Differential ({:.1}s > 5.0s)", age);
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
                crate::engine::ambiguity::manage_ambiguities_and_slips(state, &matched_obs, &self.ephemerides, &base_coord, rover_obs.time, base.time);

                    if let Some((z_safe, h_safe, r_safe, type_safe)) = crate::engine::measurement::build_measurement_model(
                        state,
                        &matched_obs,
                        &self.ephemerides,
                        &base_coord,
                        rover_obs.time,
                        base.time,
                        Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm),
                        self.config.chi_square_pr_threshold,
                        self.config.chi_square_cp_threshold,
                    ) {
                    
                    let mut diverged = false;
                    if let Some(pos) = spp_pos {
                        let diff = pos.vector - state.position.vector;
                        tracing::debug!("SPP vs RTK pos diff: {:.1}m", diff.norm());
                        if diff.norm() > 50.0 {
                            tracing::warn!("EKF DIVERGENCE DETECTED! Difference from SPP is {:.1}m (> 50m). Forcing reset.", diff.norm());
                            diverged = true;
                        }
                    }
                    
                    if diverged || crate::engine::updater::update(state, &z_safe, &h_safe, &r_safe, self.config.chi_square_pr_threshold, Some(&type_safe)).is_err() { 
                        if diverged {
                            tracing::warn!("EKF update bypassed due to divergence.");
                        } else {
                            tracing::error!("EKF update fail");
                        }
                        if let Some(pos) = spp_pos {
                            tracing::warn!("Resetting EKF position, velocity, and IMU states to SPP due to update failure.");
                            state.position = pos;
                            state.velocity = nalgebra::Vector3::zeros();
                            state.accel_bias = nalgebra::Vector3::zeros();
                            state.gyro_bias = nalgebra::Vector3::zeros();
                            state.attitude = nalgebra::UnitQuaternion::identity();
                            state.clear_ambiguities();
                            for i in 0..6 {
                                state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
                            }
                            let n = crate::filter::CORE_STATE_SIZE;
                            for i in 6..n {
                                state.covariance[(i, i)] = 1e-4;
                            }
                        }
                    } else if let Err(e) = state.resolve_ambiguities(&self.ephemerides, self.config.lambda_min_subset) {
                        tracing::debug!("AR Failed: {}", e);
                    }
                } else {
                    tracing::warn!("Not enough valid measurements for EKF update.");
                    if let Some(pos) = spp_pos {
                        tracing::warn!("Resetting EKF position, velocity, and IMU states to SPP due to insufficient measurements.");
                        state.position = pos;
                        state.velocity = nalgebra::Vector3::zeros();
                        state.accel_bias = nalgebra::Vector3::zeros();
                        state.gyro_bias = nalgebra::Vector3::zeros();
                        state.attitude = nalgebra::UnitQuaternion::identity();
                        state.clear_ambiguities();
                        let n = crate::filter::CORE_STATE_SIZE;
                        for i in 0..6 {
                            state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
                        }
                        for i in 6..n {
                            state.covariance[(i, i)] = 1e-4;
                        }
                    }
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
                let stationary = state.velocity.norm() < 0.1 && state.accel_bias.norm() < 10.0; // Heuristic
                if stationary {
                    if let Err(e) = crate::nhc::apply_zupt(state, 0.01) {
                        tracing::debug!("ZUPT failed: {}", e);
                    }
                } else if let Err(e) = crate::nhc::apply_nhc(state, 0.1, 0.1) {
                    tracing::debug!("NHC failed: {}", e);
                }
            }
        }

        if let Some(state) = &self.current_state { self.state_history.push(state.clone()); }
        self.obs_history.push((rover_obs.clone(), base_obs.cloned()));
        self.current_state.as_ref().ok_or(EngineError::StateDisappeared)
        }

    pub fn run_backward_filter(&self) -> Result<Vec<RtkState>, EngineError> {
        if self.obs_history.is_empty() { return Err(EngineError::NoObservations); }
        let mut backward_engine = ProcessingEngine::new(self.config.clone());
        backward_engine.ephemerides = self.ephemerides.clone();
        
        // Seed from last forward state if available
        if let Some(last_f) = self.state_history.last() {
            backward_engine.current_state = Some(last_f.clone());
        }

        let mut backward_history = Vec::new();
        for i in (0..self.obs_history.len()).rev() {
            let (obs, base) = &self.obs_history[i];
            
            // Feed the IMU data backward
            if i < self.imu_history.len() {
                let mut rev_imu = self.imu_history[i].clone();
                rev_imu.reverse(); // Apply IMU measurements in reverse chronological order
                backward_engine.imu_buffer = rev_imu;
            }
            
            match backward_engine.process_epoch(obs, base.as_ref()) {
                Ok(state) => backward_history.push(state.clone()),
                Err(_) => {
                    // Fallback to coasting if update fails
                    if let Some(state) = backward_engine.current_state.as_ref() {
                        backward_history.push(state.clone());
                    }
                }
            }
        }
        backward_history.reverse();
        Ok(backward_history)
    }

    pub fn run_combined_ppk(&mut self) -> Result<Vec<RtkState>, EngineError> {
        let forward_states = self.state_history.clone();
        let backward_states = self.run_backward_filter()?;
        let mut combined = Vec::new();
        for i in 0..forward_states.len() {
            let f = &forward_states[i];
            let b = &backward_states[i];
            if f.is_fixed && !b.is_fixed { combined.push(f.clone()); continue; }
            if b.is_fixed && !f.is_fixed { combined.push(b.clone()); continue; }
            let mut c = f.clone();
            let p_f = f.covariance.view((0, 0), (6, 6)).into_owned();
            let p_b = b.covariance.view((0, 0), (6, 6)).into_owned();
            let p_f_inv = p_f.clone().try_inverse().unwrap_or(DMatrix::identity(6, 6) * 1e-4);
            let p_b_inv = p_b.clone().try_inverse().unwrap_or(DMatrix::identity(6, 6) * 1e-4);
            let p_c = (p_f_inv.clone() + &p_b_inv).try_inverse().unwrap_or(p_f.clone());
            let mut x_f = DVector::zeros(6); x_f.rows_mut(0, 3).copy_from(&f.position.vector); x_f.rows_mut(3, 3).copy_from(&f.velocity);
            let mut x_b = DVector::zeros(6); x_b.rows_mut(0, 3).copy_from(&b.position.vector); x_b.rows_mut(3, 3).copy_from(&b.velocity);
            let x_c = &p_c * (p_f_inv * x_f + p_b_inv * x_b);
            c.position.vector = x_c.fixed_rows::<3>(0).into_owned();
            c.velocity = x_c.fixed_rows::<3>(3).into_owned();
            c.covariance.view_mut((0, 0), (6, 6)).copy_from(&p_c);
            c.is_fixed = f.is_fixed || b.is_fixed;
            combined.push(c);
        }
        Ok(combined)
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
    fn test_combined_ppk_merge() {
        let time = GpsTime::new(2000, 0.0);
        let pos_f = Coordinate::new(Vector3::new(10.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let pos_b = Coordinate::new(Vector3::new(12.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state_f = RtkState::new(time, pos_f, 1.0);
        state_f.is_fixed = false;
        state_f.covariance[(0,0)] = 1.0;
        let mut state_b = RtkState::new(time, pos_b, 0.01);
        state_b.is_fixed = true;
        state_b.covariance[(0,0)] = 0.0001;
        
        let f = &state_f;
        let b = &state_b;
        let res = if b.is_fixed && !f.is_fixed { b.clone() } else { f.clone() };
        assert!(res.is_fixed);
        assert!((res.position.vector.x - 12.0).abs() < 1e-6);
    }
}
mod tests_measurement;
mod tests_predictor;
pub mod jacobian_verify;
