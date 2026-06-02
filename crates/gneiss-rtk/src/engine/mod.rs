pub mod matcher;
pub mod predictor;
pub mod updater;
pub mod measurement;
pub mod ppp;

use nalgebra::{Vector3, DMatrix, DVector};
use gneiss_core::obs::EpochObs;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::SatelliteId;
use crate::filter::RtkState;
use crate::engine::matcher::match_observations;

/// Defines the core positioning mode
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub enum EngineMode {
    Spp,
    SppIns,
    Rtk,
    RtkIns,
    Ppp,
    PppIns,
}

impl Default for EngineMode {
    fn default() -> Self {
        EngineMode::Rtk
    }
}

/// Configuration for the RTK processing engine.
#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
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
    
    // Tuning Parameters
    pub raim_pseudorange_outlier_m: f64,
    pub chi_square_pr_threshold: f64,
    pub chi_square_cp_threshold: f64,
    pub nominal_snr_dbhz: f64,
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
            enable_nhc: true,
            enable_backward_smoothing: false,
            lambda_min_ratio: 3.0,
            lambda_min_subset: 5,
            raim_pseudorange_outlier_m: 25.0,
            chi_square_pr_threshold: 15.0,
            chi_square_cp_threshold: 1000000.0,
            nominal_snr_dbhz: 40.0,
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
            let enable_imu = matches!(self.config.mode, EngineMode::SppIns | EngineMode::RtkIns | EngineMode::PppIns);
            crate::engine::predictor::predict(state, dt, if enable_imu { 0.1 } else { 10.0 }, &self.imu_buffer);
        }
        self.imu_history.push(self.imu_buffer.clone());
        self.imu_buffer.clear();
    }

    pub fn process_epoch(&mut self, rover_obs: &EpochObs, base_obs: Option<&EpochObs>) -> Result<&RtkState, EngineError> {
        let err = match self.config.mode {
            EngineMode::Spp | EngineMode::SppIns => self.process_spp(rover_obs).err(),
            EngineMode::Rtk | EngineMode::RtkIns => self.process_rtk(rover_obs, base_obs).err(),
            EngineMode::Ppp | EngineMode::PppIns => crate::engine::ppp::process_ppp(self, rover_obs).err(),
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

    pub fn process_spp(&mut self, rover_obs: &EpochObs) -> Result<&RtkState, EngineError> {
        let mut spp_pos = None;
        let mut spp_cdt = 0.0;
        if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, None, &crate::spp::SppConfig::default(), None) {
            spp_pos = Some(spp_res.position);
            spp_cdt = spp_res.cdt;
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

            if let Err(e) = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat) {
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
        if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, None, &crate::spp::SppConfig::default(), None) {
            spp_pos = Some(spp_res.position);
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
            
            let matched_obs = match_observations(rover_obs, base);
            
            if state.epoch_count < 5 || state.epoch_count % 100 == 0 {
                tracing::info!("Epoch {}: Matched {} satellites", state.epoch_count, matched_obs.len());
            }

            if matched_obs.len() >= 5 {
                for (r, b) in &matched_obs {
                    let mut slip = false;
                    let prev_lock = *state.locktimes.get(&(r.sat, 1)).unwrap_or(&0);
                    let mut new_lock = prev_lock + 1;
                    
                    if let Some(r_lock) = r.locktime {
                        if r_lock == 0 {
                            slip = true;
                            new_lock = 0;
                        } else if r_lock < prev_lock {
                            slip = true;
                            new_lock = r_lock;
                        } else {
                            new_lock = r_lock;
                        }
                    }

                    if new_lock == 0 && state.locktimes.contains_key(&(r.sat, 1)) {
                        slip = true;
                    }

                    state.locktimes.insert((r.sat, 1), new_lock);
                    state.locktimes.insert((r.sat, 2), new_lock);
                    
                    let r_freq_num = self.ephemerides.iter().find(|e| e.sat() == r.sat).map(|e| e.freq_num()).unwrap_or(0);
                    let b_freq_num = self.ephemerides.iter().find(|e| e.sat() == b.sat).map(|e| e.freq_num()).unwrap_or(0);
                    
                    let (r_f1, r_f2) = gneiss_core::signal::satellite_frequencies(r.sat, r_freq_num);
                    let (b_f1, b_f2) = gneiss_core::signal::satellite_frequencies(b.sat, b_freq_num);

                    let mut slip_l1 = slip;
                    let mut slip_l2 = slip;
                    if *state.reject_counts.get(&(r.sat, 1)).unwrap_or(&0) > 3 { slip_l1 = true; }
                    if *state.reject_counts.get(&(r.sat, 2)).unwrap_or(&0) > 3 { slip_l2 = true; }

                    if !state.ambiguity_keys.contains(&(r.sat, 1)) || slip_l1 || slip_l2 {
                        if slip_l1 || slip_l2 { 
                            state.remove_ambiguity(r.sat, 1); 
                            state.remove_ambiguity(r.sat, 2);
                            state.reject_counts.insert((r.sat, 1), 0);
                            state.reject_counts.insert((r.sat, 2), 0);
                        }
                        
                        if let (Some(r_cp1), Some(b_cp1)) = (r.cp_l1, b.cp_l1) {
                            let lam_r1 = 299792458.0 / r_f1;
                            let lam_b1 = 299792458.0 / b_f1;
                            let cp_l1_rov = r_cp1 * lam_r1;
                            let cp_l1_base = b_cp1 * lam_b1;

                            let mut initialized = false;
                            if state.covariance[(0,0)] < 0.1 {
                                // Try to find an anchor satellite to cancel the clock bias
                                for (anchor_r, anchor_b) in matched_obs.iter() {
                                    if anchor_r.sat == r.sat { continue; }
                                    if let Some(anchor_idx) = state.ambiguity_keys.iter().position(|&(s, f)| s == anchor_r.sat && f == 1) {
                                        if state.covariance[(15 + anchor_idx, 15 + anchor_idx)] < 0.05 {
                                            if let (Some(ar_cp), Some(ab_cp)) = (anchor_r.cp_l1, anchor_b.cp_l1) {
                                                let (ar_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == anchor_r.sat).unwrap().position(rover_obs.time);
                                                let (ab_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == anchor_b.sat).unwrap().position(base.time);
                                                let ar_dist_rov = (state.position.vector - ar_sat_vec).norm();
                                                let ar_dist_base = (base_coord.vector - ab_sat_vec).norm();
                                                
                                                let a_lam_r = 299792458.0 / r_f1;
                                                let a_lam_b = 299792458.0 / b_f1;
                                                let a_cp_rov = ar_cp * a_lam_r;
                                                let a_cp_base = ab_cp * a_lam_b;
                                                
                                                let anchor_sd = state.ambiguities[anchor_idx];
                                                let b_clock_rov = a_cp_rov - ar_dist_rov - anchor_sd;
                                                let b_clock_base = a_cp_base - ar_dist_base; // Base has no ambiguity in SD, assuming SD = rov - base

                                                let (r_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == r.sat).unwrap().position(rover_obs.time);
                                                let (b_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == b.sat).unwrap().position(base.time);
                                                let dist_rov = (state.position.vector - r_sat_vec).norm();
                                                let dist_base = (base_coord.vector - b_sat_vec).norm();
                                                
                                                let initial_est_l1 = (cp_l1_rov - dist_rov - b_clock_rov) - (cp_l1_base - dist_base - b_clock_base);
                                                state.add_ambiguity(r.sat, 1, initial_est_l1, 100.0);
                                                initialized = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                            
                            if !initialized {
                                let initial_est_l1 = (cp_l1_rov - r.pr_l1) - (cp_l1_base - b.pr_l1);
                                state.add_ambiguity(r.sat, 1, initial_est_l1, 10000.0);
                            }
                        }
                        
                        if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2)) = (r.pr_l2, r.cp_l2, b.pr_l2, b.cp_l2) {
                            let lam_r2 = 299792458.0 / r_f2;
                            let lam_b2 = 299792458.0 / b_f2;
                            let cp_l2_rov = r_cp2 * lam_r2;
                            let cp_l2_base = b_cp2 * lam_b2;
                            
                            let mut initialized = false;
                            if state.covariance[(0,0)] < 0.1 {
                                for (anchor_r, anchor_b) in matched_obs.iter() {
                                    if anchor_r.sat == r.sat { continue; }
                                    if let Some(anchor_idx) = state.ambiguity_keys.iter().position(|&(s, f)| s == anchor_r.sat && f == 2) {
                                        if state.covariance[(15 + anchor_idx, 15 + anchor_idx)] < 0.05 {
                                            if let (Some(ar_cp), Some(ab_cp)) = (anchor_r.cp_l2, anchor_b.cp_l2) {
                                                let (ar_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == anchor_r.sat).unwrap().position(rover_obs.time);
                                                let (ab_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == anchor_b.sat).unwrap().position(base.time);
                                                let ar_dist_rov = (state.position.vector - ar_sat_vec).norm();
                                                let ar_dist_base = (base_coord.vector - ab_sat_vec).norm();
                                                
                                                let a_lam_r = 299792458.0 / r_f2;
                                                let a_lam_b = 299792458.0 / b_f2;
                                                let a_cp_rov = ar_cp * a_lam_r;
                                                let a_cp_base = ab_cp * a_lam_b;
                                                
                                                let anchor_sd = state.ambiguities[anchor_idx];
                                                let b_clock_rov = a_cp_rov - ar_dist_rov - anchor_sd;
                                                let b_clock_base = a_cp_base - ar_dist_base;

                                                let (r_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == r.sat).unwrap().position(rover_obs.time);
                                                let (b_sat_vec, _, _, _) = self.ephemerides.iter().find(|e| e.sat() == b.sat).unwrap().position(base.time);
                                                let dist_rov = (state.position.vector - r_sat_vec).norm();
                                                let dist_base = (base_coord.vector - b_sat_vec).norm();
                                                
                                                let initial_est_l2 = (cp_l2_rov - dist_rov - b_clock_rov) - (cp_l2_base - dist_base - b_clock_base);
                                                state.add_ambiguity(r.sat, 2, initial_est_l2, 100.0);
                                                initialized = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                            
                            if !initialized {
                                let initial_est_l2 = (cp_l2_rov - r_pr2) - (cp_l2_base - b_pr2);
                                state.add_ambiguity(r.sat, 2, initial_est_l2, 10000.0);
                            }
                        }
                    }
                    if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2), Some(r_cp1), Some(b_cp1)) = (r.pr_l2, r.cp_l2, b.pr_l2, b.cp_l2, r.cp_l1, b.cp_l1) {
                        let r_cp1_m = r_cp1 * (299792458.0 / r_f1);
                        let r_cp2_m = r_cp2 * (299792458.0 / r_f2);
                        let b_cp1_m = b_cp1 * (299792458.0 / b_f1);
                        let b_cp2_m = b_cp2 * (299792458.0 / b_f2);
                        let mw_sd = crate::combinations::melbourne_wubbena(r_cp1_m, r_cp2_m, r.pr_l1, r_pr2, r_f1, r_f2) - 
                                    crate::combinations::melbourne_wubbena(b_cp1_m, b_cp2_m, b.pr_l1, b_pr2, b_f1, b_f2);
                        state.update_mw(r.sat, mw_sd / crate::combinations::lambda_wl(r_f1, r_f2));
                    }
                }

                if let Some((z_safe, h_safe, r_safe)) = crate::engine::measurement::build_measurement_model(
                    state,
                    &matched_obs,
                    &self.ephemerides,
                    &base_coord,
                    rover_obs.time,
                    base.time,
                    Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm),
                ) {
                    if let Err(e) = crate::engine::updater::update(state, &z_safe, &h_safe, &r_safe) { 
                        tracing::error!("EKF update fail: {:?}", e); 
                    } else if let Err(e) = state.resolve_ambiguities(&self.ephemerides, self.config.lambda_min_subset) {
                        tracing::debug!("AR Failed: {}", e);
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

            if let Err(e) = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat) {
                tracing::debug!("SPP Fallback update failed: {:?}", e);
            }
        }
        
        if self.config.enable_nhc {
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

pub fn snr_scale(snr: f64) -> f64 { libm::pow(10.0, (45.0 - snr.clamp(25.0, 50.0)) / 10.0).min(100.0) }

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
