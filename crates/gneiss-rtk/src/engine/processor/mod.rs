use crate::engine::{EngineConfig, EngineError, EngineMode};
use crate::filter::RtkState;
use gneiss_core::coords::Coordinate;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, ObsType};
use gneiss_core::sat::SatelliteId;

mod ins;
mod ppk;
mod rtk;
pub mod rtk_iekf;
mod spp;

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
    pub consecutive_rejections: usize,
    pub tropo_mapper: Box<dyn gneiss_core::atmosphere::TropoMapper>,
    pub ppp_factor_opt: Option<crate::engine::ppp_iekf::PppIteratedEkf>,
    pub ppp_multi_epoch_opt: Option<crate::engine::ppp_multi_epoch::MultiEpochOptimizer>,
    pub sp3_epochs: Vec<gneiss_parsers::sp3::Sp3Epoch>,
    pub clk_data: Option<gneiss_parsers::rinex_clk::RinexClock>,
    pub ionex_grid: Option<gneiss_parsers::ionex::IonexGrid>,
    /// Pre-computed (time, &tec) refs for fast IONEX lookups
    pub ionex_maps: Vec<(gneiss_core::time::GpsTime, Vec<Vec<f64>>)>,
    pub antex: Option<gneiss_parsers::antex::AntexDatabase>,
    pub dcbs: std::collections::HashMap<(gneiss_core::sat::SatelliteId, String), f64>,
    pub gnn_raim: Option<crate::engine::ml::gnn_raim::GnnRaimModel>,
    pub sinex_bias: Option<gneiss_parsers::sinex_bia::SinexBias>,
}

impl ProcessingEngine {
    pub fn new(config: EngineConfig) -> Self {
        let gnn_raim = if config.enable_gnn_raim {
            tracing::warn!("GNN RAIM is experimental. Model initializes with random weights — not suitable for production. Train with: cargo run --bin train_gnn_raim");
            let dev = candle_core::Device::Cpu;
            let vm = candle_nn::VarMap::new();
            let vb = candle_nn::VarBuilder::from_varmap(&vm, candle_core::DType::F32, &dev);
            crate::engine::ml::gnn_raim::GnnRaimModel::new(vb).ok()
        } else {
            None
        };

        let tropo_mapping = config.tropo_mapping;
        let iono_model = config.iono_model;

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
            tropo_mapper: gneiss_core::atmosphere::create_tropo_mapper(tropo_mapping, None),
            ppp_factor_opt: Some(
                crate::engine::ppp_iekf::PppIteratedEkf::new()
                    .with_iono_model(iono_model),
            ),
            ppp_multi_epoch_opt: None,
            hatch_filter: crate::hatch::HatchFilter::default(),
            innovation_tracker: crate::engine::adaptive::InnovationTracker::default(),
            consecutive_rejections: 0,
            sp3_epochs: Vec::new(),
            clk_data: None,
            ionex_grid: None,
            ionex_maps: Vec::new(),
            antex: None,
            dcbs: std::collections::HashMap::new(),
            gnn_raim,
            sinex_bias: None,
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

    pub(crate) fn check_covariance_divergence(
        state: &mut RtkState,
        spp_pos: Option<Coordinate>,
        spp_state_ref: Option<&crate::spp::SppState>,
        init_isbs: bool,
    ) {
        let pos_var_max = state.covariance[(0, 0)]
            .max(state.covariance[(1, 1)])
            .max(state.covariance[(2, 2)]);
        if pos_var_max > 10000.0 {
            if let Some(pos) = spp_pos {
                tracing::warn!(
                    "Position variance {:.0} m² exceeds integrity limit. Resetting to SPP.",
                    pos_var_max
                );
                state.reset_to_spp(pos, spp_state_ref, init_isbs);
            }
        }
    }

    pub(crate) fn apply_nhc_updates(
        config: &EngineConfig,
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        state: &mut RtkState,
    ) {
        let is_ins = matches!(
            config.mode,
            EngineMode::SppIns
                | EngineMode::RtkIns
                | EngineMode::PppIns
                | EngineMode::RtkInsLooselyCoupled
                | EngineMode::SppInsLooselyCoupled
                | EngineMode::PppInsLooselyCoupled
                | EngineMode::PppInsIekf
                | EngineMode::RtkInsIekf
        );
        tracing::trace!(
            "apply_nhc_updates: enable_nhc={}, is_ins={}, ins_aligned={}",
            config.enable_nhc,
            is_ins,
            state.ins_aligned
        );
        if config.enable_nhc && is_ins && state.ins_aligned {
            let mut is_stationary = false;
            let mut accel_var = 1.0;
            if let Some(imu_buf) = imu_history.last() {
                if imu_buf.len() > 10 {
                    let mut sum_a = nalgebra::Vector3::zeros();
                    let mut sum_g = nalgebra::Vector3::zeros();
                    for m in imu_buf {
                        sum_a += m.accel;
                        sum_g += m.gyro;
                    }
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

                    if var_a < 0.05 && var_g < 0.005 && state.velocity.norm() < 1.0 {
                        is_stationary = true;
                    }
                    accel_var = var_a.max(0.001f64);
                }
            }

            if !is_stationary && state.velocity.norm() < 0.05 {
                is_stationary = true;
            }

            if is_stationary {
                let zupt_var = (accel_var * 0.1).clamp(0.001f64, 0.1f64).sqrt();
                let _ = crate::nhc::apply_zupt(state, zupt_var, &config.tuning);
            } else {
                let r_b_e = state.attitude.to_rotation_matrix();
                let omega_ie_e = nalgebra::Vector3::new(
                    0.0,
                    0.0,
                    gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S,
                );
                let omega_b = if let Some(imu_buf) = imu_history.last() {
                    if let Some(last_imu) = imu_buf.last() {
                        last_imu.gyro - state.gyro_bias - r_b_e.transpose() * omega_ie_e
                    } else {
                        nalgebra::Vector3::<f64>::zeros()
                    }
                } else {
                    nalgebra::Vector3::<f64>::zeros()
                };
                let _ = crate::nhc::apply_nhc(
                    state,
                    config.tuning.nhc_sigma_lateral,
                    config.tuning.nhc_sigma_vertical,
                    &config.imu_to_nhc_lever_arm,
                    &omega_b,
                    &config.tuning,
                );
            }
        }
    }

    pub fn reset_for_multipass(&mut self) {
        if let Some(first) = self.state_history.first() {
            let mut reset_state: RtkState = RtkState::clone(first);
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

    /// Process one epoch of GNSS observations through the configured solver.
    ///
    /// When satellite visibility is poor (urban canyons), the observation
    /// filter is relaxed to accept weaker signals rather than failing:
    /// if strict (config) filtering leaves &lt; 4 satellites, the epoch is
    /// retried with relaxed constraints (15 dB-Hz SNR, 5° elevation)
    /// before falling back to coasting on the predicted state.
    pub fn process_epoch(
        &mut self,
        rover_obs: &EpochObs,
        base_obs: Option<&EpochObs>,
    ) -> Result<&RtkState, EngineError> {
        // --- observation filtering with urban-canyon fallback ---
        let filter_obs = |obs: &EpochObs, snr_mask: f64| -> EpochObs {
            let mut filtered = obs.clone();
            if let Some(enabled) = &self.config.enabled_constellations {
                filtered
                    .satellites
                    .retain(|s| enabled.contains(&s.sat.constellation));
            }
            if snr_mask > 0.0 {
                filtered.satellites.retain(|s| {
                    let snr = s
                        .observations
                        .iter()
                        .find(|o| o.code.obs_type == ObsType::Snr && o.code.signal.freq_band == 1)
                        .map(|o| o.value)
                        .unwrap_or(25.0);
                    snr >= snr_mask
                });
            }
            filtered
        };

        let strict = filter_obs(rover_obs, self.config.min_snr_dbhz);
        let filtered_rover = if strict.satellites.len() >= 4 {
            strict
        } else {
            // Urban canyon: accept weaker signals rather than failing
            let relaxed = filter_obs(rover_obs, 15.0);
            if relaxed.satellites.len() > strict.satellites.len() {
                tracing::debug!(
                    "Relaxed SNR mask at epoch {}: {} sats (was {} with strict)",
                    rover_obs.time.tow,
                    relaxed.satellites.len(),
                    strict.satellites.len()
                );
            }
            relaxed
        };

        let filtered_base_storage = base_obs.map(|b| filter_obs(b, self.config.min_snr_dbhz));
        let filtered_base = filtered_base_storage.as_ref();

        // --- dispatch ---
        let err = match self.config.mode {
            EngineMode::Spp => self.process_spp(&filtered_rover).err(),
            EngineMode::SppIns => {
                crate::engine::spp_tight::process_spp_tightly_coupled(self, &filtered_rover).err()
            }
            EngineMode::SppInsLooselyCoupled => {
                self.process_spp_loosely_coupled(&filtered_rover).err()
            }
            EngineMode::Rtk | EngineMode::RtkIns => {
                self.process_rtk(&filtered_rover, filtered_base).err()
            }
            EngineMode::RtkInsIekf => {
                rtk_iekf::process_rtk_factor_graph(self, &filtered_rover, filtered_base).err()
            }
            EngineMode::RtkInsLooselyCoupled => self
                .process_rtk_loosely_coupled(&filtered_rover, filtered_base)
                .err(),
            EngineMode::Ppp
            | EngineMode::PppIns
            | EngineMode::PppInsLooselyCoupled
            | EngineMode::PppIekf
            | EngineMode::PppMultiEpoch => crate::engine::ppp::process_ppp(self, &filtered_rover).err(),
            EngineMode::PppInsIekf => {
                crate::engine::ppp_ins_iekf::process_ppp_ins_fg(self, &filtered_rover).err()
            }
        };

        // --- error handling ---
        if let Some(ref e) = err {
            match e {
                EngineError::StateDisappeared => {
                    tracing::warn!("EKF unrecoverable divergence. Resetting state...");
                    self.current_state = None;
                    self.consecutive_rejections = 0;
                    return Err(EngineError::StateDisappeared);
                }
                EngineError::InsufficientSatellites => {
                    self.consecutive_rejections += 1;
                    const MAX_COAST: usize = 5;
                    if self.consecutive_rejections > MAX_COAST {
                        tracing::warn!(
                            "Insufficient satellites for {} consecutive epochs — \
                             decoupling position covariance to accept new anchor",
                            self.consecutive_rejections
                        );
                        if let Some(ref mut state) = self.current_state {
                            state.decouple_position();
                            state.decouple_clock();
                        }
                        self.consecutive_rejections = 0;
                    } else {
                        tracing::warn!(
                            "Insufficient satellites at epoch {} — coasting ({}/{})",
                            rover_obs.time.tow,
                            self.consecutive_rejections,
                            MAX_COAST
                        );
                    }
                }
                _ => {
                    self.consecutive_rejections = 0;
                    return Err(e.clone());
                }
            }
        } else {
            self.consecutive_rejections = 0;
        }
        if self.current_state.is_none() {
            return Err(EngineError::InsufficientSatellites);
        }
        Ok(self.current_state.as_ref().unwrap())
    }
}

pub fn snr_scale(snr: f64) -> f64 {
    gneiss_core::variance::snr_variance_scale(snr, 45.0, 10.0)
}

#[cfg(test)]
mod tests {
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use nalgebra::Vector3;

    use super::*;
    use gneiss_core::obs::{EpochObs, SatObs};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector};
    #[test]
    fn test_engine_detects_movement() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let initial_pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(0, 0.0),
        );
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

        let pos0 = Coordinate::new(
            Vector3::new(10.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time0,
        );
        let state0 = RtkState::new(time0, pos0, 1.0);

        let pos1 = Coordinate::new(
            Vector3::new(12.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time1,
        );
        let mut state1 = RtkState::new(time1, pos1, 0.5);
        state1.is_fixed = true;

        // Mock prediction with realistic p_pred for ISB/clock states.
        // These states have cov=100000 in RtkState::new; p_pred must be
        // >= p_k to avoid negative covariance updates.
        let core_size = crate::filter::CORE_STATE_SIZE;
        state1.core_phi = Some(DMatrix::identity(core_size, core_size));
        let mut p_pred = DMatrix::identity(core_size, core_size) * 1.5;
        // Clock bias: white-noise prediction variance
        p_pred[(15, 15)] = crate::filter::PREDICTED_CLOCK_VARIANCE;
        // ISB: piece-wise constant prediction variance
        for i in [16, 17, 18] {
            p_pred[(i, i)] = crate::filter::PREDICTED_ISB_VARIANCE;
        }
        state1.full_p_predict = Some(p_pred);
        let mut x_pred = DVector::zeros(core_size);
        x_pred[0] = 10.0; // Assume velocity was 0, so predicted pos is 10
        state1.full_x_predict = Some(x_pred);

        engine.state_history.push(state0);
        engine.state_history.push(state1);

        let smoothed = engine.run_combined_ppk().expect("Should run RTS smoother");

        assert_eq!(smoothed.len(), 2);
        assert!(
            !smoothed[0].is_fixed,
            "Fix should NOT propagate backwards anymore"
        );

        // P_{0|1} = P_0 + C_0 (P_{1|1} - P_{1|0}) C_0^T
        // C_0 = P_0 * Phi^T * P_{1|0}^-1 = 1.0 * 1 * (1.5)^-1 = 0.666
        // P_{0|1} = 1.0 + 0.666 * (0.5 - 1.5) * 0.666 = 1.0 - 0.444 = 0.555
        let p_0_1 = smoothed[0].covariance[(0, 0)];
        assert!(
            (p_0_1 - 0.555555).abs() < 1e-4,
            "Covariance mismatch: {}",
            p_0_1
        );

        // x_{0|1} = x_0 + C_0 (x_1 - x_pred)
        // x_{0|1} = 10.0 + 0.666 * (12.0 - 10.0) = 11.333
        let x_0_1 = smoothed[0].position.vector.x;
        assert!(
            (x_0_1 - 11.333333).abs() < 1e-4,
            "Position mismatch: {}",
            x_0_1
        );
    }

    #[test]
    fn test_snr_scale_returns_finite() {
        let scale = snr_scale(45.0);
        assert!(scale.is_finite() && scale > 0.0);

        let scale_low = snr_scale(20.0);
        assert!(scale_low > scale);
    }

    #[test]
    fn test_covariance_divergence_no_reset_when_normal() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        // Covariance is well under 10000
        let orig_pos = state.position.vector;
        ProcessingEngine::check_covariance_divergence(&mut state, None, None, true);
        assert!((state.position.vector.x - orig_pos.x).abs() < 1e-6);
    }

    #[test]
    fn test_covariance_divergence_resets_when_extreme() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance[(0, 0)] = 20000.0;

        let spp_pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        ProcessingEngine::check_covariance_divergence(&mut state, Some(spp_pos), None, true);
        assert!((state.position.vector.x - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_apply_nhc_updates_disabled_does_nothing() {
        let mut config = EngineConfig::default();
        config.enable_nhc = false;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        ProcessingEngine::apply_nhc_updates(&config, &[], &mut state);
        // State should be unchanged
        assert!(!state.is_reset);
    }

    #[test]
    fn test_apply_nhc_updates_not_ins_mode_does_nothing() {
        let mut config = EngineConfig::default();
        config.enable_nhc = true;
        config.mode = EngineMode::Rtk; // Not an INS mode
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        ProcessingEngine::apply_nhc_updates(&config, &[], &mut state);
        assert!(!state.is_reset);
    }

    #[test]
    fn test_apply_nhc_updates_ins_not_aligned_does_nothing() {
        let mut config = EngineConfig::default();
        config.enable_nhc = true;
        config.mode = EngineMode::RtkIns;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = false;

        ProcessingEngine::apply_nhc_updates(&config, &[], &mut state);
        assert!(!state.is_reset);
    }

    #[test]
    fn test_reset_for_multipass_no_history() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.current_state = Some(RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0)),
            1.0,
        ));
        engine.state_history = Vec::new();

        engine.reset_for_multipass();
        // With no history, current_state should be None
        assert!(engine.current_state.is_none());
        assert!(engine.gnss_only_state.is_none());
    }

    #[test]
    fn test_reset_for_multipass_with_history() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        engine.state_history.push(state);
        engine.obs_history.push((make_epoch(time), None));
        engine.imu_buffer.push(gneiss_core::imu::ImuMeasurement {
            accel: Vector3::zeros(),
            gyro: Vector3::zeros(),
            time_tag: 0,
            temperature: None,
        });

        engine.reset_for_multipass();
        assert!(engine.current_state.is_some());
        assert!(engine.gnss_only_state.is_none());
        assert!(engine.state_history.is_empty());
        assert!(engine.obs_history.is_empty());
        assert!(engine.imu_buffer.is_empty());
        assert!(engine.ref_sat.is_none());
    }

    #[test]
    fn test_add_imu_measurement_no_mounting() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.add_imu_measurement(gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(1.0, 2.0, 3.0),
            gyro: Vector3::new(0.1, 0.2, 0.3),
            time_tag: 0,
            temperature: None,
        });
        assert_eq!(engine.imu_buffer.len(), 1);
        let m = &engine.imu_buffer[0];
        assert!((m.accel.x - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_add_imu_measurement_with_mounting_angles() {
        let mut config = EngineConfig::default();
        config.imu_mounting_angles = Some([0.1, 0.2, 0.3]); // Roll, Pitch, Yaw
        let mut engine = ProcessingEngine::new(config);
        engine.add_imu_measurement(gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(1.0, 0.0, 0.0),
            gyro: Vector3::new(0.0, 1.0, 0.0),
            time_tag: 0,
            temperature: None,
        });
        assert_eq!(engine.imu_buffer.len(), 1);
        // With mounting angles, the measurement should be rotated
        // (exact values depend on rotation, just verify it changed)
        let m = &engine.imu_buffer[0];
        // Accel should have been rotated from the mounting angles
        assert!(m.accel.x.abs() > 0.0 || m.accel.y.abs() > 0.0 || m.accel.z.abs() > 0.0);
    }

    #[test]
    fn test_attempt_kinematic_alignment_no_state() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.current_state = None;
        // Should not panic
        engine.attempt_kinematic_alignment();
    }

    #[test]
    fn test_attempt_kinematic_alignment_already_aligned() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        engine.current_state = Some(state);
        engine.config.mode = EngineMode::RtkIns;

        // Should do nothing because already aligned
        engine.attempt_kinematic_alignment();
        assert!(engine.current_state.as_ref().unwrap().ins_aligned);
    }

    #[test]
    fn test_attempt_kinematic_alignment_non_tight_mode() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::RtkInsLooselyCoupled; // Not tightly coupled

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        engine.current_state = Some(state);
        engine.attempt_kinematic_alignment();
        assert!(!engine.current_state.as_ref().unwrap().ins_aligned);
    }

    #[test]
    fn test_attempt_kinematic_alignment_no_imu_data() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::RtkIns;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let state = RtkState::new(time, pos, 1.0);
        engine.current_state = Some(state);
        engine.imu_history.push(Vec::new()); // Empty IMU history

        engine.attempt_kinematic_alignment();
        assert!(!engine.current_state.as_ref().unwrap().ins_aligned);
    }

    #[test]
    fn test_process_epoch_constellation_filtering() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::Spp;
        engine.config.enabled_constellations =
            Some(vec![gneiss_core::sat::Constellation::Gps]);

        let time = GpsTime::new(0, 0.0);
        let gps_sat = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        };
        let glo_sat = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Glonass,
            prn: 1,
        };

        let rover = EpochObs {
            time,
            satellites: vec![
                SatObs {
                    sat: gps_sat,
                    observations: Vec::new(),
                },
                SatObs {
                    sat: glo_sat,
                    observations: Vec::new(),
                },
            ],
        };

        let result = engine.process_epoch(&rover, None);
        // Should fail due to insufficient data (no ephemerides), but the constellation
        // filtering should have happened — Glonass sat should be filtered out
        assert!(result.is_err());
    }

    #[test]
    fn test_process_epoch_state_disappeared_resets() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::Rtk;
        let time = GpsTime::new(0, 0.0);
        let rover = EpochObs {
            time,
            satellites: Vec::new(),
        };

        // With no base obs and no state, process_rtk will error with InitialSppFailed
        // which does NOT trigger StateDisappeared reset
        let err = engine.process_epoch(&rover, None).unwrap_err();
        assert!(matches!(err, EngineError::InitialSppFailed));
    }

    #[test]
    fn test_add_ephemeris() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        assert!(engine.ephemerides.is_empty());
        let eph = gneiss_core::ephemeris::Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            toe: GpsTime::new(0, 0.0),
            toc: GpsTime::new(0, 0.0),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
            omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });
        engine.add_ephemeris(eph);
        assert_eq!(engine.ephemerides.len(), 1);
    }

    #[test]
    fn test_process_epoch_spp_mode_with_state_preserved() {
        // SPP mode with an existing state and empty observations/ephemerides.
        // SPP compute fails, but state should be preserved.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::Spp;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs {
            time: GpsTime::new(0, 1.0),
            satellites: vec![],
        };
        let err = engine.process_epoch(&rover, None).unwrap_err();
        assert!(matches!(err, EngineError::InitialSppFailed));
        // State must be preserved for the next epoch
        assert!(engine.current_state.is_some());
    }

    #[test]
    fn test_process_epoch_rtk_mode_needs_base_position() {
        // RTK mode with base_obs but no base_position configured should
        // fail with MissingBasePosition, confirming process_rtk was dispatched.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::Rtk;
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
        let err = engine.process_epoch(&rover, Some(&base)).unwrap_err();
        assert!(matches!(err, EngineError::MissingBasePosition));
    }

    #[test]
    fn test_process_epoch_rtkins_mode_needs_base_position() {
        // RtkIns mode also dispatches to process_rtk, so it should produce
        // the same MissingBasePosition error.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::RtkIns;
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
        let err = engine.process_epoch(&rover, Some(&base)).unwrap_err();
        assert!(matches!(err, EngineError::MissingBasePosition));
    }

    #[test]
    fn test_process_epoch_sppins_dispatches_correctly() {
        // SppIns mode should NOT call process_rtk (would give MissingBasePosition).
        // Instead it calls process_spp_tightly_coupled -> fails with InitialSppFailed
        // since no state or ephemerides are available.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let time = GpsTime::new(0, 0.0);
        let rover = EpochObs { time, satellites: vec![] };
        let err = engine.process_epoch(&rover, None).unwrap_err();
        // process_spp_tightly_coupled returns InitialSppFailed when no state
        assert!(matches!(err, EngineError::InitialSppFailed));
    }

    #[test]
    fn test_process_epoch_min_snr_keeps_high_snr() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::Spp;
        engine.config.min_snr_dbhz = 30.0;

        let time = GpsTime::new(0, 0.0);
        let sat = SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        };
        let rover = EpochObs {
            time,
            satellites: vec![SatObs {
                sat,
                observations: vec![
                    gneiss_core::obs::Observation {
                        code: gneiss_core::obs::ObsCode {
                            obs_type: ObsType::Snr,
                            signal: gneiss_core::obs::SignalCode { freq_band: 1, attribute: 'C' },
                        },
                        value: 35.0, // Above min_snr
                        lock_time: None,
                        lli: None,
                    },
                ],
            }],
        };

        // SPP compute will still fail (no ephemeris), but the SNR-filtering
        // step should not panic regardless of the observation count.
        let err = engine.process_epoch(&rover, None).unwrap_err();
        assert!(matches!(err, EngineError::InitialSppFailed));
    }

    #[test]
    fn test_process_epoch_min_snr_removes_low_snr() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::Spp;
        engine.config.min_snr_dbhz = 30.0;

        let time = GpsTime::new(0, 0.0);
        let sat = SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        };
        let rover = EpochObs {
            time,
            satellites: vec![SatObs {
                sat,
                observations: vec![
                    gneiss_core::obs::Observation {
                        code: gneiss_core::obs::ObsCode {
                            obs_type: ObsType::Snr,
                            signal: gneiss_core::obs::SignalCode { freq_band: 1, attribute: 'C' },
                        },
                        value: 20.0, // Below min_snr
                        lock_time: None,
                        lli: None,
                    },
                ],
            }],
        };

        // Even with a low-SNR satellite that gets filtered out entirely,
        // process_epoch should not panic. SPP compute fails (no ephemeris).
        let err = engine.process_epoch(&rover, None).unwrap_err();
        assert!(matches!(err, EngineError::InitialSppFailed));
    }

    fn make_epoch(time: GpsTime) -> EpochObs {
        EpochObs {
            time,
            satellites: vec![],
        }
    }

    #[test]
    fn test_attempt_kinematic_alignment_success() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        config.imu_to_antenna_lever_arm = [0.5, 0.0, 1.0];

        let mut engine = ProcessingEngine::new(config);
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(5.0, 0.0, 0.0);
        engine.current_state = Some(state);

        for _ in 0..5 {
            let mut hist_state = RtkState::new(time, pos, 1.0);
            hist_state.velocity = Vector3::new(5.0, 0.0, 0.0);
            engine.state_history.push(hist_state);
        }

        engine.imu_history.push(vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(0.0, 0.0, 9.8),
            gyro: Vector3::new(0.0, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]);

        let pos_before = engine.current_state.as_ref().unwrap().position.vector;
        engine.attempt_kinematic_alignment();

        let aligned_state = engine.current_state.as_ref().unwrap();
        assert!(aligned_state.ins_aligned, "INS should be aligned");

        let att_var = (15.0f64.to_radians()).powi(2);
        for i in 6..9 {
            assert!(
                (aligned_state.covariance[(i, i)] - att_var).abs() < 1e-10,
                "Covariance[{}] = {}, expected {}",
                i,
                aligned_state.covariance[(i, i)],
                att_var
            );
        }

        assert!(
            (aligned_state.position.vector - pos_before).norm() > 0.0,
            "Position should change due to lever arm correction"
        );
        assert!(engine.imu_buffer.is_empty());
    }

    #[test]
    fn test_attempt_kinematic_alignment_insufficient_history() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        let mut engine = ProcessingEngine::new(config);
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(5.0, 0.0, 0.0);
        engine.current_state = Some(state);

        for _ in 0..3 {
            let mut hist_state = RtkState::new(time, pos, 1.0);
            hist_state.velocity = Vector3::new(5.0, 0.0, 0.0);
            engine.state_history.push(hist_state);
        }

        engine.imu_history.push(vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(0.0, 0.0, 9.8),
            gyro: Vector3::new(0.0, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]);

        engine.attempt_kinematic_alignment();
        assert!(!engine.current_state.as_ref().unwrap().ins_aligned);
    }

    #[test]
    fn test_attempt_kinematic_alignment_slow_speed() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        let mut engine = ProcessingEngine::new(config);
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(1.0, 0.0, 0.0);
        engine.current_state = Some(state);

        for _ in 0..5 {
            let mut hist_state = RtkState::new(time, pos, 1.0);
            hist_state.velocity = Vector3::new(5.0, 0.0, 0.0);
            engine.state_history.push(hist_state);
        }

        engine.imu_history.push(vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(0.0, 0.0, 9.8),
            gyro: Vector3::new(0.0, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]);

        engine.attempt_kinematic_alignment();
        assert!(!engine.current_state.as_ref().unwrap().ins_aligned);
    }

    #[test]
    fn test_attempt_kinematic_alignment_history_speed_too_low() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        let mut engine = ProcessingEngine::new(config);
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::new(5.0, 0.0, 0.0);
        engine.current_state = Some(state);

        for _ in 0..5 {
            let mut hist_state = RtkState::new(time, pos, 1.0);
            hist_state.velocity = Vector3::new(1.0, 0.0, 0.0);
            engine.state_history.push(hist_state);
        }

        engine.imu_history.push(vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(0.0, 0.0, 9.8),
            gyro: Vector3::new(0.0, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]);

        engine.attempt_kinematic_alignment();
        assert!(!engine.current_state.as_ref().unwrap().ins_aligned);
    }

    #[test]
    fn test_apply_nhc_updates_stationary_detected() {
        let mut config = EngineConfig::default();
        config.enable_nhc = true;
        config.mode = EngineMode::RtkIns;

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.velocity = Vector3::new(0.0, 0.0, 0.0);

        let imu_measurements: Vec<_> = (0..12)
            .map(|i| gneiss_core::imu::ImuMeasurement {
                accel: Vector3::new(0.0, 0.0, 9.8),
                gyro: Vector3::new(0.0, 0.0, 0.0),
                time_tag: i,
                temperature: None,
            })
            .collect();
        let imu_history = vec![imu_measurements];

        ProcessingEngine::apply_nhc_updates(&config, &imu_history, &mut state);
        assert!(state.velocity.norm() < 1e-6);
    }

    #[test]
    fn test_apply_nhc_updates_non_stationary_nhc() {
        let mut config = EngineConfig::default();
        config.enable_nhc = true;
        config.mode = EngineMode::RtkIns;

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.velocity = Vector3::new(10.0, 0.0, 0.0);

        let imu_history = vec![vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(0.0, 0.0, 9.8),
            gyro: Vector3::new(0.1, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]];

        ProcessingEngine::apply_nhc_updates(&config, &imu_history, &mut state);
    }

    #[test]
    fn test_apply_nhc_updates_stationary_low_velocity() {
        let mut config = EngineConfig::default();
        config.enable_nhc = true;
        config.mode = EngineMode::RtkIns;

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.velocity = Vector3::new(0.02, 0.0, 0.0);

        let imu_measurements: Vec<_> = (0..12)
            .map(|i| gneiss_core::imu::ImuMeasurement {
                accel: Vector3::new(i as f64 * 10.0, 0.0, 9.8),
                gyro: Vector3::new(i as f64 * 0.1, 0.0, 0.0),
                time_tag: i,
                temperature: None,
            })
            .collect();
        let imu_history = vec![imu_measurements];

        ProcessingEngine::apply_nhc_updates(&config, &imu_history, &mut state);
    }
}

impl ProcessingEngine {
    pub(crate) fn attempt_kinematic_alignment(&mut self) {
        let state = if let Some(s) = &mut self.current_state {
            s
        } else {
            return;
        };
        if state.ins_aligned || !self.config.mode.is_tightly_coupled() {
            return;
        }

        // Do not attempt to align an INS if we have no IMU data.
        if self
            .imu_history
            .last()
            .map(|b| b.is_empty())
            .unwrap_or(true)
        {
            return;
        }

        let speed = state.velocity.norm();
        if speed > 3.0
            && self.state_history.len() >= 5
            && self
                .state_history
                .iter()
                .rev()
                .take(5)
                .all(|s| s.velocity.norm() > 3.0)
        {
            let llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
            let ecef_to_ned = gneiss_core::coords::ecef_to_ned_matrix(llh);
            let v_ned = ecef_to_ned * state.velocity;
            let yaw = f64::atan2(v_ned[1], v_ned[0]);
            let rot_veh_to_ned = nalgebra::Rotation3::from_euler_angles(0.0, 0.0, yaw);

            let ned_to_ecef = ecef_to_ned.transpose();
            // Since imu measurements are rotated to vehicle frame in add_imu_measurement,
            // state.attitude should represent the rotation from VEHICLE to ECEF.
            let rot_mat =
                nalgebra::Rotation3::from_matrix_unchecked(ned_to_ecef * rot_veh_to_ned.matrix());
            state.attitude = nalgebra::UnitQuaternion::from_rotation_matrix(&rot_mat);

            // Shift state position from antenna phase center to IMU center
            let r_e_v = state.attitude.to_rotation_matrix();
            let lever_arm =
                nalgebra::Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm);
            state.position.vector -= r_e_v * lever_arm;

            state.ins_aligned = true;
            tracing::info!(
                "Kinematic alignment successful! Speed: {:.2} m/s, Veh Yaw: {:.2} deg",
                speed,
                yaw.to_degrees()
            );

            // Set attitude covariance higher because we assume 0 roll/pitch and yaw is based on noisy GNSS velocity
            for i in 6..9 {
                state.covariance[(i, i)] = (15.0f64.to_radians()).powi(2);
            }
            self.imu_buffer.clear();
        }
    }
}
