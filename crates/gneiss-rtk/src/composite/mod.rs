//! Unified Composite Navigation Integration.
//!
//! Provides modular interfaces composing the 15-state Error-State Kalman Filter
//! (ESKF) with:
//! 1. **Integer PPP-AR (`tc_ppp`)**: Tightly-Coupled PPP/INS pipeline,
//!    coupling un-differenced carrier-phase and pseudorange observations with
//!    attitude/lever-arm Jacobians and integer ambiguity resolution.
//! 2. **Network RTK VRS (`tc_rtk`)**: Tightly-Coupled Network RTK/INS
//!    pipeline, coupling localized VRS synthesized observables and
//!    double-difference carrier-phase innovations with the 15-state ESKF.
//! 3. **Unified Composite Multi-Mode Engine (`UnifiedCompositeEngine`)**:
//!    Seamless mode-switching between TC-RTK and TC-PPP with hysteresis and
//!    inertial bias/covariance preservation.

pub mod tc_ppp;
pub mod tc_rtk;

use nalgebra::{UnitQuaternion, Vector3};

use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;

pub use crate::estimators::eskf::types::EngineError;
pub use crate::estimators::eskf::{EskfState, Matrix15, Vector15};

pub use tc_ppp::{FloatAmbiguity, TcPppConfig, TightlyCoupledPppIns};
pub use tc_rtk::{DoubleDiffAmbiguity, TcRtkConfig, TightlyCoupledNetworkRtkIns};

/// Type alias for epoch observation data.
pub type EpochObservation = EpochObs;

/// Operational composite navigation mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompositeMode {
    /// Tightly-coupled PPP/INS (base-station free).
    TightlyCoupledPpp,
    /// Tightly-coupled Network RTK/INS (VRS assisted).
    TightlyCoupledRtk,
    /// Inertial dead-reckoning during GNSS outage.
    DeadReckoning,
}

/// Consolidated navigation solution produced by composite engines.
#[derive(Debug, Clone, PartialEq)]
pub struct NavSolution {
    /// Solution epoch time.
    pub time: GpsTime,
    /// ECEF position (meters).
    pub pos_ecef: Vector3<f64>,
    /// ECEF velocity (meters/second).
    pub vel_ecef: Vector3<f64>,
    /// Body-to-ECEF attitude quaternion.
    pub attitude: UnitQuaternion<f64>,
    /// Estimated accelerometer bias (m/s²).
    pub accel_bias: Vector3<f64>,
    /// Estimated gyroscope bias (rad/s).
    pub gyro_bias: Vector3<f64>,
    /// 15x15 error state covariance matrix.
    pub cov: Matrix15<f64>,
    /// Number of satellites utilized in epoch solution.
    pub num_satellites: usize,
    /// Whether integer ambiguities were resolved.
    pub is_fixed: bool,
    /// Active composite navigation mode.
    pub mode: CompositeMode,
    /// Ambiguity resolution validation ratio test statistic.
    pub ratio: Option<f64>,
}

/// CORS network station observation packet for a single epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct StationEpoch {
    /// Station alphanumeric identifier (e.g. "P181").
    pub station_id: String,
    /// Approximate or surveyed ECEF coordinates (meters).
    pub station_pos: Vector3<f64>,
    /// Raw observation data for this epoch.
    pub obs: EpochObs,
}

/// Unified multi-mode coordinator managing seamless transitions between
/// Tightly-Coupled Network RTK/INS and Tightly-Coupled PPP/INS.
pub struct UnifiedCompositeEngine {
    /// Active navigation mode.
    pub mode: CompositeMode,
    /// Tightly-coupled PPP/INS pipeline.
    pub tc_ppp: TightlyCoupledPppIns,
    /// Tightly-coupled Network RTK/INS pipeline.
    pub tc_rtk: TightlyCoupledNetworkRtkIns,
    /// Consecutive epoch counter for mode hold hysteresis.
    pub mode_hold_counter: usize,
    /// Minimum required epochs before switching mode back (prevents chatter).
    pub min_hold_epochs: usize,
    /// Maximum network baseline radius (meters) for RTK operation.
    pub cors_max_radius_m: f64,
}

impl UnifiedCompositeEngine {
    /// Construct a new unified multi-mode composite navigation engine.
    pub fn new(
        tc_ppp: TightlyCoupledPppIns,
        tc_rtk: TightlyCoupledNetworkRtkIns,
        cors_max_radius_m: f64,
    ) -> Self {
        Self {
            mode: CompositeMode::TightlyCoupledRtk,
            tc_ppp,
            tc_rtk,
            mode_hold_counter: 0,
            min_hold_epochs: 10,
            cors_max_radius_m,
        }
    }

    /// Retrieve active ESKF state reference.
    pub fn active_state(&self) -> &EskfState {
        match self.mode {
            CompositeMode::TightlyCoupledRtk => &self.tc_rtk.eskf,
            CompositeMode::TightlyCoupledPpp | CompositeMode::DeadReckoning => &self.tc_ppp.eskf,
        }
    }

    /// Retrieve mutable active ESKF state reference.
    pub fn active_state_mut(&mut self) -> &mut EskfState {
        match self.mode {
            CompositeMode::TightlyCoupledRtk => &mut self.tc_rtk.eskf,
            CompositeMode::TightlyCoupledPpp | CompositeMode::DeadReckoning => {
                &mut self.tc_ppp.eskf
            }
        }
    }

    /// Seamlessly transfer full 15-state ESKF estimates and covariance
    /// between pipelines while preserving sensor biases and attitude.
    pub fn switch_mode(&mut self, target_mode: CompositeMode) {
        if self.mode == target_mode {
            return;
        }
        match (self.mode, target_mode) {
            (CompositeMode::TightlyCoupledRtk, CompositeMode::TightlyCoupledPpp) => {
                Self::copy_state(&self.tc_rtk.eskf, &mut self.tc_ppp.eskf);
                self.mode = CompositeMode::TightlyCoupledPpp;
                self.mode_hold_counter = 0;
            }
            (CompositeMode::TightlyCoupledPpp, CompositeMode::TightlyCoupledRtk) => {
                Self::copy_state(&self.tc_ppp.eskf, &mut self.tc_rtk.eskf);
                self.mode = CompositeMode::TightlyCoupledRtk;
                self.mode_hold_counter = 0;
            }
            (_, target) => {
                self.mode = target;
                self.mode_hold_counter = 0;
            }
        }
    }

    /// Process incoming IMU and GNSS observations with autonomous mode selection.
    pub fn process_epoch(
        &mut self,
        imu_samples: &[gneiss_core::imu::ImuMeasurement],
        rover_obs: &EpochObservation,
        cors_obs: Option<&[StationEpoch]>,
    ) -> Result<NavSolution, EngineError> {
        self.mode_hold_counter = self.mode_hold_counter.saturating_add(1);
        let can_use_rtk = cors_obs.is_some_and(|stations| !stations.is_empty());

        self.evaluate_mode_transition(can_use_rtk);
        self.execute_active_epoch(imu_samples, rover_obs, cors_obs)
    }

    fn evaluate_mode_transition(&mut self, can_use_rtk: bool) {
        if self.mode_hold_counter < self.min_hold_epochs {
            return;
        }
        if self.mode == CompositeMode::TightlyCoupledRtk && !can_use_rtk {
            self.switch_mode(CompositeMode::TightlyCoupledPpp);
        } else if self.mode == CompositeMode::TightlyCoupledPpp && can_use_rtk {
            self.switch_mode(CompositeMode::TightlyCoupledRtk);
        }
    }

    fn execute_active_epoch(
        &mut self,
        imu_samples: &[gneiss_core::imu::ImuMeasurement],
        rover_obs: &EpochObservation,
        cors_obs: Option<&[StationEpoch]>,
    ) -> Result<NavSolution, EngineError> {
        match self.mode {
            CompositeMode::TightlyCoupledRtk => {
                let stations = cors_obs.unwrap_or(&[]);
                self.tc_rtk.process_epoch(imu_samples, rover_obs, stations)
            }
            CompositeMode::TightlyCoupledPpp | CompositeMode::DeadReckoning => {
                self.tc_ppp.process_epoch(imu_samples, rover_obs)
            }
        }
    }

    fn copy_state(src: &EskfState, dst: &mut EskfState) {
        dst.pos_ecef = src.pos_ecef;
        dst.vel_ecef = src.vel_ecef;
        dst.attitude = src.attitude;
        dst.accel_bias = src.accel_bias;
        dst.gyro_bias = src.gyro_bias;
        dst.cov = 0.5 * (src.cov + src.cov.transpose());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::post_process::vrs::VrsSynthesizer;

    #[test]
    fn test_mode_switch_preserves_state_and_biases() {
        let pos = Vector3::new(100.0, 200.0, 300.0);
        let vel = Vector3::new(1.0, 2.0, 3.0);
        let att = UnitQuaternion::identity();
        let mut eskf_init = EskfState::new(pos, vel, att);
        eskf_init.accel_bias = Vector3::new(0.01, -0.02, 0.03);
        eskf_init.gyro_bias = Vector3::new(0.001, -0.002, 0.001);

        let ppp = TightlyCoupledPppIns::new(eskf_init.clone(), TcPppConfig::default());
        let synth = VrsSynthesizer::new("TEST", pos);
        let rtk = TightlyCoupledNetworkRtkIns::new(eskf_init, synth, TcRtkConfig::default());

        let mut engine = UnifiedCompositeEngine::new(ppp, rtk, 50_000.0);
        assert_eq!(engine.mode, CompositeMode::TightlyCoupledRtk);

        engine.switch_mode(CompositeMode::TightlyCoupledPpp);
        assert_eq!(engine.mode, CompositeMode::TightlyCoupledPpp);

        let ppp_state = &engine.tc_ppp.eskf;
        assert_eq!(ppp_state.accel_bias, Vector3::new(0.01, -0.02, 0.03));
        assert_eq!(ppp_state.gyro_bias, Vector3::new(0.001, -0.002, 0.001));
        assert_eq!(ppp_state.pos_ecef, pos);
    }

    #[test]
    fn test_hysteresis_prevents_rapid_chatter() {
        let eskf = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
        let ppp = TightlyCoupledPppIns::new(eskf.clone(), TcPppConfig::default());
        let synth = VrsSynthesizer::new("TEST", Vector3::zeros());
        let rtk = TightlyCoupledNetworkRtkIns::new(eskf, synth, TcRtkConfig::default());
        let mut engine = UnifiedCompositeEngine::new(ppp, rtk, 50_000.0);

        engine.evaluate_mode_transition(false);
        assert_eq!(engine.mode, CompositeMode::TightlyCoupledRtk);

        engine.mode_hold_counter = 10;
        engine.evaluate_mode_transition(false);
        assert_eq!(engine.mode, CompositeMode::TightlyCoupledPpp);
    }

    #[test]
    fn test_unified_composite_engine_active_state_access() {
        let pos = Vector3::new(50.0, 60.0, 70.0);
        let eskf = EskfState::new(pos, Vector3::zeros(), UnitQuaternion::identity());
        let ppp = TightlyCoupledPppIns::new(eskf.clone(), TcPppConfig::default());
        let synth = VrsSynthesizer::new("TEST", pos);
        let rtk = TightlyCoupledNetworkRtkIns::new(eskf, synth, TcRtkConfig::default());
        let mut engine = UnifiedCompositeEngine::new(ppp, rtk, 50_000.0);

        assert_eq!(engine.active_state().pos_ecef, pos);
        engine.active_state_mut().pos_ecef += Vector3::new(1.0, 1.0, 1.0);
        assert_eq!(engine.active_state().pos_ecef, pos + Vector3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn test_unified_composite_engine_process_epoch() {
        let pos = Vector3::new(10.0, 20.0, 30.0);
        let eskf = EskfState::new(pos, Vector3::zeros(), UnitQuaternion::identity());
        let ppp = TightlyCoupledPppIns::new(eskf.clone(), TcPppConfig::default());
        let synth = VrsSynthesizer::new("TEST", pos);
        let rtk = TightlyCoupledNetworkRtkIns::new(eskf, synth, TcRtkConfig::default());
        let mut engine = UnifiedCompositeEngine::new(ppp, rtk, 50_000.0);

        let imu = vec![gneiss_core::imu::ImuMeasurement::new(1000, Vector3::new(0.0, 0.0, -9.81), Vector3::zeros())];
        let obs = EpochObs { time: GpsTime::new(2200, 10.0), satellites: vec![] };
        let sol = engine.process_epoch(&imu, &obs, None).expect("process_epoch failed");
        assert_eq!(sol.mode, CompositeMode::DeadReckoning);
    }
}
