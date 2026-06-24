use crate::engine::types::{DynamicsModel, EngineMode, IonosphereModel};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(default)]
pub struct AutoTuneConfig {
    pub enabled: bool,
    pub max_sigma_ab: f64,
    pub min_sigma_ab: f64,
    pub max_sigma_gb: f64,
    pub min_sigma_gb: f64,
}

impl Default for AutoTuneConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            max_sigma_ab: 1e-2,
            min_sigma_ab: 1e-6,
            max_sigma_gb: 1e-3,
            min_sigma_gb: 1e-7,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(default)]
pub struct EkfTuningConfig {
    // Measurement Base Variances
    pub pr_base_var: f64,
    pub cp_base_var: f64,
    pub dop_base_var: f64,

    // SNR Variance Model: a^2 + b^2 / 10^(CN0/10)
    pub snr_a: f64,
    pub snr_b: f64,

    // IMU Process Noise
    pub sigma_v: f64,   // Velocity Random Walk
    pub sigma_phi: f64, // Angular Random Walk
    pub sigma_ab: f64,  // Accel bias instability
    pub sigma_gb: f64,  // Gyro bias instability

    // Outlier Thresholds
    pub loosely_coupled_mahalanobis_sq: f64,
    pub phase_outlier_ratio_thresh: f64,
    pub doppler_outlier_ratio_mult: f64,
    pub pr_abs_thresh: f64,
    pub cp_abs_thresh: f64,
    pub dop_abs_thresh: f64,

    // Huber Estimator Thresholds
    pub huber_threshold_loosely: f64,
    pub huber_threshold_tightly: f64,

    // Ambiguity Resolution
    pub min_ar_success_rate: f64,

    pub ekf_max_iterations: usize,

    // FGO robust thresholds
    pub fgo_pr_robust_threshold: f64,
    pub fgo_cp_robust_threshold: f64,
    pub fgo_dop_robust_threshold: f64,

    // Auto-tuning Constraints
    pub auto_tune: AutoTuneConfig,

    pub nhc_sigma_lateral: f64,
    pub nhc_sigma_vertical: f64,
}

impl Default for EkfTuningConfig {
    fn default() -> Self {
        Self {
            pr_base_var: 1.0,
            cp_base_var: 9e-6,
            dop_base_var: 1.0,
            snr_a: 1.0,
            snr_b: 150.0,
            sigma_v: 0.1,
            sigma_phi: 0.01,
            sigma_ab: 1e-4,
            sigma_gb: 1e-5,
            loosely_coupled_mahalanobis_sq: 1000.0,
            phase_outlier_ratio_thresh: 5.0,
            doppler_outlier_ratio_mult: 2.0,
            pr_abs_thresh: 50.0,
            cp_abs_thresh: 0.10,
            dop_abs_thresh: 2.0,
            huber_threshold_loosely: 10.0,
            huber_threshold_tightly: 3.0,
            min_ar_success_rate: 0.999,
            ekf_max_iterations: 3,
            fgo_pr_robust_threshold: 10.0,
            fgo_cp_robust_threshold: 3.0,
            fgo_dop_robust_threshold: 5.0,
            auto_tune: Default::default(),
            nhc_sigma_lateral: 2.0,
            nhc_sigma_vertical: 2.0,
        }
    }
}

/// Configuration for the RTK processing engine.
#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
#[serde(default)]
pub struct EngineConfig {
    pub mode: EngineMode,
    pub initial_position: Option<[f64; 3]>,
    pub base_position: Option<[f64; 3]>,
    pub receiver_antenna_type: Option<String>,
    pub base_datum_transform: Option<gneiss_geodesy::helmert::HelmertParams>,
    pub imu_to_antenna_lever_arm: [f64; 3],
    pub imu_mounting_angles: Option<[f64; 3]>, // [Roll, Pitch, Yaw] in radians
    pub imu_to_nhc_lever_arm: [f64; 3],        // [x, y, z] from IMU to NHC point in body frame
    pub enable_nhc: bool,
    pub enable_backward_smoothing: bool,
    pub lambda_min_ratio: f64,
    pub lambda_min_subset: usize,
    pub enabled_constellations: Option<Vec<gneiss_core::sat::Constellation>>,

    // Tuning Parameters
    pub raim_pseudorange_outlier_m: f64,
    pub chi_square_pr_threshold: f64,
    pub chi_square_cp_threshold: f64,
    pub phase_windup_enabled: bool,
    pub min_snr_dbhz: f64,
    /// Satellite elevation mask in degrees.
    pub elevation_mask_deg: f64,
    pub dynamics_model: DynamicsModel,
    /// Auto-detect dynamics from observations (overrides dynamics_model).
    pub auto_detect_dynamics: bool,
    pub doppler_slip_threshold_cycles: f64,
    pub max_reject_count: usize,
    pub max_base_age_s: f64,
    pub spp_consistency_threshold_m: f64,
    pub initial_ambiguity_variance: f64,
    pub ar_min_epoch_count: u32,
    pub ar_min_lock: u32,
    pub ar_ffrt_prob: f64,

    // Solver Mode Flags
    /// Ionosphere model selection.
    pub iono_model: IonosphereModel,
    /// Enable horizontal troposphere gradients (adds 2 state params).
    pub enable_tropo_gradients: bool,
    /// Enable integer ambiguity resolution.
    pub enable_ar: bool,

    // Process Noise
    /// Receiver clock bias process noise (m²/s).  For white-noise clock
    /// (phi=0), this is the per-epoch reset variance; for random-walk
    /// clock (phi=1), this is the continuous process noise.
    pub process_noise_cb: f64,
    /// Receiver clock drift process noise (m²/s³).
    pub process_noise_cd: f64,
    /// Inter-system bias process noise (m²/s).  ISBs are modelled as
    /// piece-wise constants with small random-walk drift (~0.3 m/hr
    /// for stable receivers per CODE/WHU analysis).
    pub process_noise_isb: f64,
    pub process_noise_zwd: f64,
    pub process_noise_iono: f64,
    /// Ambiguity float process noise (m²/s).  RALPH value 1e-4 allows
    /// slow re-convergence after cycle slips.  Was 1e-8 (too stiff).
    pub process_noise_amb_float: f64,
    pub process_noise_amb_fixed: f64,

    pub uduc_ar: bool,

    /// Tropospheric mapping function selection.
    pub tropo_mapping: gneiss_core::atmosphere::TropoMapping,

    // GNN RAIM
    pub enable_gnn_raim: bool,
    pub export_gnn_dataset_path: Option<String>,

    // External Tuning configuration
    pub tuning: crate::engine::config::EkfTuningConfig,
}

impl Default for EngineConfig {
    fn default() -> Self {
        Self {
            mode: EngineMode::Rtk,
            initial_position: None,
            base_position: None,
            receiver_antenna_type: None,
            base_datum_transform: None,
            imu_to_antenna_lever_arm: [0.0; 3],
            imu_mounting_angles: None,
            imu_to_nhc_lever_arm: [0.0; 3],
            enable_nhc: false,
            enable_backward_smoothing: false,
            lambda_min_ratio: 1.5,
            lambda_min_subset: 5,
            enabled_constellations: None,
            raim_pseudorange_outlier_m: 25.0,
            chi_square_pr_threshold: 3.0,
            chi_square_cp_threshold: 1e6,
            phase_windup_enabled: true,
            min_snr_dbhz: 25.0,
            elevation_mask_deg: 5.0,
            dynamics_model: DynamicsModel::Static,
            auto_detect_dynamics: false, // TODO: implement auto-detection from velocity estimates
            doppler_slip_threshold_cycles: 5.0,
            max_reject_count: 3,
            max_base_age_s: 5.0,
            spp_consistency_threshold_m: 15.0,
            initial_ambiguity_variance: 10000.0,
            ar_min_epoch_count: 5,
            ar_min_lock: 3,
            ar_ffrt_prob: 0.001,
            iono_model: IonosphereModel::default(),
            enable_tropo_gradients: false,
            enable_ar: false,
            process_noise_cb: 1.0,    // σ=1 m/s for TCXO random walk
            process_noise_cd: 10.0,   // σ=0.55 m/s per epoch (RALPH: was 1e4)
            process_noise_isb: 0.1,   // ~0.3 m/hr random walk
            process_noise_zwd: 1e-8,
            process_noise_iono: 1e-6,
            process_noise_amb_float: 1e-4,   // allows ambiguity re-convergence (RALPH: was 1e-8)
            process_noise_amb_fixed: 1e-12,
            tuning: Default::default(),
            uduc_ar: false,
            tropo_mapping: gneiss_core::atmosphere::TropoMapping::default(),
            enable_gnn_raim: false,
            export_gnn_dataset_path: None,
        }
    }
}
