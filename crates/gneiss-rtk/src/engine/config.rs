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

    // IMU Process Noise
    pub sigma_v: f64,    // Velocity Random Walk
    pub sigma_phi: f64,  // Angular Random Walk
    pub sigma_ab: f64,   // Accel bias instability
    pub sigma_gb: f64,   // Gyro bias instability

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

    // Auto-tuning Constraints
    pub auto_tune: AutoTuneConfig,
}

impl Default for EkfTuningConfig {
    fn default() -> Self {
        Self {
            pr_base_var: 0.09,
            cp_base_var: 9e-6,
            dop_base_var: 0.1,
            sigma_v: 0.01,
            sigma_phi: 0.001,
            sigma_ab: 1e-4,
            sigma_gb: 1e-5,
            loosely_coupled_mahalanobis_sq: 250.0,
            phase_outlier_ratio_thresh: 5.0,
            doppler_outlier_ratio_mult: 2.0,
            pr_abs_thresh: 40.0,
            cp_abs_thresh: 1.0,
            dop_abs_thresh: 15.0,
            huber_threshold_loosely: 10.0,
            huber_threshold_tightly: 3.0,
            auto_tune: Default::default(),
        }
    }
}
