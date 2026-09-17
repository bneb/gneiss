use serde::{Deserialize, Serialize};

// ---- Type-state configuration: impossible states are unrepresentable ----

/// Top-level engine configuration.  The mode tag determines which variant
/// is deserialized, and each variant carries only its relevant fields.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "mode")]
pub enum EngineConfig {
    #[serde(rename = "spp")]
    Spp(SppConfig),
    #[serde(rename = "ppp")]
    Ppp(PppConfig),
    #[serde(rename = "rtk")]
    Rtk(RtkConfig),
    #[serde(rename = "rtk-ins")]
    RtkIns(RtkInsConfig),
    #[serde(rename = "ppp-ins")]
    PppIns(PppInsConfig),
}

impl Default for EngineConfig {
    fn default() -> Self {
        Self::Rtk(RtkConfig::default())
    }
}

// ---- SPP ----

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct SppConfig {
    pub elevation_mask_deg: f64,
    pub min_snr_dbhz: f64,
    pub initial_position: Option<[f64; 3]>,
}

impl Default for SppConfig {
    fn default() -> Self {
        Self {
            elevation_mask_deg: 15.0,
            min_snr_dbhz: 25.0,
            initial_position: None,
        }
    }
}

// ---- PPP ----

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct PppConfig {
    pub elevation_mask_deg: f64,
    pub min_snr_dbhz: f64,
    pub initial_position: Option<[f64; 3]>,
    pub window_size: usize,
    pub ar: Option<ArConfig>,
    pub is_kinematic: bool,
    pub enable_glonass: bool,
    pub enable_galileo: bool,
    pub initial_pos_sigma_m: Option<f64>,
}

impl Default for PppConfig {
    fn default() -> Self {
        Self {
            elevation_mask_deg: 15.0,
            min_snr_dbhz: 25.0,
            initial_position: None,
            window_size: 10,
            ar: None,
            is_kinematic: false,
            enable_glonass: false,
            enable_galileo: false,
            initial_pos_sigma_m: None,
        }
    }
}

// ---- RTK ----

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtkConfig {
    /// Base station position (ECEF meters).  Required — deserialization
    /// fails if absent.
    pub base_position: [f64; 3],

    /// Maximum age of base observations (seconds).
    #[serde(default = "default_max_base_age")]
    pub max_base_age_s: f64,

    /// Initial rover position (from RINEX header or survey mark).
    #[serde(default)]
    pub initial_position: Option<[f64; 3]>,

    /// Elevation mask (degrees).
    #[serde(default = "default_elevation_mask")]
    pub elevation_mask_deg: f64,

    /// Minimum signal-to-noise ratio (dB-Hz).
    #[serde(default = "default_min_snr")]
    pub min_snr_dbhz: f64,

    /// Sliding window size (epochs).
    #[serde(default = "default_window_size")]
    pub window_size: usize,

    /// Ambiguity resolution configuration.  None = float-only.
    #[serde(default)]
    pub ar: Option<ArConfig>,
}

impl Default for RtkConfig {
    fn default() -> Self {
        Self {
            base_position: [0.0; 3],
            max_base_age_s: 5.0,
            initial_position: None,
            elevation_mask_deg: 15.0,
            min_snr_dbhz: 25.0,
            window_size: 10,
            ar: None,
        }
    }
}

// ---- RTK-INS ----

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtkInsConfig {
    #[serde(flatten)]
    pub rtk: RtkConfig,
    /// IMU configuration.  Required — deserialization fails if absent.
    pub imu: ImuConfig,
}

// ---- PPP-INS ----

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PppInsConfig {
    #[serde(flatten)]
    pub ppp: PppConfig,
    /// IMU configuration.  Required.
    pub imu: ImuConfig,
}

// ---- IMU ----

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ImuConfig {
    /// IMU-to-antenna lever arm in body frame (meters).
    pub lever_arm: [f64; 3],
    /// IMU mounting angles [roll, pitch, yaw] in radians.
    pub mounting_angles: [f64; 3],
    /// Enable non-holonomic constraints.
    pub enable_nhc: bool,
    /// NHC reference point in body frame (meters from IMU).
    pub nhc_lever_arm: [f64; 3],
    /// IMU process noise tuning.
    #[serde(default)]
    pub tuning: ImuTuning,
}

impl Default for ImuConfig {
    fn default() -> Self {
        Self {
            lever_arm: [0.0; 3],
            mounting_angles: [0.0; 3],
            enable_nhc: false,
            nhc_lever_arm: [0.0; 3],
            tuning: ImuTuning::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ImuTuning {
    /// Velocity random walk (m/s/√s).
    pub sigma_v: f64,
    /// Angular random walk (rad/√s).
    pub sigma_phi: f64,
    /// Accel bias instability (m/s²).
    pub sigma_ab: f64,
    /// Gyro bias instability (rad/s).
    pub sigma_gb: f64,
}

impl Default for ImuTuning {
    fn default() -> Self {
        Self {
            sigma_v: 0.1,
            sigma_phi: 0.01,
            sigma_ab: 1e-4,
            sigma_gb: 1e-5,
        }
    }
}

// ---- Ambiguity Resolution ----

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ArConfig {
    /// Minimum LAMBDA ratio to accept a fix.
    #[serde(default = "default_lambda_ratio")]
    pub min_ratio: f64,
    /// Minimum number of satellites in the ambiguity subset.
    pub min_subsets: usize,
    /// Target false-fix probability for FFRT.
    pub ffrt_prob: f64,
}

impl Default for ArConfig {
    fn default() -> Self {
        Self {
            min_ratio: 2.5,
            min_subsets: 5,
            ffrt_prob: 0.001,
        }
    }
}

// ---- Default value helpers ----

const fn default_max_base_age() -> f64 { 5.0 }
const fn default_elevation_mask() -> f64 { 15.0 }
const fn default_min_snr() -> f64 { 25.0 }
const fn default_window_size() -> usize { 10 }
const fn default_lambda_ratio() -> f64 { 2.5 }

// ---- Tests ----

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rtk_config_requires_base_position() {
        let json = r#"{"mode": "rtk", "base_position": [1.0, 2.0, 3.0]}"#;
        let cfg: EngineConfig = serde_json::from_str(json).unwrap();
        match cfg {
            EngineConfig::Rtk(rtk) => {
                assert_eq!(rtk.base_position, [1.0, 2.0, 3.0]);
                assert_eq!(rtk.window_size, 10);
            }
            _ => panic!("expected Rtk"),
        }
    }

    #[test]
    fn rtk_ins_config_requires_imu() {
        let json = r#"{"mode": "rtk-ins", "base_position": [1.0, 2.0, 3.0], "imu": {"lever_arm": [0.5, 0.0, 1.0]}}"#;
        let cfg: EngineConfig = serde_json::from_str(json).unwrap();
        match cfg {
            EngineConfig::RtkIns(rtk_ins) => {
                assert_eq!(rtk_ins.rtk.base_position, [1.0, 2.0, 3.0]);
                assert_eq!(rtk_ins.imu.lever_arm, [0.5, 0.0, 1.0]);
            }
            _ => panic!("expected RtkIns"),
        }
    }

    #[test]
    fn ppp_config_defaults() {
        let json = r#"{"mode": "ppp"}"#;
        let cfg: EngineConfig = serde_json::from_str(json).unwrap();
        match cfg {
            EngineConfig::Ppp(ppp) => {
                assert_eq!(ppp.window_size, 10);
                assert!(ppp.ar.is_none());
            }
            _ => panic!("expected Ppp"),
        }
    }

    #[test]
    fn ar_config_default_ratio() {
        let ar = ArConfig::default();
        assert_eq!(ar.min_ratio, 2.5);
        assert_eq!(ar.min_subsets, 5);
    }

    #[test]
    fn imu_tuning_defaults() {
        let tuning = ImuTuning::default();
        assert_eq!(tuning.sigma_v, 0.1);
        assert_eq!(tuning.sigma_phi, 0.01);
    }
}
