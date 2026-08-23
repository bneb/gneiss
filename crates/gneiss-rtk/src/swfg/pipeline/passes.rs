//! Built-in observation correction passes for the factor graph pipeline.

use nalgebra::Vector3;

/// A raw GNSS observation before any corrections are applied.
/// Correction passes mutate fields in-place as they are applied.
#[derive(Debug, Clone)]
pub struct RawObservation {
    pub satellite: u16,
    pub constellation_id: u8,
    pub pr_l1: f64,
    pub pr_l2: Option<f64>,
    pub cp_l1: Option<f64>,
    pub cp_l1_lli: Option<u8>,
    pub cp_l2: Option<f64>,
    pub doppler: f64,
    pub snr_dbhz: f64,
    pub sat_pos_ecef: Vector3<f64>,
    pub sat_vel_ecef: Vector3<f64>,
    /// Satellite clock bias (meters). Set by clock correction pass.
    pub sat_clock_m: f64,
    pub f1: f64,
    pub f2: f64,
    pub freq_num: i8,
    pub elevation_rad: f64,
    pub azimuth_rad: f64,
    // Fields set by correction passes:
    pub tropo_dry_m: f64,
    pub tropo_map_wet: f64,
    pub iono_l1_m: f64,
    pub variance_m2: f64,
    pub cp_variance_m2: f64,
}

/// A corrected observation ready to be converted into a factor.
#[derive(Debug, Clone)]
pub struct CorrectedObservation {
    pub satellite: u16,
    pub constellation_id: u8,
    pub pr_l1: f64,
    pub pr_l2: Option<f64>,
    pub cp_l1: Option<f64>,
    pub cp_l1_lli: Option<u8>,
    pub cp_l2: Option<f64>,
    pub doppler: f64,
    pub snr_dbhz: f64,
    pub sat_pos_ecef: Vector3<f64>,
    pub sat_clock_m: f64,
    pub f1: f64,
    pub f2: f64,
    pub freq_num: i8,
    pub elevation_rad: f64,
    /// Troposphere dry delay (meters). Set by tropo pass.
    pub tropo_dry_m: f64,
    /// Troposphere wet mapping function value.
    pub tropo_map_wet: f64,
    /// Ionosphere L1 delay (meters). Set by iono pass.
    pub iono_l1_m: f64,
    /// Measurement variance (meters²). Set by variance model pass.
    pub variance_m2: f64,
    /// Carrier phase variance (meters²).
    pub cp_variance_m2: f64,
}

/// Receiver state needed by correction passes.
#[derive(Debug, Clone)]
pub struct ReceiverState {
    /// Rover position (ECEF, meters).
    pub position_ecef: Vector3<f64>,
    /// Receiver clock bias per constellation (meters).
    pub clock_bias_m: Vec<f64>,
    /// Troposphere zenith wet delay (meters).
    pub zwd_m: f64,
    /// GLONASS IFB slope (m/freq_num).
    pub ifb_glo: f64,
    /// Receiver LLH (radians, radians, meters).
    pub llh_rad: Vector3<f64>,
    /// GPS time for ionosphere model evaluation.
    pub time: gneiss_core::time::GpsTime,
}

/// A composable correction applied to a raw observation.
pub trait CorrectionPass: std::fmt::Debug {
    /// Apply this correction to a single observation.
    fn apply(&self, obs: &mut RawObservation, rx_state: &ReceiverState);

    /// Human-readable name for diagnostics.
    fn name(&self) -> &str;
}

/// Satellite clock correction from broadcast ephemeris.
#[derive(Debug, Clone, Default)]
pub struct BroadcastClockCorrection;

impl CorrectionPass for BroadcastClockCorrection {
    fn apply(&self, _obs: &mut RawObservation, _rx: &ReceiverState) {}

    fn name(&self) -> &str {
        "broadcast_clock"
    }
}

/// Satellite clock correction from precise products (SP3/CLK).
#[derive(Debug, Clone)]
pub struct PreciseClockCorrection {
    pub clock_biases_m: Vec<f64>,
    pub constellation_id: u8,
}

impl CorrectionPass for PreciseClockCorrection {
    fn apply(&self, obs: &mut RawObservation, _rx: &ReceiverState) {
        if obs.constellation_id != self.constellation_id {
            return;
        }
        let prn = obs.satellite as usize;
        if prn < self.clock_biases_m.len() {
            obs.sat_clock_m = self.clock_biases_m[prn];
        }
    }

    fn name(&self) -> &str {
        "precise_clock"
    }
}

/// Saastamoinen troposphere model + mapping.
#[derive(Debug, Clone, Default)]
pub struct SaastamoinenTropo;

impl CorrectionPass for SaastamoinenTropo {
    fn apply(&self, obs: &mut RawObservation, rx: &ReceiverState) {
        let h_m = rx.llh_rad.z;
        let p_hpa = 1013.25 * (1.0 - 2.2557e-5 * h_m).powi(5);
        let zdry_m = 0.002277 * p_hpa;
        let m_wet = 1.001 / (0.002001 + f64::sin(obs.elevation_rad).powi(2)).sqrt();

        obs.tropo_dry_m = zdry_m * m_wet;
        obs.tropo_map_wet = m_wet;
    }

    fn name(&self) -> &str {
        "saastamoinen_tropo"
    }
}

/// Klobuchar ionosphere model.
#[derive(Debug, Clone, Default)]
pub struct KlobucharIono {
    pub alpha: [f64; 4],
    pub beta: [f64; 4],
}

impl CorrectionPass for KlobucharIono {
    fn apply(&self, obs: &mut RawObservation, rx: &ReceiverState) {
        use gneiss_core::atmosphere::AtmosphereModel;
        let params = gneiss_core::atmosphere::KlobucharParams {
            alpha: self.alpha,
            beta: self.beta,
        };
        obs.iono_l1_m = AtmosphereModel::iono_klobuchar(
            &params,
            rx.llh_rad,
            obs.azimuth_rad,
            obs.elevation_rad,
            rx.time,
        );
    }

    fn name(&self) -> &str {
        "klobuchar_iono"
    }
}

/// C/N0-based measurement variance model.
#[derive(Debug, Clone)]
pub struct SnrVarianceModel {
    pub pr_base_var: f64,
    pub cp_base_var: f64,
    pub snr_a: f64,
    pub snr_b: f64,
}

impl Default for SnrVarianceModel {
    fn default() -> Self {
        Self {
            pr_base_var: 0.5,
            cp_base_var: 9e-6,
            snr_a: 1.0,
            snr_b: 150.0,
        }
    }
}

impl CorrectionPass for SnrVarianceModel {
    fn apply(&self, obs: &mut RawObservation, _rx: &ReceiverState) {
        let snr = obs.snr_dbhz;
        let cn0_factor = self.snr_a.powi(2) + self.snr_b.powi(2) / 10.0_f64.powf(snr / 10.0);
        let el_factor = 1.0 / f64::sin(obs.elevation_rad).max(0.1);
        obs.variance_m2 = self.pr_base_var * cn0_factor * el_factor;
        obs.cp_variance_m2 = self.cp_base_var * el_factor;
    }

    fn name(&self) -> &str {
        "snr_variance"
    }
}
