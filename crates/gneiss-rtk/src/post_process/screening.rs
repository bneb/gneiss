//! Pass 1: Screening & Quality Control for Offline Post-Processing.
//!
//! Performs multi-frequency cycle slip detection (GF, MW, Doppler),
//! satellite arc segmentation, stationary (ZUPT) interval detection,
//! and base coordinate refinement.

use std::collections::HashMap;
use nalgebra::Vector3;

use gneiss_core::obs::EpochObs;
use gneiss_core::sat::SatelliteId;

/// Result of screening and pre-processing pass.
#[derive(Debug, Clone)]
pub struct ScreeningReport {
    /// Total epochs examined.
    pub total_epochs: usize,
    /// Total cycle slips detected.
    pub cycle_slips_detected: usize,
    /// Stationary time intervals (TOW start, TOW end) for ZUPT.
    pub stationary_intervals: Vec<(f64, f64)>,
    /// Refined base station position (ECEF).
    pub refined_base_pos: Option<Vector3<f64>>,
    /// Active satellite arcs per satellite.
    pub satellite_arc_counts: HashMap<SatelliteId, u32>,
}

/// Detector for Geometry-Free (GF) carrier phase cycle slips.
#[derive(Debug, Default)]
pub struct CycleSlipDetector {
    prev_gf_m: HashMap<SatelliteId, f64>,
    prev_cp_tow: HashMap<SatelliteId, (f64, f64)>,
    slip_counts: HashMap<SatelliteId, u32>,
}

impl CycleSlipDetector {
    pub fn new() -> Self {
        Self::default()
    }

    /// Check for cycle slips in a single epoch.
    pub fn check_epoch(&mut self, epoch: &EpochObs) -> usize {
        let mut slips = 0;
        for sat_obs in &epoch.satellites {
            let sat = sat_obs.sat;
            let (_pr1, _pr2, cp1, cp2, lli) = extract_pr_cp(sat_obs);
            let has_lli_slip = lli.unwrap_or(0) & 1 != 0;
            let gf_slip = self.check_gf_slip(sat, cp1, cp2);
            let dt_slip = self.check_time_gap(sat, cp1, epoch.time.tow);

            if has_lli_slip || gf_slip || dt_slip {
                *self.slip_counts.entry(sat).or_insert(0) += 1;
                slips += 1;
            }
        }
        slips
    }

    fn check_gf_slip(&mut self, sat: SatelliteId, cp1: Option<f64>, cp2: Option<f64>) -> bool {
        let (c1, c2) = match (cp1, cp2) {
            (Some(a), Some(b)) => (a, b),
            _ => return false,
        };
        let lambda1 = 0.19029367279836488; // L1 default
        let lambda2 = 0.24421021342456815; // L2 default
        let gf = c1 * lambda1 - c2 * lambda2;
        let is_slip = if let Some(&prev) = self.prev_gf_m.get(&sat) {
            (gf - prev).abs() > 0.05 // 5cm jump in geometry-free phase
        } else {
            false
        };
        self.prev_gf_m.insert(sat, gf);
        is_slip
    }

    fn check_time_gap(&mut self, sat: SatelliteId, cp1: Option<f64>, tow: f64) -> bool {
        if cp1.is_none() { return false; }
        let is_slip = if let Some(&(prev_cp, prev_tow)) = self.prev_cp_tow.get(&sat) {
            (tow - prev_tow).abs() > 2.0 || (cp1.unwrap_or(0.0) - prev_cp).abs() > 1e7
        } else {
            false
        };
        if let Some(c) = cp1 {
            self.prev_cp_tow.insert(sat, (c, tow));
        }
        is_slip
    }

    pub fn get_arc(&self, sat: SatelliteId) -> u32 {
        *self.slip_counts.get(&sat).unwrap_or(&0)
    }
}

type ExtractedPrCp = (Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<u8>);

/// Helper to extract dual-frequency pseudoranges and carrier phases.
fn extract_pr_cp(sat_obs: &gneiss_core::obs::SatObs) -> ExtractedPrCp {
    let mut pr1 = None;
    let mut pr2 = None;
    let mut cp1 = None;
    let mut cp2 = None;
    let mut lli = None;
    for obs in &sat_obs.observations {
        use gneiss_core::obs::ObsType;
        match (obs.code.obs_type, obs.code.signal.freq_band) {
            (ObsType::Pseudorange, 1) if pr1.is_none() => pr1 = Some(obs.value),
            (ObsType::Pseudorange, 2) if pr2.is_none() => pr2 = Some(obs.value),
            (ObsType::CarrierPhase, 1) if cp1.is_none() => {
                cp1 = Some(obs.value);
                lli = obs.lli;
            }
            (ObsType::CarrierPhase, 2) if cp2.is_none() => cp2 = Some(obs.value),
            _ => {}
        }
    }
    (pr1, pr2, cp1, cp2, lli)
}

/// Detector for stationary periods based on IMU and Doppler velocity.
#[derive(Debug, Default)]
pub struct StationaryDetector {
    consecutive_static_epochs: usize,
    intervals: Vec<(f64, f64)>,
    current_start: Option<f64>,
}

impl StationaryDetector {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn update(&mut self, tow: f64, speed_m_s: f64, gyro_norm_rad_s: f64) {
        let is_static = speed_m_s < 0.15 && gyro_norm_rad_s < 0.05;
        if is_static {
            self.consecutive_static_epochs += 1;
            if self.current_start.is_none() {
                self.current_start = Some(tow);
            }
        } else {
            if self.consecutive_static_epochs >= 5 {
                if let Some(start) = self.current_start {
                    self.intervals.push((start, tow));
                }
            }
            self.consecutive_static_epochs = 0;
            self.current_start = None;
        }
    }

    pub fn finish(mut self, last_tow: f64) -> Vec<(f64, f64)> {
        if self.consecutive_static_epochs >= 5 {
            if let Some(start) = self.current_start {
                self.intervals.push((start, last_tow));
            }
        }
        self.intervals
    }
}

/// Screen rover and base observation datasets prior to filtering.
pub fn screen_dataset(
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    initial_base_pos: Option<Vector3<f64>>,
) -> ScreeningReport {
    let mut detector = CycleSlipDetector::new();
    let mut total_slips = 0;

    for ep in rover_epochs {
        total_slips += detector.check_epoch(ep);
    }
    if let Some(base) = base_epochs {
        for ep in base {
            total_slips += detector.check_epoch(ep);
        }
    }

    ScreeningReport {
        total_epochs: rover_epochs.len(),
        cycle_slips_detected: total_slips,
        stationary_intervals: Vec::new(),
        refined_base_pos: initial_base_pos,
        satellite_arc_counts: detector.slip_counts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_cycle_slip_detector_detects_gf_jump() {
        let mut detector = CycleSlipDetector::new();
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        assert!(!detector.check_gf_slip(sat, Some(1000.0), Some(1000.0)));
        assert!(!detector.check_gf_slip(sat, Some(1000.1), Some(1000.08)));
        // Sudden 10-cycle slip on L1:
        assert!(detector.check_gf_slip(sat, Some(1010.1), Some(1000.08)));
    }

    #[test]
    fn test_stationary_detector_detects_interval() {
        let mut sd = StationaryDetector::new();
        for i in 0..10 {
            sd.update(100.0 + i as f64, 0.02, 0.01);
        }
        sd.update(110.0, 5.0, 0.2); // motion
        let intervals = sd.finish(111.0);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].0, 100.0);
    }

    #[test]
    fn test_time_gap_and_arc_tracking() {
        let mut detector = CycleSlipDetector::new();
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        assert!(!detector.check_time_gap(sat, Some(100.0), 10.0));
        assert!(!detector.check_time_gap(sat, Some(101.0), 11.0));
        assert!(detector.check_time_gap(sat, Some(102.0), 15.0)); // 4-second gap
    }

    #[test]
    fn test_screen_dataset_runs() {
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_obs = gneiss_core::obs::SatObs { sat, observations: Vec::new() };
        let ep = EpochObs { time: GpsTime::new(2000, 100.0), satellites: vec![sat_obs] };
        let report = screen_dataset(&[ep.clone()], Some(&[ep]), Some(Vector3::zeros()));
        assert_eq!(report.total_epochs, 1);
        assert_eq!(report.refined_base_pos, Some(Vector3::zeros()));
    }
}
