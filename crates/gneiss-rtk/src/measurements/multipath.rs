//! Code multipath estimation using the MP (Multipath) combination.
//!
//! The MP combination uses dual-frequency carrier phase to isolate code
//! multipath on pseudorange measurements:
//!
//!   MP1 = P1 - L1 - β*(L1 - L2)    where β = 2/(α-1), α = (f1/f2)²
//!
//! This cancels geometry, clocks, troposphere, and ionosphere (first-order).
//! What remains is code multipath on P1 (meter-level) plus carrier-phase
//! multipath (mm-level, negligible) plus a constant integer-ambiguity term.
//!
//! By low-pass filtering MP1 over a sliding window, we estimate the code
//! multipath bias and subtract it from P1 before the EKF update. This is
//! the same technique used by NovAtel, Leica, and other high-end receivers
//! to achieve clean code measurements.

use gneiss_core::obs::{EpochObs, ObsType, SignalCode};
use gneiss_core::sat::SatelliteId;
use std::collections::HashMap;

/// State for one satellite/frequency pair in the multipath estimator.
#[derive(Debug, Clone)]
struct MpState {
    /// Running mean of MP values (converges to ambiguity constant ~10,000m)
    mean: f64,
    /// Number of samples accumulated
    count: usize,
    /// Maximum window size for the slow mean
    max_window: usize,
    /// Current estimated multipath deviation from mean (meters)
    bias: f64,
}

impl MpState {
    fn new(max_window: usize) -> Self {
        Self { mean: 0.0, count: 0, max_window, bias: 0.0 }
    }

    fn push(&mut self, mp_value: f64) {
        if self.count == 0 {
            // Initialize mean to first MP value (captures ambiguity)
            self.mean = mp_value;
            self.count = 1;
        } else {
            // Slow EMA to track the ambiguity baseline.  The mean should
            // converge to the ambiguity constant + mean(multipath).  Since
            // multipath has near-zero mean over diverse geometry, the mean
            // ≈ ambiguity.  The deviation (mp - mean) is the multipath.
            let alpha = 1.0 / (self.max_window as f64).min(self.count as f64 + 1.0);
            self.mean += alpha * (mp_value - self.mean);
            self.count += 1;
            if self.count > self.max_window * 10 {
                self.count = self.max_window; // prevent unbounded growth
            }
        }
        // Bias = deviation from mean = estimated multipath (meters)
        // Apply after collecting enough samples (>10) for a reasonable mean
        if self.count > 10 {
            self.bias = mp_value - self.mean;
        }
    }

    fn reset(&mut self) {
        self.mean = 0.0;
        self.count = 0;
        self.bias = 0.0;
    }
}

/// Code multipath estimator using dual-frequency MP combination.
///
/// Estimates code multipath on L1 and L2 pseudorange using carrier phase
/// from both frequencies. Requires dual-frequency observations.
pub struct MultipathEstimator {
    /// Per-satellite state keyed by (sat, freq_band)
    states: HashMap<(SatelliteId, u8), MpState>,
    /// Sliding window size (epochs). 200 epochs ≈ 40s at 5Hz — long enough
    /// to average code multipath but short enough to track changes.
    window_size: usize,
    /// Cycle slip detection threshold (meters). If a new MP value differs
    /// from the running mean by more than this, the state is reset.
    slip_threshold_m: f64,
}

impl MultipathEstimator {
    pub fn new(window_size: usize) -> Self {
        Self {
            states: HashMap::new(),
            window_size,
            slip_threshold_m: 50.0, // conservative: reset on >50m jump
        }
    }

    /// Apply multipath correction to an EpochObs in-place.
    /// Modifies pseudorange values by subtracting the estimated multipath bias.
    /// Returns the number of satellites/frequencies corrected.
    pub fn correct_epoch(&mut self, epoch: &mut EpochObs) -> usize {
        let mut corrected = 0usize;

        // First pass: find L1 and L2 pseudorange and carrier phase for each satellite
        for sat_obs in &mut epoch.satellites {
            let pr_l1_idx = sat_obs.observations.iter().position(|o| {
                o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 1
            });
            let cp_l1_idx = sat_obs.observations.iter().position(|o| {
                o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1
            });

            if let (Some(pi1), Some(ci1)) = (pr_l1_idx, cp_l1_idx) {
                let pr_l2_idx = sat_obs.observations.iter().position(|o| {
                    o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 2
                });
                let cp_l2_idx = sat_obs.observations.iter().position(|o| {
                    o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 2
                });

                if let (Some(pi2), Some(ci2)) = (pr_l2_idx, cp_l2_idx) {
                    let p1 = sat_obs.observations[pi1].value;
                    let l1 = sat_obs.observations[ci1].value;
                    let p2 = sat_obs.observations[pi2].value;
                    let l2 = sat_obs.observations[ci2].value;

                    // Compute β = 2/(α-1) where α = (f1/f2)²
                    let (f1, f2) = gneiss_core::signal::satellite_frequencies(
                        sat_obs.sat,
                        0, // CDMA-only metric: for GLONASS the f1/f2 ratio is 9/7 for EVERY channel, so nominal k=0 yields identical alpha/beta/gamma (verified over k in [-12,12])
                    );
                    let alpha = (f1 / f2).powi(2);
                    let beta = 2.0 / (alpha - 1.0);

                    // MP1 = P1 - L1 - β*(L1 - L2)
                    let mp1 = p1 - l1 - beta * (l1 - l2);

                    let key1 = (sat_obs.sat, 1u8);
                    let state1 = self.states.entry(key1).or_insert_with(|| {
                        MpState::new(self.window_size)
                    });

                    // Cycle slip check: if MP jumps > threshold, reset
                    if state1.count > 10
                        && (mp1 - state1.bias).abs() > self.slip_threshold_m
                    {
                        state1.reset();
                    }
                    state1.push(mp1);

                    // Correct P1 by subtracting estimated multipath bias
                    sat_obs.observations[pi1].value -= state1.bias;
                    corrected += 1;

                    // MP2 = P2 - L2 + (2*alpha/(alpha-1))*(L1 - L2)  [approximation]
                    // Actually: MP2 = P2 - L2 - (2*alpha/(alpha-1))*(L1 - L2)
                    let gamma = 2.0 * alpha / (alpha - 1.0);
                    let mp2 = p2 - l2 - gamma * (l1 - l2);

                    let key2 = (sat_obs.sat, 2u8);
                    let state2 = self.states.entry(key2).or_insert_with(|| {
                        MpState::new(self.window_size)
                    });

                    if state2.count > 10
                        && (mp2 - state2.bias).abs() > self.slip_threshold_m
                    {
                        state2.reset();
                    }
                    state2.push(mp2);

                    // Correct P2
                    sat_obs.observations[pi2].value -= state2.bias;
                    corrected += 1;
                }
            }
        }

        corrected
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::obs::{ObsCode, Observation, SatObs};
    use gneiss_core::sat::Constellation;

    #[test]
    fn test_mp_combination_smoothes_noise() {
        let mut est = MultipathEstimator::new(100);
        let sig1 = SignalCode { freq_band: 1, attribute: 'C' };
        let sig2 = SignalCode { freq_band: 2, attribute: 'C' };
        let pr1_code = ObsCode { obs_type: ObsType::Pseudorange, signal: sig1 };
        let cp1_code = ObsCode { obs_type: ObsType::CarrierPhase, signal: sig1 };
        let pr2_code = ObsCode { obs_type: ObsType::Pseudorange, signal: sig2 };
        let cp2_code = ObsCode { obs_type: ObsType::CarrierPhase, signal: sig2 };
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };

        // Feed 50 epochs with 5m code multipath bias -> should converge
        for i in 0..150 {
            let p1 = 20_000_000.0 + 5.0 + (i as f64 % 3.0 - 1.5); // 5m bias + noise
            let l1 = 20_000_000.0 + (i as f64 * 1.0_f64); // geometric + clock drift
            let p2 = 20_000_000.0 + 5.5 + (i as f64 % 2.5 - 1.25); // 5.5m bias + noise
            let l2 = 20_000_000.0 + (i as f64 * 1.0_f64); // same geometric + clock drift
            let mut epoch = EpochObs {
                time: gneiss_core::time::GpsTime::new(0, i as f64),
                satellites: vec![SatObs {
                    sat,
                    observations: vec![
                        Observation { code: pr1_code, value: p1, lock_time: None, lli: None },
                        Observation { code: cp1_code, value: l1, lock_time: None, lli: None },
                        Observation { code: pr2_code, value: p2, lock_time: None, lli: None },
                        Observation { code: cp2_code, value: l2, lock_time: None, lli: None },
                    ],
                }],
            };
            est.correct_epoch(&mut epoch);
            let corrected_p1 = epoch.satellites[0].observations[0].value;
            // After convergence, corrected P1 should be close to (p1 - 5.0)
            if i > 50 {
                assert!((corrected_p1 - (20_000_000.0 + (i as f64 % 3.0 - 1.5))).abs() < 2.0,
                    "epoch {}: corrected P1 = {:.3}, expected near {:.3}",
                    i, corrected_p1, 20_000_000.0 + (i as f64 % 3.0 - 1.5));
            }
        }
    }
}
