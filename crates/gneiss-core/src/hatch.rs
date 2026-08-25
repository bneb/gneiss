//! Carrier-smoothed-code (Hatch) filter.
//!
//! Reduces pseudorange noise by √N using carrier phase as a precise
//! increment between epochs. Standard technique in all commercial
//! receivers; improves float ambiguity initialisation quality.

use crate::sat::SatelliteId;
use alloc::vec::Vec;

/// Per-satellite Hatch filter state.
#[derive(Debug, Clone)]
struct HatchState {
    smoothed_code_m: f64,
    prev_phase_cyc: f64,
}

/// Multi-satellite Hatch filter for one receiver.
/// Multi-satellite Hatch filter state (linear scan; ~32 sats typical).
#[derive(Debug, Default)]
pub struct HatchFilter {
    alpha: f64,
    sats: Vec<(SatelliteId, HatchState)>,
}

impl HatchFilter {
    /// Create with smoothing factor alpha = 1/N (N = effective window).
    pub fn new(window_epochs: usize) -> Self {
        Self {
            alpha: 1.0 / window_epochs.max(1) as f64,
            sats: Vec::new(),
        }
    }

    /// Apply filter to one satellite observation.
    ///
    /// Returns the smoothed pseudorange. Falls back to raw code on first
    /// observation or after a cycle slip (phase discontinuity > threshold).
    pub fn apply(
        &mut self,
        sat: SatelliteId,
        raw_code_m: f64,
        phase_cyc: f64,
        lambda_m: f64,
    ) -> f64 {
        const SLIP_THRESHOLD_CYC: f64 = 10.0;
        if let Some(idx) = self.sats.iter().position(|(s, _)| *s == sat) {
            let st = &mut self.sats[idx].1;
            let d_phase = phase_cyc - st.prev_phase_cyc;
            if d_phase.abs() > SLIP_THRESHOLD_CYC {
                // Cycle slip: reseed from raw code.
                *st = HatchState { smoothed_code_m: raw_code_m, prev_phase_cyc: phase_cyc };
                return raw_code_m;
            }
            let predicted = st.smoothed_code_m + d_phase * lambda_m;
            let smoothed = self.alpha * raw_code_m + (1.0 - self.alpha) * predicted;
            st.smoothed_code_m = smoothed;
            st.prev_phase_cyc = phase_cyc;
            smoothed
        } else {
            self.sats.push((sat, HatchState {
                smoothed_code_m: raw_code_m,
                prev_phase_cyc: phase_cyc,
            }));
            raw_code_m
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sat::Constellation;

    fn mk_sat(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn }
    }

    #[test]
    fn test_first_observation_returns_raw() {
        let mut hf = HatchFilter::new(10);
        let s = mk_sat(1);
        assert_eq!(hf.apply(s, 20_000_000.0, 100_000.0, 0.19), 20_000_000.0);
    }

    #[test]
    fn test_noise_reduction_converges() {
        let mut hf = HatchFilter::new(100); // α = 0.01
        let s = mk_sat(5);
        let lambda = 0.190_293_672_798_364_64; // L1
        // Physically consistent simulation: range advances with the phase.
        let mut true_range = 20_000_000.0;
        let mut phase = 100_000.0;
        let _ = hf.apply(s, true_range, phase, lambda);
        let mut last = true_range;
        for i in 1..200 {
            true_range += lambda; // consistent with +1 cycle phase advance
            phase += 1.0;
            let noise = if i % 2 == 0 { 3.0 } else { -3.0 };
            last = hf.apply(s, true_range + noise, phase, lambda);
        }
        assert!(
            (last - true_range).abs() < 0.5,
            "smoothed {} vs true {}, residual {:.3}",
            last, true_range, (last - true_range).abs()
        );
    }

    #[test]
    fn test_cycle_slip_reseeds() {
        let mut hf = HatchFilter::new(10);
        let s = mk_sat(3);
        let lambda = 0.19;
        let _ = hf.apply(s, 20_000_000.0, 50_000.0, lambda);
        // Large phase jump = cycle slip → returns raw code
        let result = hf.apply(s, 20_000_100.0, 60_000.0, lambda);
        assert_eq!(result, 20_000_100.0, "slip should reseed from raw");
    }
}
