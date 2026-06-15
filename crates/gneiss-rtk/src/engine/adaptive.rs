use std::collections::HashMap;
use gneiss_core::sat::SatelliteId;

/// Exponential smoothing constant for NIS tracking.
/// α = 0.95 gives ~20-epoch effective window.
const NIS_SMOOTHING_ALPHA: f64 = 0.95;

/// Minimum scale factor — never deflate below static R.
const MIN_SCALE: f64 = 1.0;

/// Maximum scale factor — cap to prevent filter stalling.
const MAX_SCALE: f64 = 100.0;

/// Per-satellite, per-signal normalized innovation squared (NIS) tracker.
/// Computes adaptive R scaling based on observed innovation statistics.
#[derive(Debug, Clone)]
pub struct InnovationTracker {
    /// Exponentially weighted NIS average per (satellite, frequency).
    nis_avg: HashMap<(SatelliteId, u8), f64>,
}

impl InnovationTracker {
    pub fn new() -> Self {
        Self { nis_avg: HashMap::new() }
    }

    /// Update the NIS tracker with a new innovation and return the adaptive R scale.
    ///
    /// `innovation`: the measurement residual (z - h*x)
    /// `predicted_var`: S_ii = H*P*H' + R (the innovation covariance diagonal)
    ///
    /// Returns a scale factor >= 1.0 to multiply into the static R value.
    pub fn update_and_scale(
        &mut self, sat: SatelliteId, freq: u8,
        innovation: f64, predicted_var: f64,
    ) -> f64 {
        let nis = compute_nis(innovation, predicted_var);
        let key = (sat, freq);
        let avg = self.nis_avg.entry(key).or_insert(1.0);
        *avg = NIS_SMOOTHING_ALPHA * *avg + (1.0 - NIS_SMOOTHING_ALPHA) * nis;
        avg.clamp(MIN_SCALE, MAX_SCALE)
    }

    /// Query current scale factor without updating.
    pub fn current_scale(&self, sat: SatelliteId, freq: u8) -> f64 {
        self.nis_avg.get(&(sat, freq)).copied().unwrap_or(1.0).clamp(MIN_SCALE, MAX_SCALE)
    }

    /// Prune satellites not seen for many epochs.
    pub fn prune(&mut self, active_sats: &[(SatelliteId, u8)]) {
        self.nis_avg.retain(|k, _| active_sats.contains(k));
    }
}

impl Default for InnovationTracker {
    fn default() -> Self { Self::new() }
}

/// Compute normalized innovation squared: v² / S
fn compute_nis(innovation: f64, predicted_var: f64) -> f64 {
    if predicted_var <= 0.0 { return 1.0; }
    (innovation * innovation) / predicted_var
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;

    fn test_sat() -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn: 1 }
    }

    #[test]
    fn nominal_innovations_yield_unit_scale() {
        let mut tracker = InnovationTracker::new();
        let sat = test_sat();
        // innovation = 1.0, predicted_var = 1.0 → NIS = 1.0
        for _ in 0..50 {
            let scale = tracker.update_and_scale(sat, 1, 1.0, 1.0);
            assert!((scale - 1.0).abs() < 0.1, "Nominal NIS should give scale ≈ 1.0, got {}", scale);
        }
    }

    #[test]
    fn large_innovations_inflate_scale() {
        let mut tracker = InnovationTracker::new();
        let sat = test_sat();
        // innovation = 10.0, predicted_var = 1.0 → NIS = 100.0
        for _ in 0..50 {
            tracker.update_and_scale(sat, 1, 10.0, 1.0);
        }
        let scale = tracker.current_scale(sat, 1);
        assert!(scale > 10.0, "Large NIS should inflate scale well above 1, got {}", scale);
    }

    #[test]
    fn scale_decays_after_clean_measurements() {
        let mut tracker = InnovationTracker::new();
        let sat = test_sat();
        // Pump up the NIS with bad measurements
        for _ in 0..20 {
            tracker.update_and_scale(sat, 1, 10.0, 1.0);
        }
        let peak_scale = tracker.current_scale(sat, 1);
        assert!(peak_scale > 5.0);

        // Feed clean measurements and verify decay
        for _ in 0..100 {
            tracker.update_and_scale(sat, 1, 1.0, 1.0);
        }
        let decayed_scale = tracker.current_scale(sat, 1);
        assert!(decayed_scale < peak_scale, "Scale should decay with clean data");
        assert!(decayed_scale < 3.0, "After 100 clean epochs, scale should be near 1, got {}", decayed_scale);
    }

    #[test]
    fn scale_clamped_to_max() {
        let mut tracker = InnovationTracker::new();
        let sat = test_sat();
        // Enormous innovation
        for _ in 0..100 {
            tracker.update_and_scale(sat, 1, 1000.0, 1.0);
        }
        let scale = tracker.current_scale(sat, 1);
        assert!(scale <= MAX_SCALE, "Scale must be clamped to MAX_SCALE");
        assert_eq!(scale, MAX_SCALE);
    }

    #[test]
    fn unknown_satellite_returns_unit_scale() {
        let tracker = InnovationTracker::new();
        let sat = test_sat();
        assert_eq!(tracker.current_scale(sat, 1), 1.0);
    }

    #[test]
    fn independent_satellite_tracking() {
        let mut tracker = InnovationTracker::new();
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 5 };
        // Only sat1 gets bad measurements
        for _ in 0..50 {
            tracker.update_and_scale(sat1, 1, 10.0, 1.0);
            tracker.update_and_scale(sat2, 1, 1.0, 1.0);
        }
        assert!(tracker.current_scale(sat1, 1) > 10.0);
        assert!((tracker.current_scale(sat2, 1) - 1.0).abs() < 0.1);
    }
}
