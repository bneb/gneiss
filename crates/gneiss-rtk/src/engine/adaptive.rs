use gneiss_core::sat::SatelliteId;
use std::collections::HashMap;

/// Exponential smoothing constant for NIS and SNR tracking.
const SMOOTHING_ALPHA: f64 = 0.95;

/// Minimum scale factor — never deflate below static R.
const MIN_SCALE: f64 = 1.0;

/// Maximum scale factor — cap to prevent filter stalling and masking outliers.
const MAX_SCALE: f64 = 3.0;

/// Threshold for SNR variance to trigger urban canyon penalty
const SNR_VAR_THRESHOLD: f64 = 10.0;

/// Per-satellite, per-signal tracking for NIS and SNR.
/// Computes adaptive R scaling based on innovation and signal statistics.
#[derive(Debug, Clone)]
pub struct InnovationTracker {
    /// Exponentially weighted NIS average per (satellite, frequency).
    nis_avg: HashMap<(SatelliteId, u8), f64>,
    /// Exponentially weighted SNR average
    snr_avg: HashMap<(SatelliteId, u8), f64>,
    /// Exponentially weighted SNR variance
    snr_var: HashMap<(SatelliteId, u8), f64>,
}

impl InnovationTracker {
    pub fn new() -> Self {
        Self {
            nis_avg: HashMap::new(),
            snr_avg: HashMap::new(),
            snr_var: HashMap::new(),
        }
    }

    /// Update the tracker with a new innovation and SNR, return the adaptive R scale.
    ///
    /// `innovation`: the measurement residual (z - h*x)
    /// `predicted_var`: S_ii = H*P*H' + R
    /// `snr`: Current Signal-to-Noise Ratio (dBHz)
    ///
    /// Returns a scale factor >= 1.0 to multiply into the static R value.
    pub fn update_and_scale(
        &mut self,
        sat: SatelliteId,
        freq: u8,
        innovation: f64,
        predicted_var: f64,
        snr: Option<f64>,
    ) -> f64 {
        let nis = compute_nis(innovation, predicted_var);
        let key = (sat, freq);

        let avg = self.nis_avg.entry(key).or_insert(1.0);
        *avg = SMOOTHING_ALPHA * *avg + (1.0 - SMOOTHING_ALPHA) * nis;

        let mut snr_penalty = 1.0;
        if let Some(snr_val) = snr {
            let s_avg = self.snr_avg.entry(key).or_insert(snr_val);
            let s_var = self.snr_var.entry(key).or_insert(0.0);

            let diff = snr_val - *s_avg;
            *s_avg = SMOOTHING_ALPHA * *s_avg + (1.0 - SMOOTHING_ALPHA) * snr_val;
            *s_var = SMOOTHING_ALPHA * *s_var + (1.0 - SMOOTHING_ALPHA) * (diff * diff);

            if *s_var > SNR_VAR_THRESHOLD {
                snr_penalty = 1.0 + (*s_var - SNR_VAR_THRESHOLD) * 0.2;
            }
        }

        (*avg * snr_penalty).clamp(MIN_SCALE, MAX_SCALE)
    }

    /// Query current scale factor without updating.
    pub fn current_scale(&self, sat: SatelliteId, freq: u8) -> f64 {
        self.nis_avg
            .get(&(sat, freq))
            .copied()
            .unwrap_or(1.0)
            .clamp(MIN_SCALE, MAX_SCALE)
    }

    /// Return the maximum SNR variance across all tracked satellites.
    /// Useful for global urban canyon detection.
    pub fn max_snr_variance(&self) -> f64 {
        self.snr_var.values().copied().fold(0.0, f64::max)
    }

    pub fn is_stalling(&self) -> bool {
        let avg_scale = self.nis_avg.values().sum::<f64>() / self.nis_avg.len().max(1) as f64;
        avg_scale > 10.0
    }

    pub fn get_total_nis(&self) -> f64 {
        self.nis_avg.values().sum::<f64>()
    }

    /// Prune satellites not seen for many epochs.
    pub fn prune(&mut self, active_sats: &[(SatelliteId, u8)]) {
        self.nis_avg.retain(|k, _| active_sats.contains(k));
        self.snr_avg.retain(|k, _| active_sats.contains(k));
        self.snr_var.retain(|k, _| active_sats.contains(k));
    }
}

impl Default for InnovationTracker {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute normalized innovation squared: v² / S
fn compute_nis(innovation: f64, predicted_var: f64) -> f64 {
    if predicted_var <= 0.0 {
        return 1.0;
    }
    (innovation * innovation) / predicted_var
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;

    fn test_sat() -> SatelliteId {
        SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        }
    }

    #[test]
    fn nominal_innovations_yield_unit_scale() {
        let mut tracker = InnovationTracker::new();
        let sat = test_sat();
        for _ in 0..50 {
            let scale = tracker.update_and_scale(sat, 1, 1.0, 1.0, Some(45.0));
            assert!((scale - 1.0).abs() < 0.1);
        }
    }

    #[test]
    fn large_innovations_inflate_scale() {
        let mut tracker = InnovationTracker::new();
        let sat = test_sat();
        for _ in 0..50 {
            tracker.update_and_scale(sat, 1, 10.0, 1.0, Some(45.0));
        }
        assert!(tracker.current_scale(sat, 1) > 2.5);
    }

    #[test]
    fn snr_variance_inflates_scale() {
        let mut tracker = InnovationTracker::new();
        let sat = test_sat();
        // Constant SNR (no penalty)
        let mut scale1 = 1.0;
        for _ in 0..50 {
            scale1 = tracker.update_and_scale(sat, 1, 1.0, 1.0, Some(40.0));
        }

        let mut tracker_var = InnovationTracker::new();
        // High variance SNR (fading/multipath)
        let mut scale2 = 1.0;
        for i in 0..50 {
            let snr = if i % 2 == 0 { 40.0 } else { 20.0 };
            scale2 = tracker_var.update_and_scale(sat, 1, 1.0, 1.0, Some(snr));
        }

        assert!(
            scale2 > scale1 * 2.0,
            "High SNR variance should heavily penalize scale: scale2={}, scale1={}",
            scale2,
            scale1
        );
        assert!(tracker_var.max_snr_variance() > 15.0);
    }
}
