//! Multi-base RTK position combiner.
//!
//! When multiple base stations are available, each baseline produces an
//! independent RTK position fix. The combiner merges them via weighted
//! averaging, where shorter baselines (less ionospheric decorrelation)
//! receive higher weight.
//!
//! Weights: w_i = 1 / (d_i² + σ₀²)
//! where d_i = baseline distance in meters, σ₀ = base uncertainty (~1m)

use nalgebra::Vector3;

/// A position fix from a single base station.
#[derive(Clone, Debug)]
pub struct BaseFix {
    /// ECEF position from this baseline
    pub position: Vector3<f64>,
    /// Baseline distance in meters (shorter = more reliable)
    pub baseline_distance_m: f64,
    /// Position quality indicator (0-1, 1=best, e.g. AR fix ratio)
    pub quality: f64,
}

/// Multi-base RTK position combiner.
pub struct MultiBaseCombiner {
    fixes: Vec<BaseFix>,
    /// Base uncertainty floor (meters)
    sigma_0: f64,
}

impl MultiBaseCombiner {
    pub fn new(sigma_0: f64) -> Self {
        Self { fixes: Vec::new(), sigma_0 }
    }

    /// Add a position fix from one base station.
    pub fn add_fix(&mut self, position: Vector3<f64>, baseline_m: f64, quality: f64) {
        self.fixes.push(BaseFix { position, baseline_distance_m: baseline_m, quality });
    }

    /// Number of contributing base stations.
    pub fn num_bases(&self) -> usize { self.fixes.len() }

    /// Weight for a fix: higher for closer bases and higher quality.
    fn weight(&self, fix: &BaseFix) -> f64 {
        let d = fix.baseline_distance_m.max(100.0); // min 100m to prevent blowup
        fix.quality / (d * d + self.sigma_0 * self.sigma_0)
    }

    /// Combined position via weighted average.
    pub fn weighted_position(&self) -> Option<Vector3<f64>> {
        if self.fixes.is_empty() { return None; }
        if self.fixes.len() == 1 { return Some(self.fixes[0].position); }

        let total_weight: f64 = self.fixes.iter().map(|f| self.weight(f)).sum();
        if total_weight <= 0.0 { return None; }

        let mut pos = Vector3::zeros();
        for fix in &self.fixes {
            let w = self.weight(fix) / total_weight;
            pos += fix.position * w;
        }
        Some(pos)
    }

    /// Robust position: drop the worst fix and re-combine if ≥3 bases.
    pub fn robust_position(&self) -> Option<Vector3<f64>> {
        if self.fixes.len() >= 3 {
            // Compute distances from weighted average, drop farthest
            let avg = self.weighted_position()?;
            let mut dists: Vec<(usize, f64)> = self.fixes.iter().enumerate()
                .map(|(i, f)| (i, (f.position - avg).norm()))
                .collect();
            dists.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            let worst_idx = dists[0].0;

            let mut trimmed = MultiBaseCombiner::new(self.sigma_0);
            for (i, fix) in self.fixes.iter().enumerate() {
                if i != worst_idx { trimmed.add_fix(fix.position, fix.baseline_distance_m, fix.quality); }
            }
            trimmed.weighted_position()
        } else {
            self.weighted_position()
        }
    }

    /// Expected accuracy improvement factor from N bases.
    /// sqrt(1/N) improvement from averaging independent estimates.
    pub fn expected_improvement(&self) -> f64 {
        let n = self.fixes.len().max(1) as f64;
        1.0 / n.sqrt()
    }
}

impl Default for MultiBaseCombiner {
    fn default() -> Self { Self::new(1.0) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_base_returns_direct() {
        let mut c = MultiBaseCombiner::new(1.0);
        c.add_fix(Vector3::new(1.0, 2.0, 3.0), 5000.0, 1.0);
        let pos = c.weighted_position().unwrap();
        assert!((pos - Vector3::new(1.0, 2.0, 3.0)).norm() < 1e-9);
    }

    #[test]
    fn test_closer_base_gets_higher_weight() {
        let mut c = MultiBaseCombiner::new(1.0);
        // Close base (5km): should dominate
        c.add_fix(Vector3::new(1.0, 0.0, 0.0), 5_000.0, 1.0);
        // Far base (50km): should have little influence
        c.add_fix(Vector3::new(10.0, 0.0, 0.0), 50_000.0, 1.0);

        let pos = c.weighted_position().unwrap();
        // Close base weight = 1/25e6 ≈ 4e-8, Far base weight = 1/2.5e9 ≈ 4e-10
        // Ratio ~100:1. Expected position near 1.0, not 5.5
        assert!(pos.x < 1.5, "Close base should dominate: got x={:.3}", pos.x);
    }

    #[test]
    fn test_two_bases_improve_over_single() {
        let truth = Vector3::new(1_000_000.0, 2_000_000.0, 3_000_000.0);

        // Simple LCG for deterministic noise
        let mut seed: u64 = 99;
        let mut noise = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 12) as f64 / (u64::MAX >> 12) as f64 - 0.5) * 2.0
        };

        // Base at 10km with σ=0.5m noise
        let mut c = MultiBaseCombiner::new(1.0);
        c.add_fix(truth + Vector3::new(noise(), noise(), noise()) * 0.5, 10_000.0, 1.0);
        let err_1 = (c.weighted_position().unwrap() - truth).norm();

        // Add second base at 15km
        c.add_fix(truth + Vector3::new(noise(), noise(), noise()) * 0.5, 15_000.0, 1.0);
        let err_2 = (c.weighted_position().unwrap() - truth).norm();

        // Two bases should give better or equal accuracy
        assert!(err_2 <= err_1 * 1.1,
            "Two bases ({:.3}m) should not be worse than one ({:.3}m)", err_2, err_1);
    }

    #[test]
    fn test_robust_drops_outlier() {
        let mut c = MultiBaseCombiner::new(1.0);
        // Two good bases (10km, 12km)
        c.add_fix(Vector3::new(0.0, 0.0, 0.0), 10_000.0, 1.0);
        c.add_fix(Vector3::new(0.1, -0.1, 0.0), 12_000.0, 1.0);
        // Outlier: far base (50km) with bad position — should get low weight
        c.add_fix(Vector3::new(50.0, 50.0, -50.0), 50_000.0, 0.5);

        let robust = c.robust_position().unwrap();
        let weighted = c.weighted_position().unwrap();

        // Robust should be closer to truth (0,0,0) than weighted
        let r_err = robust.norm();
        let w_err = weighted.norm();
        assert!(r_err < w_err,
            "Robust ({:.3}m) should beat weighted ({:.3}m) with outlier", r_err, w_err);
    }

    #[test]
    fn test_expected_improvement() {
        let mut c = MultiBaseCombiner::new(1.0);
        assert!((c.expected_improvement() - 1.0).abs() < 1e-9);
        c.add_fix(Vector3::zeros(), 1.0, 1.0);
        c.add_fix(Vector3::zeros(), 1.0, 1.0);
        c.add_fix(Vector3::zeros(), 1.0, 1.0);
        c.add_fix(Vector3::zeros(), 1.0, 1.0);
        assert!((c.expected_improvement() - 0.5).abs() < 0.01);
    }
}
