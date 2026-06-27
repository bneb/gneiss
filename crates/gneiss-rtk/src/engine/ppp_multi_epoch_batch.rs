//! Multi-epoch batch PPP solver.
//!
//! For static receivers, position is constant across epochs. We solve
//! each epoch independently then combine via inverse-covariance-weighted
//! average. This is simpler and more robust than a full batch SVD, and
//! avoids the conditioning issues of stacking 100+ measurements.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

/// Weighted position estimate with covariance.
#[derive(Clone, Debug)]
pub struct WeightedPosition {
    pub pos: Vector3<f64>,
    pub cov: Matrix3<f64>, // 3×3 position covariance
    pub weight: f64,       // trace-based weight for averaging
}

/// Multi-epoch position combiner. Takes single-epoch position estimates
/// and combines them via inverse-covariance-weighted averaging.
pub struct MultiEpochCombiner {
    positions: Vec<WeightedPosition>,
}

impl MultiEpochCombiner {
    pub fn new() -> Self {
        Self { positions: Vec::new() }
    }

    /// Add a single-epoch position estimate with its covariance.
    pub fn add_epoch(&mut self, pos: Vector3<f64>, cov_3x3: Matrix3<f64>) {
        let weight = 1.0 / (cov_3x3[(0,0)] + cov_3x3[(1,1)] + cov_3x3[(2,2)]).max(1e-6);
        self.positions.push(WeightedPosition { pos, cov: cov_3x3, weight });
    }

    /// Compute the weighted average position from all epochs.
    pub fn weighted_average(&self) -> Option<Vector3<f64>> {
        if self.positions.is_empty() { return None; }

        let total_weight: f64 = self.positions.iter().map(|p| p.weight).sum();
        if total_weight <= 0.0 { return None; }

        let mut avg = Vector3::zeros();
        for p in &self.positions {
            avg += p.pos * (p.weight / total_weight);
        }
        Some(avg)
    }

    /// Compute median position (robust to outliers).
    pub fn median_position(&self) -> Option<Vector3<f64>> {
        if self.positions.is_empty() { return None; }
        let n = self.positions.len();

        let mut xs: Vec<f64> = self.positions.iter().map(|p| p.pos.x).collect();
        let mut ys: Vec<f64> = self.positions.iter().map(|p| p.pos.y).collect();
        let mut zs: Vec<f64> = self.positions.iter().map(|p| p.pos.z).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
        zs.sort_by(|a, b| a.partial_cmp(b).unwrap());

        Some(Vector3::new(xs[n/2], ys[n/2], zs[n/2]))
    }
}

impl Default for MultiEpochCombiner {
    fn default() -> Self { Self::new() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weighted_average_two_epochs() {
        let mut combiner = MultiEpochCombiner::new();
        // Add a precise epoch
        combiner.add_epoch(
            Vector3::new(1.0, 0.0, 0.0),
            Matrix3::identity() * 0.01, // σ=0.1m, weight=333
        );
        // Add a noisy epoch
        combiner.add_epoch(
            Vector3::new(2.0, 0.0, 0.0),
            Matrix3::identity() * 1.0, // σ=1m, weight=0.33
        );

        let avg = combiner.weighted_average().unwrap();
        // Precise epoch dominates: should be near 1.0, not 1.5
        assert!((avg.x - 1.0).abs() < 0.15,
            "Weighted avg should favor precise epoch: got {:.3}", avg.x);
    }

    #[test]
    fn test_weighted_average_converges_to_truth() {
        let truth = Vector3::new(1_000_000.0, 2_000_000.0, 3_000_000.0);
        let mut combiner = MultiEpochCombiner::new();

        // Simple LCG for deterministic "noise"
        let mut seed: u64 = 42;
        let mut noise = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 12) as f64 / (u64::MAX >> 12) as f64 - 0.5) * 2.0
        };

        // Add 20 epochs with decreasing noise (simulating convergence)
        for i in 0..20 {
            let sigma = 5.0 / (1.0 + i as f64 * 0.3); // σ starts at 5m, drops to ~0.8m
            let err = Vector3::new(noise(), noise(), noise()) * sigma;
            combiner.add_epoch(
                truth + err,
                Matrix3::identity() * (sigma * sigma),
            );
        }

        let avg = combiner.weighted_average().unwrap();
        let err = (avg - truth).norm();

        // With 20 epochs of weighted averaging, error should be well under 1m
        assert!(err < 0.5, "20-epoch weighted avg should be sub-0.5m: got {:.3}m", err);
    }

    #[test]
    fn test_median_position_resists_outliers() {
        let truth = Vector3::new(0.0, 0.0, 0.0);
        let mut combiner = MultiEpochCombiner::new();

        // Add 9 good epochs
        for _ in 0..9 {
            combiner.add_epoch(
                Vector3::new(0.0, 0.0, 0.0),
                Matrix3::identity() * 0.01,
            );
        }
        // Add 1 outlier
        combiner.add_epoch(
            Vector3::new(100.0, 100.0, 100.0),
            Matrix3::identity() * 0.01, // same weight — would corrupt weighted avg
        );

        let med = combiner.median_position().unwrap();
        assert_eq!(med, Vector3::zeros(),
            "Median should reject outlier: got ({:.1}, {:.1}, {:.1})", med.x, med.y, med.z);
    }
}
