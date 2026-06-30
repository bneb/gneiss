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

// ---------------------------------------------------------------------------
// Cross-Base AR Validator
// ---------------------------------------------------------------------------

/// Result of cross-base AR consensus check.
#[derive(Clone, Debug)]
pub struct CrossBaseResult {
    /// Whether the cross-base consensus check passed.
    pub accepted: bool,
    /// Indices of bases whose AR fixes agree (within threshold of each other).
    pub agreed_indices: Vec<usize>,
    /// Consensus position (weighted average of agreeing fixes), if any.
    pub consensus_position: Option<Vector3<f64>>,
    /// Maximum pairwise disagreement among all fix-producing bases (meters).
    pub max_disagreement_m: f64,
}

/// Cross-base AR validator for multi-base RTK.
///
/// When ≥2 bases produce AR fixes, this validator checks that their position
/// solutions agree within a configurable threshold. Disagreement indicates
/// that different bases likely picked different (wrong) NL integers due to
/// different code multipath at each base. In that case, all fixes are rejected
/// and the solution stays in float mode.
///
/// If bases agree, both have likely found correct integers — the
/// weighted-average fix is accepted.
pub struct CrossBaseArValidator {
    /// Maximum allowable 3D position difference between two bases' AR fixes (meters).
    /// Default 0.3m — about 1.5× NL wavelength, tight enough to catch wrong integer fixes.
    agreement_threshold_m: f64,
}

impl CrossBaseArValidator {
    /// Create a new validator with the given agreement threshold.
    pub fn new(agreement_threshold_m: f64) -> Self {
        Self {
            agreement_threshold_m,
        }
    }

    /// Validate AR fixes from multiple bases.
    ///
    /// `fixes` is a list of `(position, is_fixed, baseline_distance_m)` tuples,
    /// one per base. Returns the consensus result.
    ///
    /// Logic:
    /// - < 2 AR-fixed bases → accept whatever we have (no cross-validation possible)
    /// - ≥ 2 AR-fixed bases → check pairwise agreement
    ///   - If the two closest fixes agree within threshold → accept, return
    ///     consensus of all agreeing fixes
    ///   - If no pair agrees within threshold → reject all fixes
    pub fn validate(
        &self,
        fixes: &[(Vector3<f64>, bool, f64)],
    ) -> CrossBaseResult {
        // Collect indices of bases that produced AR fixes
        let fixed_indices: Vec<usize> = fixes
            .iter()
            .enumerate()
            .filter(|(_, (_, is_fixed, _))| *is_fixed)
            .map(|(i, _)| i)
            .collect();

        if fixed_indices.len() < 2 {
            // Not enough fixed bases to cross-validate.
            // Return the single fix (if any) as accepted.
            let consensus = if fixed_indices.len() == 1 {
                Some(fixes[fixed_indices[0]].0)
            } else {
                None
            };
            return CrossBaseResult {
                accepted: true,
                agreed_indices: fixed_indices,
                consensus_position: consensus,
                max_disagreement_m: 0.0,
            };
        }

        // Compute pairwise distances among fixed positions
        let n_fixed = fixed_indices.len();
        let mut max_disagreement = 0.0_f64;
        let mut min_distance = f64::MAX;
        let mut min_pair: (usize, usize) = (0, 0);

        for i in 0..n_fixed {
            for j in (i + 1)..n_fixed {
                let dist = (fixes[fixed_indices[i]].0 - fixes[fixed_indices[j]].0).norm();
                max_disagreement = max_disagreement.max(dist);
                if dist < min_distance {
                    min_distance = dist;
                    min_pair = (fixed_indices[i], fixed_indices[j]);
                }
            }
        }

        if min_distance > self.agreement_threshold_m {
            // Closest pair still disagrees → all fixes are suspect
            return CrossBaseResult {
                accepted: false,
                agreed_indices: Vec::new(),
                consensus_position: None,
                max_disagreement_m: max_disagreement,
            };
        }

        // At least the closest pair agrees. Find all fixes that agree with
        // either of the agreeing pair (transitive closure).
        let mut agreed: Vec<usize> = vec![min_pair.0, min_pair.1];
        for &idx in &fixed_indices {
            if agreed.contains(&idx) {
                continue;
            }
            // Check if this fix agrees with any already-agreed fix
            let agrees = agreed.iter().any(|&a| {
                (fixes[idx].0 - fixes[a].0).norm() <= self.agreement_threshold_m
            });
            if agrees {
                agreed.push(idx);
            }
        }

        // Weighted average of agreeing fixes (weight = 1/baseline_distance)
        let mut weighted_pos = Vector3::zeros();
        let mut total_weight = 0.0_f64;
        for &idx in &agreed {
            let baseline = fixes[idx].2.max(100.0); // min 100m
            let w = 1.0 / (baseline * baseline);
            weighted_pos += fixes[idx].0 * w;
            total_weight += w;
        }

        let consensus = if total_weight > 0.0 {
            Some(weighted_pos / total_weight)
        } else {
            Some(fixes[agreed[0]].0)
        };

        CrossBaseResult {
            accepted: true,
            agreed_indices: agreed,
            consensus_position: consensus,
            max_disagreement_m: max_disagreement,
        }
    }
}

impl Default for CrossBaseArValidator {
    fn default() -> Self {
        Self::new(0.3)
    }
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

    // -----------------------------------------------------------------------
    // CrossBaseArValidator tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_cross_base_less_than_two_fixes_accepted() {
        let v = CrossBaseArValidator::new(0.3);
        // Only one base with a fix → accepted (nothing to cross-validate)
        let fixes = vec![
            (Vector3::new(1.0, 0.0, 0.0), true, 5000.0),
            (Vector3::new(5.0, 0.0, 0.0), false, 8000.0), // float, not fixed
        ];
        let result = v.validate(&fixes);
        assert!(result.accepted);
        assert_eq!(result.agreed_indices, vec![0]);
        assert!(result.consensus_position.is_some());
    }

    #[test]
    fn test_cross_base_two_agreeing_fixes_accepted() {
        let v = CrossBaseArValidator::new(0.3);
        let fixes = vec![
            (Vector3::new(1.0, 0.0, 0.0), true, 5000.0),
            (Vector3::new(1.15, 0.0, 0.0), true, 8000.0), // 0.15m away → agrees
        ];
        let result = v.validate(&fixes);
        assert!(result.accepted);
        assert_eq!(result.agreed_indices.len(), 2);
        assert!(result.consensus_position.is_some());
        // Consensus should be between the two positions
        let cp = result.consensus_position.unwrap();
        assert!(cp.x > 1.0 && cp.x < 1.15);
    }

    #[test]
    fn test_cross_base_two_disagreeing_fixes_rejected() {
        let v = CrossBaseArValidator::new(0.3);
        let fixes = vec![
            (Vector3::new(1.0, 0.0, 0.0), true, 5000.0),
            (Vector3::new(2.0, 0.0, 0.0), true, 8000.0), // 1.0m away → disagrees
        ];
        let result = v.validate(&fixes);
        assert!(!result.accepted);
        assert!(result.agreed_indices.is_empty());
        assert!(result.consensus_position.is_none());
        assert!(result.max_disagreement_m > 0.9);
    }

    #[test]
    fn test_cross_base_three_two_agree_one_disagree() {
        let v = CrossBaseArValidator::new(0.3);
        let fixes = vec![
            (Vector3::new(1.0, 0.0, 0.0), true, 5000.0),
            (Vector3::new(1.1, 0.0, 0.0), true, 8000.0),  // agrees with #0
            (Vector3::new(5.0, 0.0, 0.0), true, 12000.0), // 4m away → disagrees
        ];
        let result = v.validate(&fixes);
        assert!(result.accepted);
        assert_eq!(result.agreed_indices.len(), 2);
        assert!(!result.agreed_indices.contains(&2)); // index 2 excluded
        assert!(result.max_disagreement_m > 3.9);
    }

    #[test]
    fn test_cross_base_no_fixes_accepted() {
        let v = CrossBaseArValidator::new(0.3);
        let fixes = vec![
            (Vector3::new(1.0, 0.0, 0.0), false, 5000.0),
            (Vector3::new(2.0, 0.0, 0.0), false, 8000.0),
        ];
        let result = v.validate(&fixes);
        assert!(result.accepted); // nothing to reject
        assert!(result.agreed_indices.is_empty());
        assert!(result.consensus_position.is_none());
    }

    #[test]
    fn test_cross_base_default_threshold() {
        let v = CrossBaseArValidator::default();
        // Default threshold is 0.3m
        let fixes = vec![
            (Vector3::new(0.0, 0.0, 0.0), true, 5000.0),
            (Vector3::new(0.25, 0.0, 0.0), true, 8000.0), // 0.25m < 0.3m
        ];
        let result = v.validate(&fixes);
        assert!(result.accepted);
    }

    #[test]
    fn test_cross_base_at_threshold_boundary() {
        // Exactly at threshold should be rejected (> not >=)
        let v = CrossBaseArValidator::new(0.3);
        let fixes = vec![
            (Vector3::new(0.0, 0.0, 0.0), true, 5000.0),
            (Vector3::new(0.3, 0.0, 0.0), true, 8000.0), // exactly 0.3m
        ];
        let result = v.validate(&fixes);
        // 0.3 is NOT > 0.3, so it should be accepted
        assert!(result.accepted);
    }

    #[test]
    fn test_cross_base_barely_over_threshold() {
        let v = CrossBaseArValidator::new(0.3);
        let fixes = vec![
            (Vector3::new(0.0, 0.0, 0.0), true, 5000.0),
            (Vector3::new(0.3000001, 0.0, 0.0), true, 8000.0),
        ];
        let result = v.validate(&fixes);
        assert!(!result.accepted);
    }
}
