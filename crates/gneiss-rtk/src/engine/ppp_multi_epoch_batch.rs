//! Multi-epoch batch position solver for static PPP.
//!
//! For static receivers, carrier phase measurements from AR-fixed satellites
//! are mm-precision unbiased range measurements. Stacking these across epochs
//! with diverse satellite geometry yields significantly more precise position
//! than the single-epoch IEKF, which is limited by pseudorange (~5m) anchoring.
//!
//! This solver:
//! 1. Accumulates "pseudo-range" measurements derived from IF carrier phase
//!    (cp_if - N_if_fixed - tropo - clock) across epochs
//! 2. Jointly solves for a single static position using weighted least squares
//! 3. Estimates per-satellite biases to absorb residual errors from the IEKF
//!    (fractional-cycle IF ambiguity, ZWD, satellite orbit/clock, etc.)
//! 4. Uses Schur complement to efficiently eliminate biases from the 3D position solve
//! 5. Applies the refined position to the output state

use nalgebra::{Matrix3, Vector3};

use gneiss_core::sat::SatelliteId;

/// A single CP-derived range measurement for batch position solving.
#[derive(Clone, Debug)]
struct CpRangeMeasurement {
    /// Satellite identifier for per-satellite bias estimation
    sat_id: SatelliteId,
    /// Satellite position at time of measurement (ECEF, meters)
    sat_pos: Vector3<f64>,
    /// Measured range derived from IF carrier phase:
    ///   range = cp_if - N_if_fixed - tropo_dry - map_wet * zwd - clock
    measured_range: f64,
    /// Measurement weight (1/variance). CP noise is ~1cm, so weight ≈ 10000.
    weight: f64,
}

/// Batch position solver for static PPP with per-satellite bias estimation.
///
/// Accumulates CP-derived range measurements across epochs and solves for
/// the static receiver position using Gauss-Newton least squares.  Per-satellite
/// bias terms absorb residual IEKF errors (fractional-cycle IF ambiguity, ZWD,
/// satellite orbit/clock) and are eliminated via Schur complement for an
/// unbiased 3D position estimate.
pub struct StaticPositionBatchSolver {
    measurements: Vec<CpRangeMeasurement>,
    /// Current best position estimate (updated after each solve)
    position: Option<Vector3<f64>>,
    /// Position covariance from the last successful solve
    position_cov: Option<Matrix3<f64>>,
    /// Minimum number of measurements before attempting a solve
    min_measurements: usize,
    /// Maximum number of stored measurements (oldest evicted first)
    max_measurements: usize,
    /// L2 regularization strength for per-satellite biases (1/m²).
    /// Larger values = stronger shrinkage toward zero.
    /// Default 400 corresponds to σ_bias ≈ 5cm prior.
    bias_lambda: f64,
}

impl StaticPositionBatchSolver {
    pub fn new() -> Self {
        Self {
            measurements: Vec::new(),
            position: None,
            position_cov: None,
            min_measurements: 6,
            max_measurements: 1500, // ~150 epochs × 10 sats
            bias_lambda: 400.0,   // σ_bias ≈ 5cm prior
        }
    }

    /// Add a CP-derived range measurement for one AR-fixed satellite.
    ///
    /// # Arguments
    /// * `sat_id` - Satellite identifier (for per-satellite bias estimation)
    /// * `sat_pos` - Satellite ECEF position at measurement time
    /// * `cp_if` - Ionosphere-free carrier phase measurement (meters)
    /// * `n_if_fixed` - AR-fixed IF ambiguity (meters)
    /// * `tropo_dry` - Dry troposphere delay (meters)
    /// * `map_wet` - Wet troposphere mapping function
    /// * `zwd` - Zenith wet delay estimate (meters)
    /// * `clock` - Receiver clock bias estimate (meters)
    /// * `elevation` - Satellite elevation angle (radians)
    pub fn add_measurement(
        &mut self,
        sat_id: SatelliteId,
        sat_pos: Vector3<f64>,
        cp_if: f64,
        n_if_fixed: f64,
        tropo_dry: f64,
        map_wet: f64,
        zwd: f64,
        clock: f64,
        elevation: f64,
    ) {
        let measured_range = cp_if - n_if_fixed - tropo_dry - map_wet * zwd - clock;
        // Weight: σ=1cm for CP, elevated by sin(elevation) for low-elevation noise
        let sin_el = elevation.sin().max(0.1);
        let variance = 0.0001 / (sin_el * sin_el); // σ=1cm at zenith
        let weight = 1.0 / variance.max(1e-9);

        self.measurements.push(CpRangeMeasurement {
            sat_id,
            sat_pos,
            measured_range,
            weight,
        });

        // Evict oldest if over capacity
        while self.measurements.len() > self.max_measurements {
            self.measurements.remove(0);
        }
    }

    /// Solve for static position using all accumulated measurements.
    ///
    /// Uses Gauss-Newton nonlinear least squares with per-satellite bias
    /// estimation.  Biases are eliminated via Schur complement, yielding
    /// a 3D position estimate that is robust to systematic errors in the
    /// IF ambiguity / ZWD / satellite orbit inputs.
    ///
    /// Returns `(refined_position, position_covariance_3x3)` on success.
    pub fn solve(&mut self, initial_pos: Vector3<f64>) -> Option<(Vector3<f64>, Matrix3<f64>)> {
        if self.measurements.len() < self.min_measurements {
            return None;
        }

        let mut pos = self.position.unwrap_or(initial_pos);
        let max_iter = 5;
        let convergence_tol = 0.001; // 1mm

        // --- Pre-scan: map satellite IDs to contiguous indices ---------------
        let mut sat_to_idx: std::collections::HashMap<SatelliteId, usize> =
            std::collections::HashMap::new();
        for m in &self.measurements {
            let len = sat_to_idx.len();
            sat_to_idx.entry(m.sat_id).or_insert(len);
        }
        let n_sats = sat_to_idx.len();
        // Per-satellite accumulators for the Schur complement
        let mut sat_sum_w: Vec<f64> = vec![0.0f64; n_sats];
        let mut sat_sum_we: Vec<nalgebra::Vector3<f64>> = vec![nalgebra::Vector3::<f64>::zeros(); n_sats];
        let mut sat_sum_wr: Vec<f64> = vec![0.0f64; n_sats];

        for _iter in 0..max_iter {
            // Normal equation blocks
            let mut h_pp = Matrix3::zeros(); // position block (3×3)
            let mut g_p = Vector3::zeros(); // position RHS (3×1)

            // Reset per-satellite accumulators each iteration (residual changes)
            for i in 0..n_sats {
                sat_sum_w[i] = 0.0;
                sat_sum_we[i] = Vector3::zeros();
                sat_sum_wr[i] = 0.0;
            }

            for m in &self.measurements {
                let diff = pos - m.sat_pos;
                let dist = diff.norm();
                if dist < 1.0 {
                    continue;
                }
                // Line-of-sight unit vector (receiver → satellite)
                let e = diff / dist;
                let predicted = dist;
                let residual = m.measured_range - predicted;
                let w = m.weight;
                let si = sat_to_idx[&m.sat_id];

                // Accumulate position normal equations
                for i in 0..3 {
                    g_p[i] += w * e[i] * residual;
                    for j in 0..3 {
                        h_pp[(i, j)] += w * e[i] * e[j];
                    }
                }

                // Accumulate per-satellite statistics for Schur complement
                sat_sum_w[si] += w;
                for i in 0..3 {
                    sat_sum_we[si][i] += w * e[i];
                }
                sat_sum_wr[si] += w * residual;
            }

            // --- Schur complement: eliminate per-satellite biases -----------
            // The full system is:
            //   [H_pp  H_pb] [dx] = [g_p]
            //   [H_bp  H_bb] [db]   [g_b]
            //
            // H_bb is diagonal with H_bb[i,i] = sat_sum_w[i] + lambda (reg.)
            // H_pb[:,i] = sat_sum_we[i]
            // g_b[i]    = sat_sum_wr[i]
            //
            // Reduced normal equation (position only):
            //   H_red = H_pp - H_pb * (H_bb+λI)^{-1} * H_bp
            //   g_red = g_p  - H_pb * (H_bb+λI)^{-1} * g_b
            let mut h_red = h_pp;
            let mut g_red = g_p;

            for i in 0..n_sats {
                let denom = sat_sum_w[i] + self.bias_lambda;
                if denom < 1e-12 {
                    continue;
                }
                let inv = 1.0 / denom;
                let we = sat_sum_we[i];
                let wr = sat_sum_wr[i];

                // H_red -= inv * we * we^T
                for r in 0..3 {
                    for c in 0..3 {
                        h_red[(r, c)] -= inv * we[r] * we[c];
                    }
                }
                // g_red -= inv * we * wr
                for r in 0..3 {
                    g_red[r] -= inv * we[r] * wr;
                }
            }

            // --- Solve reduced 3×3 system ----------------------------------
            if let Some(h_inv) = try_invert_3x3(&h_red) {
                let dx = h_inv * g_red;
                pos += dx;
                if dx.norm() < convergence_tol {
                    self.position_cov = Some(h_inv);
                    break;
                }
            } else {
                return None;
            }
        }

        self.position = Some(pos);
        let cov = self.position_cov.unwrap_or_else(|| Matrix3::identity() * 1000.0);
        Some((pos, cov))
    }

    /// Number of accumulated measurements.
    pub fn num_measurements(&self) -> usize {
        self.measurements.len()
    }

    /// Clear all accumulated measurements (e.g., after cycle slip or mode change).
    pub fn clear(&mut self) {
        self.measurements.clear();
        self.position = None;
        self.position_cov = None;
    }
}

impl Default for StaticPositionBatchSolver {
    fn default() -> Self {
        Self::new()
    }
}

/// Try to invert a 3×3 matrix. Returns None if singular.
fn try_invert_3x3(m: &Matrix3<f64>) -> Option<Matrix3<f64>> {
    let det = m[(0, 0)] * (m[(1, 1)] * m[(2, 2)] - m[(1, 2)] * m[(2, 1)])
        - m[(0, 1)] * (m[(1, 0)] * m[(2, 2)] - m[(1, 2)] * m[(2, 0)])
        + m[(0, 2)] * (m[(1, 0)] * m[(2, 1)] - m[(1, 1)] * m[(2, 0)]);

    if det.abs() < 1e-15 {
        return None;
    }

    let inv_det = 1.0 / det;
    let mut inv = Matrix3::zeros();

    inv[(0, 0)] = (m[(1, 1)] * m[(2, 2)] - m[(1, 2)] * m[(2, 1)]) * inv_det;
    inv[(0, 1)] = (m[(0, 2)] * m[(2, 1)] - m[(0, 1)] * m[(2, 2)]) * inv_det;
    inv[(0, 2)] = (m[(0, 1)] * m[(1, 2)] - m[(0, 2)] * m[(1, 1)]) * inv_det;
    inv[(1, 0)] = (m[(1, 2)] * m[(2, 0)] - m[(1, 0)] * m[(2, 2)]) * inv_det;
    inv[(1, 1)] = (m[(0, 0)] * m[(2, 2)] - m[(0, 2)] * m[(2, 0)]) * inv_det;
    inv[(1, 2)] = (m[(0, 2)] * m[(1, 0)] - m[(0, 0)] * m[(1, 2)]) * inv_det;
    inv[(2, 0)] = (m[(1, 0)] * m[(2, 1)] - m[(1, 1)] * m[(2, 0)]) * inv_det;
    inv[(2, 1)] = (m[(0, 1)] * m[(2, 0)] - m[(0, 0)] * m[(2, 1)]) * inv_det;
    inv[(2, 2)] = (m[(0, 0)] * m[(1, 1)] - m[(0, 1)] * m[(1, 0)]) * inv_det;

    Some(inv)
}

// Keep the old MultiEpochCombiner for backward compatibility but deprecated.
/// Multi-epoch position combiner (DEPRECATED: use StaticPositionBatchSolver).
///
/// Simple weighted average of per-epoch positions. This does not account for
/// satellite geometry and is inferior to the batch solver for static receivers.
pub struct MultiEpochCombiner {
    positions: Vec<WeightedPosition>,
}

/// Weighted position estimate with covariance.
#[derive(Clone, Debug)]
pub struct WeightedPosition {
    pub pos: Vector3<f64>,
    pub cov: Matrix3<f64>,
    pub weight: f64,
}

impl MultiEpochCombiner {
    pub fn new() -> Self {
        Self { positions: Vec::new() }
    }

    pub fn add_epoch(&mut self, pos: Vector3<f64>, cov_3x3: Matrix3<f64>) {
        let weight = 1.0 / (cov_3x3[(0,0)] + cov_3x3[(1,1)] + cov_3x3[(2,2)]).max(1e-6);
        self.positions.push(WeightedPosition { pos, cov: cov_3x3, weight });
    }

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
    use gneiss_core::sat::Constellation;

    #[test]
    fn test_batch_solver_empty() {
        let mut solver = StaticPositionBatchSolver::new();
        assert_eq!(solver.solve(Vector3::new(1.0, 0.0, 0.0)), None);
    }

    #[test]
    fn test_batch_solver_single_sat_converges() {
        let mut solver = StaticPositionBatchSolver::new();
        let truth = Vector3::new(1_000_000.0, 2_000_000.0, 3_000_000.0);
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let sat3 = SatelliteId { constellation: Constellation::Gps, prn: 3 };

        // Add measurements from satellites at different positions
        // Satellite 1: directly overhead (z-direction)
        for _ in 0..5 {
            let sat = truth + Vector3::new(0.0, 0.0, 20_000_000.0);
            let dist = (sat - truth).norm();
            solver.add_measurement(sat1, sat, dist, 0.0, 0.0, 0.0, 0.0, 0.0, std::f64::consts::FRAC_PI_2);
        }
        // Satellite 2: east direction
        for _ in 0..5 {
            let sat = truth + Vector3::new(20_000_000.0, 0.0, 0.0);
            let dist = (sat - truth).norm();
            solver.add_measurement(sat2, sat, dist, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5);
        }
        // Satellite 3: north direction
        for _ in 0..5 {
            let sat = truth + Vector3::new(0.0, 20_000_000.0, 0.0);
            let dist = (sat - truth).norm();
            solver.add_measurement(sat3, sat, dist, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5);
        }

        let result = solver.solve(Vector3::new(1_000_100.0, 2_000_000.0, 3_000_000.0)); // 100m off in x
        assert!(result.is_some());
        let (refined, _cov) = result.unwrap();
        let err = (refined - truth).norm();
        assert!(err < 0.01, "Batch solver should converge to mm precision, got {:.3}m", err);
    }

    #[test]
    fn test_batch_solver_noisy_measurements() {
        let mut solver = StaticPositionBatchSolver::new();
        let truth = Vector3::new(1_000_000.0, 2_000_000.0, 3_000_000.0);
        let sat_ids: Vec<SatelliteId> = (0..7)
            .map(|i| SatelliteId { constellation: Constellation::Gps, prn: i + 1 })
            .collect();

        // Add many measurements with noise in different geometries
        for i in 0..50 {
            let az = (i as f64) * std::f64::consts::TAU / 50.0_f64;
            let el: f64 = 0.5; // 30 degrees
            let range = 20_000_000.0;
            let sat = truth + Vector3::new(
                range * el.cos() * az.cos(),
                range * el.cos() * az.sin(),
                range * el.sin(),
            );
            let dist = (sat - truth).norm();
            // Add ±1cm noise
            let noise = ((i % 7) as f64 - 3.0) * 0.005;
            let sat_id = sat_ids[(i % 7)];
            solver.add_measurement(sat_id, sat, dist + noise, 0.0, 0.0, 0.0, 0.0, 0.0, el);
        }

        // Start from 5m away
        let result = solver.solve(truth + Vector3::new(3.0, 4.0, 0.0));
        assert!(result.is_some());
        let (refined, _cov) = result.unwrap();
        let err = (refined - truth).norm();
        // With 50 measurements at 1cm noise each, position should be sub-cm
        assert!(err < 0.02, "Batch solver with 50 measurements should be sub-2cm, got {:.3}m", err);
    }

    #[test]
    fn test_weighted_average_two_epochs() {
        let mut combiner = MultiEpochCombiner::new();
        combiner.add_epoch(
            Vector3::new(1.0, 0.0, 0.0),
            Matrix3::identity() * 0.01,
        );
        combiner.add_epoch(
            Vector3::new(2.0, 0.0, 0.0),
            Matrix3::identity() * 1.0,
        );
        let avg = combiner.weighted_average().unwrap();
        assert!((avg.x - 1.0).abs() < 0.15,
            "Weighted avg should favor precise epoch: got {:.3}", avg.x);
    }

    #[test]
    fn test_invert_3x3_identity() {
        let m = Matrix3::identity();
        let inv = try_invert_3x3(&m).unwrap();
        assert!((inv - m).norm() < 1e-10);
    }
}
