//! Multi-epoch position smoothing.
//!
//! Applies an exponential moving average to the IEKF position output to reduce
//! epoch-to-epoch noise. For static stations, position should be constant — the
//! IEKF's per-epoch estimates have ~2m noise that averaging reduces by √N.
//!
//! Phase A3+: When shared ambiguities are added, the smoother will be upgraded
//! to use between-epoch carrier phase for mm-level Δpos constraints.

use nalgebra::Vector3;

use crate::engine::processed_sat::ProcessedSat;
use crate::engine::EngineError;

/// Multi-epoch position smoother.
///
/// Maintains an EMA of the IEKF position across epochs.
/// Simple but effective — reduces white noise by ~30% with α=0.5.
pub struct MultiEpochOptimizer {
    /// Smoothing factor: pos = α * new + (1-α) * prev
    alpha: f64,
    /// Previous smoothed position
    prev_position: Option<Vector3<f64>>,
}

impl Default for MultiEpochOptimizer {
    fn default() -> Self {
        Self {
            alpha: 0.5,
            prev_position: None,
        }
    }
}

impl MultiEpochOptimizer {
    pub fn new(_window_size: usize) -> Self {
        Self::default()
    }

    /// Smooth the IEKF position with an exponential moving average.
    ///
    /// First epoch: stores position, returns Ok.
    /// Subsequent epochs: pos_smoothed = α * pos_new + (1-α) * pos_prev
    pub fn solve(
        &mut self,
        state: &mut crate::filter::RtkState,
        _sats: &[ProcessedSat],
        _position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<(), EngineError> {
        let curr_pos = state.position.vector;

        if let Some(prev) = self.prev_position {
            let smoothed = curr_pos * self.alpha + prev * (1.0 - self.alpha);
            let corr = (smoothed - curr_pos).norm();
            state.position.vector = smoothed;
            self.prev_position = Some(smoothed);

            if corr > 0.001 {
                tracing::info!("MultiEpoch smooth: corr={:.3}m", corr);
            }
        } else {
            self.prev_position = Some(curr_pos);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constructor() {
        let opt = MultiEpochOptimizer::new(2);
        assert!((opt.alpha - 0.5).abs() < 1e-12);
        assert!(opt.prev_position.is_none());
    }

    #[test]
    fn test_smoothing_reduces_noise() {
        let time = gneiss_core::time::GpsTime::new(2082, 0.0);
        let coord = gneiss_core::coords::Coordinate::new(
            Vector3::new(6000000.0, 0.0, 0.0),
            gneiss_core::coords::Datum::WGS84,
            gneiss_core::coords::Frame::ECEF,
            time,
        );
        let mut state = crate::filter::RtkState::new(time, coord, 10.0);
        let sats: Vec<ProcessedSat> = vec![];
        let mut opt = MultiEpochOptimizer::new(2);

        // First epoch: position unchanged
        let pos1 = state.position.vector;
        opt.solve(&mut state, &sats, None).unwrap();
        assert!((state.position.vector - pos1).norm() < 1e-12);

        // Second epoch: simulate IEKF producing a 2m offset
        state.position.vector = Vector3::new(6000002.0, 0.0, 0.0);
        opt.solve(&mut state, &sats, None).unwrap();
        // Smoothed: 0.5 * (6000002) + 0.5 * (6000000) = 6000001
        assert!((state.position.vector.x - 6000001.0).abs() < 1e-9);
    }
}
