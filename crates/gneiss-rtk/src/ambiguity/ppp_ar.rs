//! Integer Ambiguity Resolution for Precise Point Positioning (PPP-AR).
//!
//! Implements Decoupled Clock / Fractional Phase Bias (FCB/OSB) Wide-Lane & Narrow-Lane
//! integer ambiguity resolution:
//!
//! 1. **Wide-Lane (WL) AR**: Resolved per-satellite from Melbourne-Wübbena (MW) LC:
//!    \hat{N}_{WL} = \text{round}\left( \frac{\Phi_{MW} - B_{WL}}{\lambda_{WL}} \right)
//!
//! 2. **Narrow-Lane (NL) AR**: Resolved via LAMBDA over the Ionosphere-Free (IF) LC
//!    conditioned on fixed Wide-Lane integers:
//!    L_{IF} - \frac{c f_2}{f_1^2 - f_2^2} \lambda_{WL} N_{WL} = \frac{c}{f_1 + f_2} N_{NL}

use nalgebra::{DMatrix, DVector};
use crate::ambiguity::lambda::resolve_lambda;

/// Fractional cycle bias correction for a single satellite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SatellitePhaseBiases {
    /// Wide-lane phase bias (cycles).
    pub bias_wl_cycles: f64,
    /// Narrow-lane phase bias (cycles).
    pub bias_nl_cycles: f64,
}

/// Wide-lane ambiguity state for a tracked satellite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WideLaneCandidate {
    /// Satellite index / identifier.
    pub sat_idx: usize,
    /// Smoothed Melbourne-Wübbena measurement in cycles.
    pub mw_cycles: f64,
    /// Standard deviation of MW measurement in cycles.
    pub mw_std_cycles: f64,
    /// Satellite wide-lane phase bias in cycles.
    pub bias_wl_cycles: f64,
}

/// Result of wide-lane fixing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FixedWideLane {
    pub sat_idx: usize,
    pub n_wl: i32,
    pub fractional_residual: f64,
}

/// Result of PPP Integer Ambiguity Resolution.
#[derive(Debug, Clone, PartialEq)]
pub struct PppArResult {
    /// Fixed wide-lane integer ambiguities.
    pub fixed_wl: Vec<FixedWideLane>,
    /// Fixed narrow-lane integer ambiguities.
    pub fixed_nl: Vec<(usize, i32)>,
    /// Ratio test statistic from LAMBDA.
    pub lambda_ratio: f64,
    /// Whether full integer fix succeeded.
    pub is_fixed: bool,
}

/// PPP-AR Decoupled Clock / Fractional Phase Bias Solver.
pub struct PppArSolver;

impl PppArSolver {
    /// Resolve Wide-Lane integer ambiguities from smoothed Melbourne-Wübbena measurements.
    pub fn fix_wide_lane(
        candidates: &[WideLaneCandidate],
        max_fractional_error: f64,
        max_sigma: f64,
    ) -> Vec<FixedWideLane> {
        let mut fixed = Vec::new();

        for c in candidates {
            if c.mw_std_cycles > max_sigma {
                continue;
            }

            let corrected_mw = c.mw_cycles - c.bias_wl_cycles;
            let n_rounded = libm::round(corrected_mw) as i32;
            let frac_err = libm::fabs(corrected_mw - (n_rounded as f64));

            if frac_err <= max_fractional_error {
                fixed.push(FixedWideLane {
                    sat_idx: c.sat_idx,
                    n_wl: n_rounded,
                    fractional_residual: frac_err,
                });
            }
        }

        fixed
    }

    /// Resolve Narrow-Lane integer ambiguities using LAMBDA conditioned on fixed Wide-Lane integers.
    pub fn fix_narrow_lane(
        float_nl_ambiguities: &DVector<f64>,
        cov_nl: &DMatrix<f64>,
        fixed_wl: &[FixedWideLane],
        min_ratio: f64,
    ) -> Option<PppArResult> {
        let n_amb = float_nl_ambiguities.len();
        if n_amb < 4 || fixed_wl.len() < 4 {
            return None;
        }

        let lambda_res = resolve_lambda(float_nl_ambiguities, cov_nl).ok()?;

        if lambda_res.ratio < min_ratio {
            return Some(PppArResult {
                fixed_wl: fixed_wl.to_vec(),
                fixed_nl: Vec::new(),
                lambda_ratio: lambda_res.ratio,
                is_fixed: false,
            });
        }

        let fixed_nl = (0..n_amb)
            .map(|i| (fixed_wl[i].sat_idx, lambda_res.best_integers[i] as i32))
            .collect();

        Some(PppArResult {
            fixed_wl: fixed_wl.to_vec(),
            fixed_nl,
            lambda_ratio: lambda_res.ratio,
            is_fixed: true,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wide_lane_rounding_with_fractional_biases() {
        let candidates = vec![
            WideLaneCandidate {
                sat_idx: 1,
                mw_cycles: 12.23,
                mw_std_cycles: 0.08,
                bias_wl_cycles: 0.25, // 12.23 - 0.25 = 11.98 -> 12
            },
            WideLaneCandidate {
                sat_idx: 2,
                mw_cycles: -5.72,
                mw_std_cycles: 0.05,
                bias_wl_cycles: 0.30, // -5.72 - 0.30 = -6.02 -> -6
            },
            WideLaneCandidate {
                sat_idx: 3,
                mw_cycles: 8.50, // ambiguous
                mw_std_cycles: 0.10,
                bias_wl_cycles: 0.0,
            },
        ];

        let fixed = PppArSolver::fix_wide_lane(&candidates, 0.15, 0.20);
        assert_eq!(fixed.len(), 2);
        assert_eq!(fixed[0].sat_idx, 1);
        assert_eq!(fixed[0].n_wl, 12);
        assert_eq!(fixed[1].sat_idx, 2);
        assert_eq!(fixed[1].n_wl, -6);
    }
}
