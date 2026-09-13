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
        let is_fixed = lambda_res.ratio >= min_ratio;
        let fixed_nl = if is_fixed {
            (0..n_amb)
                .map(|i| (fixed_wl[i].sat_idx, lambda_res.best_integers[i] as i32))
                .collect()
        } else {
            Vec::new()
        };
        Some(PppArResult {
            fixed_wl: fixed_wl.to_vec(),
            fixed_nl,
            lambda_ratio: lambda_res.ratio,
            is_fixed,
        })
    }

    /// Resolve Between-Satellite Single-Difference (SD) Wide-Lane integer ambiguities.
    pub fn fix_sd_wide_lane(
        ref_sat_idx: usize,
        ref_mw_cycles: f64,
        ref_bias_wl_cycles: f64,
        candidates: &[WideLaneCandidate],
        max_fractional_error: f64,
    ) -> Vec<FixedWideLane> {
        let mut fixed = Vec::new();
        let ref_corrected_mw = ref_mw_cycles - ref_bias_wl_cycles;
        for c in candidates {
            if c.sat_idx == ref_sat_idx {
                continue;
            }
            let c_corrected_mw = c.mw_cycles - c.bias_wl_cycles;
            let sd_mw = c_corrected_mw - ref_corrected_mw;
            let n_rounded = libm::round(sd_mw) as i32;
            let frac_err = libm::fabs(sd_mw - (n_rounded as f64));
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

    /// Resolve single-differenced ambiguities using LAMBDA search.
    pub fn fix_single_diff_ambiguities(
        &self,
        float_ambiguities: &[f64],
        cov: &DMatrix<f64>,
        wavelengths: &[f64],
    ) -> Result<Vec<i32>, crate::estimators::eskf::types::EngineError> {
        let n = float_ambiguities.len();
        if n < 4 || cov.nrows() != n || cov.ncols() != n || wavelengths.len() != n {
            return Err(crate::estimators::eskf::types::EngineError::InvalidMeasurement(
                "Need at least 4 ambiguities with matching cov and wavelengths".into(),
            ));
        }

        let mut float_cycles = DVector::zeros(n);
        let mut cov_cycles = DMatrix::zeros(n, n);
        for i in 0..n {
            let wl_i = if wavelengths[i] > 0.0 { wavelengths[i] } else { 1.0 };
            float_cycles[i] = float_ambiguities[i] / wl_i;
            for j in 0..n {
                let wl_j = if wavelengths[j] > 0.0 { wavelengths[j] } else { 1.0 };
                cov_cycles[(i, j)] = cov[(i, j)] / (wl_i * wl_j);
            }
        }

        let lambda_res = resolve_lambda(&float_cycles, &cov_cycles)
            .map_err(|e| crate::estimators::eskf::types::EngineError::Internal(format!("LAMBDA failed: {:?}", e)))?;
        if lambda_res.ratio < 2.0 {
            return Err(crate::estimators::eskf::types::EngineError::Internal("Ratio test failed".into()));
        }
        Ok(lambda_res.best_integers.iter().map(|&v| v as i32).collect())
    }

    /// Form single-differenced float vector and covariance against reference satellite at index 0.
    pub fn form_single_differences(
        float_amb: &[f64],
        cov: &DMatrix<f64>,
    ) -> (DVector<f64>, DMatrix<f64>, DMatrix<f64>) {
        let m = float_amb.len();
        let n = m - 1;
        let mut d = DMatrix::zeros(n, m);
        for k in 0..n {
            d[(k, 0)] = -1.0;
            d[(k, k + 1)] = 1.0;
        }
        let float_vec = DVector::from_row_slice(float_amb);
        let sd_float = &d * &float_vec;
        let sd_cov = &d * cov * d.transpose();
        (sd_float, sd_cov, d)
    }

    /// Back-substitute fixed single differences into un-differenced ambiguity states.
    pub fn backsubstitute_sd_fix(
        float_amb: &[f64],
        cov: &DMatrix<f64>,
        d: &DMatrix<f64>,
        sd_cov: &DMatrix<f64>,
        sd_float: &DVector<f64>,
        fixed_sd: &DVector<f64>,
    ) -> Option<DVector<f64>> {
        let chol = sd_cov.clone().cholesky()?;
        let delta = fixed_sd - sd_float;
        let y = chol.solve(&delta);
        let float_vec = DVector::from_row_slice(float_amb);
        let q_dt = cov * d.transpose();
        Some(float_vec + q_dt * y)
    }

    fn form_nl_system(
        sd_float: &DVector<f64>,
        sd_cov: &DMatrix<f64>,
        f1: f64,
        f2: f64,
        fixed_n_wl: &[i32],
    ) -> (DVector<f64>, DMatrix<f64>, f64, f64) {
        let n = sd_float.len();
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let lambda_nl = c / (f1 + f2);
        let alpha = c * f2 / (f1 * f1 - f2 * f2);
        let mut float_nl = DVector::zeros(n);
        let mut cov_nl = DMatrix::zeros(n, n);
        for k in 0..n {
            float_nl[k] = (sd_float[k] - alpha * (fixed_n_wl[k] as f64)) / lambda_nl;
        }
        for i in 0..n {
            for j in 0..n {
                cov_nl[(i, j)] = sd_cov[(i, j)] / (lambda_nl * lambda_nl);
            }
        }
        (float_nl, cov_nl, lambda_nl, alpha)
    }

    /// Execute single-differenced Narrow-Lane LAMBDA AR conditioned on fixed Wide-Lane integers.
    pub fn resolve_sd_with_fixed_wl(
        float_amb: &[f64],
        cov: &DMatrix<f64>,
        f1: f64,
        f2: f64,
        fixed_n_wl: &[i32],
        min_ratio: f64,
    ) -> Option<(DVector<f64>, f64)> {
        if float_amb.len() < 4 || fixed_n_wl.len() != float_amb.len() - 1 {
            return None;
        }
        let (sd_float, sd_cov, d) = Self::form_single_differences(float_amb, cov);
        let (float_nl, cov_nl, lambda_nl, alpha) =
            Self::form_nl_system(&sd_float, &sd_cov, f1, f2, fixed_n_wl);
        let lambda_res = resolve_lambda(&float_nl, &cov_nl).ok()?;
        if lambda_res.ratio < min_ratio {
            return None;
        }
        let best_ints = lambda_res.best_integers.as_slice();
        let mut fixed_sd_m = DVector::zeros(best_ints.len());
        for k in 0..best_ints.len() {
            fixed_sd_m[k] = lambda_nl * best_ints[k] + alpha * (fixed_n_wl[k] as f64);
        }
        let fixed_undiff = Self::backsubstitute_sd_fix(
            float_amb, cov, &d, &sd_cov, &sd_float, &fixed_sd_m,
        )?;
        Some((fixed_undiff, lambda_res.ratio))
    }

    /// Execute single-differenced Wide-Lane / Narrow-Lane cascade LAMBDA AR.
    pub fn resolve_sd_cascade(
        float_amb: &[f64],
        cov: &DMatrix<f64>,
        f1: f64,
        f2: f64,
        wl_biases: &[f64],
        min_ratio: f64,
    ) -> Option<(DVector<f64>, f64)> {
        let n = if float_amb.is_empty() { 0 } else { float_amb.len() - 1 };
        let fixed_wl: Vec<i32> = wl_biases.iter().take(n).map(|&b| libm::round(b) as i32).collect();
        Self::resolve_sd_with_fixed_wl(float_amb, cov, f1, f2, &fixed_wl, min_ratio)
    }
}

/// Multi-epoch Melbourne-Wübbena phase-code combination accumulator for PPP-AR.
#[derive(Debug, Clone, Default)]
pub struct PppMwTracker {
    arcs: std::collections::HashMap<(u8, u16), (f64, usize, u32)>,
}

impl PppMwTracker {
    pub fn new() -> Self {
        Self { arcs: std::collections::HashMap::new() }
    }

    pub fn update(&mut self, sat_key: (u8, u16), mw_cycles: f64, epoch: u32, slip: bool) {
        let entry = self.arcs.entry(sat_key).or_insert((mw_cycles, 0, epoch));
        if slip || epoch > entry.2 + 2 {
            *entry = (mw_cycles, 1, epoch);
        } else {
            entry.1 += 1;
            entry.0 += (mw_cycles - entry.0) / (entry.1 as f64);
            entry.2 = epoch;
        }
    }

    pub fn get_smoothed(&self, sat_key: (u8, u16)) -> Option<(f64, usize)> {
        let &(mean, count, _) = self.arcs.get(&sat_key)?;
        Some((mean, count))
    }

    pub fn fix_sd_wide_lane_subset(
        &self,
        ref_sat: (u8, u16),
        candidates: &[(u8, u16)],
        min_count: usize,
        max_frac_err: f64,
    ) -> Vec<((u8, u16), i32)> {
        let (ref_mw, _ref_count) = match self.get_smoothed(ref_sat) {
            Some(v) if v.1 >= min_count => v,
            _ => return Vec::new(),
        };
        let mut fixed = Vec::new();
        for &sat in candidates {
            if let Some((mw, count)) = self.get_smoothed(sat) {
                if count >= min_count {
                    let sd_mw = mw - ref_mw;
                    let rounded = libm::round(sd_mw) as i32;
                    let frac = (sd_mw - rounded as f64).abs();
                    if std::env::var("PPP_AR_DEBUG").is_ok() {
                        println!("DEBUG MW: sat={:?} mw={:.3} ref={:?} ref_mw={:.3} sd={:.3} frac={:.3}", sat, mw, ref_sat, ref_mw, sd_mw, frac);
                    }
                    if frac <= max_frac_err {
                        fixed.push((sat, rounded));
                    }
                }
            }
        }
        fixed
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

    #[test]
    fn test_between_satellite_single_difference_wide_lane_fixing() {
        let candidates = vec![
            WideLaneCandidate {
                sat_idx: 10,
                mw_cycles: 105.42,
                mw_std_cycles: 0.05,
                bias_wl_cycles: 0.40, // corrected = 105.02
            },
            WideLaneCandidate {
                sat_idx: 14,
                mw_cycles: 82.15,
                mw_std_cycles: 0.06,
                bias_wl_cycles: 0.13, // corrected = 82.02
            },
        ];

        // Reference satellite: corrected = 50.00
        let sd_fixed = PppArSolver::fix_sd_wide_lane(1, 50.35, 0.35, &candidates, 0.10);
        assert_eq!(sd_fixed.len(), 2);
        assert_eq!(sd_fixed[0].sat_idx, 10);
        assert_eq!(sd_fixed[0].n_wl, 55); // 105.02 - 50.00 = 55.02 -> 55
        assert_eq!(sd_fixed[1].sat_idx, 14);
        assert_eq!(sd_fixed[1].n_wl, 32); // 82.02 - 50.00 = 32.02 -> 32
    }

    #[test]
    fn test_fix_single_diff_ambiguities() {
        let solver = PppArSolver;
        let float_amb = vec![1.01, 2.02, -3.01, 5.02];
        let mut cov = DMatrix::zeros(4, 4);
        for i in 0..4 { cov[(i, i)] = 0.0001; }
        let wavelengths = vec![1.0, 1.0, 1.0, 1.0];
        let fixed = solver.fix_single_diff_ambiguities(&float_amb, &cov, &wavelengths).unwrap();
        assert_eq!(fixed, vec![1, 2, -3, 5]);
    }

    #[test]
    fn test_sd_form_and_backsubstitute() {
        let float_amb = vec![10.0, 12.01, 15.02, 7.99, 18.01];
        let mut cov = DMatrix::zeros(5, 5);
        for i in 0..5 { cov[(i, i)] = 0.01; }
        let (sd_float, sd_cov, d) = PppArSolver::form_single_differences(&float_amb, &cov);
        assert_eq!(sd_float.len(), 4);
        assert_eq!(sd_cov.nrows(), 4);

        let fixed_sd = DVector::from_row_slice(&[2.0, 5.0, -2.0, 8.0]);
        let undiff = PppArSolver::backsubstitute_sd_fix(&float_amb, &cov, &d, &sd_cov, &sd_float, &fixed_sd)
            .expect("backsubstitution must invert positive definite cov");
        let resulting_sd = &d * &undiff;
        for k in 0..4 {
            assert!((resulting_sd[k] - fixed_sd[k]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_ppp_mw_tracker_and_fixed_wl() {
        let mut tracker = PppMwTracker::new();
        for ep in 0..15 {
            tracker.update((0, 1), 50.0 + 0.05 * (ep as f64 % 2.0), ep, false);
            tracker.update((0, 2), 70.0 - 0.04 * (ep as f64 % 2.0), ep, false);
            tracker.update((0, 3), 85.0 + 0.02 * (ep as f64 % 2.0), ep, false);
        }
        let fixed = tracker.fix_sd_wide_lane_subset((0, 1), &[(0, 2), (0, 3)], 10, 0.35);
        assert_eq!(fixed.len(), 2);
        assert_eq!(fixed[0].1, 20); // 70 - 50 = 20
        assert_eq!(fixed[1].1, 35); // 85 - 50 = 35
    }

    #[test]
    fn test_resolve_sd_with_fixed_wl() {
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let lambda_nl = c / (f1 + f2);
        let alpha = c * f2 / (f1 * f1 - f2 * f2);

        let wl_integers = [3, -2, 5, 1];
        let nl_integers = [12, -8, 15, 7];
        let n = wl_integers.len();
        let m = n + 1;

        let mut expected_sd = DVector::zeros(n);
        for k in 0..n {
            expected_sd[k] = lambda_nl * (nl_integers[k] as f64) + alpha * (wl_integers[k] as f64);
        }

        let mut float_amb = vec![10.0; m];
        for k in 0..n {
            float_amb[k + 1] = float_amb[0] + expected_sd[k] + 0.001;
        }

        let mut cov = DMatrix::zeros(m, m);
        for i in 0..m { cov[(i, i)] = 0.0001; }

        let (fixed_undiff, ratio) = PppArSolver::resolve_sd_with_fixed_wl(
            &float_amb, &cov, f1, f2, &wl_integers, 2.0,
        ).expect("LAMBDA AR must succeed with clean ambiguities");

        assert!(ratio >= 2.0);
        for k in 0..n {
            let actual_sd = fixed_undiff[k + 1] - fixed_undiff[0];
            assert!((actual_sd - expected_sd[k]).abs() < 1e-4);
        }
    }
}

