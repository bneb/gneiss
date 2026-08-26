//! Sidereal phase primitives: conversion, bin folding, and the uniform-null
//! structure diagnostic. Pure `(tow, value)` sample math with no trajectory
//! or truth dependency — the production residual-stacking path lives here.

/// Mean sidereal day length in seconds (GPS repeat ground track).
pub const SIDEREAL_PERIOD_S: f64 = 86_164.0;
/// Default phase-bin count: 86164 s / 240 ~= 359 s (~6 min) per bin.
pub const DEFAULT_BINS: usize = 240;
/// Bin structure is claimed when the uniform-null p-value drops below this.
pub const STRUCTURE_P_THRESHOLD: f64 = 0.01;

/// One scalar observation tied to a time-of-week.
#[derive(Debug, Clone, Copy)]
pub struct Sample {
    pub tow_s: f64,
    pub value: f64,
}

/// Fractional sidereal phase in `[0, 1)`; wraps at every sidereal day.
///
/// NOT invariant under GPS week rollover: a week (604800 s) is 7 solar
/// days but ~7.018 sidereal days, so the phase slips ~27.5 min across a
/// week boundary (verified by test).
#[must_use]
pub fn sidereal_phase(tow_s: f64) -> f64 {
    tow_s.rem_euclid(SIDEREAL_PERIOD_S) / SIDEREAL_PERIOD_S
}

/// Bin index for a phase; out-of-domain input folds into bin 0, exact 1.0
/// clamps into the last bin.
#[must_use]
pub fn phase_bin(phase: f64, n_bins: usize) -> usize {
    if n_bins == 0 || !phase.is_finite() || phase < 0.0 {
        return 0;
    }
    ((phase * n_bins as f64).floor() as usize).min(n_bins - 1)
}

/// Running-mean/std accumulator (Welford) for one phase bin.
#[derive(Debug, Clone, Default)]
pub struct BinStat {
    pub count: u32,
    pub mean: f64,
    m2: f64,
}

impl BinStat {
    #[must_use]
    pub fn std(&self) -> f64 {
        if self.count > 1 {
            (self.m2 / (self.count - 1) as f64).sqrt()
        } else {
            0.0
        }
    }

    fn push(&mut self, v: f64) {
        self.count += 1;
        let d = v - self.mean;
        self.mean += d / self.count as f64;
        self.m2 += d * (v - self.mean);
    }
}

/// Fold samples onto sidereal phase bins. Non-finite input is skipped.
#[must_use]
pub fn fold(samples: &[Sample], n_bins: usize) -> Vec<BinStat> {
    let mut bins = vec![BinStat::default(); n_bins];
    for s in samples {
        if !s.tow_s.is_finite() || !s.value.is_finite() {
            continue;
        }
        let idx = phase_bin(sidereal_phase(s.tow_s), n_bins);
        bins[idx].push(s.value);
    }
    bins
}

/// One-way ANOVA-style uniformity verdict for folded bins: `chi2 =
/// SS_between / s2_within` is approximately chi-square with `dof = k-1`
/// under the no-structure null (exact for Gaussian noise up to the pooled
/// variance estimate; with thousands of epochs the approximation error is
/// negligible).
#[derive(Debug, Clone)]
pub struct Diagnostics {
    pub chi2: f64,
    pub dof: usize,
    pub p_value: f64,
}

impl Diagnostics {
    #[must_use]
    pub fn is_structured(&self) -> bool {
        self.dof > 0 && self.p_value.is_finite() && self.p_value < STRUCTURE_P_THRESHOLD
    }
}

#[must_use]
pub fn diagnose(bins: &[BinStat]) -> Diagnostics {
    let occupied: Vec<&BinStat> = bins.iter().filter(|b| b.count > 0).collect();
    let k = occupied.len();
    let n: u32 = bins.iter().map(|b| b.count).sum();
    if k < 2 || n <= k as u32 {
        return Diagnostics { chi2: 0.0, dof: 0, p_value: 1.0 };
    }
    let nf = n as f64;
    let grand = occupied.iter().map(|b| b.count as f64 * b.mean).sum::<f64>() / nf;
    let ss_b: f64 = occupied
        .iter()
        .map(|b| {
            let d = b.mean - grand;
            b.count as f64 * d * d
        })
        .sum();
    let ss_w: f64 = bins.iter().map(|b| b.m2).sum();
    let s2_w = ss_w / (n - k as u32) as f64;
    let (chi2, p) = anova_statistic(ss_b, s2_w, k - 1);
    Diagnostics { chi2, dof: k - 1, p_value: p }
}

fn anova_statistic(ss_b: f64, s2_w: f64, dof: usize) -> (f64, f64) {
    if s2_w <= 0.0 {
        return if ss_b <= 0.0 { (0.0, 1.0) } else { (f64::INFINITY, 0.0) };
    }
    let chi2 = ss_b / s2_w;
    (chi2, chi2_sf(chi2, dof))
}

/// Upper-tail survival function of chi-square with `dof` DOF.
#[must_use]
pub fn chi2_sf(x: f64, dof: usize) -> f64 {
    match (dof, x) {
        (0, _) => 1.0,
        (_, x) if !x.is_finite() => if x.is_nan() { 1.0 } else { 0.0 },
        (_, x) if x <= 0.0 => 1.0,
        (_, x) => gamma_q(dof as f64 / 2.0, x / 2.0),
    }
}

const GAMMA_ITMAX: usize = 300;
const GAMMA_EPS: f64 = 3.0e-12;
const GAMMA_FPMIN: f64 = 1.0e-300;

/// Regularized upper incomplete gamma Q(a, x) (series / continued fraction).
fn gamma_q(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 1.0;
    }
    if x < a + 1.0 {
        1.0 - gamma_series_p(a, x)
    } else {
        gamma_cf_q(a, x)
    }
}

fn gamma_series_p(a: f64, x: f64) -> f64 {
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..GAMMA_ITMAX {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * GAMMA_EPS {
            break;
        }
    }
    sum * (-x + a * x.ln() - libm::lgamma(a)).exp()
}

fn gamma_cf_q(a: f64, x: f64) -> f64 {
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / GAMMA_FPMIN;
    let mut d = 1.0 / b.max(GAMMA_FPMIN);
    let mut h = d;
    for i in 1..=GAMMA_ITMAX {
        let fi = i as f64;
        let an = -fi * (fi - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < GAMMA_FPMIN {
            d = GAMMA_FPMIN;
        }
        c = b + an / c;
        if c.abs() < GAMMA_FPMIN {
            c = GAMMA_FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < GAMMA_EPS {
            break;
        }
    }
    (-x + a * x.ln() - libm::lgamma(a)).exp() * h
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    /// Deterministic LCG + Box-Muller so structure tests never flake.
    pub(super) struct Lcg(u64);
    impl Lcg {
        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 11) as f64) / (1u64 << 53) as f64
        }
        fn normal(&mut self) -> f64 {
            let u1 = self.uniform().max(1e-12);
            let u2 = self.uniform();
            (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
        }
    }

    fn synthetic(n: usize, tow0: f64, amp: f64, sigma: f64, seed: u64) -> Vec<Sample> {
        let mut rng = Lcg(seed);
        (0..n)
            .map(|i| {
                let tow = tow0 + i as f64 * 30.0;
                let v = amp * (2.0 * PI * sidereal_phase(tow)).sin() + sigma * rng.normal();
                Sample { tow_s: tow, value: v }
            })
            .collect()
    }

    #[test]
    fn sidereal_phase_wraps_at_day_boundary() {
        // Just below one sidereal day -> phase just below 1.
        assert!((sidereal_phase(86_163.0) - 86_163.0 / 86_164.0).abs() < 1e-12);
        // Crossing the boundary wraps to a small positive phase.
        assert!(sidereal_phase(86_165.0) < 1e-4);
        // Arbitrary TOW beyond the period stays in [0, 1).
        for tow in [0.0, 100.0, 43_082.0, 86_164.0, 90_000.0, 604_800.0] {
            let p = sidereal_phase(tow);
            assert!((0.0..1.0).contains(&p), "phase {p} out of range for tow {tow}");
        }
    }

    #[test]
    fn sidereal_day_multiples_preserve_phase_week_boundary_shifts() {
        // Whole sidereal days bring the satellite geometry back -> same
        // phase (this is the entire basis of repeat-day stacking).
        for tow in [0.0, 100.0, 43_082.0, 86_163.5] {
            let base = sidereal_phase(tow);
            for k in 1..4 {
                let p = sidereal_phase(tow + f64::from(k) * SIDEREAL_PERIOD_S);
                assert!(
                    (p - base).abs() < 1e-9,
                    "day multiple changed phase: tow={tow} k={k} {base}->{p}"
                );
            }
        }
        // A GPS week (604800 s) is 7 SOLAR days but NOT an integer number
        // of sidereal days: the phase slips by 1652 s (~27.5 min) across a
        // week boundary. Expected, not a bug — sessions crossing week
        // boundaries see a one-off geometry offset.
        let expected_slip = 1_652.0 / SIDEREAL_PERIOD_S;
        assert!((sidereal_phase(604_800.0) - expected_slip).abs() < 1e-12);
    }

    #[test]
    fn day_boundary_bins_are_adjacent_across_wrap() {
        let b_before = phase_bin(sidereal_phase(86_160.0), 240);
        let b_after = phase_bin(sidereal_phase(86_165.0), 240);
        assert_eq!(b_before, 239);
        assert_eq!(b_after, 0);
    }

    #[test]
    fn phase_bin_edge_cases() {
        assert_eq!(phase_bin(0.0, 240), 0);
        assert_eq!(phase_bin(0.999_999, 240), 239);
        // Exactly 1.0 clamps into the last bin instead of overflowing.
        assert_eq!(phase_bin(1.0, 240), 239);
        assert_eq!(phase_bin(f64::NAN, 240), 0);
        assert_eq!(phase_bin(-0.5, 240), 0);
        assert_eq!(phase_bin(0.5, 0), 0);
    }

    #[test]
    fn chi2_sf_matches_known_critical_values() {
        // chi2(1 dof): p=0.05 at 3.8415.
        assert!((chi2_sf(3.841_459, 1) - 0.05).abs() < 1e-3);
        // chi2(2 dof): p=0.05 at 5.9915.
        assert!((chi2_sf(5.991_465, 2) - 0.05).abs() < 1e-3);
        assert!((chi2_sf(0.0, 4) - 1.0).abs() < 1e-12);
        // Monotonically decreasing in x.
        assert!(chi2_sf(10.0, 10) > chi2_sf(30.0, 10));
        assert!(chi2_sf(1_000.0, 10) < 1e-12);
    }

    #[test]
    fn fold_counts_means_and_std() {
        // Bins are ~359 s wide (86164/240): pick tows landing in bins
        // 0, 11, 24, plus one that wraps back into bin 0.
        let s = vec![
            Sample { tow_s: 0.0, value: 1.0 },      // bin 0
            Sample { tow_s: 4_000.0, value: 3.0 },  // bin 11
            Sample { tow_s: 8_700.0, value: 5.0 },  // bin 24
            Sample { tow_s: 86_170.0, value: 7.0 }, // wraps -> bin 0
        ];
        let bins = fold(&s, 240);
        assert_eq!(bins.len(), 240);
        let total: u32 = bins.iter().map(|b| b.count).sum();
        assert_eq!(total, 4);
        assert_eq!(bins[0].count, 2);
        assert!((bins[0].mean - 4.0).abs() < 1e-12);
        assert!((bins[0].std() - 3.0 * std::f64::consts::SQRT_2).abs() < 1e-9);
        assert_eq!(bins[11].count, 1);
        assert!((bins[11].mean - 3.0).abs() < 1e-12);
        assert_eq!(bins[24].count, 1);
        assert!((bins[24].mean - 5.0).abs() < 1e-12);
        // Non-finite samples are skipped, not propagated.
        let bad = vec![
            Sample { tow_s: 0.0, value: f64::NAN },
            Sample { tow_s: 1.0, value: 2.0 },
        ];
        let bins = fold(&bad, 240);
        assert_eq!(bins[0].count, 1);
        assert!((bins[0].mean - 2.0).abs() < 1e-12);
    }

    #[test]
    fn diagnose_detects_injected_sidereal_sinusoid() {
        let s = synthetic(2880, 3_600.0, 0.15, 0.02, 0xDEADBEEF);
        let diag = diagnose(&fold(&s, 240));
        assert!(diag.is_structured(), "expected structured, got {diag:?}");
        assert!(diag.p_value < 1e-6, "p too large: {}", diag.p_value);
        assert!(diag.chi2 / diag.dof as f64 > 2.0);
    }

    #[test]
    fn diagnose_pure_noise_reports_unstructured() {
        let s = synthetic(2880, 3_600.0, 0.0, 0.02, 0x12345678);
        let diag = diagnose(&fold(&s, 240));
        assert!(!diag.is_structured(), "noise flagged structured: {diag:?}");
        assert!(diag.p_value > 0.01);
        assert!((diag.chi2 / diag.dof as f64 - 1.0).abs() < 0.35);
    }

    #[test]
    fn diagnose_empty_or_degenerate_is_unstructured() {
        let d0 = diagnose(&fold(&[], 240));
        assert!(!d0.is_structured());
        assert_eq!(d0.dof, 0);
        let single = diagnose(&fold(&[Sample { tow_s: 5.0, value: 1.0 }], 240));
        assert!(!single.is_structured());
    }

    #[test]
    fn structure_threshold_is_strictly_less_than() {
        let mk = |p: f64| Diagnostics { chi2: 42.0, dof: 10, p_value: p };
        assert!(mk(STRUCTURE_P_THRESHOLD - 1e-12).is_structured());
        // Exactly AT the threshold does NOT count (strict comparison).
        assert!(!mk(STRUCTURE_P_THRESHOLD).is_structured());
        assert!(!mk(STRUCTURE_P_THRESHOLD + 1e-3).is_structured());
        // Non-finite p never claims structure.
        assert!(!mk(f64::NAN).is_structured());
        assert!(!Diagnostics { chi2: 1.0, dof: 0, p_value: 0.0 }.is_structured());
    }
}
