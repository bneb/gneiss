//! Distributional comparison for benchmark evaluation.
//!
//! [`metrics`] describes a single error distribution. This module answers the
//! two questions a single distribution cannot:
//!
//! 1. [`wilcoxon_signed_rank`] — given the *same* epochs evaluated two ways
//!    (e.g. ambiguity resolution on vs off), is one systematically better?
//!    This is the test that settles a p50-improves / p95-degrades trade: paired,
//!    nonparametric, and with no normality assumption on metre-scale errors.
//! 2. [`weibull_mle`] — how heavy is the tail? A fitted Weibull puts a single
//!    interpretable number on tail shape, so "how bad is the bad 1%" becomes a
//!    shape parameter rather than a second-guessed percentile.
//!
//! Both report an effect size, not just significance. A p-value says whether a
//! difference is detectable; the probability of superiority says how large it
//! is, which is what decides whether a change is worth keeping.

use alloc::vec::Vec;

/// Outcome of a paired Wilcoxon signed-rank test.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WilcoxonResult {
    /// Number of usable pairs (differences that are non-zero).
    pub n: usize,
    /// Pairs dropped because the two solutions agreed exactly. These carry no
    /// direction and are excluded rather than assigned to either group.
    pub n_zero_diff: usize,
    /// Sum of ranks where `a` exceeded `b`.
    pub w_plus: f64,
    /// Sum of ranks where `a` fell short of `b`.
    pub w_minus: f64,
    /// Two-sided p-value.
    pub p_value: f64,
    /// `P(a < b)`, treating smaller as better: the probability that a randomly
    /// drawn pair favours `a`. 0.5 is indistinguishable, 1.0 means `a` always wins.
    pub probability_of_superiority: f64,
}

impl WilcoxonResult {
    /// Sign of the paired effect: `+1` if `a` tends to be smaller (better).
    pub fn direction(&self) -> i32 {
        // Smaller error is better, so `a` wins when its differences are
        // negative, i.e. when the negative ranks dominate.
        match self.w_minus.partial_cmp(&self.w_plus) {
            Some(core::cmp::Ordering::Greater) => 1,
            Some(core::cmp::Ordering::Less) => -1,
            _ => 0,
        }
    }
}

/// Rank the absolute differences, assigning average ranks to ties.
fn average_ranks(abs_d: &[f64]) -> Vec<f64> {
    let n = abs_d.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| abs_d[a].partial_cmp(&abs_d[b]).unwrap_or(core::cmp::Ordering::Equal));
    let mut ranks = alloc::vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && abs_d[order[j + 1]] == abs_d[order[i]] {
            j += 1;
        }
        let avg = ((i + 1) + (j + 1)) as f64 / 2.0;
        for k in i..=j {
            ranks[order[k]] = avg;
        }
        i = j + 1;
    }
    ranks
}

/// Exact two-sided p-value by DP over the `W+` distribution. `O(n^2)`.
fn exact_w_plus_p(w_plus: f64, n: usize) -> f64 {
    // counts[w] = number of sign assignments with W+ == w, over ranks 1..=n.
    let total_u = n * (n + 1) / 2;
    let total = total_u as f64;
    let mut counts = alloc::vec![0f64; total_u + 1];
    counts[0] = 1.0;
    for r in 1..=n {
        let mut next = alloc::vec![0f64; total_u + 1];
        for (w, &c) in counts.iter().enumerate() {
            if c > 0.0 {
                next[w + r] += c; // rank assigned to b
                next[w] += c; // rank assigned to a
            }
        }
        counts = next;
    }
    // Two-sided: as far from the expectation as the observed statistic is.
    // Centring on the maximum instead would count every outcome when the
    // observed W+ is 0, returning p = 1 for a unanimous sample.
    let mean = total / 2.0;
    let obs_dev = (w_plus - mean).abs();
    let tail: f64 = counts
        .iter()
        .enumerate()
        .filter(|(w, _)| (*w as f64 - mean).abs() >= obs_dev - 1e-9)
        .map(|(_, &c)| c)
        .sum();
    (tail / (2f64).powi(n as i32)).min(1.0)
}

/// Normal approximation with tie and continuity correction.
fn normal_approx_p(w_plus: f64, w_minus: f64, n: usize, tie_sum: f64) -> f64 {
    let n_f = n as f64;
    // Var(W+) = [n(n+1)(2n+1) - tie_sum] / 24
    let var = (n_f * (n_f + 1.0) * (2.0 * n_f + 1.0) - tie_sum) / 24.0;
    if var <= 0.0 {
        return 1.0;
    }
    let diff = w_plus - w_minus;
    let cont = if diff >= 0.0 { 1.0 } else { -1.0 };
    let z = ((diff.abs() - cont) / var.sqrt()).max(0.0);
    // Two-sided normal tail via erfc.
    let p = libm::erfc(z / core::f64::consts::SQRT_2);
    p.min(1.0)
}

/// Largest n for which the exact DP is used.
const EXACT_MAX_N: usize = 20;

/// Paired Wilcoxon signed-rank test. Smaller error is treated as better.
///
/// `a` and `b` must describe the *same* epochs in the same order; pairing is the
/// entire point. Exact for `n <= 20`, normal approximation with tie and
/// continuity correction beyond that.
pub fn wilcoxon_signed_rank(a: &[f64], b: &[f64]) -> Result<WilcoxonResult, &'static str> {
    if a.len() != b.len() {
        return Err("paired samples must have equal length");
    }
    if a.len() < 2 {
        return Err("need at least 2 pairs");
    }
    let mut abs_d = Vec::new();
    let mut sign_pos = Vec::new();
    for (x, y) in a.iter().zip(b.iter()) {
        if !x.is_finite() || !y.is_finite() {
            return Err("samples must be finite");
        }
        let d = x - y;
        if d == 0.0 {
            continue;
        }
        abs_d.push(d.abs());
        sign_pos.push(d > 0.0);
    }
    let n = abs_d.len();
    if n == 0 {
        return Ok(WilcoxonResult {
            n: 0,
            n_zero_diff: a.len(),
            w_plus: 0.0,
            w_minus: 0.0,
            p_value: 1.0,
            probability_of_superiority: 0.5,
        });
    }
    let ranks = average_ranks(&abs_d);
    let mut w_plus = 0.0;
    let mut w_minus = 0.0;
    for (r, &pos) in ranks.iter().zip(sign_pos.iter()) {
        if pos {
            w_plus += r;
        } else {
            w_minus += r;
        }
    }
    // Tie correction: sum t^3 - t over groups of equal |d|.
    let mut tie_sum = 0.0;
    {
        let mut sorted = abs_d.clone();
        sorted.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));
        let mut i = 0;
        while i < n {
            let mut j = i;
            while j + 1 < n && sorted[j + 1] == sorted[i] {
                j += 1;
            }
            let t = (j + 1 - i) as f64;
            tie_sum += t * t * t - t;
            i = j + 1;
        }
    }
    let p_value = if n <= EXACT_MAX_N {
        exact_w_plus_p(w_plus, n)
    } else {
        normal_approx_p(w_plus, w_minus, n, tie_sum)
    };
    // P(a < b) from the rank sums (no ties, no zeros survive here).
    let w_total = (n * (n + 1) / 2) as f64;
    let pos = (w_minus - w_plus + w_total) / (2.0 * w_total);
    Ok(WilcoxonResult {
        n,
        n_zero_diff: a.len() - n,
        w_plus,
        w_minus,
        p_value,
        probability_of_superiority: pos.clamp(0.0, 1.0),
    })
}

/// Maximum-likelihood two-parameter Weibull fit (`shape`, `scale`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WeibullFit {
    /// Shape parameter `k`. Larger means the tail decays faster.
    pub shape: f64,
    /// Scale parameter `lambda` (the 63.2nd-percentile magnitude).
    pub scale: f64,
    /// Samples actually used (non-positive values are dropped).
    pub n: usize,
    /// Magnitude exceeded with probability 1e-3, or `None` if the fit failed.
    pub p999: Option<f64>,
}

impl WeibullFit {
    /// Quantile of the fitted distribution.
    pub fn quantile(&self, p: f64) -> f64 {
        self.scale * (-libm::log(1.0 - p.clamp(0.0, 0.9999999))).powf(1.0 / self.shape)
    }
}

/// Fit `F(x) = 1 - exp(-(x/lambda)^k)` by maximum likelihood.
///
/// Non-positive samples are dropped: a zero magnitude is degenerate for a
/// Weibull and an error of exactly 0 means "no error", not "infinitely fast".
pub fn weibull_mle(sample: &[f64]) -> Result<WeibullFit, &'static str> {
    let x: Vec<f64> = sample.iter().copied().filter(|v| *v > 0.0 && v.is_finite()).collect();
    let n = x.len();
    if n < 2 {
        return Err("need at least 2 positive samples for a Weibull fit");
    }
    let lnx: Vec<f64> = x.iter().map(|v| libm::log(*v)).collect();
    // Moment estimators give a robust starting point; the fixed point refines it.
    let n_f = n as f64;
    // Profiled log-likelihood with lambda eliminated:
    //   g(k) = n ln k - n k ln lambda(k) + (k-1) sum ln x - n
    //   lambda(k) = (sum x^k / n)^(1/k)
    // Maximised by ternary search in ln k. A closed-form fixed point exists
    // (k = n / (sum x^k ln x - xbar sum x^k)) but it oscillates and collapses
    // to ~0 on heavy-tailed data, and a hand-derived Newton derivative is easy
    // to get subtly wrong. The search is cheap, bounded, and always converges.
    let s1: f64 = lnx.iter().sum();
    let prof = |k: f64| -> f64 {
        let b: f64 = x.iter().map(|v| v.powf(k)).sum::<f64>() / n_f;
        if b <= 0.0 || !b.is_finite() {
            return f64::NEG_INFINITY;
        }
        let ln_lambda = libm::log(b) / k;
        n_f * libm::log(k) - n_f * k * ln_lambda + (k - 1.0) * s1 - n_f
    };
    let (mut lo, mut hi) = (core::f64::consts::LN_2 * 0.05, core::f64::consts::LN_2 * 100.0);
    for _ in 0..200 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if prof(m1.exp()) < prof(m2.exp()) {
            lo = m1;
        } else {
            hi = m2;
        }
    }
    let k = ((lo + hi) / 2.0).exp();
    let xk: Vec<f64> = x.iter().map(|v| v.powf(k)).collect();
    let mean_xk: f64 = xk.iter().sum::<f64>() / n as f64;
    if mean_xk <= 0.0 {
        return Err("degenerate sample");
    }
    let lambda = mean_xk.powf(1.0 / k);
    if !lambda.is_finite() || lambda <= 0.0 || !k.is_finite() || k <= 0.0 {
        return Err("Weibull fit did not converge to a valid parameter pair");
    }
    Ok(WeibullFit { shape: k, scale: lambda, n, p999: None })
}

/// Fit and attach the 99.9th percentile in one step.
pub fn weibull_with_tail(sample: &[f64]) -> Result<WeibullFit, &'static str> {
    let mut fit = weibull_mle(sample)?;
    fit.p999 = Some(fit.quantile(0.999));
    Ok(fit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn wilcoxon_detects_a_consistent_improvement() {
        let a = vec![1.0, 1.1, 0.9, 1.05, 0.95, 1.02];
        let b = vec![2.0, 2.1, 1.9, 2.05, 1.95, 2.02];
        let r = wilcoxon_signed_rank(&a, &b).expect("fit");
        assert_eq!(r.n, 6);
        assert_eq!(r.w_plus, 0.0, "a is smaller everywhere");
        assert_eq!(r.w_minus, 21.0, "ranks 1..=6 all go to the negative side");
        assert!(r.p_value < 0.05, "p = {}", r.p_value);
        assert_eq!(r.direction(), 1, "a is smaller, so a wins");
    }

    #[test]
    fn wilcoxon_is_null_when_the_samples_match() {
        let x = vec![0.3, 1.2, 0.7, 2.1, 0.4, 1.7, 0.9, 1.1];
        let r = wilcoxon_signed_rank(&x, &x).expect("fit");
        assert_eq!(r.n, 0, "all differences are zero");
        assert_eq!(r.n_zero_diff, 8);
        assert_eq!(r.p_value, 1.0);
    }

    #[test]
    fn wilcoxon_rejects_mismatched_lengths() {
        assert!(wilcoxon_signed_rank(&[1.0, 2.0], &[1.0]).is_err());
    }

    #[test]
    fn wilcoxon_probability_of_superiority_brackets_half() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![1.5, 2.5, 3.5, 4.5];
        let r = wilcoxon_signed_rank(&a, &b).expect("fit");
        assert!(r.probability_of_superiority > 0.5, "a is uniformly smaller");
    }

    #[test]
    fn weibull_recovers_known_parameters() {
        // Deterministic quantile grid of a Weibull(k=2, lambda=3).
        let k = 2.0f64;
        let lambda = 3.0f64;
        let sample: Vec<f64> = (1..=2000)
            .map(|i| {
                let p = i as f64 / 2001.0;
                lambda * (-libm::log(1.0 - p)).powf(1.0 / k)
            })
            .collect();
        let fit = weibull_with_tail(&sample).expect("fit");
        assert!((fit.shape - k).abs() < 0.05, "shape {} vs {}", fit.shape, k);
        assert!((fit.scale - lambda).abs() < 0.05, "scale {} vs {}", fit.scale, lambda);
        let p999 = fit.p999.expect("tail");
        let true_p999 = lambda * (-libm::log(1.0 - 0.999)).powf(1.0 / k);
        assert!((p999 - true_p999).abs() < 0.1);
    }

    #[test]
    fn weibull_dropped_zeros_are_reported() {
        let sample = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let fit = weibull_mle(&sample).expect("fit");
        assert_eq!(fit.n, 6, "the zero is dropped");
    }

    #[test]
    fn weibull_shape_distinguishes_heavy_from_light_tails() {
        // Heavy tail: many small values with a few very large ones.
        let mut heavy: Vec<f64> = (1..=200).map(|i| (i as f64).powf(2.0)).collect();
        heavy.extend((1..=20).map(|i| i as f64 * 500.0));
        let light: Vec<f64> = (1..=200).map(|i| (i as f64).powf(0.5)).collect();
        let h = weibull_mle(&heavy).expect("heavy");
        let l = weibull_mle(&light).expect("light");
        assert!(
            h.shape < l.shape,
            "heavy tail should fit a smaller shape: {} vs {}",
            h.shape,
            l.shape
        );
    }
}
