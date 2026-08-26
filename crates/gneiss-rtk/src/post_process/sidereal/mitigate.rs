//! Causal sidereal-phase corrections: per-bin means fitted from the FIRST
//! HALF of a session only, applied downstream to second-half epochs.

use super::phase::{fold, phase_bin, sidereal_phase, Sample};

/// Default minimum first-half observations before a bin earns a correction.
pub const MIN_BIN_COUNT_DEFAULT: u32 = 2;

/// Per-bin mean corrections fitted from the FIRST HALF of a session only.
///
/// Bins with fewer than `min_bin_count` first-half observations stay at
/// zero correction. Note the repeat-day caveat: a 359 s bin collects ~12
/// consecutive samples per sweep, so a session shorter than two sidereal
/// days never yields a true multi-repeat average — and because the sweep
/// advances monotonically, second-half epochs mostly land on bins the
/// first half never visited. Net effect: bin-mean mitigation is
/// structurally inert below ~2 sweeps (see the mod.rs gating test),
/// while diagnostics remain valid for any length.
#[derive(Debug, Clone, Default)]
pub struct Corrections {
    n_bins: usize,
    values: Vec<f64>,
    active: usize,
}

#[must_use]
pub fn fit_first_half(samples: &[Sample], n_bins: usize, min_bin_count: u32) -> Corrections {
    let mid = samples.len() / 2;
    let folded = fold(&samples[..mid], n_bins);
    let mut values = vec![0.0; n_bins];
    let mut active = 0usize;
    for (i, b) in folded.iter().enumerate() {
        if b.count >= min_bin_count && b.mean.is_finite() {
            values[i] = b.mean;
            active += 1;
        }
    }
    Corrections { n_bins, values, active }
}

impl Corrections {
    #[must_use]
    pub fn value_at(&self, tow_s: f64) -> f64 {
        if self.n_bins == 0 {
            return 0.0;
        }
        self.values[phase_bin(sidereal_phase(tow_s), self.n_bins)]
    }

    /// RMS over active (gated-in) bins only.
    #[must_use]
    pub fn rms(&self) -> f64 {
        if self.active == 0 {
            return 0.0;
        }
        (self.values.iter().map(|v| v * v).sum::<f64>() / self.active as f64).sqrt()
    }

    #[must_use]
    pub fn active_bins(&self) -> usize {
        self.active
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::post_process::sidereal::testing::{synthetic_samples as synthetic, Lcg};

    /// Two-plus sidereal days of 30 s epochs: the first half then contains
    /// a full phase sweep, so every bin earns >=2 repeats and corrections
    /// activate. (One solar day = 2880 epochs is NOT enough — see the
    /// mod.rs single-day gating test.)
    fn multi_day_synthetic(amp: f64, sigma: f64, seed: u64) -> Vec<Sample> {
        synthetic(6_000, 3_600.0, amp, sigma, seed)
    }

    #[test]
    fn mitigation_removes_at_least_80pct_of_sinusoid_second_half() {
        let s = multi_day_synthetic(0.15, 0.02, 0xFEEDFACE);
        let mid = s.len() / 2;
        let corr = fit_first_half(&s, 240, 2);
        assert!(corr.active_bins() > 200, "too few active bins");
        let sigma = 0.02_f64;
        let sys_amp = |vals: &[f64]| {
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let var = vals.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / vals.len() as f64;
            (var - sigma * sigma).max(0.0).sqrt()
        };
        let raw_second: Vec<f64> = s[mid..].iter().map(|x| x.value).collect();
        let fixed_second: Vec<f64> =
            s[mid..].iter().map(|x| x.value - corr.value_at(x.tow_s)).collect();
        let amp_before = sys_amp(&raw_second);
        let amp_after = sys_amp(&fixed_second);
        assert!(amp_before > 0.10, "sinusoid missing before mitigation");
        assert!(
            amp_after <= 0.20 * amp_before,
            "removal insufficient: before={amp_before:.4} after={amp_after:.4}"
        );
    }

    #[test]
    fn mitigation_on_pure_noise_changes_little_and_stays_small() {
        let s = multi_day_synthetic(0.0, 0.02, 0x0BADC0DE);
        let mid = s.len() / 2;
        let corr = fit_first_half(&s, 240, 2);
        // Corrections are frozen first-half bin noise (~sigma/sqrt(n/bin)).
        assert!(corr.rms() < 0.05, "corrections too large: {}", corr.rms());
        let rms =
            |vals: &[f64]| (vals.iter().map(|v| v * v).sum::<f64>() / vals.len() as f64).sqrt();
        let raw: Vec<f64> = s[mid..].iter().map(|x| x.value).collect();
        let fixed: Vec<f64> = s[mid..].iter().map(|x| x.value - corr.value_at(x.tow_s)).collect();
        assert!(
            rms(&fixed) <= 1.10 * rms(&raw),
            "mitigation degraded noise: {} -> {}",
            rms(&raw),
            rms(&fixed)
        );
    }

    #[test]
    fn short_session_below_one_sidereal_day_gates_corrections() {
        // 400 samples = 3.33 h << 1 sidereal day; only bins swept during
        // those hours hold samples at all.
        let s = synthetic(400, 0.0, 0.15, 0.02, 7);
        let corr = fit_first_half(&s, 240, 2);
        assert!(
            corr.active_bins() < 120,
            "expected sparse coverage, got {}",
            corr.active_bins()
        );
        // Degenerate guard: empty session yields inert corrections.
        let none = fit_first_half(&[], 240, 2);
        assert_eq!(none.active_bins(), 0);
        assert_eq!(none.value_at(12_345.0), 0.0);
    }

    #[test]
    fn sinusoid_phase_recovery_is_shift_robust() {
        // A cosine (phase offset pi/2) is recovered equally well — guards
        // against an implementation secretly assuming sin-only structure.
        let mut rng = Lcg(5);
        let s: Vec<Sample> = (0..6_000)
            .map(|i| {
                let tow = 3_600.0 + i as f64 * 30.0;
                let v = 0.12 * (2.0 * PI * sidereal_phase(tow) + PI / 2.0).cos()
                    + 0.02 * rng.normal();
                Sample { tow_s: tow, value: v }
            })
            .collect();
        let corr = fit_first_half(&s, 240, 2);
        assert!(corr.rms() > 0.05 && corr.rms() < 0.20);
    }
}
