//! Sidereal-time stacking multipath filter (post-processing option).
//!
//! Multipath repeats with the *sidereal* day (86164 s), not the solar day:
//! the GPS satellite ground-track geometry recurs once per sidereal period,
//! so reflection-induced errors trace the same function of sidereal phase
//! across repeats. This module folds a time series onto sidereal phase bins
//! to (a) detect that periodicity diagnostically and (b) mitigate it by
//! subtracting per-bin mean corrections fitted strictly causally from the
//! first half of a session and applied to the second half only — no
//! look-ahead, no future data touches any epoch's correction.
//!
//! Layout:
//! - [`phase`]: product-grade primitives — phase conversion, folding,
//!   uniform-null diagnostics. Operate on plain `(tow, value)` samples and
//!   can be pointed at per-satellite residuals (the production path).
//! - [`mitigate`]: causal per-bin mean corrections from the session's
//!   first half.
//! - [`report`]: diagnostic dump structures, one-line summaries and CSV.
//! - This file: benchmark-grade trajectory glue — folds truth-referenced
//!   E/N/U errors (valid where monument truth exists, as in the eval
//!   harness) and shifts smoothed positions by the fitted correction.
//!
//! Interpretation caveat for single-sweep sessions: folding ONE solar day
//! onto sidereal phase cannot distinguish multipath periodicity from any
//! slow time-of-day systematics (iono/tropo drift maps monotonically onto
//! phase). A "STRUCTURED" verdict on <2 sweeps is therefore suggestive,
//! not conclusive; repeat-day consistency needs >=2 sweeps of data.
//!
//! Gate: off by default; enabled with `GNEISS_SIDEREAL=1`.

pub mod mitigate;
pub mod phase;
pub mod report;

#[cfg(test)]
pub(crate) mod testing;

use gneiss_core::coords::{ecef_to_llh, ecef_to_ned_matrix};
use nalgebra::Vector3;

use super::{percentile, SmoothedEpoch};
use phase::{Sample, SIDEREAL_PERIOD_S};
use report::{ChannelDump, SiderealReport};

pub use mitigate::{fit_first_half, Corrections, MIN_BIN_COUNT_DEFAULT};
pub use phase::{
    chi2_sf, diagnose, fold, phase_bin, sidereal_phase, BinStat, Diagnostics, DEFAULT_BINS,
    STRUCTURE_P_THRESHOLD,
};
pub use report::{write_diag_csv, CsvWritten};

/// Signed N/E/U error (m) of an ECEF position against ECEF truth.
#[must_use]
pub fn enu_error(pos: Vector3<f64>, truth: Vector3<f64>) -> Vector3<f64> {
    let llh = ecef_to_llh(truth);
    let ned = ecef_to_ned_matrix(llh) * (pos - truth);
    Vector3::new(ned.x, ned.y, -ned.z)
}

fn round_tow(tow: f64) -> u32 {
    tow.round() as u32
}

/// Horizontal/vertical error statistics over a trajectory slice.
#[derive(Debug, Clone, Copy, Default)]
pub struct HalfMetrics {
    pub n: usize,
    pub h_p50: f64,
    pub h_p95: f64,
    pub v_p50: f64,
    pub v_rms: f64,
}

#[must_use]
pub fn half_metrics(
    epochs: &[SmoothedEpoch],
    truth_at: &impl Fn(u32) -> Option<Vector3<f64>>,
) -> HalfMetrics {
    let mut h = Vec::new();
    let mut v = Vec::new();
    for ep in epochs {
        let Some(t) = truth_at(round_tow(ep.time.tow)) else { continue };
        let e = enu_error(ep.position_ecef, t);
        h.push((e.x * e.x + e.y * e.y).sqrt());
        v.push(e.z);
    }
    h.sort_by(|a, b| a.total_cmp(b));
    v.sort_by(|a, b| a.total_cmp(b));
    HalfMetrics {
        n: h.len(),
        h_p50: percentile(&h, 0.50),
        h_p95: percentile(&h, 0.95),
        v_p50: percentile(&v, 0.50),
        v_rms: (v.iter().map(|x| x * x).sum::<f64>() / v.len().max(1) as f64).sqrt(),
    }
}

fn channel_samples(
    traj: &[SmoothedEpoch],
    truth_at: &impl Fn(u32) -> Option<Vector3<f64>>,
) -> [Vec<Sample>; 3] {
    let mut sets = [Vec::new(), Vec::new(), Vec::new()];
    for ep in traj {
        let Some(t) = truth_at(round_tow(ep.time.tow)) else { continue };
        let e = enu_error(ep.position_ecef, t);
        let tow = ep.time.tow;
        sets[0].push(Sample { tow_s: tow, value: e.x });
        sets[1].push(Sample { tow_s: tow, value: e.y });
        sets[2].push(Sample { tow_s: tow, value: e.z });
    }
    sets
}

/// Shift second-half positions opposite their fitted sidereal-phase error.
/// The NED-frame correction (n, e, -u) is rotated back to ECEF per epoch.
fn shift_second_half(
    epochs: &mut [SmoothedEpoch],
    corrs: &[Corrections],
    truth_at: &impl Fn(u32) -> Option<Vector3<f64>>,
) {
    for ep in epochs {
        let Some(t) = truth_at(round_tow(ep.time.tow)) else { continue };
        let r_ned = ecef_to_ned_matrix(ecef_to_llh(t));
        let d_ned = Vector3::new(
            corrs[0].value_at(ep.time.tow),
            corrs[1].value_at(ep.time.tow),
            -corrs[2].value_at(ep.time.tow),
        );
        ep.position_ecef -= r_ned.transpose() * d_ned;
    }
}

fn fit_channels(
    sample_sets: &[Vec<Sample>; 3],
    n_bins: usize,
    min_bin_count: u32,
) -> (Vec<ChannelDump>, Vec<Corrections>) {
    let names = ["north", "east", "up"];
    let mut channels = Vec::with_capacity(3);
    let mut corrs = Vec::with_capacity(3);
    for (name, samples) in names.iter().zip(sample_sets.iter()) {
        let folded = fold(samples, n_bins);
        let corr = fit_first_half(samples, n_bins, min_bin_count);
        channels.push(ChannelDump {
            channel: (*name).to_string(),
            diag: diagnose(&folded),
            corr_rms_m: corr.rms(),
            active_bins: corr.active_bins(),
            bins: folded,
        });
        corrs.push(corr);
    }
    (channels, corrs)
}

/// Fold + diagnose + causally mitigate a smoothed trajectory in place.
///
/// The session is split chronologically at `traj.len()/2`; corrections are
/// fitted ONLY from first-half errors and applied ONLY to second-half
/// epochs, so no epoch ever consumes future information.
pub fn apply_to_trajectory(
    mut traj: Vec<SmoothedEpoch>,
    truth_at: impl Fn(u32) -> Option<Vector3<f64>>,
    n_bins: usize,
    min_bin_count: u32,
) -> (Vec<SmoothedEpoch>, SiderealReport) {
    let mid = traj.len() / 2;
    let before = half_metrics(&traj[mid..], &truth_at);
    let sample_sets = channel_samples(&traj, &truth_at);
    let (channels, corrs) = fit_channels(&sample_sets, n_bins, min_bin_count);
    shift_second_half(&mut traj[mid..], &corrs, &truth_at);
    let after = half_metrics(&traj[mid..], &truth_at);
    let report = SiderealReport {
        period_s: SIDEREAL_PERIOD_S,
        n_bins,
        split_index: mid,
        samples: traj.len(),
        channels,
        before,
        after,
    };
    (traj, report)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::post_process::sidereal::testing::{mk_epoch, truth_at, TRUTH};

    #[test]
    fn enu_error_axes_at_equator_prime_meridian() {
        // enu_error returns (x=north, y=east, z=up).
        // At lat=0/lon=0: up=+X, east=+Y, north=+Z.
        let e = enu_error(Vector3::new(TRUTH.x + 10.0, TRUTH.y, TRUTH.z), TRUTH);
        assert!((e.z - 10.0).abs() < 1e-6, "up {e:?}");
        assert!(e.x.abs() < 1e-6 && e.y.abs() < 1e-6);
        let e = enu_error(Vector3::new(TRUTH.x, TRUTH.y + 5.0, TRUTH.z), TRUTH);
        assert!((e.y - 5.0).abs() < 1e-6, "east {e:?}");
        let e = enu_error(Vector3::new(TRUTH.x, TRUTH.y, TRUTH.z + 3.0), TRUTH);
        assert!((e.x - 3.0).abs() < 1e-6, "north {e:?}");
    }

    #[test]
    fn half_metrics_percentiles_sane() {
        let mut traj: Vec<SmoothedEpoch> = Vec::new();
        for i in 0..100 {
            let tow = i as f64 * 30.0;
            // +Y is east here, +Z is north: h = sqrt(1^2+2^2), up = 0.
            traj.push(mk_epoch(tow, TRUTH + Vector3::new(0.0, 1.0, 2.0)));
        }
        let m = half_metrics(&traj, &truth_at);
        assert_eq!(m.n, 100);
        assert!((m.h_p50 - 5.0_f64.sqrt()).abs() < 1e-9);
        assert_eq!(m.h_p95, m.h_p50);
        assert!(m.v_p50.abs() < 1e-12);
        assert!(m.v_rms.abs() < 1e-12);
        // Empty slice -> zeroed metrics, no panic.
        let none = half_metrics(&[], &truth_at);
        assert_eq!(none.n, 0);
    }

    #[test]
    fn apply_to_trajectory_first_half_untouched_second_half_improved() {
        // Inject a sidereal sinusoid along EAST (+ECEF Y at this site).
        // 6000 epochs ~= 2.09 sidereal days so the first half holds a full
        // phase sweep (~4 samples/bin) and corrections activate.
        let mut traj: Vec<SmoothedEpoch> = Vec::new();
        let mut rng = testing::Lcg(42);
        for i in 0..6_000 {
            let tow = 3_600.0 + i as f64 * 30.0;
            let err = 0.5 * (2.0 * PI * sidereal_phase(tow)).sin() + 0.02 * rng.normal();
            traj.push(mk_epoch(tow, TRUTH + Vector3::new(0.0, err, 0.0)));
        }
        let original = traj.clone();
        let (fixed, report) = apply_to_trajectory(traj, truth_at, 240, 2);
        let mid = original.len() / 2;
        // Causality: first half must be bit-identical.
        for (a, b) in original[..mid].iter().zip(fixed[..mid].iter()) {
            assert_eq!(a.position_ecef, b.position_ecef, "first half was modified");
        }
        // Second half moved toward truth for the bulk of epochs.
        let improved = fixed[mid..]
            .iter()
            .zip(original[mid..].iter())
            .filter(|(f, o)| {
                (f.position_ecef - TRUTH).norm() < (o.position_ecef - TRUTH).norm()
            })
            .count();
        assert!(improved > (original.len() - mid) * 9 / 10, "improved only {improved}");
        assert!(report.after.h_p50 < report.before.h_p50);
        assert!(report.after.v_p50.abs() < 0.10);
        assert!(report.channels[1].diag.is_structured());
    }

    #[test]
    fn single_solar_day_session_corrections_confined_to_swept_bins() {
        // THE sub-sidereal-day red team, encoded: one solar day of 30 s
        // epochs sweeps ~half the phase circle ONCE. Each swept bin still
        // collects ~12 consecutive samples (bin ~= 359 s), but those are a
        // SINGLE repeat. Corrections therefore exist only for swept bins,
        // and because the second half keeps sweeping NEW phases, almost no
        // second-half epoch lands on a bin with first-half coverage — the
        // mitigation is structurally inert on sub-2-sweep sessions while
        // the diagnostic still runs.
        use super::testing::synthetic_samples as synthetic;
        let s = synthetic(2_880, 3_600.0, 0.30, 0.005, 11);
        let corr = fit_first_half(&s, 240, 2);
        assert!(
            (100..140).contains(&corr.active_bins()),
            "expected ~one sweep of bins, got {}",
            corr.active_bins()
        );
        // A phase never visited by the first half gets NO correction.
        let unswept_tow = 0.8 * SIDEREAL_PERIOD_S; // phase 0.8, swept only in 2nd half
        assert_eq!(corr.value_at(unswept_tow), 0.0);
        single_day_trajectory_only_overlap_moves(corr.active_bins());
    }

    /// End-to-end half of the single-day red team: only epochs whose BIN
    /// was active in the first half may move (the half-split falls inside
    /// one bin, and the session's 1.0027 sweeps re-enter the start bins at
    /// the very end).
    fn single_day_trajectory_only_overlap_moves(swept: usize) {
        use super::testing::{mk_epoch, truth_at, TRUTH};
        let mut traj: Vec<SmoothedEpoch> = Vec::new();
        for i in 0..2_880 {
            let tow = 3_600.0 + i as f64 * 30.0;
            let err = 0.30 * (2.0 * PI * sidereal_phase(tow)).sin();
            traj.push(mk_epoch(tow, TRUTH + Vector3::new(0.0, err, 0.0)));
        }
        let s_tows: Vec<f64> = (0..2_880).map(|i| 3_600.0 + i as f64 * 30.0).collect();
        let first_lo = phase_bin(sidereal_phase(s_tows[0]), 240);
        let first_hi = phase_bin(sidereal_phase(s_tows[1439]), 240);
        let original = traj.clone();
        let (fixed, report) = apply_to_trajectory(traj, truth_at, 240, 2);
        let mut changed = 0usize;
        for (a, b) in original.iter().zip(fixed.iter()) {
            if a.position_ecef != b.position_ecef {
                changed += 1;
                let bin = phase_bin(sidereal_phase(a.time.tow), 240);
                assert!(
                    bin >= first_lo && bin <= first_hi,
                    "moved epoch outside first-half bins: phase={}",
                    sidereal_phase(a.time.tow)
                );
            }
        }
        assert!(changed <= original.len() / 100, "too many moved: {changed}");
        assert!(report.channels.iter().all(|c| c.active_bins == swept));
        // With (almost) no corrections applied, pre/post metrics coincide.
        assert_eq!(report.before.h_p50, report.after.h_p50);
    }
}
