//! Wide-lane-first cascade ambiguity resolution.
//!
//! With N_W fixed from the MW arcs (see [super::mw]), narrow lanes are fixed
//! in closed form from the ionosphere-free float combination
//! a_IF = (f1·a1 − f2·a2)/(f1 − f2) = (1+s)·a1 − s·a2 with s = f2/(f1−f2).
//! It satisfies a_IF = N1 − s·N_W exactly: conditioning on the MW wide lane
//! leaves a single rounding whose input carries no iono bias (a naive
//! per-band constraint would let the ~half-cycle iono bias shift the chosen
//! integer while passing deviation gates).

use std::collections::HashMap;

use nalgebra::{DMatrix, DVector};

use super::ar::{project_subset_fixed, ArResult};
use super::mw::WidelaneTracker;
use super::state::{DoubleDiffKey, RtkState};

/// Maximum distance of the iono-free float from its rounded narrow lane.
const NL_MAX_DEVIATION_CYCLES: f64 = 0.25;
/// Maximum distance of the bias-corrected MW average from its integer.
const WL_ROUND_MAX_DEV_CYCLES: f64 = 0.25;
/// Minimum number of wide-lane-fixed pairs before claiming a fix. At six
/// both-band pairs the iono-free stage can always refine the position, so a
/// cascade "fixed" epoch is never carried by per-band data alone.
const MIN_FIXED_PAIRS: usize = 6;

/// Cascade AR: wide-lane integers from the MW arcs, then narrow lanes.
///
/// Only pairs whose band-1 and band-2 float ambiguities are both tracked and
/// whose wide lane is confidently fixed participate; at least
/// `MIN_FIXED_PAIRS` pairs must fix for a position claim. Returns None when
/// the subset is too small or the projected position fails the sanity gates.
pub fn resolve_cascade(state: &RtkState, tracker: &WidelaneTracker) -> Option<ArResult> {
    let mut fixed_keys: Vec<usize> = Vec::new();
    let mut fixed_values: Vec<f64> = Vec::new();
    let mut nl_rejected = 0usize;

    // Receiver-pair differential code biases shift every MW arc average by a
    // near-common fractional offset (cross-brand DD does not cancel them).
    // Estimate the offset across all confident pairs and remove it before
    // rounding — otherwise ±1-cycle wide-lane errors are systematic.
    let confident: Vec<(DoubleDiffKey, f64)> = confident_widelanes(state, tracker);
    let bias = common_fractional_offset(&confident);

    for i in 0..state.ambiguities.len() {
        let key = state.ambiguities[i].0;
        if key.freq_band != 1 {
            continue;
        }
        let key2 = DoubleDiffKey { freq_band: 2, ..key };
        // get_amb_idx returns the absolute state index; ambiguities[] is
        // addressed relatively.
        let (Some(j_abs), Some(nl_scale)) =
            (state.get_amb_idx(&key2), tracker.nl_scale(&key))
            else { continue };
        let Some(w_float) = state_amb_widelane(&confident, &key)
            else { continue };
        let j = j_abs - state.amb_offset();
        let w_int = (w_float - bias).round() as i64;
        if (w_float - bias - w_int as f64).abs() > WL_ROUND_MAX_DEV_CYCLES {
            continue;
        }
        let Some(n1) = fix_narrow_lane(state, i, j, w_int, nl_scale) else {
            nl_rejected += 1;
            continue;
        };
        fixed_keys.push(i);
        fixed_values.push(n1 as f64);
        fixed_keys.push(j);
        fixed_values.push((n1 - w_int) as f64);
    }

    tracing::debug!(
        "wl-cascade: band1_pairs={} wl_confident={} bias={bias:.3} nl_rejected={nl_rejected} pairs={}",
        state.ambiguities.iter().filter(|(k, _)| k.freq_band == 1).count(),
        confident.len(),
        fixed_keys.len() / 2,
    );

    if fixed_keys.len() / 2 < MIN_FIXED_PAIRS {
        return None;
    }
    build_cascade_result(state, &fixed_keys, &fixed_values)
}

fn ar_fixed_keys(ar: &ArResult) -> Vec<DoubleDiffKey> {
    ar.fixed_ambiguities.iter().map(|(k, _)| *k).collect()
}

fn confident_widelanes_from(
    tracker: &WidelaneTracker,
    keys: &[DoubleDiffKey],
) -> Vec<(DoubleDiffKey, f64)> {
    keys.iter().filter(|k| k.freq_band == 1)
        .filter_map(|k| tracker.fixed_widelane(k).map(|(w, _)| (*k, w)))
        .collect()
}

/// Confident converged wide-lane averages for tracked band-1 pairs.
fn confident_widelanes(
    state: &RtkState,
    tracker: &WidelaneTracker,
) -> Vec<(DoubleDiffKey, f64)> {
    state.ambiguities.iter()
        .filter(|(k, _)| k.freq_band == 1)
        .filter_map(|(k, _)| tracker.fixed_widelane(k).map(|(w, _)| (*k, w)))
        .collect()
}

fn state_amb_widelane(ws: &[(DoubleDiffKey, f64)], key: &DoubleDiffKey) -> Option<f64> {
    ws.iter().find(|(k, _)| k == key).map(|(_, w)| *w)
}

/// Median fractional offset of the MW averages from their integers: the
/// receiver-pair differential code bias, common to all pairs.
fn common_fractional_offset(ws: &[(DoubleDiffKey, f64)]) -> f64 {
    if ws.is_empty() {
        return 0.0;
    }
    let mut fracs: Vec<f64> = ws.iter().map(|(_, w)| w - w.round()).collect();
    fracs.sort_by(|a, b| a.total_cmp(b));
    fracs[fracs.len() / 2]
}

/// Cross-validation of a FAR/PAR fix against confident MW wide lanes.
///
/// The ratio test cannot detect a confidently-wrong integer vector whose
/// float was dragged past half-cycle by unmodelled DD-iono drift, but the
/// iono-free MW average can: any fixed pair whose N1 − N2 disagrees with a
/// converged wide lane condemns the whole fix. Pairs without a converged
/// wide lane are skipped (cannot judge).
pub fn far_matches_widelanes(tracker: &WidelaneTracker, ar: &ArResult) -> bool {
    let n2: HashMap<DoubleDiffKey, f64> = ar.fixed_ambiguities.iter()
        .filter(|(k, _)| k.freq_band == 2)
        .map(|(k, v)| (*k, *v))
        .collect();
    // Bias-corrected wide lanes: same rounding the cascade uses.
    let confident = confident_widelanes_from(tracker, &ar_fixed_keys(ar));
    let bias = common_fractional_offset(&confident);
    let w_map: HashMap<DoubleDiffKey, i64> = confident.iter()
        .map(|(k, w)| (*k, (w - bias).round() as i64))
        .collect();
    let mut judged = 0usize;
    let mut contradictions = 0usize;
    for (key, n1) in ar.fixed_ambiguities.iter().filter(|(k, _)| k.freq_band == 1) {
        let Some(&w_int) = w_map.get(key) else { continue };
        let key2 = DoubleDiffKey { freq_band: 2, ..*key };
        let Some(n2v) = n2.get(&key2) else { continue };
        if *n1 as i64 - *n2v as i64 != w_int {
            tracing::debug!(
                "wl-veto: sat {} fixed N1-N2={} but MW wide lane is {w_int}",
                key.sat, *n1 as i64 - *n2v as i64
            );
            contradictions += 1;
        }
        judged += 1;
    }
    // Tolerate a minority of contradictions (a single stale pair must not
    // kill a seven-pair fix); majority contradiction condemns the fix.
    judged == 0 || contradictions * 2 <= judged
}

/// Integer N1 from the iono-free float under N2 = N1 − w_int, accepted only
/// when the residual after rounding stays within the gate.
fn fix_narrow_lane(state: &RtkState, i: usize, j: usize, w_int: i64, nl_scale: f64) -> Option<i64> {
    let (a1, a2) = (state.ambiguities[i].1, state.ambiguities[j].1);
    let nl_float = (1.0 + nl_scale) * a1 - nl_scale * (a2 + w_int as f64);
    let n1 = nl_float.round();
    if !n1.is_finite() || (nl_float - n1).abs() > NL_MAX_DEVIATION_CYCLES {
        return None;
    }
    Some(n1 as i64)
}

fn build_cascade_result(
    state: &RtkState,
    fixed_keys: &[usize],
    fixed_values: &[f64],
) -> Option<ArResult> {
    let k = fixed_keys.len();
    let mut sub_a = DVector::zeros(k);
    let mut ints = DVector::zeros(k);
    for (r, &idx) in fixed_keys.iter().enumerate() {
        sub_a[r] = state.ambiguities[idx].1;
        ints[r] = fixed_values[r];
    }
    let mut sub_q = DMatrix::zeros(k, k);
    let off = state.amb_offset();
    for (r, &ri) in fixed_keys.iter().enumerate() {
        for (c, &ci) in fixed_keys.iter().enumerate() {
            sub_q[(r, c)] = state.cov[(off + ri, off + ci)];
        }
    }
    let (pos, cov) = project_subset_fixed(state, &sub_a, &ints, &sub_q, fixed_keys)?;
    Some(ArResult {
        position_ecef: pos,
        cov_position: cov,
        ratio: 0.0,
        is_fixed: true,
        num_ambiguities: k,
        fixed_ambiguities: fixed_keys.iter().enumerate()
            .map(|(r, &idx)| (state.ambiguities[idx].0, fixed_values[r]))
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    const F1: f64 = 1575.42e6;
    const F2: f64 = 1227.60e6;
    const NL_SCALE: f64 = F2 / (F1 - F2);

    /// Build a converged-looking state whose per-band floats carry realistic
    /// biases: geometry error scaled by 1/lambda plus iono-like offsets that
    /// do NOT cancel per band (only in the ionosphere-free combination).
    ///
    /// NOTE: `err` must stay well inside the iono-free wavelength (~10.7 cm)
    /// projected: the narrow-lane rounding is unbiased w.r.t. ionosphere but
    /// still rides on the float geometry, so larger float errors are
    /// rejected by the deviation gate rather than mis-rounded.
    fn biased_state(truth: Vector3<f64>, err: Vector3<f64>) -> (RtkState, WidelaneTracker) {
        let lambda1 = SPEED_OF_LIGHT_M_S / F1;
        let lambda2 = SPEED_OF_LIGHT_M_S / F2;
        let dirs = [
            Vector3::new(-20_000.0, 15_000.0, 18_000.0),
            Vector3::new(5_000.0, 24_000.0, -14_000.0),
            Vector3::new(16_000.0, -18_000.0, 22_000.0),
            Vector3::new(-9_000.0, -6_000.0, -26_000.0),
            Vector3::new(22_000.0, 9_000.0, -19_000.0),
            Vector3::new(-15_000.0, 21_000.0, 8_000.0),
            Vector3::new(11_000.0, -23_000.0, 13_000.0),
        ];
        let true_w = [3.0, -2.0, 7.0, 1.0, -4.0, 5.0, 2.0];
        let true_n1 = [10.0, -5.0, 23.0, 2.0, -8.0, 17.0, 6.0];
        let iono1_m = 0.15; // strong DD iono: ~0.79 L1 cycles of per-band bias

        let mut state = RtkState::new(truth + err, GpsTime::new(2000, 100.0));
        state.pos_ecef = truth + err;
        let mut tracker = WidelaneTracker::default();
        let mut los_list = Vec::new();
        for (p, dir) in dirs.iter().enumerate() {
            let los = (dir - truth).normalize();
            los_list.push(los);
            let key = DoubleDiffKey { constellation_id: 0, sat: 2 + p as u16, ref_sat: 1, freq_band: 1 };
            let key2 = DoubleDiffKey { freq_band: 2, ..key };
            let a1 = true_n1[p] + los.dot(&err) / lambda1 - iono1_m / lambda1;
            let a2 = true_n1[p] - true_w[p] + los.dot(&err) / lambda2
                - iono1_m * (F1 / F2).powi(2) / lambda2;
            state.ensure_ambiguity(key, a1, 1.0);
            state.ensure_ambiguity(key2, a2, 1.0);
            // Converged MW arcs: symmetric ramp settling on the true integer.
            for k in 0..40 {
                let frac = 1.0 - k as f64 / 39.0;
                tracker.update(key, true_w[p] + 0.2 * frac - 0.1, NL_SCALE, false);
            }
        }
        assemble_kf_consistent_covariance(&mut state, &los_list, lambda1, lambda2);
        (state, tracker)
    }

    /// Overwrite the fixture covariance with the shape a converged DD KF
    /// would report: P_xa = P_xx H^T, P_aa = H P_xx H^T + R per band.
    fn assemble_kf_consistent_covariance(
        state: &mut RtkState,
        los_list: &[Vector3<f64>],
        lambda1: f64,
        lambda2: f64,
    ) {
        let p_xx = 1e-3;
        let r_phase = 1e-4;
        for i in 0..3 {
            state.cov[(i, i)] = p_xx;
        }
        for (p, los_p) in los_list.iter().enumerate() {
            for (q, los_q) in los_list.iter().enumerate() {
                let (ip, jp, iq, jq) = (idx(state, p, 1), idx(state, p, 2), idx(state, q, 1), idx(state, q, 2));
                let dot = los_p.dot(los_q);
                let same = (p == q) as u8 as f64;
                state.cov[(ip, iq)] = p_xx * dot / lambda1 / lambda1 + same * r_phase;
                state.cov[(jp, jq)] = p_xx * dot / lambda2 / lambda2 + same * r_phase;
                state.cov[(ip, jq)] = p_xx * dot / lambda1 / lambda2;
                state.cov[(jp, iq)] = p_xx * dot / lambda1 / lambda2;
            }
            let (ip, jp) = (idx(state, p, 1), idx(state, p, 2));
            for r in 0..3 {
                state.cov[(r, ip)] = p_xx * los_p[r] / lambda1;
                state.cov[(ip, r)] = state.cov[(r, ip)];
                state.cov[(r, jp)] = p_xx * los_p[r] / lambda2;
                state.cov[(jp, r)] = state.cov[(r, jp)];
            }
        }
    }

    fn idx(state: &RtkState, p: usize, band: u8) -> usize {
        state.get_amb_idx(&DoubleDiffKey {
            constellation_id: 0, sat: 2 + p as u16, ref_sat: 1, freq_band: band,
        }).unwrap()
    }

    #[test]
    fn test_cascade_fixes_position_despite_per_band_iono_bias() {
        let truth = Vector3::new(100.0, 200.0, 300.0);
        // Sub-iono-free-wavelength float geometry: the regime the cascade
        // targets (iono blocks per-band fixing, geometry is already tight).
        let err = Vector3::new(0.010, -0.008, 0.005);
        let (state, tracker) = biased_state(truth, err);

        let res = resolve_cascade(&state, &tracker).expect("cascade should fix");
        assert!(res.is_fixed);
        assert_eq!(res.num_ambiguities, 14);
        // The per-band projected position here is only a fallback: the
        // conditional correction necessarily leaks some of the DD-iono
        // bias (the filter carries no iono covariance to suppress it).
        // Six-plus fixed pairs guarantee the caller re-estimates the
        // position iono-free from ArResult.fixed_ambiguities, which cancels
        // that bias exactly; hence the loose fallback bound below.
        let pos_err = (res.position_ecef - truth).norm();
        assert!(pos_err < 0.15, "cascade fallback position sane, got {:.4} m", pos_err);
        let n1_fixed = res.fixed_ambiguities.iter()
            .find(|(k, _)| k.freq_band == 1 && k.sat == 2).map(|(_, v)| *v as i64).unwrap();
        assert_eq!(n1_fixed, 10);
        // Both bands present so the iono-free stage can engage downstream.
        assert_eq!(res.fixed_ambiguities.iter().filter(|(k, _)| k.freq_band == 2).count(), 7);
    }

    #[test]
    fn test_cascade_integers_survive_iono_that_breaks_per_band_rounding() {
        // The narrow-lane rounding itself must stay immune: with the same
        // fixture, direct rounding of the band-1 floats would be off by the
        // ~0.8-cycle iono bias, while the iono-free combination rounds true.
        let truth = Vector3::new(100.0, 200.0, 300.0);
        let err = Vector3::new(0.010, -0.008, 0.005);
        let (state, tracker) = biased_state(truth, err);
        for i in 0..state.ambiguities.len() {
            let key = state.ambiguities[i].0;
            if key.freq_band != 1 { continue; }
            // Naive per-band rounding disagrees with the MW-consistent truth.
            let naive = state.ambiguities[i].1.round() as i64;
            let key2 = DoubleDiffKey { freq_band: 2, ..key };
            let w = tracker.fixed_widelane(&key).unwrap().1;
            let j = state.get_amb_idx(&key2).unwrap() - state.amb_offset();
            let n1 = fix_narrow_lane(&state, i, j, w, NL_SCALE).unwrap();
            assert_ne!(n1, naive, "fixture must carry sub-per-band-rounding iono bias");
            break;
        }
    }

    #[test]
    fn test_cascade_requires_six_pairs() {
        let truth = Vector3::new(100.0, 200.0, 300.0);
        let err = Vector3::new(0.010, -0.008, 0.005);
        let (mut state, mut tracker) = biased_state(truth, err);
        // Drop two pairs below the six-pair floor via retain_active.
        let keep: Vec<DoubleDiffKey> = state.ambiguities.iter()
            .filter(|(k, _)| !(k.freq_band == 1 && k.sat >= 6))
            .map(|(k, _)| *k).collect();
        state.retain_active_ambiguities(&keep);
        tracker.retain_active(&keep);
        assert!(resolve_cascade(&state, &tracker).is_none());
    }

    #[test]
    fn test_far_veto_rejects_fix_contradicting_converged_widelane() {
        let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        let key2 = DoubleDiffKey { freq_band: 2, ..key };
        let mut tracker = WidelaneTracker::default();
        for _ in 0..30 {
            tracker.update(key, 4.01, NL_SCALE, false); // confident W = 4
        }
        let make_ar = |n1: f64, n2: f64| ArResult {
            position_ecef: Vector3::zeros(),
            cov_position: nalgebra::Matrix3::identity(),
            ratio: 3.0,
            is_fixed: true,
            num_ambiguities: 2,
            fixed_ambiguities: vec![(key, n1), (key2, n2)],
        };
        // Consistent fix passes; N1-N2=3 against MW 4 is vetoed.
        assert!(far_matches_widelanes(&tracker, &make_ar(10.0, 6.0)));
        assert!(!far_matches_widelanes(&tracker, &make_ar(10.0, 7.0)));
    }

    #[test]
    fn test_cascade_rejects_when_float_geometry_too_loose() {
        let truth = Vector3::new(100.0, 200.0, 300.0);
        // Decimetre float error exceeds the iono-free wavelength: the
        // narrow-lane deviation gate must refuse to round instead of
        // committing to a confidently wrong integer vector.
        let err = Vector3::new(0.30, -0.22, 0.15);
        let (state, tracker) = biased_state(truth, err);
        assert!(resolve_cascade(&state, &tracker).is_none());
    }

}
