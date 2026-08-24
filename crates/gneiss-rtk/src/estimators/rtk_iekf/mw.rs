//! Melbourne–Wübbena observable construction and per-pair arc tracking.
//!
//! The double-difference MW combination is geometry-free and first-order
//! ionosphere-free:
//!
//!   mw      = (∇Δφ1 − ∇Δφ2) − ∇ΔR_N / λ_W          [cycles]
//!   ∇ΔR_N   = (f1·∇ΔP1 + f2·∇ΔP2)/(f1 + f2)         [m]
//!   λ_W     = c/(f1 − f2)                           [m]
//!
//! so its arc average converges to the integer wide-lane ambiguity
//! N_W = N1 − N2 even at long baselines where the L1 DD ionosphere bias
//! (>= 9 cm at 50 km, about half an L1 cycle) makes direct per-band fixing
//! unreliable.

use std::collections::HashMap;

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::obs::SatObs;
use gneiss_core::sat::SatelliteId;

use super::state::DoubleDiffKey;

/// Minimum epochs in an MW arc average before it may be declared fixed.
pub const MIN_TRACK_EPOCHS: u32 = 20;
/// Maximum standard error of the arc average for a wide-lane fix attempt.
pub const MAX_SIGMA_MEAN_CYCLES: f64 = 0.12;
/// Maximum distance of the arc average from the nearest integer.
pub const MAX_DEVIATION_CYCLES: f64 = 0.25;
/// Single-epoch MW innovation beyond this many cycles restarts the arc
/// (catches base-side slips that rover-only detectors cannot see).
const SLIP_INNOVATION_CYCLES: f64 = 1.0;
/// Epochs before the innovation slip gate arms (mean still unreliable).
const INNOVATION_ARM_EPOCHS: u32 = 5;

/// Running mean/variance of the DD MW observable for one pair (Welford).
#[derive(Debug, Clone, Copy)]
pub struct MwTrack {
    mean: f64,
    m2: f64,
    n: u32,
}

impl MwTrack {
    fn new(first: f64) -> Self {
        Self { mean: first, m2: 0.0, n: 1 }
    }

    fn push(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        self.m2 += delta * (x - self.mean);
    }

    pub fn mean(&self) -> f64 {
        self.mean
    }

    pub fn count(&self) -> u32 {
        self.n
    }

    /// Sample standard deviation; zero while the track has a single epoch.
    pub fn sigma(&self) -> f64 {
        if self.n < 2 { 0.0 } else { (self.m2 / (self.n - 1) as f64).sqrt() }
    }

    /// Standard error of the arc mean (the rounding-relevant precision).
    fn sigma_mean(&self) -> f64 {
        self.sigma() / (self.n.max(1) as f64).sqrt()
    }

    /// Absorb one epoch; a large innovation restarts the arc (cycle slip).
    fn absorb(&mut self, x: f64) {
        if self.n >= INNOVATION_ARM_EPOCHS && (x - self.mean).abs() > SLIP_INNOVATION_CYCLES {
            *self = MwTrack::new(x);
        } else {
            self.push(x);
        }
    }
}

/// Per-pair arc averages of the double-difference MW observable plus the
/// narrow-lane scale factor s = f2/(f1-f2) captured when the pair was formed.
///
/// Keys are the band-1 DD keys; a reference-satellite switch changes the key
/// and therefore starts a fresh arc automatically.
#[derive(Debug, Default, Clone)]
pub struct WidelaneTracker {
    tracks: HashMap<DoubleDiffKey, MwTrack>,
    nl_scales: HashMap<DoubleDiffKey, f64>,
    /// Network-solved satellite wide-lane UPDs (cycles, sum-zero): applied
    /// before rounding so MW arc means round to the true integer despite
    /// per-satellite biases.
    pub sat_upd: Option<HashMap<u16, f64>>,
}

impl WidelaneTracker {
    /// Absorb one DD MW observation; `slip` forces an arc restart.
    pub fn update(&mut self, key: DoubleDiffKey, mw_cycles: f64, nl_scale: f64, slip: bool) {
        let next = match self.tracks.get_mut(&key) {
            Some(track) if !slip => {
                track.absorb(mw_cycles);
                *track
            }
            _ => MwTrack::new(mw_cycles),
        };
        self.tracks.insert(key, next);
        self.nl_scales.insert(key, nl_scale);
    }

    /// Narrow-lane scale factor captured for a pair, if any.
    pub fn nl_scale(&self, key: &DoubleDiffKey) -> Option<f64> {
        self.nl_scales.get(key).copied()
    }

    /// Drop tracks whose pair disappeared from the active DD set.
    pub fn retain_active(&mut self, active_keys: &[DoubleDiffKey]) {
        self.tracks.retain(|k, _| active_keys.contains(k));
    }

    /// Forget one pair's arc (cycle-slip evidence): the next observation
    /// starts a fresh average.
    pub fn reset_pair(&mut self, key: &DoubleDiffKey) {
        self.tracks.remove(key);
    }

    /// Converged arc means `(mean, count)` per tracked pair.
    pub fn arc_means(&self) -> HashMap<DoubleDiffKey, (f64, u32)> {
        self.tracks.iter()
            .filter(|(_, t)| t.count() >= MIN_TRACK_EPOCHS)
            .map(|(k, t)| (*k, (t.mean(), t.count())))
            .collect()
    }

    /// Fixed `(w_float, w_int)` once the arc average converged near an integer.
    pub fn fixed_widelane(&self, key: &DoubleDiffKey) -> Option<(f64, i64)> {
        let track = self.tracks.get(key)?;
        if track.count() < MIN_TRACK_EPOCHS || track.sigma_mean() > MAX_SIGMA_MEAN_CYCLES {
            tracing::debug!(
                "wl-track {:?}: n={} sigma_mean={:.3} (gates n>={MIN_TRACK_EPOCHS}, sigma<={MAX_SIGMA_MEAN_CYCLES})",
                key.sat, track.count(), track.sigma_mean()
            );
            return None;
        }
        let mut w = track.mean();
        if let Some(upd) = self.sat_upd.as_ref() {
            if let (Some(&us), Some(&ur)) = (upd.get(&key.sat), upd.get(&key.ref_sat)) {
                w -= us - ur;
            }
        }
        let w_int = w.round();
        if (w - w_int).abs() > MAX_DEVIATION_CYCLES {
            tracing::debug!("wl-track {} deviates after UPD correction: w={w:.3}", key.sat);
            return None;
        }
        Some((w, w_int as i64))
    }
}

/// DD MW observable in cycles from one satellite pair, order
/// `[rov_sat, rov_ref, bas_sat, bas_ref]` (matching repo differencing).
/// Returns `(mw_cycles, nl_scale)` with nl_scale = f2/(f1-f2).
pub fn mw_dd_cycles(
    f1_hz: f64,
    f2_hz: f64,
    phi1: [f64; 4],
    phi2: [f64; 4],
    p1_m: [f64; 4],
    p2_m: [f64; 4],
) -> Option<(f64, f64)> {
    if f1_hz <= 0.0 || f2_hz <= 0.0 || (f1_hz - f2_hz).abs() < 1e6 {
        return None;
    }
    let dd = |v: [f64; 4]| (v[0] - v[1]) - (v[2] - v[3]);
    let lambda_wl = SPEED_OF_LIGHT_M_S / (f1_hz - f2_hz);
    let rn = (f1_hz * dd(p1_m) + f2_hz * dd(p2_m)) / (f1_hz + f2_hz);
    Some(((dd(phi1) - dd(phi2)) - rn / lambda_wl, f2_hz / (f1_hz - f2_hz)))
}


#[cfg(test)]
mod tests {
    use super::*;

    const F1: f64 = 1575.42e6;
    const F2: f64 = 1227.60e6;
    const NL_SCALE: f64 = F2 / (F1 - F2);

    #[test]
    fn test_mw_recovers_widelane_through_iono_and_geometry() {
        let lambda_wl = SPEED_OF_LIGHT_M_S / (F1 - F2);
        let lambda1 = SPEED_OF_LIGHT_M_S / F1;
        let lambda2 = SPEED_OF_LIGHT_M_S / F2;
        for &(geom_m, iono1_m, n1, n2) in &[
            (1200.5, 0.0, 3.0, -1.0),
            (25000.75, 0.10, 12.0, 5.0),
            (-8000.25, 0.45, -4.0, 9.0),
        ] {
            let iono2_m = iono1_m * (F1 / F2).powi(2);
            let dd_phi1 = geom_m / lambda1 + n1 - iono1_m / lambda1;
            let dd_phi2 = geom_m / lambda2 + n2 - iono2_m / lambda2;
            let dd_p1 = geom_m + iono1_m;
            let dd_p2 = geom_m + iono2_m;
            let rn = (F1 * dd_p1 + F2 * dd_p2) / (F1 + F2);
            let mw = (dd_phi1 - dd_phi2) - rn / lambda_wl;
            assert!((mw - (n1 - n2)).abs() < 1e-9, "MW must equal N1-N2, got {mw} vs {}", n1 - n2);
        }
    }

    #[test]
    fn test_tracker_converges_to_integer_after_noise() {
        let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        let mut tracker = WidelaneTracker::default();
        // Deterministic pseudo-noise around true W = 5; epoch spread ~0.55 cyc.
        let fracs = [0.4, -0.6, 0.2, -0.9, 0.7, -0.3, 0.9, -0.5, 0.1, -0.8];
        for k in 0..60 {
            let noisy = 5.0 + fracs[k % fracs.len()] * 0.55;
            tracker.update(key, noisy, NL_SCALE, false);
        }
        let (w, w_int) = tracker.fixed_widelane(&key).expect("should converge");
        assert_eq!(w_int, 5);
        assert!((w - 5.0).abs() <= MAX_DEVIATION_CYCLES);
    }

    #[test]
    fn test_tracker_rejects_biased_mean_and_resets_on_jump() {
        let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        let mut tracker = WidelaneTracker::default();
        for _ in 0..30 {
            tracker.update(key, 5.02, NL_SCALE, false);
        }
        assert!(tracker.fixed_widelane(&key).is_some());
        // A +2-cycle jump must reset the arc: no fix until it re-converges.
        tracker.update(key, 7.03, NL_SCALE, false);
        assert!(tracker.fixed_widelane(&key).is_none());
        for _ in 0..30 {
            tracker.update(key, 7.01, NL_SCALE, false);
        }
        assert_eq!(tracker.fixed_widelane(&key).unwrap().1, 7);
    }

}


/// Phase and code observables for one frequency band across a satellite pair
/// on both receivers, order `[rov_sat, rov_ref, bas_sat, bas_ref]`.
fn band_quad(rov_s: &SatObs, rov_ref: &SatObs, bas_s: &SatObs, bas_ref: &SatObs, band: u8) -> Option<([f64; 4], [f64; 4])> {
    Some((
        [
            rov_s.get_observable_phase(band)?,
            rov_ref.get_observable_phase(band)?,
            bas_s.get_observable_phase(band)?,
            bas_ref.get_observable_phase(band)?,
        ],
        [
            rov_s.get_observable(band)?,
            rov_ref.get_observable(band)?,
            bas_s.get_observable(band)?,
            bas_ref.get_observable(band)?,
        ],
    ))
}

/// Absorb one DD MW observation for a satellite pair straight from raw obs.
///
/// Arc resets fire on LLI flags from EITHER receiver on EITHER band — the
/// rover-only slip detector cannot see base-side slips, which would
/// otherwise shift every wide lane identically and invisibly.
#[allow(clippy::too_many_arguments)] // same obs bundle as build_single_dd_pair
pub fn update_tracker_from_obs(
    tracker: &mut WidelaneTracker,
    sat_id: SatelliteId,
    ref_sat_id: u16,
    rov_s: &SatObs,
    bas_s: &SatObs,
    rov_ref: &SatObs,
    bas_ref: &SatObs,
    glo_k: i8,
) {
    let key = DoubleDiffKey {
        constellation_id: sat_id.constellation as u8,
        sat: sat_id.prn as u16,
        ref_sat: ref_sat_id,
        freq_band: 1,
    };
    // Secondary band: GPS/GLONASS L2 (2); Galileo exports carry E5a on
    // band 5 instead of L2. Policy matches all other consumers: require
    // band-2 phase at ALL FOUR stations so a single dropped slot cannot
    // silently switch an arc between E5b- and E5a-based wide lanes.
    let b2 = if rov_s.get_observable_phase(2).is_some()
        && bas_s.get_observable_phase(2).is_some()
        && rov_ref.get_observable_phase(2).is_some()
        && bas_ref.get_observable_phase(2).is_some()
    {
        2
    } else {
        5
    };
    let f1 = gneiss_core::signal::get_frequency(sat_id, 1, glo_k);
    let f2 = gneiss_core::signal::get_frequency(sat_id, b2, glo_k);
    let band1 = band_quad(rov_s, rov_ref, bas_s, bas_ref, 1);
    let band2 = band_quad(rov_s, rov_ref, bas_s, bas_ref, b2);
    let slip = [rov_s, rov_ref, bas_s, bas_ref].iter().any(|o| {
        o.get_lli(1).is_some_and(|l| l & 1 != 0)
            || o.get_lli(b2).is_some_and(|l| l & 1 != 0)
    });
    if let (Some((phi1, code1)), Some((phi2, code2))) = (band1, band2) {
        if let Some((mw_cycles, nl_scale)) = mw_dd_cycles(f1, f2, phi1, phi2, code1, code2) {
            tracker.update(key, mw_cycles, nl_scale, slip);
        }
    }
}


// ---------------------------------------------------------------------------
// Network UPD estimation from phase-only wide-lane arc means
// ---------------------------------------------------------------------------

/// Result of the cross-base satellite wide-lane UPD decomposition.
#[derive(Debug, Clone, Default)]
pub struct NetworkUpdSolution {
    /// Satellite wide-lane UPD estimates (cycles), sum-to-zero constrained.
    pub sat_upd: HashMap<u16, f64>,
    /// Post-fit residual of every input observation (cycles).
    pub residuals: Vec<(DoubleDiffKey, f64)>,
    /// RMS of post-fit residuals across all used observations.
    pub residual_rms: f64,
}

/// Solve satellite wide-lane UPDs from per-base PHASE-only wide-lane arc
/// means.
///
/// Inputs are the converged means of `pw = dd_phi1 - dd_phi2 - geom_wl`
/// (geometry and troposphere removed with the known trajectory), keyed by
/// band-1 DD pairs. For each observation:
///
///   pw_arc_mean = N_W + u_sat - u_ref  (+ noise), all in cycles
///
/// so the fractional part carries the satellite UPD difference. The
/// integer N_W is eliminated by working with `frac = mean - round(mean)`
/// wrapped to [-0.5, 0.5); the rank defect is fixed by sum(u) = 0.
/// Observations whose wrapped fraction sits within 1e-3 of +-0.5 are
/// skipped (unresolvable half-cycle ambiguity between two hypotheses).
pub fn solve_network_upd(
    per_base_means: &[HashMap<DoubleDiffKey, f64>],
) -> NetworkUpdSolution {
    // Observations: (sat, ref, wrapped fraction) per DD pair.
    let mut obs: Vec<(DoubleDiffKey, u16, u16, f64)> = Vec::new();
    for means in per_base_means {
        for (k, m) in means {
            let frac = m - m.round();
            let wrapped = frac.rem_euclid(1.0);
            // Unresolvable half-cycle: skip rather than guess a side.
            if (wrapped - 0.5).abs() < 1e-3 {
                continue;
            }
            let signed = if wrapped > 0.5 { wrapped - 1.0 } else { wrapped };
            obs.push((*k, k.sat, k.ref_sat, signed));
        }
    }
    let mut sats: Vec<u16> = obs.iter().flat_map(|(_, s, r, _)| [*s, *r]).collect();
    sats.sort_unstable();
    sats.dedup();
    if sats.len() < 2 || obs.is_empty() {
        return NetworkUpdSolution::default();
    }

    // Eliminate the rank defect via substitution: u_last = -sum(u_free).
    // Each observation contributes c . u_free + k to the predicted
    // fraction, with c_i = d(i,sat) - d(i,ref), k = -(d(sat,last)-d(ref,last)).
    let last = *sats.last().expect("non-empty");
    let free: Vec<u16> = sats[..sats.len() - 1].to_vec();
    let fidx: HashMap<u16, usize> = free.iter().enumerate()
        .map(|(i, &v)| (v, i)).collect();
    let m = free.len();

    let mut ata = vec![vec![0.0_f64; m]; m];
    let mut atb = vec![0.0_f64; m];
    for (key, s, r, frac) in &obs {
        let mut coef = vec![0.0_f64; m];
        for &(station, sign) in &[(s, 1.0), (r, -1.0)] {
            match fidx.get(station) {
                Some(&ci) => coef[ci] += sign,
                // Eliminated satellite: u_last = -sum(u_free), so its
                // contribution lands on EVERY free coefficient.
                None => {
                    for ci in coef.iter_mut() {
                        *ci -= sign;
                    }
                }
            }
        }
        for i in 0..m {
            for j2 in 0..m {
                ata[i][j2] += coef[i] * coef[j2];
            }
            atb[i] += coef[i] * frac;
        }
        let _ = key;
    }

    // Solve the normal equations (Gaussian elimination + back-substitution).
    let mut sol_free = atb.clone();
    for col in 0..m {
        let piv = (col..m)
            .fold(col, |best, r| {
                if ata[r][col].abs() > ata[best][col].abs() { r } else { best }
            });
        if ata[piv][col].abs() < 1e-12 {
            return NetworkUpdSolution::default();
        }
        ata.swap(piv, col);
        sol_free.swap(piv, col);
        // Forward elimination only; the back-substitution pass below
        // resolves the unknowns top-down.
        let pivot_row: Vec<f64> = ata[col].clone();
        for rr in col + 1..m {
            let fct = ata[rr][col] / ata[col][col];
            if fct == 0.0 { continue; }
            for cc in col..m {
                ata[rr][cc] -= fct * pivot_row[cc];
            }
            sol_free[rr] -= fct * sol_free[col];
        }
    }
    for row in (0..m).rev() {
        let mut acc = sol_free[row];
        for cc in row + 1..m {
            acc -= ata[row][cc] * sol_free[cc];
        }
        sol_free[row] = acc / ata[row][row];
    }

    let mut sat_upd: HashMap<u16, f64> =
        free.iter().enumerate().map(|(i, &v)| (v, sol_free[i])).collect();
    sat_upd.insert(last, -sol_free.iter().sum::<f64>());
    // Exact post-fit residuals with the solved UPDs.
    let mut residuals = Vec::with_capacity(obs.len());
    let mut sq = 0.0_f64;
    for (key, s, r, frac) in &obs {
        let pred = sat_upd[s] - sat_upd[r];
        let mut e = frac - pred;
        while e >= 0.5 { e -= 1.0; }
        while e < -0.5 { e += 1.0; }
        sq += e * e;
        residuals.push((*key, e));
    }
    let residual_rms = if obs.is_empty() { 0.0 } else { (sq / obs.len() as f64).sqrt() };

    NetworkUpdSolution { sat_upd, residuals, residual_rms }
}

#[cfg(test)]
mod upd_tests {
    use super::*;

    #[test]
    fn test_network_upd_solver_recovers_satellite_offsets() {
        // True satellite UPDs (sum zero), observed from two bases over
        // pairs against ref G21, with small deterministic noise.
        let true_u: HashMap<u16, f64> = [
            (21u16, 0.0_f64),
            (2u16, 0.30),
            (5u16, -0.20),
            (12u16, -0.10),
        ]
        .into_iter()
        .collect();
        let mut per_base: Vec<HashMap<DoubleDiffKey, f64>> = Vec::new();
        for b in 0..2u32 {
            let mut means = HashMap::new();
            for (&s, &u) in true_u.iter() {
                if s == 21 {
                    continue;
                }
                let key = DoubleDiffKey { constellation_id: 0, sat: s, ref_sat: 21, freq_band: 1 };
                let noise = if b == 0 { 0.01 } else { -0.01 };
                // N_W chosen so the arc mean = N_W + frac difference.
                let n_w = 7.0 + s as f64;
                means.insert(key, n_w + u - true_u[&21] + noise);
            }
            per_base.push(means);
        }

        let sol = solve_network_upd(&per_base);
        assert!(sol.residual_rms < 0.03, "rms {}", sol.residual_rms);
        for (&s, &u) in true_u.iter() {
            let est = sol.sat_upd.get(&s).copied().unwrap_or(f64::NAN);
            assert!(
                (est - u).abs() < 0.05,
                "sat {s}: estimated {est:.3} vs true {u:.3}"
            );
        }
    }
}
