//! Opt-in AR-quality gates (`GNEISS_AR_GATE=1`), both off by default.
//!
//! * **AR elevation mask** — integer resolution runs on a higher elevation
//!   cut-off than measurement processing (cf. RTKLIB `prcopt.elmaskar`,
//!   applied inside `ddmat()`): low-elevation pairs carry the largest
//!   multipath/tropospheric residuals and are the first to break the ratio
//!   test, so they are withheld from LAMBDA while still feeding the float
//!   filter. Implemented as a *temporary restricted view*: the live filter
//!   state keeps every float, so a pair that dips below the mask resumes AR
//!   with its accumulated ambiguity once it climbs back above.
//! * **Phase-code coherency bias init** — a freshly initialised DD
//!   ambiguity seeded from code-minus-phase inherits that pair's code
//!   multipath and iono divergence, so the filter sees a step in the phase
//!   innovation the epoch a satellite rises. RTKLIB `udbias()` counteracts
//!   this by shifting established biases by their mean raw-vs-state offset
//!   ("correct phase-bias offset to enssure phase-code coherency"); here the
//!   equivalent correction is applied to the new seed: subtract the median
//!   divergence of the pairs already tracked this epoch on the same band.

use nalgebra::{DMatrix, Vector3};

use super::state::{DoubleDiffKey, RtkState};
use super::update::DoubleDiffMeasurement;

/// Median of `xs` (mean of the two central values for even length).
/// Empty input yields 0.0 — callers treat "no samples" as a no-op offset.
pub(crate) fn median(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let mut sorted = xs.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mid = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

/// Median code-minus-phase divergence among this epoch's already-tracked
/// pairs on `band`. Divergences from other bands are ignored: the offset is
/// expressed in cycles, so it does not transfer across wavelengths.
pub(crate) fn coherence_offset(band: u8, divergences: &[(u8, f64)]) -> f64 {
    let same_band: Vec<f64> = divergences.iter()
        .filter(|(b, _)| *b == band)
        .map(|(_, d)| *d)
        .collect();
    median(&same_band)
}

/// Temporary AR-only view of `state` restricted to DD pairs whose *both*
/// members clear `mask_rad` elevation above the rover. Pairs without a
/// matching measurement this epoch fail open (kept) so bookkeeping gaps can
/// never silently shrink the LAMBDA input below the ungated behaviour.
/// Never mutates `state`.
pub(crate) fn elevation_filtered_view(
    state: &RtkState,
    meas: &[DoubleDiffMeasurement],
    mask_rad: f64,
) -> RtkState {
    let rx = state.pos_ecef;
    let llh = gneiss_core::coords::ecef_to_llh(rx);
    let keep: Vec<usize> = state.ambiguities.iter().enumerate()
        .filter(|(_, (k, _))| pair_clears_mask(meas, k, llh, rx, mask_rad))
        .map(|(i, _)| i)
        .collect();
    if keep.len() == state.ambiguities.len() {
        return state.clone();
    }
    restricted_state(state, &keep)
}

/// Elevation test for one DD pair: both the satellite and the reference
/// satellite must sit at or above `mask` as seen from the rover.
fn pair_clears_mask(
    meas: &[DoubleDiffMeasurement],
    key: &DoubleDiffKey,
    llh: Vector3<f64>,
    rx: Vector3<f64>,
    mask: f64,
) -> bool {
    match meas.iter().find(|m| m.key == *key) {
        None => true,
        Some(m) => {
            let el_sat = gneiss_core::coords::az_el(llh, rx, m.sat_pos).1;
            let el_ref = gneiss_core::coords::az_el(llh, rx, m.ref_pos).1;
            el_sat >= mask && el_ref >= mask
        }
    }
}

/// Copy of `state` holding only the ambiguities at `keep` (indices into the
/// ambiguity block), with the covariance reduced to the matching rows and
/// columns. Non-ambiguity columns always survive.
fn restricted_state(state: &RtkState, keep: &[usize]) -> RtkState {
    let off = state.amb_offset();
    let mut cols: Vec<usize> = (0..off).collect();
    cols.extend(keep.iter().map(|&i| off + i));
    let dim = cols.len();
    let mut cov = DMatrix::zeros(dim, dim);
    for (r, &old_r) in cols.iter().enumerate() {
        for (c, &old_c) in cols.iter().enumerate() {
            cov[(r, c)] = state.cov[(old_r, old_c)];
        }
    }
    let mut out = state.clone();
    out.ambiguities = keep.iter().map(|&i| state.ambiguities[i]).collect();
    out.cov = cov;
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;

    const ROVER: Vector3<f64> = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);

    fn key(sat: u16, band: u8) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: band }
    }

    /// Satellite placed `range_m` from the rover along `dir` (unit-ish).
    fn meas_for(k: DoubleDiffKey, dir: Vector3<f64>, range_m: f64) -> DoubleDiffMeasurement {
        DoubleDiffMeasurement {
            key: k,
            dd_pr_m: 0.0,
            dd_cp_cycles: None,
            sat_pos: ROVER + dir.normalize() * range_m,
            ref_pos: ROVER + Vector3::y() * 2.6e7, // overhead reference
            base_pos: ROVER,
            lambda: 0.19,
            pr_var_m2: 1.0,
            cp_var_cycles2: 1.0,
            dm_wet_rov: 0.0,
            dgrad_n_rov: 0.0,
            dgrad_e_rov: 0.0,
            tide_dd_m: 0.0,
                dd_clk_m: 0.0,
            dd_pcv_m: 0.0,
        }
    }

    fn elevation(pos: Vector3<f64>) -> f64 {
        let llh = gneiss_core::coords::ecef_to_llh(ROVER);
        gneiss_core::coords::az_el(llh, ROVER, pos).1
    }

    #[test]
    fn median_handles_odd_even_unsorted_and_empty() {
        assert_eq!(median(&[]), 0.0);
        assert_eq!(median(&[2.0]), 2.0);
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(median(&[-5.0, 5.0]), 0.0);
    }

    #[test]
    fn coherence_offset_is_median_within_same_band_only() {
        let divs = [(1u8, 10.0), (2u8, 99.0), (1u8, 12.0), (1u8, 11.0), (2u8, -50.0)];
        assert_eq!(coherence_offset(1, &divs), 11.0);
        assert_eq!(coherence_offset(2, &divs), 24.5);
        assert_eq!(coherence_offset(5, &divs), 0.0, "unknown band -> no offset");
        assert_eq!(coherence_offset(1, &[]), 0.0, "no prior pairs -> no offset");
    }

    #[test]
    fn elevation_filter_keeps_only_pairs_with_both_members_above_mask() {
        let mut state = RtkState::new(ROVER, GpsTime::new(2200, 300000.0));
        state.ensure_ambiguity(key(2, 1), 10.5, 4.0); // high pair
        state.ensure_ambiguity(key(3, 1), 20.5, 9.0); // low pair
        let up = ROVER.normalize();
        // Direction orthogonal to local up -> near-horizon satellite.
        let east = Vector3::new(-up.z, 0.0, up.x).normalize();
        let meas = vec![
            meas_for(key(2, 1), up, 2.4e7),
            meas_for(key(3, 1), east, 2.4e7),
        ];
        let hi_el = elevation(meas[0].sat_pos);
        let lo_el = elevation(meas[1].sat_pos);
        assert!(hi_el > 0.5, "sanity: zenith satellite el={hi_el}");
        assert!(lo_el.abs() < 0.05, "sanity: horizon satellite el={lo_el}");

        let view = elevation_filtered_view(&state, &meas, 0.2618); // 15 deg
        assert_eq!(view.ambiguities.len(), 1, "low pair excluded");
        assert_eq!(view.ambiguities[0].0, key(2, 1));
        assert_eq!(view.get_amb_idx(&key(3, 1)), None);

        // Covariance submatrix follows the kept ambiguity.
        let off = view.amb_offset();
        assert_eq!(view.cov.nrows(), off + 1);
        assert!((view.cov[(off, off)] - 4.0).abs() < 1e-12, "kept var preserved");

        // Non-destructive: the live state is untouched.
        assert_eq!(state.ambiguities.len(), 2);
        assert_eq!(state.cov.nrows(), state.amb_offset() + 2);
    }

    #[test]
    fn elevation_filter_fails_open_without_measurement_geometry() {
        let mut state = RtkState::new(ROVER, GpsTime::new(2200, 300000.0));
        state.ensure_ambiguity(key(9, 1), 1.0, 1.0); // no meas entry for this key
        let up = ROVER.normalize();
        let meas = vec![meas_for(key(2, 1), up, 2.4e7)];
        let view = elevation_filtered_view(&state, &meas, 1.0);
        assert_eq!(view.ambiguities.len(), 1, "unassessable pair is kept");
    }

    #[test]
    fn elevation_filter_below_every_mask_returns_equivalent_state() {
        let mut state = RtkState::new(ROVER, GpsTime::new(2200, 300000.0));
        state.ensure_ambiguity(key(2, 1), 10.5, 4.0);
        state.ensure_ambiguity(key(3, 1), 20.5, 9.0);
        let up = ROVER.normalize();
        let meas = vec![
            meas_for(key(2, 1), up, 2.4e7),
            meas_for(key(3, 1), up, 2.45e7),
        ];
        let view = elevation_filtered_view(&state, &meas, 0.0);
        assert_eq!(view.ambiguities, state.ambiguities);
        assert_eq!(view.cov, state.cov);
    }
}
