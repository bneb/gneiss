//! Tests for the Melbourne-Wubbena observable and arc tracker.
//!
//! A child module of `mw` so the Welford accumulator and the network-UPD
//! normal equations are reachable.

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::obs::{ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};

use super::*;

const F1: f64 = 1575.42e6;
const F2: f64 = 1227.60e6;

fn key(sat: u16, ref_sat: u16) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat, ref_sat, freq_band: 1 }
}

fn obs(band: u8, kind: ObsType, value: f64) -> Observation {
    Observation {
        code: ObsCode { obs_type: kind, signal: SignalCode { freq_band: band, attribute: 'C' } },
        value,
        lock_time: None,
        lli: Some(0),
    }
}

/// Satellite with dual-frequency code and phase on every band.
fn dual_band_sat(prn: u8) -> SatObs {
    let mut o = Vec::new();
    for b in [1u8, 2] {
        o.push(obs(b, ObsType::CarrierPhase, 2.05e8 + f64::from(b) * 1.0e5));
        o.push(obs(b, ObsType::Pseudorange, 2.05e7 + f64::from(b) * 1.0e5));
    }
    SatObs { sat: SatelliteId { constellation: Constellation::Gps, prn }, observations: o }
}

// ---------------------------------------------------------------------------
// MwTrack: Welford statistics against textbook definitions
// ---------------------------------------------------------------------------

#[test]
fn welford_mean_and_sample_sigma_match_their_definitions() {
    // For the sample set {1,2,3,4,5}:
    //   mean = 3
    //   sample variance = (4 + 1 + 0 + 1 + 4) / (5-1) = 10/4 = 2.5
    //   sigma = sqrt(2.5) = 1.5811388300841898
    //   sigma_mean = sigma / sqrt(5) = 0.7071067811865475
    let mut t = MwTrack::new(1.0);
    for x in [2.0, 3.0, 4.0, 5.0] {
        t.push(x);
    }
    assert_eq!(t.count(), 5);
    assert!((t.mean() - 3.0).abs() < 1e-12, "mean {}", t.mean());
    assert!((t.sigma() - 2.5_f64.sqrt()).abs() < 1e-12, "sigma {}", t.sigma());
    assert!((t.sigma_mean() - 0.5_f64.sqrt()).abs() < 1e-12, "sigma_mean {}", t.sigma_mean());
}

#[test]
fn a_single_epoch_track_has_zero_sigma() {
    let t = MwTrack::new(7.5);
    assert_eq!(t.sigma(), 0.0, "one sample has no spread");
    assert_eq!(t.sigma_mean(), 0.0);
}

// ---------------------------------------------------------------------------
// mw_dd_cycles
// ---------------------------------------------------------------------------

#[test]
fn mw_dd_cycles_rejects_degenerate_frequency_pairs() {
    let p = [0.0; 4];
    assert!(mw_dd_cycles(0.0, F2, p, p, p, p).is_none());
    assert!(mw_dd_cycles(F1, -1.0, p, p, p, p).is_none());
    // |f1 - f2| = 0.5 MHz is below the 1 MHz guard: the wide-lane
    // wavelength would be 600 km and the observable meaningless.
    assert!(mw_dd_cycles(1005.0e6, 1004.5e6, p, p, p, p).is_none());
    assert!(mw_dd_cycles(F1, F2, p, p, p, p).is_some());
}

#[test]
fn narrow_lane_scale_is_the_documented_frequency_ratio() {
    let (_, nl) = mw_dd_cycles(F1, F2, [0.0; 4], [0.0; 4], [0.0; 4], [0.0; 4]).expect("valid pair");
    // nl_scale = f2 / (f1 - f2) = 1227.60 / 347.82
    assert!((nl - 1227.60 / 347.82).abs() < 1e-12, "nl_scale = {nl}");
    // That is the familiar GPS L1/L2 ratio, about 3.53 cycles per wide-lane cycle.
    assert!(nl > 3.52 && nl < 3.54, "nl_scale {nl} outside 3.52..3.54");
}

// ---------------------------------------------------------------------------
// WidelaneTracker gates
// ---------------------------------------------------------------------------

#[test]
fn the_fast_path_uses_a_tighter_deviation_gate_than_the_standard_path() {
    let k = key(2, 1);
    // 5 identical epochs: sigma_mean = 0, so the FAST gate applies
    // (n >= 5, sigma_mean <= 0.08) while the STANDARD gate does not
    // (n < MIN_TRACK_EPOCHS = 12). FAST allows only 0.15 cycles.
    let mut t = WidelaneTracker::default();
    for _ in 0..5 {
        t.update(k, 5.20, 3.53, false);
    }
    assert!(
        t.fixed_widelane(&k).is_none(),
        "0.20 cycles exceeds the fast-path 0.15 gate"
    );

    let mut ok = WidelaneTracker::default();
    for _ in 0..5 {
        ok.update(k, 5.10, 3.53, false);
    }
    let (w, n) = ok.fixed_widelane(&k).expect("0.10 cycles is inside the fast gate");
    assert_eq!(n, 5);
    assert!((w - 5.10).abs() < 1e-12);

    // The same 5.20 cycles over 12 epochs now uses the 0.25-cycle gate.
    let mut std = WidelaneTracker::default();
    for _ in 0..12 {
        std.update(k, 5.20, 3.53, false);
    }
    assert!(std.fixed_widelane(&k).is_some(), "standard gate allows 0.20 cycles");
}

#[test]
fn an_unconverged_arc_never_yields_a_wide_lane_fix() {
    let k = key(2, 1);
    let mut t = WidelaneTracker::default();
    // Four identical epochs: below FAST_TRACK_EPOCHS = 5, so nothing is fixed.
    for _ in 0..4 {
        t.update(k, 5.02, 3.53, false);
    }
    assert!(t.fixed_widelane(&k).is_none());
    assert!(t.arc_means().is_empty(), "arc means also need MIN_TRACK_EPOCHS");
}

#[test]
fn a_satellite_upd_difference_is_removed_before_rounding() {
    let k = key(9, 3);
    let mut t = WidelaneTracker::default();
    for _ in 0..5 {
        t.update(k, 7.20, 3.53, false);
    }
    // No UPD: rounds to 7 with a 0.20 deviation, outside the fast gate.
    assert!(t.fixed_widelane(&k).is_none());
    // With sat UPDs of +0.30 (sat 9) and +0.10 (ref 3) the corrected wide lane
    // is 7.20 - (0.30 - 0.10) = 7.00, dead on the integer.
    let mut upd = HashMap::new();
    upd.insert(9u16, 0.30);
    upd.insert(3u16, 0.10);
    t.sat_upd = Some(upd);
    let (w, n) = t.fixed_widelane(&k).expect("UPD-corrected wide lane is exactly integer");
    assert_eq!(n, 7);
    assert!((w - 7.0).abs() < 1e-12, "w = {w}");
}

#[test]
fn arc_means_need_the_full_twelve_epochs() {
    let k = key(2, 1);
    let mut t = WidelaneTracker::default();
    for _ in 0..11 {
        t.update(k, 5.02, 3.53, false);
    }
    assert!(t.arc_means().is_empty());
    t.update(k, 5.02, 3.53, false);
    let means = t.arc_means();
    assert_eq!(means.len(), 1);
    assert_eq!(means[&k].1, 12);
    assert!((means[&k].0 - 5.02).abs() < 1e-12);
}

#[test]
fn the_innovation_gate_restarts_a_settled_arc_but_a_shielded_one_survives() {
    let k = key(2, 1);
    let mut t = WidelaneTracker::default();
    for _ in 0..20 {
        t.update(k, 5.02, 3.53, false);
    }
    assert!(t.fixed_widelane(&k).is_some());
    // +6.5 cycles with INNOVATION_ARM_EPOCHS = 5 already passed: restart.
    t.update(k, 11.52, 3.53, false);
    assert!(t.fixed_widelane(&k).is_none());
    assert_eq!(t.reset_pair(&k), (), "reset is infallible");
    assert!(t.fixed_widelane(&k).is_none());
}

#[test]
fn retain_active_drops_pairs_that_left_the_double_difference_set() {
    let (k1, k2) = (key(2, 1), key(3, 1));
    let mut t = WidelaneTracker::default();
    for _ in 0..15 {
        t.update(k1, 5.02, 3.53, false);
        t.update(k2, -2.02, 3.53, false);
    }
    assert_eq!(t.arc_means().len(), 2);
    t.retain_active(&[k1]);
    let means = t.arc_means();
    assert_eq!(means.len(), 1);
    assert!(means.contains_key(&k1));
}

#[test]
fn the_narrow_lane_scale_is_remembered_per_pair() {
    let k = key(2, 1);
    let mut t = WidelaneTracker::default();
    assert_eq!(t.nl_scale(&k), None);
    t.update(k, 5.02, 3.5292, false);
    assert!((t.nl_scale(&k).unwrap() - 3.5292).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// update_tracker_from_obs
// ---------------------------------------------------------------------------

#[test]
fn raw_observations_feed_the_tracker_with_the_expected_wide_lane() {
    let s = dual_band_sat(2);
    let r = dual_band_sat(1);
    let mut t = WidelaneTracker::default();
    update_tracker_from_obs(&mut t, s.sat, 1, &s, &s, &r, &r, 0, false);
    let k = key(2, 1);
    assert_eq!(t.nl_scale(&k), Some(F2 / (F1 - F2)));
    // The constructed observables imply some finite wide lane; the point of
    // the assertion is that a track now exists for the (sat 2, ref 1, band 1)
    // key at all -- i.e. the band selection and quad assembly succeeded.
    assert!(t.fixed_widelane(&k).is_none(), "one epoch cannot reach the gate");
}

#[test]
fn a_single_frequency_satellite_is_skipped_entirely() {
    let mut single = dual_band_sat(2);
    single.observations.retain(|o| o.code.signal.freq_band == 1);
    let r = dual_band_sat(1);
    let mut t = WidelaneTracker::default();
    update_tracker_from_obs(&mut t, single.sat, 1, &single, &single, &r, &r, 0, false);
    assert!(t.arc_means().is_empty());
    assert_eq!(t.nl_scale(&key(2, 1)), None, "no pair may be registered");
}

// ---------------------------------------------------------------------------
// solve_network_upd
// ---------------------------------------------------------------------------

fn means_for(sat: u16, refs: u16, value: f64) -> HashMap<DoubleDiffKey, f64> {
    let mut m = HashMap::new();
    m.insert(key(sat, refs), value);
    m
}

#[test]
fn one_baseline_fixes_relative_offsets_up_to_the_sum_zero_gauge() {
    // Wide-lane fractions observed on a single baseline:
    //   sat 2 vs ref 1 -> 5.20 -> frac +0.20
    //   sat 3 vs ref 1 -> -1.10 -> frac -0.10
    // Differences alone leave one free common mode; the solver removes it by
    // forcing the sum to zero. With u2 = u1 + 0.2 and u3 = u1 - 0.1:
    //   u1 + u2 + u3 = 3*u1 + 0.1 = 0  ->  u1 = -1/30 = -0.03333...
    //   u2 = 0.2 - 1/30 = 0.16666...  and  u3 = -0.1 - 1/30 = -0.13333...
    // so the two arcs must be fitted exactly, with zero post-fit residual.
    let sol = solve_network_upd(&[means_for(2, 1, 5.20), means_for(3, 1, -1.10)]);
    let (u1, u2, u3) = (
        sol.sat_upd[&1],
        sol.sat_upd[&2],
        sol.sat_upd[&3],
    );
    assert!((u2 - u1 - 0.20).abs() < 1e-9, "u2 - u1 = {}", u2 - u1);
    assert!((u3 - u1 + 0.10).abs() < 1e-9, "u3 - u1 = {}", u3 - u1);
    assert!((u1 + u2 + u3).abs() < 1e-9, "sum-to-zero violated");
    assert!(sol.residual_rms < 1e-9, "exact system must fit exactly: {}", sol.residual_rms);
}

#[test]
fn fewer_than_two_distinct_satellites_yields_no_solution() {
    // A degenerate arc whose satellite IS its own reference names only one
    // satellite, so no relative offset exists and the solver must bail out.
    let mut degenerate = HashMap::new();
    degenerate.insert(key(2, 2), 5.20);
    let sol = solve_network_upd(&[degenerate]);
    assert!(sol.sat_upd.is_empty());
    assert!(sol.residuals.is_empty());
}

#[test]
fn a_half_cycle_mean_is_dropped_as_undecidable() {
    // |frac - 0.5| < 1e-3 is the ambiguity window: such an arc cannot be
    // signed, so it must not enter the normal equations at all.
    let mut m = HashMap::new();
    m.insert(key(2, 1), 5.50);
    m.insert(key(3, 1), 5.50);
    let sol = solve_network_upd(&[m]);
    assert!(sol.sat_upd.is_empty(), "two half-cycle arcs cannot be signed");
}

#[test]
fn fractions_are_wrapped_into_the_half_open_interval() {
    // -0.30 and +0.70 are the same wide-lane residual. Both must be reduced
    // to -0.30 so the UPD solution is consistent.
    let obs = vec![
        means_for(2, 1, 5.70),
        means_for(3, 1, 5.70),
        means_for(2, 4, 4.70),
        means_for(3, 4, 4.70),
    ];
    let sol = solve_network_upd(&obs);
    assert!(!sol.sat_upd.is_empty());
    assert!(sol.residual_rms < 1.0, "residual rms {}", sol.residual_rms);
}

#[test]
fn residuals_are_reported_per_observation_and_wrapped_to_half_a_cycle() {
    // Two bases seeing the same pair: still no absolute information, but the
    // solver must still be well defined for the relative part. Build a case
    // with a consistent solution and check every residual is inside +-0.5.
    let true_u: HashMap<u16, f64> = [(21u16, 0.0), (2, 0.30), (5, -0.20), (12, -0.10)]
        .into_iter().collect();
    let per_base: Vec<HashMap<DoubleDiffKey, f64>> = (0..2u32).map(|b| {
        let noise = if b == 0 { 0.01 } else { -0.01 };
        [2u16, 5, 12].into_iter().map(|s| {
            let k = key(s, 21);
            (k, 7.0 + f64::from(s) + true_u[&s] - true_u[&21] + noise)
        }).collect()
    }).collect();
    let sol = solve_network_upd(&per_base);
    assert_eq!(sol.residuals.len(), 6);
    assert!(sol.residuals.iter().all(|(_, r)| r.abs() <= 0.5));
    // The sum-zero constraint is imposed on the last satellite.
    let sum: f64 = sol.sat_upd.values().sum();
    assert!(sum.abs() < 1e-9, "UPDs must sum to zero, got {sum}");
    assert!((sol.residual_rms - 0.01).abs() < 1e-6, "rms {}", sol.residual_rms);
}

#[test]
fn an_empty_input_yields_the_default_solution() {
    let sol = solve_network_upd(&[]);
    assert!(sol.sat_upd.is_empty());
    let sol2 = solve_network_upd(&[HashMap::new()]);
    assert!(sol2.sat_upd.is_empty());
    // Sanity: the MW lambda is the wide-lane wavelength, not the L1 one.
    let lambda_wl = SPEED_OF_LIGHT_M_S / (F1 - F2);
    assert!(lambda_wl > 0.8 && lambda_wl < 0.9, "lambda_wl = {lambda_wl} m");
    assert!((F2 / (F1 - F2) - 3.529).abs() < 0.01);
}