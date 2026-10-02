//! Tests for the pure helpers inside the DD formation module.
//!
//! A child module of `formation` so the private differencing helpers are
//! reachable.

use gneiss_core::obs::{ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};

use super::*;

fn obs(band: u8, kind: ObsType, value: f64) -> Observation {
    Observation {
        code: ObsCode { obs_type: kind, signal: SignalCode { freq_band: band, attribute: 'C' } },
        value,
        lock_time: None,
        lli: Some(0),
    }
}

fn sat(prn: u8, cons: Constellation, records: Vec<Observation>) -> SatObs {
    SatObs { sat: SatelliteId { constellation: cons, prn }, observations: records }
}

fn gps(prn: u8, records: Vec<Observation>) -> SatObs {
    sat(prn, Constellation::Gps, records)
}

// ---------------------------------------------------------------------------
// canonical_bands_for_constellation
// ---------------------------------------------------------------------------

#[test]
fn canonical_bands_are_the_documented_per_constellation_sets() {
    assert_eq!(canonical_bands_for_constellation(Constellation::Gps), &[1, 2, 5]);
    // QZSS shares the GPS CDMA bands.
    assert_eq!(canonical_bands_for_constellation(Constellation::Qzss), &[1, 2, 5]);
    assert_eq!(canonical_bands_for_constellation(Constellation::Galileo), &[1, 5, 7]);
    assert_eq!(canonical_bands_for_constellation(Constellation::Beidou), &[1, 6, 7]);
    // Everything else (GLONASS FDMA, SBAS, Navic) gets the L1/L2 pair.
    for c in [Constellation::Glonass, Constellation::Sbas, Constellation::Navic] {
        assert_eq!(canonical_bands_for_constellation(c), &[1, 2], "{c:?}");
    }
}

#[test]
fn every_canonical_band_set_is_inside_the_rtk_band_universe() {
    for c in [
        Constellation::Gps,
        Constellation::Glonass,
        Constellation::Galileo,
        Constellation::Beidou,
        Constellation::Qzss,
    ] {
        for b in canonical_bands_for_constellation(c) {
            assert!(ALL_RTK_BANDS.contains(b), "{c:?} band {b} is not in ALL_RTK_BANDS");
        }
    }
}

// ---------------------------------------------------------------------------
// filter_constellation_sats
// ---------------------------------------------------------------------------

#[test]
fn the_gps_bucket_absorbs_qzss() {
    let info = vec![
        (SatelliteId { constellation: Constellation::Gps, prn: 6 }, Vector3::new(1.0, 0.0, 0.0)),
        (SatelliteId { constellation: Constellation::Qzss, prn: 2 }, Vector3::new(0.0, 1.0, 0.0)),
        (SatelliteId { constellation: Constellation::Galileo, prn: 11 }, Vector3::new(0.0, 0.0, 1.0)),
    ];
    let gps_bucket = filter_constellation_sats(&info, Constellation::Gps as u8);
    assert_eq!(gps_bucket.len(), 2, "GPS + QZSS share one double-difference bucket");
    let gal = filter_constellation_sats(&info, Constellation::Galileo as u8);
    assert_eq!(gal.len(), 1);
    assert_eq!(gal[0].0.prn, 11);
}

// ---------------------------------------------------------------------------
// raw_dd_observables
// ---------------------------------------------------------------------------

#[test]
fn double_differencing_is_rover_minus_base_after_satellite_minus_reference() {
    // Code (band 1):
    //   rover  sat - rover  ref = 25.0 - 24.0 =  1.0
    //   base   sat - base   ref = 25.5 - 24.3 =  1.2
    //   dd                        = 1.0 - 1.2    = -0.2 m
    // Phase (band 1):
    //   rover  1.000e8 - 9.90e7 = 1.0e6
    //   base   1.050e8 - 9.97e7 = 5.3e6
    //   dd                        = 1.0e6 - 5.3e6 = -4.3e6 cycles
    let rs = gps(2, vec![obs(1, ObsType::Pseudorange, 25.0), obs(1, ObsType::CarrierPhase, 1.000e8)]);
    let rr = gps(1, vec![obs(1, ObsType::Pseudorange, 24.0), obs(1, ObsType::CarrierPhase, 9.90e7)]);
    let bs = gps(2, vec![obs(1, ObsType::Pseudorange, 25.5), obs(1, ObsType::CarrierPhase, 1.050e8)]);
    let br = gps(1, vec![obs(1, ObsType::Pseudorange, 24.3), obs(1, ObsType::CarrierPhase, 9.97e7)]);
    let (pr, cp) = raw_dd_observables(&rs, &bs, &rr, &br, 1).expect("all four present");
    assert!((pr - (-0.2)).abs() < 1e-12, "dd_pr = {pr}");
    assert!((cp.unwrap() - (-4.3e6)).abs() < 1.0, "dd_cp = {:?}", cp);
}

#[test]
fn a_missing_code_observation_drops_the_pair_entirely() {
    let with = gps(2, vec![obs(1, ObsType::Pseudorange, 25.0)]);
    let without = gps(2, vec![obs(1, ObsType::CarrierPhase, 1.0e8)]);
    let ref_sat = gps(1, vec![obs(1, ObsType::Pseudorange, 24.0)]);
    assert!(raw_dd_observables(&with, &with, &ref_sat, &ref_sat, 1).is_some());
    assert!(raw_dd_observables(&without, &without, &ref_sat, &ref_sat, 1).is_none());
}

#[test]
fn missing_carrier_phase_keeps_the_code_pair_with_no_phase() {
    // A code-only pair is still usable: the phase component must be None, not
    // zero, so downstream gates can tell "no phase" from "phase equals zero".
    let rs = gps(2, vec![obs(1, ObsType::Pseudorange, 25.0)]);
    let rr = gps(1, vec![obs(1, ObsType::Pseudorange, 24.0)]);
    let (pr, cp) = raw_dd_observables(&rs, &rs, &rr, &rr, 1).expect("code present");
    assert!((pr - 0.0).abs() < 1e-12);
    assert!(cp.is_none());
}

#[test]
fn a_band_without_observations_yields_nothing() {
    let rs = gps(2, vec![obs(1, ObsType::Pseudorange, 25.0)]);
    assert!(raw_dd_observables(&rs, &rs, &rs, &rs, 2).is_none());
}

// ---------------------------------------------------------------------------
// select_secondary_phase_band
// ---------------------------------------------------------------------------

#[test]
fn the_secondary_band_follows_the_constellation_ladder() {
    let mk_gps = |b: u8| gps(2, vec![obs(b, ObsType::CarrierPhase, 1.0e8)]);
    assert_eq!(select_secondary_phase_band(Constellation::Gps, &mk_gps(2), &mk_gps(2)), Some(2));
    assert_eq!(select_secondary_phase_band(Constellation::Gps, &mk_gps(5), &mk_gps(5)), Some(5));
    assert_eq!(select_secondary_phase_band(Constellation::Gps, &mk_gps(7), &mk_gps(7)), None);

    let mk_bds = |b: u8| {
        gneiss_core::obs::SatObs {
            sat: SatelliteId { constellation: Constellation::Beidou, prn: 6 },
            observations: vec![obs(b, ObsType::CarrierPhase, 1.0e8)],
        }
    };
    // BeiDou secondary priority: 7 (B2I), then 6 (B3I). Never 2.
    assert_eq!(select_secondary_phase_band(Constellation::Beidou, &mk_bds(7), &mk_bds(7)), Some(7));
    assert_eq!(select_secondary_phase_band(Constellation::Beidou, &mk_bds(6), &mk_bds(6)), Some(6));
    assert_eq!(select_secondary_phase_band(Constellation::Beidou, &mk_bds(2), &mk_bds(2)), None);
}

#[test]
fn a_secondary_band_seen_by_only_one_receiver_is_refused() {
    let both = gps(2, vec![obs(2, ObsType::CarrierPhase, 1.0e8)]);
    let none = gps(2, Vec::new());
    // L2 phase exists on the rover only: it cannot form a double difference.
    assert_eq!(select_secondary_phase_band(Constellation::Gps, &both, &none), None);
}

// ---------------------------------------------------------------------------
// check_wl_slip
// ---------------------------------------------------------------------------

#[test]
fn a_wide_lane_slip_on_any_band_invalidates_the_pair() {
    let f = |pairs: &[(u8, bool)]| -> bool {
        let m: HashMap<u8, bool> = pairs.iter().copied().collect();
        GnssRtkIekf::check_wl_slip(&m)
    };
    assert!(!f(&[(1, false)]), "clean L1");
    assert!(!f(&[(1, false), (2, false), (5, false)]), "clean pair");
    assert!(f(&[(1, true)]), "L1 slip");
    assert!(f(&[(1, false), (2, true)]), "L2 slip widens the wide lane");
    assert!(f(&[(1, false), (5, true)]), "E5a slip widens the wide lane");
    assert!(f(&[(2, true)]), "slip on a band with no L1 record still counts");
}

// ---------------------------------------------------------------------------
// geom_dd_tropo
// ---------------------------------------------------------------------------

#[test]
fn a_zero_length_double_difference_is_exactly_zero() {
    // Rover sitting on the base with sat == ref degenerates the geometry:
    // (rs - rr) and base_dd both vanish.
    let p = Vector3::new(1.0e7, 2.0e7, 1.5e7);
    assert!(geom_dd_tropo(p, p, Vector3::zeros(), Vector3::zeros()).abs() < 1e-9);
}

#[test]
fn the_geometry_term_follows_the_range_derivative_with_sign() {
    // Moving the rover 1000 m along the line of sight changes
    // |sat - rover| by -1000 m, so geom_dd must drop by ~1000 m.
    // The troposphere term is identical at both points for this geometry,
    // so the delta isolates the range term.
    // A real station, a satellite 20 000 km along +x and a reference 1e9 m
    // along +y (far enough that its range derivative over the step is under
    // 1 mm). Advancing the rover 100 km toward the satellite must drop the
    // geometry term by 100 km: the range enters with a MINUS sign. The
    // troposphere term is unchanged to within a few centimetres over that
    // step, so a 2 m tolerance isolates the range contribution.
    let site = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    let sat = site + Vector3::new(2.0e7, 0.0, 0.0);
    let refs = site + Vector3::new(0.0, 1.0e9, 0.0);
    let base = site + Vector3::new(300.0, -100.0, 0.0);
    let a = geom_dd_tropo(sat, refs, base, site);
    let b = geom_dd_tropo(sat, refs, base, site + Vector3::new(1.0e5, 0.0, 0.0));
    assert!(
        (b - a + 1.0e5).abs() < 2.0,
        "geometry must track the range rate, delta = {} (want -100000)",
        b - a
    );
}

#[test]
fn beidou_phase_widelane_records_valid_wavelength_and_finite_cycles() {
    let t0 = gneiss_core::time::GpsTime::new(2000, 0.0);
    let mut iekf = GnssRtkIekf::new(Vector3::zeros(), t0, 1.0);
    let sat_id = SatelliteId { constellation: Constellation::Beidou, prn: 6 };
    let sat_pos = Vector3::new(2.0e7, 0.0, 0.0);
    let ref_pos = Vector3::new(0.0, 2.0e7, 0.0);
    let base_pos = Vector3::new(100.0, 0.0, 0.0);

    for _ in 0..12 {
        iekf.update_phase_wl(sat_id, sat_pos, 1, ref_pos, base_pos, 100.0, 80.0, 7, 0, false);
    }
    let key = DoubleDiffKey {
        constellation_id: Constellation::Beidou as u8,
        sat: 6,
        ref_sat: 1,
        freq_band: 1,
    };
    let (mean, count) = iekf.pw_tracker.arc_means().get(&key).copied().expect("arc must be tracked");
    assert_eq!(count, 12);
    assert!(mean.is_finite(), "PW mean must be finite, got {mean}");
}