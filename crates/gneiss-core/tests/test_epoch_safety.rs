//! Integration and safety tests for typed epoch systems and scale alignment.

use gneiss_core::ephemeris::BeidouEphemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::{
    BdtEpoch, Epoch, EpochKey, GlonassScale, GpsEpoch, GpsScale, GpsTime, GstScale,
    TimeDelta, UtcScale, BDT_OFFSET_NANOS, BDT_OFFSET_SECONDS, NANOS_PER_SEC,
};
use std::collections::HashMap;

#[test]
fn test_subtraction_across_week_rollover_exact() {
    let t_prev: GpsEpoch = Epoch::new(2200, 604_790.0);
    let t_curr: GpsEpoch = Epoch::new(2201, 10.0);

    let dt = t_curr - t_prev;
    assert_eq!(dt, TimeDelta::from_seconds(20.0));
    assert_eq!(dt.as_nanos(), 20 * (NANOS_PER_SEC as i64));
    assert_eq!(dt.as_millis(), 20_000);

    let dt_rev = t_prev - t_curr;
    assert_eq!(dt_rev, TimeDelta::from_seconds(-20.0));
    assert_eq!(dt_rev.as_nanos(), -20 * (NANOS_PER_SEC as i64));
}

#[test]
fn test_addition_and_subtraction_of_timedelta() {
    let t0: GpsEpoch = Epoch::new(2200, 604_790.0);
    let delta = TimeDelta::from_seconds(25.0);

    let t1 = t0 + delta;
    assert_eq!(t1.week(), 2201);
    assert_eq!(t1.tow_nanos(), 15 * NANOS_PER_SEC);

    let t_back = t1 - delta;
    assert_eq!(t_back, t0);
}

/// Clock-bias probe tolerance. `af1` is inflated 7 orders above any real
/// broadcast value, so a 14 s GPST/BDT epoch mix-up moves the clock bias by
/// `af1 * 14 = 1.4e-3` s. This tolerance sits 6 orders below that error and
/// 6 orders above the ~1e-14 s of f64 round-off in the epoch subtraction.
const CLOCK_BIAS_TOL_S: f64 = 1e-9;

/// Builds a BeiDou B1I ephemeris probe for clock-bias testing.
///
/// Two elements are deliberately synthetic:
/// - `m0 = 0` makes the fixed-point eccentric-anomaly solve return exactly
///   `0.0`, so the relativistic term `F * e * sqrt_a * sin(E)` in
///   `calc_keplerian` vanishes identically and the expected bias at TOC is
///   exactly `af0 - tgd1`, to the last bit.
/// - `af1 = -1e-4` s/s is 7 orders above a real broadcast drift. With a
///   realistic `af1 ~ -6.8e-12` the same 14 s mix-up moves the clock by only
///   ~1e-10 s, which a 1e-10 tolerance would wave through — that is why the
///   original version of this test could not detect the error it named.
///
/// All other elements are realistic broadcast values and are unconstrained by
/// the clock-bias assertion.
fn beidou_clock_probe(sat: SatelliteId, toe: GpsTime, toc: GpsTime) -> BeidouEphemeris {
    BeidouEphemeris {
        sat,
        toe,
        toc,
        af0: -2.71548e-4,
        af1: -1.0e-4,
        af2: 0.0,
        aode: 70,
        crs: 70.0,
        delta_n: 4.32e-9,
        m0: 0.0,
        cuc: -1.76e-7,
        e: 1.93e-3,
        cus: 4.60e-6,
        sqrt_a: 5153.5396,
        cic: 1.35e-7,
        omega0: -2.968,
        cis: -3.91e-8,
        i0: 0.9798,
        crc: 298.0,
        omega: -1.061,
        omega_dot: -8.07e-9,
        idot: 4.57e-11,
        tgd1: 1.2e-9,
        tgd2: 0.0,
        aodc: 70,
    }
}

/// Asserts the clock bias evaluated in GPST equals the BDT-frame expectation.
#[track_caller]
fn assert_bias_matches_bdt_toc(eph: &BeidouEphemeris, t_eval_gpst: GpsTime) {
    let (_, _, clk_err, _) = eph.position(t_eval_gpst);
    let expected = eph.af0 - eph.tgd1;
    assert!(
        (clk_err - expected).abs() < CLOCK_BIAS_TOL_S,
        "clock bias at TOC should be af0 - tgd1 = {expected:e}, got {clk_err:e} \
         (delta {:e} s); a 14 s GPST/BDT mix-up would show ~1.4e-3 s",
        (clk_err - expected).abs(),
    );
}

#[test]
fn test_bdt_offset_constant_matches_icd() {
    // Pinned independently of the implementation so the canonical constant
    // cannot silently drift from the ICD value the probe tests assert against.
    assert_eq!(BDT_OFFSET_SECONDS, 14.0);
    assert_eq!(BDT_OFFSET_NANOS, 14 * (NANOS_PER_SEC as i64));
}

#[test]
fn test_beidou_ephemeris_clock_bias_at_toc() {
    // PRN 1 is a BeiDou GEO satellite (prn <= 5), exercising the GEO branch.
    let toc_bdt = GpsTime::new(2105, 424_800.0);
    let eph = beidou_clock_probe(SatelliteId::new(Constellation::Beidou, 1), toc_bdt, toc_bdt);

    // toe/toc are BDT; callers pass GPST, which runs 14 s ahead.
    assert_bias_matches_bdt_toc(&eph, toc_bdt + 14.0);
}

#[test]
fn test_beidou_meo_clock_bias_at_toc() {
    // PRN 6 is neither <= 5 nor >= 59, exercising the non-GEO branch, which
    // takes a different Earth-rotation path than the GEO branch above.
    let toc_bdt = GpsTime::new(2105, 424_800.0);
    let eph = beidou_clock_probe(SatelliteId::new(Constellation::Beidou, 6), toc_bdt, toc_bdt);

    assert_bias_matches_bdt_toc(&eph, toc_bdt + 14.0);
}

#[test]
fn test_beidou_clock_bias_at_toc_across_week_boundary() {
    // TOC sits 5 s before the week rolls over, so the GPST->BDT subtraction has
    // to carry the week number down by one, not just the time of week.
    let toc_bdt = GpsTime::new(2105, 604_795.0);
    let eph = beidou_clock_probe(SatelliteId::new(Constellation::Beidou, 1), toc_bdt, toc_bdt);

    let t_eval_gpst = toc_bdt + 14.0;
    assert_eq!(t_eval_gpst.week, 2106, "probe must actually cross the week boundary");
    assert!((t_eval_gpst.tow - 9.0).abs() < 1e-9, "expected tow 9.0, got {}", t_eval_gpst.tow);

    assert_bias_matches_bdt_toc(&eph, t_eval_gpst);
}

#[test]
fn test_bdt_gpst_scale_offset_exact() {
    let bdt_epoch: BdtEpoch = Epoch::new(2100, 100_000.0);
    let gps_epoch = bdt_epoch.to_gpst();

    assert_eq!(gps_epoch.week(), 2100);
    assert_eq!(gps_epoch.tow_nanos(), 100_000 * NANOS_PER_SEC + (BDT_OFFSET_NANOS as u64));

    let bdt_back = gps_epoch.to_bdt();
    assert_eq!(bdt_back, bdt_epoch);
}

#[test]
fn test_gst_gpst_scale_exact() {
    let gst: Epoch<GstScale> = Epoch::new(2100, 50.0);
    let gps = gst.to_gpst();
    assert_eq!(gps.week(), gst.week());
    assert_eq!(gps.tow_nanos(), gst.tow_nanos());

    let gst_back = gps.to_gst();
    assert_eq!(gst_back, gst);
}

#[test]
fn test_glonass_gpst_scale_conversion() {
    let glo: Epoch<GlonassScale> = Epoch::new(2100, 10_800.0);
    let gps = glo.to_gpst(18);
    // GPST = GLONASST - 10800 + 18 = 18.0s
    assert_eq!(gps.tow_nanos(), 18 * NANOS_PER_SEC);

    let glo_back = gps.to_glonass(18);
    assert_eq!(glo_back, glo);
}

#[test]
fn test_utc_gpst_scale_conversion() {
    let utc: Epoch<UtcScale> = Epoch::new(2100, 100.0);
    let gps = utc.to_gpst(18);
    // GPST = UTC + 18s
    assert_eq!(gps.tow_nanos(), 118 * NANOS_PER_SEC);

    let utc_back = gps.to_utc(18);
    assert_eq!(utc_back, utc);
}

#[test]
fn test_epoch_key_hash_and_matching() {
    let mut map: HashMap<EpochKey<GpsScale>, &'static str> = HashMap::new();

    let e1: GpsEpoch = Epoch::new(2200, 604_799.5);
    let e2: GpsEpoch = Epoch::new(2201, 0.5);

    map.insert(e1.to_key(), "prev_week");
    map.insert(e2.to_key(), "next_week");

    assert_eq!(map.get(&e1.to_key()), Some(&"prev_week"));
    assert_eq!(map.get(&e2.to_key()), Some(&"next_week"));

    let diff_ms = e2.to_key().continuous_ms() - e1.to_key().continuous_ms();
    assert_eq!(diff_ms, 1000);
}

#[test]
fn test_tolerance_based_is_within() {
    let t_rover: GpsEpoch = Epoch::new(2200, 100.000);
    let t_base_close: GpsEpoch = Epoch::new(2200, 100.001); // 1 ms drift
    let t_base_far: GpsEpoch = Epoch::new(2200, 100.050); // 50 ms drift

    let tol = TimeDelta::from_millis(5);
    assert!(t_rover.is_within(t_base_close, tol));
    assert!(!t_rover.is_within(t_base_far, tol));
}

#[test]
fn test_legacy_gpstime_interop() {
    let legacy = GpsTime::new(2200, 123.456);
    let typed: GpsEpoch = legacy.into();
    assert_eq!(typed.week(), 2200);

    let roundtrip: GpsTime = typed.into();
    assert_eq!(roundtrip.week, legacy.week);
    assert!((roundtrip.tow - legacy.tow).abs() < 1e-8);
}
