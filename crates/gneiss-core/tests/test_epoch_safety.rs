//! Integration and safety tests for typed epoch systems and scale alignment.

use gneiss_core::ephemeris::BeidouEphemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::{
    BdtEpoch, BdtScale, Epoch, EpochKey, GlonassScale, GpsEpoch, GpsScale, GpsTime, GstScale,
    TimeDelta, UtcScale, BDT_OFFSET_NANOS, NANOS_PER_SEC,
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

#[test]
fn test_beidou_ephemeris_clock_bias_at_toc() {
    let sat = SatelliteId::new(Constellation::Beidou, 1);
    let toc_bdt = GpsTime::new(2105, 424_800.0);
    let toe_bdt = GpsTime::new(2105, 424_800.0);
    let af0 = -2.71548e-4;
    let af1 = -6.82121e-12;
    let af2 = 0.0;
    let tgd1 = 1.2e-9;

    let eph = BeidouEphemeris {
        sat, toe: toe_bdt, toc: toc_bdt, af0, af1, af2,
        aode: 70, crs: 70.0, delta_n: 4.32e-9, m0: 2.456,
        cuc: -1.76e-7, e: 1.93e-3, cus: 4.60e-6, sqrt_a: 5153.5396,
        cic: 1.35e-7, omega0: -2.968, cis: -3.91e-8, i0: 0.9798,
        crc: 298.0, omega: -1.061, omega_dot: -8.07e-9, idot: 4.57e-11,
        tgd1, tgd2: 0.0, aodc: 70,
    };

    // Satellite position and clock evaluated at TOC:
    // When time in GPST is (toc_bdt + 14s), t_bdt inside position is toc_bdt.
    let t_eval_gpst = GpsTime::new(toc_bdt.week, toc_bdt.tow + 14.0);
    let (_, _, clk_err, _) = eph.position(t_eval_gpst);

    // At TOC, tc = t_bdt - toc_bdt = 0.0. Clock correction is af0 - tgd1 (neglecting relativity).
    let expected_clk = af0 - tgd1;
    assert!((clk_err - expected_clk).abs() < 1e-10, "Clock bias at TOC must match af0 - tgd1 without 14s error");
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
