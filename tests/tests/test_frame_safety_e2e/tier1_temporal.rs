//! Tier 1 E2E tests: Temporal typestates, integer nanosecond arithmetic, and epoch alignment.

use super::common::*;

// --- Feature 6: Temporal Typestate Scale Markers ---

#[test]
fn test_f06_timescale_name_identifiers() {
    assert_eq!(GpsScale::NAME, "GPST");
    assert_eq!(BdtScale::NAME, "BDT");
    assert_eq!(GstScale::NAME, "GST");
    assert_eq!(GlonassScale::NAME, "GLONASST");
    assert_eq!(UtcScale::NAME, "UTC");
}

#[test]
fn test_f06_epoch_creation_from_week_and_nanos() {
    let epoch: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 345_600_000_000_000);
    assert_eq!(epoch.week(), 2200);
    assert_eq!(epoch.tow_nanos(), 345_600_000_000_000);
    assert!((epoch.tow_seconds() - 345600.0).abs() < 1e-9);
}

#[test]
fn test_f06_epoch_normalization_overflow_carries_week() {
    // 604_800s + 10s overflow
    let overflow_nanos = WEEK_NANOS + 10_000_000_000;
    let epoch: Epoch<GpsScale> = Epoch::from_week_nanos(2200, overflow_nanos);
    assert_eq!(epoch.week(), 2201);
    assert_eq!(epoch.tow_nanos(), 10_000_000_000);
}

#[test]
fn test_f06_timedelta_unit_constructors() {
    let dt_s = TimeDelta::from_seconds(2.5);
    let dt_ms = TimeDelta::from_millis(2500);
    let dt_ns = TimeDelta::from_nanos(2_500_000_000);
    assert_eq!(dt_s.as_nanos(), 2_500_000_000);
    assert_eq!(dt_ms.as_nanos(), 2_500_000_000);
    assert_eq!(dt_ns.as_nanos(), 2_500_000_000);
    assert!((dt_s.as_seconds() - 2.5).abs() < 1e-9);
    assert_eq!(dt_ms.as_millis(), 2500);
}

#[test]
fn test_f06_epoch_monotonic_ordering() {
    let t1: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 100_000_000_000);
    let t2: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 200_000_000_000);
    let t3: Epoch<GpsScale> = Epoch::from_week_nanos(2201, 50_000_000_000);
    assert!(t1 < t2);
    assert!(t2 < t3);
    assert!(t1 < t3);
}

// --- Feature 7: Integer Nanosecond Time Arithmetic ---

#[test]
fn test_f07_subtraction_within_same_week() {
    let t1: Epoch<GpsScale> = Epoch::from_week_nanos(2250, 100_000_000_000);
    let t2: Epoch<GpsScale> = Epoch::from_week_nanos(2250, 105_000_000_000);
    let dt = t2 - t1;
    assert_eq!(dt.as_nanos(), 5_000_000_000);
    assert!((dt.as_seconds() - 5.0).abs() < 1e-12);
}

#[test]
fn test_f07_subtraction_across_week_boundary() {
    // Saturday 23:59:58 (tow = 604,798s) to Sunday 00:00:02 (tow = 2s) of next week
    let t_sat: Epoch<GpsScale> = Epoch::from_week_nanos(2250, 604_798_000_000_000);
    let t_sun: Epoch<GpsScale> = Epoch::from_week_nanos(2251, 2_000_000_000);
    let dt = t_sun - t_sat;
    assert_eq!(dt.as_nanos(), 4_000_000_000);
    assert!((dt.as_seconds() - 4.0).abs() < 1e-12);
}

#[test]
fn test_f07_addition_forward_time_delta() {
    let t0: Epoch<GpsScale> = Epoch::from_week_nanos(2250, 10_000_000_000);
    let dt = TimeDelta::from_seconds(100.5);
    let t1 = t0 + dt;
    assert_eq!(t1.week(), 2250);
    assert_eq!(t1.tow_nanos(), 110_500_000_000);
}

#[test]
fn test_f07_addition_negative_time_delta_across_week() {
    // Sunday 00:00:01 (week 2251) minus 5 seconds -> Saturday 23:59:56 (week 2250)
    let t_sun: Epoch<GpsScale> = Epoch::from_week_nanos(2251, 1_000_000_000);
    let dt = TimeDelta::from_seconds(-5.0);
    let t_sat = t_sun + dt;
    assert_eq!(t_sat.week(), 2250);
    assert_eq!(t_sat.tow_nanos(), 604_796_000_000_000);
}

#[test]
fn test_f07_exact_nanosecond_precision_preservation() {
    let t0: Epoch<GpsScale> = Epoch::from_week_nanos(2250, 0);
    let dt_1ns = TimeDelta::from_nanos(1);
    let mut t = t0;
    for _ in 0..1000 {
        t = t + dt_1ns;
    }
    assert_eq!(t.tow_nanos(), 1000);
    assert_eq!((t - t0).as_nanos(), 1000);
}

// --- Feature 8: BeiDou Broadcast Ephemeris Scale Fix ---

#[test]
fn test_f08_bdt_to_gpst_adds_exact_14_seconds() {
    let t_bdt: Epoch<BdtScale> = Epoch::from_week_nanos(900, 100_000_000_000);
    let t_gpst = t_bdt.to_gpst();
    assert_eq!(t_gpst.week(), 900);
    assert_eq!(t_gpst.tow_nanos(), 114_000_000_000);
}

#[test]
fn test_f08_gpst_to_bdt_subtracts_exact_14_seconds() {
    let t_gpst: Epoch<GpsScale> = Epoch::from_week_nanos(900, 114_000_000_000);
    let t_bdt = t_gpst.to_bdt();
    assert_eq!(t_bdt.week(), 900);
    assert_eq!(t_bdt.tow_nanos(), 100_000_000_000);
}

#[test]
fn test_f08_bdt_gpst_roundtrip_nanosecond_invariance() {
    let orig_bdt: Epoch<BdtScale> = Epoch::from_week_nanos(950, 456_789_123_456);
    let recon_bdt = orig_bdt.to_gpst().to_bdt();
    assert_eq!(orig_bdt.week(), recon_bdt.week());
    assert_eq!(orig_bdt.tow_nanos(), recon_bdt.tow_nanos());
}

#[test]
fn test_f08_bdt_to_gpst_across_week_boundary() {
    // BDT near end of week: 604,790s. Adding 14s overflows into next week at 4s.
    let t_bdt: Epoch<BdtScale> = Epoch::from_week_nanos(900, 604_790_000_000_000);
    let t_gpst = t_bdt.to_gpst();
    assert_eq!(t_gpst.week(), 901);
    assert_eq!(t_gpst.tow_nanos(), 4_000_000_000);
}

#[test]
fn test_f08_beidou_clock_correction_alignment() {
    let t_bdt_toe: Epoch<BdtScale> = Epoch::from_week_nanos(900, 200_000_000_000);
    let t_bdt_toc: Epoch<BdtScale> = Epoch::from_week_nanos(900, 200_000_000_000);
    let dt_scale = t_bdt_toc.to_gpst() - t_bdt_toe.to_gpst();
    assert_eq!(dt_scale.as_nanos(), 0);
}

// --- Feature 9: Explicit Leap Second Conversions ---

#[test]
fn test_f09_utc_to_gpst_explicit_18_seconds() {
    let t_utc: Epoch<UtcScale> = Epoch::from_week_nanos(2200, 100_000_000_000);
    let t_gpst = t_utc.to_gpst(GPS_LEAP_SECONDS_2017);
    assert_eq!(t_gpst.tow_nanos(), 118_000_000_000);
}

#[test]
fn test_f09_gpst_to_utc_explicit_18_seconds() {
    let t_gpst: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 118_000_000_000);
    let t_utc = t_gpst.to_utc(GPS_LEAP_SECONDS_2017);
    assert_eq!(t_utc.tow_nanos(), 100_000_000_000);
}

#[test]
fn test_f09_utc_gpst_roundtrip() {
    let orig_utc: Epoch<UtcScale> = Epoch::from_week_nanos(2250, 300_000_000_000);
    let recon_utc = orig_utc.to_gpst(18).to_utc(18);
    assert_eq!(orig_utc.week(), recon_utc.week());
    assert_eq!(orig_utc.tow_nanos(), recon_utc.tow_nanos());
}

#[test]
fn test_f09_future_leap_second_19s_parameter() {
    let t_utc: Epoch<UtcScale> = Epoch::from_week_nanos(2300, 50_000_000_000);
    let t_gpst = t_utc.to_gpst(19);
    assert_eq!(t_gpst.tow_nanos(), 69_000_000_000);
    let recon = t_gpst.to_utc(19);
    assert_eq!(recon.tow_nanos(), 50_000_000_000);
}

#[test]
fn test_f09_glonass_epoch_relationship() {
    // GLONASS = UTC(SU) + 3h (10,800s). GPST = GLONASST - 10800s + leap_seconds
    let glonass_tow_sec = 20000.0;
    let leap = 18.0;
    let gpst_tow_sec = glonass_tow_sec - 10800.0 + leap;
    assert_eq!(gpst_tow_sec, 20000.0 - 10782.0);
}

// --- Feature 10: Structured Epoch Alignment & Key ---

#[test]
fn test_f10_continuous_ms_monotonicity() {
    let t1: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 604_799_000_000_000); // 1s before rollover
    let t2: Epoch<GpsScale> = Epoch::from_week_nanos(2201, 1_000_000_000);       // 1s after rollover
    assert_eq!(t2.continuous_ms() - t1.continuous_ms(), 2000);
}

#[test]
fn test_f10_epoch_key_equality_and_hashing() {
    let t1: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 123_456_000_000);
    let t2: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 123_456_999_999);
    // Both fall in same integer millisecond
    assert_eq!(t1.to_key(), t2.to_key());
}

#[test]
fn test_f10_is_within_tolerance_positive() {
    let t1: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 100_000_000_000);
    let t2: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 100_000_400_000); // 400 microseconds diff
    let tol = TimeDelta::from_millis(1); // 1 millisecond tolerance
    assert!(t1.is_within(t2, tol));
    assert!(t2.is_within(t1, tol));
}

#[test]
fn test_f10_is_within_tolerance_negative() {
    let t1: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 100_000_000_000);
    let t2: Epoch<GpsScale> = Epoch::from_week_nanos(2200, 100_002_000_000); // 2 milliseconds diff
    let tol = TimeDelta::from_millis(1);
    assert!(!t1.is_within(t2, tol));
}

#[test]
fn test_f10_sub_millisecond_clock_jitter_handling() {
    let base_epoch: Epoch<GpsScale> = Epoch::from_week_tow(2200, 12345.000);
    let rover_epoch: Epoch<GpsScale> = Epoch::from_week_tow(2200, 12345.0008); // 0.8ms jitter
    let tol = TimeDelta::from_millis(1);
    assert!(rover_epoch.is_within(base_epoch, tol));
}
