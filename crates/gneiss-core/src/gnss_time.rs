//! Explicit GNSS time-system types.
//!
//! [`crate::time::GpsTime`] assumes every timestamp already is GPS system time
//! (GPST). That assumption silently breaks for the other constellations:
//!
//! | System  | Native time scale           | Same instant expressed in GPST          |
//! |---------|-----------------------------|------------------------------------------|
//! | GPS     | GPST                        | identical                                |
//! | GLONASS | UTC(SU) + 3 h               | GLONASST − 3 h + leap seconds (−10782 s) |
//! | BeiDou  | BDT (epoch 2006-01-01 UTC)  | BDT + 14 s                               |
//! | Galileo | GST (GGTO A₀ ≈ 0)           | GST + 0 s (nominal)                      |
//!
//! These offsets were previously applied ad hoc at each use site
//! (`toc + (18.0 - 10800.0)`, `tow - 14.0`, …), which is easy to forget or
//! mis-sign. This module makes the time system part of the *type*: values are
//! constructed in their native scale ([`GnssTime::new`],
//! [`GnssTime::from_gpst`]) and converted explicitly ([`GnssTime::to_gpst`]) —
//! one typed home for the offsets, onto which the remaining ad-hoc call sites
//! (rinex.rs, ephemeris.rs) migrate as a follow-up.
//!
//! # Sign convention
//!
//! [`TimeSystem::gpst_offset`] returns the seconds to **add** to a reading in
//! `self` to obtain GPST. Beware the classic GLONASS trap: because GLONASST is
//! UTC(SU)+3 h, a GLONASS clock reads ~10782 s *more* than GPST at the same
//! instant, so converting GLONASS → GPST **subtracts** 10782 s.
//!
//! # References
//!
//! * GLONASS ICD v5.1: GLONASST = UTC(SU) + 3 h; UTC(SU) itself steps by leap
//!   seconds, which is why only the GLONASS conversion depends on the leap count.
//! * BDS-SIS-ICD: BDT starts 2006-01-01 00:00:00 UTC, no leap seconds added
//!   since ⇒ GPST − BDT = 14 s, fixed by design.
//! * Galileo OS-SIS-ICD: GST shares the 1999-08-22 epoch with GPST week 1024;
//!   the GGTO (A₀) is broadcast and nominally 0.

use crate::time::GpsTime;

/// Seconds in one GNSS week.
const SECONDS_IN_WEEK: f64 = 604800.0;

/// The GNSS system time scales understood by the engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeSystem {
    /// GPS system time (GPST).
    Gps,
    /// GLONASS system time = UTC(SU) + 3 h.
    Glonass,
    /// BeiDou system time (BDT), epoch 2006-01-01 00:00:00 UTC.
    Bdt,
    /// Galileo system time (GST).
    Gst,
}

impl TimeSystem {
    /// Leap seconds between the GPST epoch (1980-01-06) and UTC, i.e. how far
    /// GPST leads UTC. 18 s since 2017-01-01; bump when IERS/BIPM announce the
    /// next leap second — only the GLONASS conversion depends on it (BDT/GST
    /// offsets are fixed relative to GPST by design and immune to future
    /// leap-second insertions).
    pub const GPS_LEAP_SECONDS: i32 = 18;

    /// The fixed UTC(SU)+3 h component of GLONASS time, in seconds.
    pub const GLONASS_UTC_OFFSET: f64 = 10800.0;

    /// GPST − BDT. Fixed by the BDT epoch definition; never changes.
    pub const BDT_TO_GPST: f64 = 14.0;

    /// Nominal GST → GPST offset (GGTO A₀). Aligned by design at the shared
    /// 1999-08-22 epoch; the real GGTO is broadcast and sub-microsecond.
    pub const GST_TO_GPST: f64 = 0.0;

    /// GLONASS → GPST offset under the current leap count:
    /// GPST = GLONASST − 3 h + leap seconds.
    pub const GLONASS_TO_GPST: f64 = Self::GPS_LEAP_SECONDS as f64 - Self::GLONASS_UTC_OFFSET;

    /// Seconds to add to a reading in `self` to obtain the same instant in GPST.
    pub const fn gpst_offset(self) -> f64 {
        match self {
            TimeSystem::Gps => 0.0,
            TimeSystem::Glonass => Self::GLONASS_TO_GPST,
            TimeSystem::Bdt => Self::BDT_TO_GPST,
            TimeSystem::Gst => Self::GST_TO_GPST,
        }
    }

    /// GLONASS → GPST offset under a future/hypothetical leap-second count.
    pub fn glonass_offset_with_leap_seconds(gps_leap_seconds: i32) -> f64 {
        f64::from(gps_leap_seconds) - Self::GLONASS_UTC_OFFSET
    }
}

/// A timestamp expressed in its native GNSS system time.
///
/// The `week`/`tow` pair means different physical instants depending on
/// [`GnssTime::sys`]; convert through [`GnssTime::to_gpst`] before mixing with
/// anything that assumes GPST.
#[derive(Debug, Clone, Copy)]
pub struct GnssTime {
    /// Which system time scale `week`/`tow` are expressed in.
    pub sys: TimeSystem,
    /// Week number in the system's own counting (wraps on u32 overflow).
    pub week: u32,
    /// Time of week in seconds, normalized to `[0, 604800)` by
    /// [`GnssTime::new`] and [`GnssTime::from_gpst`].
    pub tow: f64,
}

impl GnssTime {
    /// Constructs a timestamp in any system, carrying whole weeks out of
    /// `tow` so it lands in `[0, 604800)`. Non-finite `tow` is passed through
    /// untouched (garbage in, garbage out).
    pub fn new(sys: TimeSystem, week: u32, tow: f64) -> Self {
        let (week, tow) = normalized(week, tow);
        GnssTime { sys, week, tow }
    }

    /// THE way to obtain GPST seconds: applies this system's offset,
    /// normalizing across week boundaries. Every consumer that mixes a
    /// non-GPS timestamp with GPST data must go through here.
    pub fn to_gpst(&self) -> GpsTime {
        self.to_gpst_with_leap_seconds(TimeSystem::GPS_LEAP_SECONDS)
    }

    /// Like [`GnssTime::to_gpst`] but with an explicit GPS↔UTC leap-second
    /// count (integer seconds, e.g. 18 today, 19 after the next IERS leap),
    /// affecting the GLONASS conversion only.
    pub fn to_gpst_with_leap_seconds(&self, gps_leap_seconds: i32) -> GpsTime {
        let offset = match self.sys {
            TimeSystem::Glonass => TimeSystem::glonass_offset_with_leap_seconds(gps_leap_seconds),
            sys => sys.gpst_offset(),
        };
        let tow = self.tow + offset;
        if !tow.is_finite() {
            // Poisoned input (NaN/±inf, possibly via direct field access) must
            // NOT reach GpsTime's normalization loops, which would spin
            // forever; hand the poison through instead.
            return GpsTime { week: self.week, tow };
        }
        GpsTime::new(self.week, tow)
    }

    /// Converts a GPST timestamp *into* `sys`.
    pub fn from_gpst(sys: TimeSystem, t: GpsTime) -> Self {
        Self::from_gpst_with_leap_seconds(sys, t, TimeSystem::GPS_LEAP_SECONDS)
    }

    /// Leap-aware inverse of [`GnssTime::to_gpst_with_leap_seconds`]:
    /// `from_gpst_with_leap_seconds(sys, t, n).to_gpst_with_leap_seconds(n) == t`.
    pub fn from_gpst_with_leap_seconds(sys: TimeSystem, t: GpsTime, gps_leap_seconds: i32) -> Self {
        let offset = match sys {
            TimeSystem::Glonass => TimeSystem::glonass_offset_with_leap_seconds(gps_leap_seconds),
            sys => sys.gpst_offset(),
        };
        Self::new(sys, t.week, t.tow - offset)
    }
}

/// Largest |tow| (seconds) that `normalized` will wrap: ~2^53 weeks ≈ 285
/// million years. Beyond that the value is garbage data; wrapping it would be
/// meaningless, so it passes through untouched instead of being looped on.
const MAX_NORMALIZABLE_TOW: f64 = 9.0e15;

/// Carries whole weeks out of `tow`, leaving it in `[0, SECONDS_IN_WEEK)`
/// (or untouched when non-finite / beyond [`MAX_NORMALIZABLE_TOW`]). Week
/// arithmetic wraps mod 2^32, matching `GpsTime`.
fn normalized(week: u32, tow: f64) -> (u32, f64) {
    if !tow.is_finite() || tow.abs() > MAX_NORMALIZABLE_TOW {
        return (week, tow);
    }
    // IEEE fmod is exact and O(1) for any magnitude — a subtract-one-week loop
    // is how naive kernels hang on hostile input. The remainder keeps the sign
    // of `tow`; the discarded part is always an exact whole number of weeks.
    let mut r = tow % SECONDS_IN_WEEK;
    let mut w = i64::from(week) + ((tow - r) / SECONDS_IN_WEEK) as i64;
    if r < 0.0 {
        r += SECONDS_IN_WEEK;
        w -= 1;
    }
    if r >= SECONDS_IN_WEEK {
        // Rounding can park a tiny negative remainder exactly ON the boundary
        // (e.g. tow = -1e-12 adds up to exactly 604800.0). Roll forward to
        // preserve the invariant r ∈ [0, WEEK); the represented instant is
        // unchanged to f64 resolution (half-ulp there is ~6e-11 s).
        r -= SECONDS_IN_WEEK;
        w += 1;
    }
    let r = if r == 0.0 { 0.0 } else { r }; // canonicalize -0.0 from exact multiples
    (w as u32, r)
}

// ============================================================================
// Tests
//
// Written FIRST (TDD): they were developed and verified against deliberately
// broken stubs (offsets zeroed, normalization skipped) and only then the
// implementation above was allowed to satisfy them. They pin:
//   1. each system's to_gpst() against its ICD offset,
//   2. from_gpst/to_gpst round trips,
//   3. week-boundary carries near tow = 0 and 604800,
//   4. negative / multi-week tow normalization,
//   5. bit-exact reproduction of the still-live ad-hoc corrections in
//      gneiss-parsers/src/rinex.rs:756-757 and gneiss-core/src/ephemeris.rs
//      (migrating those call sites onto this module is a tracked follow-up).
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const W: u32 = 2312;
    const SYSTEMS: [TimeSystem; 4] = [
        TimeSystem::Gps,
        TimeSystem::Glonass,
        TimeSystem::Bdt,
        TimeSystem::Gst,
    ];

    #[test]
    fn icd_offsets_are_pinned() {
        assert_eq!(TimeSystem::GPS_LEAP_SECONDS, 18);
        assert_eq!(TimeSystem::GLONASS_UTC_OFFSET, 10800.0);
        assert_eq!(TimeSystem::GLONASS_TO_GPST, 18.0 - 10800.0);
        assert_eq!(TimeSystem::BDT_TO_GPST, 14.0);
        assert_eq!(TimeSystem::GST_TO_GPST, 0.0);
        assert_eq!(TimeSystem::Gps.gpst_offset(), 0.0);
        assert_eq!(TimeSystem::Glonass.gpst_offset(), -10782.0);
        assert_eq!(TimeSystem::Bdt.gpst_offset(), 14.0);
        assert_eq!(TimeSystem::Gst.gpst_offset(), 0.0);
    }

    #[test]
    fn each_system_to_gpst_applies_icd_offset() {
        assert_eq!(
            GnssTime::new(TimeSystem::Gps, W, 100.0).to_gpst(),
            GpsTime::new(W, 100.0)
        );
        // GLONASS second 0 of week 2300 is GPST 594018 s of week 2299.
        assert_eq!(
            GnssTime::new(TimeSystem::Glonass, 2300, 0.0).to_gpst(),
            GpsTime::new(2299, 594_018.0)
        );
        assert_eq!(
            GnssTime::new(TimeSystem::Bdt, 1100, 100.0).to_gpst(),
            GpsTime::new(1100, 114.0)
        );
        assert_eq!(
            GnssTime::new(TimeSystem::Gst, 1400, 12_345.5).to_gpst(),
            GpsTime::new(1400, 12_345.5)
        );
    }

    #[test]
    fn reproduces_rinex_glonass_toc_correction() {
        // Legacy expression (gneiss-parsers/src/rinex.rs:756):
        // toc_gpst = toc_gpst + (18.0 - 10800.0)
        for raw_tow in [345_600.0, 345_600.5, 10_779.0, 0.0] {
            let toc_parsed = GpsTime::new(W, raw_tow);
            let legacy = toc_parsed + (18.0 - 10800.0);
            let typed = GnssTime::new(TimeSystem::Glonass, toc_parsed.week, toc_parsed.tow).to_gpst();
            assert_eq!(typed, legacy, "GLONASS mismatch for raw tow {raw_tow}");
        }
    }

    #[test]
    fn reproduces_rinex_beidou_toc_correction() {
        // Legacy expression (gneiss-parsers/src/rinex.rs:757):
        // toc_gpst = toc_gpst + 14.0
        for raw_tow in [0.0, 100.5, 604_793.0] {
            let toc_parsed = GpsTime::new(W, raw_tow);
            let legacy = toc_parsed + 14.0;
            let typed = GnssTime::new(TimeSystem::Bdt, toc_parsed.week, toc_parsed.tow).to_gpst();
            assert_eq!(typed, legacy, "BeiDou mismatch for raw tow {raw_tow}");
        }
        // 604793 + 14 crosses the week: first seconds of the NEXT GPST week.
        assert_eq!(
            GnssTime::new(TimeSystem::Bdt, W, 604_793.0).to_gpst(),
            GpsTime::new(W + 1, 7.0)
        );
    }

    #[test]
    fn reproduces_ephemeris_beidou_conversion() {
        // Legacy expression (crates/gneiss-core/src/ephemeris.rs,
        // BeidouEphemeris::position and position_iono_free):
        // let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
        let t = GpsTime::new(1100, 3.0);
        let legacy = GpsTime::new(t.week, t.tow - 14.0);
        let typed = GnssTime::from_gpst(TimeSystem::Bdt, t);
        assert_eq!(
            (typed.week, typed.tow),
            (legacy.week, legacy.tow),
            "from_gpst(Bdt) must reproduce the ad-hoc keplerian input"
        );
        assert_eq!(typed.to_gpst(), t);
    }

    #[test]
    fn round_trip_all_systems_is_exact_for_integral_tows() {
        let tows = [0.0, 1.0, 43_200.0, 86_399.0, 604_799.0];
        for sys in SYSTEMS {
            for tow in tows {
                let origin = GpsTime::new(W, tow);
                let native = GnssTime::from_gpst(sys, origin);
                // Shifted systems must genuinely move the representation;
                // guards against a symmetric sign bug cancelling itself.
                if matches!(sys, TimeSystem::Glonass | TimeSystem::Bdt) {
                    assert!(
                        native.tow != origin.tow || native.week != origin.week,
                        "{sys:?}: offset vanished at tow {tow}"
                    );
                }
                assert_eq!(native.to_gpst(), origin, "{sys:?} tow {tow}");
            }
        }
    }

    #[test]
    fn round_trip_fractional_tows_within_nanosecond_tolerance() {
        // (tow - off) + off is not always bit-exact in f64 (ulp near 1e4 s is
        // ~2e-12 s), so require sub-nanosecond agreement instead.
        let tows = [0.1, 123_456.789, 604_799.999_999];
        for sys in SYSTEMS {
            for tow in tows {
                let origin = GpsTime::new(W, tow);
                let back = GnssTime::from_gpst(sys, origin).to_gpst();
                assert_eq!(back.week, origin.week, "{sys:?} week drifted");
                let drift = back - origin;
                assert!(
                    drift.abs() < 1e-9,
                    "{sys:?}: round-trip drift {drift} s at tow {tow}"
                );
            }
        }
    }

    #[test]
    fn leap_second_parameter_changes_only_glonass() {
        // Hypothetical next leap second: 19 s.
        let glo = GnssTime::new(TimeSystem::Glonass, 2300, 0.0);
        assert_eq!(
            glo.to_gpst_with_leap_seconds(19),
            GpsTime::new(2299, 594_019.0)
        );
        // Default path still uses the pinned 18 s.
        assert_eq!(glo.to_gpst(), GpsTime::new(2299, 594_018.0));
        // BDT/GST are leap-independent by design.
        let bdt = GnssTime::new(TimeSystem::Bdt, 1, 10.0);
        assert_eq!(bdt.to_gpst_with_leap_seconds(37), GpsTime::new(1, 24.0));
        let gst = GnssTime::new(TimeSystem::Gst, 9, 777.0);
        assert_eq!(gst.to_gpst_with_leap_seconds(99), GpsTime::new(9, 777.0));
    }

    #[test]
    fn leap_second_parameter_round_trips() {
        let origin = GpsTime::new(W, 456_789.0);
        let native = GnssTime::from_gpst_with_leap_seconds(TimeSystem::Glonass, origin, 19);
        assert_eq!(native.to_gpst_with_leap_seconds(19), origin);
        // Converting back with a stale leap count betrays exactly 1 s.
        assert_eq!(native.to_gpst() - origin, -1.0);
    }

    #[test]
    fn week_boundary_carries_across_conversions() {
        // BDT lags GPST: early-week GPST lands at the END of the previous BDT week.
        let bdt = GnssTime::from_gpst(TimeSystem::Bdt, GpsTime::new(1200, 5.0));
        assert_eq!(
            (bdt.sys, bdt.week, bdt.tow),
            (TimeSystem::Bdt, 1199, 604_791.0)
        );
        assert_eq!(bdt.to_gpst(), GpsTime::new(1200, 5.0));

        // GLONASS leads GPST: late-week GPST spills into the NEXT GLONASS week.
        let glo = GnssTime::from_gpst(TimeSystem::Glonass, GpsTime::new(1500, 604_000.0));
        assert_eq!((glo.week, glo.tow), (1501, 9_982.0));
        assert_eq!(glo.to_gpst(), GpsTime::new(1500, 604_000.0));

        // GST shifts nothing, even exactly on the boundary.
        let gst = GnssTime::from_gpst(TimeSystem::Gst, GpsTime::new(1500, 0.0));
        assert_eq!((gst.week, gst.tow), (1500, 0.0));
    }

    #[test]
    fn conversions_crossing_week_end_land_on_next_week() {
        // 604786 s BDT + 14 s == exactly the start of the next GPST week.
        assert_eq!(
            GnssTime::new(TimeSystem::Bdt, 5, 604_786.0).to_gpst(),
            GpsTime::new(6, 0.0)
        );
        // GLONASS 3000 s of week 1200 is still week 1199 in GPST.
        assert_eq!(
            GnssTime::new(TimeSystem::Glonass, 1200, 3_000.0).to_gpst(),
            GpsTime::new(1199, 597_018.0)
        );
    }

    #[test]
    fn negative_and_multi_week_tow_is_normalized_at_construction() {
        let g = GnssTime::new(TimeSystem::Bdt, 100, -1.0);
        assert_eq!((g.week, g.tow), (99, 604_799.0));

        let g = GnssTime::new(TimeSystem::Gps, 100, -2.0 * SECONDS_IN_WEEK - 1.0);
        assert_eq!((g.week, g.tow), (97, 604_799.0));

        let g = GnssTime::new(TimeSystem::Gst, 100, 2.0 * SECONDS_IN_WEEK);
        assert_eq!((g.week, g.tow), (102, 0.0));

        let g = GnssTime::new(TimeSystem::Glonass, 100, -SECONDS_IN_WEEK);
        assert_eq!((g.week, g.tow), (99, 0.0));

        // Exactly one full week rolls over; the boundary is inclusive.
        let g = GnssTime::new(TimeSystem::Bdt, 5, SECONDS_IN_WEEK);
        assert_eq!((g.week, g.tow), (6, 0.0));
    }

    #[test]
    fn non_finite_tow_is_passed_through_without_wrapping() {
        let n = GnssTime::new(TimeSystem::Gps, 7, f64::NAN);
        assert_eq!(n.week, 7);
        assert!(n.tow.is_nan());
        assert!(n.to_gpst().tow.is_nan());

        let inf = GnssTime::new(TimeSystem::Gps, 7, f64::INFINITY);
        assert_eq!(inf.week, 7);
        assert!(inf.tow.is_infinite());
    }

    #[test]
    fn from_gpst_records_requested_system() {
        assert_eq!(
            GnssTime::from_gpst(TimeSystem::Gps, GpsTime::new(1, 1.0)).sys,
            TimeSystem::Gps
        );
        assert_eq!(
            GnssTime::from_gpst(TimeSystem::Glonass, GpsTime::new(1, 1.0)).sys,
            TimeSystem::Glonass
        );
        assert_eq!(
            GnssTime::from_gpst(TimeSystem::Bdt, GpsTime::new(1, 1.0)).sys,
            TimeSystem::Bdt
        );
        assert_eq!(
            GnssTime::from_gpst(TimeSystem::Gst, GpsTime::new(1, 1.0)).sys,
            TimeSystem::Gst
        );
    }

    #[test]
    fn seam_between_weeks_stays_in_range() {
        // A tiny negative remainder rounds onto exactly 604800.0 during the
        // borrow; the kernel must roll forward so the documented invariant
        // tow ∈ [0, 604800) survives (the instant is unchanged to f64 ulp).
        let g = GnssTime::new(TimeSystem::Gps, 7, -1e-12);
        assert_eq!((g.week, g.tow), (7, 0.0));
        assert!((0.0..SECONDS_IN_WEEK).contains(&g.tow));

        let g = GnssTime::new(TimeSystem::Glonass, 100, -SECONDS_IN_WEEK);
        assert_eq!(g.tow, 0.0);
        assert!(!g.tow.is_sign_negative(), "negative zero must be canonicalized");
    }

    #[test]
    fn absurd_magnitudes_pass_through_without_hanging() {
        for tow in [-1.0e300, 1.0e300] {
            let g = GnssTime::new(TimeSystem::Gps, 42, tow);
            assert_eq!((g.week, g.tow), (42, tow), "garbage must not be wrapped");
        }
    }

    #[test]
    fn week_counter_wraps_at_u32_edges() {
        let g = GnssTime::new(TimeSystem::Bdt, u32::MAX, SECONDS_IN_WEEK);
        assert_eq!((g.week, g.tow), (0, 0.0));

        // BDT lags: early-week GPST borrows from week 0 into u32::MAX.
        let g = GnssTime::from_gpst(TimeSystem::Bdt, GpsTime::new(0, 5.0));
        assert_eq!((g.week, g.tow), (u32::MAX, 604_791.0));

        // GLONASS leads: late-week GPST spills from u32::MAX into week 0.
        let g = GnssTime::from_gpst(TimeSystem::Glonass, GpsTime::new(u32::MAX, 604_000.0));
        assert_eq!((g.week, g.tow), (0, 9_982.0));
    }

    #[test]
    fn non_finite_tow_bypasses_normalization_on_conversion() {
        // Direct field construction bypasses GnssTime::new; to_gpst must hand
        // the poison through instead of entering GpsTime's unbounded loops.
        let g = GnssTime {
            sys: TimeSystem::Gps,
            week: 3,
            tow: f64::NEG_INFINITY,
        };
        let out = g.to_gpst();
        assert_eq!(out.week, 3);
        assert!(out.tow.is_infinite());
    }
}
