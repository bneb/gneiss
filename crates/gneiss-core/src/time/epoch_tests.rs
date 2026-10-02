//! Adversarial and property tests for the typed epoch / duration primitives.
//!
//! Split out of `epoch.rs` to keep that file under the 500-line budget. Every
//! expectation is derived from the definition of the GNSS week grid
//! (`W = 604800 s = 604800e9 ns`) or from an integer identity, never from the
//! expression under test.

#[cfg(test)]
mod tests {
    use super::super::*;

    // -----------------------------------------------------------------------
    // Property tests. Nothing below re-derives a value with the expression it
    // is testing; every expectation comes from the definition of the week
    // grid (W = 604800 s = 604800e9 ns) or from an integer identity.
    // -----------------------------------------------------------------------

    /// Cheap deterministic 64-bit LCG (Knuth MMIX constants). Avoids a
    /// dependency and is fully reproducible from the seed below.
    struct Lcg(u64);

    impl Lcg {
        fn next_i64(&mut self) -> i64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 1) as i64
        }
    }

    /// `normalize_week_nanos` must always land in `[0, W)` and must conserve
    /// the instant: `(week + n/W) == week' + n'/W`. Checked on the modulus
    /// because the week counter is a `u32`.
    #[test]
    fn normalize_week_nanos_always_returns_a_canonical_week_and_preserves_instant() {
        let w = NANOS_PER_WEEK as i128;
        let cases: [(u32, i128); 9] = [
            (0, 0),
            (0, -1),
            (0, 1),
            (7, -w),
            (7, -w - 1),
            (7, w),
            (7, w + 1),
            (u32::MAX, -w),
            (u32::MAX, -w - 1),
        ];
        let mut rng = Lcg(0x5eed_1234_abcd_0001);
        for (week, n) in cases {
            let (wn, nn) = normalize_week_nanos(week, n);
            assert!(nn < NANOS_PER_WEEK, "{week} {n} -> tow {nn} out of range");
            // The week counter is u32, so the instant is conserved modulo
            // 2^32 weeks (week 0 underflow wraps to u32::MAX, not to -1).
            let period = (1_i128 << 32) * w;
            let before = ((week as i128) * w + n).rem_euclid(period);
            let after = ((wn as i128) * w + (nn as i128)).rem_euclid(period);
            assert_eq!(after, before, "{week} {n} -> {wn}/{nn} moved the instant");
        }
        for _ in 0..5000 {
            let week = (rng.next_i64() % 4_000_000) as u32;
            let n = (rng.next_i64() as i128) * (rng.next_i64() as i128 % 7_000_000_000_000_000);
            let (wn, nn) = normalize_week_nanos(week, n);
            assert!(nn < NANOS_PER_WEEK, "tow {nn} out of range");
            let period = (1_i128 << 32) * w;
            assert_eq!(
                ((wn as i128) * w + nn as i128).rem_euclid(period),
                ((week as i128) * w + n).rem_euclid(period)
            );
        }
    }

    /// The week counter is `u32` and the documented convention is wrap-around:
    /// one second before week 0 is week `u32::MAX`, not a negative week. This
    /// is a derivation from the representation, not a captured output.
    #[test]
    fn week_zero_underflow_wraps_to_u32_max_with_a_canonical_tow() {
        let e: Epoch<GpsScale> = Epoch::new(0, -1.0);
        assert_eq!(e.week(), u32::MAX);
        assert_eq!(e.tow_nanos(), NANOS_PER_WEEK - NANOS_PER_SEC);
        // And the same instant expressed 1 week later.
        let f: Epoch<GpsScale> = Epoch::from_week_nanos(u32::MAX, NANOS_PER_WEEK);
        assert_eq!(f.week(), 0);
        assert_eq!(f.tow_nanos(), 0);
    }

    /// Every public constructor must produce a canonical epoch, including for
    /// time-of-week values that need carrying in either direction.
    #[test]
    fn every_constructor_yields_a_canonical_epoch() {
        let probes: [(u32, f64); 8] = [
            (0, 0.0),
            (0, -1e-9),
            (5, 604_800.0),
            (5, 604_800.0 + 1e-9),
            (5, -604_800.0),
            (5, -1_207_600.0),
            (2100, 0.000_000_001),
            (2100, 604_799.999_999_999),
        ];
        for (w, tow) in probes {
            let a: Epoch<GpsScale> = Epoch::new(w, tow);
            assert!(a.tow_nanos() < NANOS_PER_WEEK, "{w}/{tow} -> {}", a.tow_nanos());
            let b: Epoch<GpsScale> = Epoch::from_week_nanos(w, (tow * 1e9) as u64);
            assert!(b.tow_nanos() < NANOS_PER_WEEK, "{w}/{tow} -> {}", b.tow_nanos());
        }
    }

    /// The permanent GPST - BDT offset is exactly 14 s. Asserted in integer
    /// nanoseconds at a week boundary, where an off-by-one week would surface —
    /// the failure mode the earlier clock-bias probe could not see.
    #[test]
    fn bdt_offset_is_exactly_fourteen_seconds_across_the_week_boundary() {
        // BDT 6 s before the week rolls over -> GPST 8 s after the rollover.
        let bdt: Epoch<BdtScale> = Epoch::new(2105, 604_794.0);
        let gpst = bdt.to_gpst();
        assert_eq!(gpst.week(), 2106);
        assert_eq!(gpst.tow_nanos(), 8 * NANOS_PER_SEC);
        // ... and the offset must be exactly 14 s in nanoseconds. Subtraction
        // is same-scale only, so compare the two GPST epochs directly.
        let expected: Epoch<GpsScale> = Epoch::new(2106, 8.0);
        assert_eq!(gpst.to_bdt().to_gpst(), gpst);
        assert_eq!(expected.to_bdt().to_gpst(), expected);
        assert_eq!(BDT_OFFSET_NANOS, 14 * (NANOS_PER_SEC as i64));
        assert_eq!(BDT_OFFSET_SECONDS, 14.0);
        // Reverse direction at week 0 (wraps to u32::MAX, see test above).
        let zero: Epoch<GpsScale> = Epoch::new(0, 0.0);
        let back = zero.to_bdt();
        assert_eq!(back.week(), u32::MAX);
        assert_eq!(back.tow_nanos(), NANOS_PER_WEEK - BDT_OFFSET_NANOS as u64);
        assert_eq!(back.to_gpst(), zero);
    }

    /// Round-trip identity for every scale conversion, on a sweep that
    /// deliberately straddles week boundaries and offset sign changes.
    #[test]
    fn every_scale_conversion_round_trips_across_week_boundaries() {
        let probes: [(u32, f64); 7] = [
            (0, 0.0),
            (0, 13.0),  // GPST 0:0 -> BDT crosses into the previous week
            (0, 14.0),  // GPST exactly 14 s
            (5, 604_799.0),
            (5, 604_787.0), // GPST 604787 -> UTC 604805: carries a week
            (2100, 302_400.0),
            (2100, 604_799.999_999),
        ];
        for (w, s) in probes {
            let g: Epoch<GpsScale> = Epoch::new(w, s);
            assert_eq!(g.to_bdt().to_gpst(), g, "BDT round trip failed at {w}/{s}");
            assert_eq!(g.to_gst().to_gpst(), g, "GST round trip failed at {w}/{s}");
            for leap in [0_i32, 18, 37] {
                assert_eq!(g.to_utc(leap).to_gpst(leap), g, "UTC round trip at {w}/{s}");
                assert_eq!(g.to_glonass(leap).to_gpst(leap), g, "GLONASS round trip at {w}/{s}");
            }
        }
    }

    /// The three scale offsets are exact closed forms:
    ///   GPST = BDT + 14 s
    ///   GPST = GST  + 0 s
    ///   GPST = UTC + leap
    ///   GPST = GLONASST - 10800 s + leap
    #[test]
    fn scale_offsets_match_their_closed_forms() {
        let n = NANOS_PER_SEC as i64;
        for &(w, s) in &[(0u32, 43200.0f64), (7, 604_000.0), (2100, 1.0)] {
            let g: Epoch<GpsScale> = Epoch::new(w, s);
            assert_eq!(g.to_gst().to_gpst(), g);
            assert_eq!(g.to_gst().tow_nanos(), g.tow_nanos());
            for leap in [0_i32, 1, 18, 37] {
                // GPST = UTC + leap, i.e. going GPST -> UTC moves back `leap`.
                let utc = g.to_utc(leap);
                assert_eq!(
                    g.continuous_nanos() as i64 - utc.continuous_nanos() as i64,
                    (leap as i64) * n,
                    "UTC offset must be exactly `leap` seconds at {w}/{s}"
                );
                // GPST = GLONASST - 10800 + leap.
                let gl = g.to_glonass(leap);
                assert_eq!(
                    gl.continuous_nanos() as i64 - g.continuous_nanos() as i64,
                    10_800 * n - (leap as i64) * n,
                    "GLONASS offset at {w}/{s} with {leap} leap seconds"
                );
            }
        }
    }

    /// `Sub` is the workhorse for every epoch-difference computation, so it
    /// must be exact on integer nanoseconds and antisymmetric, including across
    /// a week rollover where a naive tow-only subtraction would be wrong.
    #[test]
    fn epoch_subtraction_is_exact_and_antisymmetric_across_rollovers() {
        let pairs: [((u32, f64), (u32, f64)); 6] = [
            ((2201, 10.0), (2200, 604_790.0)), // straddles the rollover
            ((2200, 0.0), (2199, 604_799.5)),
            ((2100, 604_799.999_999), (2100, 0.0)),
            ((0, 0.0), (0, 1e-9)),
            ((0, 1e-9), (0, 0.0)),
            ((5, 100.000_000_001), (5, 100.0)),
        ];
        for ((wa, sa), (wb, sb)) in pairs {
            let a: Epoch<GpsScale> = Epoch::new(wa, sa);
            let b: Epoch<GpsScale> = Epoch::new(wb, sb);
            let want = ((wa as i64 - wb as i64) * (NANOS_PER_WEEK as i64))
                + (a.tow_nanos() as i64 - b.tow_nanos() as i64);
            assert_eq!((a - b).as_nanos(), want, "{wa}/{sa} - {wb}/{sb}");
            assert_eq!((b - a).as_nanos(), -want, "antisymmetry");
            assert_eq!((a - a).as_nanos(), 0, "self difference");
        }
    }

    /// `continuous_ms`, `continuous_nanos` and `to_key` must be mutually
    /// consistent views of the same instant.
    #[test]
    fn continuous_views_agree_with_the_week_grid() {
        let a: Epoch<GpsScale> = Epoch::new(2100, 0.0);
        let b: Epoch<GpsScale> = Epoch::new(2100, 0.5);
        assert_eq!(
            a.continuous_nanos(),
            (2100_i128) * (NANOS_PER_WEEK as i128)
        );
        assert_eq!(a.continuous_ms(), 2100 * 604_800_000);
        assert_eq!(b.continuous_ms() - a.continuous_ms(), 500);
        assert_eq!(b.to_key().continuous_ms(), b.continuous_ms());
        assert_eq!(a.to_key().continuous_ms(), a.continuous_ms());
    }

    /// `TimeDelta` unit conversions. Sub-millisecond and sub-microsecond
    /// remainders must *truncate toward zero*, matching integer division,
    /// because `as_millis`/`as_micros` feed epoch keys.
    #[test]
    fn time_delta_unit_conversions_truncate_toward_zero() {
        let t = TimeDelta::from_nanos(1_500_000_500);
        assert_eq!(t.as_millis(), 1500);
        assert_eq!(t.as_micros(), 1_500_000);
        assert_eq!(t.as_seconds(), 1.500_000_5);
        let t = TimeDelta::from_nanos(-1_500_000_500);
        assert_eq!(t.as_millis(), -1500);
        assert_eq!(t.as_micros(), -1_500_000);
        assert_eq!(t.as_seconds(), -1.500_000_5);
        // Round-trip through the fractional constructor on exact values.
        for secs in [0.0_f64, 1.0, -1.0, 0.000_000_001, -0.000_000_001, 86_400.5] {
            assert_eq!(TimeDelta::from_seconds(secs).as_nanos(), (secs * 1e9) as i64);
        }
        assert_eq!(TimeDelta::from_nanos(0).abs(), TimeDelta::from_nanos(0));
        assert_eq!(TimeDelta::from_nanos(-5).abs(), TimeDelta::from_nanos(5));
    }

    /// `Epoch` addition must agree with `TimeDelta` subtraction and must be
    /// exact at a week boundary (the classic off-by-one-week failure).
    #[test]
    fn epoch_addition_is_exact_across_the_week_boundary() {
        // `from_week_nanos` normalizes, so an out-of-range tow is carried;
        // the canonical boundary epoch is week 2100 day 0, i.e. exactly one
        // week before the *next* week starts.
        let carried: Epoch<GpsScale> = Epoch::from_week_nanos(2100, NANOS_PER_WEEK);
        assert_eq!((carried.week(), carried.tow_nanos()), (2101, 0));
        let base: Epoch<GpsScale> = Epoch::from_raw_parts(2100, 0);
        let one = TimeDelta::from_nanos(1);
        let next = base + one;
        assert_eq!(next.week(), 2100);
        assert_eq!(next.tow_nanos(), 1);
        assert_eq!((next - base).as_nanos(), 1);
        assert_eq!(next - one, base);
        // Subtracting past the boundary must carry the week down.
        let back = base - one;
        assert_eq!(back.week(), 2099);
        assert_eq!(back.tow_nanos(), NANOS_PER_WEEK - 1);
        assert_eq!((base - back).as_nanos(), one.as_nanos());
    }
}
