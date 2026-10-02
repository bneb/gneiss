//! Time representations, scales, and high-precision epoch arithmetic.

use core::ops::{Add, Sub};

pub mod epoch;
pub mod scales;

pub use epoch::*;
pub use scales::*;

const SECONDS_IN_WEEK: f64 = 604800.0;

/// Represents GPS Time consisting of a continuous week number and time of week in seconds.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct GpsTime {
    pub week: u32,
    pub tow: f64,
}

impl GpsTime {
    /// Creates a new `GpsTime` and normalizes the time of week.
    pub fn new(week: u32, tow: f64) -> Self {
        GpsTime { week, tow }.normalize()
    }

    /// Creates a new `GpsTime` from UTC calendar date and time.
    pub fn from_calendar(
        year: i32,
        month: i32,
        day: i32,
        hour: i32,
        minute: i32,
        sec: f64,
    ) -> Self {
        let mut y = year;
        let mut m = month;
        if m <= 2 {
            y -= 1;
            m += 12;
        }

        let d = day as f64 + hour as f64 / 24.0 + minute as f64 / 1440.0 + sec / 86400.0;
        let a = libm::floor(y as f64 / 100.0);
        let b = 2.0 - a + libm::floor(a / 4.0);
        let jd = libm::floor(365.25 * (y as f64 + 4716.0))
            + libm::floor(30.6001 * (m as f64 + 1.0))
            + d
            + b
            - 1524.5;

        let diff = jd - 2444244.5; // JD of Jan 6 1980
        let week = libm::floor(diff / 7.0);
        let tow = (diff - week * 7.0) * 86400.0;

        Self::new(week as u32, tow)
    }

    /// Returns the fractional year (e.g. 2020.5) for geodetic epoch transformations.
    pub fn to_fractional_year(&self) -> f64 {
        let jd = 2444244.5 + (self.week as f64 * 7.0) + (self.tow / 86400.0);
        // J2000 epoch is JD 2451545.0 (Jan 1, 2000, 12:00)
        // 1 Julian year = 365.25 days
        2000.0 + (jd - 2451545.0) / 365.25
    }

    /// Normalizes the time so that `tow` is strictly in the range [0.0, 604800.0).
    ///
    /// Non-finite `tow` is returned unchanged: the loops below subtract a
    /// finite constant per iteration, so `tow = +-inf` would spin forever
    /// (`inf - 604800.0 == inf`). Callers that cannot tolerate a non-finite
    /// time must reject it themselves.
    #[must_use]
    pub fn normalize(mut self) -> Self {
        if !self.tow.is_finite() {
            return self;
        }
        while self.tow >= SECONDS_IN_WEEK {
            self.tow -= SECONDS_IN_WEEK;
            self.week = self.week.wrapping_add(1);
        }
        while self.tow < 0.0 {
            self.tow += SECONDS_IN_WEEK;
            self.week = self.week.wrapping_sub(1);
        }
        self
    }

    /// Convert into strongly typed `Epoch<GpsScale>`.
    pub fn to_epoch(&self) -> Epoch<GpsScale> {
        Epoch::from_week_tow(self.week, self.tow)
    }
}

impl Add<f64> for GpsTime {
    type Output = Self;

    fn add(self, seconds: f64) -> Self::Output {
        GpsTime::new(self.week, self.tow + seconds)
    }
}

impl Sub<f64> for GpsTime {
    type Output = Self;

    fn sub(self, seconds: f64) -> Self::Output {
        GpsTime::new(self.week, self.tow - seconds)
    }
}

impl Sub<GpsTime> for GpsTime {
    type Output = f64;

    /// Returns the difference in seconds between two GPS times.
    fn sub(self, other: GpsTime) -> Self::Output {
        let week_diff = self.week as i64 - other.week as i64;
        (week_diff as f64 * SECONDS_IN_WEEK) + self.tow - other.tow
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpstime_normalization() {
        let t1 = GpsTime::new(100, 604800.0);
        assert_eq!(t1.week, 101);
        assert_eq!(t1.tow, 0.0);

        let t2 = GpsTime::new(100, 604801.5);
        assert_eq!(t2.week, 101);
        assert_eq!(t2.tow, 1.5);

        let t3 = GpsTime::new(100, -1.0);
        assert_eq!(t3.week, 99);
        assert_eq!(t3.tow, 604799.0);

        let t4 = GpsTime::new(100, -604801.0);
        assert_eq!(t4.week, 98);
        assert_eq!(t4.tow, 604799.0);
    }

    #[test]
    fn test_gpstime_addition() {
        let t = GpsTime::new(100, 10.0);
        let t2 = t + 604800.0;
        assert_eq!(t2.week, 101);
        assert_eq!(t2.tow, 10.0);
    }

    #[test]
    fn test_gpstime_subtraction() {
        let t = GpsTime::new(100, 10.0);
        let t2 = t - 20.0;
        assert_eq!(t2.week, 99);
        assert_eq!(t2.tow, 604790.0);
    }

    #[test]
    fn test_gpstime_difference() {
        let t1 = GpsTime::new(101, 10.0);
        let t2 = GpsTime::new(100, 604790.0);
        let diff = t1 - t2;
        assert_eq!(diff, 20.0);
    }

    // -----------------------------------------------------------------------
    // Golden vectors. Every expected (week, tow) below is derived from the
    // GPS week/day grid, not read back from the implementation.
    //
    //   GPS week 0 began Sunday 1980-01-06 00:00:00 UTC.
    //   1980-01-06 -> 2000-01-06 = 20 yr = 20*365 + 5 leap days = 7305 days,
    //     where the leap days are Feb-29 of 1980, 1984, 1988, 1992, 1996
    //     (2000-02-29 falls after the 2000-01-06 mark). 7305 = 1043*7 + 4.
    //   Hence week 1043 day 4 = 2000-01-06 and week 1043 day 0 = 2000-01-02.
    //   2000-01-01 is 7300 days after week 0, and 7300 = 1042*7 + 6, so it is
    //     week 1042 day 6. J2000.0 adds 12 h:
    //     => week 1042, tow = 6*86400 + 43200 = 561600.
    // -----------------------------------------------------------------------

    /// (UTC civil date) -> (week, tow) from the grid derivation above.
    /// `y, mo, d` date, `h, mi, s` time of day, then the expected grid value.
    fn assert_grid(
        y: i32,
        mo: i32,
        d: i32,
        h: i32,
        mi: i32,
        s: f64,
        want: (u32, f64),
    ) {
        let t = GpsTime::from_calendar(y, mo, d, h, mi, s);
        assert_eq!(t.week, want.0, "week for {y}-{mo}-{d}");
        assert!(
            (t.tow - want.1).abs() < 1e-6,
            "tow for {y}-{mo}-{d}: got {}, want {}",
            t.tow,
            want.1
        );
    }

    #[test]
    fn calendar_golden_vectors_against_the_gps_week_grid() {
        assert_grid(1980, 1, 6, 0, 0, 0.0, (0, 0.0)); // week 0, day 0 (Sunday)
        assert_grid(1980, 1, 6, 12, 0, 0.0, (0, 43200.0)); // day 0, noon
        assert_grid(1980, 1, 13, 0, 0, 0.0, (1, 0.0)); // day 7 -> week 1
        assert_grid(2000, 1, 1, 12, 0, 0.0, (1042, 561600.0)); // J2000.0
        // 2000-01-01 00:00 is the same day, 12 h earlier.
        assert_grid(2000, 1, 1, 0, 0, 0.0, (1042, 518400.0));
        // 2024-01-01 was a Monday: week 2295 day 1.
        // 1980-01-06 -> 2024-01-06 = 44 yr = 44*365 + 11 leap days = 16071 d
        // (leap Feb-29 of 1980..2020); 16071 = 2295*7 + 6, so week 2295 day 0
        // is 2023-12-31 and day 1 is 2024-01-01.
        assert_grid(2024, 1, 1, 0, 0, 0.0, (2295, 86400.0));
        assert_grid(2024, 1, 8, 0, 0, 0.0, (2296, 86400.0)); // +7 days
        // Gregorian century rule: 2000 is a leap year (divisible by 400),
        // 1900 is not (divisible by 100 but not 400). 2000-02-29 exists and
        // 2000-03-01 must be exactly one day after 2000-02-28.
        assert_grid(2000, 2, 29, 0, 0, 0.0, (1051, 172800.0));
        assert_grid(2000, 3, 1, 0, 0, 0.0, (1051, 259200.0));
    }

    #[test]
    fn fractional_year_is_exact_at_j2000_and_linear_after() {
        let j2000 = GpsTime::from_calendar(2000, 1, 1, 12, 0, 0.0);
        assert!((j2000.to_fractional_year() - 2000.0).abs() < 1e-9);
        // 100 Gregorian years contain 25 leap days (Feb-29 of 2000..2096;
        // 2100 is NOT a leap year under the century rule), so
        // 100*365 + 25 = 36525 days = exactly one Julian century.
        let century = GpsTime::from_calendar(2100, 1, 1, 12, 0, 0.0);
        assert_eq!(century.week, 6260, "J2000 + 36525 d, week carry");
        assert_eq!(century.tow, 475200.0, "J2000 + 36525 d, tow");
        assert!(
            (century.to_fractional_year() - 2100.0).abs() < 1e-9,
            "one Julian century later must read 2100.0, got {}",
            century.to_fractional_year()
        );
        // Same thing expressed purely in the week grid: 36525 d = 5217 wk + 6 d.
        let via_grid = GpsTime::new(1042 + 5217, 561600.0 + 6.0 * 86400.0);
        assert_eq!(via_grid, century);
    }

    /// `normalize` subtracts a finite constant per iteration, so a non-finite
    /// `tow` never terminates. It must instead be returned unchanged.
    #[test]
    fn normalize_terminates_on_non_finite_time_of_week() {
        for bad in [f64::INFINITY, f64::NEG_INFINITY] {
            let t = GpsTime::new(100, bad);
            assert_eq!(t.week, 100);
            assert_eq!(t.tow, bad, "non-finite tow must pass through, not loop");
        }
        let nan = GpsTime::new(100, f64::NAN);
        assert!(nan.tow.is_nan(), "NaN tow must pass through");
    }

    /// Property: for any finite input the result is canonical, and the week
    /// arithmetic is exact for inputs within a week of a week boundary.
    #[test]
    fn normalization_is_canonical_and_preserves_the_instant() {
        // `new` replaces the time of week outright, so the reference instant is
        // (week, raw) in absolute terms, not base.tow + raw.
        let origin_week = 2000.0;
        let instant = origin_week * 604_800.0;
        let mut s = 1_i64;
        for _ in 0..2000 {
            s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            let raw = (s % 2_000_000_000) as f64 / 1000.0 - 1_000_000.0; // +/- 1e6 s
            let t = GpsTime::new(2000, raw);
            assert!(t.tow >= 0.0 && t.tow < 604_800.0, "tow {} out of range", t.tow);
            let got = t.week as f64 * 604_800.0 + t.tow;
            let want = instant + raw;
            assert!(
                (got - want).abs() < 1e-3,
                "week {} tow {} -> instant {got}, want {want}",
                t.week,
                t.tow
            );
        }
    }
}
