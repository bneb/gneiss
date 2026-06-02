use core::ops::{Add, Sub};

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
    pub fn from_calendar(year: i32, month: i32, day: i32, hour: i32, minute: i32, sec: f64) -> Self {
        let mut y = year;
        let mut m = month;
        if m <= 2 {
            y -= 1;
            m += 12;
        }
        
        let d = day as f64 + hour as f64 / 24.0 + minute as f64 / 1440.0 + sec / 86400.0;
        
        let a = libm::floor(y as f64 / 100.0);
        let b = 2.0 - a + libm::floor(a / 4.0);
        let jd = libm::floor(365.25 * (y as f64 + 4716.0)) + libm::floor(30.6001 * (m as f64 + 1.0)) + d + b - 1524.5;
        
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
    #[must_use]
    pub fn normalize(mut self) -> Self {
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
}
