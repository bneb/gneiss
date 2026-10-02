//! Typed epoch and duration primitives with integer nanosecond precision.

use core::marker::PhantomData;
use core::ops::{Add, AddAssign, Neg, Sub, SubAssign};

use super::scales::{BdtScale, GlonassScale, GpsScale, GstScale, TimeScale, UtcScale};
use super::GpsTime;

/// Exact number of seconds in one standard GNSS week.
pub const SECONDS_PER_WEEK: u64 = 604_800;

/// Number of nanoseconds in one standard second.
pub const NANOS_PER_SEC: u64 = 1_000_000_000;

/// Number of nanoseconds in one standard GNSS week: `604_800 * 10^9`.
pub const NANOS_PER_WEEK: u64 = SECONDS_PER_WEEK * NANOS_PER_SEC;

/// Permanent ICD offset between GPST and BDT in nanoseconds (GPST - BDT = 14s).
pub const BDT_OFFSET_NANOS: i64 = 14 * (NANOS_PER_SEC as i64);

/// Permanent ICD offset between GPST and BDT in fractional seconds (GPST - BDT = 14s).
/// Derived from [`BDT_OFFSET_NANOS`] so the integer and `f64` forms cannot drift apart.
pub const BDT_OFFSET_SECONDS: f64 = BDT_OFFSET_NANOS as f64 / NANOS_PER_SEC as f64;

/// GLONASS timezone offset from UTC(SU) in nanoseconds (3 hours = 10,800s).
pub const GLONASS_HOURS_NANOS: i64 = 3 * 3600 * (NANOS_PER_SEC as i64);

/// Standard leap seconds between GPST and UTC since 2017-01-01 (18 seconds).
pub const DEFAULT_LEAP_SECONDS: i32 = 18;

/// Normalizes signed nanoseconds and week into canonical `(week, tow_nanos)` with `tow_nanos < NANOS_PER_WEEK`.
#[inline]
pub(crate) fn normalize_week_nanos(week: u32, total_nanos: i128) -> (u32, u64) {
    let week_nanos = NANOS_PER_WEEK as i128;
    let extra_weeks = total_nanos.div_euclid(week_nanos);
    let rem_nanos = total_nanos.rem_euclid(week_nanos) as u64;
    let new_week = (week as i64).wrapping_add(extra_weeks as i64) as u32;
    (new_week, rem_nanos)
}

/// Signed duration represented with nanosecond precision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct TimeDelta {
    pub nanos: i64,
}

impl TimeDelta {
    /// Construct duration from integer nanoseconds.
    pub const fn from_nanos(nanos: i64) -> Self {
        Self { nanos }
    }

    /// Construct duration from integer microseconds.
    pub const fn from_micros(micros: i64) -> Self {
        Self { nanos: micros * 1_000 }
    }

    /// Construct duration from integer milliseconds.
    pub const fn from_millis(millis: i64) -> Self {
        Self { nanos: millis * 1_000_000 }
    }

    /// Construct duration from fractional seconds.
    pub fn from_seconds(secs: f64) -> Self {
        Self { nanos: libm::round(secs * 1e9) as i64 }
    }

    /// Return duration as integer nanoseconds.
    pub const fn as_nanos(&self) -> i64 {
        self.nanos
    }

    /// Return duration as integer microseconds.
    pub const fn as_micros(&self) -> i64 {
        self.nanos / 1_000
    }

    /// Return duration as integer milliseconds.
    pub const fn as_millis(&self) -> i64 {
        self.nanos / 1_000_000
    }

    /// Return duration as fractional seconds.
    pub fn as_seconds(&self) -> f64 {
        self.nanos as f64 / (NANOS_PER_SEC as f64)
    }

    /// Absolute value of duration.
    pub const fn abs(&self) -> Self {
        Self { nanos: self.nanos.abs() }
    }
}

impl Add for TimeDelta {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self { nanos: self.nanos + rhs.nanos }
    }
}

impl Sub for TimeDelta {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self { nanos: self.nanos - rhs.nanos }
    }
}

impl Neg for TimeDelta {
    type Output = Self;
    fn neg(self) -> Self {
        Self { nanos: -self.nanos }
    }
}

impl AddAssign for TimeDelta {
    fn add_assign(&mut self, rhs: Self) {
        self.nanos += rhs.nanos;
    }
}

impl SubAssign for TimeDelta {
    fn sub_assign(&mut self, rhs: Self) {
        self.nanos -= rhs.nanos;
    }
}

/// Strictly typed epoch within a specified time scale `S`.
///
/// Cross-scale subtraction fails to compile:
/// ```compile_fail
/// use gneiss_core::time::{Epoch, GpsScale, BdtScale};
/// let gps: Epoch<GpsScale> = Epoch::from_week_nanos(2000, 0);
/// let bdt: Epoch<BdtScale> = Epoch::from_week_nanos(2000, 0);
/// let _ = gps - bdt;
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Epoch<S: TimeScale> {
    week: u32,
    tow_nanos: u64,
    _scale: PhantomData<S>,
}

impl<S: TimeScale> Epoch<S> {
    /// Construct epoch from pre-normalized week and nanosecond time of week.
    pub const fn from_raw_parts(week: u32, tow_nanos: u64) -> Self {
        Self { week, tow_nanos, _scale: PhantomData }
    }

    /// Construct epoch and normalize into canonical `tow_nanos ∈ [0, 604_800_000_000_000)`.
    pub fn from_week_nanos(week: u32, tow_nanos: u64) -> Self {
        let (w, n) = normalize_week_nanos(week, tow_nanos as i128);
        Self { week: w, tow_nanos: n, _scale: PhantomData }
    }

    /// Construct epoch from week and fractional seconds of week.
    pub fn from_week_tow(week: u32, tow_sec: f64) -> Self {
        let nanos = libm::round(tow_sec * 1e9) as i128;
        let (w, n) = normalize_week_nanos(week, nanos);
        Self { week: w, tow_nanos: n, _scale: PhantomData }
    }

    /// Alias for `from_week_tow`.
    pub fn new(week: u32, tow_sec: f64) -> Self {
        Self::from_week_tow(week, tow_sec)
    }

    /// Continuous week number.
    pub const fn week(&self) -> u32 {
        self.week
    }

    /// Time of week in integer nanoseconds.
    pub const fn tow_nanos(&self) -> u64 {
        self.tow_nanos
    }

    /// Time of week in fractional seconds.
    pub fn tow_seconds(&self) -> f64 {
        self.tow_nanos as f64 / (NANOS_PER_SEC as f64)
    }

    /// Continuous millisecond timestamp since scale origin, immune to week rollover.
    pub fn continuous_ms(&self) -> u64 {
        (self.week as u64) * 604_800_000 + (self.tow_nanos / 1_000_000)
    }

    /// Continuous nanosecond timestamp since scale origin.
    pub fn continuous_nanos(&self) -> i128 {
        (self.week as i128) * (NANOS_PER_WEEK as i128) + (self.tow_nanos as i128)
    }

    /// Construct discrete collision-free hash/map key.
    pub fn to_key(&self) -> EpochKey<S> {
        EpochKey::new(self.continuous_ms())
    }

    /// Returns true if `other` is within the given `tolerance` window.
    pub fn is_within(&self, other: Self, tolerance: TimeDelta) -> bool {
        let diff = if *self >= other { *self - other } else { other - *self };
        diff <= tolerance
    }
}

impl<S: TimeScale> Sub<Epoch<S>> for Epoch<S> {
    type Output = TimeDelta;

    fn sub(self, other: Epoch<S>) -> TimeDelta {
        let week_diff = self.week as i64 - other.week as i64;
        let nanos_diff = self.tow_nanos as i64 - other.tow_nanos as i64;
        let total_nanos = week_diff * (NANOS_PER_WEEK as i64) + nanos_diff;
        TimeDelta::from_nanos(total_nanos)
    }
}

impl<S: TimeScale> Add<TimeDelta> for Epoch<S> {
    type Output = Self;

    fn add(self, delta: TimeDelta) -> Self {
        let total_nanos = self.tow_nanos as i128 + delta.nanos as i128;
        let (w, n) = normalize_week_nanos(self.week, total_nanos);
        Self { week: w, tow_nanos: n, _scale: PhantomData }
    }
}

impl<S: TimeScale> Sub<TimeDelta> for Epoch<S> {
    type Output = Self;

    fn sub(self, delta: TimeDelta) -> Self {
        self + (-delta)
    }
}

impl<S: TimeScale> AddAssign<TimeDelta> for Epoch<S> {
    fn add_assign(&mut self, delta: TimeDelta) {
        *self = *self + delta;
    }
}

impl<S: TimeScale> SubAssign<TimeDelta> for Epoch<S> {
    fn sub_assign(&mut self, delta: TimeDelta) {
        *self = *self - delta;
    }
}

/// Discrete continuous millisecond key for collision-free hash indexing across week rollovers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EpochKey<S: TimeScale> {
    continuous_ms: u64,
    _scale: PhantomData<S>,
}

impl<S: TimeScale> EpochKey<S> {
    /// Construct key from continuous milliseconds.
    pub const fn new(continuous_ms: u64) -> Self {
        Self { continuous_ms, _scale: PhantomData }
    }

    /// Access continuous millisecond counter.
    pub const fn continuous_ms(&self) -> u64 {
        self.continuous_ms
    }
}

// ---------------------------------------------------------------------------
// Explicit Cross-Scale Conversions
// ---------------------------------------------------------------------------

impl Epoch<BdtScale> {
    /// Convert BeiDou system time to GPS time (+14.0 seconds permanent offset).
    pub fn to_gpst(&self) -> Epoch<GpsScale> {
        let (w, n) = normalize_week_nanos(self.week, self.tow_nanos as i128 + BDT_OFFSET_NANOS as i128);
        Epoch { week: w, tow_nanos: n, _scale: PhantomData }
    }
}

impl Epoch<GpsScale> {
    /// Convert GPS time to BeiDou system time (-14.0 seconds permanent offset).
    pub fn to_bdt(&self) -> Epoch<BdtScale> {
        let (w, n) = normalize_week_nanos(self.week, self.tow_nanos as i128 - BDT_OFFSET_NANOS as i128);
        Epoch { week: w, tow_nanos: n, _scale: PhantomData }
    }

    /// Convert GPS time to Galileo system time (nominal 0-second offset).
    pub fn to_gst(&self) -> Epoch<GstScale> {
        Epoch { week: self.week, tow_nanos: self.tow_nanos, _scale: PhantomData }
    }

    /// Convert GPS time to GLONASS time with explicit leap seconds.
    pub fn to_glonass(&self, leap_seconds: i32) -> Epoch<GlonassScale> {
        let offset = GLONASS_HOURS_NANOS - (leap_seconds as i64) * (NANOS_PER_SEC as i64);
        let (w, n) = normalize_week_nanos(self.week, self.tow_nanos as i128 + offset as i128);
        Epoch { week: w, tow_nanos: n, _scale: PhantomData }
    }

    /// Convert GPS time to UTC with explicit leap seconds.
    pub fn to_utc(&self, leap_seconds: i32) -> Epoch<UtcScale> {
        let offset = (leap_seconds as i64) * (NANOS_PER_SEC as i64);
        let (w, n) = normalize_week_nanos(self.week, self.tow_nanos as i128 - offset as i128);
        Epoch { week: w, tow_nanos: n, _scale: PhantomData }
    }

    /// Convert to legacy `GpsTime`.
    pub fn to_gps_time(&self) -> GpsTime {
        (*self).into()
    }
}

impl Epoch<GstScale> {
    /// Convert Galileo system time to GPS time (nominal 0-second offset).
    pub fn to_gpst(&self) -> Epoch<GpsScale> {
        Epoch { week: self.week, tow_nanos: self.tow_nanos, _scale: PhantomData }
    }
}

impl Epoch<GlonassScale> {
    /// Convert GLONASS time to GPS time with explicit leap seconds.
    pub fn to_gpst(&self, leap_seconds: i32) -> Epoch<GpsScale> {
        let offset = -GLONASS_HOURS_NANOS + (leap_seconds as i64) * (NANOS_PER_SEC as i64);
        let (w, n) = normalize_week_nanos(self.week, self.tow_nanos as i128 + offset as i128);
        Epoch { week: w, tow_nanos: n, _scale: PhantomData }
    }
}

impl Epoch<UtcScale> {
    /// Convert UTC to GPS time with explicit leap seconds.
    pub fn to_gpst(&self, leap_seconds: i32) -> Epoch<GpsScale> {
        let offset = (leap_seconds as i64) * (NANOS_PER_SEC as i64);
        let (w, n) = normalize_week_nanos(self.week, self.tow_nanos as i128 + offset as i128);
        Epoch { week: w, tow_nanos: n, _scale: PhantomData }
    }
}

// ---------------------------------------------------------------------------
// Type Aliases
// ---------------------------------------------------------------------------

pub type GpsEpoch = Epoch<GpsScale>;
pub type BdtEpoch = Epoch<BdtScale>;
pub type GstEpoch = Epoch<GstScale>;
pub type GlonassEpoch = Epoch<GlonassScale>;
pub type UtcEpoch = Epoch<UtcScale>;

pub type BdtTime = Epoch<BdtScale>;
pub type GstTime = Epoch<GstScale>;
pub type GlonassTime = Epoch<GlonassScale>;
pub type UtcTime = Epoch<UtcScale>;

// ---------------------------------------------------------------------------
// Interop with legacy GpsTime
// ---------------------------------------------------------------------------

impl From<GpsTime> for Epoch<GpsScale> {
    fn from(gps: GpsTime) -> Self {
        Epoch::from_week_tow(gps.week, gps.tow)
    }
}

impl From<Epoch<GpsScale>> for GpsTime {
    fn from(epoch: Epoch<GpsScale>) -> Self {
        GpsTime::new(epoch.week(), epoch.tow_seconds())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalization_and_week_carry() {
        let e1: Epoch<GpsScale> = Epoch::from_week_nanos(100, NANOS_PER_WEEK);
        assert_eq!(e1.week(), 101);
        assert_eq!(e1.tow_nanos(), 0);

        let e2: Epoch<GpsScale> = Epoch::new(100, -1.0);
        assert_eq!(e2.week(), 99);
        assert_eq!(e2.tow_nanos(), NANOS_PER_WEEK - NANOS_PER_SEC);
    }

    #[test]
    fn test_rollover_safe_subtraction() {
        let t_curr: Epoch<GpsScale> = Epoch::new(2201, 10.0);
        let t_prev: Epoch<GpsScale> = Epoch::new(2200, 604790.0);
        let diff = t_curr - t_prev;
        assert_eq!(diff, TimeDelta::from_seconds(20.0));
        assert_eq!(diff.as_nanos(), 20 * (NANOS_PER_SEC as i64));
    }

    #[test]
    fn test_bdt_gpst_offset_roundtrip() {
        let bdt: Epoch<BdtScale> = Epoch::new(2000, 100.0);
        let gpst = bdt.to_gpst();
        assert_eq!(gpst.week(), 2000);
        assert_eq!(gpst.tow_seconds(), 114.0);

        let roundtrip = gpst.to_bdt();
        assert_eq!(roundtrip, bdt);
    }

    #[test]
    fn test_leap_seconds_utc_gpst() {
        let utc: Epoch<UtcScale> = Epoch::new(2000, 100.0);
        let gpst = utc.to_gpst(18);
        assert_eq!(gpst.tow_seconds(), 118.0);
        let utc_back = gpst.to_utc(18);
        assert_eq!(utc_back, utc);
    }

    #[test]
    fn test_is_within_tolerance() {
        let t1: Epoch<GpsScale> = Epoch::new(2000, 10.0);
        let t2: Epoch<GpsScale> = Epoch::new(2000, 10.005);
        assert!(t1.is_within(t2, TimeDelta::from_millis(10)));
        assert!(!t1.is_within(t2, TimeDelta::from_millis(2)));
    }

    #[test]
    fn test_epoch_key_continuity() {
        let e1: Epoch<GpsScale> = Epoch::new(2000, 604799.0);
        let e2: Epoch<GpsScale> = Epoch::new(2001, 1.0);
        let k1 = e1.to_key();
        let k2 = e2.to_key();
        assert_eq!(k2.continuous_ms() - k1.continuous_ms(), 2000);
        assert!(k2 > k1);
    }
}

#[cfg(test)]
#[path = "epoch_tests.rs"]
mod epoch_property;
