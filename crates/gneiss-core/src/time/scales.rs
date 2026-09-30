//! Time scale typestates for GNSS and civil time systems.

use core::fmt;

/// Trait defining a temporal scale / coordinate system.
pub trait TimeScale: Copy + Clone + PartialEq + Eq + PartialOrd + Ord + fmt::Debug + 'static {
    /// Human-readable abbreviation for the time scale.
    const NAME: &'static str;
}

/// GPS System Time (GPST). Continuous time scale with origin 1980-01-06 00:00:00 UTC.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct GpsScale;

impl TimeScale for GpsScale {
    const NAME: &'static str = "GPST";
}

/// BeiDou Navigation Satellite System Time (BDT). Continuous time scale with origin 2006-01-01 00:00:00 UTC.
/// Permanent ICD relationship: GPST - BDT = 14.0 seconds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct BdtScale;

impl TimeScale for BdtScale {
    const NAME: &'static str = "BDT";
}

/// Galileo System Time (GST). Continuous time scale with origin 1999-08-21 23:59:47 UTC.
/// Nominal relationship: GPST == GST (within nanoseconds).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct GstScale;

impl TimeScale for GstScale {
    const NAME: &'static str = "GST";
}

/// GLONASS Time (GLONASST). Referencing UTC(SU) + 3 hours, containing leap second steps.
/// GPST = GLONASST - 10800s + leap_seconds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct GlonassScale;

impl TimeScale for GlonassScale {
    const NAME: &'static str = "GLONASST";
}

/// Coordinated Universal Time (UTC). Civil time stepped by leap seconds.
/// GPST = UTC + leap_seconds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct UtcScale;

impl TimeScale for UtcScale {
    const NAME: &'static str = "UTC";
}
