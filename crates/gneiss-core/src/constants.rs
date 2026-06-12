//! Physical and geodetic constants used throughout the Gneiss navigation engine.
//!
//! Values conform to WGS84 and IERS conventions unless otherwise noted.

/// Speed of light in vacuum (m/s) — ITU-R / IAU
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// WGS84 Earth rotation rate (rad/s)
pub const EARTH_ROTATION_RATE_RAD_S: f64 = 7.292_115_146_7e-5;

/// WGS84 semi-major axis (m)
pub const WGS84_SEMI_MAJOR_AXIS_M: f64 = 6_378_137.0;

/// WGS84 gravitational parameter GM (m³/s²)
pub const WGS84_GM_M3_S2: f64 = 3.986_005e14;

/// WGS84 J2 zonal harmonic (dimensionless)
pub const WGS84_J2: f64 = 1.082_627e-3;

/// MAD-to-sigma scaling factor for normal distributions
pub const MAD_NORMAL_SCALE_FACTOR: f64 = 1.4826;

/// Seconds per GPS week
pub const SECONDS_PER_GPS_WEEK: f64 = 604_800.0;

/// Seconds per day
pub const SECONDS_PER_DAY: f64 = 86_400.0;

/// Astronomical unit (m)
pub const ASTRONOMICAL_UNIT_M: f64 = 149_597_870_700.0;

/// Kilometers to meters conversion factor
pub const KM_TO_M: f64 = 1_000.0;

/// Microseconds to seconds conversion factor
pub const MICROSECONDS_TO_SECONDS: f64 = 1e-6;
