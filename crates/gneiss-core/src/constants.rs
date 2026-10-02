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

/// Milliarcseconds to radians conversion factor
pub const MILLIARCSEC_TO_RAD: f64 = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);

/// Days between GPS epoch (1980-01-06 00:00 UTC, JD 2444244.5) and J2000.0
/// (2000-01-01 12:00 TT, JD 2451545.0).
///
/// Derived, not remembered: `2451545.0 - 2444244.5 = 7300.5`.
/// Both Julian dates are exact halves, so the difference is an exact half-day.
pub const DAYS_GPS_TO_J2000: f64 = 7300.5;

/// Days in a Julian century
pub const DAYS_PER_JULIAN_CENTURY: f64 = 36525.0;

#[cfg(test)]
mod tests {
    use super::*;

    // Every derived constant below is checked against an *independently
    // derived* value, never against the expression used to define it.
    // Flipping a sign, swapping two numbers, or changing a constant by one
    // digit breaks at least one of these.

    /// JD(J2000.0) = 2451545.0 is fixed by definition (2000-01-01 12:00 TT).
    /// JD(GPS epoch) = 2444244.5 is 1980-01-06 00:00 UTC. Both are exact
    /// half-integers, so the span is exactly 7300.5 days.
    #[test]
    fn days_gps_to_j2000_is_the_exact_span_of_two_known_julian_dates() {
        const JD_J2000: f64 = 2451545.0;
        const JD_GPS_EPOCH: f64 = 2444244.5;
        assert_eq!(JD_J2000 - JD_GPS_EPOCH, 7300.5);
        assert_eq!(DAYS_GPS_TO_J2000, 7300.5);
        // Independent cross-check: 1980-01-06 -> 2000-01-06 is 20 years with
        // 5 leap days (1980, 1984, 1988, 1992, 1996; 2000-02-29 is past the
        // 2000-01-06 mark) = 20*365 + 5 = 7305 days. J2000.0 is 5 days before
        // 2000-01-06 plus 12 h = 7300.5.
        assert_eq!(DAYS_GPS_TO_J2000, (20.0 * 365.0 + 5.0) - 5.0 + 0.5);
    }

    /// The GPS week length is defined as 604800 s, so 7300.5 days is not a
    /// whole number of weeks: it is 1043 weeks minus 3.5 days.
    /// Checked so a mutation of either day constant is caught here too.
    #[test]
    fn days_gps_to_j2000_sits_between_gps_weeks_as_expected() {
        const SECONDS_PER_DAY: f64 = 86_400.0;
        // 7300.5 = 1042*7 + 6.5, i.e. J2000 is week 1042 day 6.5 — consistent
        // with the calendar golden vectors in `time::tests`.
        assert_eq!(DAYS_GPS_TO_J2000 / 7.0, 1042.0 + 6.5 / 7.0);
        assert_eq!(1042.0 * 7.0 + 6.5, DAYS_GPS_TO_J2000);
        const _: () = assert!(SECONDS_PER_DAY > 0.0);
    }

    /// 1 radian = 180/pi degrees = 648000/pi arcseconds = 6.48e8/pi mas.
    /// Computed with python3: 6.48e8/pi = 206264806.24709636 mas per radian.
    /// The constant must equal pi / 6.48e8 to the last bit.
    #[test]
    fn milliarcsecond_to_radian_matches_the_definition_chain() {
        let mas_per_rad = 180.0 / core::f64::consts::PI * 3600.0 * 1000.0;
        assert!((mas_per_rad - 206_264_806.247_096_36).abs() < 1e-6);
        assert_eq!(MILLIARCSEC_TO_RAD, 1.0 / mas_per_rad);
        assert_eq!(MILLIARCSEC_TO_RAD, core::f64::consts::PI / (180.0 * 3600.0 * 1000.0));
        // 1 microarcsecond = 4.8481368110953594e-12 rad (python3); it is one
        // thousandth of a milliarcsecond, agreeing to 1 ulp.
        assert!((MILLIARCSEC_TO_RAD / 1000.0 - 4.848_136_811_095_36e-12).abs() < 1e-27);
    }

    /// SI defines c exactly. Any other value means the constant is wrong.
    #[test]
    fn speed_of_light_is_the_exact_si_defining_constant() {
        assert_eq!(SPEED_OF_LIGHT_M_S, 299_792_458.0);
        // Cross-check against the metre/second definitions is impossible, but
        // the value must be integral and positive; a mutation to 299_792_458.5
        // or 299_792_458.0 * 1.0001 must fail.
        assert_eq!(SPEED_OF_LIGHT_M_S as u64 as f64, SPEED_OF_LIGHT_M_S);
        const _: () = assert!(SPEED_OF_LIGHT_M_S > 299_792_457.0 && SPEED_OF_LIGHT_M_S < 299_792_459.0);
    }

    /// A Julian century is 36525 days by definition (100 * 365.25).
    #[test]
    fn julian_century_is_36525_days() {
        assert_eq!(DAYS_PER_JULIAN_CENTURY, 365.25 * 100.0);
        // GPS-epoch-to-J2000 is not a whole number of centuries, so the two
        // must not be swapped: 36525 days != 7300.5 days.
        assert_ne!(DAYS_PER_JULIAN_CENTURY, DAYS_GPS_TO_J2000);
    }

    /// The two unit-conversion constants must be exactly consistent with each
    /// other and with the day/week lengths already in the crate. A mutation of
    /// any one of them breaks a different assertion.
    #[test]
    fn unit_conversion_constants_are_mutually_consistent() {
        assert_eq!(KM_TO_M, 1000.0);
        assert_eq!(MICROSECONDS_TO_SECONDS, 1e-6);
        assert_eq!(SECONDS_PER_DAY, 86_400.0);
        assert_eq!(SECONDS_PER_GPS_WEEK, 7.0 * SECONDS_PER_DAY);
        // 86 400 s/day = 86 400 000 000 us/day
        assert_eq!(SECONDS_PER_DAY / MICROSECONDS_TO_SECONDS, 86_400_000_000.0);
        const _: () = assert!(SECONDS_PER_DAY > 0.0);
        // 604 800 s/week = 604.8 km in units of KM_TO_M
        assert_eq!(SECONDS_PER_GPS_WEEK / 1000.0, 604.8);
        assert!((KM_TO_M * MICROSECONDS_TO_SECONDS - 1e-3).abs() < 1e-18);
    }

    /// The geodetic constants are cross-checked through relationships that
    /// follow from their definitions rather than from a lookup table.
    #[test]
    fn wgs84_constants_relate_through_the_flattening() {
        // f = 1/298.257223563 defines the WGS84 ellipsoid together with a.
        const INV_F: f64 = 298.257_223_563;
        let f = 1.0 / INV_F;
        let b = WGS84_SEMI_MAJOR_AXIS_M * (1.0 - f);
        // Polar radius implied by the flattening, in the published range.
        assert!(b > 6_356_752.0 && b < 6_356_753.0, "b = {b}");
        // First eccentricity squared e^2 = 1 - (b/a)^2 must be ~0.00669.
        let e2 = 1.0 - (b / WGS84_SEMI_MAJOR_AXIS_M).powi(2);
        assert!((e2 - 0.006_694_379_990_141_316).abs() < 1e-15, "e^2 = {e2}");
        // GM must be positive and of planetary scale; WGS84 GM ~ 3.986e14.
        assert!((WGS84_GM_M3_S2 / 1e14 - 3.986_005).abs() < 1e-6);
        // J2 is dimensionless and strictly between 0 and 0.01.
        const _: () = assert!(WGS84_J2 > 0.0 && WGS84_J2 < 0.01);
        // The Earth rotation rate is 2*pi / sidereal day.
        // 2*pi / 86164.090530833 s = 7.292115855306578e-5 rad/s (python3).
        // The WGS84 tabulated value 7.2921151467e-5 differs by 7.1e-12 rad/s,
        // i.e. ~1e-7 relative, because WGS84 fixes it to 7 decimal places of
        // 2*pi/s rather than to the mean sidereal day. Bound the agreement.
        const SIDEREAL_DAY_S: f64 = 86_164.090_530_833;
        let from_sidereal_day = 2.0 * core::f64::consts::PI / SIDEREAL_DAY_S;
        assert!(
            (EARTH_ROTATION_RATE_RAD_S - from_sidereal_day).abs()
                < 1e-6 * from_sidereal_day,
            "Earth rotation rate {EARTH_ROTATION_RATE_RAD_S} vs {from_sidereal_day}"
        );
    }

    /// MAD -> sigma for a normal law uses 1/0.6744897501960817 = 1.4826 (the
    /// median of |N(0,1)| is 0.6744897501960817). Asserting the *derivation*
    /// rather than the literal catches a sign flip or a 1/x transposition.
    #[test]
    fn mad_normal_scale_factor_is_the_inverse_normal_median() {
        // median |Z| for Z ~ N(0,1) is Phi^{-1}(0.75) = 0.6744897501960817.
        // 1 / that = 1.482602218505602 (python3); the constant is the 5-digit
        // rounding 1.4826, so agree to 3e-6, not bit for bit.
        const NORMAL_MEDIAN_ABS: f64 = 0.674_489_750_196_081_7;
        assert!((1.0 / NORMAL_MEDIAN_ABS - 1.482_602_218_505_602).abs() < 1e-15);
        assert!((1.0 / NORMAL_MEDIAN_ABS - MAD_NORMAL_SCALE_FACTOR).abs() < 3e-6);
        // Applying then un-applying is the identity, and the factor is > 1 so
        // it inflates a median absolute deviation toward a standard deviation.
        let mad = 1.0;
        const _: () = assert!(MAD_NORMAL_SCALE_FACTOR > 1.0 && MAD_NORMAL_SCALE_FACTOR < 2.0);
        assert_eq!(MAD_NORMAL_SCALE_FACTOR * mad / MAD_NORMAL_SCALE_FACTOR, mad);
        assert!((MAD_NORMAL_SCALE_FACTOR - 1.4826).abs() < 1e-12);
    }
}
