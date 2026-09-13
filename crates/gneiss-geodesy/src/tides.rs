//! Solid Earth Tide displacement modeling according to IERS Conventions (2010).
//!
//! Evaluates degree-2 and degree-3 elastic deformation induced by the Moon and Sun,
//! using frequency-independent Love and Shida numbers (Step 1 of dehanttideinel.f).

use core::f64::consts::PI;
use nalgebra::Vector3;

const GM_RATIO_MOON: f64 = 0.0123000371; // M_moon / M_earth
const GM_RATIO_SUN: f64 = 332946.0482;   // M_sun / M_earth

#[inline]
fn vec_norm(v: &Vector3<f64>) -> f64 {
    libm::sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

#[inline]
fn deg2rad(deg: f64) -> f64 {
    deg * (PI / 180.0)
}

#[inline]
fn rotate_eci_to_ecef(v: &Vector3<f64>, gmst_rad: f64) -> Vector3<f64> {
    let (s, c) = libm::sincos(gmst_rad);
    Vector3::new(c * v.x + s * v.y, -s * v.x + c * v.y, v.z)
}

fn compute_sun_eci(t_cent: f64, eps: f64) -> Vector3<f64> {
    let l_sun = deg2rad(280.460 + 36000.770 * t_cent);
    let g_sun = deg2rad(357.528 + 35999.050 * t_cent);
    let lambda_sun = l_sun + deg2rad(1.915 * libm::sin(g_sun) + 0.020 * libm::sin(2.0 * g_sun));
    let r_sun = 1.495978707e11 * (1.00014 - 0.01671 * libm::cos(g_sun) - 0.00014 * libm::cos(2.0 * g_sun));
    Vector3::new(
        r_sun * libm::cos(lambda_sun),
        r_sun * libm::sin(lambda_sun) * libm::cos(eps),
        r_sun * libm::sin(lambda_sun) * libm::sin(eps),
    )
}

fn compute_moon_eci(t_cent: f64, eps: f64) -> Vector3<f64> {
    let l_moon = deg2rad(218.316 + 481267.881 * t_cent);
    let m_moon = deg2rad(134.963 + 477198.867 * t_cent);
    let lambda_moon = l_moon + deg2rad(6.289 * libm::sin(m_moon));
    let r_moon = 384400000.0 * (1.0 - 0.0549 * libm::cos(m_moon));
    Vector3::new(
        r_moon * libm::cos(lambda_moon),
        r_moon * libm::sin(lambda_moon) * libm::cos(eps),
        r_moon * libm::sin(lambda_moon) * libm::sin(eps),
    )
}

/// Compute solar and lunar approximate ECEF coordinates at a given GPS time.
pub fn solar_lunar_positions(gps_tow: f64, gps_week: u32) -> (Vector3<f64>, Vector3<f64>) {
    let t_days = (gps_week as f64) * 7.0 + gps_tow / 86400.0;
    let d_j2000 = t_days - 7300.5;
    let t_cent = d_j2000 / 36525.0;
    let eps = deg2rad(23.439291 - 0.0130042 * t_cent);

    let sun_eci = compute_sun_eci(t_cent, eps);
    let moon_eci = compute_moon_eci(t_cent, eps);

    let gmst_rad = deg2rad((280.46061837 + 360.98564736629 * d_j2000) % 360.0);

    (rotate_eci_to_ecef(&sun_eci, gmst_rad), rotate_eci_to_ecef(&moon_eci, gmst_rad))
}

#[derive(Debug, Clone, Copy)]
struct LoveNumbers {
    h2: f64,
    l2: f64,
    h3: f64,
    l3: f64,
}

/// Computes the 3D Solid Earth Tide displacement vector in ECEF frame (meters).
///
/// Implements IERS Conventions (2010) Step 1 spherical Love number elastic deformation.
pub fn solid_earth_tide(
    station_ecef: Vector3<f64>,
    gps_tow: f64,
    gps_week: u32,
) -> Vector3<f64> {
    let r_station = vec_norm(&station_ecef);
    if r_station < 1.0 {
        return Vector3::zeros();
    }
    let r_hat = station_ecef / r_station;
    let sin_lat = r_hat.z;
    let p2_lat = 0.5 * (3.0 * sin_lat * sin_lat - 1.0);

    let love = LoveNumbers {
        h2: 0.6078 - 0.0006 * p2_lat,
        l2: 0.0847 + 0.0002 * p2_lat,
        h3: 0.292,
        l3: 0.015,
    };

    let (sun_pos, moon_pos) = solar_lunar_positions(gps_tow, gps_week);

    let mut delta_r = Vector3::zeros();
    delta_r += body_tide_displacement(r_station, &r_hat, &moon_pos, GM_RATIO_MOON, love);
    delta_r += body_tide_displacement(r_station, &r_hat, &sun_pos, GM_RATIO_SUN, love);

    delta_r
}

fn body_tide_displacement(
    r_station: f64,
    r_hat: &Vector3<f64>,
    body_pos: &Vector3<f64>,
    gm_ratio: f64,
    love: LoveNumbers,
) -> Vector3<f64> {
    let r_body = vec_norm(body_pos);
    if r_body < 1.0 {
        return Vector3::zeros();
    }
    let r_body_hat = body_pos / r_body;
    let z = r_hat.dot(&r_body_hat); // Cosine of zenith angle

    // Degree 2 displacement (IERS 2010 eq. 7.5)
    let factor2 = gm_ratio * (r_station * r_station * r_station * r_station / (r_body * r_body * r_body));
    let rad2 = (1.5 * love.h2 * z * z - 0.5 * love.h2 - 3.0 * love.l2 * z * z) * r_hat;
    let tang2 = (3.0 * love.l2 * z) * r_body_hat;
    let disp2 = factor2 * (rad2 + tang2);

    // Degree 3 displacement (IERS 2010 eq. 7.6)
    let factor3 = gm_ratio * (r_station * r_station * r_station * r_station * r_station / (r_body * r_body * r_body * r_body));
    let z2 = z * z;
    let rad3 = (love.h3 * (2.5 * z * z2 - 1.5 * z) - love.l3 * (7.5 * z * z2 - 1.5 * z)) * r_hat;
    let tang3 = (love.l3 * (7.5 * z2 - 1.5)) * r_body_hat;
    let disp3 = factor3 * (rad3 + tang3);

    disp2 + disp3
}

/// 11-constituent Ocean Tide Loading (OTL) harmonic amplitudes and phases for a station.
///
/// Order of constituents (IERS Conventions 2010 Chapter 7.2):
/// M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa.
#[derive(Debug, Clone, Copy, Default)]
pub struct OtlHarmonics {
    /// Radial/Up amplitudes (meters) for the 11 constituents.
    pub amp_radial_m: [f64; 11],
    /// Radial/Up phase lags (radians) for the 11 constituents.
    pub phase_radial_rad: [f64; 11],
    /// East amplitudes (meters).
    pub amp_east_m: [f64; 11],
    /// East phase lags (radians).
    pub phase_east_rad: [f64; 11],
    /// North amplitudes (meters).
    pub amp_north_m: [f64; 11],
    /// North phase lags (radians).
    pub phase_north_rad: [f64; 11],
}

/// Angular speeds of the 11 principal tidal constituents in radians per second.
pub const OTL_FREQUENCIES_RAD_S: [f64; 11] = [
    1.405189e-4, // M2
    1.454441e-4, // S2
    1.378797e-4, // N2
    1.458423e-4, // K2
    7.292117e-5, // K1
    6.759774e-5, // O1
    7.252295e-5, // P1
    6.495854e-5, // Q1
    0.053234e-4, // Mf
    0.026392e-4, // Mm
    0.001991e-4, // Ssa
];

/// Computes the 3D Ocean Tide Loading crustal displacement in Local (East, North, Up) coordinates (meters).
pub fn ocean_tide_loading_enu(harmonics: &OtlHarmonics, gps_tow: f64) -> Vector3<f64> {
    let mut d_e = 0.0;
    let mut d_n = 0.0;
    let mut d_u = 0.0;

    for (i, &omega) in OTL_FREQUENCIES_RAD_S.iter().enumerate() {
        let arg = omega * gps_tow;
        d_u += harmonics.amp_radial_m[i] * libm::cos(arg - harmonics.phase_radial_rad[i]);
        d_e += harmonics.amp_east_m[i] * libm::cos(arg - harmonics.phase_east_rad[i]);
        d_n += harmonics.amp_north_m[i] * libm::cos(arg - harmonics.phase_north_rad[i]);
    }

    Vector3::new(d_e, d_n, d_u)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solid_earth_tide_bounded_magnitude() {
        let station = Vector3::new(6378137.0, 0.0, 0.0);
        let tide = solid_earth_tide(station, 43200.0, 2137);

        let magnitude = vec_norm(&tide);
        assert!(magnitude < 0.50, "Earth tide magnitude {:.4}m exceeded 0.50m physical bound", magnitude);
        assert!(magnitude > 0.01, "Earth tide magnitude {:.4}m should be non-zero and at centimeter level", magnitude);
    }

    #[test]
    fn test_solid_earth_tide_zero_position_safe() {
        let tide = solid_earth_tide(Vector3::zeros(), 0.0, 0);
        assert_eq!(tide, Vector3::zeros());
    }

    #[test]
    fn test_ocean_tide_loading_enu_evaluation() {
        let mut harm = OtlHarmonics::default();
        // M2 amplitude 2.5 cm vertical, 0.5 cm east
        harm.amp_radial_m[0] = 0.025;
        harm.amp_east_m[0] = 0.005;

        let disp = ocean_tide_loading_enu(&harm, 0.0);
        assert!((disp.z - 0.025).abs() < 1e-6);
        assert!((disp.x - 0.005).abs() < 1e-6);
        assert_eq!(disp.y, 0.0);
    }
}
