use crate::sun::{moon_position_ecef, sun_position_ecef};
use crate::time::GpsTime;
use nalgebra::Vector3;

/// Nominal degree 2 Love number
const H2: f64 = 0.609;
/// Nominal degree 2 Shida number
const L2: f64 = 0.085;
/// Gravitational parameter of the Sun (m^3 / s^2)
const GM_SUN: f64 = 1.32712440042E20;
/// Gravitational parameter of the Moon (m^3 / s^2)
const GM_MOON: f64 = 4.902800066E12;
/// Gravitational parameter of the Earth (m^3 / s^2)
const GM_EARTH: f64 = 3.986004415E14;

/// Calculates the Solid Earth Tides (SET) displacement vector in ECEF frame.
/// Returns the displacement (dx, dy, dz) in meters.
///
/// Implements IERS Conventions 2010, Ch. 7.1.1, Step 1 only.
///
/// TODO (Bug 14 — Tier 4): Add Step 2 corrections from IERS Conventions 2010,
/// Section 7.1.1, Eq. (7.5)–(7.6):
///   - Frequency-dependent diurnal band corrections using Earth's resonance kernel.
///   - Removal of permanent deformation (zero-frequency term) to produce tide-free
///     coordinates compatible with ITRF.
///     This introduces ≤1–2 cm systematic error in height and horizontal components
///     that is currently absorbed by the troposphere zenith wet delay state.
///
/// Reference: IERS Conventions 2010, Section 7.1.1, Tables 7.3a/7.3b.
pub fn solid_earth_tides_ecef(t: GpsTime, pos_ecef: Vector3<f64>) -> Vector3<f64> {
    let r_sun = sun_position_ecef(t);
    let r_moon = moon_position_ecef(t);

    let mut disp = Vector3::zeros();

    disp += compute_tide_contribution(pos_ecef, r_sun, GM_SUN);
    disp += compute_tide_contribution(pos_ecef, r_moon, GM_MOON);

    disp
}

/// 11 standard IERS Ocean Tide Loading constituents (M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa).
pub const OTL_CONSTITUENT_COUNT: usize = 11;

/// Standard constituent angular frequencies in rad/s.
pub const OTL_ANGULAR_FREQUENCIES_RAD_S: [f64; OTL_CONSTITUENT_COUNT] = [
    1.405189025e-4, // M2
    1.454441043e-4, // S2
    1.378796995e-4, // N2
    1.458421570e-4, // K2
    7.292115855e-5, // K1
    6.759774415e-5, // O1
    7.252294578e-5, // P1
    6.495854122e-5, // Q1
    5.3234144e-6,   // Mf
    2.639203e-6,    // Mm
    3.98213e-7,     // Ssa
];

/// 11-constituent Ocean Tide Loading parameter set for a site (amplitudes in meters, phases in radians).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OceanTideParams {
    /// Radial / Up amplitude for 11 constituents (m).
    pub amp_radial_m: [f64; OTL_CONSTITUENT_COUNT],
    /// West amplitude for 11 constituents (m).
    pub amp_west_m: [f64; OTL_CONSTITUENT_COUNT],
    /// South amplitude for 11 constituents (m).
    pub amp_south_m: [f64; OTL_CONSTITUENT_COUNT],
    /// Radial / Up phase for 11 constituents (rad).
    pub ph_radial_rad: [f64; OTL_CONSTITUENT_COUNT],
    /// West phase for 11 constituents (rad).
    pub ph_west_rad: [f64; OTL_CONSTITUENT_COUNT],
    /// South phase for 11 constituents (rad).
    pub ph_south_rad: [f64; OTL_CONSTITUENT_COUNT],
}

impl OceanTideParams {
    /// Zero parameters (no ocean tide displacement).
    #[must_use]
    pub const fn zeros() -> Self {
        Self {
            amp_radial_m: [0.0; OTL_CONSTITUENT_COUNT],
            amp_west_m: [0.0; OTL_CONSTITUENT_COUNT],
            amp_south_m: [0.0; OTL_CONSTITUENT_COUNT],
            ph_radial_rad: [0.0; OTL_CONSTITUENT_COUNT],
            ph_west_rad: [0.0; OTL_CONSTITUENT_COUNT],
            ph_south_rad: [0.0; OTL_CONSTITUENT_COUNT],
        }
    }

    /// Compute 3D ENU displacement in metres (East, North, Up) at epoch `t`.
    #[must_use]
    pub fn displacement_enu(&self, t: GpsTime) -> Vector3<f64> {
        let t_sec = t.week as f64 * 604_800.0 + t.tow;
        let mut d_up = 0.0;
        let mut d_west = 0.0;
        let mut d_south = 0.0;

        for (i, &omega) in OTL_ANGULAR_FREQUENCIES_RAD_S.iter().enumerate() {
            let arg_up = omega * t_sec - self.ph_radial_rad[i];
            let arg_w = omega * t_sec - self.ph_west_rad[i];
            let arg_s = omega * t_sec - self.ph_south_rad[i];

            d_up += self.amp_radial_m[i] * libm::cos(arg_up);
            d_west += self.amp_west_m[i] * libm::cos(arg_w);
            d_south += self.amp_south_m[i] * libm::cos(arg_s);
        }

        Vector3::new(-d_west, -d_south, d_up)
    }
}

/// Ocean Tide Loading (OTL) displacement in ECEF frame.
pub fn ocean_tide_loading_ecef(t: GpsTime, pos_ecef: Vector3<f64>, params: &OceanTideParams) -> Vector3<f64> {
    let enu = params.displacement_enu(t);
    crate::coords::enu_to_ecef(pos_ecef, enu)
}

fn compute_tide_contribution(
    pos_ecef: Vector3<f64>,
    r_celestial: Vector3<f64>,
    gm: f64,
) -> Vector3<f64> {
    let r_norm = pos_ecef.norm();
    if r_norm < 1e-6 {
        return Vector3::zeros();
    }
    let r_hat = pos_ecef / r_norm;

    let dist = r_celestial.norm();
    if dist < 1e-6 {
        return Vector3::zeros();
    }
    let r_celestial_hat = r_celestial / dist;

    let dot = r_celestial_hat.dot(&r_hat);

    // Scale coefficient
    let coeff = (gm / GM_EARTH) * libm::pow(r_norm, 4.0) / libm::pow(dist, 3.0);

    // Radial displacement component (h2)
    let dr = H2 * r_hat * (1.5 * dot * dot - 0.5);

    // Transverse displacement component (l2)
    let dt = 3.0 * L2 * dot * (r_celestial_hat - dot * r_hat);

    coeff * (dr + dt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::GpsTime;
    use nalgebra::Vector3;

    #[test]
    fn test_celestial_distances_at_epoch() {
        let t = GpsTime::new(2105, 0.0);
        let r_sun = crate::sun::sun_position_ecef(t);
        let r_moon = crate::sun::moon_position_ecef(t);
        // These assertions validate the celestial position functions.
        assert!(
            r_sun.norm() > 1e10 && r_sun.norm() < 2e11,
            "Sun at {} m — expected ~1.5e11",
            r_sun.norm()
        );
        assert!(
            r_moon.norm() > 3e8 && r_moon.norm() < 5e8,
            "Moon at {} m — expected ~3.8e8",
            r_moon.norm()
        );
    }

    #[test]
    fn test_set_displacement_magnitude_realistic() {
        // At any epoch, SET displacement should be 1-50 cm for a mid-latitude station.
        let t = GpsTime::new(2105, 0.0); // arbitrary epoch
        // P224 Sibley Volcanic: approximate ECEF
        let pos = Vector3::new(-2688201.0, -4265643.0, 3893778.0);
        let disp = solid_earth_tides_ecef(t, pos);
        let norm = disp.norm();
        assert!(
            norm > 0.01 && norm < 0.60,
            "SET displacement {:.2} cm outside expected [1, 60] cm range",
            norm * 100.0
        );
    }

    #[test]
    fn test_set_varies_with_epoch() {
        // Displacement must change over hours (tidal period ~12.4 h)
        let pos = Vector3::new(-2688201.0, -4265643.0, 3893778.0);
        let d0 = solid_earth_tides_ecef(GpsTime::new(2105, 0.0), pos);
        let d6 = solid_earth_tides_ecef(GpsTime::new(2105, 21600.0), pos); // +6h
        let d12 = solid_earth_tides_ecef(GpsTime::new(2105, 43200.0), pos); // +12h
        assert!(
            (d0 - d6).norm() > 0.001 || (d0 - d12).norm() > 0.001,
            "SET displacement should vary over 6-12 h periods"
        );
    }

    #[test]
    fn test_set_differential_between_nearby_stations_small() {
        // Two stations 15 km apart should have very similar SET displacement.
        // Differential is what matters for DD; expect < 5 mm at this distance.
        let t = GpsTime::new(2105, 0.0);
        let p224 = Vector3::new(-2688201.0, -4265643.0, 3893778.0);
        // P181 is ~15km NW
        let p181 = Vector3::new(-2697941.0, -4255089.0, 3898009.0);
        let set_p224 = solid_earth_tides_ecef(t, p224);
        let set_p181 = solid_earth_tides_ecef(t, p181);
        let diff = (set_p224 - set_p181).norm();
        assert!(
            diff < 0.005,
            "SET differential between 15-km stations should be < 5 mm, got {:.3} mm",
            diff * 1000.0
        );
    }

    #[test]
    fn test_ocean_tide_loading_zeros_and_m2_response() {
        let t0 = GpsTime::new(0, 0.0);
        let pos = Vector3::new(-2688201.0, -4265643.0, 3893778.0);
        let p_zero = OceanTideParams::zeros();
        assert_eq!(p_zero.displacement_enu(t0), Vector3::zeros());
        assert_eq!(ocean_tide_loading_ecef(t0, pos, &p_zero), Vector3::zeros());

        let mut p_m2 = OceanTideParams::zeros();
        p_m2.amp_radial_m[0] = 0.03; // 3 cm M2 vertical
        p_m2.amp_west_m[0] = 0.01;   // 1 cm M2 west
        let enu0 = p_m2.displacement_enu(t0);
        assert!((enu0.z - 0.03).abs() < 1e-12);
        assert!((enu0.x - (-0.01)).abs() < 1e-12);
        let ecef0 = ocean_tide_loading_ecef(t0, pos, &p_m2);
        assert!((ecef0.norm() - enu0.norm()).abs() < 1e-9);

        // Arbitrary epoch: displacement bounded by amplitude sum
        let t_arb = GpsTime::new(2105, 12345.0);
        let enu_arb = p_m2.displacement_enu(t_arb);
        assert!(enu_arb.z.abs() <= 0.03 + 1e-12);
        assert!(enu_arb.x.abs() <= 0.01 + 1e-12);
    }
}
