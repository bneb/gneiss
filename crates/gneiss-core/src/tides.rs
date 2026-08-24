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

/// TODO (Bug 13 — Tier 4): Ocean Tide Loading (OTL) corrections.
///
/// In coastal areas OTL displacements can reach 5–10 cm at diurnal/semi-diurnal
/// frequencies. Standard PPP engines parse BLQ files (produced by the IERS OTL
/// provider at <http://holt.oso.chalmers.se/loading/>) to obtain 11 tidal-constituent
/// amplitude/phase vectors for each station, then sum:
///
///   d_OTL(t) = Σ_k A_k cos(χ_k(t) + φ_k - u_k)
///
/// where A_k is the 3D amplitude vector, χ_k the astronomical argument,
/// φ_k the constituent phase, and u_k the ANTEX correction.
///
/// Implementation requirements:
///   1. Parse BLQ file → per-station `OceanTideParams` struct (11 rows × 6 cols).
///   2. Compute tidal arguments χ_k from IERS astronomical tables (Doodson numbers).
///   3. Evaluate and accumulate 11-constituent sum at each epoch.
///   4. Apply to receiver position vector before forming line-of-sight geometry.
///
/// This function is a placeholder that returns zero until BLQ support is added.
#[allow(unused_variables)]
pub fn ocean_tide_loading_ecef(_t: GpsTime, _pos_ecef: Vector3<f64>) -> Vector3<f64> {
    // TODO: parse BLQ, compute tidal arguments, evaluate 11-constituent sum.
    Vector3::zeros()
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

    fn test_set_displacement_magnitude_realistic() {
        // At any epoch, SET displacement should be 10-50 cm for a mid-latitude station.
        let t = GpsTime::new(2105, 0.0); // arbitrary epoch
        // P224 Sibley Volcanic: approximate ECEF
        let pos = Vector3::new(-2688201.0, -4265643.0, 3893778.0);
        let disp = solid_earth_tides_ecef(t, pos);
        let norm = disp.norm();
        assert!(
            norm > 0.05 && norm < 0.60,
            "SET displacement {} m outside expected [5, 60] cm range",
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
}
