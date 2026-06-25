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

    /// Bug 14 placeholder: Step 2 SET frequency-dependent corrections.
    /// When implemented, applying Step 2 to the IERS benchmark position should
    /// reduce the height residual to <1 mm vs. Step-1-only output.
    #[test]
    #[ignore = "TODO(Bug14): IERS SET Step 2 not yet implemented"]
    fn test_solid_earth_tide_step2_diurnal_band() {
        // IERS Conventions 2010, example station at lat=45°, lon=0°, h=100m.
        // Expected Step-2 correction at J2000.0 epoch ≈ (dx, dy, dz) [mm].
        // Acceptance: each component within 1 mm of IERS tabulated values.
        let _t = GpsTime::new(1042, 0.0); // approximate J2000.0
        let _pos = Vector3::new(4_517_590.0, 0.0, 4_487_348.0); // ~45°N ECEF
        // When implemented: assert (step1+step2 - iers_ref).norm() < 0.001
        todo!("Implement IERS SET Step 2 corrections")
    }

    /// Bug 13 placeholder: Ocean Tide Loading correction.
    /// When implemented, OTL displacement at a coastal station must be non-zero
    /// and vary with epoch (at least 1 mm amplitude for any major tidal species).
    #[test]
    #[ignore = "TODO(Bug13): Ocean Tide Loading not yet implemented (requires BLQ parser)"]
    fn test_ocean_tide_loading_nonzero_coastal() {
        let t = GpsTime::new(2000, 43200.0);
        // Odaiba, Tokyo Bay — well-known coastal benchmark (lat≈35.6°N, lon≈139.8°E)
        let pos = Vector3::new(-3_960_000.0, 3_360_000.0, 3_690_000.0);
        let otl = ocean_tide_loading_ecef(t, pos);
        // Once implemented: norm should be > 0.001 m (1 mm) for a coastal site
        assert!(
            otl.norm() > 0.001,
            "OTL displacement should be > 1 mm at coastal site, got {:.3} m",
            otl.norm()
        );
    }
}
