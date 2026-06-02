use nalgebra::Vector3;
use crate::time::GpsTime;
use crate::sun::{sun_position_ecef, moon_position_ecef};

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
/// Equatorial radius of the Earth (m)
const R_EARTH: f64 = 6378137.0;

/// Calculates the Solid Earth Tides (SET) displacement vector in ECEF frame.
/// Returns the displacement (dx, dy, dz) in meters.
pub fn solid_earth_tides_ecef(t: GpsTime, pos_ecef: Vector3<f64>) -> Vector3<f64> {
    let r_sun = sun_position_ecef(t);
    let r_moon = moon_position_ecef(t);
    
    let mut disp = Vector3::zeros();
    
    disp += compute_tide_contribution(pos_ecef, r_sun, GM_SUN);
    disp += compute_tide_contribution(pos_ecef, r_moon, GM_MOON);
    
    disp
}

fn compute_tide_contribution(pos_ecef: Vector3<f64>, r_celestial: Vector3<f64>, gm: f64) -> Vector3<f64> {
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
    let coeff = (gm / GM_EARTH) * libm::pow(R_EARTH, 4.0) / libm::pow(dist, 3.0);
    
    // Radial displacement component (h2)
    let dr = H2 * r_hat * (1.5 * dot * dot - 0.5);
    
    // Transverse displacement component (l2)
    let dt = 3.0 * L2 * dot * (r_celestial_hat - dot * r_hat);
    
    coeff * (dr + dt)
}
