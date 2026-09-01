//! Relativistic clock dilation and Shapiro gravitational time delay.
//!
//! Implements IERS Conventions (2010) Chapter 10 & 11 relativistic corrections.

use nalgebra::Vector3;

const SPEED_OF_LIGHT: f64 = 299_792_458.0; // m/s
const GM_EARTH: f64 = 3.986004418e14;       // Earth gravitational constant (m^3/s^2)

/// Computes the periodic relativistic satellite clock correction (meters).
///
/// Implements IERS Conventions (2010) eq. 10.1:
///   \Delta \rho_{\text{rel}} = -\frac{2 \mathbf{r}^s \cdot \mathbf{v}^s}{c}
///
/// Note: Multiply by 1/c to obtain seconds (\Delta t_{\text{rel}}).
#[inline]
pub fn periodic_relativistic_range_correction(
    sat_pos_ecef: &Vector3<f64>,
    sat_vel_ecef: &Vector3<f64>,
) -> f64 {
    let r_dot_v = sat_pos_ecef.dot(sat_vel_ecef);
    -2.0 * r_dot_v / SPEED_OF_LIGHT
}

/// Computes the gravitational Shapiro time delay in range units (meters).
///
/// Implements IERS Conventions (2010) eq. 11.2:
///   \Delta \rho_{\text{grav}} = \frac{2 G M_E}{c^2} \ln\left( \frac{r_s + r_r + \rho}{r_s + r_r - \rho} \right)
pub fn gravitational_shapiro_delay(
    sat_pos_ecef: &Vector3<f64>,
    rx_pos_ecef: &Vector3<f64>,
) -> f64 {
    let r_s = libm::sqrt(sat_pos_ecef.x * sat_pos_ecef.x + sat_pos_ecef.y * sat_pos_ecef.y + sat_pos_ecef.z * sat_pos_ecef.z);
    let r_r = libm::sqrt(rx_pos_ecef.x * rx_pos_ecef.x + rx_pos_ecef.y * rx_pos_ecef.y + rx_pos_ecef.z * rx_pos_ecef.z);
    let los = sat_pos_ecef - rx_pos_ecef;
    let rho = libm::sqrt(los.x * los.x + los.y * los.y + los.z * los.z);

    if r_s < 1.0 || r_r < 1.0 || (r_s + r_r - rho) <= 0.0 {
        return 0.0;
    }

    let factor = (2.0 * GM_EARTH) / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    let ratio = (r_s + r_r + rho) / (r_s + r_r - rho);
    factor * libm::log(ratio)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_periodic_relativistic_range_correction_textbook_vector() {
        // Known test case: r = 26560 km with radial velocity vr = 50.0 m/s
        let pos = Vector3::new(26_560_000.0, 0.0, 0.0);
        let vel = Vector3::new(50.0, 3_874.0, 0.0);

        let delta_rho = periodic_relativistic_range_correction(&pos, &vel);
        let delta_t_ns = (delta_rho / SPEED_OF_LIGHT) * 1e9;

        // Exact analytical: -2 * (26560000 * 50) / c = -8.85946 m (-29.552 ns)
        assert!((delta_rho - (-8.859463)).abs() < 1e-4, "delta_rho was {:.6}m", delta_rho);
        assert!((delta_t_ns - (-29.552)).abs() < 1e-2, "delta_t was {:.3}ns", delta_t_ns);
    }

    #[test]
    fn test_gravitational_shapiro_delay_zenith_satellite() {
        // Station on equator, satellite at zenith on GPS orbit (26560 km)
        let rx = Vector3::new(6_378_137.0, 0.0, 0.0);
        let sat = Vector3::new(26_560_000.0, 0.0, 0.0);

        let delay = gravitational_shapiro_delay(&sat, &rx);

        // Analytical textbook value: 2*GM/c^2 * ln(r_s/r_r) = 8.870056 mm * ln(4.164225) = 12.653 mm
        assert!((delay - 0.012653).abs() < 1e-5, "Shapiro delay was {:.6}m", delay);
    }
}
