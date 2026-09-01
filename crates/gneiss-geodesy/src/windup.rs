//! Continuous Phase Windup Modeling for Circularly Polarized GNSS Transmissions.
//!
//! Implements the Wu et al. (1993) phase windup formulation for RHCP satellite transmissions
//! with continuous cycle wrap tracking across passes.

use core::f64::consts::PI;
use nalgebra::Vector3;

#[inline]
fn vec_norm(v: &Vector3<f64>) -> f64 {
    libm::sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

/// Continuous phase windup tracker across epochs.
#[derive(Debug, Clone, Copy, Default)]
pub struct PhaseWindupTracker {
    pub prev_windup_rad: f64,
    pub is_initialized: bool,
}

impl PhaseWindupTracker {
    pub fn new() -> Self {
        Self {
            prev_windup_rad: 0.0,
            is_initialized: false,
        }
    }

    /// Computes the continuous phase windup correction in radians.
    ///
    /// Parameters:
    /// - `sat_pos`: Satellite ECEF position
    /// - `sun_pos`: Sun ECEF position
    /// - `rx_pos`: Receiver ECEF position
    /// - `rx_up`: Receiver local Up unit vector
    /// - `rx_north`: Receiver local North unit vector
    /// - `rx_east`: Receiver local East unit vector
    pub fn update(
        &mut self,
        sat_pos: &Vector3<f64>,
        sun_pos: &Vector3<f64>,
        rx_pos: &Vector3<f64>,
        _rx_up: &Vector3<f64>,
        rx_north: &Vector3<f64>,
        rx_east: &Vector3<f64>,
    ) -> f64 {
        let los = sat_pos - rx_pos;
        let rho = vec_norm(&los);
        if rho < 1.0 {
            return self.prev_windup_rad;
        }
        let k_hat = los / rho; // Unit vector from receiver to satellite

        // Satellite body unit vectors
        let e_z = -sat_pos / vec_norm(sat_pos);
        let e_sun = sun_pos / vec_norm(sun_pos);
        let e_y_unnorm = e_z.cross(&e_sun);
        let e_y_len = vec_norm(&e_y_unnorm);
        if e_y_len < 1e-6 {
            return self.prev_windup_rad;
        }
        let e_y = e_y_unnorm / e_y_len;
        let e_x = e_y.cross(&e_z);

        // Transmitting effective dipole vector D' (Wu et al. 1993 eq. 8)
        let d_prime = e_x - k_hat * k_hat.dot(&e_x) - k_hat.cross(&e_y);

        // Receiving effective dipole vector D (Wu et al. 1993 eq. 9)
        let e_rx_x = rx_north; // Local north
        let e_rx_y = rx_east;  // Local east
        let d_rx = e_rx_x - k_hat * k_hat.dot(e_rx_x) + k_hat.cross(e_rx_y);

        let d_prime_len = vec_norm(&d_prime);
        let d_rx_len = vec_norm(&d_rx);
        if d_prime_len < 1e-6 || d_rx_len < 1e-6 {
            return self.prev_windup_rad;
        }

        let cos_phi = (d_prime.dot(&d_rx) / (d_prime_len * d_rx_len)).clamp(-1.0, 1.0);
        let cross_k = d_prime.cross(&d_rx);
        let sign = if k_hat.dot(&cross_k) >= 0.0 { 1.0 } else { -1.0 };
        let fractional_phi = sign * libm::acos(cos_phi);

        if !self.is_initialized {
            self.prev_windup_rad = fractional_phi;
            self.is_initialized = true;
            return fractional_phi;
        }

        // Continuous wrap integer tracking: \Delta \phi = \delta \phi + 2\pi * round((\phi_{prev} - \delta \phi)/2\pi)
        let delta = self.prev_windup_rad - fractional_phi;
        let two_pi = 2.0 * PI;
        let n_wrap = libm::round(delta / two_pi);
        let continuous_phi = fractional_phi + two_pi * n_wrap;

        self.prev_windup_rad = continuous_phi;
        continuous_phi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phase_windup_continuous_smoothness() {
        let mut tracker = PhaseWindupTracker::new();

        let sat = Vector3::new(0.0, 26_560_000.0, 0.0);
        let sun = Vector3::new(1.495e11, 0.0, 0.0);
        let rx = Vector3::new(6_378_137.0, 0.0, 0.0);
        let up = Vector3::new(1.0, 0.0, 0.0);
        let north = Vector3::new(0.0, 0.0, 1.0);
        let east = Vector3::new(0.0, 1.0, 0.0);

        let w1 = tracker.update(&sat, &sun, &rx, &up, &north, &east);
        assert!(!w1.is_nan());

        // Slight motion should maintain continuous angle without 2pi step
        let sat2 = Vector3::new(100.0, 26_560_000.0, 0.0);
        let w2 = tracker.update(&sat2, &sun, &rx, &up, &north, &east);
        assert!((w2 - w1).abs() < 0.1, "Windup step was discontinuous: {:.4} -> {:.4}", w1, w2);
    }
}
