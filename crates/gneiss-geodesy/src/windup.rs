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

    /// Computes continuous phase windup correction in radians with optional receiver attitude.
    pub fn update_with_attitude(
        &mut self,
        sat_pos: &Vector3<f64>,
        sun_pos: &Vector3<f64>,
        rx_pos: &Vector3<f64>,
        rx_north: &Vector3<f64>,
        rx_east: &Vector3<f64>,
        r_body_to_ecef: Option<&nalgebra::Matrix3<f64>>,
    ) -> f64 {
        let los = sat_pos - rx_pos;
        let rho = vec_norm(&los);
        if rho < 1.0 {
            return self.prev_windup_rad;
        }
        let k_hat = los / rho;

        let e_z = -sat_pos / vec_norm(sat_pos);
        let e_sun = sun_pos / vec_norm(sun_pos);
        let e_y_unnorm = e_z.cross(&e_sun);
        let e_y_len = vec_norm(&e_y_unnorm);
        if e_y_len < 1e-6 {
            return self.prev_windup_rad;
        }
        let e_y = e_y_unnorm / e_y_len;
        let e_x = e_y.cross(&e_z);
        let d_prime = e_x - k_hat * k_hat.dot(&e_x) - k_hat.cross(&e_y);

        let (e_rx_x, e_rx_y) = if let Some(r_b2e) = r_body_to_ecef {
            (r_b2e * rx_north, r_b2e * rx_east)
        } else {
            (*rx_north, *rx_east)
        };
        let d_rx = e_rx_x - k_hat * k_hat.dot(&e_rx_x) + k_hat.cross(&e_rx_y);

        let (dp_len, dr_len) = (vec_norm(&d_prime), vec_norm(&d_rx));
        if dp_len < 1e-6 || dr_len < 1e-6 {
            return self.prev_windup_rad;
        }

        let cos_phi = (d_prime.dot(&d_rx) / (dp_len * dr_len)).clamp(-1.0, 1.0);
        let sign = if k_hat.dot(&d_prime.cross(&d_rx)) >= 0.0 { 1.0 } else { -1.0 };
        let frac_phi = sign * libm::acos(cos_phi);

        if !self.is_initialized {
            self.prev_windup_rad = frac_phi;
            self.is_initialized = true;
            return frac_phi;
        }

        let delta = self.prev_windup_rad - frac_phi;
        let n_wrap = libm::round(delta / (2.0 * PI));
        let continuous_phi = frac_phi + 2.0 * PI * n_wrap;
        self.prev_windup_rad = continuous_phi;
        continuous_phi
    }

    /// Computes continuous phase windup correction in radians for a static/topocentric receiver.
    pub fn update(
        &mut self,
        sat_pos: &Vector3<f64>,
        sun_pos: &Vector3<f64>,
        rx_pos: &Vector3<f64>,
        _rx_up: &Vector3<f64>,
        rx_north: &Vector3<f64>,
        rx_east: &Vector3<f64>,
    ) -> f64 {
        self.update_with_attitude(sat_pos, sun_pos, rx_pos, rx_north, rx_east, None)
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

    #[test]
    fn test_phase_windup_receiver_attitude_rotation() {
        let mut tracker = PhaseWindupTracker::new();
        let sat = Vector3::new(10_000_000.0, 15_000_000.0, 20_000_000.0);
        let sun = Vector3::new(1.495e11, 0.0, 0.0);
        let rx = Vector3::new(0.0, 0.0, 6_378_137.0);
        let up = Vector3::new(0.0, 0.0, 1.0);
        let north = Vector3::new(1.0, 0.0, 0.0);
        let east = Vector3::new(0.0, 1.0, 0.0);

        let w0 = tracker.update(&sat, &sun, &rx, &up, &north, &east);
        assert!(!w0.is_nan());

        // 90-degree yaw rotation around up axis (z)
        let angle = core::f64::consts::FRAC_PI_2;
        let r_yaw = nalgebra::Matrix3::new(
            angle.cos(), -angle.sin(), 0.0,
            angle.sin(),  angle.cos(), 0.0,
            0.0,          0.0,         1.0,
        );

        let mut tracker_yaw = PhaseWindupTracker::new();
        let w_yaw = tracker_yaw.update_with_attitude(&sat, &sun, &rx, &north, &east, Some(&r_yaw));
        let diff = (w_yaw - w0).abs();
        assert!((diff - angle).abs() < 1e-4, "90-degree receiver yaw must induce pi/2 windup, got diff: {:.4}", diff);
    }
}
