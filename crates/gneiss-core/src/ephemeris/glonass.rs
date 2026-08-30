//! GLONASS broadcast ephemeris orbital model using 4th-order Runge-Kutta integration in PZ-90.

use crate::sat::SatelliteId;
use crate::time::GpsTime;
use nalgebra::Vector3;
use super::{J2_GLO, MU_GLO, OMEGA_E_GLO, RADIUS_GLO};

#[derive(Debug, Clone, PartialEq)]
pub struct GlonassEphemeris {
    pub sat: SatelliteId,
    pub toe: GpsTime,
    pub freq_num: i8,
    pub tau_n: f64,
    pub gamma_n: f64,
    pub delta_tau_n: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub ax: f64,
    pub ay: f64,
    pub az: f64,
}

fn glonass_derivatives(state: &[f64; 6], acc: &[f64; 3]) -> [f64; 6] {
    let r2 = state[0] * state[0] + state[1] * state[1] + state[2] * state[2];
    let r = libm::sqrt(r2);
    let r3 = r2 * r;
    let ae2 = RADIUS_GLO * RADIUS_GLO;
    let factor = 1.5 * J2_GLO * MU_GLO * ae2 / (r2 * r3);
    let z2_r2 = state[2] * state[2] / r2;

    let ax = -MU_GLO * state[0] / r3 - factor * state[0] * (1.0 - 5.0 * z2_r2)
        + OMEGA_E_GLO * OMEGA_E_GLO * state[0]
        + 2.0 * OMEGA_E_GLO * state[4]
        + acc[0];
    let ay = -MU_GLO * state[1] / r3 - factor * state[1] * (1.0 - 5.0 * z2_r2)
        + OMEGA_E_GLO * OMEGA_E_GLO * state[1]
        - 2.0 * OMEGA_E_GLO * state[3]
        + acc[1];
    let az = -MU_GLO * state[2] / r3 - factor * state[2] * (3.0 - 5.0 * z2_r2) + acc[2];

    [state[3], state[4], state[5], ax, ay, az]
}

fn rk4_step(state: &[f64; 6], acc: &[f64; 3], h: f64) -> [f64; 6] {
    let k1 = glonass_derivatives(state, acc);

    let mut s2 = [0.0; 6];
    for i in 0..6 {
        s2[i] = state[i] + 0.5 * h * k1[i];
    }
    let k2 = glonass_derivatives(&s2, acc);

    let mut s3 = [0.0; 6];
    for i in 0..6 {
        s3[i] = state[i] + 0.5 * h * k2[i];
    }
    let k3 = glonass_derivatives(&s3, acc);

    let mut s4 = [0.0; 6];
    for i in 0..6 {
        s4[i] = state[i] + h * k3[i];
    }
    let k4 = glonass_derivatives(&s4, acc);

    let mut next_state = [0.0; 6];
    for i in 0..6 {
        next_state[i] = state[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    next_state
}

impl GlonassEphemeris {
    pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
        let dt = t - self.toe;
        let mut state = [self.x, self.y, self.z, self.vx, self.vy, self.vz];
        let acc = [self.ax, self.ay, self.az];

        let step = if dt < 0.0 { -30.0 } else { 30.0 };
        let mut t_rem = dt;

        while libm::fabs(t_rem) > 1e-14 {
            let h = if libm::fabs(t_rem) < libm::fabs(step) {
                t_rem
            } else {
                step
            };
            state = rk4_step(&state, &acc, h);
            t_rem -= h;
        }

        let clk_err = self.tau_n + self.gamma_n * dt;
        let clk_drift = self.gamma_n;

        (
            Vector3::new(state[0], state[1], state[2]),
            Vector3::new(state[3], state[4], state[5]),
            clk_err,
            clk_drift,
        )
    }
}
