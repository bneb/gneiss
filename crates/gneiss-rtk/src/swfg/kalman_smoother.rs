//! GNSS Kalman Filter with Rauch-Tung-Striebel (RTS) smoother.
//!
//! State vector (8D):
//!   x[0..3] — position ECEF (m)
//!   x[3..6] — velocity ECEF (m/s)
//!   x[6]    — receiver clock bias (m)
//!   x[7]    — receiver clock drift (m/s)

use gneiss_core::time::GpsTime;
use nalgebra::{DMatrix, DVector, Vector3};

/// Per-epoch snapshot stored during the forward EKF pass.
#[derive(Clone)]
pub struct EkfSnapshot {
    pub time: GpsTime,
    pub x_predicted: DVector<f64>,
    pub p_predicted: DMatrix<f64>,
    pub x_updated: DVector<f64>,
    pub p_updated: DMatrix<f64>,
    pub n_obs: usize,
}

/// Smoothed state output from the RTS backward pass.
#[derive(Clone)]
pub struct SmoothedState {
    pub time: GpsTime,
    pub position_ecef: Vector3<f64>,
    pub velocity_ecef: Vector3<f64>,
    pub clock_bias_m: f64,
    pub p_smoothed: DMatrix<f64>,
    pub n_obs: usize,
}

/// Process noise tuning parameters.
#[derive(Debug, Clone)]
pub struct ProcessNoise {
    /// Velocity random walk (m²/s³). Controls how fast position
    /// uncertainty grows during GNSS outages.
    pub q_vel: f64,
    /// Clock bias random walk (m²/s).
    pub q_clk: f64,
    /// Clock drift random walk (m²/s³).
    pub q_drift: f64,
}

impl Default for ProcessNoise {
    fn default() -> Self {
        Self {
            q_vel: 1.0,
            q_clk: 100.0,
            q_drift: 1.0,
        }
    }
}

const STATE_DIM: usize = 8;

/// EKF + RTS smoother for GNSS-only post-processing.
pub struct GnssKalmanSmoother {
    pub x: DVector<f64>,
    pub p: DMatrix<f64>,
    pub history: Vec<EkfSnapshot>,
    pub process_noise: ProcessNoise,
    last_time: Option<GpsTime>,
}

impl GnssKalmanSmoother {
    pub fn new(initial_pos: Vector3<f64>, pn: ProcessNoise) -> Self {
        let mut x = DVector::zeros(STATE_DIM);
        x[0] = initial_pos.x;
        x[1] = initial_pos.y;
        x[2] = initial_pos.z;

        // Initial covariance: large for position/clock, moderate for velocity
        let mut p = DMatrix::zeros(STATE_DIM, STATE_DIM);
        for i in 0..3 { p[(i, i)] = 100.0 * 100.0; }      // 100m position
        for i in 3..6 { p[(i, i)] = 10.0 * 10.0; }         // 10 m/s velocity
        p[(6, 6)] = 1e6;                                     // clock bias
        p[(7, 7)] = 1e4;                                     // clock drift

        Self {
            x,
            p,
            history: Vec::new(),
            process_noise: pn,
            last_time: None,
        }
    }

    /// Build the state transition matrix F for a given dt.
    fn build_f(dt: f64) -> DMatrix<f64> {
        let mut f = DMatrix::identity(STATE_DIM, STATE_DIM);
        // Position += velocity * dt
        for i in 0..3 { f[(i, i + 3)] = dt; }
        // Clock bias += drift * dt
        f[(6, 7)] = dt;
        f
    }

    /// Build the process noise matrix Q for a given dt.
    fn build_q(&self, dt: f64) -> DMatrix<f64> {
        let mut q = DMatrix::zeros(STATE_DIM, STATE_DIM);
        let pn = &self.process_noise;

        // Velocity-driven position noise (integrated white noise on accel):
        //   q_pos = q_vel * dt³/3,  cross = q_vel * dt²/2
        let dt2 = dt * dt;
        let dt3 = dt2 * dt;
        for i in 0..3 {
            q[(i, i)] = pn.q_vel * dt3 / 3.0;
            q[(i, i + 3)] = pn.q_vel * dt2 / 2.0;
            q[(i + 3, i)] = pn.q_vel * dt2 / 2.0;
            q[(i + 3, i + 3)] = pn.q_vel * dt;
        }
        // Clock noise
        q[(6, 6)] = pn.q_clk * dt;
        q[(7, 7)] = pn.q_drift * dt;

        q
    }

    /// EKF predict step. Propagates state and covariance forward by dt.
    pub fn predict(&mut self, time: GpsTime) {
        let dt = self.last_time.map_or(1.0, |t| time.tow - t.tow);
        let dt = dt.max(0.01); // guard against zero/negative

        let f = Self::build_f(dt);
        let q = self.build_q(dt);

        self.x = &f * &self.x;
        self.p = &f * &self.p * f.transpose() + q;

        self.last_time = Some(time);
    }

    /// EKF update step with pseudorange observations.
    ///
    /// `observations`: Vec of (sat_pos_ecef, pseudorange_corrected, variance_m2)
    /// where pseudorange_corrected = raw_pr - sat_clock + tropo + iono
    pub fn update(
        &mut self,
        time: GpsTime,
        observations: &[(Vector3<f64>, f64, f64)],
    ) {
        let n_obs = observations.len();

        // Store predicted state before update
        let x_pred = self.x.clone();
        let p_pred = self.p.clone();

        if n_obs >= 4 {
            let mut h = DMatrix::zeros(n_obs, STATE_DIM);
            let mut z = DVector::zeros(n_obs);
            let mut r = DMatrix::zeros(n_obs, n_obs);

            let rx_pos = Vector3::new(self.x[0], self.x[1], self.x[2]);

            for (i, (sat_pos, pr_corrected, var)) in observations.iter().enumerate() {
                let diff = rx_pos - sat_pos;
                let range = diff.norm().max(1.0);
                let e = diff / range; // unit vector rx → sat

                // Predicted measurement: geometric range + clock bias
                let h_pred = range + self.x[6];

                // Innovation
                z[i] = pr_corrected - h_pred;

                // Jacobian row: d(range)/d(pos) = e^T, d/d(vel) = 0, d/d(clk) = 1, d/d(drift) = 0
                h[(i, 0)] = e.x;
                h[(i, 1)] = e.y;
                h[(i, 2)] = e.z;
                h[(i, 6)] = 1.0;

                r[(i, i)] = *var;
            }

            // Kalman gain
            let h_t = h.transpose();
            let s = &h * &self.p * &h_t + &r;
            if let Some(s_inv) = s.try_inverse() {
                let k = &self.p * &h_t * &s_inv;

                // State update
                self.x = &self.x + &k * &z;

                // Joseph-form covariance update: P = (I - KH) P (I - KH)^T + K R K^T
                let i_kh = DMatrix::identity(STATE_DIM, STATE_DIM) - &k * &h;
                self.p = &i_kh * &self.p * i_kh.transpose() + &k * &r * k.transpose();
            }
        }

        self.history.push(EkfSnapshot {
            time,
            x_predicted: x_pred,
            p_predicted: p_pred,
            x_updated: self.x.clone(),
            p_updated: self.p.clone(),
            n_obs,
        });
    }

    /// RTS backward smoother. Returns smoothed states for all epochs.
    pub fn smooth(&self) -> Vec<SmoothedState> {
        let n = self.history.len();
        if n == 0 { return Vec::new(); }

        let mut smoothed_x = vec![DVector::zeros(STATE_DIM); n];
        let mut smoothed_p = vec![DMatrix::zeros(STATE_DIM, STATE_DIM); n];

        // Initialize from last epoch's filtered state
        smoothed_x[n - 1] = self.history[n - 1].x_updated.clone();
        smoothed_p[n - 1] = self.history[n - 1].p_updated.clone();

        // Backward sweep
        for k in (0..n - 1).rev() {
            let p_k = &self.history[k].p_updated;
            let p_pred_next = &self.history[k + 1].p_predicted;

            // Smoother gain: G_k = P_{k|k} F^T P_{k+1|k}^{-1}
            let dt = self.history[k + 1].time.tow - self.history[k].time.tow;
            let dt = dt.max(0.01);
            let f = Self::build_f(dt);

            let p_pred_inv = match p_pred_next.clone().try_inverse() {
                Some(inv) => inv,
                None => {
                    // Fallback: use filtered state if inversion fails
                    smoothed_x[k] = self.history[k].x_updated.clone();
                    smoothed_p[k] = self.history[k].p_updated.clone();
                    continue;
                }
            };

            let g = p_k * f.transpose() * &p_pred_inv;

            // Smoothed state
            let x_k = &self.history[k].x_updated;
            let x_pred_next = &self.history[k + 1].x_predicted;
            smoothed_x[k] = x_k + &g * (&smoothed_x[k + 1] - x_pred_next);

            // Smoothed covariance
            smoothed_p[k] = p_k + &g * (&smoothed_p[k + 1] - p_pred_next) * g.transpose();
        }

        // Convert to output
        smoothed_x
            .into_iter()
            .zip(smoothed_p)
            .enumerate()
            .map(|(i, (x, p))| SmoothedState {
                time: self.history[i].time,
                position_ecef: Vector3::new(x[0], x[1], x[2]),
                velocity_ecef: Vector3::new(x[3], x[4], x[5]),
                clock_bias_m: x[6],
                p_smoothed: p,
                n_obs: self.history[i].n_obs,
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_smoother() -> GnssKalmanSmoother {
        let pos = Vector3::new(1e6, 2e6, 3e6);
        GnssKalmanSmoother::new(pos, ProcessNoise::default())
    }

    #[test]
    fn test_predict_constant_velocity() {
        let mut kf = make_smoother();
        // Set velocity to (1, 2, 3) m/s
        kf.x[3] = 1.0;
        kf.x[4] = 2.0;
        kf.x[5] = 3.0;
        kf.last_time = Some(GpsTime::new(2000, 100.0));

        kf.predict(GpsTime::new(2000, 101.0)); // dt = 1s

        assert!((kf.x[0] - 1e6 - 1.0).abs() < 1e-10);
        assert!((kf.x[1] - 2e6 - 2.0).abs() < 1e-10);
        assert!((kf.x[2] - 3e6 - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_predict_covariance_growth() {
        let mut kf = make_smoother();
        let p_before = kf.p.clone();
        kf.predict(GpsTime::new(2000, 101.0));

        // Position covariance should grow due to process noise
        for i in 0..3 {
            assert!(kf.p[(i, i)] > p_before[(i, i)]);
        }
    }

    #[test]
    fn test_update_reduces_covariance() {
        let mut kf = make_smoother();
        kf.predict(GpsTime::new(2000, 100.0));

        let p_before = kf.p.clone();

        // Simulate 6 satellites at known positions
        let sats: Vec<(Vector3<f64>, f64, f64)> = (0..6)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::FRAC_PI_3;
                let sat_pos = Vector3::new(
                    1e6 + 2e7 * angle.cos(),
                    2e6 + 2e7 * angle.sin(),
                    3e6 + 1e7,
                );
                let range = (sat_pos - Vector3::new(1e6, 2e6, 3e6)).norm();
                (sat_pos, range, 4.0) // 2m sigma
            })
            .collect();

        kf.update(GpsTime::new(2000, 100.0), &sats);

        // Position covariance should decrease after update
        for i in 0..3 {
            assert!(
                kf.p[(i, i)] < p_before[(i, i)],
                "P[{i},{i}] should decrease: {} < {}",
                kf.p[(i, i)],
                p_before[(i, i)]
            );
        }
    }

    #[test]
    fn test_joseph_form_symmetry() {
        let mut kf = make_smoother();
        kf.predict(GpsTime::new(2000, 100.0));

        let sats: Vec<(Vector3<f64>, f64, f64)> = (0..6)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::FRAC_PI_3;
                let sat_pos = Vector3::new(
                    1e6 + 2e7 * angle.cos(),
                    2e6 + 2e7 * angle.sin(),
                    3e6 + 1e7,
                );
                let range = (sat_pos - Vector3::new(1e6, 2e6, 3e6)).norm();
                (sat_pos, range, 4.0)
            })
            .collect();

        kf.update(GpsTime::new(2000, 100.0), &sats);

        // P should remain symmetric
        let diff = (&kf.p - kf.p.transpose()).norm();
        assert!(diff < 1e-10, "P not symmetric: asymmetry = {diff}");
    }

    #[test]
    fn test_rts_smooth_reduces_uncertainty() {
        let mut kf = make_smoother();
        let _t0 = GpsTime::new(2000, 100.0);

        // Run 5 epochs with observations
        for i in 0..5 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            kf.predict(t);

            let sats: Vec<(Vector3<f64>, f64, f64)> = (0..6)
                .map(|j| {
                    let angle = j as f64 * std::f64::consts::FRAC_PI_3;
                    let sat_pos = Vector3::new(
                        1e6 + 2e7 * angle.cos(),
                        2e6 + 2e7 * angle.sin(),
                        3e6 + 1e7,
                    );
                    let range = (sat_pos - Vector3::new(1e6, 2e6, 3e6)).norm();
                    (sat_pos, range, 4.0)
                })
                .collect();

            kf.update(t, &sats);
        }

        let smoothed = kf.smooth();
        assert_eq!(smoothed.len(), 5);

        // First epoch's smoothed covariance should be smaller than filtered
        // (RTS uses future information)
        let filtered_p0_trace = kf.history[0].p_updated.trace();
        let smoothed_p0_trace = smoothed[0].p_smoothed.trace();
        assert!(
            smoothed_p0_trace <= filtered_p0_trace + 1e-6,
            "Smoothed trace {} should be <= filtered trace {}",
            smoothed_p0_trace,
            filtered_p0_trace
        );
    }
}
