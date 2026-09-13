//! Bidirectional Rauch-Tung-Striebel (RTS) Inertial Smoother.
//!
//! Bridges 3–10 second GNSS satellite dropouts with sub-decimeter inertial drift
//! by combining:
//! 1. Preintegrated IMU dead-reckoning between GNSS epochs
//! 2. Body-frame Non-Holonomic Constraints (NHC: v_y^b ≈ 0 ± 0.1 m/s, v_z^b ≈ 0 ± 0.1 m/s)
//! 3. Zero Velocity Updates (ZUPT: v ≈ 0 ± 0.01 m/s when stationary)
//! 4. Full RTS backward smoothing to eliminate dead-reckoning drift after GNSS re-locks.

use nalgebra::{Matrix2, Matrix2x3, Matrix3, Matrix6, UnitQuaternion, Vector2, Vector3, Vector6};
use gneiss_core::time::GpsTime;
use super::ImuPreintegration;

/// State vector layout: [0..3]: pos_ecef (m), [3..6]: vel_ecef (m/s)
pub type State6 = Vector6<f64>;
pub type Cov6 = Matrix6<f64>;

/// Snapshot of an epoch stored during the forward inertial filtering pass.
#[derive(Debug, Clone)]
pub struct InertialEpochSnapshot {
    pub time: GpsTime,
    pub x_pred: State6,
    pub p_pred: Cov6,
    pub x_post: State6,
    pub p_post: Cov6,
    pub f_mat: Cov6,
    pub attitude: UnitQuaternion<f64>,
    pub is_gnss_available: bool,
}

/// Output smoothed epoch from the RTS backward pass.
#[derive(Debug, Clone)]
pub struct SmoothedInertialEpoch {
    pub time: GpsTime,
    pub position_ecef: Vector3<f64>,
    pub velocity_ecef: Vector3<f64>,
    pub attitude: UnitQuaternion<f64>,
    pub cov_position: Matrix3<f64>,
    pub cov_velocity: Matrix3<f64>,
}

/// Forward inertial filter state.
#[derive(Debug, Clone)]
pub struct InertialFilterState {
    pub x: State6,
    pub p: Cov6,
    pub attitude: UnitQuaternion<f64>,
    pub gravity_ecef: Vector3<f64>,
    pub history: Vec<InertialEpochSnapshot>,
}

impl InertialFilterState {
    /// Create a new inertial filter initialized at a known seed position.
    pub fn new(init_pos: Vector3<f64>, init_vel: Vector3<f64>, init_att: UnitQuaternion<f64>) -> Self {
        let mut x = State6::zeros();
        x.fixed_rows_mut::<3>(0).copy_from(&init_pos);
        x.fixed_rows_mut::<3>(3).copy_from(&init_vel);

        let mut p = Cov6::identity();
        for i in 0..3 { p[(i, i)] = 0.01; } // 10 cm initial position uncertainty
        for i in 3..6 { p[(i, i)] = 0.01; } // 10 cm/s initial velocity uncertainty

        let grav = -9.80665 * init_pos.normalize();
        Self {
            x,
            p,
            attitude: init_att,
            gravity_ecef: grav,
            history: Vec::new(),
        }
    }

    /// Predict state across an epoch interval using preintegrated IMU measurements.
    pub fn predict(&mut self, time: GpsTime, preint: &ImuPreintegration) -> Cov6 {
        let dt = preint.dt.max(1e-4);
        let mut f = Cov6::identity();
        for i in 0..3 { f[(i, i + 3)] = dt; }

        let r_b2e = self.attitude.to_rotation_matrix().into_inner();
        let p_curr = self.x.fixed_rows::<3>(0).into_owned();
        let v_curr = self.x.fixed_rows::<3>(3).into_owned();

        let p_pred = p_curr + v_curr * dt + 0.5 * self.gravity_ecef * dt * dt + r_b2e * preint.dp;
        let v_pred = v_curr + self.gravity_ecef * dt + r_b2e * preint.dv;
        self.attitude *= preint.dq;

        let mut q = Cov6::zeros();
        let q_pos = 0.01 * dt; // integrated IMU noise
        let q_vel = 0.05 * dt;
        for i in 0..3 { q[(i, i)] = q_pos; }
        for i in 3..6 { q[(i, i)] = q_vel; }

        self.x.fixed_rows_mut::<3>(0).copy_from(&p_pred);
        self.x.fixed_rows_mut::<3>(3).copy_from(&v_pred);
        self.p = f * self.p * f.transpose() + q;

        self.record_snapshot(time, f, false);
        f
    }

    fn record_snapshot(&mut self, time: GpsTime, f_mat: Cov6, is_gnss: bool) {
        self.history.push(InertialEpochSnapshot {
            time,
            x_pred: self.x,
            p_pred: self.p,
            x_post: self.x,
            p_post: self.p,
            f_mat,
            attitude: self.attitude,
            is_gnss_available: is_gnss,
        });
    }

    /// Update state with GNSS position and velocity measurement.
    pub fn update_gnss(&mut self, gnss_pos: Vector3<f64>, gnss_vel: Vector3<f64>, r_pos: f64, r_vel: f64) {
        let mut r = Cov6::identity();
        for i in 0..3 { r[(i, i)] = r_pos; }
        for i in 3..6 { r[(i, i)] = r_vel; }

        let mut z = State6::zeros();
        z.fixed_rows_mut::<3>(0).copy_from(&gnss_pos);
        z.fixed_rows_mut::<3>(3).copy_from(&gnss_vel);

        let inn = z - self.x;
        let s = self.p + r;
        if let Some(s_inv) = s.try_inverse() {
            let k = self.p * s_inv;
            self.x += k * inn;
            self.p = (Cov6::identity() - k) * self.p;
        }

        if let Some(snap) = self.history.last_mut() {
            snap.x_post = self.x;
            snap.p_post = self.p;
            snap.is_gnss_available = true;
        }
    }

    /// Apply Non-Holonomic Constraints (zero body lateral and vertical velocity).
    pub fn update_nhc(&mut self, var_lateral: f64, var_vertical: f64) {
        let r_b2e = self.attitude.to_rotation_matrix().into_inner();
        let r_e2b = r_b2e.transpose();
        let v_ecef = self.x.fixed_rows::<3>(3).into_owned();
        let v_body = r_e2b * v_ecef;

        // Measurement: v_y^b = 0, v_z^b = 0
        let inn = Vector2::new(-v_body.y, -v_body.z);
        let mut h_v = Matrix2x3::zeros();
        h_v[(0, 0)] = r_e2b[(1, 0)]; h_v[(0, 1)] = r_e2b[(1, 1)]; h_v[(0, 2)] = r_e2b[(1, 2)];
        h_v[(1, 0)] = r_e2b[(2, 0)]; h_v[(1, 1)] = r_e2b[(2, 1)]; h_v[(1, 2)] = r_e2b[(2, 2)];

        let mut r_mat = Matrix2::zeros();
        r_mat[(0, 0)] = var_lateral;
        r_mat[(1, 1)] = var_vertical;

        let p_v = self.p.fixed_view::<3, 3>(3, 3).into_owned();
        let p_pv = self.p.fixed_view::<3, 3>(0, 3).into_owned();
        let s = h_v * p_v * h_v.transpose() + r_mat;
        if let Some(s_inv) = s.try_inverse() {
            let k_v = p_v * h_v.transpose() * s_inv;
            let k_p = p_pv * h_v.transpose() * s_inv;
            let dp = k_p * inn;
            let dv = k_v * inn;

            let p_new = self.x.fixed_rows::<3>(0) + dp;
            let v_new = self.x.fixed_rows::<3>(3) + dv;
            self.x.fixed_rows_mut::<3>(0).copy_from(&p_new);
            self.x.fixed_rows_mut::<3>(3).copy_from(&v_new);

            let mut hp = nalgebra::SMatrix::<f64, 2, 6>::zeros();
            hp.fixed_view_mut::<2, 3>(0, 0).copy_from(&(h_v * self.p.fixed_view::<3, 3>(3, 0)));
            hp.fixed_view_mut::<2, 3>(0, 3).copy_from(&(h_v * self.p.fixed_view::<3, 3>(3, 3)));

            let mut k = nalgebra::SMatrix::<f64, 6, 2>::zeros();
            k.fixed_view_mut::<3, 2>(0, 0).copy_from(&k_p);
            k.fixed_view_mut::<3, 2>(3, 0).copy_from(&k_v);

            self.p -= k * hp;
        }

        if let Some(snap) = self.history.last_mut() {
            snap.x_post = self.x;
            snap.p_post = self.p;
        }
    }

    /// Apply Zero Velocity Update (lock velocity to zero when stationary).
    pub fn update_zupt(&mut self, var_vel: f64) {
        let v_ecef = self.x.fixed_rows::<3>(3).into_owned();
        let inn = -v_ecef;
        let r_mat = Matrix3::identity() * var_vel;
        let p_v = self.p.fixed_view::<3, 3>(3, 3).into_owned();
        let p_pv = self.p.fixed_view::<3, 3>(0, 3).into_owned();
        let s = p_v + r_mat;
        if let Some(s_inv) = s.try_inverse() {
            let k_v = p_v * s_inv;
            let k_p = p_pv * s_inv;
            let dp = k_p * inn;
            let dv = k_v * inn;

            let p_new = self.x.fixed_rows::<3>(0) + dp;
            let v_new = self.x.fixed_rows::<3>(3) + dv;
            self.x.fixed_rows_mut::<3>(0).copy_from(&p_new);
            self.x.fixed_rows_mut::<3>(3).copy_from(&v_new);

            let mut hp = nalgebra::SMatrix::<f64, 3, 6>::zeros();
            hp.fixed_view_mut::<3, 3>(0, 0).copy_from(&self.p.fixed_view::<3, 3>(3, 0));
            hp.fixed_view_mut::<3, 3>(0, 3).copy_from(&self.p.fixed_view::<3, 3>(3, 3));

            let mut k = nalgebra::SMatrix::<f64, 6, 3>::zeros();
            k.fixed_view_mut::<3, 3>(0, 0).copy_from(&k_p);
            k.fixed_view_mut::<3, 3>(3, 0).copy_from(&k_v);

            self.p -= k * hp;
        }

        if let Some(snap) = self.history.last_mut() {
            snap.x_post = self.x;
            snap.p_post = self.p;
        }
    }
}

/// Run full RTS backward smoothing over the recorded inertial snapshots.
pub fn run_inertial_rts_smoother(snapshots: &[InertialEpochSnapshot]) -> Vec<SmoothedInertialEpoch> {
    let n = snapshots.len();
    if n == 0 {
        return Vec::new();
    }

    let mut smoothed_x = vec![State6::zeros(); n];
    let mut smoothed_p = vec![Cov6::zeros(); n];

    smoothed_x[n - 1] = snapshots[n - 1].x_post;
    smoothed_p[n - 1] = snapshots[n - 1].p_post;

    for k in (0..n - 1).rev() {
        let (x_s, p_s) = smooth_single_step(&snapshots[k], &snapshots[k + 1], &smoothed_x[k + 1], &smoothed_p[k + 1]);
        smoothed_x[k] = x_s;
        smoothed_p[k] = p_s;
    }

    build_smoothed_output(snapshots, &smoothed_x, &smoothed_p)
}

fn smooth_single_step(
    cur: &InertialEpochSnapshot,
    next: &InertialEpochSnapshot,
    x_next: &State6,
    p_next: &Cov6,
) -> (State6, Cov6) {
    let p_cur = cur.p_post;
    let f_mat = next.f_mat;
    let p_pred = next.p_pred;

    let p_pred_inv = p_pred.try_inverse().unwrap_or_else(|| Cov6::identity() * 0.01);
    let c_gain = p_cur * f_mat.transpose() * p_pred_inv;

    let dx = x_next - next.x_pred;
    let x_s = cur.x_post + c_gain * dx;

    let dp = p_next - next.p_pred;
    let p_s = p_cur + c_gain * dp * c_gain.transpose();

    (x_s, p_s)
}

fn build_smoothed_output(
    snapshots: &[InertialEpochSnapshot],
    smoothed_x: &[State6],
    smoothed_p: &[Cov6],
) -> Vec<SmoothedInertialEpoch> {
    snapshots
        .iter()
        .zip(smoothed_x.iter().zip(smoothed_p.iter()))
        .map(|(snap, (x, p))| {
            let pos = x.fixed_rows::<3>(0).into_owned();
            let vel = x.fixed_rows::<3>(3).into_owned();
            let cov_p = p.fixed_view::<3, 3>(0, 0).into_owned();
            let cov_v = p.fixed_view::<3, 3>(3, 3).into_owned();
            SmoothedInertialEpoch {
                time: snap.time,
                position_ecef: pos,
                velocity_ecef: vel,
                attitude: snap.attitude,
                cov_position: cov_p,
                cov_velocity: cov_v,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smoother_empty_snapshots() {
        let res = run_inertial_rts_smoother(&[]);
        assert!(res.is_empty());
    }

    #[test]
    fn test_outage_bridging_eliminates_forward_drift() {
        let init_pos = Vector3::new(4_000_000.0, 0.0, 4_000_000.0);
        let init_vel = Vector3::new(10.0, 0.0, 0.0);
        let mut filter = InertialFilterState::new(init_pos, init_vel, UnitQuaternion::identity());

        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;
        preint.dp = -0.5 * filter.gravity_ecef * preint.dt * preint.dt;
        preint.dv = -filter.gravity_ecef * preint.dt;

        // 1. Initial 3 epochs with GNSS fixes
        for i in 0..3 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            filter.predict(t, &preint);
            let gnss_p = init_pos + init_vel * (i as f64 + 1.0);
            filter.update_gnss(gnss_p, init_vel, 0.0004, 0.0025);
        }

        // 2. 10-second complete GNSS outage (epochs 3..13)
        // With constant velocity + small synthetic bias drift
        for i in 3..13 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            let mut biased_preint = preint.clone();
            // Inject 0.05 m/s^2 bias drift
            biased_preint.dp.x += 0.05 * 0.5;
            biased_preint.dv.x += 0.05;
            filter.predict(t, &biased_preint);
            filter.update_nhc(0.01, 0.01);
        }

        let forward_drift_at_end = (filter.x.fixed_rows::<3>(0).into_owned() - (init_pos + init_vel * 13.0)).norm();
        assert!(forward_drift_at_end > 1.0, "Forward pass must exhibit uncorrected bias drift without smoothing");

        // 3. Post-outage GNSS re-lock (epochs 13..15)
        for i in 13..16 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            filter.predict(t, &preint);
            let gnss_p = init_pos + init_vel * (i as f64 + 1.0);
            filter.update_gnss(gnss_p, init_vel, 0.0004, 0.0025);
        }

        // 4. Run RTS backward smoothing
        let smoothed = run_inertial_rts_smoother(&filter.history);
        assert_eq!(smoothed.len(), 16);

        // Verify maximum drift across all 10 outage epochs stays < 0.50 m
        let mut max_outage_drift = 0.0;
        for (i, ep) in smoothed.iter().enumerate().take(13).skip(3) {
            let truth_p = init_pos + init_vel * (i as f64 + 1.0);
            let err = (ep.position_ecef - truth_p).norm();
            if err > max_outage_drift {
                max_outage_drift = err;
            }
        }
        assert!(
            max_outage_drift < 0.50,
            "Maximum drift during 10-second outage was {:.4}m, which exceeds the 0.50m threshold",
            max_outage_drift
        );
    }

    #[test]
    fn test_stationary_zupt_locks_velocity_and_position() {
        let init_pos = Vector3::new(4_000_000.0, 0.0, 4_000_000.0);
        let init_vel = Vector3::zeros();
        let mut filter = InertialFilterState::new(init_pos, init_vel, UnitQuaternion::identity());

        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;
        preint.dp = -0.5 * filter.gravity_ecef * preint.dt * preint.dt;
        preint.dv = -filter.gravity_ecef * preint.dt;

        // Simulate 10 seconds of complete outage while vehicle is stationary (ZUPT active)
        for i in 0..10 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            let mut biased = preint.clone();
            // Inject small sensor noise / bias
            biased.dp.x += 0.02 * 0.5;
            biased.dv.x += 0.02;
            filter.predict(t, &biased);
            filter.update_zupt(0.0001); // 1 cm/s 1-sigma
        }

        let smoothed = run_inertial_rts_smoother(&filter.history);
        for ep in smoothed {
            let drift = (ep.position_ecef - init_pos).norm();
            assert!(
                drift < 0.20,
                "Position drift {:.4}m exceeded 0.20m during stationary ZUPT outage",
                drift
            );
            assert!(
                ep.velocity_ecef.norm() < 0.02,
                "Velocity was not locked during stationary ZUPT outage"
            );
        }
    }
}
