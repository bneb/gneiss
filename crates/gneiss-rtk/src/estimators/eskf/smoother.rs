use gneiss_core::time::GpsTime;
use nalgebra::UnitQuaternion;

use super::types::{clamp_vector, EngineError, EskfSnapshot, EskfState, Matrix15, Vector15};

#[derive(Debug, Clone, Default)]
pub struct EskfSmoother {
    snapshots: Vec<EskfSnapshot>,
}

impl EskfSmoother {
    pub fn new() -> Self {
        Self {
            snapshots: Vec::new(),
        }
    }

    pub fn push(&mut self, snapshot: EskfSnapshot) {
        self.snapshots.push(snapshot);
    }

    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    pub fn clear(&mut self) {
        self.snapshots.clear();
    }

    pub fn snapshots(&self) -> &[EskfSnapshot] {
        &self.snapshots
    }

    pub fn smooth(&self) -> Result<Vec<EskfState>, EngineError> {
        if self.snapshots.is_empty() {
            return Err(EngineError::EmptyHistory);
        }
        let n = self.snapshots.len();
        let mut smoothed = vec![self.snapshots[n - 1].state_post.clone(); n];

        for k in (0..n - 1).rev() {
            let next_snap = &self.snapshots[k + 1];
            let curr_snap = &self.snapshots[k];

            let c_k = compute_smoother_gain(
                &curr_snap.state_post.cov,
                &next_snap.phi,
                &next_snap.state_pred.cov,
            )?;
            let dx_next = compute_state_discrepancy(&smoothed[k + 1], &next_snap.state_pred);

            smoothed[k] = apply_smoother_correction(
                &curr_snap.state_post,
                &c_k,
                &dx_next,
                &next_snap.state_pred.cov,
                &smoothed[k + 1].cov,
            );
        }
        Ok(smoothed)
    }

    pub fn smooth_with_time(&self) -> Result<Vec<(GpsTime, EskfState)>, EngineError> {
        let states = self.smooth()?;
        let epochs = self
            .snapshots
            .iter()
            .zip(states)
            .map(|(snap, state)| (snap.time, state))
            .collect();
        Ok(epochs)
    }
}

fn compute_smoother_gain(
    p_post: &Matrix15<f64>,
    phi_next: &Matrix15<f64>,
    p_pred_next: &Matrix15<f64>,
) -> Result<Matrix15<f64>, EngineError> {
    let mut reg = *p_pred_next;
    for i in 0..15 {
        reg[(i, i)] = reg[(i, i)] * 1.0001 + 1e-12;
    }
    let p_pred_inv = reg
        .try_inverse()
        .ok_or(EngineError::InversionError)?;
    Ok(p_post * phi_next.transpose() * p_pred_inv)
}

fn compute_state_discrepancy(
    smoothed_next: &EskfState,
    pred_next: &EskfState,
) -> Vector15<f64> {
    let mut dx = Vector15::zeros();
    dx.fixed_rows_mut::<3>(0).copy_from(&(smoothed_next.pos_ecef - pred_next.pos_ecef));
    dx.fixed_rows_mut::<3>(3).copy_from(&(smoothed_next.vel_ecef - pred_next.vel_ecef));

    let dq = smoothed_next.attitude * pred_next.attitude.inverse();
    dx.fixed_rows_mut::<3>(6).copy_from(&dq.scaled_axis());
    dx.fixed_rows_mut::<3>(9).copy_from(&(smoothed_next.accel_bias - pred_next.accel_bias));
    dx.fixed_rows_mut::<3>(12).copy_from(&(smoothed_next.gyro_bias - pred_next.gyro_bias));
    dx
}

fn apply_smoother_correction(
    post_k: &EskfState,
    c_k: &Matrix15<f64>,
    dx_next: &Vector15<f64>,
    p_pred_next: &Matrix15<f64>,
    p_smooth_next: &Matrix15<f64>,
) -> EskfState {
    let mut s_k = post_k.clone();
    let dx_k = c_k * dx_next;

    s_k.pos_ecef += dx_k.fixed_rows::<3>(0);
    s_k.vel_ecef += dx_k.fixed_rows::<3>(3);

    let dq = UnitQuaternion::from_scaled_axis(dx_k.fixed_rows::<3>(6).into_owned());
    s_k.attitude = UnitQuaternion::new_normalize((dq * s_k.attitude).into_inner());

    s_k.accel_bias += dx_k.fixed_rows::<3>(9);
    s_k.gyro_bias += dx_k.fixed_rows::<3>(12);
    clamp_vector(&mut s_k.accel_bias, 0.5);
    clamp_vector(&mut s_k.gyro_bias, 0.05);

    let dp = p_smooth_next - p_pred_next;
    let new_cov = post_k.cov + c_k * dp * c_k.transpose();
    s_k.cov = 0.5 * (new_cov + new_cov.transpose());
    s_k
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_smoother_empty_returns_error() {
        let smoother = EskfSmoother::new();
        assert_eq!(smoother.smooth(), Err(EngineError::EmptyHistory));
    }

    #[test]
    fn test_smoother_single_epoch() {
        let mut smoother = EskfSmoother::new();
        let state = EskfState::new(Vector3::new(1.0, 2.0, 3.0), Vector3::zeros(), UnitQuaternion::identity());
        let snap = EskfSnapshot {
            time: GpsTime::new(2000, 100.0),
            state_pred: state.clone(),
            state_post: state.clone(),
            phi: Matrix15::identity(),
            is_gnss_available: true,
        };
        smoother.push(snap);
        let smoothed = smoother.smooth().expect("smoothing failed");
        assert_eq!(smoothed.len(), 1);
        assert_eq!(smoothed[0].pos_ecef, state.pos_ecef);
    }

    fn make_snapshot(time_s: f64, p_pred: f64, p_post: f64, gnss: bool) -> EskfSnapshot {
        let mut s_pred = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
        s_pred.cov[(0, 0)] = p_pred;
        let mut s_post = s_pred.clone();
        s_post.cov[(0, 0)] = p_post;
        EskfSnapshot {
            time: GpsTime::new(2000, time_s),
            state_pred: s_pred,
            state_post: s_post,
            phi: Matrix15::identity(),
            is_gnss_available: gnss,
        }
    }

    #[test]
    fn test_smoother_reduces_intermediate_covariance() {
        let mut smoother = EskfSmoother::new();
        smoother.push(make_snapshot(0.0, 1.0, 1.0, true));
        smoother.push(make_snapshot(1.0, 5.0, 5.0, false));
        smoother.push(make_snapshot(2.0, 6.0, 0.05, true));

        let smoothed = smoother.smooth().expect("smooth failed");
        assert!(smoothed[1].cov[(0, 0)] < 5.0, "Smoother must reduce covariance during outage");
    }

    /// Hand evaluation of one Rauch-Tung-Striebel backward pass.
    ///
    /// With two snapshots and phi = I the smoother computes, for k = 0,
    ///   C_0      = P_post_0 * P_pred_1^-1
    ///   dx_1     = x_post_1 - x_pred_1            (already smoothed next state)
    ///   x_sm_0   = x_post_0 + C_0 * dx_1
    ///   P_sm_0   = P_post_0 + C_0 * (P_sm_1 - P_pred_1) * C_0^T
    /// The predictor inverse is ridge-regularised as
    ///   (P_pred_1 * 1.0001 + 1e-12 I)^-1  (diagonal only).
    fn two_snapshot_fixture() -> EskfSmoother {
        let mut s0_pred = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
        s0_pred.cov[(0, 0)] = 5.0;
        let mut s0_post = s0_pred.clone();
        s0_post.cov[(0, 0)] = 2.0;
        let mut s1_pred = EskfState::new(Vector3::new(4.0, 0.0, 0.0), Vector3::zeros(), UnitQuaternion::identity());
        s1_pred.cov[(0, 0)] = 4.0;
        let mut s1_post = s1_pred.clone();
        s1_post.cov[(0, 0)] = 1.0;
        s1_post.pos_ecef = Vector3::new(10.0, 0.0, 0.0);
        let mut smoother = EskfSmoother::new();
        smoother.push(EskfSnapshot {
            time: GpsTime::new(2000, 0.0),
            state_pred: s0_pred,
            state_post: s0_post,
            phi: Matrix15::identity(),
            is_gnss_available: true,
        });
        smoother.push(EskfSnapshot {
            time: GpsTime::new(2000, 1.0),
            state_pred: s1_pred,
            state_post: s1_post,
            phi: Matrix15::identity(),
            is_gnss_available: false,
        });
        smoother
    }

    #[test]
    fn smoother_bookkeeping_is_ordered_and_clearable() {
        let mut s = EskfSmoother::new();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
        assert!(s.snapshots().is_empty());
        let two = two_snapshot_fixture();
        let first = two.snapshots()[0].time;
        s.push(EskfSnapshot {
            time: GpsTime::new(1999, 500.0),
            state_pred: EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity()),
            state_post: EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity()),
            phi: Matrix15::identity(),
            is_gnss_available: true,
        });
        let _ = first;
        s.push(EskfSnapshot {
            time: GpsTime::new(1999, 501.0),
            state_pred: EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity()),
            state_post: EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity()),
            phi: Matrix15::identity(),
            is_gnss_available: true,
        });
        assert_eq!(s.len(), 2);
        assert!(!s.is_empty());
        // snapshots() must preserve push order so epoch k lines up with state k.
        assert!((s.snapshots()[0].time.tow - 500.0).abs() < 1e-12);
        assert!((s.snapshots()[1].time.tow - 501.0).abs() < 1e-12);
        s.clear();
        assert_eq!(s.len(), 0);
        assert!(s.is_empty());
        assert!(s.snapshots().is_empty());
        assert_eq!(s.smooth(), Err(EngineError::EmptyHistory));
    }

    #[test]
    fn rts_backward_pass_matches_the_closed_form() {
        // dx_1 = 10 - 4 = 6 m. C_0 = 2 / (4 * 1.0001) = 0.5 / 1.0001.
        // x_sm_0 = 0 + 6 * 0.5 / 1.0001 = 3 / 1.0001 = 2.999700030
        // P_sm_0 = 2 + (0.5/1.0001)^2 * (1 - 4) = 2 - 3*(0.5/1.0001)^2 = 1.250149978
        let smoother = two_snapshot_fixture();
        let sm = smoother.smooth().expect("smooth");
        let c = 0.5 / 1.0001;
        assert!((sm[0].pos_ecef.x - 6.0 * c).abs() < 1e-9, "x = {}", sm[0].pos_ecef.x);
        assert!((sm[0].pos_ecef.y).abs() < 1e-15 && (sm[0].pos_ecef.z).abs() < 1e-15);
        assert!((sm[0].cov[(0, 0)] - (2.0 - 3.0 * c * c)).abs() < 1e-9,
            "P = {}", sm[0].cov[(0, 0)]);
        // The last epoch is the anchor: the loop never rewrites it.
        assert_eq!(sm[1].pos_ecef, Vector3::new(10.0, 0.0, 0.0));
        assert!((sm[1].cov[(0, 0)] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn smoothing_pulls_the_earlier_state_toward_the_later_discrepancy() {
        // Qualitative: a forward gap of 6 m must move the earlier estimate
        // forward by strictly less than the full gap (Kalman gain < 1) and
        // strictly more than zero.
        let smoother = two_snapshot_fixture();
        let sm = smoother.smooth().expect("smooth");
        assert!(sm[0].pos_ecef.x > 0.0);
        assert!(sm[0].pos_ecef.x < 6.0, "gain must be < 1, got {}", sm[0].pos_ecef.x);
        // Uncertainty must shrink or stay equal relative to the prior.
        assert!(sm[0].cov[(0, 0)] < 5.0);
    }

    #[test]
    fn smoother_keeps_covariance_symmetric_through_the_backward_pass() {
        let mut smoother = EskfSmoother::new();
        for k in 0..5 {
            let mut s_pred = EskfState::new(
                Vector3::new(k as f64, 0.0, 0.0),
                Vector3::zeros(),
                UnitQuaternion::identity(),
            );
            s_pred.cov[(0, 0)] = 4.0 + k as f64;
            s_pred.cov[(0, 1)] = 0.5;
            s_pred.cov[(1, 0)] = 0.5;
            let mut s_post = s_pred.clone();
            s_post.cov[(0, 0)] = 0.5 + k as f64;
            smoother.push(EskfSnapshot {
                time: GpsTime::new(2000, k as f64),
                state_pred: s_pred,
                state_post: s_post,
                phi: Matrix15::identity(),
                is_gnss_available: k % 2 == 0,
            });
        }
        let sm = smoother.smooth().expect("smooth");
        for (k, s) in sm.iter().enumerate() {
            assert!((s.cov - s.cov.transpose()).norm() < 1e-12, "epoch {k} covariance not symmetric");
            assert!(s.cov[(0, 0)] > 0.0, "epoch {k} lost positive variance");
        }
    }

    #[test]
    fn smooth_with_time_pairs_each_epoch_with_its_state() {
        let smoother = two_snapshot_fixture();
        let states = smoother.smooth().expect("smooth");
        let timed = smoother.smooth_with_time().expect("smooth_with_time");
        assert_eq!(timed.len(), states.len());
        for (i, (t, s)) in timed.iter().enumerate() {
            assert!((t.tow - smoother.snapshots()[i].time.tow).abs() < 1e-12);
            assert!((t.week - smoother.snapshots()[i].time.week) == 0);
            assert_eq!(s.pos_ecef, states[i].pos_ecef);
        }
        assert_eq!(timed[0].0, GpsTime::new(2000, 0.0));
        assert_eq!(timed[1].0, GpsTime::new(2000, 1.0));
        // An empty smoother must fail the same way through this entry point.
        assert_eq!(EskfSmoother::new().smooth_with_time(), Err(EngineError::EmptyHistory));
    }
}
