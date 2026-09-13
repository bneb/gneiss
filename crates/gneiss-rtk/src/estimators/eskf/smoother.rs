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
}
