//! Rauch-Tung-Striebel (RTS) Backward Smoother for RTK.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use gneiss_core::time::GpsTime;
use crate::post_process::combiner::SmoothedEpoch;

/// Per-epoch snapshot saved during the forward IEKF pass.
#[derive(Debug, Clone)]
pub struct IekfSnapshot {
    pub time: GpsTime,
    pub x_pred: DVector<f64>,
    pub p_pred: DMatrix<f64>,
    pub x_post: DVector<f64>,
    pub p_post: DMatrix<f64>,
    pub f_mat: DMatrix<f64>,
    pub is_fixed: bool,
    pub n_sats: usize,
    pub quality: u8,
    /// DoubleDiffKey list matching the ambiguity columns of x_post
    /// (starting at amb_offset). Empty when the state has no ambiguities;
    /// only populated when `GnssRtkIekf.track_ambiguity_keys` is set.
    pub amb_keys: Vec<crate::estimators::rtk_iekf::state::DoubleDiffKey>,
}

/// Run full RTS backward smoothing over all stored snapshots.
pub fn run_rts_smoother(snapshots: &[IekfSnapshot]) -> Vec<SmoothedEpoch> {
    let n = snapshots.len();
    if n == 0 {
        return Vec::new();
    }

    let mut smoothed: Vec<(DVector<f64>, DMatrix<f64>)> = Vec::with_capacity(n);
    smoothed.resize(n, (DVector::zeros(0), DMatrix::zeros(0, 0)));
    smoothed[n - 1] = (snapshots[n - 1].x_post.clone(), snapshots[n - 1].p_post.clone());

    for k in (0..n - 1).rev() {
        let (x_next, p_next) = &smoothed[k + 1];
        smoothed[k] = smooth_single_step(&snapshots[k], &snapshots[k + 1], x_next, p_next);
    }

    build_smoothed_epochs(snapshots, &smoothed)
}

fn smooth_single_step(
    cur: &IekfSnapshot,
    next: &IekfSnapshot,
    x_next: &DVector<f64>,
    p_next: &DMatrix<f64>,
) -> (DVector<f64>, DMatrix<f64>) {
    let common_dim = cur.x_post.len().min(6).min(next.x_pred.len().min(6));
    let p_cur = cur.p_post.view_range(0..common_dim, 0..common_dim).clone_owned();
    let f_blk = next.f_mat.view_range(0..common_dim, 0..common_dim).clone_owned();
    let p_pred = next.p_pred.view_range(0..common_dim, 0..common_dim).clone_owned();

    let p_pred_inv = p_pred.try_inverse().unwrap_or_else(|| DMatrix::identity(common_dim, common_dim) * 0.01);
    let c_gain = &p_cur * f_blk.transpose() * &p_pred_inv;

    let dx = x_next.rows(0, common_dim) - next.x_pred.rows(0, common_dim);
    let x_s = cur.x_post.rows(0, common_dim) + &c_gain * dx;

    let dp = p_next.view_range(0..common_dim, 0..common_dim) - next.p_pred.view_range(0..common_dim, 0..common_dim);
    let p_s = &p_cur + &c_gain * dp * c_gain.transpose();

    (x_s, p_s)
}

fn build_smoothed_epochs(
    snapshots: &[IekfSnapshot],
    smoothed: &[(DVector<f64>, DMatrix<f64>)],
) -> Vec<SmoothedEpoch> {
    let mut out = Vec::with_capacity(snapshots.len());
    for (snap, (x_s, p_s)) in snapshots.iter().zip(smoothed.iter()) {
        out.push(build_single_smoothed_epoch(snap, x_s, p_s));
    }
    out
}

fn build_single_smoothed_epoch(
    snap: &IekfSnapshot,
    x_s: &DVector<f64>,
    p_s: &DMatrix<f64>,
) -> SmoothedEpoch {
    let pos = Vector3::new(x_s[0], x_s[1], x_s[2]);
    let vel = if x_s.len() >= 6 {
        Some(Vector3::new(x_s[3], x_s[4], x_s[5]))
    } else {
        None
    };

    let mut cov_pos = Matrix3::zeros();
    for r in 0..3 {
        for c in 0..3 {
            cov_pos[(r, c)] = p_s[(r, c)];
        }
    }

    let (std_e, std_n, std_u) = compute_enu_stds(pos, cov_pos);

    SmoothedEpoch {
        time: snap.time,
        position_ecef: pos,
        velocity_ecef: vel,
        attitude: None,
        cov_position: cov_pos,
        std_east: std_e,
        std_north: std_n,
        std_up: std_u,
        separation_3d: 0.0,
        quality: snap.quality,
        n_satellites: snap.n_sats,
    }
}

fn compute_enu_stds(pos_ecef: Vector3<f64>, cov_ecef: Matrix3<f64>) -> (f64, f64, f64) {
    let llh = gneiss_core::coords::ecef_to_llh(pos_ecef);
    let r_enu = gneiss_core::coords::ecef_to_ned_matrix(llh);
    let cov_enu = r_enu * cov_ecef * r_enu.transpose();
    (
        cov_enu[(1, 1)].max(0.0).sqrt(),
        cov_enu[(0, 0)].max(0.0).sqrt(),
        cov_enu[(2, 2)].max(0.0).sqrt(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rts_smoother_empty() {
        let res = run_rts_smoother(&[]);
        assert!(res.is_empty());
    }

    #[test]
    fn test_rts_smoother_single_epoch() {
        let snap = IekfSnapshot {
            time: GpsTime::new(2000, 100.0),
            x_pred: DVector::zeros(6),
            p_pred: DMatrix::identity(6, 6),
            x_post: DVector::zeros(6),
            p_post: DMatrix::identity(6, 6) * 0.5,
            f_mat: DMatrix::identity(6, 6),
            is_fixed: true,
            n_sats: 8,
            quality: 1,
            amb_keys: Vec::new(),
        };
        let res = run_rts_smoother(&[snap]);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].quality, 1);
    }
}
