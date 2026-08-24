//! Pass 4: Optimal Bidirectional Fusion and Covariance Intersection Smoother.
//!
//! Fuses forward and backward filter runs into an optimal, continuous,
//! post-processed trajectory with rigorous error bounds.

use std::collections::BTreeMap;
use nalgebra::{Matrix3, Vector3};

use gneiss_core::time::GpsTime;
use crate::post_process::forward::FilteredEpoch;

/// A fully smoothed post-processed epoch.
#[derive(Debug, Clone)]
pub struct SmoothedEpoch {
    pub time: GpsTime,
    pub position_ecef: Vector3<f64>,
    pub velocity_ecef: Option<Vector3<f64>>,
    pub attitude: Option<nalgebra::UnitQuaternion<f64>>,
    pub cov_position: Matrix3<f64>,
    pub std_east: f64,
    pub std_north: f64,
    pub std_up: f64,
    pub separation_3d: f64,
    pub quality: u8,
    pub n_satellites: usize,
}

/// Combine forward and backward trajectories using covariance intersection.
///
/// `strict_disagreement` enables the long-baseline honesty rule: when the
/// two passes disagree by more than [`STRICT_DISAGREE_M`], neither fixed
/// claim is trusted and the product is downgraded to float quality.
pub fn combine_trajectories(
    forward: &[FilteredEpoch],
    backward: &BTreeMap<u64, FilteredEpoch>,
    strict_disagreement: bool,
) -> Vec<SmoothedEpoch> {
    let mut smoothed = Vec::with_capacity(forward.len());

    for fwd in forward {
        let tow_ms = (fwd.time.tow * 1000.0).round() as u64;
        let bwd_opt = backward.get(&tow_ms);

        let epoch_res = match bwd_opt {
            Some(bwd) => combine_bidirectional_epoch(fwd, bwd, strict_disagreement),
            None => single_pass_epoch(fwd),
        };
        smoothed.push(epoch_res);
    }
    smoothed
}

/// Passes disagreeing beyond this cannot both be right for a static
/// monument: whichever side wins, a fixed claim would be dishonest.
pub const STRICT_DISAGREE_M: f64 = 0.50;

/// Helper to fuse forward and backward estimates for a single epoch.
fn combine_bidirectional_epoch(fwd: &FilteredEpoch, bwd: &FilteredEpoch, strict: bool) -> SmoothedEpoch {
    let sep = (fwd.position_ecef - bwd.position_ecef).norm();

    // Strict mode (long-baseline path): large disagreement means at least
    // one pass is wrong. Fuse as float with honest quality instead of
    // blessing either side's fixed claim.
    if strict && sep > STRICT_DISAGREE_M {
        let q_merged = fwd.quality.min(bwd.quality);
        let (pos, cov, _) = fuse_covariances(fwd, bwd, q_merged);
        let (std_e, std_n, std_u) = compute_enu_stds(pos, cov);
        return SmoothedEpoch {
            time: fwd.time,
            position_ecef: pos,
            velocity_ecef: match (&fwd.velocity_ecef, &bwd.velocity_ecef) {
                (Some(vf), Some(vb)) => Some(0.5 * (vf + vb)),
                (Some(vf), None) => Some(*vf),
                (None, Some(vb)) => Some(*vb),
                (None, None) => None,
            },
            attitude: None,
            cov_position: cov,
            std_east: std_e,
            std_north: std_n,
            std_up: std_u,
            separation_3d: sep,
            quality: 2,
            n_satellites: fwd.n_satellites.max(bwd.n_satellites),
        };
    }

    let (pos, cov, q) = if fwd.is_fixed && !bwd.is_fixed {
        if sep < 0.50 {
            fuse_covariances(fwd, bwd, 1)
        } else {
            (fwd.position_ecef, fwd.cov_position, 1)
        }
    } else if bwd.is_fixed && !fwd.is_fixed {
        if sep < 0.50 {
            fuse_covariances(fwd, bwd, 1)
        } else {
            (bwd.position_ecef, bwd.cov_position, 1)
        }
    } else if fwd.is_fixed && bwd.is_fixed {
        if sep < 0.20 {
            fuse_covariances(fwd, bwd, 1)
        } else if fwd.cov_position.trace() <= bwd.cov_position.trace() {
            (fwd.position_ecef, fwd.cov_position, 1)
        } else {
            (bwd.position_ecef, bwd.cov_position, 1)
        }
    } else {
        let q_merged = fwd.quality.min(bwd.quality);
        if sep < 10.0 {
            fuse_covariances(fwd, bwd, q_merged)
        } else if fwd.cov_position.trace() <= bwd.cov_position.trace() {
            (fwd.position_ecef, fwd.cov_position, q_merged)
        } else {
            (bwd.position_ecef, bwd.cov_position, q_merged)
        }
    };

    let (std_e, std_n, std_u) = compute_enu_stds(pos, cov);
    let vel = match (fwd.velocity_ecef, bwd.velocity_ecef) {
        (Some(vf), Some(vb)) => Some(0.5 * vf + 0.5 * vb),
        (Some(vf), None) => Some(vf),
        (None, Some(vb)) => Some(vb),
        (None, None) => None,
    };

    SmoothedEpoch {
        time: fwd.time,
        position_ecef: pos,
        velocity_ecef: vel,
        attitude: None,
        cov_position: cov,
        std_east: std_e,
        std_north: std_n,
        std_up: std_u,
        separation_3d: sep,
        quality: q,
        n_satellites: fwd.n_satellites.max(bwd.n_satellites),
    }
}

/// Helper for epochs present in only one filter direction.
fn single_pass_epoch(fwd: &FilteredEpoch) -> SmoothedEpoch {
    let (std_e, std_n, std_u) = compute_enu_stds(fwd.position_ecef, fwd.cov_position);
    SmoothedEpoch {
        time: fwd.time,
        position_ecef: fwd.position_ecef,
        velocity_ecef: fwd.velocity_ecef,
        attitude: None,
        cov_position: fwd.cov_position,
        std_east: std_e,
        std_north: std_n,
        std_up: std_u,
        separation_3d: 0.0,
        quality: fwd.quality,
        n_satellites: fwd.n_satellites,
    }
}

/// Performs optimal inverse-covariance weighting between two estimates.
fn fuse_covariances(fwd: &FilteredEpoch, bwd: &FilteredEpoch, q: u8) -> (Vector3<f64>, Matrix3<f64>, u8) {
    let inv_fwd = fwd.cov_position.try_inverse().unwrap_or_else(|| Matrix3::identity() * 100.0);
    let inv_bwd = bwd.cov_position.try_inverse().unwrap_or_else(|| Matrix3::identity() * 100.0);
    let inv_sum = inv_fwd + inv_bwd;
    let cov = inv_sum.try_inverse().unwrap_or_else(|| fwd.cov_position * 0.5);

    let pos = cov * (inv_fwd * fwd.position_ecef + inv_bwd * bwd.position_ecef);
    (pos, cov, q)
}

/// Converts ECEF position covariance to ENU standard deviations.
fn compute_enu_stds(pos_ecef: Vector3<f64>, cov_ecef: Matrix3<f64>) -> (f64, f64, f64) {
    let llh = gneiss_core::coords::ecef_to_llh(pos_ecef);
    let r_enu_ecef = gneiss_core::coords::ecef_to_ned_matrix(llh);
    let cov_enu = r_enu_ecef * cov_ecef * r_enu_ecef.transpose();

    let std_e = cov_enu[(1, 1)].max(0.0).sqrt();
    let std_n = cov_enu[(0, 0)].max(0.0).sqrt();
    let std_u = cov_enu[(2, 2)].max(0.0).sqrt();
    (std_e, std_n, std_u)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fuse_covariances_reduces_uncertainty() {
        let p1 = Vector3::new(100.0, 0.0, 0.0);
        let p2 = Vector3::new(100.2, 0.0, 0.0);
        let cov1 = Matrix3::identity() * 0.04; // 20cm
        let cov2 = Matrix3::identity() * 0.04;

        let ep1 = FilteredEpoch {
            time: GpsTime::new(2000, 100.0), position_ecef: p1, velocity_ecef: None,
            attitude: None, cov_position: cov1, n_satellites: 8, quality: 2, is_fixed: false,
        };
        let ep2 = FilteredEpoch {
            time: GpsTime::new(2000, 100.0), position_ecef: p2, velocity_ecef: None,
            attitude: None, cov_position: cov2, n_satellites: 8, quality: 2, is_fixed: false,
        };

        let (pos, cov, _) = fuse_covariances(&ep1, &ep2, 2);
        assert!((pos.x - 100.1).abs() < 1e-6);
        assert!(cov[(0, 0)] < cov1[(0, 0)]);
    }
}
