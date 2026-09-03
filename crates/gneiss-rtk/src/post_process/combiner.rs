//! Pass 4: Optimal Bidirectional Fusion and Covariance Intersection Smoother.
//!
//! Fuses forward and backward filter runs into an optimal, continuous,
//! post-processed trajectory with rigorous error bounds.

use std::collections::BTreeMap;
use nalgebra::{Matrix3, Vector3};

use gneiss_core::time::GpsTime;
use crate::post_process::dynamics;
use crate::post_process::dynamics::ProcessingDynamics;
use crate::post_process::iekf_pass::FilteredEpoch;

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
/// two passes disagree by more than the profile's disagreement limit,
/// neither fixed claim is trusted and the product is downgraded to float
/// quality. The limit is [`STRICT_DISAGREE_M`] for static monuments and
/// sigma-scaled (see [`dynamics::kinematic_sep_limit_m`]) for kinematic
/// rovers, where absolute metres would flag legitimate motion as fraud.
pub fn combine_trajectories(
    forward: &[FilteredEpoch],
    backward: &BTreeMap<u64, FilteredEpoch>,
    strict_disagreement: bool,
    prof: ProcessingDynamics,
) -> Vec<SmoothedEpoch> {
    let mut smoothed = Vec::with_capacity(forward.len());

    for fwd in forward {
        let tow_ms = (fwd.time.tow * 1000.0).round() as u64;
        let bwd_opt = backward.get(&tow_ms);

        let epoch_res = match bwd_opt {
            Some(bwd) => combine_bidirectional_epoch(fwd, bwd, strict_disagreement, prof),
            None => single_pass_epoch(fwd),
        };
        smoothed.push(epoch_res);
    }
    smoothed
}

/// Passes disagreeing beyond this cannot both be right for a static
/// monument: whichever side wins, a fixed claim would be dishonest.
pub const STRICT_DISAGREE_M: f64 = 0.50;

/// Formal (1-sigma) position magnitude of one pass, from the trace of
/// its ECEF covariance: sqrt(trace/3) of an isotropic-equivalent error.
fn formal_sigma_m(cov: &Matrix3<f64>) -> f64 {
    (cov.trace().max(0.0) / 3.0).sqrt()
}

/// Profile-dependent separation limits for one epoch pair.
struct SepLimits {
    /// Honesty gate: beyond this, fixed claims are downgraded.
    strict_m: f64,
    /// Fuse window when both passes claim fixes.
    both_fixed_fuse_m: f64,
}

impl SepLimits {
    fn for_profile(prof: ProcessingDynamics, fwd: &FilteredEpoch, bwd: &FilteredEpoch) -> Self {
        match prof {
            ProcessingDynamics::Static => SepLimits {
                strict_m: STRICT_DISAGREE_M,
                both_fixed_fuse_m: 0.20,
            },
            ProcessingDynamics::Kinematic => {
                let (sf, sb) = (formal_sigma_m(&fwd.cov_position), formal_sigma_m(&bwd.cov_position));
                // Floors equal the audited static constants, so the
                // sigma-scaled rule only ever WIDENS the validated
                // tolerances (measured: formal sigmas are >10x
                // optimistic vs realized fwd/bwd disagreement).
                SepLimits {
                    strict_m: dynamics::kinematic_sep_limit_m(
                        dynamics::KIN_DISAGREE_K_SIGMA, sf, sb,
                        dynamics::KIN_THRESHOLD_FLOOR_M, dynamics::KIN_THRESHOLD_CAP_M,
                    ),
                    both_fixed_fuse_m: dynamics::kinematic_sep_limit_m(
                        dynamics::KIN_FUSE_BOTH_FIXED_K_SIGMA, sf, sb,
                        dynamics::KIN_BOTH_FIXED_FLOOR_M, dynamics::KIN_BOTH_FIXED_CAP_M,
                    ),
                }
            }
        }
    }
}

/// Helper to fuse forward and backward estimates for a single epoch.
fn combine_bidirectional_epoch(
    fwd: &FilteredEpoch,
    bwd: &FilteredEpoch,
    strict: bool,
    prof: ProcessingDynamics,
) -> SmoothedEpoch {
    let sep = (fwd.position_ecef - bwd.position_ecef).norm();
    let limits = SepLimits::for_profile(prof, fwd, bwd);

    let (pos, cov, mut q) = if fwd.is_fixed && !bwd.is_fixed {
        (fwd.position_ecef, fwd.cov_position, 1)
    } else if bwd.is_fixed && !fwd.is_fixed {
        (bwd.position_ecef, bwd.cov_position, 1)
    } else if fwd.is_fixed && bwd.is_fixed {
        if sep < limits.both_fixed_fuse_m {
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

    // Long-baseline honesty: when both passes claim fixed integers but disagree beyond
    // the profile threshold, neither can be trusted as fixed.
    // When only one pass is fixed, allow the fixed solution unless gross divergence (>2.0m).
    let both_fixed = fwd.is_fixed && bwd.is_fixed;
    let limit_m = if both_fixed { limits.strict_m } else { limits.strict_m.max(2.0) };
    if strict && sep > limit_m && q == 1 {
        q = 2;
    }

    let (std_e, std_n, std_u) = gneiss_core::coords::ecef_cov_to_enu_std(pos, cov);
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
    let (std_e, std_n, std_u) = gneiss_core::coords::ecef_cov_to_enu_std(fwd.position_ecef, fwd.cov_position);
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

    /// Synthetic fixed-claim epoch pair at `sep_m` separation with both
    /// passes reporting isotropic formal sigma `sigma_m`.
    fn fixed_pair(sep_m: f64, sigma_m: f64) -> (FilteredEpoch, FilteredEpoch) {
        let mk = |x: f64| FilteredEpoch {
            time: GpsTime::new(2000, 100.0),
            position_ecef: Vector3::new(x, 0.0, 0.0),
            velocity_ecef: None,
            attitude: None,
            cov_position: Matrix3::identity() * (sigma_m * sigma_m),
            n_satellites: 8,
            quality: 1,
            is_fixed: true,
        };
        (mk(0.0), mk(sep_m))
    }

    fn combined_quality(sep_m: f64, sigma_m: f64, prof: ProcessingDynamics) -> u8 {
        let (fwd, bwd) = fixed_pair(sep_m, sigma_m);
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, bwd);
        combine_trajectories(&[fwd], &map, true, prof).remove(0).quality
    }

    #[test]
    fn static_profile_keeps_legacy_absolute_rule() {
        // Legacy behaviour: 0.60 m disagreement downgrades a static fix.
        assert_eq!(combined_quality(0.60, 0.01, ProcessingDynamics::Static), 2);
        // Just inside stays fixed.
        assert_eq!(combined_quality(0.40, 0.01, ProcessingDynamics::Static), 1);
    }

    #[test]
    fn kinematic_threshold_scales_with_reported_sigma() {
        // Same 0.60 m separation as the failing static case, but the
        // passes report formal sigma 0.20 m: limit = 6*0.20 = 1.2 m,
        // so honest agreement tracking motion stays fixed.
        assert_eq!(combined_quality(0.60, 0.20, ProcessingDynamics::Kinematic), 1);
        // Scale-up property: larger reported sigma, wider tolerance.
        assert_eq!(combined_quality(1.00, 0.30, ProcessingDynamics::Kinematic), 1);
        // Tiny sigmas bind at the audited 0.50 m floor: the kinematic
        // rule is never stricter than the validated static bound.
        assert_eq!(combined_quality(0.45, 1e-4, ProcessingDynamics::Kinematic), 1);
        // Just beyond the floor it still bites even with tiny sigmas.
        assert_eq!(combined_quality(0.55, 1e-4, ProcessingDynamics::Kinematic), 2);
        // Huge sigmas hit the cap (10 m): divergence cannot excuse fraud.
        assert_eq!(combined_quality(11.0, 1e3, ProcessingDynamics::Kinematic), 2);
    }

    #[test]
    fn kinematic_downgrades_where_static_would_and_more() {
        // For identical epochs, kinematic must never be STRICTER than
        // static beyond the floor: its limit is >= floor and grows with
        // sigma, so any sep the static rule forgives is forgiven too
        // whenever sigma >= ~8.3 cm (floor/6).
        for sigma in [0.10_f64, 0.15, 0.30, 1.0] {
            for sep in [0.05_f64, 0.45, 0.49] {
                let qs = combined_quality(sep, sigma, ProcessingDynamics::Static);
                let qk = combined_quality(sep, sigma, ProcessingDynamics::Kinematic);
                assert!(
                    qk <= qs,
                    "kinematic (sigma={sigma}) must not downgrade sep={sep} that static forgives"
                );
            }
        }
    }

    #[test]
    fn strict_off_never_downgrades() {
        let (fwd, bwd) = fixed_pair(50.0, 0.01);
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, bwd);
        let mut out = combine_trajectories(&[fwd], &map, false, ProcessingDynamics::Static);
        assert_eq!(out.remove(0).quality, 1);
    }
}
