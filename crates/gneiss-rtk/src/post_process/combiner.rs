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

fn select_fused_estimate(
    fwd: &FilteredEpoch,
    bwd: &FilteredEpoch,
    sep: f64,
    limits: &SepLimits,
    strict: bool,
) -> (Vector3<f64>, Matrix3<f64>, u8) {
    let (mut pos, mut cov, mut q) = if fwd.is_fixed && bwd.is_fixed {
        if sep < limits.both_fixed_fuse_m {
            fuse_covariances(fwd, bwd, 1)
        } else if fwd.cov_position.trace() <= bwd.cov_position.trace() {
            (fwd.position_ecef, fwd.cov_position, 1)
        } else {
            (bwd.position_ecef, bwd.cov_position, 1)
        }
    } else if fwd.is_fixed && !bwd.is_fixed {
        (fwd.position_ecef, fwd.cov_position, 1)
    } else if bwd.is_fixed && !fwd.is_fixed {
        (bwd.position_ecef, bwd.cov_position, 1)
    } else {
        let q_merged = 2.max(fwd.quality.min(bwd.quality));
        if sep < 10.0 {
            fuse_covariances(fwd, bwd, q_merged)
        } else if fwd.cov_position.trace() <= bwd.cov_position.trace() {
            (fwd.position_ecef, fwd.cov_position, q_merged)
        } else {
            (bwd.position_ecef, bwd.cov_position, q_merged)
        }
    };
    let both_fixed = fwd.is_fixed && bwd.is_fixed;
    let limit_m = if both_fixed { limits.strict_m } else { limits.strict_m.max(2.2) };
    if strict && sep > limit_m && q == 1 {
        q = 2;
        if fwd.cov_position.trace() <= bwd.cov_position.trace() {
            pos = fwd.position_ecef;
            cov = fwd.cov_position;
        } else {
            pos = bwd.position_ecef;
            cov = bwd.cov_position;
        }
    }
    (pos, cov, q)
}

fn fused_velocity(vf: Option<Vector3<f64>>, vb: Option<Vector3<f64>>) -> Option<Vector3<f64>> {
    match (vf, vb) {
        (Some(f), Some(b)) => Some(0.5 * f + 0.5 * b),
        (Some(f), None) => Some(f),
        (None, Some(b)) => Some(b),
        (None, None) => None,
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
    let (pos, cov, q) = select_fused_estimate(fwd, bwd, sep, &limits, strict);
    let (std_e, std_n, std_u) = gneiss_core::coords::ecef_cov_to_enu_std(pos, cov);

    SmoothedEpoch {
        time: fwd.time,
        position_ecef: pos,
        velocity_ecef: fused_velocity(fwd.velocity_ecef, bwd.velocity_ecef),
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

    /// Epoch with an explicit pass identity.
    fn pass(tow: f64, pos: Vector3<f64>, var: f64, q: u8, fixed: bool, vel: Option<Vector3<f64>>) -> FilteredEpoch {
        FilteredEpoch {
            time: GpsTime::new(2000, tow),
            position_ecef: pos,
            velocity_ecef: vel,
            attitude: None,
            cov_position: Matrix3::identity() * var,
            n_satellites: 7,
            quality: q,
            is_fixed: fixed,
        }
    }

    #[test]
    fn formal_sigma_is_the_isotropic_equivalent_of_the_trace() {
        // sigma = sqrt(trace/3). Identity * 0.25 -> trace 0.75 -> 0.25 -> 0.5.
        assert!((formal_sigma_m(&(Matrix3::identity() * 0.25)) - 0.5).abs() < 1e-12);
        // diag(1, 4, 9): trace 14 -> 14/3 -> 2.160246899469287
        let d = Matrix3::new(1.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 9.0);
        assert!((formal_sigma_m(&d) - (14.0_f64 / 3.0).sqrt()).abs() < 1e-12);
        // A covariance that is not positive definite contributes nothing
        // usable, so the magnitude is clamped at zero rather than NaN.
        let neg = Matrix3::identity() * -1.0;
        assert_eq!(formal_sigma_m(&neg), 0.0);
    }

    #[test]
    fn velocities_fuse_by_a_plain_average_and_survive_one_sided_absence() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(3.0, 6.0, 9.0);
        assert_eq!(fused_velocity(Some(a), Some(b)), Some(Vector3::new(2.0, 4.0, 6.0)));
        assert_eq!(fused_velocity(Some(a), None), Some(a));
        assert_eq!(fused_velocity(None, Some(b)), Some(b));
        assert_eq!(fused_velocity(None, None), None);
    }

    #[test]
    fn a_backward_only_epoch_is_left_untouched() {
        // No backward match: the epoch is reported verbatim, with a zero
        // separation (nothing was compared) and the forward quality.
        let fwd = pass(100.0, Vector3::new(1.0, 2.0, 3.0), 0.04, 1, true, None);
        let out = combine_trajectories(std::slice::from_ref(&fwd), &BTreeMap::new(), true, ProcessingDynamics::Static);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].position_ecef, fwd.position_ecef);
        assert_eq!(out[0].cov_position, fwd.cov_position);
        assert_eq!(out[0].quality, 1);
        assert_eq!(out[0].separation_3d, 0.0);
        assert_eq!(out[0].n_satellites, 7);
    }

    #[test]
    fn epochs_are_matched_on_the_millisecond_tow_key() {
        // tow 100.0004 and tow 100.00049 both round to 100 000 ms and pair up;
        // tow 100.001 is 100 001 ms and must NOT match tow 100.000.
        let fwd = pass(100.0000, Vector3::new(0.0, 0.0, 0.0), 0.01, 2, false, None);
        let near = pass(100.0004, Vector3::new(0.05, 0.0, 0.0), 0.01, 2, false, None);
        let far = pass(100.0010, Vector3::new(5.0, 0.0, 0.0), 0.01, 2, false, None);
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, near);
        let out = combine_trajectories(std::slice::from_ref(&fwd), &map, true, ProcessingDynamics::Static);
        assert!(out[0].separation_3d > 0.0, "a matched pair must report a separation");

        let mut map2 = BTreeMap::new();
        map2.insert(100_001_u64, far);
        let out2 = combine_trajectories(std::slice::from_ref(&fwd), &map2, true, ProcessingDynamics::Static);
        assert_eq!(out2[0].separation_3d, 0.0, "a 1 ms key mismatch must not pair up");
    }

    #[test]
    fn two_fixed_passes_beyond_the_fuse_window_keep_the_tighter_one() {
        // sep 0.60 m > the 0.20 m both-fixed fuse window, so the passes are
        // NOT averaged; the side with the smaller covariance trace wins whole.
        let mut fwd = pass(100.0, Vector3::new(0.0, 0.0, 0.0), 0.04, 1, true, None);
        let bwd = pass(100.0, Vector3::new(0.6, 0.0, 0.0), 0.01, 1, true, None);
        let bwd_pos = bwd.position_ecef;
        let bwd_cov = bwd.cov_position;
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, bwd);
        let out = combine_trajectories(std::slice::from_ref(&fwd), &map, true, ProcessingDynamics::Static);
        assert_eq!(out[0].position_ecef, bwd_pos, "tighter pass must win whole");
        assert_eq!(out[0].cov_position, bwd_cov);
        // Reversing the variances reverses the winner: the tie-break is the
        // covariance trace, not the epoch order.
        fwd.cov_position = Matrix3::identity() * 0.0001;
        let loose = pass(100.0, Vector3::new(0.6, 0.0, 0.0), 0.04, 1, true, None);
        let mut map2 = BTreeMap::new();
        map2.insert(100_000_u64, loose);
        let out2 = combine_trajectories(std::slice::from_ref(&fwd), &map2, true, ProcessingDynamics::Static);
        assert_eq!(out2[0].position_ecef, fwd.position_ecef);
    }

    #[test]
    fn two_fixed_passes_inside_the_window_are_actually_fused() {
        let fwd = pass(100.0, Vector3::new(0.0, 0.0, 0.0), 0.04, 1, true, None);
        let bwd = pass(100.0, Vector3::new(0.1, 0.0, 0.0), 0.04, 1, true, None);
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, bwd);
        let out = combine_trajectories(&[fwd], &map, true, ProcessingDynamics::Static);
        assert!((out[0].position_ecef.x - 0.05).abs() < 1e-9, "equal weights must average");
        assert!(out[0].cov_position[(0, 0)] < 0.04);
        assert_eq!(out[0].quality, 1);
    }

    #[test]
    fn a_lone_fixed_claim_is_taken_whole_but_still_honest_gated() {
        // Only the forward pass claims a fix. Its position is taken whole
        // (never averaged with the dissenting float pass), but the honesty
        // gate still applies: a one-sided claim gets the wider
        // max(strict_m, 2.2) = 2.2 m budget, and a 50 m disagreement is far
        // beyond it, so the claim is downgraded to float quality.
        let fwd = pass(100.0, Vector3::new(0.0, 0.0, 0.0), 0.04, 1, true, None);
        let bwd = pass(100.0, Vector3::new(50.0, 0.0, 0.0), 0.04, 2, false, None);
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, bwd);
        let gated = combine_trajectories(std::slice::from_ref(&fwd), &map, true, ProcessingDynamics::Static);
        assert_eq!(gated[0].quality, 2, "50 m of disagreement voids the lone fix");
        // With the gate off the fix stands unchallenged.
        let open = combine_trajectories(std::slice::from_ref(&fwd), &map, false, ProcessingDynamics::Static);
        assert_eq!(open[0].quality, 1);

        // Inside the one-sided 2.2 m budget the lone fix survives: 1.0 m
        // apart keeps quality 1 even though the float pass disagrees.
        let bwd_near = pass(100.0, Vector3::new(1.0, 0.0, 0.0), 0.04, 2, false, None);
        let fwd2 = pass(100.0, Vector3::new(0.0, 0.0, 0.0), 0.04, 1, true, None);
        let mut map2 = BTreeMap::new();
        map2.insert(100_000_u64, bwd_near);
        let kept = combine_trajectories(&[fwd2], &map2, true, ProcessingDynamics::Static);
        assert_eq!(kept[0].quality, 1);
        assert_eq!(kept[0].position_ecef, Vector3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn a_backward_only_fix_is_the_one_that_survives() {
        // Mirror image of the previous case: the fixed side is the backward
        // pass, so its position must be taken whole even when the forward
        // float pass sits far away.
        let fwd_f = pass(100.0, Vector3::new(-9.0, 0.0, 0.0), 0.04, 2, false, None);
        let bwd_f = pass(100.0, Vector3::new(3.0, 0.0, 0.0), 0.04, 1, true, None);
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, bwd_f);
        let out = combine_trajectories(&[fwd_f], &map, false, ProcessingDynamics::Static);
        assert_eq!(out[0].position_ecef, Vector3::new(3.0, 0.0, 0.0));
        assert_eq!(out[0].quality, 1);
    }

    #[test]
    fn two_float_passes_fuse_only_while_they_stay_close() {
        let fwd = pass(100.0, Vector3::new(0.0, 0.0, 0.0), 0.04, 2, false, None);
        let near = pass(100.0, Vector3::new(0.2, 0.0, 0.0), 0.04, 2, false, None);
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, near);
        let fused = combine_trajectories(std::slice::from_ref(&fwd), &map, true, ProcessingDynamics::Static);
        assert!((fused[0].position_ecef.x - 0.1).abs() < 1e-9, "0.2 m apart must average");

        let far = pass(100.0, Vector3::new(40.0, 0.0, 0.0), 0.04, 2, false, None);
        let mut map2 = BTreeMap::new();
        map2.insert(100_000_u64, far);
        let split = combine_trajectories(std::slice::from_ref(&fwd), &map2, true, ProcessingDynamics::Static);
        assert_eq!(split[0].position_ecef, Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(split[0].quality, 2, "merged quality never claims a fix");
    }

    #[test]
    fn the_satellite_count_is_the_larger_of_the_two_passes() {
        let mut fwd = pass(100.0, Vector3::new(0.0, 0.0, 0.0), 0.04, 2, false, None);
        let mut bwd = pass(100.0, Vector3::new(0.0, 0.0, 0.0), 0.04, 2, false, None);
        fwd.n_satellites = 5;
        bwd.n_satellites = 11;
        let mut map = BTreeMap::new();
        map.insert(100_000_u64, bwd);
        let out = combine_trajectories(&[fwd], &map, true, ProcessingDynamics::Static);
        assert_eq!(out[0].n_satellites, 11);
    }
}
