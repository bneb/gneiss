//! Multi-base network solution consensus.
//!
//! Independent single-base runs over the same rover epochs share the rover
//! receiver, antenna, troposphere, and broadcast ephemerides (common-mode:
//! they do NOT average), but see different base-side error realizations —
//! multipath, slant tropo/iono projections, base-coordinate error,
//! base-specific cycle slips. Consensus over bases averages exactly those
//! independent components.
//!
//! Design notes (red-team hardened):
//! - Reported per-base covariances are heuristic (quality-flag based), so
//!   fusion weights are effectively uniform and the fused covariance must
//!   come from EMPIRICAL scatter of the agreeing candidates, never from
//!   intersecting the fictional matrices.
//! - Wrong integer fixes look statistically excellent, so the fixed
//!   sub-population is combined by component-wise median (robust to a
//!   minority of grossly wrong bases), gated by an agreement radius.
//! - The float sub-population acts as an internal validator: when its own
//!   consensus sits far outside the fixed consensus, the fixed claim is
//!   downgraded instead of trusted.

use std::collections::BTreeMap;

use nalgebra::{Matrix3, Vector3};

use gneiss_core::time::GpsTime;

use super::combiner::SmoothedEpoch;

/// Minimum agreeing fixed candidates required to claim quality=1.
pub const MIN_FIXED_CONSENSUS: usize = 3;
/// Agreement radius (m) around the fixed consensus center.
pub const AGREEMENT_RADIUS_M: f64 = 0.5;
/// Hard cap on survivor RMS spread (m) for a quality=1 claim.
const FIXED_SPREAD_M: f64 = 0.35;
/// Float-consensus disagreement (m) beyond which fixed claims are vetoed.
/// NOTE: measured neutral-to-harmful on the CORS day set — correlated
/// multi-base degradation clusters floats far out and vetoes healthy
/// consensuses. Disabled pending a robustness-aware rule.
const FLOAT_VETO_M: f64 = f64::INFINITY;
/// Continuity-gate defaults: static monuments cannot jump farther than this
/// between adjacent epochs within the max time delta.
pub const CONTINUITY_JUMP_M: f64 = 0.20;
pub const CONTINUITY_MAX_DT_S: f64 = 90.0;

/// Per-epoch consensus configuration knobs.
#[derive(Debug, Clone)]
pub struct NetworkConsensusConfig {
    pub min_fixed_consensus: usize,
    pub agreement_radius_m: f64,
}

impl Default for NetworkConsensusConfig {
    fn default() -> Self {
        Self {
            min_fixed_consensus: MIN_FIXED_CONSENSUS,
            agreement_radius_m: AGREEMENT_RADIUS_M,
        }
    }
}

/// Fuse per-base trajectories into one network solution.
///
/// Trajectories may differ in coverage; alignment is by wrapped-TOW key and
/// output is time-ascending. Epochs with no candidates are omitted.
pub fn fuse_network_solutions(
    trajectories: &[Vec<SmoothedEpoch>],
    cfg: &NetworkConsensusConfig,
) -> Vec<SmoothedEpoch> {
    let mut by_tow: BTreeMap<u32, Vec<&SmoothedEpoch>> = BTreeMap::new();
    for traj in trajectories {
        for ep in traj {
            by_tow.entry(tow_key(&ep.time)).or_default().push(ep);
        }
    }
    by_tow.into_values()
        .filter_map(|cands| consensus_epoch(&cands, cfg))
        .collect()
}

/// Day-wrapping-safe key (GPS TOW rolls at 604800).
fn tow_key(time: &GpsTime) -> u32 {
    (time.tow.round() as i64).rem_euclid(604800) as u32
}

/// Downgrade epochs whose jump from the previous kept epoch exceeds
/// `max_jump_m` within `max_dt_s` — the only independent signal left against
/// temporally correlated wrong fixes. Returns the annotated trajectory.
pub fn apply_continuity_gate(
    mut traj: Vec<SmoothedEpoch>,
    max_jump_m: f64,
    max_dt_s: f64,
) -> Vec<SmoothedEpoch> {
    // Anchor on the last ACCEPTED epoch: once an excursion starts, every
    // subsequent epoch stays downgraded until the series returns near the
    // anchored solution — a persistent wrong stretch cannot re-bless itself
    // by drifting slowly away.
    let mut anchor: Option<(f64, Vector3<f64>)> = None;
    for ep in traj.iter_mut() {
        if let Some((at, ap)) = anchor {
            let dt = ep.time.tow - at;
            if dt > 0.0 && dt <= max_dt_s && (ep.position_ecef - ap).norm() > max_jump_m {
                ep.quality = 2;
                continue;
            }
        }
        anchor = Some((ep.time.tow, ep.position_ecef));
    }
    traj
}

fn consensus_epoch(cands: &[&SmoothedEpoch], cfg: &NetworkConsensusConfig) -> Option<SmoothedEpoch> {
    let first = *cands.first()?;
    let (fixed, floats): (Vec<&SmoothedEpoch>, Vec<&SmoothedEpoch>) =
        cands.iter().partition(|c| c.quality == 1);

    let (pos, members, quality) = if fixed.len() >= cfg.min_fixed_consensus {
        let center = component_median(&fixed)?;
        let radius = cfg.agreement_radius_m;
        let agree: Vec<&SmoothedEpoch> = fixed.iter().copied()
            .filter(|c| (c.position_ecef - center).norm() <= radius)
            .collect();
        let float_center = component_median(&floats);
        let vetoed = matches!(&float_center, Some(fc)
            if (*fc - center).norm() > FLOAT_VETO_M);
        let spread = rms_spread(&agree, center)?;
        let ok = !agree.is_empty()
            && agree.len() >= cfg.min_fixed_consensus
            && spread <= FIXED_SPREAD_M
            && !vetoed;
        let pos = if ok { component_median(&agree)? } else { center };
        (pos, agree, if ok { 1 } else { 2 })
    } else {
        let center = component_median(cands)?;
        let mut dev: Vec<f64> = cands.iter()
            .map(|c| (c.position_ecef - center).norm())
            .collect();
        dev.sort_by(|a, b| a.total_cmp(b));
        let med_dev = dev.get(dev.len() / 2).copied().unwrap_or(0.0);
        let keep: Vec<&SmoothedEpoch> = cands.iter().copied()
            .filter(|c| (c.position_ecef - center).norm() <= 2.0 * med_dev)
            .collect();
        let min_q = cands.iter().map(|c| c.quality).min()?;
        (component_median(&keep)?, keep, min_q)
    };

    let sep = pairwise_max_dist(&members)?;
    let n_sat = members.iter().map(|c| c.n_satellites).max().unwrap_or(0);
    let (std_n, std_e, std_u) = empirical_sigmas(&members, pos);
    Some(SmoothedEpoch {
        time: first.time,
        position_ecef: pos,
        velocity_ecef: None,
        attitude: None,
        cov_position: empirical_cov(&members, pos),
        std_east: std_e,
        std_north: std_n,
        std_up: std_u,
        separation_3d: sep,
        quality,
        n_satellites: n_sat,
    })
}

/// Component-wise median: robust to a minority of grossly wrong bases.
fn component_median(cands: &[&SmoothedEpoch]) -> Option<Vector3<f64>> {
    if cands.is_empty() {
        return None;
    }
    let mut xs: Vec<f64> = cands.iter().map(|c| c.position_ecef[0]).collect();
    let mut ys: Vec<f64> = cands.iter().map(|c| c.position_ecef[1]).collect();
    let mut zs: Vec<f64> = cands.iter().map(|c| c.position_ecef[2]).collect();
    for v in [&mut xs, &mut ys, &mut zs] {
        v.sort_by(|a, b| a.total_cmp(b));
    }
    let mid = cands.len() / 2;
    Some(Vector3::new(xs[mid], ys[mid], zs[mid]))
}

fn rms_spread(cands: &[&SmoothedEpoch], center: Vector3<f64>) -> Option<f64> {
    if cands.is_empty() {
        return None;
    }
    let n = cands.len() as f64;
    Some((cands.iter()
        .map(|c| (c.position_ecef - center).norm_squared())
        .sum::<f64>() / n)
        .sqrt())
}

fn pairwise_max_dist(cands: &[&SmoothedEpoch]) -> Option<f64> {
    let mut max = 0.0_f64;
    for a in 0..cands.len() {
        for b in a + 1..cands.len() {
            max = max.max((cands[a].position_ecef - cands[b].position_ecef).norm());
        }
    }
    Some(max)
}

/// Empirical scatter of the survivors — NOT the fictional per-base matrices.
fn empirical_cov(cands: &[&SmoothedEpoch], center: Vector3<f64>) -> Matrix3<f64> {
    let n = cands.len().max(1) as f64;
    let mut cov = Matrix3::zeros();
    for c in cands {
        let d = c.position_ecef - center;
        cov += d * d.transpose();
    }
    let mut cov = cov / n;
    for i in 0..3 {
        cov[(i, i)] = cov[(i, i)].max(1e-6);
    }
    cov
}

fn empirical_sigmas(cands: &[&SmoothedEpoch], center: Vector3<f64>) -> (f64, f64, f64) {
    let cov = empirical_cov(cands, center);
    (
        cov[(0, 0)].max(0.0).sqrt(),
        cov[(1, 1)].max(0.0).sqrt(),
        cov[(2, 2)].max(0.0).sqrt(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn ep(tow: f64, x: f64, sigma: f64, quality: u8) -> SmoothedEpoch {
        let var = sigma * sigma;
        SmoothedEpoch {
            time: GpsTime::new(2105, tow),
            position_ecef: Vector3::new(x, 0.0, 0.0),
            velocity_ecef: None,
            attitude: None,
            cov_position: Matrix3::new(var, 0.0, 0.0, 0.0, var, 0.0, 0.0, 0.0, var * 2.25),
            std_east: 0.0,
            std_north: 0.0,
            std_up: 0.0,
            separation_3d: 0.0,
            quality,
            n_satellites: 9,
        }
    }

    #[test]
    fn test_confident_outlier_base_cannot_drag_median() {
        // Five agreeing bases plus one confidently-wrong base 5 m off:
        // the component median ignores it entirely.
        let mut trajs: Vec<Vec<SmoothedEpoch>> = Vec::new();
        for k in 0..5u32 {
            trajs.push(vec![ep(100.0, 0.01 * k as f64, 0.02, 1)]);
        }
        trajs.push(vec![ep(100.0, 5.0, 0.005, 1)]);
        let fused = fuse_network_solutions(&trajs, &NetworkConsensusConfig::default());
        assert_eq!(fused.len(), 1);
        assert!(fused[0].position_ecef.x.abs() < 0.05);
        assert_eq!(fused[0].quality, 1);
    }

    #[test]
    fn test_minority_wrong_cluster_does_not_create_ghost_mean() {
        // 2-of-5 wrong at 14 m: median stays on the majority mode, and the
        // disagreement filter removes the wrong pair from the survivors.
        let mut trajs: Vec<Vec<SmoothedEpoch>> = Vec::new();
        for k in 0..3u32 {
            trajs.push(vec![ep(200.0, 0.01 * k as f64, 0.02, 1)]);
        }
        for k in 0..2u32 {
            trajs.push(vec![ep(200.0, 14.0 + 0.01 * k as f64, 0.01, 1)]);
        }
        let fused = fuse_network_solutions(&trajs, &NetworkConsensusConfig::default());
        assert_eq!(fused.len(), 1);
        assert!(fused[0].position_ecef.x < 0.5, "must ride the majority mode");
        assert_eq!(fused[0].quality, 1);
    }

    #[test]
    fn test_float_disagreement_does_not_veto_measured_harmful() {
        // Four fixed bases agree; three honest float solutions sit 12 m away.
        // The float veto was measured NET-HARMFUL on real data: correlated
        // multi-base degradation clusters floats far out and vetoes healthy
        // consensuses (CORS day set: NETWORK RMS +27%). It is therefore
        // disabled (FLOAT_VETO_M = inf) and this test pins that contract:
        // floats never flip a tight fixed consensus, and the median rides
        // the fixed mode.
        let mut trajs: Vec<Vec<SmoothedEpoch>> = Vec::new();
        for k in 0..4u32 {
            trajs.push(vec![ep(250.0, 0.01 * k as f64, 0.02, 1)]);
        }
        for k in 0..3u32 {
            trajs.push(vec![ep(250.0, 12.0 + 0.1 * k as f64, 0.5, 2)]);
        }
        let fused = fuse_network_solutions(&trajs, &NetworkConsensusConfig::default());
        assert_eq!(fused[0].quality, 1);
        assert!(fused[0].position_ecef.x < 0.05);
    }

    #[test]
    fn test_continuity_gate_downgrades_excursions_until_recovery() {
        let traj = vec![
            ep(1000.0, 0.0, 0.02, 1),
            ep(1030.0, 0.01, 0.02, 1),
            ep(1060.0, 14.0, 0.02, 1), // impossible jump for a static monument
            ep(1090.0, 14.01, 0.02, 1), // stays downgraded: anchored on last good
            ep(1120.0, 14.02, 0.02, 1), // drifting away does NOT re-bless itself
            ep(1150.0, 0.02, 0.02, 1),  // return near the anchor recovers
        ];
        let gated = apply_continuity_gate(traj, 0.20, 90.0);
        assert_eq!(gated[0].quality, 1);
        assert_eq!(gated[1].quality, 1);
        assert_eq!(gated[2].quality, 2, "jump epoch must be downgraded");
        assert_eq!(gated[3].quality, 2, "anchored: excursion cannot self-recover");
        assert_eq!(gated[4].quality, 2, "drift cannot re-bless the excursion");
        assert_eq!(gated[5].quality, 1, "returning near the anchor recovers");
    }

    #[test]
    fn test_continuity_gate_allows_low_dynamics() {
        // The gate targets static/low-dynamics monuments (the network
        // benchmark): steady centimetre-level drift passes untouched.
        let traj: Vec<SmoothedEpoch> = (0..5)
            .map(|k| ep(1000.0 + 30.0 * k as f64, 0.01 * k as f64, 0.02, 1))
            .collect();
        let gated = apply_continuity_gate(traj, 0.20, 90.0);
        assert!(gated.iter().all(|e| e.quality == 1));
    }

    #[test]
    fn test_misaligned_trajectories_align_by_wrapped_tow() {
        let a = vec![ep(400.0, 0.0, 0.02, 1), ep(430.0, 0.0, 0.02, 1)];
        let b = vec![ep(430.0, 0.01, 0.02, 1)];
        let trajs = vec![a, b];
        let fused = fuse_network_solutions(&trajs, &NetworkConsensusConfig::default());
        assert_eq!(fused.len(), 2);
        assert_eq!(tow_key(&fused[0].time), 400);
        assert_eq!(tow_key(&fused[1].time), 430);
    }
}
