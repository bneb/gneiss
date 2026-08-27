//! Factor graph builders for undifferenced (SPP/PPP) and double-differenced (RTK) epochs.

use std::collections::{BTreeSet, HashMap};
use nalgebra::Vector3;

use crate::swfg::engine::accumulator::{DdPseudorangeAccumulator, DdSatKey};
use crate::swfg::pipeline::dd_factors::{
    DdCarrierPhaseFactor, DdDopplerFactor, DdPseudorangeFactor, WidelaneConstraintFactor,
};
use crate::swfg::pipeline::passes::CorrectedObservation;
use crate::swfg::pipeline::{
    build_carrier_phase_factor, build_pseudorange_factor, RawObservation,
};
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableKind};

/// Measurements collected for post-fit residual checking.
pub struct CpMeasurementRecord {
    pub sat: u16,
    pub var_amb: VariableId,
    pub dd_cp_m: f64,
    pub base_dd: f64,
    pub sat_pos: Vector3<f64>,
    pub ref_pos: Vector3<f64>,
    pub lambda: f64,
}

/// Context for building RTK double-differenced factors.
pub struct RtkFactorContext<'a> {
    pub epoch: u32,
    pub pose_id: VariableId,
    pub prev_pose_id: Option<VariableId>,
    pub init_pos: Vector3<f64>,
    pub prev_position: Option<Vector3<f64>>,
    pub dt_sec: f64,
    pub base_pos: Vector3<f64>,
    pub corrected: &'a [CorrectedObservation],
    pub base_raw: &'a [RawObservation],
}

/// Selects reference satellite per constellation with hysteresis.
pub fn select_ref_satellite(
    const_obs: &[&CorrectedObservation],
    cid: u8,
    ref_sat_map: &HashMap<u8, u16>,
) -> u16 {
    let best_obs = match const_obs
        .iter()
        .max_by(|a, b| a.elevation_rad.total_cmp(&b.elevation_rad))
    {
        Some(b) => b,
        None => return 0,
    };
    if let Some(&prev_ref) = ref_sat_map.get(&cid) {
        let prev_ref_obs = const_obs.iter().find(|o| o.satellite == prev_ref);
        match prev_ref_obs {
            Some(prev) if prev.elevation_rad >= 0.26 => prev_ref, // 15° — keep
            Some(_) if best_obs.elevation_rad < 0.52 => prev_ref, // no candidate > 30°
            _ => best_obs.satellite,
        }
    } else {
        best_obs.satellite
    }
}

/// Build undifferenced PR and CP factors for SPP/PPP mode.
pub fn build_undifferenced_factors(
    solver: &mut SlidingWindowSolver,
    corrected: &[CorrectedObservation],
    epoch: u32,
    pose_id: VariableId,
    zwd_id: Option<VariableId>,
) {
    for obs in corrected {
        let clock_id = solver
            .graph
            .variables
            .iter()
            .find(|(_, n)| {
                matches!(
                    n.kind,
                    VariableKind::ClockBias {
                        epoch: e,
                        constellation_id: c,
                    } if e == epoch && c == obs.constellation_id
                )
            })
            .map(|(id, _)| *id);

        // Lazy find-or-create: a GLONASS observation reaching this point
        // has already survived ephemeris matching, so the variable this
        // returns is guaranteed to get a factor below. See
        // SlidingWindowSolver::ensure_ifb_glonass.
        let var_ifb = if obs.constellation_id == 1 {
            Some(solver.ensure_ifb_glonass())
        } else {
            None
        };
        let pr_factor = build_pseudorange_factor(
            obs, epoch, pose_id, clock_id, zwd_id, var_ifb,
        );
        solver.graph.add_factor(pr_factor);

        if let Some(cp_l1) = obs.cp_l1 {
            let amb_id = solver.ensure_ambiguity(obs.satellite, 1);
            let lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / obs.f1.max(1.0);
            let current_val = solver
                .graph
                .variables
                .get(&amb_id)
                .map_or(0.0, |n| n.value[0]);
            if current_val == 0.0 {
                let float_amb = (cp_l1 * lambda - obs.pr_l1) / lambda;
                solver.graph.set_value(amb_id, &[float_amb]);
                let amb_prior = crate::swfg::factor::PriorFactor::new(
                    amb_id,
                    nalgebra::DVector::from_element(1, float_amb),
                    4.0,
                );
                solver.graph.add_factor(Box::new(amb_prior));
            }
            let mut cp_obs = obs.clone();
            cp_obs.cp_variance_m2 = 2.5e-3;
            let cp_factor = build_carrier_phase_factor(
                &cp_obs, epoch, pose_id, clock_id, zwd_id, amb_id, var_ifb, 0.0,
            );
            solver.graph.add_factor(cp_factor);
        }
    }
}

/// Helper to build L1 and L2 carrier phase factors for a single satellite.
#[allow(clippy::too_many_arguments)]
fn build_sat_carrier_phase(
    solver: &mut SlidingWindowSolver,
    obs: &CorrectedObservation,
    ref_obs: &CorrectedObservation,
    base_sat: &RawObservation,
    base_ref: &RawObservation,
    ctx: &RtkFactorContext<'_>,
    arc: u32,
    mw_acc: &mut DdPseudorangeAccumulator,
    cp_records: &mut Vec<CpMeasurementRecord>,
) {
    let (cp_rov_s, cp_rov_r) = match (obs.cp_l1, ref_obs.cp_l1) {
        (Some(s), Some(r)) => (s, r),
        _ => return,
    };
    let (cp_bas_s, cp_bas_r) = match (base_sat.cp_l1, base_ref.cp_l1) {
        (Some(s), Some(r)) => (s, r),
        _ => return,
    };

    let l1_lambda_sat = gneiss_core::constants::SPEED_OF_LIGHT_M_S / obs.f1;
    let l1_lambda_ref = gneiss_core::constants::SPEED_OF_LIGHT_M_S / ref_obs.f1;
    let dd_cp_m = (cp_rov_s * l1_lambda_sat - cp_rov_r * l1_lambda_ref)
        - (cp_bas_s * l1_lambda_sat - cp_bas_r * l1_lambda_ref);

    let cur_pos = match ctx.prev_position {
        Some(prev) if (prev - ctx.init_pos).norm() < 30.0 => prev,
        _ => ctx.init_pos,
    };
    let rover_sat_range = (obs.sat_pos_ecef - cur_pos).norm();
    let rover_ref_range = (ref_obs.sat_pos_ecef - cur_pos).norm();
    let base_range_sat = (obs.sat_pos_ecef - ctx.base_pos).norm();
    let base_range_ref = (ref_obs.sat_pos_ecef - ctx.base_pos).norm();
    let base_dd_range = base_range_sat - base_range_ref;

    let expected_dd_cp_m = (rover_sat_range - rover_ref_range) - base_dd_range;
    let initial_amb = (dd_cp_m - expected_dd_cp_m) / l1_lambda_sat;

    let (var_amb, is_new) = solver.ensure_dd_ambiguity(
        obs.constellation_id, obs.satellite, ref_obs.satellite, 1, arc, initial_amb,
    );

    let dd_cp_factor = DdCarrierPhaseFactor {
        var_pose: ctx.pose_id,
        var_amb,
        dd_cp_obs_m: dd_cp_m,
        sat_pos: obs.sat_pos_ecef,
        ref_pos: ref_obs.sat_pos_ecef,
        base_pos: ctx.base_pos,
        base_dd_range,
        lambda: l1_lambda_sat,
        variance_m2: 2.0 * obs.cp_variance_m2,
        elevation_rad: obs.elevation_rad,
        ref_elevation_rad: ref_obs.elevation_rad,
        is_new_amb: is_new,
        variables: vec![ctx.pose_id, var_amb],
    };
    solver.graph.add_factor(Box::new(dd_cp_factor));

    build_sat_l2_carrier_phase(
        solver, obs, ref_obs, base_sat, base_ref, ctx, arc,
        var_amb, dd_cp_m, base_dd_range, expected_dd_cp_m, mw_acc, cp_records,
    );
}

/// Helper to build L2 carrier phase factor and MW widelane constraint.
#[allow(clippy::too_many_arguments)]
fn build_sat_l2_carrier_phase(
    solver: &mut SlidingWindowSolver,
    obs: &CorrectedObservation,
    ref_obs: &CorrectedObservation,
    base_sat: &RawObservation,
    base_ref: &RawObservation,
    ctx: &RtkFactorContext<'_>,
    arc: u32,
    var_amb1: VariableId,
    dd_cp_m: f64,
    base_dd_range: f64,
    expected_dd_cp_m: f64,
    mw_acc: &mut DdPseudorangeAccumulator,
    cp_records: &mut Vec<CpMeasurementRecord>,
) {
    let (cp_rov_s2, cp_rov_r2) = match (obs.cp_l2, ref_obs.cp_l2) {
        (Some(s), Some(r)) => (s, r),
        _ => return,
    };
    let (cp_bas_s2, cp_bas_r2) = match (base_sat.cp_l2, base_ref.cp_l2) {
        (Some(s), Some(r)) => (s, r),
        _ => return,
    };

    let l2_lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / obs.f2;
    let dd_cp_m2 = ((cp_rov_s2 - cp_rov_r2) - (cp_bas_s2 - cp_bas_r2)) * l2_lambda;
    let initial_amb2 = (dd_cp_m2 - expected_dd_cp_m) / l2_lambda;

    let (var_amb2, is_new2) = solver.ensure_dd_ambiguity(
        obs.constellation_id, obs.satellite, ref_obs.satellite, 2, arc, initial_amb2,
    );

    let dd_cp_factor2 = DdCarrierPhaseFactor {
        var_pose: ctx.pose_id,
        var_amb: var_amb2,
        dd_cp_obs_m: dd_cp_m2,
        sat_pos: obs.sat_pos_ecef,
        ref_pos: ref_obs.sat_pos_ecef,
        base_pos: ctx.base_pos,
        base_dd_range,
        lambda: l2_lambda,
        variance_m2: 2.0 * obs.cp_variance_m2,
        elevation_rad: obs.elevation_rad,
        ref_elevation_rad: ref_obs.elevation_rad,
        is_new_amb: is_new2,
        variables: vec![ctx.pose_id, var_amb2],
    };
    solver.graph.add_factor(Box::new(dd_cp_factor2));

    // MW Widelane constraint
    if let (Some(pr_s2), Some(pr_r2), Some(b_s2), Some(b_r2)) = (
        obs.pr_l2, ref_obs.pr_l2, base_sat.pr_l2, base_ref.pr_l2,
    ) {
        if (obs.f1 - ref_obs.f1).abs() < 1.0 {
            let dd_pr1 = (obs.pr_l1 - ref_obs.pr_l1) - (base_sat.pr_l1 - base_ref.pr_l1);
            let dd_pr2 = (pr_s2 - pr_r2) - (b_s2 - b_r2);
            let mw_dd = crate::measurements::combinations::melbourne_wubbena(
                dd_cp_m, dd_cp_m2, dd_pr1, dd_pr2, obs.f1, obs.f2,
            );
            let lambda_wl = crate::measurements::combinations::lambda_wl(obs.f1, obs.f2);
            let n_wl = mw_dd / lambda_wl;
            let mw_key = DdSatKey {
                sat: obs.satellite,
                ref_sat: ref_obs.satellite,
                frequency: 0,
            };
            mw_acc.add_observation(mw_key, n_wl);
            if let Some((mean_wl, std_dev, count)) = mw_acc.get_stats(&mw_key) {
                if count >= 3 && std_dev < 0.25 {
                    let wl_fac = WidelaneConstraintFactor {
                        var_amb1,
                        var_amb2,
                        fixed_n_wl: mean_wl.round(),
                        variance: 1e-6,
                        variables: vec![var_amb1, var_amb2],
                    };
                    solver.graph.add_factor(Box::new(wl_fac));
                }
            }
        }
    }

    cp_records.push(CpMeasurementRecord {
        sat: obs.satellite,
        var_amb: var_amb2,
        dd_cp_m: dd_cp_m2,
        base_dd: base_dd_range,
        sat_pos: obs.sat_pos_ecef,
        ref_pos: ref_obs.sat_pos_ecef,
        lambda: l2_lambda,
    });
}

/// Builds all double-differenced factors for an RTK epoch.
#[allow(clippy::too_many_arguments)]
pub fn build_rtk_dd_factors(
    solver: &mut SlidingWindowSolver,
    ctx: &RtkFactorContext<'_>,
    ref_sat_map: &mut HashMap<u8, u16>,
    slip_counts: &mut HashMap<u16, u32>,
    mw_acc: &mut DdPseudorangeAccumulator,
) -> Vec<CpMeasurementRecord> {
    let mut cp_records = Vec::new();
    let constellations: BTreeSet<u8> = ctx.corrected.iter().map(|o| o.constellation_id).collect();

    for cid in constellations {
        let const_obs: Vec<_> = ctx.corrected.iter().filter(|o| o.constellation_id == cid).collect();
        if const_obs.len() < 2 { continue; }

        let ref_sat_id = select_ref_satellite(&const_obs, cid, ref_sat_map);
        if let Some(&prev_ref) = ref_sat_map.get(&cid) {
            if prev_ref != ref_sat_id {
                let stale_ids: Vec<_> = solver.graph.variables.iter()
                    .filter(|(_, n)| matches!(n.kind, VariableKind::DdAmbiguity { constellation_id, .. } if constellation_id == cid))
                    .map(|(id, _)| *id).collect();
                for id in stale_ids { solver.graph.remove_variable(id); }
            }
        }
        ref_sat_map.insert(cid, ref_sat_id);

        let ref_obs = match const_obs.iter().find(|o| o.satellite == ref_sat_id) {
            Some(o) => *o,
            None => continue,
        };
        let base_ref = match ctx.base_raw.iter().find(|b| b.satellite == ref_obs.satellite && b.constellation_id == cid && b.pr_l1 > 0.0) {
            Some(b) => b,
            None => continue,
        };

        for obs in &const_obs {
            if obs.satellite == ref_obs.satellite { continue; }
            let base_sat = match ctx.base_raw.iter().find(|b| b.satellite == obs.satellite && b.constellation_id == cid && b.pr_l1 > 0.0) {
                Some(b) => b,
                None => continue,
            };

            let dd_pr = (obs.pr_l1 - ref_obs.pr_l1) - (base_sat.pr_l1 - base_ref.pr_l1);
            let base_dd_range = (obs.sat_pos_ecef - ctx.base_pos).norm() - (ref_obs.sat_pos_ecef - ctx.base_pos).norm();

            let dd_pr_factor = DdPseudorangeFactor {
                var_pose: ctx.pose_id,
                dd_pr_obs: dd_pr,
                sat_pos: obs.sat_pos_ecef,
                ref_pos: ref_obs.sat_pos_ecef,
                base_pos: ctx.base_pos,
                base_dd_range,
                variance_m2: 2.0 * obs.variance_m2,
                elevation_rad: obs.elevation_rad,
                ref_elevation_rad: ref_obs.elevation_rad,
                variables: vec![ctx.pose_id],
            };
            solver.graph.add_factor(Box::new(dd_pr_factor));

            // Doppler factor
            if let Some(prev_p) = ctx.prev_pose_id {
                if solver.graph.variables.contains_key(&prev_p) && ctx.dt_sec <= 2.0 && obs.doppler.abs() > 1.0 && ref_obs.doppler.abs() > 1.0 {
                    let l1_lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / obs.f1;
                    let rov_dop = -(obs.doppler - ref_obs.doppler) * l1_lambda;
                    let bas_dop = -(base_sat.doppler - base_ref.doppler) * l1_lambda;
                    let dop_factor = DdDopplerFactor {
                        var_pose_prev: prev_p,
                        var_pose_curr: ctx.pose_id,
                        dd_doppler_m_s: rov_dop - bas_dop,
                        dt: ctx.dt_sec,
                        sat_pos: obs.sat_pos_ecef,
                        ref_pos: ref_obs.sat_pos_ecef,
                        variance_m2: 0.05,
                        variables: vec![prev_p, ctx.pose_id],
                    };
                    solver.graph.add_factor(Box::new(dop_factor));
                }
            }

            let has_slip = |lli: Option<u8>| lli.unwrap_or(0) & 1 != 0;
            if has_slip(obs.cp_l1_lli) || has_slip(ref_obs.cp_l1_lli) || has_slip(base_sat.cp_l1_lli) || has_slip(base_ref.cp_l1_lli) {
                *slip_counts.entry(obs.satellite).or_insert(0) += 1;
            }
            let arc = *slip_counts.get(&obs.satellite).unwrap_or(&0);

            build_sat_carrier_phase(
                solver, obs, ref_obs, base_sat, base_ref, ctx, arc,
                mw_acc, &mut cp_records,
            );
        }
    }
    cp_records
}
