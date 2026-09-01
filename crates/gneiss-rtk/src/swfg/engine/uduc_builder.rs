//! Undifferenced Uncombined (UDUC) factor graph builder for PPP and PPP-AR modes.

use std::collections::HashMap;
use nalgebra::Vector3;

use crate::swfg::pipeline::factors::uduc::{
    GeodeticNormalizations, SlantIonoRandomWalkFactor, UducCarrierPhaseFactor, UducPseudorangeFactor,
};
use crate::swfg::pipeline::passes::CorrectedObservation;
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableKind};

pub struct UducFactorContext<'a> {
    pub epoch: u32,
    pub pose_id: VariableId,
    pub zwd_id: Option<VariableId>,
    pub prev_epoch: Option<u32>,
    pub dt_sec: f64,
    pub norm: GeodeticNormalizations,
    pub slip_counts: &'a mut HashMap<u16, u32>,
    pub windup_trackers: &'a mut HashMap<u16, gneiss_geodesy::windup::PhaseWindupTracker>,
    pub rover_time: gneiss_core::time::GpsTime,
    pub rx_pos: Vector3<f64>,
}

/// Build Undifferenced Uncombined (UDUC) factors for PPP-AR mode with slant ionosphere states.
pub fn build_uduc_factors(
    solver: &mut SlidingWindowSolver,
    corrected: &[CorrectedObservation],
    ctx: &mut UducFactorContext<'_>,
) {
    let (sun_pos, _) = gneiss_geodesy::tides::solar_lunar_positions(ctx.rover_time.tow, ctx.rover_time.week);
    let ref_llh = gneiss_core::coords::ecef_to_llh(ctx.rx_pos);
    let sin_lat = ref_llh.x.sin();
    let cos_lat = ref_llh.x.cos();
    let sin_lon = ref_llh.y.sin();
    let cos_lon = ref_llh.y.cos();

    let rx_up = Vector3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat);
    let rx_north = Vector3::new(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat);
    let rx_east = Vector3::new(-sin_lon, cos_lon, 0.0);

    for obs in corrected {
        let clock_id = solver
            .graph
            .variables
            .iter()
            .find(|(_, n)| {
                matches!(
                    n.kind,
                    VariableKind::ClockBias { epoch: e, constellation_id: c } if e == ctx.epoch && c == obs.constellation_id
                )
            })
            .map(|(id, _)| *id);

        let c_id = match clock_id {
            Some(c) => c,
            None => continue,
        };

        // Ensure per-satellite slant ionosphere state for this epoch
        let iono_kind = VariableKind::IonosphereSlant { epoch: ctx.epoch, satellite: obs.satellite };
        let iono_id = solver.graph.add_variable(iono_kind);
        solver.graph.set_value(iono_id, &[obs.iono_l1_m]);

        if let Some(pe) = ctx.prev_epoch {
            let prev_iono_kind = VariableKind::IonosphereSlant { epoch: pe, satellite: obs.satellite };
            if let Some((prev_id, _)) = solver.graph.variables.iter().find(|(_, n)| n.kind == prev_iono_kind) {
                let rw_factor = SlantIonoRandomWalkFactor::new(*prev_id, iono_id, ctx.dt_sec, 4.0e-4);
                solver.graph.add_factor(Box::new(rw_factor));
            }
        }

        let gamma = (obs.f1 / obs.f2.max(1.0)).powi(2);

        // L1 Pseudorange factor
        let pr_factor = UducPseudorangeFactor::new(
            obs.clone(), ctx.epoch, obs.pr_l1, ctx.pose_id, c_id, ctx.zwd_id, Some(iono_id), 1.0, ctx.norm,
        );
        solver.graph.add_factor(Box::new(pr_factor));

        // L2 Pseudorange factor
        if let Some(pr_l2) = obs.pr_l2 {
            let pr2_factor = UducPseudorangeFactor::new(
                obs.clone(), ctx.epoch, pr_l2, ctx.pose_id, c_id, ctx.zwd_id, Some(iono_id), gamma, ctx.norm,
            );
            solver.graph.add_factor(Box::new(pr2_factor));
        }

        let tracker = ctx.windup_trackers.entry(obs.satellite).or_default();
        let windup_rad = tracker.update(&obs.sat_pos_ecef, &sun_pos, &ctx.rx_pos, &rx_up, &rx_north, &rx_east);

        if obs.cp_l1_lli.unwrap_or(0) & 1 != 0 {
            *ctx.slip_counts.entry(obs.satellite).or_insert(0) += 1;
        }
        let arc = *ctx.slip_counts.entry(obs.satellite).or_insert(0);

        let sin_el = obs.elevation_rad.sin().max(0.1);
        let cp_var = (0.003 / sin_el).powi(2);

        // L1 Carrier phase factor
        if let Some(cp_l1) = obs.cp_l1 {
            let amb_id = solver.ensure_ambiguity(obs.satellite, 1, arc);
            let lambda1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / obs.f1.max(1.0);
            let windup_m1 = (windup_rad / (2.0 * std::f64::consts::PI)) * lambda1;

            let current_val = solver.graph.variables.get(&amb_id).map_or(0.0, |n| n.value[0]);
            if current_val == 0.0 {
                let float_amb = (cp_l1 * lambda1 - windup_m1 - obs.pr_l1) / lambda1;
                solver.graph.set_value(amb_id, &[float_amb]);
                let amb_prior = crate::swfg::factor::PriorFactor::new(
                    amb_id, nalgebra::DVector::from_element(1, float_amb), 10_000.0,
                );
                solver.graph.add_factor(Box::new(amb_prior));
            }
            let mut cp_obs = obs.clone();
            cp_obs.cp_variance_m2 = cp_var;
            let cp_factor = UducCarrierPhaseFactor::new(
                cp_obs, ctx.epoch, cp_l1, ctx.pose_id, c_id, ctx.zwd_id, Some(iono_id), amb_id, lambda1, 1.0, windup_m1, ctx.norm,
            );
            solver.graph.add_factor(Box::new(cp_factor));
        }

        // L2 Carrier phase factor
        if let Some(cp_l2) = obs.cp_l2 {
            let amb2_id = solver.ensure_ambiguity(obs.satellite, 2, arc);
            let lambda2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / obs.f2.max(1.0);
            let windup_m2 = (windup_rad / (2.0 * std::f64::consts::PI)) * lambda2;

            let current_val = solver.graph.variables.get(&amb2_id).map_or(0.0, |n| n.value[0]);
            let pr2_m = obs.pr_l2.unwrap_or(obs.pr_l1);
            if current_val == 0.0 {
                let float_amb2 = (cp_l2 * lambda2 - windup_m2 - pr2_m) / lambda2;
                solver.graph.set_value(amb2_id, &[float_amb2]);
                let amb2_prior = crate::swfg::factor::PriorFactor::new(
                    amb2_id, nalgebra::DVector::from_element(1, float_amb2), 10_000.0,
                );
                solver.graph.add_factor(Box::new(amb2_prior));
            }
            let mut cp2_obs = obs.clone();
            cp2_obs.cp_variance_m2 = cp_var;
            let cp2_factor = UducCarrierPhaseFactor::new(
                cp2_obs, ctx.epoch, cp_l2, ctx.pose_id, c_id, ctx.zwd_id, Some(iono_id), amb2_id, lambda2, gamma, windup_m2, ctx.norm,
            );
            solver.graph.add_factor(Box::new(cp2_factor));
        }
    }
}
