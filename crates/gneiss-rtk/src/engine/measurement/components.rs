use crate::engine::measurement::types::{
    DdCarrierPhaseParams, DdContext, DdMeasurementContext,
};
use crate::engine::measurement::{
    compute_dd_carrier_phase, compute_dd_doppler, compute_dd_pseudorange,
};
use crate::engine::measurement::EkfGeometryContext;
use crate::engine::measurement_math;
use crate::filter::RtkState;

pub struct DdComponents {
    pub comp_pr_dd: f64,
    pub iono_dd_l1: f64,
    pub iono_dd_l2: f64,
    pub var_factor: f64,
    pub ref_var_factor: f64,
    pub h_zwd: f64,
}

pub(crate) fn compute_dd_components(
    state: &RtkState,
    geom: &EkfGeometryContext,
    env: &crate::engine::measurement::types::MeasurementEnvironment,
    ctx: &DdContext,
) -> DdComponents {
    let (tropo_dd, iono_dd_l1, iono_dd_l2) = measurement_math::compute_atmospheric_delays(
        state.time,
        geom.pos_apc,
        geom.base_coord_vec,
        ctx.sat_state.rov_pos,
        ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos,
        ctx.ref_state.bas_pos,
        ctx.sat_state.f1,
        ctx.sat_state.f2,
        ctx.ref_state.f1,
        ctx.ref_state.f2,
        env.klobuchar_params.as_ref(),
    );

    let base_llh = gneiss_core::coords::ecef_to_llh(geom.base_coord_vec);
    let rov_llh = gneiss_core::coords::ecef_to_llh(geom.pos_apc);
    let (_, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, ctx.sat_state.rov_pos);
    let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, ctx.ref_state.rov_pos);
    let (_, el_bas_sat) =
        gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, ctx.sat_state.bas_pos);
    let (_, el_bas_ref) =
        gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, ctx.ref_state.bas_pos);

    let baseline_dist = (geom.pos_apc - geom.base_coord_vec).norm();
    let (var_factor, ref_var_factor) = measurement_math::compute_variance_factors(
        &measurement_math::VarianceFactors {
            snr_rov_sat: ctx.rov_sat.snr,
            snr_rov_ref: ctx.rov_ref.snr,
            el_rov_sat,
            el_rov_ref,
            el_bas_sat,
            el_bas_ref,
            snr_a: env.tuning.snr_a,
            snr_b: env.tuning.snr_b,
            gnn_var_sat: env.gnn_variances.get(&ctx.rov_sat.sat).copied(),
            gnn_var_ref: env.gnn_variances.get(&ctx.rov_ref.sat).copied(),
            baseline_distance_m: baseline_dist,
        },
    );
    let (h_zwd, zwd_dd) = measurement_math::compute_zwd_mapping(el_rov_sat, el_rov_ref, state.zwd);

    let comp_pr_dd = measurement_math::compute_geometric_dd(
        geom.pos_apc,
        geom.base_coord_vec,
        ctx.sat_state.rov_pos,
        ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos,
        ctx.ref_state.bas_pos,
    ) + tropo_dd
        + zwd_dd;

    DdComponents {
        comp_pr_dd,
        iono_dd_l1,
        iono_dd_l2,
        var_factor,
        ref_var_factor,
        h_zwd,
    }
}

pub(crate) fn generate_measurement_updates(
    state: &RtkState,
    mctx: &DdMeasurementContext,
    geom: &EkfGeometryContext,
    ref_idx_l1: Option<usize>,
    ref_idx_l2: Option<usize>,
    updates: &mut crate::engine::measurement::types::EkfUpdates,
) {
    let sat = mctx.ctx.rov_sat.sat;
    for u in compute_dd_pseudorange(mctx) {
        updates.push(u, sat);
    }

    let p = DdCarrierPhaseParams {
        is_fixed: state.is_fixed,
        ambiguities: &state.ambiguities,
        sat_idx_l1: state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == sat && f == 1),
        ref_idx_l1,
        sat_idx_l2: state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == sat && f == 2),
        ref_idx_l2,
        iono_idx_sat: state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == sat && f == 3),
        iono_idx_ref: state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == mctx.ctx.rov_ref.sat && f == 3),
        iono_state_vals: &state.ambiguities,
        cp_base_var: mctx.env.tuning.cp_base_var,
    };

    for u in compute_dd_carrier_phase(mctx, &p) {
        updates.push(u, sat);
    }

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(u) = compute_dd_doppler(
        mctx.ctx,
        geom.pos_apc,
        geom.base_coord_vec,
        mctx.geom.h_r,
        &r_b_e_rot,
        &state.velocity,
        &mctx.env.omega_b,
        &mctx.env.lever_arm,
        mctx.env.tuning.dop_base_var,
        geom.state_size,
        mctx.comps.var_factor,
        mctx.comps.ref_var_factor,
    ) {
        updates.push(u, sat);
    }
}
