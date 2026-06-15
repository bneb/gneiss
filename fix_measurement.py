import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    code = f.read()

b1_old = """struct UpdateGeometry {
    comp_dd: f64,
    h_r: Vector3<f64>,
    h_att: Vector3<f64>,
    h_zwd: f64,
    state_size: usize,
}

struct VarianceWeights {
    val: f64,
    ref_val: f64,
}"""
b1_new = """pub struct UpdateGeometry {
    pub comp_dd: f64,
    pub h_r: Vector3<f64>,
    pub h_att: Vector3<f64>,
    pub h_zwd: f64,
    pub state_size: usize,
}

pub struct VarianceWeights {
    pub val: f64,
    pub ref_val: f64,
}

pub struct DdMeasurementContext<'a> {
    pub ctx: &'a DdContext<'a>,
    pub geom: &'a UpdateGeometry,
    pub comps: &'a DdComponents,
    pub env: &'a MeasurementEnvironment<'a>,
}

pub struct DdCarrierPhaseParams<'a> {
    pub is_fixed: bool,
    pub ambiguities: &'a [f64],
    pub sat_idx_l1: Option<usize>,
    pub ref_idx_l1: Option<usize>,
    pub sat_idx_l2: Option<usize>,
    pub ref_idx_l2: Option<usize>,
    pub cp_base_var: f64,
}"""

b2_old = """pub fn compute_dd_pseudorange(
    ctx: &DdContext, comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, state_size: usize,
    var_factor: f64, env: &MeasurementEnvironment,
    h_zwd: f64, ref_var_factor: f64
) -> Vec<SingleUpdate> {
    let mut updates = Vec::new();
    let geom = UpdateGeometry { comp_dd: comp_pr_dd, h_r, h_att, h_zwd, state_size };
    let var = VarianceWeights { val: env.tuning.pr_base_var * var_factor, ref_val: env.tuning.pr_base_var * ref_var_factor };

    if [ctx.rov_sat.pr_l1, ctx.base_sat.pr_l1, ctx.rov_ref.pr_l1, ctx.ref_base.pr_l1].iter().all(|&x| x > 0.0) {
        updates.push(compute_pseudorange_update(
            [ctx.rov_sat.pr_l1, ctx.rov_ref.pr_l1, ctx.base_sat.pr_l1, ctx.ref_base.pr_l1], iono_dd_l1, &geom, &var
        ));
    }

    if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [ctx.rov_ref.pr_l2, ctx.rov_sat.pr_l2, ctx.ref_base.pr_l2, ctx.base_sat.pr_l2] {
        updates.push(compute_pseudorange_update([rs2, rr2, bs2, br2], iono_dd_l2, &geom, &var));
    }
    updates
}"""
b2_new = """pub fn compute_dd_pseudorange(mctx: &DdMeasurementContext) -> Vec<SingleUpdate> {
    let mut updates = Vec::new();
    let var = VarianceWeights { 
        val: mctx.env.tuning.pr_base_var * mctx.comps.var_factor, 
        ref_val: mctx.env.tuning.pr_base_var * mctx.comps.ref_var_factor 
    };

    if [mctx.ctx.rov_sat.pr_l1, mctx.ctx.base_sat.pr_l1, mctx.ctx.rov_ref.pr_l1, mctx.ctx.ref_base.pr_l1].iter().all(|&x| x > 0.0) {
        updates.push(compute_pseudorange_update(
            [mctx.ctx.rov_sat.pr_l1, mctx.ctx.rov_ref.pr_l1, mctx.ctx.base_sat.pr_l1, mctx.ctx.ref_base.pr_l1], mctx.comps.iono_dd_l1, mctx.geom, &var
        ));
    }

    if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [mctx.ctx.rov_ref.pr_l2, mctx.ctx.rov_sat.pr_l2, mctx.ctx.ref_base.pr_l2, mctx.ctx.base_sat.pr_l2] {
        updates.push(compute_pseudorange_update([rs2, rr2, bs2, br2], mctx.comps.iono_dd_l2, mctx.geom, &var));
    }
    updates
}"""

b3_old = """pub fn compute_dd_carrier_phase(
    ctx: &DdContext, is_fixed: bool, ambiguities: &[f64],
    sat_idx_l1: Option<usize>, ref_idx_l1: Option<usize>, sat_idx_l2: Option<usize>, ref_idx_l2: Option<usize>,
    comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64, h_r: Vector3<f64>, h_att: Vector3<f64>,
    state_size: usize, var_factor: f64, ref_var_factor: f64, cp_base_var: f64, h_zwd: f64
) -> Vec<SingleUpdate> {
    let mut updates = Vec::new();
    let geom = UpdateGeometry { comp_dd: comp_pr_dd, h_r, h_att, h_zwd, state_size };
    let var = VarianceWeights {
        val: if is_fixed { 1e-6 * var_factor } else { cp_base_var * var_factor },
        ref_val: if is_fixed { 1e-6 * ref_var_factor } else { cp_base_var * ref_var_factor }
    };

    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l1, ref_idx_l1) {
        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = [ctx.rov_ref.cp_l1, ctx.rov_sat.cp_l1, ctx.ref_base.cp_l1, ctx.base_sat.cp_l1] {
            updates.push(compute_carrier_phase_update(
                [rs1, rr1, bs1, br1], [ctx.sat_state.f1, ctx.ref_state.f1], [sat_idx, ref_idx], ambiguities,
                iono_dd_l1, &geom, &var, 1
            ));
        }
    }

    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l2, ref_idx_l2) {
        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [ctx.rov_ref.cp_l2, ctx.rov_sat.cp_l2, ctx.ref_base.cp_l2, ctx.base_sat.cp_l2] {
            updates.push(compute_carrier_phase_update(
                [rs2, rr2, bs2, br2], [ctx.sat_state.f2, ctx.ref_state.f2], [sat_idx, ref_idx], ambiguities,
                iono_dd_l2, &geom, &var, 2
            ));
        }
    }
    updates
}"""
b3_new = """pub fn compute_dd_carrier_phase(mctx: &DdMeasurementContext, p: &DdCarrierPhaseParams) -> Vec<SingleUpdate> {
    let mut updates = Vec::new();
    let var = VarianceWeights {
        val: if p.is_fixed { 1e-6 * mctx.comps.var_factor } else { p.cp_base_var * mctx.comps.var_factor },
        ref_val: if p.is_fixed { 1e-6 * mctx.comps.ref_var_factor } else { p.cp_base_var * mctx.comps.ref_var_factor }
    };

    if let (Some(sat_idx), Some(ref_idx)) = (p.sat_idx_l1, p.ref_idx_l1) {
        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = [mctx.ctx.rov_ref.cp_l1, mctx.ctx.rov_sat.cp_l1, mctx.ctx.ref_base.cp_l1, mctx.ctx.base_sat.cp_l1] {
            updates.push(compute_carrier_phase_update(
                [rs1, rr1, bs1, br1], [mctx.ctx.sat_state.f1, mctx.ctx.ref_state.f1], [sat_idx, ref_idx], p.ambiguities,
                mctx.comps.iono_dd_l1, mctx.geom, &var, 1
            ));
        }
    }

    if let (Some(sat_idx), Some(ref_idx)) = (p.sat_idx_l2, p.ref_idx_l2) {
        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [mctx.ctx.rov_ref.cp_l2, mctx.ctx.rov_sat.cp_l2, mctx.ctx.ref_base.cp_l2, mctx.ctx.base_sat.cp_l2] {
            updates.push(compute_carrier_phase_update(
                [rs2, rr2, bs2, br2], [mctx.ctx.sat_state.f2, mctx.ctx.ref_state.f2], [sat_idx, ref_idx], p.ambiguities,
                mctx.comps.iono_dd_l2, mctx.geom, &var, 2
            ));
        }
    }
    updates
}"""

b4_old = """    let DdComponents { comp_pr_dd, iono_dd_l1, iono_dd_l2, var_factor, ref_var_factor, h_zwd } = 
        compute_dd_components(state, geom, env, &ctx);

    generate_measurement_updates(state, &ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom, var_factor, env, h_zwd, ref_var_factor, ref_idx_l1, ref_idx_l2, updates);"""
b4_new = """    let comps = compute_dd_components(state, geom, env, &ctx);
    let ugeom = UpdateGeometry { comp_dd: comps.comp_pr_dd, h_r, h_att, h_zwd: comps.h_zwd, state_size: geom.state_size };
    let mctx = DdMeasurementContext { ctx: &ctx, geom: &ugeom, comps: &comps, env };

    generate_measurement_updates(state, &mctx, geom, ref_idx_l1, ref_idx_l2, updates);"""

b5_old = """#[allow(clippy::too_many_arguments)]
fn generate_measurement_updates(
    state: &mut RtkState, ctx: &DdContext, comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, geom: &EkfGeometryContext, var_factor: f64,
    env: &MeasurementEnvironment, h_zwd: f64, ref_var_factor: f64,
    ref_idx_l1: Option<usize>, ref_idx_l2: Option<usize>, updates: &mut EkfUpdates
) {
    let sat = ctx.rov_sat.sat;
    for u in compute_dd_pseudorange(ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size, var_factor, env, h_zwd, ref_var_factor) {
        updates.push(u, sat);
    }

    let sat_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 1);
    let sat_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 2);
    for u in compute_dd_carrier_phase(ctx, state.is_fixed, &state.ambiguities, sat_idx_l1, ref_idx_l1, sat_idx_l2, ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size, var_factor, ref_var_factor, env.tuning.cp_base_var, h_zwd) {
        updates.push(u, sat);
    }

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(u) = compute_dd_doppler(ctx, geom.pos_apc, geom.base_coord_vec, h_r, &r_b_e_rot, &state.velocity, &env.omega_b, &env.lever_arm, env.tuning.dop_base_var, geom.state_size, var_factor, ref_var_factor) {
        updates.push(u, sat);
    }
}"""
b5_new = """fn generate_measurement_updates(
    state: &RtkState, mctx: &DdMeasurementContext, geom: &EkfGeometryContext,
    ref_idx_l1: Option<usize>, ref_idx_l2: Option<usize>, updates: &mut EkfUpdates
) {
    let sat = mctx.ctx.rov_sat.sat;
    for u in compute_dd_pseudorange(mctx) {
        updates.push(u, sat);
    }

    let p = DdCarrierPhaseParams {
        is_fixed: state.is_fixed,
        ambiguities: &state.ambiguities,
        sat_idx_l1: state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 1),
        ref_idx_l1,
        sat_idx_l2: state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 2),
        ref_idx_l2,
        cp_base_var: mctx.env.tuning.cp_base_var,
    };
    
    for u in compute_dd_carrier_phase(mctx, &p) {
        updates.push(u, sat);
    }

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(u) = compute_dd_doppler(mctx.ctx, geom.pos_apc, geom.base_coord_vec, mctx.geom.h_r, &r_b_e_rot, &state.velocity, &mctx.env.omega_b, &mctx.env.lever_arm, mctx.env.tuning.dop_base_var, geom.state_size, mctx.comps.var_factor, mctx.comps.ref_var_factor) {
        updates.push(u, sat);
    }
}"""

b6_old = """        let updates = compute_dd_pseudorange(
            &ctx, 0.0, 0.0, 0.0, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros(), 22, 1.0, &env, 0.0, 1.0
        );"""
b6_new = """        let ugeom = crate::engine::measurement::UpdateGeometry { comp_dd: 0.0, h_r: Vector3::new(1.0, 0.0, 0.0), h_att: Vector3::zeros(), h_zwd: 0.0, state_size: 22 };
        let comps = crate::engine::measurement::DdComponents { comp_pr_dd: 0.0, iono_dd_l1: 0.0, iono_dd_l2: 0.0, var_factor: 1.0, ref_var_factor: 1.0, h_zwd: 0.0 };
        let mctx = crate::engine::measurement::DdMeasurementContext { ctx: &ctx, geom: &ugeom, comps: &comps, env: &env };
        let updates = compute_dd_pseudorange(&mctx);"""

b7_old = """        let updates = compute_dd_carrier_phase(
            &ctx, false, &state.ambiguities, Some(0), Some(1), Some(2), Some(3),
            0.0, 0.0, 0.0, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros(), 22, 1.0, 1.0, 0.01, 0.0
        );"""
b7_new = """        let ugeom = crate::engine::measurement::UpdateGeometry { comp_dd: 0.0, h_r: Vector3::new(1.0, 0.0, 0.0), h_att: Vector3::zeros(), h_zwd: 0.0, state_size: 22 };
        let comps = crate::engine::measurement::DdComponents { comp_pr_dd: 0.0, iono_dd_l1: 0.0, iono_dd_l2: 0.0, var_factor: 1.0, ref_var_factor: 1.0, h_zwd: 0.0 };
        let mctx = crate::engine::measurement::DdMeasurementContext { ctx: &ctx, geom: &ugeom, comps: &comps, env: &env };
        let p = crate::engine::measurement::DdCarrierPhaseParams { is_fixed: false, ambiguities: &state.ambiguities, sat_idx_l1: Some(0), ref_idx_l1: Some(1), sat_idx_l2: Some(2), ref_idx_l2: Some(3), cp_base_var: 0.01 };
        let updates = compute_dd_carrier_phase(&mctx, &p);"""

code = code.replace(b1_old, b1_new)
code = code.replace(b2_old, b2_new)
code = code.replace(b3_old, b3_new)
code = code.replace(b4_old, b4_new)
code = code.replace(b5_old, b5_new)
code = code.replace(b6_old, b6_new)
code = code.replace(b7_old, b7_new)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(code)

print("done")
