import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_pr = """pub fn compute_dd_pseudorange(
    ctx: &DdContext,
    comp_pr_dd: f64,
    iono_dd_l1: f64,
    iono_dd_l2: f64,
    h_r: Vector3<f64>,
    h_att: Vector3<f64>,
    state_size: usize,
    var_factor: f64,
    env: &MeasurementEnvironment,
 h_zwd: f64, ref_var_factor: f64) -> Vec<(f64, Vec<f64>, f64, u8, f64)> {
    let mut updates = Vec::new();
    let pr_base_var = env.tuning.pr_base_var;

    if [ctx.rov_sat.pr_l1, ctx.base_sat.pr_l1, ctx.rov_ref.pr_l1, ctx.ref_base.pr_l1].iter().all(|&x| x > 0.0) {
        let pr_dd = (ctx.rov_sat.pr_l1 - ctx.rov_ref.pr_l1) - (ctx.base_sat.pr_l1 - ctx.ref_base.pr_l1);
        let mut h_pr1 = vec![0.0; state_size];
        h_pr1[0] = h_r.x; h_pr1[1] = h_r.y; h_pr1[2] = h_r.z;
        h_pr1[6] = h_att.x; h_pr1[7] = h_att.y; h_pr1[8] = h_att.z;
        if state_size > 20 { h_pr1[20] = h_zwd; }
        updates.push((pr_dd - (comp_pr_dd + iono_dd_l1), h_pr1, pr_base_var * var_factor, 0, pr_base_var * ref_var_factor));
    }

    if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [ctx.rov_ref.pr_l2, ctx.rov_sat.pr_l2, ctx.ref_base.pr_l2, ctx.base_sat.pr_l2] {
        let pr_dd_l2 = (rs2 - rr2) - (bs2 - br2);
        let mut h_pr2 = vec![0.0; state_size];
        h_pr2[0] = h_r.x; h_pr2[1] = h_r.y; h_pr2[2] = h_r.z;
        h_pr2[6] = h_att.x; h_pr2[7] = h_att.y; h_pr2[8] = h_att.z;
        if state_size > 20 { h_pr2[20] = h_zwd; }
        updates.push((pr_dd_l2 - (comp_pr_dd + iono_dd_l2), h_pr2, pr_base_var * var_factor, 0, pr_base_var * ref_var_factor));
    }

    updates
}"""

new_pr = """fn compute_pseudorange_update(
    rs: f64, rr: f64, bs: f64, br: f64, comp_pr_dd: f64, iono_dd: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, h_zwd: f64,
    state_size: usize, r_val: f64, r_ref_val: f64
) -> (f64, Vec<f64>, f64, u8, f64) {
    let pr_dd = (rs - rr) - (bs - br);
    let mut h_pr = vec![0.0; state_size];
    h_pr[0] = h_r.x; h_pr[1] = h_r.y; h_pr[2] = h_r.z;
    h_pr[6] = h_att.x; h_pr[7] = h_att.y; h_pr[8] = h_att.z;
    if state_size > 20 { h_pr[20] = h_zwd; }
    (pr_dd - (comp_pr_dd + iono_dd), h_pr, r_val, 0, r_ref_val)
}

#[allow(clippy::too_many_arguments)]
pub fn compute_dd_pseudorange(
    ctx: &DdContext, comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, state_size: usize,
    var_factor: f64, env: &MeasurementEnvironment,
    h_zwd: f64, ref_var_factor: f64
) -> Vec<(f64, Vec<f64>, f64, u8, f64)> {
    let mut updates = Vec::new();
    let r_val = env.tuning.pr_base_var * var_factor;
    let r_ref = env.tuning.pr_base_var * ref_var_factor;

    if [ctx.rov_sat.pr_l1, ctx.base_sat.pr_l1, ctx.rov_ref.pr_l1, ctx.ref_base.pr_l1].iter().all(|&x| x > 0.0) {
        updates.push(compute_pseudorange_update(
            ctx.rov_sat.pr_l1, ctx.rov_ref.pr_l1, ctx.base_sat.pr_l1, ctx.ref_base.pr_l1,
            comp_pr_dd, iono_dd_l1, h_r, h_att, h_zwd, state_size, r_val, r_ref
        ));
    }

    if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [ctx.rov_ref.pr_l2, ctx.rov_sat.pr_l2, ctx.ref_base.pr_l2, ctx.base_sat.pr_l2] {
        updates.push(compute_pseudorange_update(
            rs2, rr2, bs2, br2, comp_pr_dd, iono_dd_l2, h_r, h_att, h_zwd, state_size, r_val, r_ref
        ));
    }
    updates
}"""

content = content.replace(old_pr, new_pr)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
