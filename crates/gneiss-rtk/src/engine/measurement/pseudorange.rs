use crate::engine::measurement::types::{SingleUpdate, UpdateGeometry, VarianceWeights};
use crate::engine::measurement::DdMeasurementContext;

fn compute_pseudorange_update(
    pr: [f64; 4],
    iono_dd: f64,
    geom: &UpdateGeometry,
    var: &VarianceWeights,
) -> SingleUpdate {
    let pr_dd = (pr[0] - pr[1]) - (pr[2] - pr[3]);
    let mut h_pr = vec![0.0; geom.state_size];
    h_pr[0] = geom.h_r.x;
    h_pr[1] = geom.h_r.y;
    h_pr[2] = geom.h_r.z;
    h_pr[6] = geom.h_att.x;
    h_pr[7] = geom.h_att.y;
    h_pr[8] = geom.h_att.z;
    if geom.state_size > 20 {
        h_pr[20] = geom.h_zwd;
    }
    SingleUpdate {
        z: pr_dd - (geom.comp_dd + iono_dd),
        h: h_pr,
        r: var.val,
        type_code: 0,
        r_ref: var.ref_val,
    }
}

/// Computes the L1 pseudorange double-difference update if all four L1
/// pseudoranges are valid (positive).
fn compute_l1_pseudorange(
    mctx: &DdMeasurementContext,
    var: &VarianceWeights,
) -> Option<SingleUpdate> {
    let valid = [
        mctx.ctx.rov_sat.pr_l1,
        mctx.ctx.base_sat.pr_l1,
        mctx.ctx.rov_ref.pr_l1,
        mctx.ctx.ref_base.pr_l1,
    ]
    .iter()
    .all(|&x| x > 0.0);
    if valid {
        Some(compute_pseudorange_update(
            [
                mctx.ctx.rov_sat.pr_l1,
                mctx.ctx.rov_ref.pr_l1,
                mctx.ctx.base_sat.pr_l1,
                mctx.ctx.ref_base.pr_l1,
            ],
            mctx.comps.iono_dd_l1,
            mctx.geom,
            var,
        ))
    } else {
        None
    }
}

/// Computes the L2 pseudorange double-difference update if all four L2
/// pseudoranges are available.
fn compute_l2_pseudorange(
    mctx: &DdMeasurementContext,
    var: &VarianceWeights,
) -> Option<SingleUpdate> {
    if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [
        mctx.ctx.rov_ref.pr_l2,
        mctx.ctx.rov_sat.pr_l2,
        mctx.ctx.ref_base.pr_l2,
        mctx.ctx.base_sat.pr_l2,
    ] {
        Some(compute_pseudorange_update(
            [rs2, rr2, bs2, br2],
            mctx.comps.iono_dd_l2,
            mctx.geom,
            var,
        ))
    } else {
        None
    }
}

pub fn compute_dd_pseudorange(mctx: &DdMeasurementContext) -> Vec<SingleUpdate> {
    let var = VarianceWeights {
        val: mctx.env.tuning.pr_base_var * mctx.comps.var_factor,
        ref_val: mctx.env.tuning.pr_base_var * mctx.comps.ref_var_factor,
    };
    let mut updates = Vec::new();
    if let Some(u) = compute_l1_pseudorange(mctx, &var) {
        updates.push(u);
    }
    if let Some(u) = compute_l2_pseudorange(mctx, &var) {
        updates.push(u);
    }
    updates
}
