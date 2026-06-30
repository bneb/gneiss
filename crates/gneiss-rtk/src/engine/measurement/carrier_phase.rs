use crate::engine::measurement::types::{SingleUpdate, UpdateGeometry, VarianceWeights, DdCarrierPhaseParams};
use crate::engine::measurement::DdMeasurementContext;

#[allow(clippy::too_many_arguments)]
fn compute_carrier_phase_update(
    cp: [f64; 4],
    f: [f64; 2],
    idx: [usize; 2],
    ambiguities: &[f64],
    iono_dd: f64,
    geom: &UpdateGeometry,
    var: &VarianceWeights,
    freq_idx: u8,
    iono_idx: Option<[usize; 2]>,  // [sat_iono_idx, ref_iono_idx]
) -> SingleUpdate {
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let lam_sat = c / f[0];
    let lam_ref = c / f[1];
    let cp_dd = (cp[0] * lam_sat - cp[1] * lam_ref) - (cp[2] * lam_sat - cp[3] * lam_ref);
    let n_dd = ambiguities[idx[0]] - ambiguities[idx[1]];

    // Ionosphere state contribution to predicted DD
    let iono_scale = if freq_idx == 1 { 1.0 } else { (f[0] / f[1]).powi(2) };
    let iono_state_dd = if let Some([si, ri]) = iono_idx {
        iono_scale * (ambiguities[si] - ambiguities[ri])
    } else {
        0.0
    };

    let mut h_cp = vec![0.0; geom.state_size];
    h_cp[0] = geom.h_r.x;
    h_cp[1] = geom.h_r.y;
    h_cp[2] = geom.h_r.z;
    h_cp[6] = geom.h_att.x;
    h_cp[7] = geom.h_att.y;
    h_cp[8] = geom.h_att.z;
    h_cp[crate::filter::CORE_STATE_SIZE + idx[0]] = 1.0;
    h_cp[crate::filter::CORE_STATE_SIZE + idx[1]] = -1.0;
    // Ionosphere state Jacobian: ∂z/∂iono_sat = +scale, ∂z/∂iono_ref = -scale
    if let Some([si, ri]) = iono_idx {
        h_cp[crate::filter::CORE_STATE_SIZE + si] = iono_scale;
        h_cp[crate::filter::CORE_STATE_SIZE + ri] = -iono_scale;
    }
    if geom.state_size > 20 {
        h_cp[20] = geom.h_zwd;
    }

    SingleUpdate {
        // Innovation includes predicted iono state contribution
        z: cp_dd - (geom.comp_dd - iono_dd - iono_state_dd + n_dd),
        h: h_cp,
        r: var.val,
        type_code: freq_idx,
        r_ref: var.ref_val,
    }
}

pub fn compute_dd_carrier_phase(
    mctx: &DdMeasurementContext,
    p: &DdCarrierPhaseParams,
) -> Vec<SingleUpdate> {
    let mut updates = Vec::new();
    let var = VarianceWeights {
        val: if p.is_fixed {
            1e-6 * mctx.comps.var_factor
        } else {
            p.cp_base_var * mctx.comps.var_factor
        },
        ref_val: if p.is_fixed {
            1e-6 * mctx.comps.ref_var_factor
        } else {
            p.cp_base_var * mctx.comps.ref_var_factor
        },
    };

    let iono_pair = match (p.iono_idx_sat, p.iono_idx_ref) {
        (Some(si), Some(ri)) => Some([si, ri]),
        _ => None,
    };

    if let (Some(sat_idx), Some(ref_idx)) = (p.sat_idx_l1, p.ref_idx_l1) {
        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = [
            mctx.ctx.rov_ref.cp_l1,
            mctx.ctx.rov_sat.cp_l1,
            mctx.ctx.ref_base.cp_l1,
            mctx.ctx.base_sat.cp_l1,
        ] {
            updates.push(compute_carrier_phase_update(
                [rs1, rr1, bs1, br1],
                [mctx.ctx.sat_state.f1, mctx.ctx.ref_state.f1],
                [sat_idx, ref_idx],
                p.ambiguities,
                mctx.comps.iono_dd_l1,
                mctx.geom,
                &var,
                1,
                iono_pair,
            ));
        }
    }

    if let (Some(sat_idx), Some(ref_idx)) = (p.sat_idx_l2, p.ref_idx_l2) {
        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [
            mctx.ctx.rov_ref.cp_l2,
            mctx.ctx.rov_sat.cp_l2,
            mctx.ctx.ref_base.cp_l2,
            mctx.ctx.base_sat.cp_l2,
        ] {
            updates.push(compute_carrier_phase_update(
                [rs2, rr2, bs2, br2],
                [mctx.ctx.sat_state.f2, mctx.ctx.ref_state.f2],
                [sat_idx, ref_idx],
                p.ambiguities,
                mctx.comps.iono_dd_l2,
                mctx.geom,
                &var,
                2,
                iono_pair,
            ));
        }
    }
    updates
}
