import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_cp = """pub fn compute_dd_carrier_phase(
    ctx: &DdContext,
    is_fixed: bool,
    ambiguities: &[f64],
    sat_idx_l1: Option<usize>,
    ref_idx_l1: Option<usize>,
    sat_idx_l2: Option<usize>,
    ref_idx_l2: Option<usize>,
    comp_pr_dd: f64,
    iono_dd_l1: f64,
    iono_dd_l2: f64,
    h_r: Vector3<f64>,
    h_att: Vector3<f64>,
    state_size: usize,
    var_factor: f64,
    ref_var_factor: f64,
    cp_base_var: f64,
    h_zwd: f64) -> Vec<(f64, Vec<f64>, f64, u8, f64)> {
    let mut updates = Vec::new();
    let r_val = if is_fixed { 1e-6 * var_factor } else { cp_base_var * var_factor };
    let r_ref_val = if is_fixed { 1e-6 * ref_var_factor } else { cp_base_var * ref_var_factor };
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l1, ref_idx_l1) {
        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = [ctx.rov_ref.cp_l1, ctx.rov_sat.cp_l1, ctx.ref_base.cp_l1, ctx.base_sat.cp_l1] {
            let lam_ref_1 = c / ctx.ref_state.f1;
            let lam_sat_1 = c / ctx.sat_state.f1;
            let cp_dd_l1 = (rs1 * lam_sat_1 - rr1 * lam_ref_1) - (bs1 * lam_sat_1 - br1 * lam_ref_1);
            let n_dd_l1 = ambiguities[sat_idx] - ambiguities[ref_idx];
            
            let mut h_cp1 = vec![0.0; state_size];
            h_cp1[0] = h_r.x; h_cp1[1] = h_r.y; h_cp1[2] = h_r.z;
            h_cp1[6] = h_att.x; h_cp1[7] = h_att.y; h_cp1[8] = h_att.z;
            h_cp1[crate::filter::CORE_STATE_SIZE + sat_idx] = 1.0; 
            h_cp1[crate::filter::CORE_STATE_SIZE + ref_idx] = -1.0;
            
            if state_size > 20 { h_cp1[20] = h_zwd; }
            updates.push((cp_dd_l1 - (comp_pr_dd - iono_dd_l1 + n_dd_l1), h_cp1, r_val, 1, r_ref_val));
        }
    }

    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l2, ref_idx_l2) {
        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [ctx.rov_ref.cp_l2, ctx.rov_sat.cp_l2, ctx.ref_base.cp_l2, ctx.base_sat.cp_l2] {
            let lam_ref_2 = c / ctx.ref_state.f2;
            let lam_sat_2 = c / ctx.sat_state.f2;
            let cp_dd_l2 = (rs2 * lam_sat_2 - rr2 * lam_ref_2) - (bs2 * lam_sat_2 - br2 * lam_ref_2);
            let n_dd_l2 = ambiguities[sat_idx] - ambiguities[ref_idx];
            
            let mut h_cp2 = vec![0.0; state_size];
            h_cp2[0] = h_r.x; h_cp2[1] = h_r.y; h_cp2[2] = h_r.z;
            h_cp2[6] = h_att.x; h_cp2[7] = h_att.y; h_cp2[8] = h_att.z;
            h_cp2[crate::filter::CORE_STATE_SIZE + sat_idx] = 1.0; 
            h_cp2[crate::filter::CORE_STATE_SIZE + ref_idx] = -1.0;
            
            if state_size > 20 { h_cp2[20] = h_zwd; }
            updates.push((cp_dd_l2 - (comp_pr_dd - iono_dd_l2 + n_dd_l2), h_cp2, r_val, 2, r_ref_val));
        }
    }

    updates
}"""

new_cp = """#[allow(clippy::too_many_arguments)]
fn compute_carrier_phase_update(
    rs: f64, rr: f64, bs: f64, br: f64, f_sat: f64, f_ref: f64,
    sat_idx: usize, ref_idx: usize, ambiguities: &[f64],
    comp_pr_dd: f64, iono_dd: f64, h_r: Vector3<f64>, h_att: Vector3<f64>,
    h_zwd: f64, state_size: usize, r_val: f64, r_ref_val: f64, freq_idx: u8
) -> (f64, Vec<f64>, f64, u8, f64) {
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let lam_ref = c / f_ref;
    let lam_sat = c / f_sat;
    let cp_dd = (rs * lam_sat - rr * lam_ref) - (bs * lam_sat - br * lam_ref);
    let n_dd = ambiguities[sat_idx] - ambiguities[ref_idx];
    
    let mut h_cp = vec![0.0; state_size];
    h_cp[0] = h_r.x; h_cp[1] = h_r.y; h_cp[2] = h_r.z;
    h_cp[6] = h_att.x; h_cp[7] = h_att.y; h_cp[8] = h_att.z;
    h_cp[crate::filter::CORE_STATE_SIZE + sat_idx] = 1.0; 
    h_cp[crate::filter::CORE_STATE_SIZE + ref_idx] = -1.0;
    if state_size > 20 { h_cp[20] = h_zwd; }
    
    (cp_dd - (comp_pr_dd - iono_dd + n_dd), h_cp, r_val, freq_idx, r_ref_val)
}

#[allow(clippy::too_many_arguments)]
pub fn compute_dd_carrier_phase(
    ctx: &DdContext, is_fixed: bool, ambiguities: &[f64],
    sat_idx_l1: Option<usize>, ref_idx_l1: Option<usize>,
    sat_idx_l2: Option<usize>, ref_idx_l2: Option<usize>,
    comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, state_size: usize,
    var_factor: f64, ref_var_factor: f64, cp_base_var: f64, h_zwd: f64
) -> Vec<(f64, Vec<f64>, f64, u8, f64)> {
    let mut updates = Vec::new();
    let r_val = if is_fixed { 1e-6 * var_factor } else { cp_base_var * var_factor };
    let r_ref_val = if is_fixed { 1e-6 * ref_var_factor } else { cp_base_var * ref_var_factor };

    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l1, ref_idx_l1) {
        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = [ctx.rov_ref.cp_l1, ctx.rov_sat.cp_l1, ctx.ref_base.cp_l1, ctx.base_sat.cp_l1] {
            updates.push(compute_carrier_phase_update(
                rs1, rr1, bs1, br1, ctx.sat_state.f1, ctx.ref_state.f1, sat_idx, ref_idx, ambiguities,
                comp_pr_dd, iono_dd_l1, h_r, h_att, h_zwd, state_size, r_val, r_ref_val, 1
            ));
        }
    }

    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l2, ref_idx_l2) {
        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [ctx.rov_ref.cp_l2, ctx.rov_sat.cp_l2, ctx.ref_base.cp_l2, ctx.base_sat.cp_l2] {
            updates.push(compute_carrier_phase_update(
                rs2, rr2, bs2, br2, ctx.sat_state.f2, ctx.ref_state.f2, sat_idx, ref_idx, ambiguities,
                comp_pr_dd, iono_dd_l2, h_r, h_att, h_zwd, state_size, r_val, r_ref_val, 2
            ));
        }
    }
    updates
}"""

content = content.replace(old_cp, new_cp)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
