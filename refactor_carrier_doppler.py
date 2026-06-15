import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Refactor compute_dd_carrier_phase
old_cp = """pub fn compute_dd_carrier_phase(
    state: &RtkState,
    ctx: &DdContext,
    ref_idx_l1: Option<usize>,
    ref_idx_l2: Option<usize>,
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
    let cp_base_var = env.tuning.cp_base_var;
    let r_val = if state.is_fixed { 1e-6 * var_factor } else { cp_base_var * var_factor };
    let r_ref_val = if state.is_fixed { 1e-6 * ref_var_factor } else { cp_base_var * ref_var_factor };
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

    let sat_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 1);
    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l1, ref_idx_l1) {
        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = [ctx.rov_ref.cp_l1, ctx.rov_sat.cp_l1, ctx.ref_base.cp_l1, ctx.base_sat.cp_l1] {
            let lam_ref_1 = c / ctx.ref_state.f1;
            let lam_sat_1 = c / ctx.sat_state.f1;
            let cp_dd_l1 = (rs1 * lam_sat_1 - rr1 * lam_ref_1) - (bs1 * lam_sat_1 - br1 * lam_ref_1);
            let n_dd_l1 = state.ambiguities[sat_idx] - state.ambiguities[ref_idx];
            
            let mut h_cp1 = vec![0.0; state_size];
            h_cp1[0] = h_r.x; h_cp1[1] = h_r.y; h_cp1[2] = h_r.z;
            h_cp1[6] = h_att.x; h_cp1[7] = h_att.y; h_cp1[8] = h_att.z;
            h_cp1[crate::filter::CORE_STATE_SIZE + sat_idx] = 1.0; 
            h_cp1[crate::filter::CORE_STATE_SIZE + ref_idx] = -1.0;
            
            if state_size > 20 { h_cp1[20] = h_zwd; }
            updates.push((cp_dd_l1 - (comp_pr_dd - iono_dd_l1 + n_dd_l1), h_cp1, r_val, 1, r_ref_val));
        }
    }

    let sat_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 2);
    if let (Some(sat_idx), Some(ref_idx)) = (sat_idx_l2, ref_idx_l2) {
        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = [ctx.rov_ref.cp_l2, ctx.rov_sat.cp_l2, ctx.ref_base.cp_l2, ctx.base_sat.cp_l2] {
            let lam_ref_2 = c / ctx.ref_state.f2;
            let lam_sat_2 = c / ctx.sat_state.f2;
            let cp_dd_l2 = (rs2 * lam_sat_2 - rr2 * lam_ref_2) - (bs2 * lam_sat_2 - br2 * lam_ref_2);
            let n_dd_l2 = state.ambiguities[sat_idx] - state.ambiguities[ref_idx];
            
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

new_cp = """pub fn compute_dd_carrier_phase(
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

content = content.replace(old_cp, new_cp)


old_dop = """pub fn compute_dd_doppler(
    state: &RtkState,
    ctx: &DdContext,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    h_r: Vector3<f64>,
    state_size: usize,
    var_factor: f64,
    env: &MeasurementEnvironment,
 ref_var_factor: f64) -> Option<(f64, Vec<f64>, f64, u8, f64)> {
    let dop_valid = [ctx.rov_sat.doppler, ctx.rov_ref.doppler, ctx.base_sat.doppler, ctx.ref_base.doppler].iter().all(|&x| x != 0.0);
    if dop_valid {
        let lam_sat_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / ctx.sat_state.f1;
        let lam_ref_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / ctx.ref_state.f1;
        
        let e_sat_rov = (ctx.sat_state.rov_pos - pos_apc).normalize();
        let e_ref_rov = (ctx.ref_state.rov_pos - pos_apc).normalize();
        let e_sat_bas = (ctx.sat_state.bas_pos - base_coord_vec).normalize();
        let e_ref_bas = (ctx.ref_state.bas_pos - base_coord_vec).normalize();
        
        let r_b_e = state.attitude.to_rotation_matrix();
        let v_ant = state.velocity + r_b_e * env.omega_b.cross(&env.lever_arm);
        
        // dr/dt = e_los · (v_sat - v_rcv)
        let rr_rov_sat = e_sat_rov.dot(&(ctx.sat_state.rov_vel - v_ant));
        let rr_rov_ref = e_ref_rov.dot(&(ctx.ref_state.rov_vel - v_ant));
        let rr_bas_sat = e_sat_bas.dot(&(ctx.sat_state.bas_vel));
        let rr_bas_ref = e_ref_bas.dot(&(ctx.ref_state.bas_vel));
        
        let predicted_dd_rr = (rr_rov_sat - rr_rov_ref) - (rr_bas_sat - rr_bas_ref);
        
        let obs_rov_sat = -ctx.rov_sat.doppler * lam_sat_1;
        let obs_rov_ref = -ctx.rov_ref.doppler * lam_ref_1;
        let obs_bas_sat = -ctx.base_sat.doppler * lam_sat_1;
        let obs_bas_ref = -ctx.ref_base.doppler * lam_ref_1;
        
        let observed_dd_rr = (obs_rov_sat - obs_rov_ref) - (obs_bas_sat - obs_bas_ref);
        let innov = observed_dd_rr - predicted_dd_rr;
        
        let dop_base_var = env.tuning.dop_base_var;
        tracing::trace!("Doppler Innov. var={:.3} innov={:.3} obs={:.3} pred={:.3}", dop_base_var, innov, observed_dd_rr, predicted_dd_rr);
        
        let mut h_dop = vec![0.0; state_size]; 
        h_dop[3] = h_r.x; h_dop[4] = h_r.y; h_dop[5] = h_r.z;
        
        let h_dop_att = doppler_attitude_jacobian(r_b_e.matrix(), &env.omega_b, &env.lever_arm, &h_r);
        h_dop[6] = h_dop_att.x; h_dop[7] = h_dop_att.y; h_dop[8] = h_dop_att.z;
        
        let h_dop_bg = r_b_e.matrix() * env.lever_arm.cross_matrix();
        let h_dop_bg = h_r.transpose() * h_dop_bg;
        h_dop[12] = h_dop_bg[0]; h_dop[13] = h_dop_bg[1]; h_dop[14] = h_dop_bg[2];
        
        let dop_base_var = env.tuning.dop_base_var;
        return Some((innov, h_dop, dop_base_var * var_factor, 3, dop_base_var * ref_var_factor));
    }
    None
}"""


new_dop = """pub fn compute_dd_doppler(
    ctx: &DdContext,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    h_r: Vector3<f64>,
    r_b_e: &nalgebra::Rotation3<f64>,
    velocity: &Vector3<f64>,
    omega_b: &Vector3<f64>,
    lever_arm: &Vector3<f64>,
    dop_base_var: f64,
    state_size: usize,
    var_factor: f64,
    ref_var_factor: f64) -> Option<(f64, Vec<f64>, f64, u8, f64)> {
    let dop_valid = [ctx.rov_sat.doppler, ctx.rov_ref.doppler, ctx.base_sat.doppler, ctx.ref_base.doppler].iter().all(|&x| x != 0.0);
    if dop_valid {
        let lam_sat_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / ctx.sat_state.f1;
        let lam_ref_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / ctx.ref_state.f1;
        
        let e_sat_rov = (ctx.sat_state.rov_pos - pos_apc).normalize();
        let e_ref_rov = (ctx.ref_state.rov_pos - pos_apc).normalize();
        let e_sat_bas = (ctx.sat_state.bas_pos - base_coord_vec).normalize();
        let e_ref_bas = (ctx.ref_state.bas_pos - base_coord_vec).normalize();
        
        let v_ant = velocity + r_b_e * omega_b.cross(lever_arm);
        
        // dr/dt = e_los · (v_sat - v_rcv)
        let rr_rov_sat = e_sat_rov.dot(&(ctx.sat_state.rov_vel - v_ant));
        let rr_rov_ref = e_ref_rov.dot(&(ctx.ref_state.rov_vel - v_ant));
        let rr_bas_sat = e_sat_bas.dot(&(ctx.sat_state.bas_vel));
        let rr_bas_ref = e_ref_bas.dot(&(ctx.ref_state.bas_vel));
        
        let predicted_dd_rr = (rr_rov_sat - rr_rov_ref) - (rr_bas_sat - rr_bas_ref);
        
        let obs_rov_sat = -ctx.rov_sat.doppler * lam_sat_1;
        let obs_rov_ref = -ctx.rov_ref.doppler * lam_ref_1;
        let obs_bas_sat = -ctx.base_sat.doppler * lam_sat_1;
        let obs_bas_ref = -ctx.ref_base.doppler * lam_ref_1;
        
        let observed_dd_rr = (obs_rov_sat - obs_rov_ref) - (obs_bas_sat - obs_bas_ref);
        let innov = observed_dd_rr - predicted_dd_rr;
        
        tracing::trace!("Doppler Innov. var={:.3} innov={:.3} obs={:.3} pred={:.3}", dop_base_var, innov, observed_dd_rr, predicted_dd_rr);
        
        let mut h_dop = vec![0.0; state_size]; 
        h_dop[3] = h_r.x; h_dop[4] = h_r.y; h_dop[5] = h_r.z;
        
        let h_dop_att = doppler_attitude_jacobian(r_b_e.matrix(), omega_b, lever_arm, &h_r);
        h_dop[6] = h_dop_att.x; h_dop[7] = h_dop_att.y; h_dop[8] = h_dop_att.z;
        
        let h_dop_bg = r_b_e.matrix() * lever_arm.cross_matrix();
        let h_dop_bg = h_r.transpose() * h_dop_bg;
        h_dop[12] = h_dop_bg[0]; h_dop[13] = h_dop_bg[1]; h_dop[14] = h_dop_bg[2];
        
        return Some((innov, h_dop, dop_base_var * var_factor, 3, dop_base_var * ref_var_factor));
    }
    None
}"""

content = content.replace(old_dop, new_dop)


# Now update calls in generate_measurement_updates
old_gen = """    for update in compute_dd_carrier_phase(
        state, ctx, ref_idx_l1, ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size,
        var_factor, env, h_zwd, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }

    if let Some(update) = compute_dd_doppler(
        state, ctx, geom.pos_apc, geom.base_coord_vec, h_r, geom.state_size,
        var_factor, env, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }"""


new_gen = """    let sat_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 1);
    let sat_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 2);

    for update in compute_dd_carrier_phase(
        ctx, state.is_fixed, &state.ambiguities, sat_idx_l1, ref_idx_l1, sat_idx_l2, ref_idx_l2,
        comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size,
        var_factor, ref_var_factor, env.tuning.cp_base_var, h_zwd
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(update) = compute_dd_doppler(
        ctx, geom.pos_apc, geom.base_coord_vec, h_r, &r_b_e_rot,
        &state.velocity, &env.omega_b, &env.lever_arm, env.tuning.dop_base_var,
        geom.state_size, var_factor, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }"""

content = content.replace(old_gen, new_gen)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
