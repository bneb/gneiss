use crate::engine::measurement::types::{DdContext, SingleUpdate};
use nalgebra::Vector3;
use crate::engine::measurement_math;

fn compute_doppler_innovation(
    ctx: &DdContext,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    v_ant: Vector3<f64>,
) -> f64 {
    let lam_sat_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / ctx.sat_state.f1;
    let lam_ref_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / ctx.ref_state.f1;

    let e_sat_rov = (ctx.sat_state.rov_pos - pos_apc).normalize();
    let e_ref_rov = (ctx.ref_state.rov_pos - pos_apc).normalize();
    let e_sat_bas = (ctx.sat_state.bas_pos - base_coord_vec).normalize();
    let e_ref_bas = (ctx.ref_state.bas_pos - base_coord_vec).normalize();

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
    observed_dd_rr - predicted_dd_rr
}

#[allow(clippy::too_many_arguments)]
pub fn compute_dd_doppler(
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
    ref_var_factor: f64,
) -> Option<SingleUpdate> {
    let dop_valid = [
        ctx.rov_sat.doppler,
        ctx.rov_ref.doppler,
        ctx.base_sat.doppler,
        ctx.ref_base.doppler,
    ]
    .iter()
    .all(|&x| x != 0.0);
    if dop_valid {
        let v_ant = velocity + r_b_e * omega_b.cross(lever_arm);
        let innov = compute_doppler_innovation(ctx, pos_apc, base_coord_vec, v_ant);

        tracing::trace!("Doppler Innov. var={:.3} innov={:.3}", dop_base_var, innov);

        let mut h_dop = vec![0.0; state_size];
        h_dop[3] = h_r.x;
        h_dop[4] = h_r.y;
        h_dop[5] = h_r.z;

        let h_dop_att = measurement_math::doppler_attitude_jacobian(r_b_e.matrix(), omega_b, lever_arm, &h_r);
        h_dop[6] = h_dop_att.x;
        h_dop[7] = h_dop_att.y;
        h_dop[8] = h_dop_att.z;

        let h_dop_bg = r_b_e.matrix() * lever_arm.cross_matrix();
        let h_dop_bg = h_r.transpose() * h_dop_bg;
        h_dop[12] = h_dop_bg[0];
        h_dop[13] = h_dop_bg[1];
        h_dop[14] = h_dop_bg[2];

        return Some(SingleUpdate {
            z: innov,
            h: h_dop,
            r: dop_base_var * var_factor,
            type_code: 3,
            r_ref: dop_base_var * ref_var_factor,
        });
    }
    None
}
