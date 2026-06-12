use nalgebra::{DMatrix, DVector, Vector3};
use gneiss_core::coords::Coordinate;
use gneiss_core::time::GpsTime;
use gneiss_core::ephemeris::Ephemeris;
use crate::filter::{RtkState, DdObservation};

pub struct MeasurementEnvironment<'a> {
    pub ephemerides: &'a [Ephemeris],
    pub base_coord: &'a Coordinate,
    pub base_time: GpsTime,
    pub lever_arm: Vector3<f64>,
    pub omega_b: Vector3<f64>,
    pub tuning: &'a crate::engine::config::EkfTuningConfig,
}

pub struct SatState {
    pub rov_pos: Vector3<f64>,
    pub rov_vel: Vector3<f64>,
    pub bas_pos: Vector3<f64>,
    pub bas_vel: Vector3<f64>,
    pub f1: f64,
    pub f2: f64,
}

pub struct DdContext<'a> {
    pub rov_sat: &'a mut DdObservation,
    pub base_sat: &'a mut DdObservation,
    pub rov_ref: &'a mut DdObservation,
    pub ref_base: &'a mut DdObservation,
    pub sat_state: &'a SatState,
    pub ref_state: &'a SatState,
}


pub fn get_sat_state(eph: &Ephemeris, pr: f64, t_rx: GpsTime, rx_pos: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let tau_pr = pr / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let t_tx_nom = GpsTime::new(t_rx.week, t_rx.tow - tau_pr);
    let (_, _, dt_s, _) = eph.position(t_tx_nom);
    let t_tx_true = GpsTime::new(t_rx.week, t_rx.tow - tau_pr - dt_s);
    let (raw_vec, raw_vel, _, _) = eph.position(t_tx_true);
    
    let mut sat_pos = raw_vec;
    let mut sat_vel = raw_vel;
    for _ in 0..2 {
        let geometric_range = (sat_pos - rx_pos).norm();
        let true_tau = geometric_range / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let theta = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * true_tau;
        let cos_t = f64::cos(theta);
        let sin_t = f64::sin(theta);
        sat_pos = nalgebra::Vector3::new(
            raw_vec.x * cos_t + raw_vec.y * sin_t,
            -raw_vec.x * sin_t + raw_vec.y * cos_t,
            raw_vec.z
        );
        sat_vel = nalgebra::Vector3::new(
            raw_vel.x * cos_t + raw_vel.y * sin_t,
            -raw_vel.x * sin_t + raw_vel.y * cos_t,
            raw_vel.z
        );
    }
    (sat_pos, sat_vel)
}

#[allow(clippy::too_many_arguments)]
pub fn compute_atmospheric_delays(
    state_time: GpsTime,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    sat_vec_rov: Vector3<f64>,
    ref_sat_vec_rov: Vector3<f64>,
    sat_vec_bas: Vector3<f64>,
    ref_sat_vec_bas: Vector3<f64>,
    sat_f1: f64, sat_f2: f64,
    ref_f1: f64, ref_f2: f64,
) -> (f64, f64, f64) {
    let tropo_params = gneiss_core::atmosphere::TropoParams::default();
    // Ideally we should pass klobuchar from EngineConfig, but we don't have access to Engine here.
    // For RTK, ionosphere cancels out over short baselines anyway.
    let iono_params = gneiss_core::atmosphere::KlobucharParams::default();
    
    let base_llh = gneiss_core::coords::ecef_to_llh(base_coord_vec);
    let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
    
    let (az_rov_sat, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_vec_rov);
    let (az_rov_ref, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, pos_apc, ref_sat_vec_rov);
    let (az_bas_sat, el_bas_sat) = gneiss_core::coords::az_el(base_llh, base_coord_vec, sat_vec_bas);
    let (az_bas_ref, el_bas_ref) = gneiss_core::coords::az_el(base_llh, base_coord_vec, ref_sat_vec_bas);
    
    let tropo_rov_sat = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, rov_llh, el_rov_sat);
    let tropo_rov_ref = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, rov_llh, el_rov_ref);
    let tropo_bas_sat = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, base_llh, el_bas_sat);
    let tropo_bas_ref = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, base_llh, el_bas_ref);
    let tropo_dd = (tropo_rov_sat - tropo_rov_ref) - (tropo_bas_sat - tropo_bas_ref);

    let iono_rov_sat = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, rov_llh, az_rov_sat, el_rov_sat, state_time);
    let iono_rov_ref = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, rov_llh, az_rov_ref, el_rov_ref, state_time);
    let iono_bas_sat = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, base_llh, az_bas_sat, el_bas_sat, state_time);
    let iono_bas_ref = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, base_llh, az_bas_ref, el_bas_ref, state_time);
    let iono_dd_l1 = (iono_rov_sat - iono_rov_ref) - (iono_bas_sat - iono_bas_ref);

    let f_ratio_sat_l2 = (sat_f1 / sat_f2).powi(2);
    let f_ratio_ref_l2 = (ref_f1 / ref_f2).powi(2);
    let iono_dd_l2 = (iono_rov_sat * f_ratio_sat_l2 - iono_rov_ref * f_ratio_ref_l2) - (iono_bas_sat * f_ratio_sat_l2 - iono_bas_ref * f_ratio_ref_l2);

    (tropo_dd, iono_dd_l1, iono_dd_l2)
}

pub fn apply_phase_windup(
    state: &mut RtkState,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    ctx: &mut DdContext,
) {
    let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);

    let prev_w_sat = *state.windup.get(&ctx.rov_sat.sat).unwrap_or(&0.0);
    let prev_w_ref = *state.windup.get(&ctx.rov_ref.sat).unwrap_or(&0.0);
    let w_sat = gneiss_core::windup::phase_windup(ctx.sat_state.rov_pos, sun_pos, pos_apc, prev_w_sat);
    let w_ref = gneiss_core::windup::phase_windup(ctx.ref_state.rov_pos, sun_pos, pos_apc, prev_w_ref);
    state.windup.insert(ctx.rov_sat.sat, w_sat);
    state.windup.insert(ctx.rov_ref.sat, w_ref);

    if let Some(cp) = &mut ctx.rov_sat.cp_l1 { *cp += w_sat; }
    if let Some(cp2) = &mut ctx.rov_sat.cp_l2 { *cp2 += w_sat; }
    if let Some(cp) = &mut ctx.rov_ref.cp_l1 { *cp += w_ref; }
    if let Some(cp2) = &mut ctx.rov_ref.cp_l2 { *cp2 += w_ref; }

    let prev_w_bas_sat = *state.windup.get(&ctx.base_sat.sat).unwrap_or(&0.0);
    let prev_w_bas_ref = *state.windup.get(&ctx.ref_base.sat).unwrap_or(&0.0);
    let w_bas_sat = gneiss_core::windup::phase_windup(ctx.sat_state.bas_pos, sun_pos, base_coord_vec, prev_w_bas_sat);
    let w_bas_ref = gneiss_core::windup::phase_windup(ctx.ref_state.bas_pos, sun_pos, base_coord_vec, prev_w_bas_ref);
    state.windup.insert(ctx.base_sat.sat, w_bas_sat);
    state.windup.insert(ctx.ref_base.sat, w_bas_ref);

    if let Some(cp) = &mut ctx.base_sat.cp_l1 { *cp += w_bas_sat; }
    if let Some(cp2) = &mut ctx.base_sat.cp_l2 { *cp2 += w_bas_sat; }
    if let Some(cp) = &mut ctx.ref_base.cp_l1 { *cp += w_bas_ref; }
    if let Some(cp2) = &mut ctx.ref_base.cp_l2 { *cp2 += w_bas_ref; }
}

#[allow(clippy::too_many_arguments)]
pub fn compute_dd_pseudorange(
    ctx: &DdContext,
    comp_pr_dd: f64,
    iono_dd_l1: f64,
    iono_dd_l2: f64,
    h_r: Vector3<f64>,
    h_att: Vector3<f64>,
    state_size: usize,
    var_factor: f64,
    env: &MeasurementEnvironment,
) -> Vec<(f64, Vec<f64>, f64, u8)> {
    let mut updates = Vec::new();
    let pr_base_var = env.tuning.pr_base_var;

    let pr1_valid = [ctx.rov_sat.pr_l1, ctx.base_sat.pr_l1, ctx.rov_ref.pr_l1, ctx.ref_base.pr_l1].iter().all(|&x| x > 0.0);
    if pr1_valid {
        let pr_dd = (ctx.rov_sat.pr_l1 - ctx.rov_ref.pr_l1) - (ctx.base_sat.pr_l1 - ctx.ref_base.pr_l1);
        let mut h_pr1 = vec![0.0; state_size];
        h_pr1[0] = h_r.x; h_pr1[1] = h_r.y; h_pr1[2] = h_r.z;
        h_pr1[6] = h_att.x; h_pr1[7] = h_att.y; h_pr1[8] = h_att.z;
        updates.push((pr_dd - (comp_pr_dd + iono_dd_l1), h_pr1, pr_base_var * var_factor, 0));
    }

    let pr2_vals = [ctx.rov_ref.pr_l2, ctx.rov_sat.pr_l2, ctx.ref_base.pr_l2, ctx.base_sat.pr_l2];
    if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = pr2_vals {
        let pr_dd_l2 = (rs2 - rr2) - (bs2 - br2);
        let mut h_pr2 = vec![0.0; state_size];
        h_pr2[0] = h_r.x; h_pr2[1] = h_r.y; h_pr2[2] = h_r.z;
        h_pr2[6] = h_att.x; h_pr2[7] = h_att.y; h_pr2[8] = h_att.z;
        updates.push((pr_dd_l2 - (comp_pr_dd + iono_dd_l2), h_pr2, pr_base_var * var_factor, 0));
    }

    updates
}

#[allow(clippy::too_many_arguments)]
pub fn compute_dd_carrier_phase(
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
) -> Vec<(f64, Vec<f64>, f64, u8)> {
    let mut updates = Vec::new();
    let cp_base_var = env.tuning.cp_base_var;
    let r_val = if state.is_fixed { 1e-6 * var_factor } else { cp_base_var * var_factor };
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

    if let (Some(sat_idx), Some(ref_idx)) = (state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 1), ref_idx_l1) {
        let cp1_vals = [ctx.rov_ref.cp_l1, ctx.rov_sat.cp_l1, ctx.ref_base.cp_l1, ctx.base_sat.cp_l1];
        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = cp1_vals {
            let lam_ref_1 = c / ctx.ref_state.f1;
            let lam_sat_1 = c / ctx.sat_state.f1;
            let cp_dd_l1 = (rs1 * lam_sat_1 - rr1 * lam_ref_1) - (bs1 * lam_sat_1 - br1 * lam_ref_1);
            let n_dd_l1 = state.ambiguities[sat_idx] - state.ambiguities[ref_idx];
            
            let mut h_cp1 = vec![0.0; state_size];
            h_cp1[0] = h_r.x; h_cp1[1] = h_r.y; h_cp1[2] = h_r.z;
            h_cp1[6] = h_att.x; h_cp1[7] = h_att.y; h_cp1[8] = h_att.z;
            h_cp1[crate::filter::CORE_STATE_SIZE + sat_idx] = 1.0; 
            h_cp1[crate::filter::CORE_STATE_SIZE + ref_idx] = -1.0;
            
            tracing::debug!("CP1 DD computed for {:?}", ctx.rov_sat.sat);
            updates.push((cp_dd_l1 - (comp_pr_dd - iono_dd_l1 + n_dd_l1), h_cp1, r_val, 1));
        } else {
            tracing::debug!("Missing CP1 obs for {:?}! rr={:?} rs={:?} br={:?} bs={:?}", ctx.rov_sat.sat, ctx.rov_ref.cp_l1, ctx.rov_sat.cp_l1, ctx.ref_base.cp_l1, ctx.base_sat.cp_l1);
        }
    } else {
        tracing::debug!("Missing ambiguity for {:?}! sat_idx={:?} ref_idx={:?}", ctx.rov_sat.sat, state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 1), ref_idx_l1);
    }

    if let (Some(sat_idx), Some(ref_idx)) = (state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 2), ref_idx_l2) {
        let cp2_vals = [ctx.rov_ref.cp_l2, ctx.rov_sat.cp_l2, ctx.ref_base.cp_l2, ctx.base_sat.cp_l2];
        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = cp2_vals {
            let lam_ref_2 = c / ctx.ref_state.f2;
            let lam_sat_2 = c / ctx.sat_state.f2;
            let cp_dd_l2 = (rs2 * lam_sat_2 - rr2 * lam_ref_2) - (bs2 * lam_sat_2 - br2 * lam_ref_2);
            let n_dd_l2 = state.ambiguities[sat_idx] - state.ambiguities[ref_idx];
            
            let mut h_cp2 = vec![0.0; state_size];
            h_cp2[0] = h_r.x; h_cp2[1] = h_r.y; h_cp2[2] = h_r.z;
            h_cp2[6] = h_att.x; h_cp2[7] = h_att.y; h_cp2[8] = h_att.z;
            h_cp2[crate::filter::CORE_STATE_SIZE + sat_idx] = 1.0; 
            h_cp2[crate::filter::CORE_STATE_SIZE + ref_idx] = -1.0;
            
            updates.push((cp_dd_l2 - (comp_pr_dd - iono_dd_l2 + n_dd_l2), h_cp2, r_val, 2));
        }
    }

    updates
}

#[allow(clippy::too_many_arguments)]
pub fn compute_dd_doppler(
    state: &RtkState,
    ctx: &DdContext,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    h_r: Vector3<f64>,
    state_size: usize,
    var_factor: f64,
    env: &MeasurementEnvironment,
) -> Option<(f64, Vec<f64>, f64, u8)> {
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
        
        let a = r_b_e * env.omega_b.cross(&env.lever_arm);
        let h_dop_att = a.cross(&h_r);
        h_dop[6] = h_dop_att.x; h_dop[7] = h_dop_att.y; h_dop[8] = h_dop_att.z;
        
        let h_dop_bg = r_b_e.matrix() * env.lever_arm.cross_matrix();
        let h_dop_bg = h_r.transpose() * h_dop_bg;
        h_dop[12] = h_dop_bg[0]; h_dop[13] = h_dop_bg[1]; h_dop[14] = h_dop_bg[2];
        
        let dop_base_var = env.tuning.dop_base_var;
        return Some((innov, h_dop, dop_base_var * var_factor, 3));
    }
    None
}

#[allow(clippy::type_complexity)]
pub fn compute_innovations(
    state: &mut RtkState,
    group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation,
    ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<(Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<(gneiss_core::sat::SatelliteId, u8)>)> {
    let mut z_vals = Vec::new();
    let mut h_rows = Vec::new();
    let mut r_vals = Vec::new();
    let mut meas_type = Vec::new();

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    
    let r_b_e = state.attitude.to_rotation_matrix();
    let pos_apc = state.position.vector + r_b_e * env.lever_arm;

    let set_rov = gneiss_core::tides::solid_earth_tides_ecef(state.time, pos_apc);
    let set_bas = gneiss_core::tides::solid_earth_tides_ecef(state.time, env.base_coord.vector);
    
    let pos_apc = pos_apc + set_rov;
    let base_coord_vec = env.base_coord.vector + set_bas;

    let time_tow = state.time.tow;
    let find_eph = |sat| {
        env.ephemerides.iter()
            .filter(|e| e.sat() == sat)
            .min_by(|a, b| {
                let da = (a.toe().tow - time_tow).abs();
                let db = (b.toe().tow - time_tow).abs();
                da.partial_cmp(&db).unwrap()
            })
    };

    let ref_eph = find_eph(ref_rover_orig.sat)?;
    let (ref_sat_vec_rov, ref_sat_vel_rov) = get_sat_state(ref_eph, ref_rover_orig.pr_l1, state.time, pos_apc);
    let (ref_sat_vec_bas, ref_sat_vel_bas) = get_sat_state(ref_eph, ref_base_orig.pr_l1, env.base_time, base_coord_vec);
    let (ref_f1, ref_f2) = gneiss_core::signal::satellite_frequencies(ref_rover_orig.sat, ref_eph.freq_num());

    let ref_state = SatState {
        rov_pos: ref_sat_vec_rov,
        rov_vel: ref_sat_vel_rov,
        bas_pos: ref_sat_vec_bas,
        bas_vel: ref_sat_vel_bas,
        f1: ref_f1,
        f2: ref_f2,
    };

    let e_ref_rov = (ref_sat_vec_rov - pos_apc).normalize();
    
    let (ref_idx_l1, ref_idx_l2) = {
        let idx1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 1)?;
        let idx2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 2);
        (idx1, idx2)
    };

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_eph(rover_sat_orig.sat) {
            
            let (sat_vec_rov, sat_vel_rov) = get_sat_state(sat_eph, rover_sat_orig.pr_l1, state.time, pos_apc);
            let (sat_vec_bas, sat_vel_bas) = get_sat_state(sat_eph, base_sat_orig.pr_l1, env.base_time, base_coord_vec);
            let (sat_f1, sat_f2) = gneiss_core::signal::satellite_frequencies(rover_sat_orig.sat, sat_eph.freq_num());

            let sat_state = SatState {
                rov_pos: sat_vec_rov,
                rov_vel: sat_vel_rov,
                bas_pos: sat_vec_bas,
                bas_vel: sat_vel_bas,
                f1: sat_f1,
                f2: sat_f2,
            };

            let e_sat_rov = (sat_vec_rov - pos_apc).normalize();
            let h_r = e_ref_rov - e_sat_rov; 
            
            let lever_ecef = r_b_e * env.lever_arm;
            let h_att = lever_ecef.cross(&h_r);

            let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(
                state.time, pos_apc, base_coord_vec, sat_vec_rov, ref_sat_vec_rov, sat_vec_bas, ref_sat_vec_bas,
                sat_f1, sat_f2, ref_f1, ref_f2
            );

            let mut rov_sat = rover_sat_orig.clone();
            let mut bas_sat = base_sat_orig.clone();
            let mut rov_ref = ref_rover_orig.clone();
            let mut bas_ref = ref_base_orig.clone();

            let mut ctx = DdContext {
                rov_sat: &mut rov_sat,
                base_sat: &mut bas_sat,
                rov_ref: &mut rov_ref,
                ref_base: &mut bas_ref,
                sat_state: &sat_state,
                ref_state: &ref_state,
            };

            apply_phase_windup(
                state, pos_apc, base_coord_vec, &mut ctx
            );

            let base_llh = gneiss_core::coords::ecef_to_llh(base_coord_vec);
            let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
            let (_, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_vec_rov);
            let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, pos_apc, ref_sat_vec_rov);
            let (_, el_bas_sat) = gneiss_core::coords::az_el(base_llh, base_coord_vec, sat_vec_bas);
            let (_, el_bas_ref) = gneiss_core::coords::az_el(base_llh, base_coord_vec, ref_sat_vec_bas);

            let var_factor = gneiss_core::variance::observation_variance(ctx.rov_sat.snr, el_rov_sat, 45.0)
                           + gneiss_core::variance::observation_variance(ctx.rov_ref.snr, el_rov_ref, 45.0)
                           + gneiss_core::variance::elevation_variance_scale(el_bas_sat)
                           + gneiss_core::variance::elevation_variance_scale(el_bas_ref);

            let comp_pr_dd = ( (pos_apc - sat_vec_rov).norm() - (pos_apc - ref_sat_vec_rov).norm() ) - ( (base_coord_vec - sat_vec_bas).norm() - (base_coord_vec - ref_sat_vec_bas).norm() ) + tropo_dd;
            
            for update in compute_dd_pseudorange(&ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, state_size, var_factor, env) {
                z_vals.push(update.0); h_rows.push(update.1); r_vals.push(update.2); meas_type.push((ctx.rov_sat.sat, update.3));
            }

            for update in compute_dd_carrier_phase(state, &ctx, Some(ref_idx_l1), ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, state_size, var_factor, env) {
                z_vals.push(update.0); h_rows.push(update.1); r_vals.push(update.2); meas_type.push((ctx.rov_sat.sat, update.3));
            }
            
            if let Some(update) = compute_dd_doppler(state, &ctx, pos_apc, base_coord_vec, h_r, state_size, var_factor, env) {
                z_vals.push(update.0); h_rows.push(update.1); r_vals.push(update.2); meas_type.push((ctx.rov_sat.sat, update.3));
            }
        }
    }
    
    Some((z_vals, h_rows, r_vals, meas_type))
}

#[allow(clippy::type_complexity)]
pub fn build_measurement_model(
    state: &mut RtkState,
    matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8)>)> {
    let mut z_all = Vec::new();
    let mut h_all = Vec::new();
    let mut r_all = Vec::new();
    let mut type_all = Vec::new();

    use std::collections::HashMap;
    use gneiss_core::sat::Constellation;
    let mut const_groups: HashMap<Constellation, Vec<(DdObservation, DdObservation)>> = HashMap::new();
    
    for obs in matched_obs {
        const_groups.entry(obs.0.sat.constellation).or_default().push(obs.clone());
    }

    for (_, group) in const_groups {
        if group.len() < 2 { continue; } 

        let mut best_score = -1.0;
        let mut ref_idx = 0;
        let rov_llh_ref = gneiss_core::coords::ecef_to_llh(state.position.vector);
        let time_tow = state.time.tow;
        for (i, (r, _)) in group.iter().enumerate() {
            if let Some(eph) = env.ephemerides.iter().filter(|e| e.sat() == r.sat).min_by(|a, b| {
                let da = (a.toe().tow - time_tow).abs();
                let db = (b.toe().tow - time_tow).abs();
                da.partial_cmp(&db).unwrap()
            }) {
                let tau = r.pr_l1 / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
                let t_tx = gneiss_core::time::GpsTime::new(state.time.week, state.time.tow - tau);
                let (sat_pos, _, _, _) = eph.position(t_tx);
                let (_, el) = gneiss_core::coords::az_el(rov_llh_ref, state.position.vector, sat_pos);
                
                let score = if r.cp_l1.is_some() { el + 100.0 } else { el };
                
                if score > best_score {
                    best_score = score;
                    ref_idx = i;
                }
            }
        }
        
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some((z, h, r, mt)) = compute_innovations(state, &group_clone, &ref_rover, &ref_base, env) {
            z_all.extend(z);
            h_all.extend(h);
            r_all.extend(r);
            type_all.extend(mt);
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    let mut safe_indices = Vec::new();
    
    tracing::trace!("Pre-filter z_all len: {}", z_all.len());
    
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size { h_row[(0, c)] = h_all[i][c]; }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;
        
        let threshold = match type_all[i].1 { 
            0 => chi_square_pr_threshold * chi_square_pr_threshold,  
            1 | 2 => chi_square_cp_threshold * chi_square_cp_threshold,  
            3 => chi_square_pr_threshold * 1000.0, // Doppler relax
            _ => chi_square_pr_threshold * chi_square_pr_threshold   
        };
        
        if chi2 <= threshold { 
            safe_indices.push(i); 
            if z_all[i].abs() > 100000.0 {
                tracing::error!("MASSIVE Z PASSED PRE-FILTER! type: {}, z: {:.1}, chi2: {:.1}, thresh: {:.1}, s_ii: {:.1}, P_pos: {:.1}", type_all[i].1, z_all[i], chi2, threshold, s_ii, state.covariance[(0,0)]);
            }
            if type_all[i].1 == 1 || type_all[i].1 == 2 {
                for c in crate::filter::CORE_STATE_SIZE..state_size {
                    if h_row[(0, c)] > 0.5 {
                        let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
                        state.reject_counts.insert(key, 0);
                    }
                }
            }
        } else {
            tracing::debug!("Rejected meas type {} with inn: {:.3}, chi2: {:.1}, threshold: {:.1}, s_ii: {:.1}", type_all[i].1, z_all[i], chi2, threshold, s_ii);
            if type_all[i].1 == 1 || type_all[i].1 == 2 {
                for c in crate::filter::CORE_STATE_SIZE..state_size {
                    if h_row[(0, c)] > 0.5 {
                        let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
                        let count = *state.reject_counts.get(&key).unwrap_or(&0) + 1;
                        state.reject_counts.insert(key, count);
                    }
                }
            }
        }
    }

    if safe_indices.len() >= 4 {
        let mut z_vec = DVector::zeros(safe_indices.len());
        let mut h_mat = DMatrix::zeros(safe_indices.len(), state_size);
        let mut r_mat = DMatrix::zeros(safe_indices.len(), safe_indices.len());
        let mut t_vec = Vec::new();

        for (new_i, &old_i) in safe_indices.iter().enumerate() {
            z_vec[new_i] = z_all[old_i];
            for c in 0..state_size { h_mat[(new_i, c)] = h_all[old_i][c]; }
            r_mat[(new_i, new_i)] = r_all[old_i];
            t_vec.push(type_all[old_i]);
        }
        Some((z_vec, h_mat, r_mat, t_vec))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    
    use crate::filter::{RtkState, DdObservation};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use gneiss_core::sat::{SatelliteId, Constellation};
    use gneiss_core::ephemeris::Ephemeris;
    use nalgebra::Vector3;

    #[test]
    fn test_measurement_model_against_rtklib_golden_data() {
        let time = GpsTime::new(2137, 422922.0);
        let mut state = RtkState::new(time, Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time), 10.0);
        state.velocity = Vector3::new(10.0, -5.0, 2.0);
        
        let ref_sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let rov_sat1 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let rov_sat2 = SatelliteId { constellation: Constellation::Gps, prn: 3 };
        
        state.add_ambiguity(ref_sat, 1, 5.0, 100.0);
        state.add_ambiguity(rov_sat1, 1, 10.0, 100.0);
        state.add_ambiguity(rov_sat2, 1, 15.0, 100.0);
        
        state.windup.insert(ref_sat, 0.0);
        state.windup.insert(rov_sat1, 0.0);
        state.windup.insert(rov_sat2, 0.0);

        let ref_rover = DdObservation { sat: ref_sat, pr_l1: 20000000.0, pr_l2: Some(20000001.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 100.0, snr: 45.0, locktime: Some(100) };
        let ref_base = DdObservation { sat: ref_sat, pr_l1: 20005000.0, pr_l2: Some(20005001.0), cp_l1: Some(100020000.0), cp_l2: Some(80016000.0), doppler: 10.0, snr: 45.0, locktime: Some(100) };
        
        let rov1_rover = DdObservation { sat: rov_sat1, pr_l1: 21000000.0, pr_l2: Some(21000001.0), cp_l1: Some(105000000.0), cp_l2: Some(84000000.0), doppler: -50.0, snr: 45.0, locktime: Some(100) };
        let rov1_base = DdObservation { sat: rov_sat1, pr_l1: 21005000.0, pr_l2: Some(21005001.0), cp_l1: Some(105020000.0), cp_l2: Some(84016000.0), doppler: 10.0, snr: 45.0, locktime: Some(100) };

        let rov2_rover = DdObservation { sat: rov_sat2, pr_l1: 22000000.0, pr_l2: Some(22000001.0), cp_l1: Some(110000000.0), cp_l2: Some(88000000.0), doppler: -20.0, snr: 45.0, locktime: Some(100) };
        let rov2_base = DdObservation { sat: rov_sat2, pr_l1: 22005000.0, pr_l2: Some(22005001.0), cp_l1: Some(110020000.0), cp_l2: Some(88016000.0), doppler: 10.0, snr: 45.0, locktime: Some(100) };

        let matched_obs = vec![(rov1_rover, rov1_base), (rov2_rover, rov2_base)];

        let eph_ref = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: ref_sat, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });

        let eph_rov1 = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: rov_sat1, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 1.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.5, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });
        
        let eph_rov2 = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: rov_sat2, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 2.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 1.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });

        let ephemerides = vec![eph_ref, eph_rov1, eph_rov2];
        let base_coord = Coordinate::new(Vector3::new(1005.0, 2005.0, 3005.0), Datum::WGS84, Frame::ECEF, time);

        // We explicitly use compute_innovations to avoid the Mahalanobis chi2 filter rejecting dummy data
        let config = crate::engine::EngineConfig::default();
        let env = super::super::measurement::MeasurementEnvironment {
            ephemerides: &ephemerides,
            base_coord: &base_coord,
            base_time: base_coord.epoch,
            lever_arm: Vector3::zeros(),
            omega_b: Vector3::zeros(),
            tuning: &config.tuning,
        };
        let (z, _h, r, _) = super::super::measurement::compute_innovations(&mut state, &matched_obs, &ref_rover, &ref_base, &env).unwrap();

        println!("Z: {:?}", z);
        
        // Lock in the golden Z vector (updated for iterative ecef_to_llh refinement)
        assert!((z[0] - 0.70573).abs() < 1e-3, "z[0]={}", z[0]);
        assert!((z[1] - 0.70623).abs() < 1e-3, "z[1]={}", z[1]);
        assert!((z[2] - -4.29582).abs() < 1e-3, "z[2]={}", z[2]);
        assert!((z[3] - 19.36016).abs() < 1e-3, "z[3]={}", z[3]);
        assert!((z[4] - 7.72370).abs() < 1e-3, "z[4]={}", z[4]);
        assert!((z[5] - 7.72444).abs() < 1e-3, "z[5]={}", z[5]);
        assert!((z[6] - -2.27859).abs() < 1e-3, "z[6]={}", z[6]);
        assert!((z[7] - 16.01613).abs() < 1e-3, "z[7]={}", z[7]);

        // Lock in the golden R diagonal
        assert!(r[0] >= 16.0); // Now scales with elevation
        assert!(r[1] >= 16.0);
        assert!(r[2] >= 0.0001);
        assert!(r[3] >= 0.1); // Doppler var
        assert!(r[4] >= 16.0);
        assert!(r[5] >= 16.0);
        assert!(r[6] >= 0.0001);
        assert!(r[7] >= 0.1); // Doppler var
    }

    #[test]
    fn test_get_sat_state() {
        use crate::engine::measurement::get_sat_state;
        let time = GpsTime::new(2137, 422922.0);
        let rx_pos = Vector3::new(1000.0, 2000.0, 3000.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 1.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });

        let pr = 20000000.0;
        let (pos, vel) = get_sat_state(&eph, pr, time, rx_pos);
        
        // Assert non-zero output
        assert!(pos.norm() > 10000.0);
        assert!(vel.norm() > 10.0);
        
        let (pos0, _vel0) = get_sat_state(&eph, 0.0, time, rx_pos);
        assert!((pos.x - pos0.x).abs() > 0.0);
    }
}
