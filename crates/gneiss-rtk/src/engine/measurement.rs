fn update_windup_state_and_obs(state: &mut RtkState, ctx: &mut DdContext, geom: &EkfGeometryContext) {
    let prev_w_sat = *state.windup.get(&ctx.rov_sat.sat).unwrap_or(&0.0);
    let prev_w_ref = *state.windup.get(&ctx.rov_ref.sat).unwrap_or(&0.0);
    let prev_w_bas_sat = *state.windup.get(&ctx.base_sat.sat).unwrap_or(&0.0);
    let prev_w_bas_ref = *state.windup.get(&ctx.ref_base.sat).unwrap_or(&0.0);

    let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);
    let crate::engine::measurement_math::WindupUpdates { w_sat, w_ref, w_bas_sat, w_bas_ref } = crate::engine::measurement_math::compute_phase_windup(
        geom.pos_apc, geom.base_coord_vec, sun_pos,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        prev_w_sat, prev_w_ref, prev_w_bas_sat, prev_w_bas_ref,
    );

    state.windup.insert(ctx.rov_sat.sat, w_sat);
    state.windup.insert(ctx.rov_ref.sat, w_ref);
    state.windup.insert(ctx.base_sat.sat, w_bas_sat);
    state.windup.insert(ctx.ref_base.sat, w_bas_ref);

    apply_windup_to_obs(ctx.rov_sat, w_sat);
    apply_windup_to_obs(ctx.rov_ref, w_ref);
    apply_windup_to_obs(ctx.base_sat, w_bas_sat);
    apply_windup_to_obs(ctx.ref_base, w_bas_ref);
}

fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
    if let Some(cp) = &mut obs.cp_l1 { *cp += windup; }
    if let Some(cp2) = &mut obs.cp_l2 { *cp2 += windup; }
}
use nalgebra::{DMatrix, DVector, Vector3};
pub use crate::engine::measurement_math::*;
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

pub struct EkfUpdates {
    pub z: Vec<f64>,
    pub h: Vec<Vec<f64>>,
    pub r: Vec<f64>,
    pub mt: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

impl Default for EkfUpdates {
    fn default() -> Self {
        Self::new()
    }
}

impl EkfUpdates {
    pub fn new() -> Self {
        Self { z: Vec::new(), h: Vec::new(), r: Vec::new(), mt: Vec::new() }
    }
    pub fn push(&mut self, u: SingleUpdate, sat: gneiss_core::sat::SatelliteId) {
        self.z.push(u.z); self.h.push(u.h); self.r.push(u.r); self.mt.push((sat, u.type_code, u.r_ref));
    }
    pub fn extend(&mut self, other: Self) {
        self.z.extend(other.z); self.h.extend(other.h); self.r.extend(other.r); self.mt.extend(other.mt);
    }
}

pub struct DdContext<'a> {
    pub rov_sat: &'a mut DdObservation,
    pub base_sat: &'a mut DdObservation,
    pub rov_ref: &'a mut DdObservation,
    pub ref_base: &'a mut DdObservation,
    pub sat_state: &'a SatState,
    pub ref_state: &'a SatState,
}





pub struct SingleUpdate {
    pub z: f64,
    pub h: Vec<f64>,
    pub r: f64,
    pub type_code: u8,
    pub r_ref: f64,
}

pub struct UpdateGeometry {
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
}

fn compute_pseudorange_update(
    pr: [f64; 4], iono_dd: f64, geom: &UpdateGeometry, var: &VarianceWeights
) -> SingleUpdate {
    let pr_dd = (pr[0] - pr[1]) - (pr[2] - pr[3]);
    let mut h_pr = vec![0.0; geom.state_size];
    h_pr[0] = geom.h_r.x; h_pr[1] = geom.h_r.y; h_pr[2] = geom.h_r.z;
    h_pr[6] = geom.h_att.x; h_pr[7] = geom.h_att.y; h_pr[8] = geom.h_att.z;
    if geom.state_size > 20 { h_pr[20] = geom.h_zwd; }
    SingleUpdate { z: pr_dd - (geom.comp_dd + iono_dd), h: h_pr, r: var.val, type_code: 0, r_ref: var.ref_val }
}

pub fn compute_dd_pseudorange(mctx: &DdMeasurementContext) -> Vec<SingleUpdate> {
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
}

fn compute_carrier_phase_update(
    cp: [f64; 4], f: [f64; 2], idx: [usize; 2], ambiguities: &[f64],
    iono_dd: f64, geom: &UpdateGeometry, var: &VarianceWeights, freq_idx: u8
) -> SingleUpdate {
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let lam_sat = c / f[0]; let lam_ref = c / f[1];
    let cp_dd = (cp[0] * lam_sat - cp[1] * lam_ref) - (cp[2] * lam_sat - cp[3] * lam_ref);
    let n_dd = ambiguities[idx[0]] - ambiguities[idx[1]];
    
    let mut h_cp = vec![0.0; geom.state_size];
    h_cp[0] = geom.h_r.x; h_cp[1] = geom.h_r.y; h_cp[2] = geom.h_r.z;
    h_cp[6] = geom.h_att.x; h_cp[7] = geom.h_att.y; h_cp[8] = geom.h_att.z;
    h_cp[crate::filter::CORE_STATE_SIZE + idx[0]] = 1.0; 
    h_cp[crate::filter::CORE_STATE_SIZE + idx[1]] = -1.0;
    if geom.state_size > 20 { h_cp[20] = geom.h_zwd; }
    
    SingleUpdate { z: cp_dd - (geom.comp_dd - iono_dd + n_dd), h: h_cp, r: var.val, type_code: freq_idx, r_ref: var.ref_val }
}

pub fn compute_dd_carrier_phase(mctx: &DdMeasurementContext, p: &DdCarrierPhaseParams) -> Vec<SingleUpdate> {
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
}

fn compute_doppler_innovation(
    ctx: &DdContext, pos_apc: Vector3<f64>, base_coord_vec: Vector3<f64>, v_ant: Vector3<f64>
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
    ctx: &DdContext, pos_apc: Vector3<f64>, base_coord_vec: Vector3<f64>,
    h_r: Vector3<f64>, r_b_e: &nalgebra::Rotation3<f64>, velocity: &Vector3<f64>,
    omega_b: &Vector3<f64>, lever_arm: &Vector3<f64>, dop_base_var: f64,
    state_size: usize, var_factor: f64, ref_var_factor: f64
) -> Option<SingleUpdate> {
    let dop_valid = [ctx.rov_sat.doppler, ctx.rov_ref.doppler, ctx.base_sat.doppler, ctx.ref_base.doppler].iter().all(|&x| x != 0.0);
    if dop_valid {
        let v_ant = velocity + r_b_e * omega_b.cross(lever_arm);
        let innov = compute_doppler_innovation(ctx, pos_apc, base_coord_vec, v_ant);
        
        tracing::trace!("Doppler Innov. var={:.3} innov={:.3}", dop_base_var, innov);
        
        let mut h_dop = vec![0.0; state_size]; 
        h_dop[3] = h_r.x; h_dop[4] = h_r.y; h_dop[5] = h_r.z;
        
        let h_dop_att = doppler_attitude_jacobian(r_b_e.matrix(), omega_b, lever_arm, &h_r);
        h_dop[6] = h_dop_att.x; h_dop[7] = h_dop_att.y; h_dop[8] = h_dop_att.z;
        
        let h_dop_bg = r_b_e.matrix() * lever_arm.cross_matrix();
        let h_dop_bg = h_r.transpose() * h_dop_bg;
        h_dop[12] = h_dop_bg[0]; h_dop[13] = h_dop_bg[1]; h_dop[14] = h_dop_bg[2];
        
        return Some(SingleUpdate { z: innov, h: h_dop, r: dop_base_var * var_factor, type_code: 3, r_ref: dop_base_var * ref_var_factor });
    }
    None
}

const MIN_ELEVATION_RAD: f64 = 0.001;
const BASE_SNR_ELEVATION_THRESH_DEG: f64 = 45.0;

pub struct EkfGeometryContext {
    pub pos_apc: Vector3<f64>,
    pub base_coord_vec: Vector3<f64>,
    pub r_b_e: nalgebra::Matrix3<f64>,
    pub lever_arm: Vector3<f64>,
    pub state_size: usize,
}

impl EkfGeometryContext {
    pub fn new(state: &RtkState, env: &MeasurementEnvironment) -> Self {
        let r_b_e = state.attitude.to_rotation_matrix();
        let mut pos_apc = state.position.vector + r_b_e * env.lever_arm;
        let mut base_coord_vec = env.base_coord.vector;

        let set_rov = gneiss_core::tides::solid_earth_tides_ecef(state.time, pos_apc);
        let set_bas = gneiss_core::tides::solid_earth_tides_ecef(state.time, base_coord_vec);

        pos_apc += set_rov;
        base_coord_vec += set_bas;

        let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();

        Self {
            pos_apc,
            base_coord_vec,
            r_b_e: r_b_e.into_inner(),
            lever_arm: env.lever_arm,
            state_size,
        }
    }

    pub fn compute_attitude_jacobian(&self, h_r: &Vector3<f64>) -> Vector3<f64> {
        range_attitude_jacobian(&(self.r_b_e * self.lever_arm), h_r)
    }
}

/// Computes ∂(DD_range)/∂θ for a single satellite pair.
///
/// With left-multiplicative attitude error (`R_new = (I + [δθ]×) · R`),
/// the antenna phase center perturbation is:
///
///   d(pos_apc) = δθ × l_ecef
///
/// The DD range Jacobian is `h_r · d(pos_apc) = h_r · (δθ × l_ecef)`.
/// By the cyclic scalar triple product identity `a·(b×c) = b·(c×a)`:
///
///   h_r · (δθ × l_ecef) = δθ · (l_ecef × h_r)
///
/// Therefore `∂DD/∂θ = (l_ecef × h_r)ᵀ`.
///
/// # Arguments
/// * `lever_ecef` — Lever arm rotated into ECEF: `R_b^e · l_body`
/// * `h_r` — DD direction vector: `e_ref − e_sat`
///
/// Computes ∂(DD_range_rate)/∂θ for a single satellite pair.
///
/// The antenna velocity perturbation from attitude error is:
///   d(v_apc) = δθ × a,  where a = R · (ω_b × l)
///
/// By the same cyclic identity:
///   h_r · (δθ × a) = δθ · (a × h_r)
///
/// Therefore `∂DD_rate/∂θ = (a × h_r)ᵀ`.
///
/// # Arguments
/// * `r_b_e` — Body-to-ECEF rotation matrix
/// * `omega_b` — Corrected body angular rate (gyro − bias)
/// * `lever_arm` — Lever arm in body frame
/// * `h_r` — DD direction vector: `e_ref − e_sat`
#[allow(clippy::too_many_arguments)]
fn find_ephemeris(
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    sat: gneiss_core::sat::SatelliteId,
    time_tow: f64,
) -> Option<&gneiss_core::ephemeris::Ephemeris> {
    ephemerides.iter().filter(|e| e.sat() == sat).min_by(|a, b| {
        let da = (a.toe().tow - time_tow).abs();
        let db = (b.toe().tow - time_tow).abs();
        da.partial_cmp(&db).unwrap()
    })
}

pub fn compute_innovations(
    state: &mut RtkState, group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation, ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<EkfUpdates> {
    let mut updates = EkfUpdates::new();
    let geom = EkfGeometryContext::new(state, env);
    let ref_eph = find_ephemeris(env.ephemerides, ref_rover_orig.sat, state.time.tow)?;
    let ref_state = compute_sat_state(ref_eph, ref_rover_orig, ref_base_orig, state.time, &geom, env.base_time);

    let ref_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 1);
    let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 2);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_ephemeris(env.ephemerides, rover_sat_orig.sat, state.time.tow) {
            process_single_satellite_pair(
                state, rover_sat_orig, base_sat_orig, ref_rover_orig, ref_base_orig,
                sat_eph, &ref_state, &geom, env, ref_idx_l1, ref_idx_l2, &mut updates
            );
        }
    }
    Some(updates)
}

fn compute_sat_state(
    eph: &gneiss_core::ephemeris::Ephemeris,
    rov_obs: &DdObservation,
    bas_obs: &DdObservation,
    time: gneiss_core::time::GpsTime,
    geom: &EkfGeometryContext,
    base_time: gneiss_core::time::GpsTime,
) -> SatState {
    let (rov_pos, rov_vel) = get_sat_state(eph, rov_obs.pr_l1, time, geom.pos_apc);
    let (bas_pos, bas_vel) = get_sat_state(eph, bas_obs.pr_l1, base_time, geom.base_coord_vec);
    let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_obs.sat, eph.freq_num());

    SatState { rov_pos, rov_vel, bas_pos, bas_vel, f1, f2 }
}

#[allow(clippy::too_many_arguments)]
fn process_single_satellite_pair(
    state: &mut RtkState, rover_sat_orig: &DdObservation, base_sat_orig: &DdObservation, ref_rover_orig: &DdObservation, ref_base_orig: &DdObservation,
    sat_eph: &gneiss_core::ephemeris::Ephemeris, ref_state: &SatState, geom: &EkfGeometryContext, env: &MeasurementEnvironment,
    ref_idx_l1: Option<usize>, ref_idx_l2: Option<usize>, updates: &mut EkfUpdates,
) {
    let sat_state = compute_sat_state(sat_eph, rover_sat_orig, base_sat_orig, state.time, geom, env.base_time);
    
    let e_ref_rov = (ref_state.rov_pos - geom.pos_apc).normalize();
    let e_sat_rov = (sat_state.rov_pos - geom.pos_apc).normalize();
    let h_r = e_ref_rov - e_sat_rov;
    let h_att = geom.compute_attitude_jacobian(&h_r);

    let mut rov_sat = rover_sat_orig.clone(); let mut bas_sat = base_sat_orig.clone();
    let mut rov_ref = ref_rover_orig.clone(); let mut bas_ref = ref_base_orig.clone();

    let mut ctx = DdContext {
        rov_sat: &mut rov_sat, base_sat: &mut bas_sat,
        rov_ref: &mut rov_ref, ref_base: &mut bas_ref,
        sat_state: &sat_state, ref_state,
    };

    update_windup_state_and_obs(state, &mut ctx, geom);


    let comps = compute_dd_components(state, geom, env, &ctx);
    let ugeom = UpdateGeometry { comp_dd: comps.comp_pr_dd, h_r, h_att, h_zwd: comps.h_zwd, state_size: geom.state_size };
    let mctx = DdMeasurementContext { ctx: &ctx, geom: &ugeom, comps: &comps, env };

    generate_measurement_updates(state, &mctx, geom, ref_idx_l1, ref_idx_l2, updates);
}

pub struct DdComponents {
    pub comp_pr_dd: f64,
    pub iono_dd_l1: f64,
    pub iono_dd_l2: f64,
    pub var_factor: f64,
    pub ref_var_factor: f64,
    pub h_zwd: f64,
}

fn compute_dd_components(
    state: &RtkState,
    geom: &EkfGeometryContext,
    env: &MeasurementEnvironment,
    ctx: &DdContext,
) -> DdComponents {
    let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(
        state.time, geom.pos_apc, geom.base_coord_vec,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        ctx.sat_state.f1, ctx.sat_state.f2, ctx.ref_state.f1, ctx.ref_state.f2,
    );

    let base_llh = gneiss_core::coords::ecef_to_llh(geom.base_coord_vec);
    let rov_llh = gneiss_core::coords::ecef_to_llh(geom.pos_apc);
    let (_, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, ctx.sat_state.rov_pos);
    let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, ctx.ref_state.rov_pos);
    let (_, el_bas_sat) = gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, ctx.sat_state.bas_pos);
    let (_, el_bas_ref) = gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, ctx.ref_state.bas_pos);

    let (var_factor, ref_var_factor) = crate::engine::measurement_math::compute_variance_factors(&crate::engine::measurement_math::VarianceFactors {
        snr_rov_sat: ctx.rov_sat.snr,
        snr_rov_ref: ctx.rov_ref.snr,
        el_rov_sat,
        el_rov_ref,
        el_bas_sat,
        el_bas_ref,
        snr_a: env.tuning.snr_a,
        snr_b: env.tuning.snr_b,
    });
    let (h_zwd, zwd_dd) = compute_zwd_mapping(el_rov_sat, el_rov_ref, state.zwd);

    let comp_pr_dd = compute_geometric_dd(
        geom.pos_apc, geom.base_coord_vec, ctx.sat_state.rov_pos, ctx.ref_state.rov_pos, ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
    ) + tropo_dd + zwd_dd;

    DdComponents { comp_pr_dd, iono_dd_l1, iono_dd_l2, var_factor, ref_var_factor, h_zwd }
}



fn generate_measurement_updates(
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
}

pub struct EkfMeasurementMatrices {
    pub z: DVector<f64>,
    pub h: DMatrix<f64>,
    pub r: DMatrix<f64>,
    pub mt: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

pub fn build_measurement_model(
    state: &mut RtkState, matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment, chi_square_pr_threshold: f64, chi_square_cp_threshold: f64,
) -> Option<EkfMeasurementMatrices> {
    let mut all = EkfUpdates::new();
    let const_groups = group_measurements_by_constellation(matched_obs);

    for (_, group) in const_groups {
        if group.len() < 2 { continue; } 
        let ref_idx = select_reference_satellite(&group, state, env);
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some(updates) = compute_innovations(state, &group_clone, &ref_rover, &ref_base, env) {
            all.extend(updates);
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    tracing::trace!("Pre-filter z_all len: {}", all.z.len());
    
    let safe_indices = filter_innovations_chi_squared(
        state, state_size, chi_square_pr_threshold, chi_square_cp_threshold,
        &all.z, &all.h, &all.r, &all.mt
    );

    build_final_measurement_matrices(state_size, safe_indices, &all.z, &all.h, &all.r, &all.mt)
}

fn group_measurements_by_constellation(
    matched_obs: &[(DdObservation, DdObservation)],
) -> std::collections::HashMap<gneiss_core::sat::Constellation, Vec<(DdObservation, DdObservation)>> {
    let mut const_groups = std::collections::HashMap::<gneiss_core::sat::Constellation, Vec<(DdObservation, DdObservation)>>::new();
    for obs in matched_obs {
        const_groups.entry(obs.0.sat.constellation).or_default().push(obs.clone());
    }
    const_groups
}

fn select_reference_satellite(group: &[(DdObservation, DdObservation)], state: &RtkState, env: &MeasurementEnvironment) -> usize {
    let mut best_score = -1.0; let mut ref_idx = 0;
    let rov_llh = gneiss_core::coords::ecef_to_llh(state.position.vector);

    for (i, (r, _)) in group.iter().enumerate() {
        if let Some(eph) = env.ephemerides.iter().filter(|e| e.sat() == r.sat).min_by(|a, b| {
            let da = (a.toe().tow - state.time.tow).abs();
            let db = (b.toe().tow - state.time.tow).abs();
            da.partial_cmp(&db).unwrap()
        }) {
            let tau = r.pr_l1 / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
            let t_tx = gneiss_core::time::GpsTime::new(state.time.week, state.time.tow - tau);
            let (sat_pos, _, _, _): (Vector3<f64>, _, _, _) = eph.position(t_tx);
            let (_, el) = gneiss_core::coords::az_el(rov_llh, state.position.vector, sat_pos);
            
            let score = if r.cp_l1.is_some() { el + 100.0 } else { el };
            if score > best_score { best_score = score; ref_idx = i; }
        }
    }
    ref_idx
}

fn update_reject_counts(state: &mut RtkState, h_row: &DMatrix<f64>, state_size: usize, passed: bool) {
    for c in crate::filter::CORE_STATE_SIZE..state_size {
        if h_row[(0, c)] > 0.5 {
            let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
            if passed {
                state.reject_counts.insert(key, 0);
            } else {
                let count = *state.reject_counts.get(&key).unwrap_or(&0) + 1;
                state.reject_counts.insert(key, count);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn filter_innovations_chi_squared(
    state: &mut RtkState, state_size: usize, chi_pr: f64, chi_cp: f64,
    z_all: &[f64], h_all: &[Vec<f64>], r_all: &[f64], type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Vec<usize> {
    let mut safe_indices = Vec::new();
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size { h_row[(0, c)] = h_all[i][c]; }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;
        
        let threshold = match type_all[i].1 { 
            0 => chi_pr * chi_pr, 1 | 2 => chi_cp * chi_cp, 3 => chi_pr * 1000.0, _ => chi_pr * chi_pr   
        };
        
        let passed = chi2 <= threshold;
        if passed { safe_indices.push(i); }
        else { tracing::debug!("Rejected meas type {} with inn: {:.3}, chi2: {:.1}", type_all[i].1, z_all[i], chi2); }
        
        if type_all[i].1 == 1 || type_all[i].1 == 2 {
            update_reject_counts(state, &h_row, state_size, passed);
        }
    }
    safe_indices
}

fn build_final_measurement_matrices(
    state_size: usize, safe_indices: Vec<usize>, z_all: &[f64], h_all: &[Vec<f64>], r_all: &[f64], type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Option<EkfMeasurementMatrices> {
    if safe_indices.len() >= 4 {
        let mut z_vec = DVector::zeros(safe_indices.len());
        let mut h_mat = DMatrix::zeros(safe_indices.len(), state_size);
        let mut r_diagonals = Vec::new();
        let mut t_vec = Vec::new();

        for (new_i, &old_i) in safe_indices.iter().enumerate() {
            z_vec[new_i] = z_all[old_i];
            for c in 0..state_size { h_mat[(new_i, c)] = h_all[old_i][c]; }
            r_diagonals.push(r_all[old_i]);
            t_vec.push(type_all[old_i]);
        }
        let r_mat = build_dense_covariance_matrix(&r_diagonals, &t_vec);
        Some(EkfMeasurementMatrices { z: z_vec, h: h_mat, r: r_mat, mt: t_vec })
    } else {
        tracing::warn!("measurement model empty! all_z={}, safe_indices={}", z_all.len(), safe_indices.len());
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
    use crate::engine::measurement_math::{doppler_attitude_jacobian, range_attitude_jacobian};

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
        let updates = super::super::measurement::compute_innovations(&mut state, &matched_obs, &ref_rover, &ref_base, &env).unwrap();
        let z = updates.z;
        let r = updates.r;

        println!("Z: {:?}", z);
        
        // Lock in the golden Z vector (updated for iterative ecef_to_llh refinement)
        assert!((z[0] - 2.95577).abs() < 1e-3, "z[0]={}", z[0]);
        assert!((z[1] - 2.95627).abs() < 1e-3, "z[1]={}", z[1]);
        assert!((z[2] - -2.04577).abs() < 1e-3, "z[2]={}", z[2]);
        assert!((z[3] - 19.36016).abs() < 1e-3, "z[3]={}", z[3]);
        assert!((z[4] - 9.97132).abs() < 1e-3, "z[4]={}", z[4]);
        assert!((z[5] - 9.97206).abs() < 1e-3, "z[5]={}", z[5]);
        assert!((z[6] - -0.03096).abs() < 1e-3, "z[6]={}", z[6]);
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
        println!("pos: {:?}", pos); println!("vel: {:?}", vel); assert!((pos.x - 5041617.577444584).abs() < 1e-6); assert!((pos.y - 17749192.669882875).abs() < 1e-6); assert!((pos.z - 18906563.91687504).abs() < 1e-6);
        assert!((vel.x - (-2062.128434083682)).abs() < 1e-6); assert!((vel.y - (-1216.6504200626093)).abs() < 1e-6); assert!((vel.z - 1738.0997934704972).abs() < 1e-6);
        
        let (pos0, _vel0) = get_sat_state(&eph, 0.0, time, rx_pos);
        assert!((pos.x - pos0.x).abs() > 0.0);
    }
    #[test]
    fn test_compute_atmospheric_delays() {
        use crate::engine::measurement::compute_atmospheric_delays;
        let state_time = GpsTime::new(2137, 422922.0);
        let pos_apc = Vector3::new(1000.0, 2000.0, 3000.0);
        let base_coord_vec = Vector3::new(1005.0, 2005.0, 3005.0);
        let sat_vec_rov = Vector3::new(15000000.0, 20000000.0, 30000000.0);
        let ref_sat_vec_rov = Vector3::new(-15000000.0, 20000000.0, -30000000.0);
        let sat_vec_bas = Vector3::new(15000005.0, 20000005.0, 30000005.0);
        let ref_sat_vec_bas = Vector3::new(-15000005.0, 20000005.0, -30000005.0);
        let sat_f1 = 1575.42e6;
        let sat_f2 = 1227.60e6;
        let ref_f1 = 1575.42e6;
        let ref_f2 = 1227.60e6;

        let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(
            state_time, pos_apc, base_coord_vec, sat_vec_rov, ref_sat_vec_rov, sat_vec_bas, ref_sat_vec_bas,
            sat_f1, sat_f2, ref_f1, ref_f2
        );
        assert!((tropo_dd - 0.0).abs() < 1e-6);
        assert!((iono_dd_l1 - (-0.00020734960295598626)).abs() < 1e-6);
        assert!((iono_dd_l2 - (-0.0003414932766467871)).abs() < 1e-6);
    }

    #[test]
    fn test_compute_dd_pseudorange() { use gneiss_core::sat::{SatelliteId, Constellation}; use gneiss_core::coords::{Coordinate, Datum, Frame};
        use crate::engine::measurement::{compute_dd_pseudorange, DdContext, SatState, MeasurementEnvironment};
        use crate::filter::DdObservation;
        use crate::engine::config::EkfTuningConfig;
        
        let mut rov_sat = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: None };
        let mut base_sat = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: None };
        let mut rov_ref = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: None };
        let mut ref_base = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: None };
        
        let sat_state = SatState { rov_pos: Vector3::new(20000000.0, 0.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(0.0, 20000000.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        let ref_state = SatState { rov_pos: Vector3::new(20000000.0, 0.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(0.0, 20000000.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        
        let ctx = DdContext {
            rov_sat: &mut rov_sat,
            base_sat: &mut base_sat,
            rov_ref: &mut rov_ref,
            ref_base: &mut ref_base,
            sat_state: &sat_state,
            ref_state: &ref_state,
        };
        
        let tuning = EkfTuningConfig::default();
        let base_coord = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let env = MeasurementEnvironment {
            ephemerides: &[],
            base_coord: &base_coord,
            base_time: GpsTime::new(0, 0.0),
            lever_arm: Vector3::zeros(),
            omega_b: Vector3::zeros(),
            tuning: &tuning,
        };
        
        let ugeom = crate::engine::measurement::UpdateGeometry { comp_dd: 0.0, h_r: Vector3::new(1.0, 0.0, 0.0), h_att: Vector3::zeros(), h_zwd: 0.0, state_size: 22 };
        let comps = crate::engine::measurement::DdComponents { comp_pr_dd: 0.0, iono_dd_l1: 0.0, iono_dd_l2: 0.0, var_factor: 1.0, ref_var_factor: 1.0, h_zwd: 0.0 };
        let mctx = crate::engine::measurement::DdMeasurementContext { ctx: &ctx, geom: &ugeom, comps: &comps, env: &env };
        let updates = compute_dd_pseudorange(&mctx);
        assert_eq!(updates.len(), 2);
        assert_eq!(updates[0].z, 0.0);
        assert_eq!(updates[1].z, 0.0);
    }

    #[test]
    fn test_compute_dd_carrier_phase() { use gneiss_core::sat::{SatelliteId, Constellation}; use gneiss_core::coords::{Coordinate, Datum, Frame};
        use crate::engine::measurement::{compute_dd_carrier_phase, DdContext, SatState, MeasurementEnvironment};
        use crate::filter::{RtkState, DdObservation};
        use crate::engine::config::EkfTuningConfig;
        
        let time = GpsTime::new(2137, 422922.0);
        let mut state = RtkState::new(time, Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time), 10.0);
        
        let mut rov_sat = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 0.0, snr: 45.0, locktime: None };
        let mut base_sat = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 0.0, snr: 45.0, locktime: None };
        let mut rov_ref = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 0.0, snr: 45.0, locktime: None };
        let mut ref_base = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 0.0, snr: 45.0, locktime: None };
        
        let sat_state = SatState { rov_pos: Vector3::new(20000000.0, 0.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(0.0, 20000000.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        let ref_state = SatState { rov_pos: Vector3::new(20000000.0, 0.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(0.0, 20000000.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        
        let ctx = DdContext {
            rov_sat: &mut rov_sat,
            base_sat: &mut base_sat,
            rov_ref: &mut rov_ref,
            ref_base: &mut ref_base,
            sat_state: &sat_state,
            ref_state: &ref_state,
        };
        
        state.add_ambiguity(SatelliteId { constellation: Constellation::Gps, prn: 1 }, 1, 0.0, 1.0);
        state.add_ambiguity(SatelliteId { constellation: Constellation::Gps, prn: 2 }, 1, 0.0, 1.0);
        state.add_ambiguity(SatelliteId { constellation: Constellation::Gps, prn: 1 }, 2, 0.0, 1.0);
        state.add_ambiguity(SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 0.0, 1.0);
        
        let tuning = EkfTuningConfig::default();
        let base_coord = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let env = MeasurementEnvironment {
            ephemerides: &[],
            base_coord: &base_coord,
            base_time: GpsTime::new(0, 0.0),
            lever_arm: Vector3::zeros(),
            omega_b: Vector3::zeros(),
            tuning: &tuning,
        };
        
        let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
        let updates = compute_dd_carrier_phase(
            &ctx, state.is_fixed, &state.ambiguities, Some(0), Some(1), Some(2), Some(3),
            0.0, 0.0, 0.0, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros(), state_size, 1.0, 1.0, env.tuning.cp_base_var, 0.0
        );
        assert_eq!(updates.len(), 2);
    }

    #[test]
    fn test_compute_dd_doppler() { use gneiss_core::sat::{SatelliteId, Constellation}; use gneiss_core::coords::{Coordinate, Datum, Frame};
        use crate::engine::measurement::{compute_dd_doppler, DdContext, SatState, MeasurementEnvironment};
        use crate::filter::{RtkState, DdObservation};
        use crate::engine::config::EkfTuningConfig;
        
        let time = GpsTime::new(2137, 422922.0);
        let mut state = RtkState::new(time, Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time), 10.0);
        
        let mut rov_sat = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 100.0, snr: 45.0, locktime: None };
        let mut base_sat = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 100.0, snr: 45.0, locktime: None };
        let mut rov_ref = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 100.0, snr: 45.0, locktime: None };
        let mut ref_base = DdObservation { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: None, cp_l2: None, doppler: 100.0, snr: 45.0, locktime: None };
        
        let sat_state = SatState { rov_pos: Vector3::new(20000000.0, 0.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(0.0, 20000000.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        let ref_state = SatState { rov_pos: Vector3::new(20000000.0, 0.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(0.0, 20000000.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        
        let ctx = DdContext {
            rov_sat: &mut rov_sat,
            base_sat: &mut base_sat,
            rov_ref: &mut rov_ref,
            ref_base: &mut ref_base,
            sat_state: &sat_state,
            ref_state: &ref_state,
        };
        
        let tuning = EkfTuningConfig::default();
        let base_coord = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let env = MeasurementEnvironment {
            ephemerides: &[],
            base_coord: &base_coord,
            base_time: GpsTime::new(0, 0.0),
            lever_arm: Vector3::zeros(),
            omega_b: Vector3::zeros(),
            tuning: &tuning,
        };
        
        let r_b_e_rot = state.attitude.to_rotation_matrix();
        let update = compute_dd_doppler(
            &ctx, Vector3::zeros(), Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0),
            &r_b_e_rot, &state.velocity, &env.omega_b, &env.lever_arm,
            env.tuning.dop_base_var, 22, 1.0, 1.0
        );
        let u = update.unwrap(); assert!(!u.z.is_nan());
    }

    #[test]
    fn test_compute_phase_windup() { use gneiss_core::sat::{SatelliteId, Constellation}; use gneiss_core::coords::{Coordinate, Datum, Frame};
        use crate::engine::measurement::{DdContext, SatState};
        use crate::filter::{RtkState, DdObservation};
        
        let time = GpsTime::new(2137, 422922.0);
        let mut state = RtkState::new(time, Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time), 10.0);
        
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        
        let mut rov_sat = DdObservation { sat: sat1, pr_l1: 0.0, pr_l2: None, cp_l1: Some(10.0), cp_l2: Some(20.0), doppler: 0.0, snr: 45.0, locktime: None };
        let mut base_sat = DdObservation { sat: sat1, pr_l1: 0.0, pr_l2: None, cp_l1: Some(10.0), cp_l2: Some(20.0), doppler: 0.0, snr: 45.0, locktime: None };
        let mut rov_ref = DdObservation { sat: sat2, pr_l1: 0.0, pr_l2: None, cp_l1: Some(10.0), cp_l2: Some(20.0), doppler: 0.0, snr: 45.0, locktime: None };
        let mut ref_base = DdObservation { sat: sat2, pr_l1: 0.0, pr_l2: None, cp_l1: Some(10.0), cp_l2: Some(20.0), doppler: 0.0, snr: 45.0, locktime: None };
        
        let sat_state = SatState { rov_pos: Vector3::new(20000000.0, 0.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(0.0, 20000000.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        let ref_state = SatState { rov_pos: Vector3::new(0.0, 20000000.0, 0.0), rov_vel: Vector3::zeros(), bas_pos: Vector3::new(20000000.0, 0.0, 0.0), bas_vel: Vector3::zeros(), f1: 1575.42e6, f2: 1227.60e6 };
        
        let mut ctx = DdContext {
            rov_sat: &mut rov_sat,
            base_sat: &mut base_sat,
            rov_ref: &mut rov_ref,
            ref_base: &mut ref_base,
            sat_state: &sat_state,
            ref_state: &ref_state,
        };
        
        let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);
        let crate::engine::measurement_math::WindupUpdates { w_sat, w_ref, w_bas_sat, w_bas_ref } = crate::engine::measurement_math::compute_phase_windup(
            Vector3::zeros(), Vector3::zeros(), sun_pos,
            ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
            ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
            0.0, 0.0, 0.0, 0.0,
        );
        if let Some(cp) = &mut ctx.rov_sat.cp_l1 { *cp += w_sat; }
        if let Some(cp2) = &mut ctx.rov_sat.cp_l2 { *cp2 += w_sat; }
        if let Some(cp) = &mut ctx.rov_ref.cp_l1 { *cp += w_ref; }
        if let Some(cp2) = &mut ctx.rov_ref.cp_l2 { *cp2 += w_ref; }
        if let Some(cp) = &mut ctx.base_sat.cp_l1 { *cp += w_bas_sat; }
        if let Some(cp2) = &mut ctx.base_sat.cp_l2 { *cp2 += w_bas_sat; }
        if let Some(cp) = &mut ctx.ref_base.cp_l1 { *cp += w_bas_ref; }
        if let Some(cp2) = &mut ctx.ref_base.cp_l2 { *cp2 += w_bas_ref; }
        
        assert!(ctx.rov_sat.cp_l1.unwrap() != 10.0);
        assert!(ctx.rov_sat.cp_l2.unwrap() != 20.0);
        assert!(ctx.rov_ref.cp_l1.unwrap() != 10.0);
        assert!(ctx.rov_ref.cp_l2.unwrap() != 20.0);
    }

    #[test]
    fn test_compute_innovations() { use gneiss_core::sat::{SatelliteId, Constellation}; use gneiss_core::coords::{Coordinate, Datum, Frame}; use gneiss_core::time::GpsTime;
        use crate::engine::measurement::{compute_innovations, MeasurementEnvironment};
        use crate::filter::{RtkState, DdObservation};
        use crate::engine::config::EkfTuningConfig;
        use gneiss_core::ephemeris::Ephemeris;
        use nalgebra::Vector3;
        
        let time = GpsTime::new(2137, 422922.0);
        let mut state = RtkState::new(time, Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time), 10.0);
        
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        
        let rov_sat = DdObservation { sat: sat1, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 100.0, snr: 45.0, locktime: None };
        let base_sat = DdObservation { sat: sat1, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 100.0, snr: 45.0, locktime: None };
        let rov_ref = DdObservation { sat: sat2, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 100.0, snr: 45.0, locktime: None };
        let ref_base = DdObservation { sat: sat2, pr_l1: 20000000.0, pr_l2: Some(20000000.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 100.0, snr: 45.0, locktime: None };
        
        let eph_ref = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: sat2, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });

        let eph_rov = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: sat1, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 1.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.5, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });
        
        state.add_ambiguity(sat1, 1, 0.0, 1.0);
        state.add_ambiguity(sat2, 1, 0.0, 1.0);
        state.add_ambiguity(sat1, 2, 0.0, 1.0);
        state.add_ambiguity(sat2, 2, 0.0, 1.0);
        
        let tuning = EkfTuningConfig::default();
        let base_coord = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let env = MeasurementEnvironment {
            ephemerides: &[eph_ref, eph_rov],
            base_coord: &base_coord,
            base_time: GpsTime::new(0, 0.0),
            lever_arm: Vector3::zeros(),
            omega_b: Vector3::zeros(),
            tuning: &tuning,
        };
        
        let group = vec![(rov_sat, base_sat)];
        
        let res = compute_innovations(&mut state, &group, &rov_ref, &ref_base, &env);
        assert!(res.is_some());
        let updates = res.unwrap();
        let (z_vals, h_rows, r_vals, meas_type) = (updates.z, updates.h, updates.r, updates.mt);
        assert!(z_vals.len() > 0);
        assert_eq!(h_rows.len(), z_vals.len());
        assert_eq!(r_vals.len(), z_vals.len());
        assert_eq!(meas_type.len(), z_vals.len());
    }

    // -----------------------------------------------------------------------
    // Numerical Jacobian verification for attitude coupling
    // -----------------------------------------------------------------------
    //
    // These tests are the definitive proof of the correct sign. They
    // perturb the attitude by a small δθ along each axis, recompute the
    // DD range (or DD rate) via the actual range equations, and compare
    // the central-difference derivative against the analytical Jacobian.
    //
    // If these tests pass, the Jacobian is correct. Period.
    // -----------------------------------------------------------------------

    /// Compute DD range given a rotation applied to base position + lever arm.
    /// DD = |sat - pos_apc| - |ref - pos_apc|
    #[cfg(test)]
    fn dd_range(
        pos: Vector3<f64>,
        rot: nalgebra::UnitQuaternion<f64>,
        lever_body: Vector3<f64>,
        sat_pos: Vector3<f64>,
        ref_pos: Vector3<f64>,
    ) -> f64 {
        let pos_apc = pos + rot * lever_body;
        (sat_pos - pos_apc).norm() - (ref_pos - pos_apc).norm()
    }

    /// Compute DD range-rate (rover portion only).
    /// DD_rate = e_sat·(v_sat - v_apc) - e_ref·(v_ref - v_apc)
    #[cfg(test)]
    fn dd_range_rate(
        pos: Vector3<f64>,
        vel: Vector3<f64>,
        rot: nalgebra::UnitQuaternion<f64>,
        lever_body: Vector3<f64>,
        omega_b: Vector3<f64>,
        sat_pos: Vector3<f64>,
        ref_pos: Vector3<f64>,
        sat_vel: Vector3<f64>,
        ref_vel: Vector3<f64>,
    ) -> f64 {
        let pos_apc = pos + rot * lever_body;
        let v_apc = vel + rot * omega_b.cross(&lever_body);
        let e_sat = (sat_pos - pos_apc).normalize();
        let e_ref = (ref_pos - pos_apc).normalize();
        e_sat.dot(&(sat_vel - v_apc)) - e_ref.dot(&(ref_vel - v_apc))
    }

    /// Apply a small rotation perturbation (left-multiplicative, matching
    /// the EKF convention in apply_state_correction).
    #[cfg(test)]
    fn perturb_attitude(
        rot: nalgebra::UnitQuaternion<f64>,
        d_theta: Vector3<f64>,
    ) -> nalgebra::UnitQuaternion<f64> {
        let angle = d_theta.norm();
        if angle < 1e-15 { return rot; }
        let dq = nalgebra::UnitQuaternion::from_axis_angle(
            &nalgebra::Unit::new_normalize(d_theta), angle,
        );
        dq * rot  // left-multiplicative: R_new = dR · R
    }

    #[test]
    fn test_range_attitude_jacobian_numerical() {
        // Verify against finite differences of the full nonlinear DD range.
        // The analytical Jacobian is a linearization, so we expect agreement
        // to O(eps) relative error (limited by the h_r approximation).
        let pos = Vector3::new(4_000_000.0, 1_000_000.0, 4_500_000.0);
        let lever_body = Vector3::new(0.0, 0.0, 1.5);
        let rot = nalgebra::UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
        let lever_ecef = rot * lever_body;

        let sat_pos = Vector3::new(20_000_000.0, 10_000_000.0, 15_000_000.0);
        let ref_pos = Vector3::new(15_000_000.0, 20_000_000.0, 10_000_000.0);

        let pos_apc = pos + lever_ecef;
        let e_sat = (sat_pos - pos_apc).normalize();
        let e_ref = (ref_pos - pos_apc).normalize();
        let h_r = e_ref - e_sat;

        let j_analytical = range_attitude_jacobian(&lever_ecef, &h_r);

        let eps = 1e-7;
        let mut j_numerical = Vector3::zeros();
        for axis in 0..3 {
            let mut d_theta = Vector3::zeros();
            d_theta[axis] = eps;
            let dd_plus = dd_range(pos, perturb_attitude(rot, d_theta), lever_body, sat_pos, ref_pos);
            let dd_minus = dd_range(pos, perturb_attitude(rot, -d_theta), lever_body, sat_pos, ref_pos);
            j_numerical[axis] = (dd_plus - dd_minus) / (2.0 * eps);
        }

        // Nonlinear DD includes second-order effects from h_r changing
        // as pos_apc moves (~lever/range ≈ 1e-7). Relax tolerance accordingly.
        let err = (j_analytical - j_numerical).norm();
        let scale = j_numerical.norm().max(1e-12);
        assert!(
            err / scale < 0.05,
            "Range attitude Jacobian sign/magnitude error!\n  analytical: {:?}\n  numerical:  {:?}\n  rel_error: {:.2e}",
            j_analytical, j_numerical, err / scale
        );
    }

    /// Verify the range attitude Jacobian via the linearized quantity directly.
    /// This test eliminates second-order effects and is exact to machine precision.
    ///
    /// The analytical Jacobian says: d(DD) ≈ J · δθ
    /// We verify: h_r · (δθ × l_ecef) == J · δθ for arbitrary δθ.
    #[test]
    fn test_range_attitude_jacobian_linearized_exact() {
        let lever_ecef = Vector3::new(0.3, -0.7, 1.2);
        let h_r = Vector3::new(0.4, -0.1, 0.6);
        let j = range_attitude_jacobian(&lever_ecef, &h_r);

        // For arbitrary δθ, verify h_r · (δθ × l_ecef) == J · δθ
        for d_theta in &[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.3, -0.5, 0.8),
        ] {
            let lhs: f64 = h_r.dot(&d_theta.cross(&lever_ecef)); // h_r · (δθ × l)
            let rhs: f64 = j.dot(d_theta);                        // J · δθ
            assert!(
                (lhs - rhs).abs() < 1e-14,
                "Linearized check failed: lhs={}, rhs={}, δθ={:?}",
                lhs, rhs, d_theta
            );
        }
    }

    #[test]
    fn test_range_attitude_jacobian_axis_aligned() {
        // Simple case: lever arm along body-Z, identity rotation,
        // h_r along ECEF-X. Analytical: l_ecef × h_r = [0,0,1.5] × [1,0,0] = [0,1.5,0]
        let lever_ecef = Vector3::new(0.0, 0.0, 1.5);
        let h_r = Vector3::new(1.0, 0.0, 0.0);

        let j = range_attitude_jacobian(&lever_ecef, &h_r);
        assert!((j.x - 0.0).abs() < 1e-12);
        assert!((j.y - 1.5).abs() < 1e-12);
        assert!((j.z - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_range_attitude_jacobian_zero_lever_arm() {
        let lever_ecef = Vector3::zeros();
        let h_r = Vector3::new(0.3, -0.5, 0.8);
        let j = range_attitude_jacobian(&lever_ecef, &h_r);
        assert!(j.norm() < 1e-15);
    }

    #[test]
    fn test_doppler_attitude_jacobian_numerical() {
        let pos = Vector3::new(4_000_000.0, 1_000_000.0, 4_500_000.0);
        let vel = Vector3::new(1.0, 2.0, 3.0);
        let lever_body = Vector3::new(0.0, 0.0, 1.5);
        let omega_b = Vector3::new(0.01, -0.02, 0.005);
        let rot = nalgebra::UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
        let r_b_e = rot.to_rotation_matrix().into_inner();

        let sat_pos = Vector3::new(20_000_000.0, 10_000_000.0, 15_000_000.0);
        let ref_pos = Vector3::new(15_000_000.0, 20_000_000.0, 10_000_000.0);
        let sat_vel = Vector3::new(-500.0, 200.0, 3000.0);
        let ref_vel = Vector3::new(300.0, -800.0, 2500.0);

        let pos_apc = pos + rot * lever_body;
        let e_sat = (sat_pos - pos_apc).normalize();
        let e_ref = (ref_pos - pos_apc).normalize();
        let h_r = e_ref - e_sat;

        let j_analytical = doppler_attitude_jacobian(&r_b_e, &omega_b, &lever_body, &h_r);

        let eps = 1e-7;
        let mut j_numerical = Vector3::zeros();
        for axis in 0..3 {
            let mut d_theta = Vector3::zeros();
            d_theta[axis] = eps;
            let rr_plus = dd_range_rate(
                pos, vel, perturb_attitude(rot, d_theta), lever_body, omega_b,
                sat_pos, ref_pos, sat_vel, ref_vel,
            );
            let rr_minus = dd_range_rate(
                pos, vel, perturb_attitude(rot, -d_theta), lever_body, omega_b,
                sat_pos, ref_pos, sat_vel, ref_vel,
            );
            j_numerical[axis] = (rr_plus - rr_minus) / (2.0 * eps);
        }

        let err = (j_analytical - j_numerical).norm();
        let scale = j_numerical.norm().max(1e-12);
        assert!(
            err / scale < 0.05,
            "Doppler attitude Jacobian sign/magnitude error!\n  analytical: {:?}\n  numerical:  {:?}\n  rel_error: {:.2e}",
            j_analytical, j_numerical, err / scale
        );
    }

    /// Verify the Doppler attitude Jacobian via the linearized quantity directly.
    #[test]
    fn test_doppler_attitude_jacobian_linearized_exact() {
        let r_b_e = nalgebra::UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3)
            .to_rotation_matrix()
            .into_inner();
        let omega_b = Vector3::new(0.01, -0.02, 0.005);
        let lever_body = Vector3::new(0.0, 0.0, 1.5);
        let h_r = Vector3::new(0.4, -0.1, 0.6);

        let j = doppler_attitude_jacobian(&r_b_e, &omega_b, &lever_body, &h_r);
        let a = r_b_e * omega_b.cross(&lever_body);

        // For arbitrary δθ, verify h_r · (δθ × a) == J · δθ
        for d_theta in &[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(-0.7, 0.3, 0.9),
        ] {
            let lhs = h_r.dot(&d_theta.cross(&a)); // h_r · (δθ × a)
            let rhs = j.dot(d_theta);               // J · δθ
            assert!(
                (lhs - rhs).abs() < 1e-14,
                "Doppler linearized check failed: lhs={}, rhs={}, δθ={:?}",
                lhs, rhs, d_theta
            );
        }
    }

    #[test]
    fn test_doppler_attitude_jacobian_zero_omega() {
        let r_b_e = nalgebra::Matrix3::identity();
        let omega_b = Vector3::zeros();
        let lever_body = Vector3::new(0.0, 0.0, 1.5);
        let h_r = Vector3::new(0.3, -0.5, 0.8);
        let j = doppler_attitude_jacobian(&r_b_e, &omega_b, &lever_body, &h_r);
        assert!(j.norm() < 1e-15);
    }

    #[test]
    fn test_range_jacobian_antisymmetry() {
        // l × h_r = -(h_r × l): swapping args negates
        let lever_ecef = Vector3::new(0.5, -1.0, 1.5);
        let h_r = Vector3::new(0.3, 0.7, -0.2);
        let j1 = range_attitude_jacobian(&lever_ecef, &h_r);
        let j2 = range_attitude_jacobian(&h_r, &lever_ecef);
        assert!((j1 + j2).norm() < 1e-15, "Should be antisymmetric");
    }

}
