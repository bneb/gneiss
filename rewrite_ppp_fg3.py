import re

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    orig = f.read()

# Replace log_ppp_convergence
t_log = """    fn log_ppp_convergence(state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, x_pred: &DVector<f64>, p_pred: &DMatrix<f64>, solver: &PppFactorGraph) {"""
r_log = """    fn log_ppp_convergence(state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, x_pred: &DVector<f64>, p_pred: &DMatrix<f64>, solver: &PppFactorGraph) {"""

# The solve function
t_solve = """    pub fn solve(&self, state: &mut RtkState, sats: &[ProcessedSat]) -> Result<(), EngineError> {"""

# Let's break solve into solve_inner and solve.
new_solve = """    pub fn solve(&self, state: &mut RtkState, sats: &[ProcessedSat]) -> Result<(), EngineError> {
        let mut outer_iter = 0;
        loop {
            if self.solve_inner(state, sats)? { break; }
            outer_iter += 1;
            if outer_iter > 3 { tracing::warn!("Max outlier rejection iterations reached."); break; }
        }
        if state.epoch_count > 10 {
            if let Err(e) = self.resolve_cascade_ar(state, sats) { tracing::info!("Cascade AR did not fix: {:?}", e); }
        }
        Ok(())
    }

    fn solve_inner(&self, state: &mut RtkState, sats: &[ProcessedSat]) -> Result<bool, EngineError> {
        let x_pred = extract_state_vector(state);
        let p_pred = state.covariance.clone();
        state.full_x_predict = Some(x_pred.clone());
        state.full_p_predict = Some(p_pred.clone());
        let mut x_i = x_pred.clone();
        let p_inv = invert_matrix(&p_pred).ok_or(EngineError::StateDisappeared)?;

        for _iter in 0..self.max_iterations {
            if let Some(dx) = self.compute_iteration_dx(state, sats, &x_i, &x_pred, &p_inv, _iter)? {
                x_i = &x_i + &dx;
                if dx.norm() < self.convergence_threshold { break; }
            } else { break; }
        }

        if let Some(sat) = self.find_worst_outlier(state, sats, &x_i) {
            tracing::warn!("PPP FG Outlier Detected for {:?}. Removing ambiguity and retrying.", sat);
            for i in 0..4 { state.remove_ambiguity(sat, i); }
            return Ok(false);
        }

        let final_p = self.compute_final_covariance(state, sats, &x_i, &p_pred, &p_inv);
        apply_state_vector(state, &x_i, final_p);
        log_ppp_convergence(state, sats, &x_i, &x_pred, &p_pred, self);
        Ok(true)
    }

    fn find_worst_outlier(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>) -> Option<gneiss_core::sat::SatelliteId> {
        let final_meas = self.build_measurements(state, sats, x_i, self.max_iterations);
        let mut worst_sat = None;
        let mut max_norm = 15.0;
        for m in &final_meas {
            if m.is_phase {
                let norm = m.res.abs() / m.raw_var.sqrt();
                if norm > max_norm { max_norm = norm; worst_sat = m.sat; }
            }
        }
        worst_sat
    }"""

start_solve = orig.find("    pub fn solve(&self, state: &mut RtkState, sats: &[ProcessedSat]) -> Result<(), EngineError> {")
end_solve = orig.find("    pub fn resolve_cascade_ar(&self, state: &mut RtkState, sats: &[ProcessedSat]) -> Result<(), &'static str> {")

if start_solve != -1 and end_solve != -1:
    orig = orig[:start_solve] + new_solve + "\n\n" + orig[end_solve:]

# Now resolve_widelane_ar
t_rwl = """    fn resolve_widelane_ar(&self, state: &RtkState, subset: &[((gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64), (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64))], x: &DVector<f64>) -> Result<(DVector<f64>, DMatrix<f64>), &'static str> {"""
new_rwl = """    fn resolve_widelane_ar(&self, state: &RtkState, subset: &[((gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64), (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64))], x: &DVector<f64>) -> Result<(DVector<f64>, DMatrix<f64>), &'static str> {
        let mut d_wl_full = DMatrix::zeros(subset.len(), state.covariance.nrows());
        for (i, (c, ref_sat)) in subset.iter().enumerate() {
            d_wl_full[(i, CORE_STATE_SIZE + c.1)] = 1.0 / c.4;
            d_wl_full[(i, CORE_STATE_SIZE + c.2)] = -1.0 / c.5;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.1)] = -1.0 / ref_sat.4;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.2)] = 1.0 / ref_sat.5;
        }

        let a_wl_full = &d_wl_full * x;
        let q_wl_full = &d_wl_full * &state.covariance * d_wl_full.transpose();
        
        let keep_indices: Vec<usize> = (0..q_wl_full.nrows()).filter(|&i| q_wl_full[(i, i)].sqrt() < 0.15).collect();
        if keep_indices.len() < 4 { return Err("Insufficient well-converged Widelane ambiguities"); }

        let mut d_wl = DMatrix::zeros(keep_indices.len(), state.covariance.nrows());
        for (i, &idx) in keep_indices.iter().enumerate() {
            for j in 0..state.covariance.nrows() { d_wl[(i, j)] = d_wl_full[(idx, j)]; }
        }
        
        let a_wl = &d_wl * x;
        let q_wl = &d_wl * &state.covariance * d_wl.transpose();
        let res_wl = crate::ambiguity::lambda::resolve_lambda(&a_wl, &q_wl).map_err(|_| "WL LAMBDA Failed")?;
        
        if res_wl.ratio < 2.0 || res_wl.success_rate < 0.99 { return Err("WL ratio test failed"); }

        let s_inv = q_wl.try_inverse().ok_or("WL Cov Inversion failed")?;
        let k_wl = &state.covariance * d_wl.transpose() * s_inv;
        let dx_wl = &k_wl * (res_wl.best_integers - a_wl);
        let p_wl = crate::math::covariance::apply_joseph_covariance_update(&state.covariance, &k_wl, &d_wl, &DMatrix::zeros(keep_indices.len(), keep_indices.len()));
        Ok((x + dx_wl, p_wl))
    }"""

start_rwl = orig.find(t_rwl)
end_rwl = orig.find("    fn resolve_narrowlane_ar")
if start_rwl != -1 and end_rwl != -1:
    orig = orig[:start_rwl] + new_rwl + "\n\n" + orig[end_rwl:]

# Now build_measurements
t_bm = """    fn build_measurements(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, iter: usize) -> Vec<FgMeasurement> {"""
new_bm = """    fn build_measurements(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, iter: usize) -> Vec<FgMeasurement> {
        let mut meas = Vec::new();
        let tide_offset = gneiss_core::tides::solid_earth_tides_ecef(state.time, state.position.vector);
        let rcv_pos = Vector3::new(x_i[0], x_i[1], x_i[2]) + tide_offset;
        let ztd = if x_i.len() > 20 && !x_i[20].is_nan() && x_i[20] != 0.0 { x_i[20] } else { state.zwd };

        for sat in sats {
            let geometric_dist = (sat.sat_pos_rot - rcv_pos).norm();
            let dist = geometric_dist - sat.pcv_correction;
            let los = (sat.sat_pos_rot - rcv_pos) / geometric_dist;
            let isb = Self::extract_isb(x_i, sat.sat_obs.sat.constellation);
            let expected_base = dist + x_i[15] + isb - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet;

            if !self.push_sat_meas(&mut meas, state, sat, x_i, iter, &los, expected_base, dist, isb) { continue; }
            if sat.doppler != 0.0 { self.push_doppler_measurement(&mut meas, sat, x_i, &los); }
        }
        meas
    }

    fn push_sat_meas(&self, meas: &mut Vec<FgMeasurement>, state: &RtkState, sat: &ProcessedSat, x_i: &DVector<f64>, iter: usize, los: &Vector3<f64>, expected_base: f64, dist: f64, isb: f64) -> bool {
        let expected_pr = if sat.is_iono_free { expected_base } else if sat.cp1.is_some() && sat.cp2.is_some() && sat.p2.is_some() {
            expected_base + find_amb_idx(state, sat.sat_obs.sat, 3).map(|i| x_i[CORE_STATE_SIZE + i]).unwrap_or(sat.iono_delay)
        } else { expected_base + sat.iono_delay };
        
        let res_pr = sat.p1 - expected_pr;
        if res_pr.abs() > 100.0 { return false; }

        if !sat.is_iono_free && sat.cp1.is_some() && sat.cp2.is_some() && sat.p2.is_some() {
            self.push_uduc_measurements(meas, state, sat, x_i, iter, los, expected_base, dist, isb);
        } else {
            self.push_pr_measurement(meas, state, sat, x_i, iter, los, expected_base, dist, isb);
            if let Some(cp1) = sat.cp1 { if cp1 != 0.0 { self.push_cp_measurement(meas, state, sat, x_i, iter, los, expected_base, dist, cp1); } }
        }
        true
    }"""

start_bm = orig.find(t_bm)
end_bm = orig.find("    fn extract_isb")
if start_bm != -1 and end_bm != -1:
    orig = orig[:start_bm] + new_bm + "\n\n" + orig[end_bm:]

# push_uduc_measurements
t_pum = """    fn push_uduc_measurements(&self, meas: &mut Vec<FgMeasurement>, state: &RtkState, sat: &ProcessedSat, x_i: &DVector<f64>, _iter: usize, los: &Vector3<f64>, expected_base: f64, _dist: f64, _isb: f64) {"""
new_pum = """    fn push_uduc_measurements(&self, meas: &mut Vec<FgMeasurement>, state: &RtkState, sat: &ProcessedSat, x_i: &DVector<f64>, _iter: usize, los: &Vector3<f64>, expected_base: f64, _dist: f64, _isb: f64) {
        let i1_idx = find_amb_idx(state, sat.sat_obs.sat, 3).map(|idx| CORE_STATE_SIZE + idx);
        let n1_idx = find_amb_idx(state, sat.sat_obs.sat, 1).map(|idx| CORE_STATE_SIZE + idx);
        let n2_idx = find_amb_idx(state, sat.sat_obs.sat, 2).map(|idx| CORE_STATE_SIZE + idx);
        let (i1, n1, n2) = (i1_idx.map(|idx| x_i[idx]).unwrap_or(0.0), n1_idx.map(|idx| x_i[idx]).unwrap_or(0.0), n2_idx.map(|idx| x_i[idx]).unwrap_or(0.0));
        let gamma = (sat.f1 * sat.f1) / (sat.f2 * sat.f2);

        let var_p1 = PSEUDORANGE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        let res_p1 = sat.p1 - (expected_base + i1);
        meas.push(FgMeasurement { res: res_p1, h_row: build_h_row_uduc(los, sat.map_wet, i1_idx, 1.0, None, x_i.len(), sat.sat_obs.sat.constellation), weight: var_p1 / apply_huber(res_p1, var_p1, self.huber_k), raw_var: var_p1, is_phase: false, sat: Some(sat.sat_obs.sat) });

        let res_p2 = sat.p2.unwrap() - (expected_base + gamma * i1);
        meas.push(FgMeasurement { res: res_p2, h_row: build_h_row_uduc(los, sat.map_wet, i1_idx, gamma, None, x_i.len(), sat.sat_obs.sat.constellation), weight: (var_p1 * 1.5) / apply_huber(res_p2, var_p1 * 1.5, self.huber_k), raw_var: var_p1 * 1.5, is_phase: false, sat: Some(sat.sat_obs.sat) });

        let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
        let var_l1 = 0.0001 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        
        let res_l1 = (sat.cp1.unwrap() + windup) * sat.lam1 - (expected_base - i1 + n1);
        meas.push(FgMeasurement { res: res_l1, h_row: build_h_row_uduc(los, sat.map_wet, i1_idx, -1.0, n1_idx, x_i.len(), sat.sat_obs.sat.constellation), weight: var_l1 / apply_huber(res_l1, var_l1, self.huber_k), raw_var: var_l1, is_phase: true, sat: Some(sat.sat_obs.sat) });

        let res_l2 = (sat.cp2.unwrap() + windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
        meas.push(FgMeasurement { res: res_l2, h_row: build_h_row_uduc(los, sat.map_wet, i1_idx, -gamma, n2_idx, x_i.len(), sat.sat_obs.sat.constellation), weight: (var_l1 * 1.5) / apply_huber(res_l2, var_l1 * 1.5, self.huber_k), raw_var: var_l1 * 1.5, is_phase: true, sat: Some(sat.sat_obs.sat) });
    }"""

start_pum = orig.find(t_pum)
end_pum = orig.find("}\n\nfn log_ppp_convergence")
if start_pum != -1 and end_pum != -1:
    orig = orig[:start_pum] + new_pum + orig[end_pum:]

# log_ppp_convergence
t_lpc = """fn log_ppp_convergence(state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, x_pred: &DVector<f64>, p_pred: &DMatrix<f64>, solver: &PppFactorGraph) {"""
new_lpc = """fn log_ppp_convergence(state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, x_pred: &DVector<f64>, p_pred: &DMatrix<f64>, solver: &PppFactorGraph) {
    let p_amb = if p_pred.nrows() > 21 { p_pred[(21, 21)] } else { 0.0 };
    tracing::info!("PPP Epoch: pos=[{:.2}, {:.2}, {:.2}], dx_norm={:.4}", x_i[0], x_i[1], x_i[2], (x_i.clone() - x_pred.clone()).norm());
    
    let meas = solver.build_measurements(state, sats, x_i, solver.max_iterations - 1);
    let (mut sum_pr, mut count_pr, mut sum_rr, mut count_rr) = (0.0, 0, 0.0, 0);
    for m in &meas {
        if m.h_row.len() == x_i.len() {
            if m.is_phase { continue; }
            if m.weight > 0.05 { sum_pr += m.res.abs(); count_pr += 1; } else { sum_rr += m.res.abs(); count_rr += 1; }
        }
    }
    if state.epoch_count.is_multiple_of(100) { tracing::trace!("Epoch {}: Mean PR Res = {:.3} m, Mean RR Res = {:.3} m/s", state.epoch_count, sum_pr / count_pr.max(1) as f64, sum_rr / count_rr.max(1) as f64); }
}"""

start_lpc = orig.find(t_lpc)
end_lpc = orig.find("fn build_weight_matrix")
if start_lpc != -1 and end_lpc != -1:
    orig = orig[:start_lpc] + new_lpc + "\n\n" + orig[end_lpc:]

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(orig)
