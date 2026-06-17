import re

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    orig = f.read()

# I will define new versions of solve, resolve_widelane_ar, build_measurements, push_uduc_measurements, log_ppp_convergence

# For solve:
solve_str = """    pub fn solve(&mut self, state: &mut RtkState, sats: &[ProcessedSat]) -> Result<(), EngineError> {
        if sats.is_empty() { return Ok(()); }
        let num_sats = sats.len();
        self.num_sats_used = num_sats;
        
        for _ in 0..self.config.max_iterations {
            self.iteration += 1;
            let (dx, x_vec, ambs) = self.compute_iteration_dx(state, sats)?;
            self.apply_state_vector(state, &dx, &ambs);
            
            if dx.norm() < self.config.convergence_threshold {
                self.converged = true;
                break;
            }
        }
        
        let (_, x_vec, ambs) = self.compute_iteration_dx(state, sats)?;
        if self.config.uduc_ar { self.resolve_cascade_ar(state, sats, &x_vec, &ambs)?; }
        self.compute_final_covariance(state, sats);
        self.log_ppp_convergence(state, sats);
        Ok(())
    }"""

# For log_ppp_convergence
log_str = """    fn log_ppp_convergence(&self, state: &RtkState, sats: &[ProcessedSat]) {
        let (mut pr, mut cp) = (0, 0);
        for sat in sats {
            if sat.p1.is_some() { pr += 1; }
            if sat.cp1.is_some() { cp += 1; }
        }
        tracing::info!("PPP converged: {}, iter: {}, sats: {}, PR: {}, CP: {}", self.converged, self.iteration, sats.len(), pr, cp);
        
        let (mut fx, mut fl) = (0, 0);
        for (_, lk) in &state.locktimes { if *lk > 10 { fx += 1; } }
        for (_, lk) in &state.locktimes { if *lk > 0 { fl += 1; } }
        tracing::debug!("PPP ambiguities: {} fixed, {} floating", fx, fl);
        
        let mut pdop = 0.0;
        let mut hdop = 0.0;
        let mut vdop = 0.0;
        if let Some(cov) = &state.covariance {
            pdop = (cov[(0, 0)] + cov[(1, 1)] + cov[(2, 2)]).sqrt();
            hdop = (cov[(0, 0)] + cov[(1, 1)]).sqrt();
            vdop = cov[(2, 2)].sqrt();
        }
        tracing::debug!("PPP DOPs: PDOP={:.2}, HDOP={:.2}, VDOP={:.2}", pdop, hdop, vdop);
    }"""

# For push_uduc_measurements
push_uduc_str = """    fn push_uduc_measurements(&mut self, state: &RtkState, sat: &ProcessedSat, pcv: f64) {
        if let Some(p1) = sat.p1 {
            let i1 = state.ambiguities.get(&(sat.sat_obs.sat, 3)).map(|a| a.value).unwrap_or(0.0);
            let inn = p1 + pcv - (sat.dist + state.rcv_clk_bias - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet - i1);
            self.meas.push(inn);
            self.weights.push(crate::engine::ppp_fg::snr_scale(sat.snr, sat.el) / 100.0);
            self.h_rows.push(crate::engine::ppp_fg::build_h_row_uduc(state, sat, 3, 0.0, -1.0));
        }
        if let Some(p2) = sat.p2 {
            let i1 = state.ambiguities.get(&(sat.sat_obs.sat, 3)).map(|a| a.value).unwrap_or(0.0);
            let gamma = (sat.f1 * sat.f1) / (sat.f2 * sat.f2);
            let inn = p2 + pcv - (sat.dist + state.rcv_clk_bias - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet - gamma * i1);
            self.meas.push(inn);
            self.weights.push(crate::engine::ppp_fg::snr_scale(sat.snr, sat.el) / 100.0);
            self.h_rows.push(crate::engine::ppp_fg::build_h_row_uduc(state, sat, 3, 0.0, -gamma));
        }
        if let Some(cp1) = sat.cp1 {
            let n1 = state.ambiguities.get(&(sat.sat_obs.sat, 1)).map(|a| a.value).unwrap_or(0.0);
            let i1 = state.ambiguities.get(&(sat.sat_obs.sat, 3)).map(|a| a.value).unwrap_or(0.0);
            let inn = cp1 * sat.lam1 + pcv - (sat.dist + state.rcv_clk_bias - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet - i1 + n1);
            self.meas.push(inn);
            self.weights.push(crate::engine::ppp_fg::snr_scale(sat.snr, sat.el) / 0.01);
            self.h_rows.push(crate::engine::ppp_fg::build_h_row_uduc(state, sat, 1, 1.0, -1.0));
        }
        if let Some(cp2) = sat.cp2 {
            let n2 = state.ambiguities.get(&(sat.sat_obs.sat, 2)).map(|a| a.value).unwrap_or(0.0);
            let i1 = state.ambiguities.get(&(sat.sat_obs.sat, 3)).map(|a| a.value).unwrap_or(0.0);
            let gamma = (sat.f1 * sat.f1) / (sat.f2 * sat.f2);
            let inn = cp2 * sat.lam2 + pcv - (sat.dist + state.rcv_clk_bias - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet - gamma * i1 + n2);
            self.meas.push(inn);
            self.weights.push(crate::engine::ppp_fg::snr_scale(sat.snr, sat.el) / 0.01);
            self.h_rows.push(crate::engine::ppp_fg::build_h_row_uduc(state, sat, 2, 1.0, -gamma));
        }
    }"""

