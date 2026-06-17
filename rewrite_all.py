import re

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    ppp = f.read()

# Fix process_ppp
ppp = ppp.replace("""pub fn process_ppp<'a>(engine: &'a mut ProcessingEngine, rover_obs: &'a EpochObs) -> Result<&'a RtkState, EngineError> {
    if !valid_pos(engine) {
        return engine.process_spp(rover_obs);
    }
    
    let dt = rover_obs.time.tow - engine.current_state.as_ref().unwrap().time.tow;
    engine.predict_state(dt);
    let state = engine.current_state.as_mut().unwrap();
    state.time = rover_obs.time;
    state.position.epoch = rover_obs.time;

    let sats = build_sats(engine, rover_obs);
    if sats.is_empty() {
        return Err(EngineError::InsufficientSatellites);
    }

    let state = engine.current_state.as_mut().unwrap();
    update_phase_ambiguities(state, &sats, rover_obs.time);
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
    
    let fg = PppFactorGraph::new();
    fg.solve(state, &sats)?;
    
    state.epoch_count += 1;
    
    // We need to clone state BEFORE taking a reference again
    let final_state = engine.current_state.as_ref().unwrap().clone();
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));
    
    Ok(engine.current_state.as_ref().unwrap())
}""", """pub fn process_ppp<'a>(engine: &'a mut ProcessingEngine, rover_obs: &'a EpochObs) -> Result<&'a RtkState, EngineError> {
    if !valid_pos(engine) { return engine.process_spp(rover_obs); }
    let dt = rover_obs.time.tow - engine.current_state.as_ref().unwrap().time.tow;
    engine.predict_state(dt);
    let state = engine.current_state.as_mut().unwrap();
    state.time = rover_obs.time; state.position.epoch = rover_obs.time;
    let sats = build_sats(engine, rover_obs);
    if sats.is_empty() { return Err(EngineError::InsufficientSatellites); }
    let state = engine.current_state.as_mut().unwrap();
    update_phase_ambiguities(state, &sats, rover_obs.time);
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
    PppFactorGraph::new().solve(state, &sats)?;
    state.epoch_count += 1;
    let final_state = engine.current_state.as_ref().unwrap().clone();
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));
    Ok(engine.current_state.as_ref().unwrap())
}""")

# Fix update_phase_ambiguities
ppp = ppp.replace("""fn update_phase_ambiguities(state: &mut RtkState, sats: &[ProcessedSat], t: gneiss_core::time::GpsTime) {
    for sat in sats {
        if let Some(cp1) = sat.cp1 {
            if cp1 == 0.0 { continue; }
            let wup = gneiss_core::windup::phase_windup(sat.sat_pos_rot, gneiss_core::sun::sun_position_ecef(t), sat.rcv_pos_ecef, *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0));
            state.windup.insert(sat.sat_obs.sat, wup);
            
            let l_meas = if sat.is_iono_free && sat.cp2.is_some() { crate::engine::ppp_math::compute_iono_free((cp1 + wup) * sat.lam1, (sat.cp2.unwrap() + wup) * sat.lam2, sat.f1, sat.f2) } else { (cp1 + wup) * sat.lam1 };
            let prev = *state.locktimes.get(&(sat.sat_obs.sat, 1)).unwrap_or(&0);
            let (slip, new_lk) = crate::engine::ppp_math::detect_cycle_slip(sat.sat_obs, prev as u32);
            state.locktimes.insert((sat.sat_obs.sat, 1), new_lk as u16);
            
            if slip {
                for i in 0..4 { state.remove_ambiguity(sat.sat_obs.sat, i); }
            }
            
            let isb = match sat.sat_obs.sat.constellation { Constellation::Glonass => state.isb_glo, Constellation::Galileo => state.isb_gal, Constellation::Beidou => state.isb_bds, _ => 0.0 };
            let expected_base = sat.dist + state.rcv_clk_bias + isb - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet;
            
            if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
                add_uduc_ambiguities(state, sat, cp1, wup, expected_base);
            } else {
                let exp = if sat.is_iono_free && sat.cp2.is_some() { expected_base } else { expected_base - sat.iono_delay };
                if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) { state.add_ambiguity(sat.sat_obs.sat, 0, l_meas - exp, 10000.0); }
                state.last_observed.insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
            }
        }
    }
}""", """fn update_phase_ambiguities(state: &mut RtkState, sats: &[ProcessedSat], t: gneiss_core::time::GpsTime) {
    for sat in sats.iter().filter(|s| s.cp1.unwrap_or(0.0) != 0.0) {
        let cp1 = sat.cp1.unwrap();
        let wup = gneiss_core::windup::phase_windup(sat.sat_pos_rot, gneiss_core::sun::sun_position_ecef(t), sat.rcv_pos_ecef, *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0));
        state.windup.insert(sat.sat_obs.sat, wup);
        let l_meas = if sat.is_iono_free && sat.cp2.is_some() { crate::engine::ppp_math::compute_iono_free((cp1 + wup) * sat.lam1, (sat.cp2.unwrap() + wup) * sat.lam2, sat.f1, sat.f2) } else { (cp1 + wup) * sat.lam1 };
        let prev = *state.locktimes.get(&(sat.sat_obs.sat, 1)).unwrap_or(&0);
        let (slip, new_lk) = crate::engine::ppp_math::detect_cycle_slip(sat.sat_obs, prev as u32);
        state.locktimes.insert((sat.sat_obs.sat, 1), new_lk as u16);
        if slip { for i in 0..4 { state.remove_ambiguity(sat.sat_obs.sat, i); } }
        let isb = match sat.sat_obs.sat.constellation { Constellation::Glonass => state.isb_glo, Constellation::Galileo => state.isb_gal, Constellation::Beidou => state.isb_bds, _ => 0.0 };
        let expected_base = sat.dist + state.rcv_clk_bias + isb - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet;
        if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
            add_uduc_ambiguities(state, sat, cp1, wup, expected_base);
        } else {
            let exp = if sat.is_iono_free && sat.cp2.is_some() { expected_base } else { expected_base - sat.iono_delay };
            if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) { state.add_ambiguity(sat.sat_obs.sat, 0, l_meas - exp, 10000.0); }
            state.last_observed.insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
        }
    }
}""")

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(ppp)

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    fg = f.read()

fg = fg.replace("""    fn resolve_widelane_ar(&self, state: &mut RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>) -> bool {
        let mut wl_meas = Vec::new();
        let mut n_wl = 0;
        
        for sat in sats {
            if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
                if let Some(n1_idx) = find_amb_idx(state, sat.sat_obs.sat, 1) {
                    if let Some(n2_idx) = find_amb_idx(state, sat.sat_obs.sat, 2) {
                        let i_n1 = n1_idx + CORE_STATE_SIZE;
                        let i_n2 = n2_idx + CORE_STATE_SIZE;
                        let n_w = x_i[i_n1] - x_i[i_n2];
                        let var_nw = state.covariance[(i_n1, i_n1)] + state.covariance[(i_n2, i_n2)] - 2.0 * state.covariance[(i_n1, i_n2)];
                        
                        if var_nw < 0.25 {
                            wl_meas.push((sat.sat_obs.sat, n_w, var_nw, n1_idx, n2_idx));
                            n_wl += 1;
                        }
                    }
                }
            }
        }
        
        if n_wl < 4 {
            tracing::info!("Cascade AR did not fix: \"Insufficient well-converged Widelane ambiguities\"");
            return false;
        }
        
        // Single difference
        let mut sd_meas = Vec::new();
        let ref_sat = wl_meas.iter().min_by(|a, b| a.2.partial_cmp(&b.2).unwrap()).unwrap();
        
        for m in &wl_meas {
            if m.0 != ref_sat.0 {
                sd_meas.push((m.0, m.1 - ref_sat.1, m.2 + ref_sat.2, m.3, m.4));
            }
        }
        false
    }""", """    fn resolve_widelane_ar(&self, state: &mut RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>) -> bool {
        let mut wl_meas = Vec::new();
        for sat in sats.iter().filter(|s| !s.is_iono_free && s.cp2.is_some() && s.p2.is_some()) {
            if let (Some(n1_idx), Some(n2_idx)) = (find_amb_idx(state, sat.sat_obs.sat, 1), find_amb_idx(state, sat.sat_obs.sat, 2)) {
                let i_n1 = n1_idx + CORE_STATE_SIZE;
                let i_n2 = n2_idx + CORE_STATE_SIZE;
                let n_w = x_i[i_n1] - x_i[i_n2];
                let var_nw = state.covariance[(i_n1, i_n1)] + state.covariance[(i_n2, i_n2)] - 2.0 * state.covariance[(i_n1, i_n2)];
                if var_nw < 0.25 { wl_meas.push((sat.sat_obs.sat, n_w, var_nw, n1_idx, n2_idx)); }
            }
        }
        if wl_meas.len() < 4 { tracing::info!("Cascade AR did not fix: \"Insufficient well-converged Widelane ambiguities\""); return false; }
        let ref_sat = wl_meas.iter().min_by(|a, b| a.2.partial_cmp(&b.2).unwrap()).unwrap();
        let mut sd_meas = Vec::new();
        for m in &wl_meas { if m.0 != ref_sat.0 { sd_meas.push((m.0, m.1 - ref_sat.1, m.2 + ref_sat.2, m.3, m.4)); } }
        false
    }""")

fg = fg.replace("""    fn push_cp_measurement(&self, meas: &mut Vec<FgMeasurement>, state: &RtkState, sat: &ProcessedSat, x_i: &DVector<f64>, iter: usize, los: &Vector3<f64>, expected_base: f64, _dist: f64, _isb: f64) {
        if sat.is_iono_free {
            if let Some(amb_idx) = find_amb_idx(state, sat.sat_obs.sat, 0) {
                let mut h_row = build_h_row(los, sat.map_wet, Some(amb_idx + CORE_STATE_SIZE), x_i.len(), sat.sat_obs.sat.constellation);
                if x_i.len() > 20 { h_row[20] = sat.map_wet; }
                let cp_res = sat.cp1.unwrap() * sat.lam1 - (expected_base + x_i[amb_idx + CORE_STATE_SIZE]);
                meas.push(FgMeasurement { res: cp_res, h_row, weight: 1.0 / (0.003 * 0.003), raw_var: 0.003 * 0.003, is_phase: true, sat: Some(sat.sat_obs.sat) });
            }
        } else if sat.cp2.is_some() && sat.p2.is_some() {
            let i_idx = find_amb_idx(state, sat.sat_obs.sat, 0);
            let n1_idx = find_amb_idx(state, sat.sat_obs.sat, 1);
            let n2_idx = find_amb_idx(state, sat.sat_obs.sat, 2);
            if i_idx.is_some() && n1_idx.is_some() && n2_idx.is_some() {
                let (i_i, n1_i, n2_i) = (i_idx.unwrap() + CORE_STATE_SIZE, n1_idx.unwrap() + CORE_STATE_SIZE, n2_idx.unwrap() + CORE_STATE_SIZE);
                let wup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
                
                let res1 = (sat.cp1.unwrap() + wup) * sat.lam1 - (expected_base - x_i[i_i] + x_i[n1_i]);
                let h1 = build_h_row_uduc(los, sat.map_wet, Some(i_i), -1.0, Some(n1_i), x_i.len(), sat.sat_obs.sat.constellation);
                meas.push(FgMeasurement { res: res1, h_row: h1, weight: 1.0 / (0.003 * 0.003), raw_var: 0.003 * 0.003, is_phase: true, sat: Some(sat.sat_obs.sat) });
                
                let gamma = (sat.f1 / sat.f2).powi(2);
                let res2 = (sat.cp2.unwrap() + wup) * sat.lam2 - (expected_base - gamma * x_i[i_i] + x_i[n2_i]);
                let h2 = build_h_row_uduc(los, sat.map_wet, Some(i_i), -gamma, Some(n2_i), x_i.len(), sat.sat_obs.sat.constellation);
                meas.push(FgMeasurement { res: res2, h_row: h2, weight: 1.0 / (0.003 * 0.003), raw_var: 0.003 * 0.003, is_phase: true, sat: Some(sat.sat_obs.sat) });
            }
        }
    }""", """    fn push_cp_measurement(&self, meas: &mut Vec<FgMeasurement>, state: &RtkState, sat: &ProcessedSat, x_i: &DVector<f64>, iter: usize, los: &Vector3<f64>, expected_base: f64, _dist: f64, _isb: f64) {
        if sat.is_iono_free {
            if let Some(a_i) = find_amb_idx(state, sat.sat_obs.sat, 0) {
                let mut h_row = build_h_row(los, sat.map_wet, Some(a_i + CORE_STATE_SIZE), x_i.len(), sat.sat_obs.sat.constellation);
                if x_i.len() > 20 { h_row[20] = sat.map_wet; }
                meas.push(FgMeasurement { res: sat.cp1.unwrap() * sat.lam1 - (expected_base + x_i[a_i + CORE_STATE_SIZE]), h_row, weight: 1.0 / 0.000009, raw_var: 0.000009, is_phase: true, sat: Some(sat.sat_obs.sat) });
            }
        } else if sat.cp2.is_some() && sat.p2.is_some() {
            if let (Some(i_idx), Some(n1_idx), Some(n2_idx)) = (find_amb_idx(state, sat.sat_obs.sat, 0), find_amb_idx(state, sat.sat_obs.sat, 1), find_amb_idx(state, sat.sat_obs.sat, 2)) {
                let (i_i, n1_i, n2_i) = (i_idx + CORE_STATE_SIZE, n1_idx + CORE_STATE_SIZE, n2_idx + CORE_STATE_SIZE);
                let wup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
                let h1 = build_h_row_uduc(los, sat.map_wet, Some(i_i), -1.0, Some(n1_i), x_i.len(), sat.sat_obs.sat.constellation);
                meas.push(FgMeasurement { res: (sat.cp1.unwrap() + wup) * sat.lam1 - (expected_base - x_i[i_i] + x_i[n1_i]), h_row: h1, weight: 1.0 / 0.000009, raw_var: 0.000009, is_phase: true, sat: Some(sat.sat_obs.sat) });
                let gamma = (sat.f1 / sat.f2).powi(2);
                let h2 = build_h_row_uduc(los, sat.map_wet, Some(i_i), -gamma, Some(n2_i), x_i.len(), sat.sat_obs.sat.constellation);
                meas.push(FgMeasurement { res: (sat.cp2.unwrap() + wup) * sat.lam2 - (expected_base - gamma * x_i[i_i] + x_i[n2_i]), h_row: h2, weight: 1.0 / 0.000009, raw_var: 0.000009, is_phase: true, sat: Some(sat.sat_obs.sat) });
            }
        }
    }""")

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(fg)

with open("crates/gneiss-rtk/src/engine/ppp_math.rs", "r") as f:
    math = f.read()

math = math.replace("""pub fn detect_cycle_slip(sat_obs: &SatObs, prev_locktime: u32) -> (bool, u32) {
    let mut slip = false;
    let mut new_lk = prev_locktime.saturating_add(1);

    for obs in &sat_obs.observations {
        if obs.code.obs_type == ObsType::CarrierPhase {
            if let Some(lli) = obs.lli {
                if lli & 1 != 0 {
                    slip = true;
                    new_lk = 1;
                }
            }
        }
    }
    
    (slip, new_lk)
}""", """pub fn detect_cycle_slip(sat_obs: &SatObs, prev_locktime: u32) -> (bool, u32) {
    let slip = sat_obs.observations.iter().any(|o| o.code.obs_type == ObsType::CarrierPhase && o.lli.unwrap_or(0) & 1 != 0);
    (slip, if slip { 1 } else { prev_locktime.saturating_add(1) })
}""")

with open("crates/gneiss-rtk/src/engine/ppp_math.rs", "w") as f:
    f.write(math)
