import sys

with open('crates/gneiss-rtk/src/filter.rs', 'r') as f:
    content = f.read()

resolve_start = content.find('    pub fn resolve_ambiguities(&self, ephemerides: &[gneiss_core::ephemeris::Ephemeris]')
prune_start = content.find('    pub fn prune_stale_ambiguities(&mut self, current_epoch: u32, threshold: u32) {')

if resolve_start == -1 or prune_start == -1:
    print("Could not find start indices")
    sys.exit(1)

pre_resolve = content[:resolve_start]
post_resolve = content[prune_start:]

new_resolve = """    pub fn resolve_ambiguities(&self, ephemerides: &[gneiss_core::ephemeris::Ephemeris], min_subset: usize, ar_min_epoch_count: u32, ar_min_lock: u32, lambda_min_ratio: f64) -> Result<(RtkState, DVector<f64>, DMatrix<f64>, f64, usize), &'static str> {
        let num_amb = self.ambiguities.len();
        tracing::debug!("resolve_ambiguities check: num_amb={} (min_subset={}), epoch_count={} (ar_min_epoch={}), ar_min_lock={}", num_amb, min_subset, self.epoch_count, ar_min_epoch_count, ar_min_lock);
        if num_amb < min_subset || self.epoch_count <= ar_min_epoch_count as usize { return Err("Insufficient data"); }
        
        let mut candidate_vars = self.find_ar_candidates(ephemerides, ar_min_lock);
        candidate_vars.retain(|c| c.3 < 10000.0);
        tracing::debug!("AR candidates: retained: {}", candidate_vars.len());
        
        if candidate_vars.len() >= min_subset {
            let max_subset = candidate_vars.len().min(24);
            let a_sd = nalgebra::DVector::from_vec(self.ambiguities.clone());
            let q_sd = self.covariance.view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb)).into_owned();
            
            for subset_size in (min_subset..=max_subset).rev() {
                if let Some(res) = self.try_resolve_subset(subset_size, &candidate_vars, ephemerides, &a_sd, &q_sd, lambda_min_ratio) {
                    return Ok(res);
                }
            }
        }
        Err("AR failed to resolve")
    }

    fn find_ar_candidates(&self, ephemerides: &[gneiss_core::ephemeris::Ephemeris], ar_min_lock: u32) -> Vec<(usize, usize, u16, f64)> {
        let num_amb = self.ambiguities.len();
        let q_sd = self.covariance.view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb)).into_owned();
        use gneiss_core::sat::Constellation;
        let constellations = [Constellation::Gps, Constellation::Galileo];
        let mut candidates = Vec::new();

        for &constell in &constellations {
            let mut best_ref_idx = None;
            let mut max_lock = 0;
            
            for i in 0..num_amb {
                let (sat, freq) = self.ambiguity_keys[i];
                if sat.constellation != constell || freq != 1 { continue; }
                let lock = *self.locktimes.get(&(sat, freq)).unwrap_or(&0);
                if lock >= ar_min_lock as u16 && lock > max_lock {
                    max_lock = lock;
                    best_ref_idx = Some(i);
                }
            }

            if let Some(ref_idx) = best_ref_idx {
                let ref_sat_id = self.ambiguity_keys[ref_idx].0;
                let l2_ref_idx = self.ambiguity_keys.iter().position(|&(s, f)| s == ref_sat_id && f == 2);

                for i in 0..num_amb {
                    if i == ref_idx || Some(i) == l2_ref_idx { continue; }
                    let (rov_sat, freq) = self.ambiguity_keys[i];
                    if rov_sat.constellation != constell { continue; }
                    let lock = *self.locktimes.get(&(rov_sat, freq)).unwrap_or(&0);
                    
                    if lock >= ar_min_lock as u16 {
                        if freq == 1 {
                            candidates.push((i, ref_idx, lock));
                        } else if let Some(r2_idx) = l2_ref_idx {
                            candidates.push((i, r2_idx, lock));
                        }
                    }
                }
            }
        }
            
        let mut candidate_vars = Vec::new();
        for &(rov, r_idx, lock) in &candidates {
            let (rov_sat_id, freq_band) = self.ambiguity_keys[rov];
            let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
            let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
            let lam = 299792458.0 / if freq_band == 1 { f1 } else { f2 };
            let q_dd = q_sd[(rov, rov)] + q_sd[(r_idx, r_idx)] - 2.0 * q_sd[(rov, r_idx)];
            let var_cycles = q_dd / (lam * lam);
            candidate_vars.push((rov, r_idx, lock, var_cycles));
        }
        
        candidate_vars.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
        candidate_vars
    }

    fn try_resolve_subset(
        &self,
        subset_size: usize,
        candidate_vars: &[(usize, usize, u16, f64)],
        ephemerides: &[gneiss_core::ephemeris::Ephemeris],
        a_sd: &DVector<f64>,
        q_sd: &DMatrix<f64>,
        lambda_min_ratio: f64
    ) -> Option<(RtkState, DVector<f64>, DMatrix<f64>, f64, usize)> {
        let num_amb = self.ambiguities.len();
        let mut d_mat = DMatrix::zeros(subset_size, num_amb);
        let mut a_cycles = DVector::zeros(subset_size);
        
        for row in 0..subset_size {
            let (rov, r_idx, _, _) = candidate_vars[row];
            let (rov_sat_id, freq_band) = self.ambiguity_keys[rov];
            let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
            let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
            let lam = 299792458.0 / if freq_band == 1 { f1 } else { f2 };

            d_mat[(row, rov)] = 1.0;
            d_mat[(row, r_idx)] = -1.0;
            a_cycles[row] = (a_sd[rov] - a_sd[r_idx]) / lam;
        }
        
        let mut q_cycles = DMatrix::zeros(subset_size, subset_size);
        for r in 0..subset_size {
            for c in 0..subset_size {
                let (rov_r, ref_r, _, _) = candidate_vars[r];
                let (rov_c, ref_c, _, _) = candidate_vars[c];
                
                let freq_r = self.ambiguity_keys[rov_r].1;
                let freq_c = self.ambiguity_keys[rov_c].1;
                
                let sat_r = self.ambiguity_keys[rov_r].0;
                let freq_num_r = ephemerides.iter().find(|e| e.sat() == sat_r).map(|e| e.freq_num()).unwrap_or(0);
                let (f1_r, f2_r) = gneiss_core::signal::satellite_frequencies(sat_r, freq_num_r);
                let lam_r = 299792458.0 / if freq_r == 1 { f1_r } else { f2_r };
                
                let sat_c = self.ambiguity_keys[rov_c].0;
                let freq_num_c = ephemerides.iter().find(|e| e.sat() == sat_c).map(|e| e.freq_num()).unwrap_or(0);
                let (f1_c, f2_c) = gneiss_core::signal::satellite_frequencies(sat_c, freq_num_c);
                let lam_c = 299792458.0 / if freq_c == 1 { f1_c } else { f2_c };
                
                let q_dd = q_sd[(rov_r, rov_c)] - q_sd[(rov_r, ref_c)] - q_sd[(ref_r, rov_c)] + q_sd[(ref_r, ref_c)];
                q_cycles[(r, c)] = q_dd / (lam_r * lam_c);
            }
        }
        
        if let Ok(res) = crate::lambda::resolve_lambda(&a_cycles, &q_cycles) {
            let dynamic_threshold = crate::ffrt::calculate_threshold(subset_size, 0.001).max(lambda_min_ratio);
            if res.ratio >= dynamic_threshold {
                tracing::info!("Multi-Const PAR Fixed! Ratio={:.2}, Ps={:.4}", res.ratio, res.success_rate);
                return self.apply_ar_correction(subset_size, candidate_vars, ephemerides, &a_cycles, &res.best_integers, res.ratio);
            }
        }
        None
    }

    fn apply_ar_correction(
        &self,
        subset_size: usize,
        candidate_vars: &[(usize, usize, u16, f64)],
        ephemerides: &[gneiss_core::ephemeris::Ephemeris],
        a_cycles: &DVector<f64>,
        best_integers: &DVector<f64>,
        ratio: f64
    ) -> Option<(RtkState, DVector<f64>, DMatrix<f64>, f64, usize)> {
        let mut da_meters = DVector::zeros(subset_size);
        for row in 0..subset_size {
            let (rov, _, _, _) = candidate_vars[row];
            let (rov_sat_id, freq_band) = self.ambiguity_keys[rov];
            let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
            let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
            let lam = 299792458.0 / if freq_band == 1 { f1 } else { f2 };
            da_meters[row] = (best_integers[row] - a_cycles[row]) * lam;
        }

        let state_size = self.covariance.nrows();
        let mut d_full = DMatrix::zeros(subset_size, state_size);
        for row in 0..subset_size {
            let (rov, r_idx, _, _) = candidate_vars[row];
            d_full[(row, CORE_STATE_SIZE + rov)] = 1.0;
            d_full[(row, CORE_STATE_SIZE + r_idx)] = -1.0;
        }

        let s = &d_full * &self.covariance * d_full.transpose();
        let s_inv = s.try_inverse()?;
        let mut k_full = &self.covariance * d_full.transpose() * &s_inv;
        
        if state_size > 15 {
            for i in 6..15 {
                for j in 0..k_full.ncols() {
                    k_full[(i, j)] = 0.0;
                }
            }
        }
        
        let dx = &k_full * &da_meters;
        let mut fixed_state = self.clone();
        fixed_state.fixed_state = None;
        fixed_state.position.vector += dx.rows(0, 3).into_owned();
        fixed_state.velocity += dx.rows(3, 3).into_owned();
        
        let d_psi = nalgebra::Vector3::new(dx[6], dx[7], dx[8]);
        let angle = d_psi.norm();
        if angle > 1e-12 {
            let dq = nalgebra::UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(d_psi / angle), angle);
            fixed_state.attitude = dq * fixed_state.attitude;
            fixed_state.attitude.renormalize();
        }

        fixed_state.accel_bias.x += dx[9];
        fixed_state.accel_bias.y += dx[10];
        fixed_state.accel_bias.z += dx[11];
        fixed_state.gyro_bias.x += dx[12];
        fixed_state.gyro_bias.y += dx[13];
        fixed_state.gyro_bias.z += dx[14];
        
        if CORE_STATE_SIZE > 15 {
            fixed_state.rcv_clk_bias += dx[15];
            fixed_state.rcv_clk_drift += dx[16];
            fixed_state.zwd += dx[17];
            fixed_state.zwd = fixed_state.zwd.max(0.0);
        }

        for i in 0..self.ambiguities.len() {
            fixed_state.ambiguities[i] += dx[CORE_STATE_SIZE + i];
        }

        let identity = DMatrix::identity(state_size, state_size);
        let i_kd = identity - &k_full * &d_full;
        let mut p_fixed = &i_kd * &self.covariance * i_kd.transpose();
        for r in 0..p_fixed.nrows() {
            for c in 0..r {
                let avg = (p_fixed[(r, c)] + p_fixed[(c, r)]) * 0.5;
                p_fixed[(r, c)] = avg;
                p_fixed[(c, r)] = avg;
            }
        }
        fixed_state.covariance = p_fixed;
        fixed_state.is_fixed = true;

        Some((fixed_state, da_meters, d_full, ratio, subset_size))
    }
"""

with open('crates/gneiss-rtk/src/filter.rs', 'w') as f:
    f.write(pre_resolve + new_resolve + post_resolve)

print("Successfully refactored resolve_ambiguities")
