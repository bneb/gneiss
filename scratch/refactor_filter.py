import re

with open("crates/gneiss-rtk/src/filter.rs", "r") as f:
    content = f.read()

# I will write a regex to find resolve_ambiguities and replace it with the new version and apply_ar_fix
replacement = """    #[allow(clippy::type_complexity)]
    pub fn resolve_ambiguities(&self, ephemerides: &[gneiss_core::ephemeris::Ephemeris], min_subset: usize, ar_min_epoch_count: u32, ar_min_lock: u32, lambda_min_ratio: f64) -> Result<(RtkState, DVector<f64>, DMatrix<f64>, f64, usize), &'static str> {
        let num_amb = self.ambiguities.len();
        if num_amb < min_subset || self.epoch_count <= ar_min_epoch_count as usize { return Err("Insufficient data"); }
        
        let candidate_vars = select_ar_candidates(self, ephemerides, ar_min_lock);
        if candidate_vars.len() < min_subset { return Err("Insufficient candidates"); }
        
        let max_subset = candidate_vars.len().min(24);
        for subset_size in (min_subset..=max_subset).rev() {
            let (d_mat_small, a_cycles, q_cycles) = build_lambda_matrices(self, &candidate_vars, subset_size, ephemerides);
            
            if let Ok(res) = crate::lambda::resolve_lambda(&a_cycles, &q_cycles) {
                let dynamic_threshold = crate::ffrt::calculate_threshold(subset_size, 0.001).max(lambda_min_ratio);
                if res.ratio >= dynamic_threshold {
                    let (fixed_state, da_meters, d_full) = self.apply_ar_fix(subset_size, &candidate_vars, &res, ephemerides)?;
                    return Ok((fixed_state, da_meters, d_full, res.ratio, subset_size));
                }
            }
        }
        Err("AR failed to resolve")
    }

    fn apply_ar_fix(&self, subset_size: usize, candidate_vars: &[(usize, usize, u16, f64)], res: &crate::lambda::LambdaResult, ephemerides: &[gneiss_core::ephemeris::Ephemeris]) -> Result<(RtkState, DVector<f64>, DMatrix<f64>), &'static str> {
        let mut da_meters = DVector::zeros(subset_size);
        let a_sd = nalgebra::DVector::from_vec(self.ambiguities.clone());
        for row in 0..subset_size {
            let (rov, r_idx, _, _) = candidate_vars[row];
            let (rov_sat_id, freq_band) = self.ambiguity_keys[rov];
            let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
            let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
            let lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_band == 1 { f1 } else { f2 };
            let a_cycle_float = (a_sd[rov] - a_sd[r_idx]) / lam;
            da_meters[row] = (res.best_integers[row] - a_cycle_float) * lam;
        }

        let state_size = self.covariance.nrows();
        let mut d_full = DMatrix::zeros(subset_size, state_size);
        for row in 0..subset_size {
            let (rov, r_idx, _, _) = candidate_vars[row];
            d_full[(row, CORE_STATE_SIZE + rov)] = 1.0;
            d_full[(row, CORE_STATE_SIZE + r_idx)] = -1.0;
        }

        let s = &d_full * &self.covariance * d_full.transpose();
        let s_inv = s.try_inverse().ok_or("Fix covariance inversion failed")?;
        let mut k_full = &self.covariance * d_full.transpose() * &s_inv;
        
        for i in 6..15 {
            for j in 0..k_full.ncols() { k_full[(i, j)] = 0.0; }
        }
        
        let dx = &k_full * &da_meters;
        let mut fixed_state = self.clone();
        crate::engine::updater::apply_state_correction(&mut fixed_state, &dx);
        let r_zero = DMatrix::zeros(subset_size, subset_size);
        fixed_state.covariance = crate::engine::updater::apply_joseph_covariance_update(&self.covariance, &k_full, &d_full, &r_zero);
        fixed_state.is_fixed = true;

        Ok((fixed_state, da_meters, d_full))
    }"""

# find resolve_ambiguities function body using regex
pattern = r"    #\[allow\(clippy::type_complexity\)\]\n    pub fn resolve_ambiguities.*?\n    }\n"
# we need to be careful with nested braces, so we will use python parsing or just index search.
import sys
start_idx = content.find("    #[allow(clippy::type_complexity)]\n    pub fn resolve_ambiguities")
if start_idx == -1:
    print("Could not find resolve_ambiguities")
    sys.exit(1)

# Find the end of the resolve_ambiguities function
end_idx = content.find("\n    pub fn prune_stale_ambiguities", start_idx)
if end_idx == -1:
    print("Could not find end of resolve_ambiguities")
    sys.exit(1)

new_content = content[:start_idx] + replacement + "\n" + content[end_idx:]

with open("crates/gneiss-rtk/src/filter.rs", "w") as f:
    f.write(new_content)

print("Refactored resolve_ambiguities in filter.rs")
