use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};
use gneiss_core::sat::SatelliteId;
use gneiss_core::coords::Coordinate;

use gneiss_core::time::GpsTime;

pub const CORE_STATE_SIZE: usize = 18;

/// Represents the state of the RTK Extended Kalman Filter (EKF).
#[derive(Debug, Clone)]
pub struct RtkState {
    pub time: GpsTime,
    pub position: Coordinate,
    pub velocity: Vector3<f64>,
    pub attitude: UnitQuaternion<f64>,
    pub accel_bias: Vector3<f64>,
    pub gyro_bias: Vector3<f64>,
    pub rcv_clk_bias: f64,
    pub rcv_clk_drift: f64,
    pub zwd: f64,
    pub gf_values: std::collections::HashMap<SatelliteId, f64>,
    
    pub ambiguities: Vec<f64>,
    pub ambiguity_keys: Vec<(SatelliteId, u8)>, // (Sat, Freq Band 1 or 2)
    pub mw_sd_ema: std::collections::HashMap<SatelliteId, f64>,
    pub mw_sd_counts: std::collections::HashMap<SatelliteId, usize>,
    pub locktimes: std::collections::HashMap<(SatelliteId, u8), u16>,
    pub last_observed: std::collections::HashMap<(SatelliteId, u8), u32>,
    pub windup: std::collections::HashMap<SatelliteId, f64>,
    pub innovation_cov: std::collections::HashMap<(SatelliteId, u8), f64>, // For IAE
    pub innovation_counts: std::collections::HashMap<(SatelliteId, u8), usize>,
    pub reject_counts: std::collections::HashMap<(SatelliteId, u8), usize>,
    pub is_fixed: bool,
    pub epoch_count: usize,
    pub covariance: DMatrix<f64>,
    
    // RTS Smoother matrices (15x15 core blocks)
    pub core_phi: Option<DMatrix<f64>>,       // State transition from k-1 to k
    pub core_p_predict: Option<DMatrix<f64>>, // Predicted covariance P_{k|k-1}
}

impl RtkState {
    pub fn new(time: GpsTime, initial_pos: Coordinate, initial_var: f64) -> Self {
        let mut cov = DMatrix::zeros(CORE_STATE_SIZE, CORE_STATE_SIZE);
        for i in 0..3 { cov[(i, i)] = initial_var; }
        for i in 3..6 { cov[(i, i)] = 100.0; }
        let att_var = (1.0f64.to_radians()).powi(2);
        for i in 6..9 { cov[(i, i)] = att_var; }
        for i in 9..12 { cov[(i, i)] = 0.01; }
        for i in 12..CORE_STATE_SIZE { cov[(i, i)] = 1e-4; }

        Self {
            time,
            position: initial_pos,
            velocity: Vector3::zeros(),
            attitude: UnitQuaternion::identity(),
            accel_bias: Vector3::zeros(),
            gyro_bias: Vector3::zeros(),
            rcv_clk_bias: 0.0,
            rcv_clk_drift: 0.0,
            zwd: 0.1,
            gf_values: std::collections::HashMap::new(),
            ambiguities: Vec::new(),
            ambiguity_keys: Vec::new(),
            mw_sd_ema: std::collections::HashMap::new(),
            mw_sd_counts: std::collections::HashMap::new(),
            locktimes: std::collections::HashMap::new(),
            last_observed: std::collections::HashMap::new(),
            windup: std::collections::HashMap::new(),
            innovation_cov: std::collections::HashMap::new(),
            innovation_counts: std::collections::HashMap::new(),
            reject_counts: std::collections::HashMap::new(),
            is_fixed: false,
            epoch_count: 0,
            covariance: cov,
            core_phi: None,
            core_p_predict: None,
        }
    }

    pub fn update_mw(&mut self, sat: SatelliteId, mw_cycles: f64) {
        let count = self.mw_sd_counts.entry(sat).or_insert(0);
        let ema = self.mw_sd_ema.entry(sat).or_insert(mw_cycles);
        let alpha = 1.0 / ((*count + 1) as f64).min(100.0);
        *ema = *ema * (1.0 - alpha) + mw_cycles * alpha;
        *count += 1;
    }

    pub fn resolve_ambiguities(&mut self, ephemerides: &[gneiss_core::ephemeris::Ephemeris], min_subset: usize) -> Result<(), &'static str> {
        self.is_fixed = false;

        let num_amb = self.ambiguities.len();
        if num_amb < 6 || self.epoch_count <= 50 { return Ok(()); }
        
        let a_sd = nalgebra::DVector::from_vec(self.ambiguities.clone());
        let q_sd = self.covariance.view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb)).into_owned();
        
        use gneiss_core::sat::Constellation;
        let constellations = [Constellation::Gps, Constellation::Galileo, Constellation::Beidou, Constellation::Qzss];
        
        let mut candidates = Vec::new();

        for &constell in &constellations {
            let mut best_ref_idx = None;
            let mut max_lock = 0;
            
            // Find best reference for this constellation (L1 freq only)
            for i in 0..num_amb {
                let (sat, freq) = self.ambiguity_keys[i];
                if sat.constellation != constell || freq != 1 { continue; }
                let lock = *self.locktimes.get(&(sat, freq)).unwrap_or(&0);
                if lock > max_lock {
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
                    
                    if lock >= 5 {
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
            
            // GLONASS uses FDMA (different wavelengths per satellite).
            // Double Differencing across different wavelengths destroys the integer property!
            if rov_sat_id.constellation == gneiss_core::sat::Constellation::Glonass {
                continue;
            }

            let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
            let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
            let lam = 299792458.0 / if freq_band == 1 { f1 } else { f2 };
            
            let q_dd = q_sd[(rov, rov)] + q_sd[(r_idx, r_idx)] - 2.0 * q_sd[(rov, r_idx)];
            let var_cycles = q_dd / (lam * lam);
            
            candidate_vars.push((rov, r_idx, lock, var_cycles));
        }
        
        // Sort primarily by variance in cycles (ascending - lower is better), fallback to lock time
        candidate_vars.sort_by(|a, b| {
            a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal)
        });

        // Prevent LAMBDA search from exploding by only passing highly converged ambiguities
        candidate_vars.retain(|c| c.3 < 1.0);

        if candidate_vars.len() >= min_subset {
            let max_subset = candidate_vars.len().min(24);
            
            for subset_size in (min_subset..=max_subset).rev() {
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
                        let (rov_r, _, _, _) = candidate_vars[r];
                        let (rov_c, _, _, _) = candidate_vars[c];
                        let freq_r = self.ambiguity_keys[rov_r].1;
                        let freq_c = self.ambiguity_keys[rov_c].1;
                        
                        let (sat_r, _) = self.ambiguity_keys[rov_r];
                        let freq_num_r = ephemerides.iter().find(|e| e.sat() == sat_r).map(|e| e.freq_num()).unwrap_or(0);
                        let (f1_r, f2_r) = gneiss_core::signal::satellite_frequencies(sat_r, freq_num_r);
                        let lam_r = 299792458.0 / if freq_r == 1 { f1_r } else { f2_r };
                        
                        let (sat_c, _) = self.ambiguity_keys[rov_c];
                        let freq_num_c = ephemerides.iter().find(|e| e.sat() == sat_c).map(|e| e.freq_num()).unwrap_or(0);
                        let (f1_c, f2_c) = gneiss_core::signal::satellite_frequencies(sat_c, freq_num_c);
                        let lam_c = 299792458.0 / if freq_c == 1 { f1_c } else { f2_c };
                        
                        let (rov_r, ref_r, _, _) = candidate_vars[r];
                        let (rov_c, ref_c, _, _) = candidate_vars[c];
                        let q_dd = q_sd[(rov_r, rov_c)] - q_sd[(rov_r, ref_c)] - q_sd[(ref_r, rov_c)] + q_sd[(ref_r, ref_c)];
                        q_cycles[(r, c)] = q_dd / (lam_r * lam_c);
                    }
                }
                
                // Removed debug print
                if let Ok(res) = crate::lambda::resolve_lambda(&a_cycles, &q_cycles) {
                    // Removed debug print
                    // Dynamic Ratio Threshold (FF-RT)
                    // We rely purely on LAMBDA Ratio Test and FF-RT for validation
                    let dynamic_threshold = crate::ffrt::calculate_threshold(subset_size, 0.001);
                    if res.ratio >= dynamic_threshold {
                        tracing::info!("Multi-Const PAR Fixed & Validated (subset {}/{})! Ratio={:.2}, Ps={:.4}", subset_size, max_subset, res.ratio, res.success_rate);

                        let mut da_meters = DVector::zeros(subset_size);
                        for row in 0..subset_size {
                            let (rov, _, _, _) = candidate_vars[row];
                            let (rov_sat_id, freq_band) = self.ambiguity_keys[rov];
                            let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
                            let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
                            let lam = 299792458.0 / if freq_band == 1 { f1 } else { f2 };
                            da_meters[row] = (res.best_integers[row] - a_cycles[row]) * lam;
                        }

                        let state_size = self.covariance.nrows();
                        let mut d_full = DMatrix::zeros(subset_size, state_size);
                        for row in 0..subset_size {
                            let (rov, r_idx, _, _) = candidate_vars[row];
                            d_full[(row, CORE_STATE_SIZE + rov)] = 1.0;
                            d_full[(row, CORE_STATE_SIZE + r_idx)] = -1.0;
                        }

                        let s = &d_full * &self.covariance * d_full.transpose();
                        let s_inv = s.clone().try_inverse().ok_or("Fix covariance inversion failed")?;
                        let k_full = &self.covariance * d_full.transpose() * &s_inv;
                        
                        let dx = &k_full * da_meters;
                        let mut fixed_state = self.clone();
                        fixed_state.position.vector += dx.rows(0, 3).into_owned();
                        fixed_state.velocity += dx.rows(3, 3).into_owned();
                        
                        for i in 0..self.ambiguities.len() {
                            fixed_state.ambiguities[i] += dx[CORE_STATE_SIZE + i];
                        }

                        // Joseph form update with R=0 to maintain positive semi-definiteness
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
                        
                        // Re-apply update to self
                        self.position = fixed_state.position;
                        self.velocity = fixed_state.velocity;
                        self.ambiguities = fixed_state.ambiguities;
                        self.covariance = fixed_state.covariance.clone();
                        self.is_fixed = true;

                        return Ok(());
                    }
                }
            }
        }

        Ok(())
    }



    pub fn prune_stale_ambiguities(&mut self, current_epoch: u32, threshold: u32) {
        let mut to_remove = Vec::new();
        for key in &self.ambiguity_keys {
            let last = *self.last_observed.get(key).unwrap_or(&0);
            if current_epoch > last && current_epoch - last > threshold {
                to_remove.push(*key);
            }
        }
        for (sat, freq) in to_remove {
            tracing::debug!("Pruning stale ambiguity for {:?} freq {}", sat, freq);
            self.remove_ambiguity(sat, freq);
            self.last_observed.remove(&(sat, freq));
            self.locktimes.remove(&(sat, freq));
        }
    }
    pub fn add_ambiguity(&mut self, sat: SatelliteId, freq: u8, initial_estimate: f64, initial_variance: f64) {
        self.ambiguities.push(initial_estimate);
        self.ambiguity_keys.push((sat, freq));

        let n_old = self.covariance.nrows();
        let n_new = n_old + 1;
        let mut new_cov = DMatrix::zeros(n_new, n_new);
        new_cov.view_mut((0, 0), (n_old, n_old)).copy_from(&self.covariance);
        new_cov[(n_old, n_old)] = initial_variance;
        self.covariance = new_cov;
    }

    pub fn remove_ambiguity(&mut self, sat: SatelliteId, freq: u8) {
        if let Some(idx) = self.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == freq) {
            self.ambiguities.remove(idx);
            self.ambiguity_keys.remove(idx);
            
            let n_old = self.covariance.nrows();
            let n_new = n_old - 1;
            let mut new_cov = DMatrix::zeros(n_new, n_new);
            let cov_idx = CORE_STATE_SIZE + idx;

            if cov_idx > 0 { new_cov.view_mut((0, 0), (cov_idx, cov_idx)).copy_from(&self.covariance.view((0, 0), (cov_idx, cov_idx))); }
            if cov_idx < n_new {
                new_cov.view_mut((0, cov_idx), (cov_idx, n_new - cov_idx)).copy_from(&self.covariance.view((0, cov_idx + 1), (cov_idx, n_old - cov_idx - 1)));
                new_cov.view_mut((cov_idx, 0), (n_new - cov_idx, cov_idx)).copy_from(&self.covariance.view((cov_idx + 1, 0), (n_old - cov_idx - 1, cov_idx)));
                new_cov.view_mut((cov_idx, cov_idx), (n_new - cov_idx, n_new - cov_idx)).copy_from(&self.covariance.view((cov_idx + 1, cov_idx + 1), (n_old - cov_idx - 1, n_old - cov_idx - 1)));
            }
            self.covariance = new_cov;
        }
    }
}

#[derive(Debug, Clone)]
pub struct DdObservation {
    pub sat: SatelliteId,
    pub pr_l1: f64, pub pr_l2: Option<f64>,
    pub cp_l1: Option<f64>, pub cp_l2: Option<f64>,
    pub doppler: f64, pub snr: f64, pub locktime: Option<u16>,
}

pub fn compute_if_combination(v1: f64, v2: Option<f64>, f1: f64, f2: f64) -> f64 {
    if let Some(v2_val) = v2 {
        let f1_2 = f1 * f1;
        let f2_2 = f2 * f2;
        (f1_2 * v1 - f2_2 * v2_val) / (f1_2 - f2_2)
    } else {
        v1
    }
}

pub fn compute_double_difference(
    rover_ref: &DdObservation, rover_sat: &DdObservation, 
    base_ref: &DdObservation, base_sat: &DdObservation,
    ref_f1: f64, ref_f2: f64, sat_f1: f64, sat_f2: f64
) -> f64 {
    let rov_ref = compute_if_combination(rover_ref.pr_l1, rover_ref.pr_l2, ref_f1, ref_f2);
    let rov_sat = compute_if_combination(rover_sat.pr_l1, rover_sat.pr_l2, sat_f1, sat_f2);
    let bas_ref = compute_if_combination(base_ref.pr_l1, base_ref.pr_l2, ref_f1, ref_f2);
    let bas_sat = compute_if_combination(base_sat.pr_l1, base_sat.pr_l2, sat_f1, sat_f2);
    (rov_sat - rov_ref) - (bas_sat - bas_ref)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;
    use gneiss_core::coords::{Datum, Frame};

    #[test]
    fn test_resolve_ambiguities_multi_constellation() {
        let time = GpsTime::new(2137, 422922.0);
        let initial_pos = Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, initial_pos, 10.0);
        state.epoch_count = 100;

        let gps_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let gps_rov1 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let gps_rov2 = SatelliteId { constellation: Constellation::Gps, prn: 3 };
        let gps_rov3 = SatelliteId { constellation: Constellation::Gps, prn: 4 };
        let gal_ref = SatelliteId { constellation: Constellation::Galileo, prn: 10 };
        let gal_rov = SatelliteId { constellation: Constellation::Galileo, prn: 11 };

        let lam = 0.19029367279836487;
        state.add_ambiguity(gps_ref, 1, 10.1 * lam, 0.0001);
        state.add_ambiguity(gps_rov1, 1, 15.15 * lam, 0.0001);
        state.add_ambiguity(gps_rov2, 1, 20.05 * lam, 0.0001);
        state.add_ambiguity(gps_rov3, 1, 25.18 * lam, 0.0001);
        state.add_ambiguity(gal_ref, 1, 30.12 * lam, 0.0001);
        state.add_ambiguity(gal_rov, 1, 36.21 * lam, 0.0001);

        for &(sat, freq) in &[(gps_ref, 1), (gps_rov1, 1), (gps_rov2, 1), (gps_rov3, 1), (gal_ref, 1), (gal_rov, 1)] {
            state.locktimes.insert((sat, freq), 100);
        }

        use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris, GalileoEphemeris};
        let ephemerides = vec![
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_ref, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_rov1, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.1, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_rov2, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.2, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_rov3, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.3, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Galileo(GalileoEphemeris { 
                sat: gal_ref, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5440.6, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, bgd_e1_e5a: 0.0,
                iod_nav: 0,
            }),
            Ephemeris::Galileo(GalileoEphemeris { 
                sat: gal_rov, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5440.6, delta_n: 0.0,
                omega0: 0.5, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, bgd_e1_e5a: 0.0,
                iod_nav: 0,
            }),
        ];

        state.resolve_ambiguities(&ephemerides, 4).expect("AR should run");

        assert!(state.is_fixed, "Should achieve fix with multi-constellation support");
        
        let idx_ref = state.ambiguity_keys.iter().position(|&(s, f)| s == gps_ref && f == 1).unwrap();
        let idx_rov = state.ambiguity_keys.iter().position(|&(s, f)| s == gps_rov1 && f == 1).unwrap();
        let dd_gps = (state.ambiguities[idx_rov] - state.ambiguities[idx_ref]) / lam;
        assert!((dd_gps.round() - 5.0).abs() < 1e-6);

        let idx_ref_gal = state.ambiguity_keys.iter().position(|&(s, f)| s == gal_ref && f == 1).unwrap();
        let idx_rov_gal = state.ambiguity_keys.iter().position(|&(s, f)| s == gal_rov && f == 1).unwrap();
        let dd_gal = (state.ambiguities[idx_rov_gal] - state.ambiguities[idx_ref_gal]) / lam;
        assert!((dd_gal.round() - 6.0).abs() < 1e-6);
    }

    #[test]
    fn test_double_difference_eliminates_clocks() {
        let sat_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_a = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let true_r_rover_ref = 20_000_000.0;
        let true_r_rover_a   = 21_000_000.0;
        let true_r_base_ref  = 20_005_000.0;
        let true_r_base_a    = 21_004_000.0;
        let rover_clk = 300.0; let base_clk = -150.0;
        let sat_ref_clk = 1000.0; let sat_a_clk = -500.0;
        let rover_ref_obs = DdObservation { sat: sat_ref, pr_l1: true_r_rover_ref + rover_clk - sat_ref_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let rover_a_obs = DdObservation { sat: sat_a, pr_l1: true_r_rover_a + rover_clk - sat_a_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let base_ref_obs = DdObservation { sat: sat_ref, pr_l1: true_r_base_ref + base_clk - sat_ref_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let base_a_obs = DdObservation { sat: sat_a, pr_l1: true_r_base_a + base_clk - sat_a_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let dd = compute_double_difference(&rover_ref_obs, &rover_a_obs, &base_ref_obs, &base_a_obs, f1, f2, f1, f2);
        assert!((dd - 1000.0).abs() < 1e-6);
    }
}
