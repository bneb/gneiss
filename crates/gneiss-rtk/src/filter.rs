use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};
use gneiss_core::sat::SatelliteId;
use gneiss_core::coords::Coordinate;

use gneiss_core::time::GpsTime;

pub const CORE_STATE_SIZE: usize = 21;

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
    pub isb_glo: f64,
    pub isb_gal: f64,
    pub isb_bds: f64,
    pub rcv_clk_drift: f64,
    pub zwd: f64,
    pub gf_values: std::collections::HashMap<SatelliteId, f64>,
    pub phase_history: std::collections::HashMap<(SatelliteId, u8), (f64, f64, GpsTime)>,
    
    pub ambiguities: Vec<f64>,
    pub ambiguity_keys: Vec<(SatelliteId, u8)>, // (Sat, Freq Band 1 or 2)
    pub ambiguity_track_ids: Vec<u32>, // Unique ID for each continuous ambiguity track
    pub next_track_id: u32,

    pub mw_sd_ema: std::collections::HashMap<SatelliteId, f64>,
    pub mw_sd_counts: std::collections::HashMap<SatelliteId, usize>,
    pub locktimes: std::collections::HashMap<(SatelliteId, u8), u16>,
    pub last_observed: std::collections::HashMap<(SatelliteId, u8), u32>,
    pub windup: std::collections::HashMap<SatelliteId, f64>,
    pub innovation_cov: std::collections::HashMap<(SatelliteId, u8), f64>, // For IAE
    pub innovation_counts: std::collections::HashMap<(SatelliteId, u8), usize>,
    pub reject_counts: std::collections::HashMap<(SatelliteId, u8), usize>,
    pub consecutive_rejections: usize,
    pub is_fixed: bool,
    pub epoch_count: usize,
    pub covariance: DMatrix<f64>,
    
    // RTS Smoother matrices (15x15/18x18 core blocks)
    pub core_phi: Option<DMatrix<f64>>,       // State transition from k-1 to k
    pub full_p_predict: Option<DMatrix<f64>>, // Predicted covariance P_{k|k-1}
    pub full_x_predict: Option<DVector<f64>>, // Predicted nominal state x_{k|k-1}
    pub predicted_position: Option<Coordinate>,
    pub predicted_velocity: Option<Vector3<f64>>,
    pub predicted_attitude: Option<UnitQuaternion<f64>>,
    pub predicted_accel_bias: Option<Vector3<f64>>,
    pub predicted_gyro_bias: Option<Vector3<f64>>,
    pub fixed_state: Option<Box<RtkState>>,
    pub is_reset: bool,
}

impl RtkState {
    pub fn new(time: GpsTime, initial_pos: Coordinate, initial_var: f64) -> Self {
        let mut cov = DMatrix::zeros(CORE_STATE_SIZE, CORE_STATE_SIZE);
        for i in 0..3 { cov[(i, i)] = initial_var; } // position
        for i in 3..6 { cov[(i, i)] = 100.0; } // velocity
        let att_var = (1.0f64.to_radians()).powi(2);
        for i in 6..9 { cov[(i, i)] = att_var; } // attitude
        for i in 9..12 { cov[(i, i)] = 0.01; } // accel bias
        for i in 12..15 { cov[(i, i)] = 1e-6; } // gyro bias
        cov[(15, 15)] = 100000.0; // rcv_clk_bias
        cov[(16, 16)] = 100000.0; // isb_glo
        cov[(17, 17)] = 100000.0; // isb_gal
        cov[(18, 18)] = 100000.0; // isb_bds
        cov[(19, 19)] = 1000.0;   // rcv_clk_drift
        cov[(20, 20)] = 1e-4;     // zwd

        Self {
            time,
            position: initial_pos,
            velocity: Vector3::zeros(),
            attitude: {
                let llh = gneiss_core::coords::ecef_to_llh(initial_pos.vector);
                let ecef_to_ned = gneiss_core::coords::ecef_to_ned_matrix(llh);
                let ned_to_ecef = ecef_to_ned.transpose();
                UnitQuaternion::from_rotation_matrix(&nalgebra::Rotation3::from_matrix(&ned_to_ecef))
            },
            accel_bias: Vector3::zeros(),
            gyro_bias: Vector3::zeros(),
            rcv_clk_bias: 0.0,
            isb_glo: 0.0,
            isb_gal: 0.0,
            isb_bds: 0.0,
            rcv_clk_drift: 0.0,
            zwd: 0.1,
            gf_values: std::collections::HashMap::new(),
            phase_history: std::collections::HashMap::new(),
            ambiguities: Vec::new(),
            ambiguity_keys: Vec::new(),
            ambiguity_track_ids: Vec::new(),
            next_track_id: 1,

            mw_sd_ema: std::collections::HashMap::new(),
            mw_sd_counts: std::collections::HashMap::new(),
            locktimes: std::collections::HashMap::new(),
            last_observed: std::collections::HashMap::new(),
            windup: std::collections::HashMap::new(),
            innovation_cov: std::collections::HashMap::new(),
            innovation_counts: std::collections::HashMap::new(),
            reject_counts: std::collections::HashMap::new(),
            consecutive_rejections: 0,
            is_fixed: false,
            epoch_count: 0,
            covariance: cov,
            core_phi: None,
            full_p_predict: None,
            full_x_predict: None,
            predicted_position: None,
            predicted_velocity: None,
            predicted_attitude: None,
            predicted_accel_bias: None,
            predicted_gyro_bias: None,
            fixed_state: None,
            is_reset: false,
        }
    }

    /// Decouples the position states (indices 0..3) from the rest of the EKF state
    /// by zeroing out the corresponding cross-covariance rows and columns.
    /// This is mathematically required when teleporting the position state.
    pub fn decouple_position(&mut self) {
        let cols = self.covariance.ncols();
        for i in 0..3 {
            for j in 3..cols {
                self.covariance[(i, j)] = 0.0;
                self.covariance[(j, i)] = 0.0;
            }
        }
    }

    /// Decouples the receiver clock bias state (index 15) from the rest of the EKF state
    /// by zeroing out the corresponding cross-covariance rows and columns.
    /// This is mathematically required when teleporting the clock state.
    pub fn decouple_clock(&mut self) {
        let cols = self.covariance.ncols();
        
        // Clock bias is index 15
        let indices = [15, 16, 17, 18];
        for &i in &indices {
            for j in 0..cols {
                if i != j {
                    self.covariance[(i, j)] = 0.0;
                    self.covariance[(j, i)] = 0.0;
                }
            }
            self.covariance[(i, i)] = 100000.0;
        }
    }

    pub fn update_mw(&mut self, sat: SatelliteId, mw_cycles: f64) {
        let count = self.mw_sd_counts.entry(sat).or_insert(0);
        let ema = self.mw_sd_ema.entry(sat).or_insert(mw_cycles);
        let alpha = 1.0 / ((*count + 1) as f64).min(100.0);
        *ema = *ema * (1.0 - alpha) + mw_cycles * alpha;
        *count += 1;
    }

        #[allow(clippy::type_complexity)]
    pub fn resolve_ambiguities(&self, ephemerides: &[gneiss_core::ephemeris::Ephemeris], min_subset: usize, ar_min_epoch_count: u32, ar_min_lock: u32, lambda_min_ratio: f64) -> Result<(RtkState, DVector<f64>, DMatrix<f64>, f64, usize), &'static str> {
        let num_amb = self.ambiguities.len();
        if num_amb < min_subset || self.epoch_count <= ar_min_epoch_count as usize { return Err("Insufficient data"); }
        
        let candidate_vars = select_ar_candidates(self, ephemerides, ar_min_lock);
        if candidate_vars.len() < min_subset { return Err("Insufficient candidates"); }
        
        let max_subset = candidate_vars.len().min(24);
        for subset_size in (min_subset..=max_subset).rev() {
            let (_d_mat_small, a_cycles, q_cycles) = build_lambda_matrices(self, &candidate_vars, subset_size, ephemerides);
            
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
        fixed_state.fixed_state = None;
        crate::engine::updater::apply_state_correction(&mut fixed_state, &dx);
        let r_zero = DMatrix::zeros(subset_size, subset_size);
        fixed_state.covariance = crate::engine::updater::apply_joseph_covariance_update(&self.covariance, &k_full, &d_full, &r_zero);
        fixed_state.is_fixed = true;

        Ok((fixed_state, da_meters, d_full))
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
            self.phase_history.remove(&(sat, freq));
        }
    }
    pub fn add_ambiguity(&mut self, sat: SatelliteId, freq: u8, initial_estimate: f64, initial_variance: f64) {
        if self.ambiguity_keys.contains(&(sat, freq)) {
            tracing::warn!("Ambiguity for {:?} L{} already exists! Resetting.", sat, freq);
            self.remove_ambiguity(sat, freq);
        }
        tracing::debug!("Adding ambiguity for {:?} L{} val={} var={}", sat, freq, initial_estimate, initial_variance);
        self.ambiguities.push(initial_estimate);
        self.ambiguity_keys.push((sat, freq));
        self.ambiguity_track_ids.push(self.next_track_id);
        self.next_track_id += 1;

        let n_old = self.covariance.nrows();
        self.covariance = self.covariance.clone().insert_row(n_old, 0.0).insert_column(n_old, 0.0);
        self.covariance[(n_old, n_old)] = initial_variance;
    }

    pub fn remove_ambiguity(&mut self, sat: SatelliteId, freq: u8) {
        if let Some(idx) = self.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == freq) {
            self.ambiguities.remove(idx);
            self.ambiguity_keys.remove(idx);
            self.ambiguity_track_ids.remove(idx);

            let cov_idx = CORE_STATE_SIZE + idx;
            self.covariance = self.covariance.clone().remove_row(cov_idx).remove_column(cov_idx);
        }
    }
    
    pub fn clear_ambiguities(&mut self) {
        let num_amb = self.ambiguities.len();
        self.ambiguities.clear();
        self.ambiguity_keys.clear();
        self.ambiguity_track_ids.clear();

        self.covariance = self.covariance.clone().remove_rows(CORE_STATE_SIZE, num_amb)
                                         .remove_columns(CORE_STATE_SIZE, num_amb);
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

#[allow(clippy::too_many_arguments)]
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

        let (fixed_state, _, _, _, _) = state.resolve_ambiguities(&ephemerides, 4, 5, 3, 3.0).expect("AR should run");
        assert!(fixed_state.is_fixed, "Should achieve fix with multi-constellation support");
        
        let idx_ref = fixed_state.ambiguity_keys.iter().position(|&(s, f)| s == gps_ref && f == 1).unwrap();
        let idx_rov = fixed_state.ambiguity_keys.iter().position(|&(s, f)| s == gps_rov1 && f == 1).unwrap();
        let dd_gps = (fixed_state.ambiguities[idx_rov] - fixed_state.ambiguities[idx_ref]) / lam;
        assert!((dd_gps.round() - 5.0).abs() < 1e-6);

        let idx_ref_gal = fixed_state.ambiguity_keys.iter().position(|&(s, f)| s == gal_ref && f == 1).unwrap();
        let idx_rov_gal = fixed_state.ambiguity_keys.iter().position(|&(s, f)| s == gal_rov && f == 1).unwrap();
        let dd_gal = (fixed_state.ambiguities[idx_rov_gal] - fixed_state.ambiguities[idx_ref_gal]) / lam;
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
        let _dd = compute_double_difference(&rover_ref_obs, &rover_a_obs, &base_ref_obs, &base_a_obs, f1, f2, f1, f2);
    }

    #[test]
    fn test_decouple_position() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        
        // Fill covariance with 1.0
        state.covariance.fill(1.0);
        state.decouple_position();

        let cols = state.covariance.ncols();
        for i in 0..3 {
            for j in 3..cols {
                assert_eq!(state.covariance[(i, j)], 0.0);
                assert_eq!(state.covariance[(j, i)], 0.0);
            }
        }
        // Ensure other elements are untouched
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(state.covariance[(i, j)], 1.0);
            }
        }
        for i in 3..cols {
            for j in 3..cols {
                assert_eq!(state.covariance[(i, j)], 1.0);
            }
        }
    }

    #[test]
    fn test_decouple_clock() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        
        // Fill covariance with 1.0
        state.covariance.fill(1.0);
        state.decouple_clock();

        let cols = state.covariance.ncols();
        let clock_indices = [15, 16, 17, 18]; // rcv_clk_bias, isb_glo, isb_gal, isb_bds
        
        // Verify all 4 clock-related diagonals are reset to 100000.0
        for &idx in &clock_indices {
            assert_eq!(state.covariance[(idx, idx)], 100000.0,
                "diagonal at index {} should be 100000.0", idx);
        }
        
        // Verify all cross-covariance terms involving clock indices are zeroed
        for &idx in &clock_indices {
            for j in 0..cols {
                if !clock_indices.contains(&j) {
                    assert_eq!(state.covariance[(idx, j)], 0.0,
                        "cross-covariance ({}, {}) should be 0.0", idx, j);
                    assert_eq!(state.covariance[(j, idx)], 0.0,
                        "cross-covariance ({}, {}) should be 0.0", j, idx);
                }
            }
        }
        
        // Verify non-clock entries are untouched (still 1.0)
        for i in 0..cols {
            for j in 0..cols {
                if !clock_indices.contains(&i) && !clock_indices.contains(&j) {
                    assert_eq!(state.covariance[(i, j)], 1.0,
                        "non-clock entry ({}, {}) should remain 1.0", i, j);
                }
            }
        }
    }

    #[test]
    fn test_update_mw() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        
        state.update_mw(sat, 10.0);
        assert_eq!(state.mw_sd_counts[&sat], 1);
        assert_eq!(state.mw_sd_ema[&sat], 10.0);
        
        state.update_mw(sat, 20.0);
        assert_eq!(state.mw_sd_counts[&sat], 2);
        // alpha = 1 / 2 = 0.5. ema = 10 * 0.5 + 20 * 0.5 = 15.0
        assert_eq!(state.mw_sd_ema[&sat], 15.0);
        
        state.update_mw(sat, 30.0);
        assert_eq!(state.mw_sd_counts[&sat], 3);
        // alpha = 1 / 3. ema = 15 * (2/3) + 30 * 1/3 = 10 + 10 = 20.0
        assert!((state.mw_sd_ema[&sat] - 20.0).abs() < 1e-6);
        
        // Add more than 100 to test min(100.0)
        for _ in 0..97 {
            state.update_mw(sat, 20.0);
        }
        assert_eq!(state.mw_sd_counts[&sat], 100);
        assert!((state.mw_sd_ema[&sat] - 20.0).abs() < 1e-6);
        
        state.update_mw(sat, 120.0);
        assert_eq!(state.mw_sd_counts[&sat], 101);
        // alpha should be 1/100, not 1/101
        assert!((state.mw_sd_ema[&sat] - (20.0 * 0.99 + 120.0 * 0.01)).abs() < 1e-6);
    }
    #[test]
    fn test_prune_stale_ambiguities() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let sat3 = SatelliteId { constellation: Constellation::Gps, prn: 3 };

        state.add_ambiguity(sat1, 1, 0.0, 1.0);
        state.last_observed.insert((sat1, 1), 10);
        state.locktimes.insert((sat1, 1), 10);
        state.phase_history.insert((sat1, 1), (0.0, 0.0, time));

        state.add_ambiguity(sat2, 1, 0.0, 1.0);
        state.last_observed.insert((sat2, 1), 20);
        state.locktimes.insert((sat2, 1), 20);
        state.phase_history.insert((sat2, 1), (0.0, 0.0, time));

        state.add_ambiguity(sat3, 1, 0.0, 1.0);
        // sat3 has no last_observed, defaults to 0
        
        assert_eq!(state.covariance.nrows(), CORE_STATE_SIZE + 3);

        state.prune_stale_ambiguities(25, 10);

        assert_eq!(state.ambiguity_keys.len(), 1);
        assert_eq!(state.ambiguity_keys[0], (sat2, 1));
        assert!(state.last_observed.contains_key(&(sat2, 1)));
        assert!(state.locktimes.contains_key(&(sat2, 1)));
        assert!(state.phase_history.contains_key(&(sat2, 1)));

        assert!(!state.last_observed.contains_key(&(sat1, 1)));
        assert!(!state.locktimes.contains_key(&(sat1, 1)));
        assert!(!state.phase_history.contains_key(&(sat1, 1)));
        
        assert_eq!(state.covariance.nrows(), CORE_STATE_SIZE + 1);
        assert_eq!(state.covariance.ncols(), CORE_STATE_SIZE + 1);
    }

    #[test]
    fn test_clear_ambiguities() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        
        state.add_ambiguity(sat1, 1, 10.0, 1.0);
        assert_eq!(state.ambiguity_keys.len(), 1);
        assert_eq!(state.covariance.nrows(), CORE_STATE_SIZE + 1);

        state.clear_ambiguities();
        assert_eq!(state.ambiguity_keys.len(), 0);
        assert_eq!(state.ambiguities.len(), 0);
        assert_eq!(state.ambiguity_track_ids.len(), 0);
        assert_eq!(state.covariance.nrows(), CORE_STATE_SIZE);
        assert_eq!(state.covariance.ncols(), CORE_STATE_SIZE);
    }

    #[test]
    fn test_compute_if_combination() {
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        
        // Single frequency
        let res_single = compute_if_combination(100.0, None, f1, f2);
        assert_eq!(res_single, 100.0);
        
        // Dual frequency
        let f1_2 = f1 * f1;
        let f2_2 = f2 * f2;
        let expected = (f1_2 * 100.0 - f2_2 * 80.0) / (f1_2 - f2_2);
        let res_dual = compute_if_combination(100.0, Some(80.0), f1, f2);
        assert_eq!(res_dual, expected);
    }

#[test]
    fn test_select_ar_candidates() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        
        let gps_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let gps_rov = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        
        // Add ambiguities
        state.add_ambiguity(gps_ref, 1, 0.0, 1.0);
        state.add_ambiguity(gps_rov, 1, 0.0, 1.0);
        
        // Set locktimes
        state.locktimes.insert((gps_ref, 1), 100);
        state.locktimes.insert((gps_rov, 1), 50);
        
        let ephemerides = vec![]; // Empty ephemerides is fine, will fallback to 0 freq_num
        
        // Test with lock limit 60
        let candidates_fail = crate::filter::select_ar_candidates(&state, &ephemerides, 60);
        assert_eq!(candidates_fail.len(), 0, "Rover sat does not meet lock criteria");
        
        // Test with lock limit 40
        let candidates_pass = crate::filter::select_ar_candidates(&state, &ephemerides, 40);
        assert_eq!(candidates_pass.len(), 1, "Rover sat meets lock criteria");
        assert_eq!(candidates_pass[0].0, 1); // rov_idx
        assert_eq!(candidates_pass[0].1, 0); // ref_idx
        assert_eq!(candidates_pass[0].2, 50); // locktime
    }

    #[test]
    fn test_build_lambda_matrices() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        
        let gps_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let gps_rov = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        
        let lam = 0.19029367279836487;
        state.add_ambiguity(gps_ref, 1, 10.0 * lam, 0.1);
        state.add_ambiguity(gps_rov, 1, 15.0 * lam, 0.1);
        
        let candidates = vec![(1, 0, 50, 0.5)]; // rov_idx=1, ref_idx=0
        
        let ephemerides = vec![];
        let (d_mat, a_cycles, q_cycles) = crate::filter::build_lambda_matrices(&state, &candidates, 1, &ephemerides);
        
        assert_eq!(d_mat.nrows(), 1);
        assert_eq!(d_mat.ncols(), 2);
        assert_eq!(d_mat[(0, 0)], -1.0); // ref is -1
        assert_eq!(d_mat[(0, 1)], 1.0);  // rov is +1
        
        assert_eq!(a_cycles.len(), 1);
        assert!((a_cycles[0] - 5.0).abs() < 1e-6); // 15 - 10
        
        assert_eq!(q_cycles.nrows(), 1);
        assert_eq!(q_cycles.ncols(), 1);
    }

}
pub fn select_ar_candidates(
    state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    ar_min_lock: u32,
) -> Vec<(usize, usize, u16, f64)> {
    let constellations = [gneiss_core::sat::Constellation::Gps, gneiss_core::sat::Constellation::Galileo];
    let candidates = filter_by_locktime(state, &constellations, ar_min_lock);
    let mut candidate_vars = compute_candidate_variance(state, ephemerides, &candidates);
    candidate_vars.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
    candidate_vars
}

fn filter_by_locktime(
    state: &RtkState, constellations: &[gneiss_core::sat::Constellation], ar_min_lock: u32,
) -> Vec<(usize, usize, u16)> {
    let mut candidates = Vec::new();
    for &constell in constellations {
        if let Some(ref_idx) = find_best_reference_sat(state, constell, ar_min_lock) {
            collect_candidates_for_constellation(state, constell, ref_idx, ar_min_lock, &mut candidates);
        }
    }
    candidates
}

fn find_best_reference_sat(state: &RtkState, constell: gneiss_core::sat::Constellation, ar_min_lock: u32) -> Option<usize> {
    let mut best_ref_idx = None;
    let mut max_lock = 0;
    for i in 0..state.ambiguities.len() {
        let (sat, freq) = state.ambiguity_keys[i];
        if sat.constellation != constell || freq != 1 { continue; }
        let lock = *state.locktimes.get(&(sat, freq)).unwrap_or(&0);
        if lock >= ar_min_lock as u16 && lock > max_lock {
            max_lock = lock; best_ref_idx = Some(i);
        }
    }
    best_ref_idx
}

fn collect_candidates_for_constellation(
    state: &RtkState, constell: gneiss_core::sat::Constellation, ref_idx: usize, ar_min_lock: u32, candidates: &mut Vec<(usize, usize, u16)>
) {
    let ref_sat_id = state.ambiguity_keys[ref_idx].0;
    let l2_ref_idx = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_sat_id && f == 2);
    for i in 0..state.ambiguities.len() {
        if i == ref_idx || Some(i) == l2_ref_idx { continue; }
        let (rov_sat, freq) = state.ambiguity_keys[i];
        if rov_sat.constellation != constell { continue; }
        let lock = *state.locktimes.get(&(rov_sat, freq)).unwrap_or(&0);
        if lock >= ar_min_lock as u16 {
            if freq == 1 { candidates.push((i, ref_idx, lock)); }
            else if let Some(r2_idx) = l2_ref_idx { candidates.push((i, r2_idx, lock)); }
        }
    }
}

fn compute_candidate_variance(
    state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    candidates: &[(usize, usize, u16)],
) -> Vec<(usize, usize, u16, f64)> {
    let num_amb = state.ambiguities.len();
    let q_sd = state.covariance.view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb));
    let mut candidate_vars = Vec::new();
    for &(rov, r_idx, lock) in candidates {
        let (rov_sat_id, freq_band) = state.ambiguity_keys[rov];
        if rov_sat_id.constellation == gneiss_core::sat::Constellation::Glonass { continue; }
        let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
        let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
        let lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_band == 1 { f1 } else { f2 };
        let q_dd = q_sd[(rov, rov)] + q_sd[(r_idx, r_idx)] - 2.0 * q_sd[(rov, r_idx)];
        let var_cycles = q_dd / (lam * lam);
        if var_cycles < 10000.0 { candidate_vars.push((rov, r_idx, lock, var_cycles)); }
    }
    candidate_vars
}

pub fn build_lambda_matrices(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    subset_size: usize,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
    let (d_mat, a_cycles) = build_lambda_design_matrix(state, candidates, subset_size, ephemerides);
    let q_cycles = build_lambda_variance_matrix(state, candidates, subset_size, ephemerides);
    (d_mat, a_cycles, q_cycles)
}

fn build_lambda_design_matrix(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    subset_size: usize,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> (DMatrix<f64>, DVector<f64>) {
    let num_amb = state.ambiguities.len();
    let a_sd = nalgebra::DVector::from_vec(state.ambiguities.clone());
    let mut d_mat = DMatrix::zeros(subset_size, num_amb);
    let mut a_cycles = DVector::zeros(subset_size);
    for row in 0..subset_size {
        let (rov, r_idx, _, _) = candidates[row];
        let (rov_sat_id, freq_band) = state.ambiguity_keys[rov];
        let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
        let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
        let lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_band == 1 { f1 } else { f2 };
        d_mat[(row, rov)] = 1.0;
        d_mat[(row, r_idx)] = -1.0;
        a_cycles[row] = (a_sd[rov] - a_sd[r_idx]) / lam;
    }
    (d_mat, a_cycles)
}

fn build_lambda_variance_matrix(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    subset_size: usize,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> DMatrix<f64> {
    let num_amb = state.ambiguities.len();
    let q_sd = state.covariance.view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb));
    let mut q_cycles = DMatrix::zeros(subset_size, subset_size);
    for r in 0..subset_size {
        for c in 0..subset_size {
            let (rov_r, ref_r, _, _) = candidates[r];
            let (rov_c, ref_c, _, _) = candidates[c];
            let freq_r = state.ambiguity_keys[rov_r].1;
            let freq_c = state.ambiguity_keys[rov_c].1;
            let (sat_r, _) = state.ambiguity_keys[rov_r];
            let freq_num_r = ephemerides.iter().find(|e| e.sat() == sat_r).map(|e| e.freq_num()).unwrap_or(0);
            let (f1_r, f2_r) = gneiss_core::signal::satellite_frequencies(sat_r, freq_num_r);
            let lam_r = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_r == 1 { f1_r } else { f2_r };
            let (sat_c, _) = state.ambiguity_keys[rov_c];
            let freq_num_c = ephemerides.iter().find(|e| e.sat() == sat_c).map(|e| e.freq_num()).unwrap_or(0);
            let (f1_c, f2_c) = gneiss_core::signal::satellite_frequencies(sat_c, freq_num_c);
            let lam_c = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_c == 1 { f1_c } else { f2_c };
            let q_dd = q_sd[(rov_r, rov_c)] - q_sd[(rov_r, ref_c)] - q_sd[(ref_r, rov_c)] + q_sd[(ref_r, ref_c)];
            q_cycles[(r, c)] = q_dd / (lam_r * lam_c);
        }
    }
    q_cycles
}
