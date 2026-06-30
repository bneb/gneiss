use gneiss_core::coords::Coordinate;
use gneiss_core::sat::SatelliteId;
use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};
use std::collections::VecDeque;

use gneiss_core::time::GpsTime;

pub const CORE_STATE_SIZE: usize = 21;

/// Ring buffer for accumulating raw DD pseudorange over a sliding window.
/// Each entry stores (dd_pr, ref_sat, rover_position) for motion compensation.
/// Maintains a running sum for O(1) mean computation.
#[derive(Debug, Clone)]
pub struct PrRingBuffer {
    buf: VecDeque<(f64, SatelliteId, nalgebra::Vector3<f64>)>,
    sum: f64,
    capacity: usize,
}

impl PrRingBuffer {
    pub fn new(capacity: usize) -> Self {
        Self { buf: VecDeque::with_capacity(capacity), sum: 0.0, capacity }
    }

    pub fn push(&mut self, dd_pr: f64, ref_sat: SatelliteId, rov_pos: nalgebra::Vector3<f64>) {
        if self.buf.len() >= self.capacity {
            self.sum -= self.buf.pop_front().map(|(v, _, _)| v).unwrap_or(0.0);
        }
        self.buf.push_back((dd_pr, ref_sat, rov_pos));
        self.sum += dd_pr;
    }

    pub fn mean(&self) -> Option<f64> {
        if self.buf.is_empty() { None }
        else { Some(self.sum / self.buf.len() as f64) }
    }

    /// Motion-compensated mean: shifts each DD PR to the reference position
    /// by subtracting the geometric DD change between the entry's rover position
    /// and the reference position. Returns (compensated_mean, ref_position).
    pub fn compensated_mean(
        &self,
        sat: gneiss_core::sat::SatelliteId,
        ref_sat: gneiss_core::sat::SatelliteId,
        ref_pos: nalgebra::Vector3<f64>,
        ephemerides: &[gneiss_core::ephemeris::Ephemeris],
        time: gneiss_core::time::GpsTime,
        base_pos: nalgebra::Vector3<f64>,
        base_time: gneiss_core::time::GpsTime,
    ) -> Option<f64> {
        if self.buf.is_empty() { return None; }
        let eph_sat = ephemerides.iter().find(|e| e.sat() == sat)?;
        let eph_ref = ephemerides.iter().find(|e| e.sat() == ref_sat)?;
        let (ref_sat_pos, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, time, ref_pos);
        let (ref_ref_pos, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, time, ref_pos);
        let (ref_bas_sat, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, base_time, base_pos);
        let (ref_bas_ref, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, base_time, base_pos);
        let ref_geom = crate::engine::measurement_math::compute_geometric_dd(ref_pos, base_pos, ref_sat_pos, ref_ref_pos, ref_bas_sat, ref_bas_ref);

        let mut comp_sum = 0.0f64;
        for &(dd_pr, _, entry_pos) in &self.buf {
            let (es, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, time, entry_pos);
            let (er, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, time, entry_pos);
            let (bs, _) = crate::engine::measurement_math::get_sat_state(eph_sat, 0.0, 0.0, base_time, base_pos);
            let (br, _) = crate::engine::measurement_math::get_sat_state(eph_ref, 0.0, 0.0, base_time, base_pos);
            let entry_geom = crate::engine::measurement_math::compute_geometric_dd(entry_pos, base_pos, es, er, bs, br);
            // Shift DD PR to reference position: dd_pr_ref = dd_pr - (entry_geom - ref_geom)
            comp_sum += dd_pr - (entry_geom - ref_geom);
        }
        Some(comp_sum / self.buf.len() as f64)
    }

    pub fn count(&self) -> usize { self.buf.len() }
}

/// Initial variance for receiver clock bias (m²).
///
/// RTKLIB uses VAR_CLK = 100² = 10,000 m² for the PPP filter, reset
/// every epoch (white-noise model).  u-blox F9P TCXO typically drifts
/// at ~1 ms/s unsteered, producing kilometer-level position errors
/// within seconds, so the clock must be estimated from measurements
/// rather than trusted.
///
/// 10,000 m² corresponds to σ = 100 m — tight enough for the filter
/// to converge quickly, loose enough not to bias position when the
/// clock estimate is wrong.
///
/// Sources: RTKLIB 2.4.3 `rtkcmn.c:udclk_ppp()`; NovAtel OEM7
/// documentation recommends 10²–10⁴ m² for TCXO receivers.
pub const INITIAL_CLOCK_BIAS_VARIANCE: f64 = 10_000.0;

/// Initial variance for inter-system biases (GLO/GAL/BDS) in m².
///
/// ISBs are directly observable from single-differenced pseudorange
/// residuals across constellations.  Commercial products (Trimble, Leica)
/// model them as piece-wise constants with σ ≈ 5–10 m initial uncertainty.
/// The RTKLIB MIS model eliminates ISB states entirely by estimating
/// separate receiver clocks per constellation.
///
/// 100 m² corresponds to σ = 10 m — well within the observable range
/// from a single epoch of pseudorange.
///
/// Sources: Tian et al. (2020), *Advances in Space Research*; RTKLIB 2.4.3
pub const INITIAL_ISB_VARIANCE: f64 = 100.0;

/// Predicted variance for clock/ISB states after one epoch of process
/// noise (initial + q·Δt).  Used when constructing synthetic prediction
/// covariances in tests; the real forward filter computes this dynamically
/// from the process noise model.
pub const PREDICTED_CLOCK_VARIANCE: f64 = 20_000.0;
/// Predicted ISB variance for tests.
pub const PREDICTED_ISB_VARIANCE: f64 = 200.0;

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
    pub ambiguity_track_ids: Vec<u32>,          // Unique ID for each continuous ambiguity track
    pub next_track_id: u32,

    pub mw_sd_ema: std::collections::HashMap<SatelliteId, f64>,
    pub mw_sd_counts: std::collections::HashMap<SatelliteId, usize>,
    /// Narrow-lane SD EMA (cycles) for NL integer validation
    pub nl_sd_ema: std::collections::HashMap<SatelliteId, f64>,
    pub nl_sd_counts: std::collections::HashMap<SatelliteId, usize>,
    /// Sliding-window DD pseudorange accumulator: (rov_sat, ref_sat) → ring buffer
    pub pr_dd_window: std::collections::HashMap<(SatelliteId, SatelliteId), PrRingBuffer>,
    pub gf_prev: std::collections::HashMap<SatelliteId, f64>,
    pub mw_prev: std::collections::HashMap<SatelliteId, f64>,
    pub locktimes: std::collections::HashMap<(SatelliteId, u8), u16>,
    pub last_observed: std::collections::HashMap<(SatelliteId, u8), u32>,
    pub windup: std::collections::HashMap<SatelliteId, f64>,
    pub current_ref_sat: std::collections::HashMap<gneiss_core::sat::Constellation, SatelliteId>,
    pub innovation_cov: std::collections::HashMap<(SatelliteId, u8), f64>, // For IAE
    pub innovation_counts: std::collections::HashMap<(SatelliteId, u8), usize>,
    pub reject_counts: std::collections::HashMap<(SatelliteId, u8), usize>,
    pub consecutive_rejections: usize,
    pub is_fixed: bool,
    pub epoch_count: usize,
    pub covariance: DMatrix<f64>,

    // RTS Smoother matrices (15x15/18x18 core blocks)
    pub core_phi: Option<DMatrix<f64>>, // State transition from k-1 to k
    pub full_p_predict: Option<DMatrix<f64>>, // Predicted covariance P_{k|k-1}
    pub full_x_predict: Option<DVector<f64>>, // Predicted nominal state x_{k|k-1}
    pub predicted_position: Option<Coordinate>,
    pub predicted_velocity: Option<Vector3<f64>>,
    pub predicted_attitude: Option<UnitQuaternion<f64>>,
    pub predicted_accel_bias: Option<Vector3<f64>>,
    pub predicted_gyro_bias: Option<Vector3<f64>>,
    pub fixed_state: Option<Box<RtkState>>,
    pub is_reset: bool,
    pub ins_aligned: bool,
}

/// Initializes the core 21x21 covariance matrix for a new EKF state.
fn init_covariance(initial_var: f64) -> DMatrix<f64> {
    let mut cov = DMatrix::zeros(CORE_STATE_SIZE, CORE_STATE_SIZE);
    for i in 0..3 {
        cov[(i, i)] = initial_var;
    } // position
    for i in 3..6 {
        cov[(i, i)] = 100.0;
    } // velocity
    let att_var = (1.0f64.to_radians()).powi(2);
    for i in 6..9 {
        cov[(i, i)] = att_var;
    } // attitude
    for i in 9..12 {
        cov[(i, i)] = 0.01;
    } // accel bias
    for i in 12..15 {
        cov[(i, i)] = 1e-6;
    } // gyro bias
    cov[(15, 15)] = INITIAL_CLOCK_BIAS_VARIANCE;
    cov[(16, 16)] = INITIAL_ISB_VARIANCE;
    cov[(17, 17)] = INITIAL_ISB_VARIANCE;
    cov[(18, 18)] = INITIAL_ISB_VARIANCE;
    cov[(19, 19)] = 1000.0; // rcv_clk_drift
    cov[(20, 20)] = 1.0; // zwd
    cov
}

/// Computes the initial attitude from the position by aligning the body frame
/// with the local NED frame (level attitude).
fn init_attitude(initial_pos: &Coordinate) -> UnitQuaternion<f64> {
    let llh = gneiss_core::coords::ecef_to_llh(initial_pos.vector);
    let ecef_to_ned = gneiss_core::coords::ecef_to_ned_matrix(llh);
    let ned_to_ecef = ecef_to_ned.transpose();
    UnitQuaternion::from_rotation_matrix(&nalgebra::Rotation3::from_matrix(&ned_to_ecef))
}

impl RtkState {
    pub fn new(time: GpsTime, initial_pos: Coordinate, initial_var: f64) -> Self {
        let cov = init_covariance(initial_var);
        Self {
            time,
            position: initial_pos,
            velocity: Vector3::zeros(),
            attitude: init_attitude(&initial_pos),
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
            nl_sd_ema: std::collections::HashMap::new(),
            nl_sd_counts: std::collections::HashMap::new(),
            pr_dd_window: std::collections::HashMap::new(),
            gf_prev: std::collections::HashMap::new(),
            mw_prev: std::collections::HashMap::new(),
            locktimes: std::collections::HashMap::new(),
            last_observed: std::collections::HashMap::new(),
            windup: std::collections::HashMap::new(),
            current_ref_sat: std::collections::HashMap::new(),
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
            ins_aligned: false,
        }
    }

    /// Decouples the position states (indices 0..3) from the rest of the EKF state
    /// by zeroing out the corresponding cross-covariance rows and columns.
    /// Inflates position variance to reflect uncertainty after coasting — this
    /// widens the adaptive PR rejection threshold so measurements can be
    /// reacquired after a divergence episode.
    /// This is mathematically required when teleporting the position state.
    pub fn decouple_position(&mut self) {
        let cols = self.covariance.ncols();
        for i in 0..3 {
            if self.covariance[(i, i)] < 900.0 {
                self.covariance[(i, i)] = 900.0;
            }
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
            self.covariance[(i, i)] = INITIAL_CLOCK_BIAS_VARIANCE;
        }
    }

    pub fn reset_to_spp(
        &mut self,
        spp_pos: Coordinate,
        spp_state: Option<&crate::spp::SppState>,
        init_isbs: bool,
    ) {
        self.position = spp_pos;
        if self.velocity.norm() > 100.0 || self.velocity.norm().is_nan() {
            self.velocity = Vector3::zeros();
        }

        if self.covariance.nrows() > 15 {
            if let Some(spp) = spp_state {
                self.rcv_clk_bias = spp.cdt;
                if init_isbs {
                    self.isb_glo = spp.cdt_glo - spp.cdt;
                    self.isb_gal = spp.cdt_gal - spp.cdt;
                    self.isb_bds = spp.cdt_bds - spp.cdt;
                } else {
                    self.isb_glo = 0.0;
                    self.isb_gal = 0.0;
                    self.isb_bds = 0.0;
                }

                tracing::debug!("SPP Reset: cdt={:.2}, bds={:.2}, glo={:.2}, init_isbs={}, isb_bds={:.2}, isb_glo={:.2}", 
                    spp.cdt, spp.cdt_bds, spp.cdt_glo, init_isbs, self.isb_bds, self.isb_glo);
            }
            if self.rcv_clk_drift.abs() > 10000.0 || self.rcv_clk_drift.is_nan() {
                self.rcv_clk_drift = 0.0;
            }
        }

        self.clear_ambiguities();

        let cols = self.covariance.ncols();
        let mut reset_indices: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 15];
        if self.covariance.nrows() > 15 && !self.ins_aligned {
            reset_indices = (0..16).collect();
        }
        for &i in &reset_indices {
            for j in 0..cols {
                if i != j && !reset_indices.contains(&j) {
                    self.covariance[(i, j)] = 0.0;
                    self.covariance[(j, i)] = 0.0;
                }
            }
        }
        for &i in &reset_indices {
            for &j in &reset_indices {
                if i != j {
                    self.covariance[(i, j)] = 0.0;
                }
            }
        }

        for i in 0..3 {
            self.covariance[(i, i)] = 100.0;
        }
        for i in 3..6 {
            self.covariance[(i, i)] = 100.0;
        }
        if self.covariance.nrows() > 15 {
            if !self.ins_aligned {
                let att_var = (1.0f64.to_radians()).powi(2);
                for i in 6..9 {
                    self.covariance[(i, i)] = att_var;
                }
                for i in 9..12 {
                    self.covariance[(i, i)] = 0.01;
                }
                for i in 12..15 {
                    self.covariance[(i, i)] = 1e-6;
                }

                self.accel_bias = Vector3::zeros();
                self.gyro_bias = Vector3::zeros();
                let llh = gneiss_core::coords::ecef_to_llh(self.position.vector);
                let ecef_to_ned = gneiss_core::coords::ecef_to_ned_matrix(llh);
                let ned_to_ecef = ecef_to_ned.transpose();
                self.attitude = UnitQuaternion::from_rotation_matrix(
                    &nalgebra::Rotation3::from_matrix(&ned_to_ecef),
                );
            }
            self.covariance[(15, 15)] = INITIAL_CLOCK_BIAS_VARIANCE;
        }

        self.is_reset = true;
        self.consecutive_rejections = 0;
        if !self.ins_aligned {
            self.ins_aligned = false;
        }
    }

    pub fn update_mw(&mut self, sat: SatelliteId, mw_cycles: f64) {
        let count = self.mw_sd_counts.entry(sat).or_insert(0);
        let ema = self.mw_sd_ema.entry(sat).or_insert(mw_cycles);
        let alpha = 1.0 / ((*count + 1) as f64).min(100.0);
        *ema = *ema * (1.0 - alpha) + mw_cycles * alpha;
        *count += 1;
    }

    /// Update the narrow-lane SD EMA with a new measurement (cycles).
    /// Used to validate LAMBDA NL integers against time-averaged code.
    pub fn update_nl(&mut self, sat: SatelliteId, nl_cycles: f64) {
        let count = self.nl_sd_counts.entry(sat).or_insert(0);
        let ema = self.nl_sd_ema.entry(sat).or_insert(nl_cycles);
        let alpha = 1.0 / ((*count + 1) as f64).min(100.0);
        *ema = *ema * (1.0 - alpha) + nl_cycles * alpha;
        *count += 1;
    }

    pub fn resolve_ambiguities(
        &self,
        ephemerides: &[gneiss_core::ephemeris::Ephemeris],
        config: &crate::engine::EngineConfig,
    ) -> Result<crate::ambiguity::AmbiguityResolutionResult, &'static str> {
        let num_amb = self.ambiguities.len();
        if num_amb < config.lambda_min_subset
            || self.epoch_count <= config.ar_min_epoch_count as usize
        {
            tracing::debug!("AR: insufficient data (amb={} need={} epoch={} need={})",
                num_amb, config.lambda_min_subset, self.epoch_count, config.ar_min_epoch_count);
            return Err("Insufficient data");
        }

        let candidate_vars = select_ar_candidates(self, ephemerides, config.ar_min_lock);
        if candidate_vars.len() < config.lambda_min_subset {
            tracing::debug!("AR: insufficient candidates (have={} need={})",
                candidate_vars.len(), config.lambda_min_subset);
            return Err("Insufficient candidates");
        }

        let max_subset = candidate_vars.len().min(24);
        let mut best_ratio = 0.0f64;
        let mut best_success_rate = 0.0f64;
        let mut best_subset_size = 0;
        for subset_size in (config.lambda_min_subset..=max_subset).rev() {
            let (_d_mat_small, a_cycles, q_cycles) =
                build_lambda_matrices(self, &candidate_vars, subset_size, ephemerides);

            if let Ok(res) = crate::lambda::resolve_lambda(&a_cycles, &q_cycles) {
                // Diagnostics: track best attempt across all subset sizes
                if res.ratio > best_ratio {
                    best_ratio = res.ratio;
                    best_success_rate = res.success_rate;
                    best_subset_size = subset_size;
                }
                let dynamic_threshold =
                    crate::ffrt::calculate_threshold(subset_size, config.ar_ffrt_prob)
                        .max(config.lambda_min_ratio);
                // Accept if it passes the ratio test AND the MW widelane validation
                if res.ratio >= dynamic_threshold
                    || res.success_rate >= config.tuning.min_ar_success_rate
                {
                    if !validate_mw_widelane(self, &candidate_vars, &res, subset_size) {
                        tracing::debug!(
                            "AR: MW WL validation failed (subset={}, ratio={:.2}, sr={:.3})",
                            subset_size, res.ratio, res.success_rate
                        );
                        continue;
                    }
                    if !validate_nl_narrowlane(self, &candidate_vars, &res, subset_size) {
                        tracing::debug!(
                            "AR: NL validation failed (subset={}, ratio={:.2}, sr={:.3})",
                            subset_size, res.ratio, res.success_rate
                        );
                        continue;
                    }
                    let fix_res =
                        self.apply_ar_fix(subset_size, &candidate_vars, &res, ephemerides)?;
                    let (fixed_state, da_meters, d_full) =
                        (fix_res.fixed_state, fix_res.z_dd, fix_res.d_full);
                    return Ok(crate::ambiguity::AmbiguityResolutionResult {
                        fixed_state,
                        z_dd: da_meters,
                        d_full,
                        ratio_test: res.ratio,
                        subset_size,
                    });
                }
            }
        }
        // Diagnostic: log why AR failed — what was the best ratio and which threshold blocked it
        let ffrt_at_best = crate::ffrt::calculate_threshold(best_subset_size, config.ar_ffrt_prob);
        let effective_threshold = ffrt_at_best.max(config.lambda_min_ratio);
        tracing::info!(
            "AR: LAMBDA failed ({} candidates, best ratio={:.2}, success_rate={:.3}, subset={}, ffrt={:.2}, threshold={:.2})",
            candidate_vars.len(), best_ratio, best_success_rate, best_subset_size, ffrt_at_best, effective_threshold
        );
        Err("AR failed to resolve")
    }

    fn apply_ar_fix(
        &self,
        subset_size: usize,
        candidate_vars: &[(usize, usize, u16, f64)],
        res: &crate::lambda::LambdaResult,
        ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    ) -> Result<crate::ambiguity::AmbiguityFixResult, &'static str> {
        let mut da_cycles = DVector::zeros(subset_size);
        let a_sd = nalgebra::DVector::from_vec(self.ambiguities.clone());
        let state_size = self.covariance.nrows();
        let mut d_full = DMatrix::zeros(subset_size, state_size);

        for row in 0..subset_size {
            let (rov, r_idx, _, _) = candidate_vars[row];
            let (rov_sat_id, freq_band) = self.ambiguity_keys[rov];
            let freq_num = ephemerides
                .iter()
                .find(|e| e.sat() == rov_sat_id)
                .map(|e| e.freq_num())
                .unwrap_or(0);
            let (f1_rov, f2_rov) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
            let lam_rov = gneiss_core::constants::SPEED_OF_LIGHT_M_S
                / if freq_band == 1 { f1_rov } else { f2_rov };

            let (ref_sat_id, _) = self.ambiguity_keys[r_idx];
            let ref_freq_num = ephemerides
                .iter()
                .find(|e| e.sat() == ref_sat_id)
                .map(|e| e.freq_num())
                .unwrap_or(0);
            let (f1_ref, f2_ref) =
                gneiss_core::signal::satellite_frequencies(ref_sat_id, ref_freq_num);
            let lam_ref = gneiss_core::constants::SPEED_OF_LIGHT_M_S
                / if freq_band == 1 { f1_ref } else { f2_ref };

            let a_cycle_float = a_sd[rov] / lam_rov - a_sd[r_idx] / lam_ref;
            da_cycles[row] = res.best_integers[row] - a_cycle_float;

            d_full[(row, CORE_STATE_SIZE + rov)] = 1.0 / lam_rov;
            d_full[(row, CORE_STATE_SIZE + r_idx)] = -1.0 / lam_ref;
        }

        let s = &d_full * &self.covariance * d_full.transpose();
        let s_inv = s.try_inverse().ok_or("Fix covariance inversion failed")?;
        let k_full = &self.covariance * d_full.transpose() * &s_inv;

        let dx = &k_full * &da_cycles;
        let mut fixed_state = self.clone();
        fixed_state.fixed_state = None;
        crate::engine::updater::apply_state_correction(&mut fixed_state, &dx);
        let r_zero = DMatrix::zeros(subset_size, subset_size);
        fixed_state.covariance = crate::math::covariance::apply_joseph_covariance_update(
            &self.covariance,
            &k_full,
            &d_full,
            &r_zero,
        );
        fixed_state.is_fixed = true;

        Ok(crate::ambiguity::AmbiguityFixResult {
            fixed_state,
            z_dd: da_cycles,
            d_full,
        })
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
    pub fn add_ambiguity(
        &mut self,
        sat: SatelliteId,
        freq: u8,
        initial_estimate: f64,
        initial_variance: f64,
    ) {
        if self.ambiguity_keys.contains(&(sat, freq)) {
            tracing::warn!(
                "Ambiguity for {:?} L{} already exists! Resetting.",
                sat,
                freq
            );
            self.remove_ambiguity(sat, freq);
        }
        tracing::debug!(
            "Adding ambiguity for {:?} L{} val={} var={}",
            sat,
            freq,
            initial_estimate,
            initial_variance
        );
        self.ambiguities.push(initial_estimate);
        self.ambiguity_keys.push((sat, freq));
        self.ambiguity_track_ids.push(self.next_track_id);
        self.next_track_id += 1;

        let n_old = self.covariance.nrows();
        self.covariance = self
            .covariance
            .clone()
            .insert_row(n_old, 0.0)
            .insert_column(n_old, 0.0);
        self.covariance[(n_old, n_old)] = initial_variance;
    }

    pub fn remove_ambiguity(&mut self, sat: SatelliteId, freq: u8) {
        if let Some(idx) = self
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == sat && f == freq)
        {
            self.ambiguities.remove(idx);
            self.ambiguity_keys.remove(idx);
            self.ambiguity_track_ids.remove(idx);

            let cov_idx = CORE_STATE_SIZE + idx;
            self.covariance = self
                .covariance
                .clone()
                .remove_row(cov_idx)
                .remove_column(cov_idx);

            self.gf_values.remove(&sat);
        }
    }

    pub fn clear_ambiguities(&mut self) {
        let num_amb = self.ambiguities.len();
        self.ambiguities.clear();
        self.ambiguity_keys.clear();
        self.ambiguity_track_ids.clear();

        self.covariance = self
            .covariance
            .clone()
            .remove_rows(CORE_STATE_SIZE, num_amb)
            .remove_columns(CORE_STATE_SIZE, num_amb);
    }
}

#[derive(Debug, Clone)]
pub struct DdObservation {
    pub sat: SatelliteId,
    pub pr_l1: f64,
    pub pr_l2: Option<f64>,
    pub cp_l1: Option<f64>,
    pub cp_l2: Option<f64>,
    pub doppler: f64,
    pub snr: f64,
    pub locktime: Option<u16>,
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
    rover_ref: &DdObservation,
    rover_sat: &DdObservation,
    base_ref: &DdObservation,
    base_sat: &DdObservation,
    ref_f1: f64,
    ref_f2: f64,
    sat_f1: f64,
    sat_f2: f64,
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
    use gneiss_core::coords::{Datum, Frame};
    use gneiss_core::sat::Constellation;

    #[test]
    fn test_resolve_ambiguities_multi_constellation() {
        let time = GpsTime::new(2137, 422922.0);
        let initial_pos = Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, initial_pos, 10.0);
        state.epoch_count = 100;

        let gps_ref = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let gps_rov1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        let gps_rov2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 3,
        };
        let gps_rov3 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 4,
        };
        let gal_ref = SatelliteId {
            constellation: Constellation::Galileo,
            prn: 10,
        };
        let gal_rov = SatelliteId {
            constellation: Constellation::Galileo,
            prn: 11,
        };

        let lam = 0.19029367279836487;
        state.add_ambiguity(gps_ref, 1, 10.1 * lam, 0.0001);
        state.add_ambiguity(gps_rov1, 1, 15.15 * lam, 0.0001);
        state.add_ambiguity(gps_rov2, 1, 20.05 * lam, 0.0001);
        state.add_ambiguity(gps_rov3, 1, 25.18 * lam, 0.0001);
        state.add_ambiguity(gal_ref, 1, 30.12 * lam, 0.0001);
        state.add_ambiguity(gal_rov, 1, 36.21 * lam, 0.0001);

        for &(sat, freq) in &[
            (gps_ref, 1),
            (gps_rov1, 1),
            (gps_rov2, 1),
            (gps_rov3, 1),
            (gal_ref, 1),
            (gal_rov, 1),
        ] {
            state.locktimes.insert((sat, freq), 100);
        }

        use gneiss_core::ephemeris::{Ephemeris, GalileoEphemeris, GpsEphemeris};
        let ephemerides = vec![
            Ephemeris::Gps(GpsEphemeris {
                sat: gps_ref,
                toe: time,
                toc: time,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.0,
                omega_dot: 0.0,
                i0: 1.0,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 0,
                iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris {
                sat: gps_rov1,
                toe: time,
                toc: time,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.1,
                omega_dot: 0.0,
                i0: 1.0,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 0,
                iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris {
                sat: gps_rov2,
                toe: time,
                toc: time,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.2,
                omega_dot: 0.0,
                i0: 1.0,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 0,
                iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris {
                sat: gps_rov3,
                toe: time,
                toc: time,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.3,
                omega_dot: 0.0,
                i0: 1.0,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 0,
                iodc: 0,
            }),
            Ephemeris::Galileo(GalileoEphemeris {
                sat: gal_ref,
                toe: time,
                toc: time,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5440.6,
                delta_n: 0.0,
                omega0: 0.0,
                omega_dot: 0.0,
                i0: 1.0,
                idot: 0.0,
                omega: 0.0,
                bgd_e1_e5a: 0.0,
                bgd_e1_e5b: 0.0,
                iod_nav: 0,
            }),
            Ephemeris::Galileo(GalileoEphemeris {
                sat: gal_rov,
                toe: time,
                toc: time,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5440.6,
                delta_n: 0.0,
                omega0: 0.5,
                omega_dot: 0.0,
                i0: 1.0,
                idot: 0.0,
                omega: 0.0,
                bgd_e1_e5a: 0.0,
                bgd_e1_e5b: 0.0,
                iod_nav: 0,
            }),
        ];

        let mut config = crate::engine::EngineConfig::default();
        config.lambda_min_subset = 4;
        config.ar_min_epoch_count = 5;
        config.ar_min_lock = 3;
        config.lambda_min_ratio = 3.0;
        config.ar_ffrt_prob = 0.001;
        config.tuning.min_ar_success_rate = 0.999;
        let res = state
            .resolve_ambiguities(&ephemerides, &config)
            .expect("AR should run");
        let fixed_state = res.fixed_state;
        assert!(
            fixed_state.is_fixed,
            "Should achieve fix with multi-constellation support"
        );

        let idx_ref = fixed_state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == gps_ref && f == 1)
            .unwrap();
        let idx_rov = fixed_state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == gps_rov1 && f == 1)
            .unwrap();
        let dd_gps = (fixed_state.ambiguities[idx_rov] - fixed_state.ambiguities[idx_ref]) / lam;
        assert!((dd_gps.round() - 5.0).abs() < 1e-6);

        let idx_ref_gal = fixed_state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == gal_ref && f == 1)
            .unwrap();
        let idx_rov_gal = fixed_state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == gal_rov && f == 1)
            .unwrap();
        let dd_gal =
            (fixed_state.ambiguities[idx_rov_gal] - fixed_state.ambiguities[idx_ref_gal]) / lam;
        assert!((dd_gal.round() - 6.0).abs() < 1e-6);
    }

    #[test]
    fn test_double_difference_eliminates_clocks() {
        let sat_ref = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat_a = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        let true_r_rover_ref = 20_000_000.0;
        let true_r_rover_a = 21_000_000.0;
        let true_r_base_ref = 20_005_000.0;
        let true_r_base_a = 21_004_000.0;
        let rover_clk = 300.0;
        let base_clk = -150.0;
        let sat_ref_clk = 1000.0;
        let sat_a_clk = -500.0;
        let rover_ref_obs = DdObservation {
            sat: sat_ref,
            pr_l1: true_r_rover_ref + rover_clk - sat_ref_clk,
            pr_l2: None,
            cp_l1: Some(0.0),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(1000),
        };
        let rover_a_obs = DdObservation {
            sat: sat_a,
            pr_l1: true_r_rover_a + rover_clk - sat_a_clk,
            pr_l2: None,
            cp_l1: Some(0.0),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(1000),
        };
        let base_ref_obs = DdObservation {
            sat: sat_ref,
            pr_l1: true_r_base_ref + base_clk - sat_ref_clk,
            pr_l2: None,
            cp_l1: Some(0.0),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(1000),
        };
        let base_a_obs = DdObservation {
            sat: sat_a,
            pr_l1: true_r_base_a + base_clk - sat_a_clk,
            pr_l2: None,
            cp_l1: Some(0.0),
            cp_l2: None,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(1000),
        };
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let dd = compute_double_difference(
            &rover_ref_obs,
            &rover_a_obs,
            &base_ref_obs,
            &base_a_obs,
            f1,
            f2,
            f1,
            f2,
        );
        assert!(
            (dd - 1000.0).abs() < 1e-6,
            "Double difference failed to eliminate clocks"
        );
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
        // Position variance inflated to >= 900, cross-terms untouched
        for i in 0..3 {
            assert!(state.covariance[(i, i)] >= 900.0);
            for j in 0..3 {
                if i != j {
                    assert_eq!(state.covariance[(i, j)], 1.0);
                }
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

        // Verify all 4 clock-related diagonals are reset to the initial variance
        for &idx in &clock_indices {
            assert_eq!(
                state.covariance[(idx, idx)],
                INITIAL_CLOCK_BIAS_VARIANCE,
                "diagonal at index {} should be {}",
                idx, INITIAL_CLOCK_BIAS_VARIANCE
            );
        }

        // Verify all cross-covariance terms involving clock indices are zeroed
        for &idx in &clock_indices {
            for j in 0..cols {
                if !clock_indices.contains(&j) {
                    assert_eq!(
                        state.covariance[(idx, j)],
                        0.0,
                        "cross-covariance ({}, {}) should be 0.0",
                        idx,
                        j
                    );
                    assert_eq!(
                        state.covariance[(j, idx)],
                        0.0,
                        "cross-covariance ({}, {}) should be 0.0",
                        j,
                        idx
                    );
                }
            }
        }

        // Verify non-clock entries are untouched (still 1.0)
        for i in 0..cols {
            for j in 0..cols {
                if !clock_indices.contains(&i) && !clock_indices.contains(&j) {
                    assert_eq!(
                        state.covariance[(i, j)],
                        1.0,
                        "non-clock entry ({}, {}) should remain 1.0",
                        i,
                        j
                    );
                }
            }
        }
    }

    #[test]
    fn test_update_mw() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };

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
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        let sat3 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 3,
        };

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
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };

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

        let gps_ref = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let gps_rov = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };

        // Add ambiguities
        state.add_ambiguity(gps_ref, 1, 0.0, 1.0);
        state.add_ambiguity(gps_rov, 1, 0.0, 1.0);

        // Set locktimes
        state.locktimes.insert((gps_ref, 1), 100);
        state.locktimes.insert((gps_rov, 1), 50);

        let ephemerides = vec![]; // Empty ephemerides is fine, will fallback to 0 freq_num

        // Test with lock limit 60
        let candidates_fail = crate::filter::select_ar_candidates(&state, &ephemerides, 60);
        assert_eq!(
            candidates_fail.len(),
            0,
            "Rover sat does not meet lock criteria"
        );

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

        let gps_ref = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let gps_rov = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };

        let lam = 0.19029367279836487;
        state.add_ambiguity(gps_ref, 1, 10.0 * lam, 0.1);
        state.add_ambiguity(gps_rov, 1, 15.0 * lam, 0.1);

        let candidates = vec![(1, 0, 50, 0.5)]; // rov_idx=1, ref_idx=0

        let ephemerides = vec![];
        let (d_mat, a_cycles, q_cycles) =
            crate::filter::build_lambda_matrices(&state, &candidates, 1, &ephemerides);

        assert_eq!(d_mat.nrows(), 1);
        assert_eq!(d_mat.ncols(), 2);
        assert!((d_mat[(0, 0)] - (-1.0 / lam)).abs() < 1e-6); // ref
        assert!((d_mat[(0, 1)] - (1.0 / lam)).abs() < 1e-6); // rov

        assert_eq!(a_cycles.len(), 1);
        assert!((a_cycles[0] - 5.0).abs() < 1e-6); // (15 / lam) - (10 / lam) = 5 / lam? No!
                                                   // Wait, a_sd is 15.0 * lam. So a_sd / lam = 15.0
                                                   // ref is 10.0 * lam. So a_sd / lam = 10.0
                                                   // 15.0 - 10.0 = 5.0

        assert_eq!(q_cycles.nrows(), 1);
        assert_eq!(q_cycles.ncols(), 1);
    }

    #[test]
    fn test_init_covariance_values() {
        let cov = init_covariance(10.0);
        assert_eq!(cov.nrows(), CORE_STATE_SIZE);
        assert_eq!(cov.ncols(), CORE_STATE_SIZE);
        for i in 0..3 {
            assert_eq!(cov[(i, i)], 10.0);
        }
        for i in 3..6 {
            assert_eq!(cov[(i, i)], 100.0);
        }
        for i in 6..9 {
            assert!((cov[(i, i)] - (1.0f64.to_radians()).powi(2)).abs() < 1e-10);
        }
        for i in 9..12 {
            assert!((cov[(i, i)] - 0.01).abs() < 1e-10);
        }
        for i in 12..15 {
            assert!((cov[(i, i)] - 1e-6).abs() < 1e-10);
        }
        assert!((cov[(15, 15)] - INITIAL_CLOCK_BIAS_VARIANCE).abs() < 1e-10);
        assert!((cov[(16, 16)] - INITIAL_ISB_VARIANCE).abs() < 1e-10);
        assert!((cov[(17, 17)] - INITIAL_ISB_VARIANCE).abs() < 1e-10);
        assert!((cov[(18, 18)] - INITIAL_ISB_VARIANCE).abs() < 1e-10);
        assert!((cov[(19, 19)] - 1000.0).abs() < 1e-10);
        assert!((cov[(20, 20)] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_remove_ambiguity() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 10.0, 1.0);
        assert_eq!(state.ambiguities.len(), 1);
        assert_eq!(state.covariance.nrows(), CORE_STATE_SIZE + 1);
        state.remove_ambiguity(sat, 1);
        assert_eq!(state.ambiguities.len(), 0);
        assert_eq!(state.covariance.nrows(), CORE_STATE_SIZE);
    }

    #[test]
    fn test_add_ambiguity_resets_existing() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 10.0, 1.0);
        assert_eq!(state.ambiguities[0], 10.0);
        // Add same ambiguity again -> should reset
        state.add_ambiguity(sat, 1, 20.0, 2.0);
        assert_eq!(state.ambiguities.len(), 1);
        assert_eq!(state.ambiguities[0], 20.0);
        assert_eq!(state.ambiguity_keys[0], (sat, 1));
    }

    #[test]
    fn test_reset_to_spp_clears_ambiguities() {
        use gneiss_core::coords::Datum;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat, 1, 10.0, 1.0);
        assert_eq!(state.ambiguities.len(), 1);

        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        state.reset_to_spp(new_pos, None, false);
        assert_eq!(state.ambiguities.len(), 0);
        assert_eq!(state.position.vector.y, 100.0);
        assert!(state.is_reset);
    }

    #[test]
    fn test_reset_to_spp_with_spp_state_and_isbs() {
        use gneiss_core::coords::Datum;
        use crate::spp::SppState;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 50.0), Datum::WGS84, Frame::ECEF, time);
        let spp_state = SppState::new(new_pos, 100.0, 110.0, 120.0, 130.0);
        state.reset_to_spp(new_pos, Some(&spp_state), true);
        assert!((state.rcv_clk_bias - 100.0).abs() < 1e-10);
        assert!((state.isb_glo - 30.0).abs() < 1e-10);
        assert!((state.isb_gal - 10.0).abs() < 1e-10);
        assert!((state.isb_bds - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_reset_to_spp_without_isbs() {
        use gneiss_core::coords::Datum;
        use crate::spp::SppState;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 50.0), Datum::WGS84, Frame::ECEF, time);
        let spp_state = SppState::new(new_pos, 100.0, 110.0, 120.0, 130.0);
        state.reset_to_spp(new_pos, Some(&spp_state), false);
        assert!((state.rcv_clk_bias - 100.0).abs() < 1e-10);
        assert!((state.isb_glo - 0.0).abs() < 1e-10);
        assert!((state.isb_gal - 0.0).abs() < 1e-10);
        assert!((state.isb_bds - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_init_attitude_produces_valid_quaternion() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let q = init_attitude(&pos);
        // Attitude should be a valid unit quaternion
        assert!((q.quaternion().w.powi(2) + q.quaternion().i.powi(2) + q.quaternion().j.powi(2) + q.quaternion().k.powi(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reset_to_spp_ins_aligned_true() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        state.ins_aligned = true;
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        state.reset_to_spp(new_pos, None, false);
        // When ins_aligned is true, attitude/IMU biases should NOT be reset
        assert!(state.ins_aligned);
        assert_eq!(state.position.vector.y, 100.0);
    }

    #[test]
    fn test_reset_to_spp_sanitizes_large_velocity() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        state.velocity = Vector3::new(500.0, 0.0, 0.0); // > 100 -> should be zeroed
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        state.reset_to_spp(new_pos, None, false);
        assert_eq!(state.velocity.norm(), 0.0);
    }

    #[test]
    fn test_reset_to_spp_sanitizes_nan_velocity() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        state.velocity = Vector3::new(f64::NAN, 0.0, 0.0);
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        state.reset_to_spp(new_pos, None, false);
        assert_eq!(state.velocity.norm(), 0.0);
    }

    #[test]
    fn test_reset_to_spp_sanitizes_large_clock_drift() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        state.rcv_clk_drift = 50000.0; // > 10000 -> should be zeroed
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        state.reset_to_spp(new_pos, None, false);
        assert_eq!(state.rcv_clk_drift, 0.0);
    }

    #[test]
    fn test_reset_to_spp_sanitizes_nan_clock_drift() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        state.rcv_clk_drift = f64::NAN;
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        state.reset_to_spp(new_pos, None, false);
        assert_eq!(state.rcv_clk_drift, 0.0);
    }

    #[test]
    fn test_reset_to_spp_covariance_reset_indices() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        // Set all covariances to 1.0 first
        state.covariance.fill(1.0);
        let new_pos = Coordinate::new(Vector3::new(6378137.0, 100.0, 50.0), Datum::WGS84, Frame::ECEF, time);
        state.reset_to_spp(new_pos, None, false);
        // Position diagonals should be reset to 100.0
        for i in 0..3 {
            assert!((state.covariance[(i, i)] - 100.0).abs() < 1e-10);
        }
        // Clock bias should be reset
        assert!((state.covariance[(15, 15)] - INITIAL_CLOCK_BIAS_VARIANCE).abs() < 1e-10);
    }

    #[test]
    fn test_find_best_reference_sat_with_reject_count() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        state.add_ambiguity(sat1, 1, 0.0, 1.0);
        state.add_ambiguity(sat2, 1, 0.0, 1.0);
        state.locktimes.insert((sat1, 1), 100);
        state.locktimes.insert((sat2, 1), 200);
        state.reject_counts.insert((sat2, 1), 1); // sat2 has rejections -> excluded

        let constell = Constellation::Gps;
        let result = find_best_reference_sat(&state, constell, 50);
        // sat1 should be chosen (sat2 has reject count > 0)
        assert_eq!(result, Some(0));
    }

    #[test]
    fn test_find_best_reference_sat_all_rejected() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.add_ambiguity(sat1, 1, 0.0, 1.0);
        state.locktimes.insert((sat1, 1), 100);
        state.reject_counts.insert((sat1, 1), 3); // rejected

        let result = find_best_reference_sat(&state, Constellation::Gps, 50);
        assert_eq!(result, None);
    }

    #[test]
    fn test_collect_candidates_l2_reference() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let sat_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_rov = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        // Add L1 and L2 for ref, L1 for rov
        state.add_ambiguity(sat_ref, 1, 0.0, 1.0);
        state.add_ambiguity(sat_ref, 2, 0.0, 1.0);
        state.add_ambiguity(sat_rov, 1, 0.0, 1.0);
        state.add_ambiguity(sat_rov, 2, 0.0, 1.0);
        // Set locktimes so all pass filter
        for &(s, f) in &[(sat_ref, 1), (sat_ref, 2), (sat_rov, 1), (sat_rov, 2)] {
            state.locktimes.insert((s, f), 100);
        }

        let mut candidates = Vec::new();
        // ref_idx = 0 (sat_ref, L1)
        collect_candidates_for_constellation(&state, Constellation::Gps, 0, 50, &mut candidates);
        // Should have: (rov L1 -> ref L1), (rov L2 -> ref L2)
        // Excludes ref itself (idx 0) and ref L2 companion (idx 1)
        assert_eq!(candidates.len(), 2);
        // Candidate 0 should be rov L1 with ref L1 (idx 2 -> idx 0)
        assert_eq!(candidates[0], (2, 0, 100));
        // Candidate 1 should be rov L2 with ref L2 (idx 3 -> idx 1 because l2_ref_idx=1)
        assert_eq!(candidates[1], (3, 1, 100));
    }

    #[test]
    fn test_collect_candidates_rejected_satellite() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let sat_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_rov = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        state.add_ambiguity(sat_ref, 1, 0.0, 1.0);
        state.add_ambiguity(sat_rov, 1, 0.0, 1.0);
        state.locktimes.insert((sat_ref, 1), 100);
        state.locktimes.insert((sat_rov, 1), 100);
        state.reject_counts.insert((sat_rov, 1), 2); // rejected

        let mut candidates = Vec::new();
        collect_candidates_for_constellation(&state, Constellation::Gps, 0, 50, &mut candidates);
        assert_eq!(candidates.len(), 0, "Rejected satellite should be excluded");
    }

    #[test]
    fn test_compute_candidate_variance_clamping() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let sat_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_rov = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        state.add_ambiguity(sat_ref, 1, 0.0, 1.0);
        state.add_ambiguity(sat_rov, 1, 0.0, 1.0);
        // Set a very high variance on the ambiguity covariance elements to trigger clamping
        let amb_idx = CORE_STATE_SIZE;
        state.covariance[(amb_idx, amb_idx)] = 1e8; // Very large variance
        state.covariance[(amb_idx + 1, amb_idx + 1)] = 1e8;
        state.covariance[(amb_idx + 1, amb_idx)] = 0.0;
        state.covariance[(amb_idx, amb_idx + 1)] = 0.0;

        let candidates = vec![(1, 0, 100)];
        let result = compute_candidate_variance(&state, &[], &candidates);
        // With empty ephemerides, lam_rov = lam_ref = LIGHT_SPEED / 1575.42e6 ~ 0.19
        // var_cycles = q_sd[1,1]/(0.19^2) + q_sd[0,0]/(0.19^2) = 2 * 1e8 / 0.0361 ~ 5.5e9 > 10000
        // So the candidate should be excluded by clamping
        assert_eq!(result.len(), 0, "High variance candidate should be excluded");
    }

    #[test]
    fn test_build_lambda_variance_matrix_structure() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        let sat_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_rov1 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let sat_rov2 = SatelliteId { constellation: Constellation::Gps, prn: 3 };
        state.add_ambiguity(sat_ref, 1, 10.0, 0.1);
        state.add_ambiguity(sat_rov1, 1, 15.0, 0.1);
        state.add_ambiguity(sat_rov2, 1, 20.0, 0.1);
        let candidates = vec![
            (1, 0, 100, 0.5),
            (2, 0, 50, 0.8),
        ];
        let q = build_lambda_variance_matrix(&state, &candidates, 2, &[]);
        assert_eq!(q.nrows(), 2);
        assert_eq!(q.ncols(), 2);
        // The matrix should be symmetric
        assert!((q[(0, 1)] - q[(1, 0)]).abs() < 1e-10);
    }

    #[test]
    fn test_filter_by_locktime_empty_constellations() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 10.0);
        // Empty constellations slice
        let result = filter_by_locktime(&state, &[], 50);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_resolve_ambiguities_insufficient_data() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 10.0);
        let config = crate::engine::EngineConfig::default();
        let result = state.resolve_ambiguities(&[], &config);
        assert!(result.is_err(), "Expected error for insufficient data");
    }

    #[test]
    fn test_resolve_ambiguities_insufficient_epochs() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 10.0);
        state.epoch_count = 2; // Below ar_min_epoch_count (usually 5)
        state.add_ambiguity(SatelliteId { constellation: Constellation::Gps, prn: 1 }, 1, 0.0, 1.0);
        state.add_ambiguity(SatelliteId { constellation: Constellation::Gps, prn: 2 }, 1, 0.0, 1.0);
        let mut config = crate::engine::EngineConfig::default();
        config.lambda_min_subset = 2;
        config.ar_min_epoch_count = 5;
        let result = state.resolve_ambiguities(&[], &config);
        assert!(result.is_err(), "Expected error for insufficient epochs");
    }
}
pub fn select_ar_candidates(
    state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    ar_min_lock: u32,
) -> Vec<(usize, usize, u16, f64)> {
    let constellations = [
        gneiss_core::sat::Constellation::Gps,
        gneiss_core::sat::Constellation::Galileo,
        gneiss_core::sat::Constellation::Beidou,
        gneiss_core::sat::Constellation::Glonass,
    ];

    // Diagnostic: log AR vs EKF reference mismatch per constellation
    for &constell in &constellations {
        let ekf_ref = state.current_ref_sat.get(&constell);
        if let Some(ar_ref_idx) = find_best_reference_sat(state, constell, ar_min_lock) {
            let ar_ref_sat = state.ambiguity_keys[ar_ref_idx].0;
            if ekf_ref != Some(&ar_ref_sat) {
                tracing::debug!(
                    "AR ref mismatch {:?}: EKF={:?} AR={:?}",
                    constell,
                    ekf_ref,
                    ar_ref_sat,
                );
            }
        }
    }

    let candidates = filter_by_locktime(state, &constellations, ar_min_lock);
    let mut candidate_vars = compute_candidate_variance(state, ephemerides, &candidates);
    candidate_vars.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));

    // Diagnostic: log top-5 candidate DD variances
    if tracing::enabled!(tracing::Level::DEBUG) && !candidate_vars.is_empty() {
        let top_n = candidate_vars.len().min(5);
        let vars_cycles: Vec<String> = candidate_vars[..top_n]
            .iter()
            .map(|(rov, _, _, var_m2)| {
                let sat = state.ambiguity_keys[*rov].0;
                let lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S
                    / gneiss_core::signal::satellite_frequencies(sat, 0).0;
                let var_cyc = var_m2 / (lam * lam);
                format!("{:?}={:.1}", sat, var_cyc)
            })
            .collect();
        tracing::debug!(
            "AR candidates: {} total, top var (cyc²): [{}]",
            candidate_vars.len(),
            vars_cycles.join(", ")
        );
    }

    candidate_vars
}

fn filter_by_locktime(
    state: &RtkState,
    constellations: &[gneiss_core::sat::Constellation],
    ar_min_lock: u32,
) -> Vec<(usize, usize, u16)> {
    let mut candidates = Vec::new();
    for &constell in constellations {
        if let Some(ref_idx) = find_best_reference_sat(state, constell, ar_min_lock) {
            collect_candidates_for_constellation(
                state,
                constell,
                ref_idx,
                ar_min_lock,
                &mut candidates,
            );
        }
    }
    candidates
}

fn find_best_reference_sat(
    state: &RtkState,
    constell: gneiss_core::sat::Constellation,
    ar_min_lock: u32,
) -> Option<usize> {
    let mut best_ref_idx = None;
    let mut max_lock = 0;
    for i in 0..state.ambiguities.len() {
        let (sat, freq) = state.ambiguity_keys[i];
        if sat.constellation != constell || freq != 1 {
            continue;
        }
        let reject_count = *state.reject_counts.get(&(sat, freq)).unwrap_or(&0);
        if reject_count > 0 {
            continue;
        }
        let lock = *state.locktimes.get(&(sat, freq)).unwrap_or(&0);
        if lock >= ar_min_lock as u16 && lock > max_lock {
            max_lock = lock;
            best_ref_idx = Some(i);
        }
    }
    best_ref_idx
}

/// Validate LAMBDA integer fix against the geometry-free Melbourne-Wübbena
/// widelane combination.  Returns false if any satellite with sufficient MW
/// data disagrees with the LAMBDA integers by more than `mw_tolerance` cycles.
fn validate_mw_widelane(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    lambda_result: &crate::ambiguity::lambda::LambdaResult,
    subset_size: usize,
) -> bool {
    let mw_tolerance: f64 = 0.3; // cycles of widelane (~86cm wavelength)
    let subset = &candidates[..subset_size];

    // Map L1 DD integers by rover satellite index for quick lookup
    // lambda_result.best_integers[i] corresponds to subset[i]
    let l1_ints: std::collections::HashMap<usize, f64> = subset
        .iter()
        .enumerate()
        .filter(|(_, &(rov, _, _, _))| state.ambiguity_keys[rov].1 == 1)
        .map(|(i, &(rov, _, _, _))| (rov, lambda_result.best_integers[i]))
        .collect();

    let l2_ints: std::collections::HashMap<usize, f64> = subset
        .iter()
        .enumerate()
        .filter(|(_, &(rov, _, _, _))| state.ambiguity_keys[rov].1 == 2)
        .map(|(i, &(rov, _, _, _))| (rov, lambda_result.best_integers[i]))
        .collect();

    for &(rov_idx, ref_idx, _, _) in subset {
        let (rov_sat, freq) = state.ambiguity_keys[rov_idx];
        if freq != 1 {
            continue; // only validate L1 candidates (they pair with L2)
        }

        // Find L2 candidate for same satellite
        if let Some(l2_rov_idx) = state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == rov_sat && f == 2)
        {
            let ref_sat = state.ambiguity_keys[ref_idx].0;

            if let (Some(&n1_dd), Some(&n2_dd)) =
                (l1_ints.get(&rov_idx), l2_ints.get(&l2_rov_idx))
            {
                let wl_lambda = n1_dd - n2_dd; // widelane from LAMBDA (cycles)

                // Wideline from MW EMA
                let mw_rov = state.mw_sd_ema.get(&rov_sat).copied().unwrap_or(0.0);
                let mw_ref = state.mw_sd_ema.get(&ref_sat).copied().unwrap_or(0.0);
                let mw_counts_rov = state.mw_sd_counts.get(&rov_sat).copied().unwrap_or(0);
                let mw_counts_ref = state.mw_sd_counts.get(&ref_sat).copied().unwrap_or(0);

                // Both rover and reference must have sufficient MW data;
                // otherwise mw_ref defaults to 0.0, introducing a systematic
                // offset equal to the reference satellite's true SD MW value
                // in EVERY DD comparison, rejecting all fixes.
                if mw_counts_rov >= 10 && mw_counts_ref >= 10 {
                    let wl_mw = mw_rov - mw_ref; // DD widelane from MW (cycles)
                    let diff = (wl_lambda - wl_mw).abs();
                    if diff > mw_tolerance {
                        tracing::debug!(
                            "MW validation failed for {:?}: LAMBDA WL={:.2} MW WL={:.2} diff={:.2} (rov_cnt={}, ref_cnt={})",
                            rov_sat, wl_lambda, wl_mw, diff, mw_counts_rov, mw_counts_ref
                        );
                        return false;
                    }
                }
            }
        }
    }
    true
}

/// Validate LAMBDA narrow-lane integers against time-averaged code-based NL.
/// Returns false if any satellite disagrees by more than `nl_tolerance` cycles.
fn validate_nl_narrowlane(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    lambda_result: &crate::ambiguity::lambda::LambdaResult,
    subset_size: usize,
) -> bool {
    let nl_tolerance: f64 = 0.5; // cycles of narrow-lane (~10.7cm wavelength)
    let subset = &candidates[..subset_size];

    // LAMBDA best_integers[i] = DD integer for subset[i] (already rover - reference)
    // Build maps: rov_idx -> DD integer
    let l1_dd_ints: std::collections::HashMap<usize, f64> = subset
        .iter()
        .enumerate()
        .filter(|(_, &(rov, _, _, _))| state.ambiguity_keys[rov].1 == 1)
        .map(|(i, &(rov, _, _, _))| (rov, lambda_result.best_integers[i]))
        .collect();

    let l2_dd_ints: std::collections::HashMap<usize, f64> = subset
        .iter()
        .enumerate()
        .filter(|(_, &(rov, _, _, _))| state.ambiguity_keys[rov].1 == 2)
        .map(|(i, &(rov, _, _, _))| (rov, lambda_result.best_integers[i]))
        .collect();

    // Use GPS frequencies for NL combination (same ratio works for all GNSS)
    let f1 = gneiss_core::signal::FREQ_GPS_L1; // 1575.42 MHz
    let f2 = gneiss_core::signal::FREQ_GPS_L2; // 1227.60 MHz

    for &(rov_idx, _, _, _) in subset {
        let (rov_sat, freq) = state.ambiguity_keys[rov_idx];
        if freq != 1 {
            continue;
        }
        // Find L2 DD integer for same satellite
        if let Some(l2_rov_idx) = state
            .ambiguity_keys
            .iter()
            .position(|&(s, f)| s == rov_sat && f == 2)
        {
            if let (Some(&n1_dd), Some(&n2_dd)) =
                (l1_dd_ints.get(&rov_idx), l2_dd_ints.get(&l2_rov_idx))
            {
                // DD narrow-lane from LAMBDA integers (cycles of NL)
                let nl_dd = (f1 * n1_dd + f2 * n2_dd) / (f1 + f2);

                // DD narrow-lane from code EMA
                let nl_rov = state.nl_sd_ema.get(&rov_sat).copied().unwrap_or(0.0);
                let nl_ref_sat = state.ambiguity_keys[subset.iter()
                    .find(|&&(r, _, _, _)| r == rov_idx)
                    .map(|&(_, r, _, _)| r)
                    .unwrap_or(0)].0;
                let nl_ref = state.nl_sd_ema.get(&nl_ref_sat).copied().unwrap_or(0.0);
                let nl_counts_rov = state.nl_sd_counts.get(&rov_sat).copied().unwrap_or(0);
                let nl_counts_ref = state.nl_sd_counts.get(&nl_ref_sat).copied().unwrap_or(0);

                if nl_counts_rov >= 10 && nl_counts_ref >= 10 {
                    let nl_ema_dd = nl_rov - nl_ref;
                    if (nl_dd - nl_ema_dd).abs() > nl_tolerance {
                        tracing::debug!(
                            "NL validation failed for {:?}: LAMBDA NL={:.2} EMA NL={:.2} diff={:.2}",
                            rov_sat, nl_dd, nl_ema_dd, (nl_dd - nl_ema_dd).abs()
                        );
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn collect_candidates_for_constellation(
    state: &RtkState,
    constell: gneiss_core::sat::Constellation,
    ref_idx: usize,
    ar_min_lock: u32,
    candidates: &mut Vec<(usize, usize, u16)>,
) {
    let ref_sat_id = state.ambiguity_keys[ref_idx].0;
    let l2_ref_idx = state
        .ambiguity_keys
        .iter()
        .position(|&(s, f)| s == ref_sat_id && f == 2);
    for i in 0..state.ambiguities.len() {
        if i == ref_idx || Some(i) == l2_ref_idx {
            continue;
        }
        let (rov_sat, freq) = state.ambiguity_keys[i];
        if rov_sat.constellation != constell {
            continue;
        }
        let reject_count = *state.reject_counts.get(&(rov_sat, freq)).unwrap_or(&0);
        if reject_count > 0 {
            continue;
        }
        let lock = *state.locktimes.get(&(rov_sat, freq)).unwrap_or(&0);
        if lock >= ar_min_lock as u16 {
            if freq == 1 {
                candidates.push((i, ref_idx, lock));
            } else if let Some(r2_idx) = l2_ref_idx {
                candidates.push((i, r2_idx, lock));
            }
        }
    }
}

fn compute_candidate_variance(
    state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    candidates: &[(usize, usize, u16)],
) -> Vec<(usize, usize, u16, f64)> {
    let num_amb = state.ambiguities.len();
    let q_sd = state
        .covariance
        .view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb));
    let mut candidate_vars = Vec::new();
    for &(rov, r_idx, lock) in candidates {
        let (rov_sat_id, freq_band) = state.ambiguity_keys[rov];
        let freq_num = ephemerides
            .iter()
            .find(|e| e.sat() == rov_sat_id)
            .map(|e| e.freq_num())
            .unwrap_or(0);
        let (f1_rov, f2_rov) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
        let lam_rov = gneiss_core::constants::SPEED_OF_LIGHT_M_S
            / if freq_band == 1 { f1_rov } else { f2_rov };

        let (ref_sat_id, _) = state.ambiguity_keys[r_idx];
        let ref_freq_num = ephemerides
            .iter()
            .find(|e| e.sat() == ref_sat_id)
            .map(|e| e.freq_num())
            .unwrap_or(0);
        let (f1_ref, f2_ref) = gneiss_core::signal::satellite_frequencies(ref_sat_id, ref_freq_num);
        let lam_ref = gneiss_core::constants::SPEED_OF_LIGHT_M_S
            / if freq_band == 1 { f1_ref } else { f2_ref };

        let var_cycles = q_sd[(rov, rov)] / (lam_rov * lam_rov)
            + q_sd[(r_idx, r_idx)] / (lam_ref * lam_ref)
            - 2.0 * q_sd[(rov, r_idx)] / (lam_rov * lam_ref);

        if var_cycles < 10000.0 {
            candidate_vars.push((rov, r_idx, lock, var_cycles));
        }
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
        let freq_num = ephemerides
            .iter()
            .find(|e| e.sat() == rov_sat_id)
            .map(|e| e.freq_num())
            .unwrap_or(0);
        let (f1_rov, f2_rov) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
        let lam_rov = gneiss_core::constants::SPEED_OF_LIGHT_M_S
            / if freq_band == 1 { f1_rov } else { f2_rov };

        let (ref_sat_id, _) = state.ambiguity_keys[r_idx];
        let ref_freq_num = ephemerides
            .iter()
            .find(|e| e.sat() == ref_sat_id)
            .map(|e| e.freq_num())
            .unwrap_or(0);
        let (f1_ref, f2_ref) = gneiss_core::signal::satellite_frequencies(ref_sat_id, ref_freq_num);
        let lam_ref = gneiss_core::constants::SPEED_OF_LIGHT_M_S
            / if freq_band == 1 { f1_ref } else { f2_ref };

        d_mat[(row, rov)] = 1.0 / lam_rov;
        d_mat[(row, r_idx)] = -1.0 / lam_ref;
        a_cycles[row] = a_sd[rov] / lam_rov - a_sd[r_idx] / lam_ref;
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
    let q_sd = state
        .covariance
        .view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb));
    let mut q_cycles = DMatrix::zeros(subset_size, subset_size);
    for r in 0..subset_size {
        for c in 0..subset_size {
            let (rov_r, ref_r, _, _) = candidates[r];
            let (rov_c, ref_c, _, _) = candidates[c];
            let freq_r = state.ambiguity_keys[rov_r].1;
            let freq_c = state.ambiguity_keys[rov_c].1;

            let (sat_rov_r, _) = state.ambiguity_keys[rov_r];
            let freq_num_rov_r = ephemerides
                .iter()
                .find(|e| e.sat() == sat_rov_r)
                .map(|e| e.freq_num())
                .unwrap_or(0);
            let (f1_rov_r, f2_rov_r) =
                gneiss_core::signal::satellite_frequencies(sat_rov_r, freq_num_rov_r);
            let lam_rov_r = gneiss_core::constants::SPEED_OF_LIGHT_M_S
                / if freq_r == 1 { f1_rov_r } else { f2_rov_r };

            let (sat_ref_r, _) = state.ambiguity_keys[ref_r];
            let freq_num_ref_r = ephemerides
                .iter()
                .find(|e| e.sat() == sat_ref_r)
                .map(|e| e.freq_num())
                .unwrap_or(0);
            let (f1_ref_r, f2_ref_r) =
                gneiss_core::signal::satellite_frequencies(sat_ref_r, freq_num_ref_r);
            let lam_ref_r = gneiss_core::constants::SPEED_OF_LIGHT_M_S
                / if freq_r == 1 { f1_ref_r } else { f2_ref_r };

            let (sat_rov_c, _) = state.ambiguity_keys[rov_c];
            let freq_num_rov_c = ephemerides
                .iter()
                .find(|e| e.sat() == sat_rov_c)
                .map(|e| e.freq_num())
                .unwrap_or(0);
            let (f1_rov_c, f2_rov_c) =
                gneiss_core::signal::satellite_frequencies(sat_rov_c, freq_num_rov_c);
            let lam_rov_c = gneiss_core::constants::SPEED_OF_LIGHT_M_S
                / if freq_c == 1 { f1_rov_c } else { f2_rov_c };

            let (sat_ref_c, _) = state.ambiguity_keys[ref_c];
            let freq_num_ref_c = ephemerides
                .iter()
                .find(|e| e.sat() == sat_ref_c)
                .map(|e| e.freq_num())
                .unwrap_or(0);
            let (f1_ref_c, f2_ref_c) =
                gneiss_core::signal::satellite_frequencies(sat_ref_c, freq_num_ref_c);
            let lam_ref_c = gneiss_core::constants::SPEED_OF_LIGHT_M_S
                / if freq_c == 1 { f1_ref_c } else { f2_ref_c };

            let val = q_sd[(rov_r, rov_c)] / (lam_rov_r * lam_rov_c)
                - q_sd[(rov_r, ref_c)] / (lam_rov_r * lam_ref_c)
                - q_sd[(ref_r, rov_c)] / (lam_ref_r * lam_rov_c)
                + q_sd[(ref_r, ref_c)] / (lam_ref_r * lam_ref_c);
            q_cycles[(r, c)] = val;
        }
    }
    q_cycles
}
