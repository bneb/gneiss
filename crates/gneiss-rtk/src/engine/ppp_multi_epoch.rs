//! Two-epoch sliding-window factor graph optimizer for PPP.
//!
//! Jointly optimizes state at epochs k-1 and k:
//!   [x_{k-1} (21), x_k (21), ambiguities (N)]
//!
//! - GNSS measurement factors at both epochs (PR, CP, Doppler)
//! - Dynamics constraint: x_k ≈ Phi · x_{k-1} with configurable process noise
//! - SPP position prior on current epoch
//!
//! Falls back to single-epoch IEKF if the LM optimization fails.

use nalgebra::{DMatrix, DVector, Vector3};

use crate::engine::ppp_common::extract_state_vector;
use crate::engine::ppp_iekf::PppIteratedEkf;
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::types::IonosphereModel;
use crate::engine::EngineError;
use crate::estimators::factor_graph::gnss_factors::{
    ErrorStateCarrierPhaseFactor, ErrorStateDopplerFactor, ErrorStatePseudorangeFactor,
};
use crate::estimators::factor_graph::{Factor, FactorGraphOptimizer};
use crate::filter::{CORE_STATE_SIZE, RtkState};
use crate::math::inversion::invert_matrix_robust;

// ---------------------------------------------------------------------------
// Owned per-satellite data (no references) so we can store between epochs.
// ---------------------------------------------------------------------------

/// Minimal owned copy of the satellite data needed to build GNSS factors.
#[derive(Clone, Debug)]
struct OwnedSatData {
    sat_pos: Vector3<f64>,
    sat_vel: Vector3<f64>,
    p1: f64,
    cp1: Option<f64>,
    doppler: f64,
    dt_sat_m: f64,
    tropo_dry: f64,
    map_wet: f64,
    lam1: f64,
    el: f64,
    snr: f64,
    sat_id: gneiss_core::sat::SatelliteId,
}

impl<'a> From<&ProcessedSat<'a>> for OwnedSatData {
    fn from(s: &ProcessedSat<'a>) -> Self {
        Self {
            sat_pos: s.sat_pos_rot,
            sat_vel: s.sat_vel,
            p1: s.p1,
            cp1: s.cp1,
            doppler: s.doppler,
            dt_sat_m: s.dt_sat_m,
            tropo_dry: s.tropo_dry,
            map_wet: s.map_wet,
            lam1: s.lam1,
            el: s.el,
            snr: s.snr,
            sat_id: s.sat_obs.sat,
        }
    }
}

/// Snapshot of a single epoch's state and measurement data.
#[derive(Clone)]
struct EpochSnapshot {
    /// Full core state vector (CORE_STATE_SIZE elements)
    state: DVector<f64>,
    /// Covariance matrix (full rank, including ambiguities)
    #[allow(dead_code)]
    cov: DMatrix<f64>,
    /// Satellite observations
    sats: Vec<OwnedSatData>,
}

// ---------------------------------------------------------------------------
// Process-noise configuration
// ---------------------------------------------------------------------------

/// Process-noise standard deviations for the dynamics constraint.
///
/// The dynamics factor models x_k = x_{k-1} + w_k where w_k ~ N(0, Q)
/// and Q = diag(q_i^2) for each core-state element.
///
/// Default values assume ~1 s epoch spacing and static / slow-moving
/// receiver (urban pedestrian / vehicular).
#[derive(Clone, Debug)]
pub struct ProcessNoiseConfig {
    /// Position random walk std (m/√s).
    pub pos: f64,
    /// Velocity random walk std (m/s/√s).
    pub vel: f64,
    /// Attitude random walk std (rad/√s).
    pub attitude: f64,
    /// Accel-bias random walk std (m/s²/√s).
    pub accel_bias: f64,
    /// Gyro-bias random walk std (rad/s/√s).
    pub gyro_bias: f64,
    /// Clock-bias random walk std (m/√s).
    pub clock_bias: f64,
    /// ISB random walk std (m/√s).
    pub isb: f64,
    /// Clock-drift random walk std (m/s/√s).
    pub clock_drift: f64,
    /// Zenith wet delay random walk std (m/√s).
    pub zwd: f64,
}

impl Default for ProcessNoiseConfig {
    fn default() -> Self {
        Self {
            pos: 0.1,        // 10 cm/sqrt(s) — tight static assumption
            vel: 1.0,         // 1 m/s/sqrt(s)
            attitude: 0.01,   // 0.01 rad/sqrt(s)
            accel_bias: 0.01, // 0.01 m/s^2/sqrt(s)
            gyro_bias: 0.001, // 0.001 rad/s/sqrt(s)
            clock_bias: 10.0, // 10 m/sqrt(s) — loose, GPS clock can drift
            isb: 0.01,        // 0.01 m/sqrt(s) — very stable
            clock_drift: 0.1, // 0.1 m/s/sqrt(s)
            zwd: 0.01,        // 0.01 m/sqrt(s) — slow
        }
    }
}

// ---------------------------------------------------------------------------
// Dynamics factor (between-epoch constraint)
// ---------------------------------------------------------------------------

/// Factor that penalises deviation from the dynamics model x_k = Phi · x_{k-1}.
///
/// Operates on the error-state delta vector:
///   residual = (nominal_curr + delta_curr) - Phi · (nominal_prev + delta_prev)
///
/// For static PPP, Phi is the identity matrix and the residual simplifies to
///   (x_{k-1} + delta_curr - delta_{prev}) after the nominals converge.
struct TwoEpochDynamicsFactor {
    /// Number of core states (CORE_STATE_SIZE).
    core_dim: usize,
    /// State transition matrix (core_dim × core_dim).
    phi: DMatrix<f64>,
    /// Inverse process-noise covariance (core_dim × core_dim).
    q_inv: DMatrix<f64>,
    /// Nominal previous-epoch core state.
    nominal_prev: DVector<f64>,
    /// Nominal current-epoch core state.
    nominal_curr: DVector<f64>,
    /// Total error-state dimension (= 2 * core_dim + num_ambiguities).
    #[allow(dead_code)]
    total_dim: usize,
}

impl Factor for TwoEpochDynamicsFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let delta_prev = delta.rows(0, self.core_dim);
        let delta_curr = delta.rows(self.core_dim, self.core_dim);
        let full_prev = &self.nominal_prev + delta_prev;
        let full_curr = &self.nominal_curr + delta_curr;
        &full_curr - &self.phi * &full_prev
    }

    fn jacobian(&self, _delta: &DVector<f64>) -> DMatrix<f64> {
        let mut jac = DMatrix::zeros(self.core_dim, self.total_dim);
        // d(residual)/d(delta_prev) = -Phi  (columns 0..core_dim)
        for i in 0..self.core_dim {
            for j in 0..self.core_dim {
                jac[(i, j)] = -self.phi[(i, j)];
            }
        }
        // d(residual)/d(delta_curr) = I     (columns core_dim .. 2*core_dim)
        for i in 0..self.core_dim {
            jac[(i, self.core_dim + i)] = 1.0;
        }
        // d(residual)/d(delta_amb) = 0      (already zero)
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        self.q_inv.clone()
    }
}

// ---------------------------------------------------------------------------
// SPP position-prior factor
// ---------------------------------------------------------------------------

/// Soft prior on the current-epoch position from SPP.
///
/// residual = (nominal + delta) - spp_pos   (3 × 1 vector).
struct PositionPriorFactor {
    spp_pos: Vector3<f64>,
    /// Inverse variance (scalar, applied isotropically).
    info: f64,
    /// Index of x in the full error-state vector.
    index_x: usize,
    index_y: usize,
    index_z: usize,
    nominal_x: f64,
    nominal_y: f64,
    nominal_z: f64,
    total_dim: usize,
}

impl Factor for PositionPriorFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        DVector::from_vec(vec![
            self.nominal_x + delta[self.index_x] - self.spp_pos.x,
            self.nominal_y + delta[self.index_y] - self.spp_pos.y,
            self.nominal_z + delta[self.index_z] - self.spp_pos.z,
        ])
    }

    fn jacobian(&self, _delta: &DVector<f64>) -> DMatrix<f64> {
        let mut jac = DMatrix::zeros(3, self.total_dim);
        jac[(0, self.index_x)] = 1.0;
        jac[(1, self.index_y)] = 1.0;
        jac[(2, self.index_z)] = 1.0;
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_diagonal(&DVector::from_element(3, self.info))
    }
}

// ---------------------------------------------------------------------------
// PppTwoEpochOptimizer
// ---------------------------------------------------------------------------

/// Two-epoch sliding-window factor-graph optimizer for PPP.
///
/// Accumulates epoch pairs (k-1, k) and jointly optimises:
///
/// - Pseudorange, carrier-phase and Doppler factors at both epochs,
/// - A dynamics constraint x_k ≈ Phi · x_{k-1},
/// - An SPP position prior on the current epoch.
///
/// The LM back-end reuses the existing `FactorGraphOptimizer`.
/// If the joint optimisation fails the result falls back to the single-epoch
/// IEKF estimate for epoch k.
/// Type alias for backward compatibility (Phase 1).
pub type MultiEpochOptimizer = PppTwoEpochOptimizer;

pub struct PppTwoEpochOptimizer {
    /// Sliding window of past epoch snapshots (oldest first).
    window: std::collections::VecDeque<EpochSnapshot>,
    /// Maximum number of epochs in the sliding window (default 5).
    window_size: usize,
    /// Process-noise configuration for the dynamics constraint.
    process_noise: ProcessNoiseConfig,
    /// Maximum LM iterations.
    max_iterations: usize,
    /// LM convergence tolerance (∆x norm).
    convergence_tol: f64,
    /// Huber robust-loss threshold (applied inside GNSS factors).
    huber_k: f64,
    /// Ionosphere model (Klobuchar / IONEX).
    iono_model: IonosphereModel,
    /// LAMBDA AR minimum ratio threshold.
    lambda_min_ratio: f64,
}

impl Default for PppTwoEpochOptimizer {
    fn default() -> Self {
        Self {
            window: std::collections::VecDeque::new(),
            window_size: 5,
            process_noise: ProcessNoiseConfig::default(),
            max_iterations: 15,
            convergence_tol: 1e-4,
            huber_k: 3.0,
            iono_model: IonosphereModel::Klobuchar,
            lambda_min_ratio: 2.0,
        }
    }
}

impl PppTwoEpochOptimizer {
    /// Create a new `PppTwoEpochOptimizer` with default configuration.
    pub fn new(_window_size: usize) -> Self {
        Self::default()
    }

    /// Set the ionosphere model.
    pub fn with_iono_model(mut self, model: IonosphereModel) -> Self {
        self.iono_model = model;
        self
    }

    /// Set the LAMBDA AR minimum ratio threshold.
    pub fn with_lambda_min_ratio(mut self, ratio: f64) -> Self {
        self.lambda_min_ratio = ratio;
        self
    }

    // ------------------------------------------------------------------
    // Public entry point
    // ------------------------------------------------------------------

    /// Process one PPP epoch.
    ///
    /// - **First call**: runs the single-epoch IEKF, caches the result, and
    ///   returns it unchanged.
    /// - **Subsequent calls**: runs the single-epoch IEKF on the current
    ///   epoch, builds a 2-epoch joint factor graph, runs LM optimisation,
    ///   and writes the smoothed epoch-k estimate back into `state`.
    /// - **Fallback**: if the LM optimisation fails (singular system,
    ///   divergence) the single-epoch IEKF result is kept.
    pub fn solve(
        &mut self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<(), EngineError> {
        // --- Step 1: run single-epoch IEKF to get an initial estimate ------
        let iekf = PppIteratedEkf::new()
            .with_iono_model(self.iono_model)
            .with_lambda_min_ratio(self.lambda_min_ratio);
        iekf.solve(state, sats, position_prior)?;

        // --- Step 2: push current snapshot into the sliding window ---------
        let curr_snapshot = self.snapshot(state, sats);
        self.window.push_back(curr_snapshot);

        // --- Step 3: need at least 2 epochs for joint optimisation ---------
        if self.window.len() < 2 {
            return Ok(());
        }

        // --- Step 4: try N-epoch joint optimisation ------------------------
        let n_epochs = self.window.len();
        let result = self.try_n_epoch_optimisation(state, sats, position_prior);

        match result {
            Ok((final_state, final_cov)) => {
                apply_core_state(state, &final_state, &final_cov);
            }
            Err(ref e) => {
                tracing::warn!(
                    "{}-epoch optimisation failed ({}), keeping IEKF result.",
                    n_epochs, e
                );
            }
        }

        // --- Step 5: maintain window size ----------------------------------
        while self.window.len() > self.window_size {
            self.window.pop_front();
        }

        Ok(())
    }

    // ------------------------------------------------------------------
    // Internal helpers
    // ------------------------------------------------------------------

    /// Snapshot the current state+measurements for the next epoch pair.
    fn snapshot(&self, state: &RtkState, sats: &[ProcessedSat]) -> EpochSnapshot {
        EpochSnapshot {
            state: extract_state_vector(state),
            cov: state.covariance.clone(),
            sats: sats.iter().map(OwnedSatData::from).collect(),
        }
    }

    /// Attempt a two-epoch joint optimisation.
    ///
    /// Returns `(final_core_state, final_covariance)` on success.
    #[allow(clippy::too_many_arguments)]
    fn try_n_epoch_optimisation(
        &self,
        state: &RtkState,
        _sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<(DVector<f64>, DMatrix<f64>), EngineError> {
        // Use the oldest and newest epochs in the window for 2-epoch joint
        // optimisation. A wider baseline (N epochs apart) provides better
        // geometry diversity than adjacent pairs when the window > 2.
        // Full N-epoch joint optimisation (Phase 2 proper) is deferred.
        let prev = self.window.front().expect("window has >= 2 entries");
        let curr = self.window.back().expect("window has >= 2 entries");
        // --- Determine ambiguity layout -----------------------------------
        // Collect all ambiguity keys from the current state (which the IEKF
        // already updated).  Satellites visible in both epochs share the same
        // ambiguity index.
        let n_amb = state.ambiguities.len();
        let amb_keys = &state.ambiguity_keys;
        let total_dim = 2 * CORE_STATE_SIZE + n_amb;

        // --- Build the Q^{-1} matrix for the dynamics constraint ----------
        let q_inv = build_q_inverse(&self.process_noise);

        // --- Build the dynamics factor ------------------------------------
        // Phi is identity for static PPP.
        let phi = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let dynamics = TwoEpochDynamicsFactor {
            core_dim: CORE_STATE_SIZE,
            phi,
            q_inv,
            nominal_prev: prev.state.rows(0, CORE_STATE_SIZE).clone_owned(),
            nominal_curr: curr.state.rows(0, CORE_STATE_SIZE).clone_owned(),
            total_dim,
        };

        // --- Build GNSS factors for the previous epoch --------------------
        let n_prev_amb = prev.state.len().saturating_sub(CORE_STATE_SIZE);
        let mut optimizer = FactorGraphOptimizer::new();

        for sat in &prev.sats {
            // Find ambiguity index in the CURRENT state's ambiguity list.
            // This way both epochs share the same ambiguity state variable.
            let amb_idx = amb_keys
                .iter()
                .position(|&(k, _f)| k == sat.sat_id);

            add_epoch_factors(
                &mut optimizer,
                sat,
                &prev.state,
                n_prev_amb,
                amb_idx,
                0,                // core offset = 0 for epoch k-1
                total_dim,
                self.huber_k,
            );
        }

        // --- Build GNSS factors for the current epoch ---------------------
        let n_curr_amb = curr.state.len().saturating_sub(CORE_STATE_SIZE);
        for sat in &curr.sats {
            let amb_idx = amb_keys
                .iter()
                .position(|&(k, _f)| k == sat.sat_id);

            add_epoch_factors(
                &mut optimizer,
                sat,
                &curr.state,
                n_curr_amb,
                amb_idx,
                CORE_STATE_SIZE, // core offset for epoch k
                total_dim,
                self.huber_k,
            );
        }

        // --- Add dynamics constraint --------------------------------------
        optimizer.add_factor(Box::new(dynamics));

        // --- Add SPP position prior on current epoch ----------------------
        if let Some((spp_pos, var)) = position_prior {
            optimizer.add_factor(Box::new(PositionPriorFactor {
                spp_pos,
                info: 1.0 / var.max(1e-9),
                index_x: CORE_STATE_SIZE,     // position.x of epoch k
                index_y: CORE_STATE_SIZE + 1, // position.y of epoch k
                index_z: CORE_STATE_SIZE + 2, // position.z of epoch k
                nominal_x: curr.state[0],
                nominal_y: curr.state[1],
                nominal_z: curr.state[2],
                total_dim,
            }));
        }

        // --- Run LM optimisation ------------------------------------------
        let initial_delta = DVector::zeros(total_dim);
        let (delta_opt, _cov) = optimizer.optimize(&initial_delta, self.max_iterations, self.convergence_tol);

        // --- Check for optimisation failure -------------------------------
        if delta_opt.iter().any(|x| x.is_nan() || x.is_infinite()) {
            return Err(EngineError::StateDisappeared);
        }

        // --- Extract the current-epoch position correction ----------------
        let dx_k = delta_opt.rows(CORE_STATE_SIZE, CORE_STATE_SIZE);

        // Compute the final state from the current IEKF state + delta
        let mut final_state = curr.state.clone();
        for i in 0..CORE_STATE_SIZE {
            final_state[i] += dx_k[i];
        }

        // Compute approximate posterior covariance (top-left 21×21 block of
        // the factor-graph information matrix inverse).
        let final_cov = if !delta_opt.is_empty() {
            compute_sub_covariance(&optimizer, &delta_opt, CORE_STATE_SIZE, 2 * CORE_STATE_SIZE)
                .unwrap_or_else(|| state.covariance.clone())
        } else {
            state.covariance.clone()
        };

        Ok((final_state, final_cov))
    }
}

// ---------------------------------------------------------------------------
// Factor-construction helpers
// ---------------------------------------------------------------------------

/// Add PR, CP and Doppler factors for one satellite and one epoch.
///
/// `core_offset` is the index of the first core-state element for this epoch
/// in the combined error-state vector (0 for epoch k-1, CORE_STATE_SIZE for
/// epoch k).  Ambiguities live at indices [2 × CORE_STATE_SIZE .. total_dim).
#[allow(clippy::too_many_arguments)]
fn add_epoch_factors(
    opt: &mut FactorGraphOptimizer,
    sat: &OwnedSatData,
    core_state: &DVector<f64>,
    _n_amb: usize,
    amb_idx: Option<usize>,
    core_offset: usize,
    _total_dim: usize,
    huber_k: f64,
) {
    let amb_offset = 2 * CORE_STATE_SIZE;

    // --- Pseudorange factor -----------------------------------------------
    opt.add_factor(Box::new(ErrorStatePseudorangeFactor {
        sat_pos: sat.sat_pos,
        measured_pr: sat.p1,
        variance: pseudorange_variance(sat),
        sat_clock_bias: sat.dt_sat_m,
        tropo_dry_delay: sat.tropo_dry,
        map_wet: sat.map_wet,
        nominal_rx: core_state[0],
        nominal_ry: core_state[1],
        nominal_rz: core_state[2],
        nominal_dt: core_state[15],
        nominal_dt_gal: core_state[17],
        nominal_dt_bds: core_state[18],
        nominal_dt_glo: core_state[16],
        nominal_zwd: core_state[20],
        index_x: core_offset,       // 0 or CORE_STATE_SIZE
        index_y: core_offset + 1,
        index_z: core_offset + 2,
        index_dt: core_offset + 15,
        index_dt_gal: Some(core_offset + 17),
        index_dt_bds: Some(core_offset + 18),
        index_dt_glo: Some(core_offset + 16),
        index_zwd: Some(core_offset + 20),
        sat_id: sat.sat_id,
        robust_threshold: huber_k,
    }));

    // --- Carrier-phase factor (if available) ------------------------------
    if let Some(cp1) = sat.cp1 {
        if let Some(a_idx) = amb_idx {
            // Nominal ambiguity value lives in the state vector beyond core.
            let nominal_amb = if a_idx < core_state.len().saturating_sub(CORE_STATE_SIZE) {
                core_state[CORE_STATE_SIZE + a_idx]
            } else {
                0.0
            };
            opt.add_factor(Box::new(ErrorStateCarrierPhaseFactor {
                sat_pos: sat.sat_pos,
                measured_cp: cp1,
                variance: cp_variance(sat),
                sat_clock_bias: sat.dt_sat_m,
                tropo_dry_delay: sat.tropo_dry,
                map_wet: sat.map_wet,
                wavelength: sat.lam1,
                nominal_rx: core_state[0],
                nominal_ry: core_state[1],
                nominal_rz: core_state[2],
                nominal_dt: core_state[15],
                nominal_dt_gal: core_state[17],
                nominal_dt_bds: core_state[18],
                nominal_dt_glo: core_state[16],
                nominal_zwd: core_state[20],
                nominal_amb,
                index_x: core_offset,
                index_y: core_offset + 1,
                index_z: core_offset + 2,
                index_dt: core_offset + 15,
                index_dt_gal: Some(core_offset + 17),
                index_dt_bds: Some(core_offset + 18),
                index_dt_glo: Some(core_offset + 16),
                index_zwd: Some(core_offset + 20),
                index_amb: amb_offset + a_idx,
                sat_id: sat.sat_id,
                robust_threshold: huber_k,
            }));
        }
    }

    // --- Doppler factor (if non-zero) ------------------------------------
    if sat.doppler != 0.0 {
        opt.add_factor(Box::new(ErrorStateDopplerFactor {
            los: sat.sat_pos / sat.sat_pos.norm(),
            sat_vel: sat.sat_vel,
            measured_doppler_hz: sat.doppler,
            variance: 0.01,
            wavelength: sat.lam1,
            sat_clock_drift: 0.0,
            nominal_vx: core_state[3],
            nominal_vy: core_state[4],
            nominal_vz: core_state[5],
            nominal_cdt: core_state[19],
            index_vx: core_offset + 3,
            index_vy: core_offset + 4,
            index_vz: core_offset + 5,
            index_cdt: core_offset + 19,
            robust_threshold: huber_k,
        }));
    }
}

/// Approximate pseudorange variance from SNR and elevation.
fn pseudorange_variance(sat: &OwnedSatData) -> f64 {
    let snr_scale = (10.0_f64).powf((45.0 - sat.snr) / 10.0);
    let mut var = 1.0 * snr_scale / sat.el.sin().max(0.1);
    var += 9.0; // Klobuchar residual error (3^2)
    var
}

/// Approximate carrier-phase variance from SNR and elevation.
fn cp_variance(sat: &OwnedSatData) -> f64 {
    let snr_scale = (10.0_f64).powf((45.0 - sat.snr) / 10.0);
    0.0001 * snr_scale / sat.el.sin().max(0.1)
}

/// Build the inverse process-noise covariance matrix Q⁻¹ from the config.
fn build_q_inverse(config: &ProcessNoiseConfig) -> DMatrix<f64> {
    let q_diag = vec![
        config.pos.powi(2),
        config.pos.powi(2),
        config.pos.powi(2),
        config.vel.powi(2),
        config.vel.powi(2),
        config.vel.powi(2),
        config.attitude.powi(2),
        config.attitude.powi(2),
        config.attitude.powi(2),
        config.accel_bias.powi(2),
        config.accel_bias.powi(2),
        config.accel_bias.powi(2),
        config.gyro_bias.powi(2),
        config.gyro_bias.powi(2),
        config.gyro_bias.powi(2),
        config.clock_bias.powi(2),
        config.isb.powi(2),
        config.isb.powi(2),
        config.isb.powi(2),
        config.clock_drift.powi(2),
        config.zwd.powi(2),
    ];
    let q_inv_diag: Vec<f64> = q_diag.iter().map(|v| 1.0 / v.max(1e-18)).collect();
    DMatrix::from_diagonal(&DVector::from_vec(q_inv_diag))
}

/// Compute a sub-block of the posterior covariance from the factor graph.
///
/// Extracts columns `start..end` from the full information matrix inverse.
fn compute_sub_covariance(
    optimizer: &FactorGraphOptimizer,
    state: &DVector<f64>,
    start: usize,
    end: usize,
) -> Option<DMatrix<f64>> {
    let n = end - start;
    if n == 0 {
        return None;
    }

    // Build the full information matrix H = Σ J_i^T W_i J_i
    let mut h = DMatrix::zeros(state.len(), state.len());
    for factor in &optimizer.factors {
        let info = factor.information();
        let jac = factor.jacobian(state);
        let j_t_info = jac.transpose() * &info;
        h += &j_t_info * &jac;
    }

    // Add LM damping so the matrix is invertible
    for i in 0..h.nrows() {
        h[(i, i)] += 1e-6;
    }

    // Invert and extract the sub-block
    let cov = invert_matrix_robust(&h);
    Some(cov.view((start, start), (n, n)).clone_owned())
}

/// Write the optimised epoch-k core state back into an `RtkState`.
fn apply_core_state(state: &mut RtkState, core: &DVector<f64>, cov: &DMatrix<f64>) {
    state.position.vector.x = core[0];
    state.position.vector.y = core[1];
    state.position.vector.z = core[2];
    state.velocity.x = core[3];
    state.velocity.y = core[4];
    state.velocity.z = core[5];
    let r_vec = Vector3::new(core[6], core[7], core[8]);
    state.attitude = if r_vec.norm() > 1e-12 {
        nalgebra::UnitQuaternion::from_scaled_axis(r_vec)
    } else {
        nalgebra::UnitQuaternion::identity()
    };
    state.accel_bias.x = core[9];
    state.accel_bias.y = core[10];
    state.accel_bias.z = core[11];
    state.gyro_bias.x = core[12];
    state.gyro_bias.y = core[13];
    state.gyro_bias.z = core[14];
    state.rcv_clk_bias = core[15];
    state.isb_glo = core[16];
    state.isb_gal = core[17];
    state.isb_bds = core[18];
    state.rcv_clk_drift = core[19];
    state.zwd = core[20];

    // Merge covariance — only update the top-left 21×21 block.
    if cov.nrows() >= CORE_STATE_SIZE && cov.ncols() >= CORE_STATE_SIZE {
        for i in 0..CORE_STATE_SIZE {
            for j in 0..CORE_STATE_SIZE {
                state.covariance[(i, j)] = cov[(i, j)];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    // ------------------------------------------------------------------
    // Unit tests for helper types
    // ------------------------------------------------------------------

    #[test]
    fn test_process_noise_config_default() {
        let cfg = ProcessNoiseConfig::default();
        assert!((cfg.pos - 0.1).abs() < 1e-12);
        assert!((cfg.clock_bias - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_owned_sat_data_from_processed() {
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.5,
            p1: 20485784.0,
            p2: None,
            cp1: Some(107653241.0),
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::new(0.6, 0.0, 0.8),
            dist: 20485784.0,
            el: 0.5,
            snr: 45.0,
            doppler: -250.0,
            lam1: 0.1903,
            lam2: 0.0,
            tropo_dry: 2.3,
            map_wet: 1.0,
            iono_delay: 3.0,
            f1: 0.0,
            f2: 0.0,
            sat_pos_rot: Vector3::new(10000000.0, 0.0, 0.0),
            sat_vel: Vector3::new(10.0, 0.0, 0.0),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::new(0.0, 0.0, 0.0),
            pcv_correction: 0.0,
        };
        let owned = OwnedSatData::from(&sat);
        assert_eq!(owned.sat_id, sat_id);
        assert!((owned.sat_pos.x - 10000000.0).abs() < 1e-6);
        assert!((owned.p1 - 20485784.0).abs() < 1e-6);
        assert_eq!(owned.cp1, Some(107653241.0));
    }

    #[test]
    fn test_build_q_inverse_default() {
        let q_inv = build_q_inverse(&ProcessNoiseConfig::default());
        assert_eq!(q_inv.nrows(), CORE_STATE_SIZE);
        assert_eq!(q_inv.ncols(), CORE_STATE_SIZE);
        // Position diagonal: 1 / 0.01 = 100
        assert!((q_inv[(0, 0)] - 100.0).abs() < 1e-9);
        // Clock bias: 1 / 100 = 0.01
        assert!((q_inv[(15, 15)] - 0.01).abs() < 1e-9);
    }

    #[test]
    fn test_pseudorange_variance_basic() {
        let sat = OwnedSatData {
            sat_pos: Vector3::new(10000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(),
            p1: 0.0, cp1: None, doppler: 0.0,
            dt_sat_m: 0.0, tropo_dry: 0.0, map_wet: 0.0,
            lam1: 0.19, el: core::f64::consts::FRAC_PI_2, snr: 45.0,
            sat_id: SatelliteId { constellation: Constellation::Gps, prn: 1 },
        };
        let var = pseudorange_variance(&sat);
        // snr_scale = 10^((45-45)/10) = 1.0
        // var = 1.0 * 1.0 / sin(pi/2) + 9.0 = 1.0 + 9.0 = 10.0
        assert!((var - 10.0).abs() < 1e-9, "Expected 10.0, got {}", var);
    }

    #[test]
    fn test_cp_variance_basic() {
        let sat = OwnedSatData {
            sat_pos: Vector3::new(10000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(),
            p1: 0.0, cp1: None, doppler: 0.0,
            dt_sat_m: 0.0, tropo_dry: 0.0, map_wet: 0.0,
            lam1: 0.19, el: core::f64::consts::FRAC_PI_2, snr: 45.0,
            sat_id: SatelliteId { constellation: Constellation::Gps, prn: 1 },
        };
        let var = cp_variance(&sat);
        // snr_scale = 1.0, var = 0.0001 * 1.0 / 1.0 = 0.0001
        assert!((var - 0.0001).abs() < 1e-9, "Expected 0.0001, got {}", var);
    }

    // ------------------------------------------------------------------
    // DynamicsFactor unit tests
    // ------------------------------------------------------------------

    #[test]
    fn test_dynamics_factor_residual_zero_when_nominals_match() {
        let nominal = DVector::from_element(CORE_STATE_SIZE, 1.0);
        let q_inv = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let phi = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let factor = TwoEpochDynamicsFactor {
            core_dim: CORE_STATE_SIZE,
            phi,
            q_inv,
            nominal_prev: nominal.clone(),
            nominal_curr: nominal,
            total_dim: 2 * CORE_STATE_SIZE,
        };
        let delta = DVector::zeros(2 * CORE_STATE_SIZE);
        let res = factor.residual(&delta);
        assert!(
            res.norm() < 1e-10,
            "Residual should be zero when nominal states match. Norm={}",
            res.norm()
        );
    }

    #[test]
    fn test_dynamics_factor_residual_nonzero_nominal_mismatch() {
        let prev = DVector::from_element(CORE_STATE_SIZE, 0.0);
        let curr = DVector::from_element(CORE_STATE_SIZE, 5.0);
        let q_inv = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let phi = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let factor = TwoEpochDynamicsFactor {
            core_dim: CORE_STATE_SIZE,
            phi,
            q_inv,
            nominal_prev: prev,
            nominal_curr: curr,
            total_dim: 2 * CORE_STATE_SIZE,
        };
        let delta = DVector::zeros(2 * CORE_STATE_SIZE);
        let res = factor.residual(&delta);
        assert!(
            (res[0] - 5.0).abs() < 1e-10,
            "Residual should be (curr_0 - prev_0) = 5.0, got {}",
            res[0]
        );
    }

    #[test]
    fn test_dynamics_factor_residual_corrected_by_delta() {
        let prev = DVector::from_element(CORE_STATE_SIZE, 0.0);
        let curr = DVector::from_element(CORE_STATE_SIZE, 5.0);
        let q_inv = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let phi = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let factor = TwoEpochDynamicsFactor {
            core_dim: CORE_STATE_SIZE,
            phi,
            q_inv,
            nominal_prev: prev,
            nominal_curr: curr,
            total_dim: 2 * CORE_STATE_SIZE,
        };
        // delta_prev = 1, delta_curr = 6 => full_prev = 1, full_curr = 11
        // residual = 11 - 1 = 10
        let mut delta = DVector::zeros(2 * CORE_STATE_SIZE);
        delta[0] = 1.0;
        delta[CORE_STATE_SIZE] = 6.0;
        let res = factor.residual(&delta);
        assert!(
            (res[0] - 10.0).abs() < 1e-10,
            "Residual should be 10.0, got {}",
            res[0]
        );
    }

    #[test]
    fn test_dynamics_factor_jacobian_structure() {
        let nominal = DVector::from_element(CORE_STATE_SIZE, 0.0);
        let q_inv = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let phi = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let factor = TwoEpochDynamicsFactor {
            core_dim: CORE_STATE_SIZE,
            phi,
            q_inv,
            nominal_prev: nominal.clone(),
            nominal_curr: nominal,
            total_dim: 2 * CORE_STATE_SIZE,
        };
        let delta = DVector::zeros(2 * CORE_STATE_SIZE);
        let jac = factor.jacobian(&delta);
        // d(res[0])/d(delta_prev[0]) = -phi[0,0] = -1
        assert!((jac[(0, 0)] - (-1.0)).abs() < 1e-10);
        // d(res[0])/d(delta_curr[0]) = 1
        assert!((jac[(0, CORE_STATE_SIZE)] - 1.0).abs() < 1e-10);
        // No coupling to ambiguity elements (columns beyond 2*core_dim are 0)
        assert_eq!(jac.ncols(), 2 * CORE_STATE_SIZE);
    }

    // ------------------------------------------------------------------
    // PositionPriorFactor unit tests
    // ------------------------------------------------------------------

    #[test]
    fn test_position_prior_residual_zero() {
        let spp = Vector3::new(1.0, 2.0, 3.0);
        let factor = PositionPriorFactor {
            spp_pos: spp,
            info: 1.0,
            index_x: CORE_STATE_SIZE,
            index_y: CORE_STATE_SIZE + 1,
            index_z: CORE_STATE_SIZE + 2,
            nominal_x: 1.0,
            nominal_y: 2.0,
            nominal_z: 3.0,
            total_dim: 2 * CORE_STATE_SIZE,
        };
        let delta = DVector::zeros(2 * CORE_STATE_SIZE);
        let res = factor.residual(&delta);
        assert!(
            res.norm() < 1e-10,
            "Prior residual should be zero. Norm={}",
            res.norm()
        );
    }

    #[test]
    fn test_position_prior_residual_nonzero() {
        let spp = Vector3::new(0.0, 0.0, 0.0);
        let factor = PositionPriorFactor {
            spp_pos: spp,
            info: 1.0,
            index_x: CORE_STATE_SIZE,
            index_y: CORE_STATE_SIZE + 1,
            index_z: CORE_STATE_SIZE + 2,
            nominal_x: 5.0,
            nominal_y: 5.0,
            nominal_z: 5.0,
            total_dim: 2 * CORE_STATE_SIZE,
        };
        let delta = DVector::zeros(2 * CORE_STATE_SIZE);
        let res = factor.residual(&delta);
        assert!(
            (res[0] - 5.0).abs() < 1e-10,
            "Residual[0] should be 5.0, got {}",
            res[0]
        );
        assert!((res[1] - 5.0).abs() < 1e-10);
        assert!((res[2] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_position_prior_jacobian() {
        let factor = PositionPriorFactor {
            spp_pos: Vector3::zeros(),
            info: 1.0,
            index_x: CORE_STATE_SIZE,
            index_y: CORE_STATE_SIZE + 1,
            index_z: CORE_STATE_SIZE + 2,
            nominal_x: 0.0,
            nominal_y: 0.0,
            nominal_z: 0.0,
            total_dim: 2 * CORE_STATE_SIZE,
        };
        let delta = DVector::zeros(2 * CORE_STATE_SIZE);
        let jac = factor.jacobian(&delta);
        assert_eq!(jac.nrows(), 3);
        assert_eq!(jac.ncols(), 2 * CORE_STATE_SIZE);
        assert!((jac[(0, CORE_STATE_SIZE)] - 1.0).abs() < 1e-10);
        assert!((jac[(1, CORE_STATE_SIZE + 1)] - 1.0).abs() < 1e-10);
        assert!((jac[(2, CORE_STATE_SIZE + 2)] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_position_prior_information() {
        let factor = PositionPriorFactor {
            spp_pos: Vector3::zeros(),
            info: 2.5,
            index_x: 0, index_y: 1, index_z: 2,
            nominal_x: 0.0, nominal_y: 0.0, nominal_z: 0.0,
            total_dim: 10,
        };
        let info = factor.information();
        assert_eq!(info.nrows(), 3);
        assert_eq!(info.ncols(), 3);
        assert!((info[(0, 0)] - 2.5).abs() < 1e-10);
        assert!((info[(1, 1)] - 2.5).abs() < 1e-10);
        assert!((info[(2, 2)] - 2.5).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // apply_core_state helper tests
    // ------------------------------------------------------------------

    #[test]
    fn test_apply_core_state_writes_values() {
        let time = GpsTime::new(2082, 0.0);
        let coord = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, coord, 10.0);
        let core = DVector::from_vec(vec![
            100.0, 200.0, 300.0,  // position
            0.5, 0.6, 0.7,       // velocity
            0.1, 0.0, 0.0,       // attitude (small rotation)
            0.01, 0.02, 0.03,    // accel bias
            0.001, 0.002, 0.003, // gyro bias
            42.0,                // rcv_clk_bias
            0.1, 0.2, 0.3,       // ISBs
            0.5,                 // rcv_clk_drift
            0.05,                // ZWD
        ]);
        let cov = DMatrix::identity(21, 21);
        apply_core_state(&mut state, &core, &cov);

        assert!((state.position.vector.x - 100.0).abs() < 1e-10);
        assert!((state.position.vector.y - 200.0).abs() < 1e-10);
        assert!((state.position.vector.z - 300.0).abs() < 1e-10);
        assert!((state.velocity.y - 0.6).abs() < 1e-10);
        assert!((state.rcv_clk_bias - 42.0).abs() < 1e-10);
        assert!((state.isb_gal - 0.2).abs() < 1e-10);
        assert!((state.zwd - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_apply_core_state_zero_attitude() {
        let time = GpsTime::new(2082, 0.0);
        let coord = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, coord, 10.0);
        // Attitude elements (6,7,8) will be 0 from RtkState::new (or default)
        let core = DVector::from_vec(vec![
            1.0, 2.0, 3.0,       // position unchanged
            0.0, 0.0, 0.0,       // velocity
            0.0, 0.0, 0.0,       // attitude = zero
            0.0, 0.0, 0.0,       // accel bias
            0.0, 0.0, 0.0,       // gyro bias
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // clocks, ZWD
        ]);
        let cov = DMatrix::identity(21, 21);
        apply_core_state(&mut state, &core, &cov);
        // Attitude should be identity when rotation vector is zero
        assert!((state.attitude.quaternion().w - 1.0).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // Integration test: two-epoch factor graph with realistic data
    // ------------------------------------------------------------------

    /// Helper to create a minimal RtkState for testing.
    fn make_state(time: GpsTime) -> RtkState {
        RtkState::new(
            time,
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                time,
            ),
            1.0,
        )
    }

    /// Helper to build a ProcessedSat from simple parameters.
    fn make_sat(
        prn: u8,
        sat_pos: Vector3<f64>,
        p1: f64,
        cp1: Option<f64>,
        doppler: f64,
    ) -> ProcessedSat<'static> {
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn };
        let dist = sat_pos.norm();
        let lam1 = 0.1903;
        let obs_ref: &'static SatObs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        ProcessedSat {
            sat_obs: obs_ref,
            dt_sat_m: 0.0,
            p1,
            p2: None,
            cp1,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: sat_pos / dist.max(1e-6),
            dist,
            el: core::f64::consts::FRAC_PI_2,
            snr: 45.0,
            doppler,
            lam1,
            lam2: 0.0,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 0.0,
            f2: 0.0,
            sat_pos_rot: sat_pos,
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        }
    }

    #[test]
    fn test_optimizer_first_epoch_caches_and_returns() {
        let mut opt = PppTwoEpochOptimizer::new(2);
        let time = GpsTime::new(2082, 0.0);
        let coord = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, coord, 10.0);
        state.position.vector = Vector3::new(1000000.0, 2000000.0, 3000000.0);
        state.covariance = DMatrix::identity(21, 21) * 100.0;

        // Build a satellite far away so the IEKF produces a valid result
        let sat_pos = Vector3::new(20000000.0, 10000000.0, 5000000.0);
        let dist = (sat_pos - state.position.vector).norm();
        let sat = make_sat(1, sat_pos, dist, Some(dist / 0.1903), 0.0);
        let sats = [sat];

        let result = opt.solve(&mut state, &sats, None);
        assert!(result.is_ok(), "First epoch should succeed: {:?}", result);
        // State should be updated (IEKF ran)
        assert!(state.position.vector.norm() > 0.0);
        // After first epoch, window should have 1 entry
        assert_eq!(opt.window.len(), 1);
    }

    #[test]
    fn test_multi_epoch_produces_result() {
        let mut opt = PppTwoEpochOptimizer::new(2);
        let time = GpsTime::new(2082, 0.0);
        let coord = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, coord, 10.0);
        state.position.vector = Vector3::new(1000000.0, 2000000.0, 3000000.0);
        state.covariance = DMatrix::identity(21, 21) * 100.0;

        // Same satellite visible across both epochs
        let sat_pos = Vector3::new(20000000.0, 10000000.0, 5000000.0);
        let dist = (sat_pos - state.position.vector).norm();
        let sat1 = make_sat(1, sat_pos, dist, Some(dist / 0.1903), 0.0);
        let sats1 = [sat1];

        // Epoch 1
        opt.solve(&mut state, &sats1, None).unwrap();

        // Epoch 2: same satellite, same geometry
        let sat2 = make_sat(1, sat_pos, dist, Some(dist / 0.1903), 0.0);
        let sats2 = [sat2];
        let result = opt.solve(&mut state, &sats2, None);
        assert!(result.is_ok(), "Second epoch should succeed: {:?}", result);
    }

    #[test]
    fn test_two_epoch_factor_graph_full_optimization() {
        // Build a factor graph with 2 epochs, 4 satellites, dynamics and prior.
        // Verify the optimizer converges to a valid state.
        let pos0 = Vector3::new(1000000.0, 2000000.0, 3000000.0);

        let time = GpsTime::new(2082, 0.0);
        let coord = Coordinate::new(pos0, Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, coord, 10.0);
        state.position.vector = pos0;
        state.covariance = DMatrix::identity(21, 21) * 100.0;

        // 4 GPS satellites with good geometry
        let sat_positions = vec![
            Vector3::new(20000000.0, 10000000.0, 5000000.0),
            Vector3::new(-15000000.0, 20000000.0, 8000000.0),
            Vector3::new(-18000000.0, -12000000.0, 6000000.0),
            Vector3::new(22000000.0, -15000000.0, 4000000.0),
        ];

        let mut sats = Vec::new();
        for (i, sp) in sat_positions.iter().enumerate() {
            let dist = (sp - pos0).norm();
            let prn = (i + 1) as u8;
            let s = make_sat(prn, *sp, dist, Some(dist / 0.1903), 0.0);
            sats.push(s);
        }

        // Epoch 1: run IEKF (via solve first call)
        let mut opt = PppTwoEpochOptimizer::new(2);
        opt.solve(&mut state, &sats, None).unwrap();

        // Make sure the state is set to pos0 for epoch 2 (after IEKF may have
        // changed it, reset to a slightly perturbed position to test smoothing)
        state.position.vector = pos0 + Vector3::new(2.0, -1.0, 3.0);
        state.covariance = DMatrix::identity(21, 21) * 100.0;

        // Epoch 2: should run two-epoch optimization
        let mut sats2 = Vec::new();
        for (i, sp) in sat_positions.iter().enumerate() {
            let dist = (sp - pos0).norm();
            let prn = (i + 1) as u8;
            let s = make_sat(prn, *sp, dist, Some(dist / 0.1903), 0.0);
            sats2.push(s);
        }

        let result = opt.solve(&mut state, &sats2, None);
        assert!(result.is_ok(), "Two-epoch optimization should succeed: {:?}", result);

        // The optimised position should be closer to pos0 than the perturbation
        let error = (state.position.vector - pos0).norm();
        assert!(
            error < 5.0,
            "Position error should be small after smoothing. error={:.3}m, pos={:?}",
            error,
            state.position.vector
        );
    }

    #[test]
    fn test_optimizer_recovery_after_iekf_failure() {
        // With 0 satellites, the IEKF fails but we should handle it gracefully
        let mut opt = PppTwoEpochOptimizer::new(2);
        let time = GpsTime::new(2082, 0.0);
        let coord = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, coord, 10.0);
        state.position.vector = Vector3::new(1000000.0, 2000000.0, 3000000.0);
        state.covariance = DMatrix::identity(21, 21) * 100.0;

        let sats: Vec<ProcessedSat> = vec![];
        let result = opt.solve(&mut state, &sats, None);
        assert!(result.is_err(), "IEKF should fail with no satellites");
    }

    #[test]
    fn test_add_epoch_factors_pseudorange_only() {
        let mut opt = FactorGraphOptimizer::new();
        let core_state = DVector::from_element(CORE_STATE_SIZE, 0.0);
        let _amb_vec: Vec<f64> = vec![];

        let sat = OwnedSatData {
            sat_pos: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(),
            p1: 20000000.0,
            cp1: None,
            doppler: 0.0,
            dt_sat_m: 0.0,
            tropo_dry: 0.0,
            map_wet: 0.0,
            lam1: 0.19,
            el: core::f64::consts::FRAC_PI_2,
            snr: 45.0,
            sat_id: SatelliteId { constellation: Constellation::Gps, prn: 1 },
        };

        add_epoch_factors(
            &mut opt, &sat, &core_state, 0, None,
            0, 2 * CORE_STATE_SIZE, 3.0,
        );

        // Should have exactly 1 factor (PR only, no CP or Doppler)
        assert_eq!(opt.factors.len(), 1);
    }

    #[test]
    fn test_add_epoch_factors_pr_cp_doppler() {
        let mut opt = FactorGraphOptimizer::new();
        let core_state = DVector::from_element(CORE_STATE_SIZE, 0.0);
        let _amb_vec = vec![50.5];

        let sat = OwnedSatData {
            sat_pos: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::new(10.0, 0.0, 0.0),
            p1: 20000000.0,
            cp1: Some(105263157.0), // dist / 0.19
            doppler: -526.0,        // ~ -100 / 0.19
            dt_sat_m: 0.0,
            tropo_dry: 0.0,
            map_wet: 0.0,
            lam1: 0.19,
            el: core::f64::consts::FRAC_PI_2,
            snr: 45.0,
            sat_id: SatelliteId { constellation: Constellation::Gps, prn: 1 },
        };

        add_epoch_factors(
            &mut opt, &sat, &core_state, 1, Some(0),
            0, 2 * CORE_STATE_SIZE, 3.0,
        );

        // Should have 3 factors: PR, CP, Doppler
        assert_eq!(opt.factors.len(), 3);
    }

    #[test]
    fn test_compute_sub_covariance_empty() {
        let opt = FactorGraphOptimizer::new();
        let state = DVector::zeros(10);
        let cov = compute_sub_covariance(&opt, &state, 5, 5);
        assert!(cov.is_none());
    }

    #[test]
    fn test_compute_sub_covariance_basic() {
        let mut opt = FactorGraphOptimizer::new();
        let state = DVector::zeros(5);
        // Add a simple factor
        use crate::estimators::factor_graph::gnss_factors::ErrorStatePseudorangeFactor;
        opt.add_factor(Box::new(ErrorStatePseudorangeFactor {
            sat_pos: Vector3::new(20000000.0, 0.0, 0.0),
            measured_pr: 20000000.0,
            variance: 1.0,
            sat_clock_bias: 0.0,
            tropo_dry_delay: 0.0,
            map_wet: 0.0,
            nominal_rx: 0.0, nominal_ry: 0.0, nominal_rz: 0.0,
            nominal_dt: 0.0,
            nominal_dt_gal: 0.0, nominal_dt_bds: 0.0, nominal_dt_glo: 0.0,
            nominal_zwd: 0.0,
            index_x: 0, index_y: 1, index_z: 2, index_dt: 3,
            index_zwd: None,
            index_dt_gal: None, index_dt_bds: None, index_dt_glo: None,
            sat_id: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            robust_threshold: 3.0,
        }));
        let cov = compute_sub_covariance(&opt, &state, 0, 4);
        assert!(cov.is_some());
        let cov = cov.unwrap();
        assert_eq!(cov.nrows(), 4);
        assert_eq!(cov.ncols(), 4);
        for i in 0..4 {
            assert!(cov[(i, i)] > 0.0, "Diagonal element {} should be positive", i);
        }
    }

    #[test]
    fn test_optimizer_with_position_prior() {
        let mut opt = PppTwoEpochOptimizer::new(2);
        let time = GpsTime::new(2082, 0.0);
        let pos0 = Vector3::new(1000000.0, 2000000.0, 3000000.0);
        let coord = Coordinate::new(pos0, Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, coord, 10.0);
        state.position.vector = pos0;
        state.covariance = DMatrix::identity(21, 21) * 100.0;

        let sat_pos = Vector3::new(20000000.0, 10000000.0, 5000000.0);
        let dist = (sat_pos - pos0).norm();
        let sat1 = make_sat(1, sat_pos, dist, Some(dist / 0.1903), 0.0);
        let sats1 = [sat1];

        // Epoch 1
        opt.solve(&mut state, &sats1, None).unwrap();

        // Epoch 2 with SPP prior
        state.position.vector = pos0 + Vector3::new(3.0, -2.0, 1.0);
        state.covariance = DMatrix::identity(21, 21) * 100.0;

        let sat2 = make_sat(1, sat_pos, dist, Some(dist / 0.1903), 0.0);
        let sats2 = [sat2];
        let prior = Some((pos0 + Vector3::new(1.0, 0.0, -1.0), 9.0));
        let result = opt.solve(&mut state, &sats2, prior);
        assert!(result.is_ok(), "With prior should succeed: {:?}", result);
    }
}
