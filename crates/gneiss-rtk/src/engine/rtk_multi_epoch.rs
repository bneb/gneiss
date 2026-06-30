//! Two-epoch RTK factor graph for AR validation.
//!
//! Jointly estimates position and DD ambiguities across epochs k-1 and k
//! with shared ambiguity states. The independent geometric constraint from
//! different satellite geometries breaks the code-multipath circularity that
//! biases single-epoch AR validation (MW, NL, PR residuals).
//!
//! The factor graph uses PRE-COMPUTED DD innovations and LOS Jacobians from
//! the EKF's H matrix, making all measurement factors linear. The only
//! non-linearity comes from the dynamics constraint coupling the two epochs.

use nalgebra::{DMatrix, DVector, Vector3};

use crate::engine::measurement::builder::EkfMeasurementMatrices;
use crate::estimators::factor_graph::gnss_factors::ErrorStateDdPseudorangeFactor;
use crate::estimators::factor_graph::gnss_factors::ErrorStateDdCarrierPhaseFactor;
use crate::estimators::factor_graph::{Factor, FactorGraphOptimizer};
use crate::filter::{CORE_STATE_SIZE, RtkState};
use gneiss_core::sat::SatelliteId;

// ---------------------------------------------------------------------------
// Factor graph state layout
// ---------------------------------------------------------------------------
// [pos_{k-1}(3), vel_{k-1}(3), clk_{k-1}(1), pos_k(3), vel_k(3), clk_k(1), ambiguities(N)]
//  0,1,2         3,4,5          6             7,8,9    10,11,12   13          14..14+N-1

const OFF_PREV_POS: usize = 0;
const OFF_PREV_VEL: usize = 3;
const OFF_PREV_CLK: usize = 6;
const OFF_CURR_POS: usize = 7;
const OFF_CURR_VEL: usize = 10;
const OFF_CURR_CLK: usize = 13;
const OFF_AMB: usize = 14;

#[allow(dead_code)]
const DIM_PER_EPOCH: usize = 7;
const DIM_BASE: usize = 14; // 2 * DIM_PER_EPOCH

/// Factor graph validation result.
pub struct FactorGraphResult {
    /// Optimized position at epoch k (ECEF, meters).
    pub pos_k: Vector3<f64>,
    /// Whether the optimization converged.
    pub converged: bool,
    /// Final residual norm (diagnostic).
    pub final_error: f64,
    /// Ambiguity corrections (delta from float values) in cycles.
    /// Index maps to state.ambiguity_keys.
    pub amb_corrections: Vec<f64>,
    /// Ambiguity keys matching amb_corrections indices.
    pub amb_keys: Vec<(SatelliteId, u8)>,
}

/// Run the two-epoch factor graph to get an independent position estimate at
/// epoch k. Uses pre-computed DD innovations from both the current EKF epoch
/// and the stored previous epoch.
///
/// Returns `None` if the factor graph fails (singular system, NaN, etc.).
pub fn run_two_epoch_factor_graph(
    state: &RtkState,
    m: &EkfMeasurementMatrices,
) -> Option<FactorGraphResult> {
    let prev = state.prev_epoch_meas.as_ref()?;
    let n_amb = state.ambiguities.len();
    let total_dim = DIM_BASE + n_amb;

    if n_amb == 0 {
        tracing::debug!("RTK FG: no ambiguities, skipping");
        return None;
    }

    let dt = (state.time.tow - prev.time).abs().max(0.1);

    let mut optimizer = FactorGraphOptimizer::new();

    // --- Measurement factors for epoch k-1 (previous) ---
    let prev_amb_map = build_ambiguity_map(&prev.ambiguity_keys, &state.ambiguity_keys);
    for meas in &prev.measurements {
        add_stored_dd_factor(
            &mut optimizer, meas, OFF_PREV_POS, OFF_AMB,
            total_dim, &prev_amb_map,
        );
    }

    // --- Measurement factors for epoch k (current) ---
    for i in 0..m.z.nrows() {
        let row_data = extract_h_row(m, i, &state.ambiguity_keys);
        add_dd_factor(
            &mut optimizer, &row_data, OFF_CURR_POS, OFF_AMB,
            total_dim, &build_identity_amb_map(n_amb),
        );
    }

    if optimizer.factors.is_empty() {
        tracing::debug!("RTK FG: no factors built, skipping");
        return None;
    }

    // --- Dynamics constraint: pos_k ≈ pos_{k-1} + vel_{k-1}·dt ---
    add_dynamics_factor(
        &mut optimizer,
        state.position.vector,
        state.velocity,
        state.rcv_clk_bias,
        prev.pos,
        prev.vel,
        prev.clk,
        dt,
        total_dim,
    );

    // --- Prior on pos_{k-1}: trust the EKF position within σ=1.0m floor ---
    let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
    let prior_sigma = (pos_var / 3.0).sqrt().max(1.0);
    let prior_info = 1.0 / (prior_sigma * prior_sigma);

    struct SimplePositionPrior {
        nominal: Vector3<f64>,
        info: f64,
        total_dim: usize,
    }

    impl Factor for SimplePositionPrior {
        fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
            DVector::from_vec(vec![
                self.nominal.x + delta[OFF_PREV_POS] - self.nominal.x,
                self.nominal.y + delta[OFF_PREV_POS + 1] - self.nominal.y,
                self.nominal.z + delta[OFF_PREV_POS + 2] - self.nominal.z,
            ])
        }
        fn jacobian(&self, _delta: &DVector<f64>) -> DMatrix<f64> {
            let mut jac = DMatrix::zeros(3, self.total_dim);
            jac[(0, OFF_PREV_POS)] = 1.0;
            jac[(1, OFF_PREV_POS + 1)] = 1.0;
            jac[(2, OFF_PREV_POS + 2)] = 1.0;
            jac
        }
        fn information(&self) -> DMatrix<f64> {
            DMatrix::from_diagonal(&DVector::from_element(3, self.info))
        }
    }

    optimizer.add_factor(Box::new(SimplePositionPrior {
        nominal: prev.pos,
        info: prior_info,
        total_dim,
    }));

    // --- Solve ---
    let initial_delta = DVector::zeros(total_dim);
    let (delta_opt, _cov) = optimizer.optimize(&initial_delta, 10, 1e-4);

    // Check for NaN/Inf
    if delta_opt.iter().any(|x| x.is_nan() || x.is_infinite()) {
        tracing::warn!("RTK FG: solution contains NaN/Inf, rejecting");
        return None;
    }

    let pos_k = Vector3::new(
        state.position.vector.x + delta_opt[OFF_CURR_POS],
        state.position.vector.y + delta_opt[OFF_CURR_POS + 1],
        state.position.vector.z + delta_opt[OFF_CURR_POS + 2],
    );

    let final_error: f64 = optimizer.factors.iter()
        .map(|f| f.residual(&delta_opt).norm_squared())
        .sum();

    let pos_delta = Vector3::new(
        delta_opt[OFF_CURR_POS],
        delta_opt[OFF_CURR_POS + 1],
        delta_opt[OFF_CURR_POS + 2],
    ).norm();

    tracing::info!(
        "RTK FG: {} factors, pos_delta={:.3}m, final_err={:.2}, pos_k=({:.1},{:.1},{:.1})",
        optimizer.factors.len(), pos_delta, final_error,
        pos_k.x, pos_k.y, pos_k.z
    );

    // Extract ambiguity corrections (delta from float values).
    // The factor graph ambiguity states are in DD cycles; the EKF stores
    // UDUC ambiguities in meters. We extract the DD corrections here
    // and convert to cycles in the caller.
    let amb_corrections: Vec<f64> = (0..n_amb)
        .map(|i| delta_opt[OFF_AMB + i])
        .collect();
    let amb_keys = state.ambiguity_keys.clone();

    Some(FactorGraphResult {
        pos_k,
        converged: pos_delta < 10.0,
        final_error,
        amb_corrections,
        amb_keys,
    })
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Extracted data for one DD measurement row from the EKF H matrix.
struct HRowData {
    sat_id: SatelliteId,
    ref_sat_id: SatelliteId,
    z: f64,
    h_pos: [f64; 3],
    variance: f64,
    is_pr: bool,
    #[allow(dead_code)]
    freq_band: u8,
    amb_idx: Option<usize>,
    ref_amb_idx: Option<usize>,
    iono_pair: Option<(usize, usize, f64)>,
}

/// Extract per-measurement data from one row of the EKF H matrix.
///
/// The EKF's H matrix has:
/// - Position columns (0,1,2) = h_r = e_ref - e_sat
/// - Ambiguity columns (CORE_STATE_SIZE + amb_idx): +1.0 for sat, -1.0 for ref
/// - Ionosphere columns: +scale for sat, -scale for ref
fn extract_h_row(
    m: &EkfMeasurementMatrices,
    row: usize,
    ambiguity_keys: &[(SatelliteId, u8)],
) -> HRowData {
    let (sat_id, type_code, _r_ref) = m.mt[row];
    let z = m.z[row];
    let variance = m.r[(row, row)].max(1e-4);

    let h_pos = [m.h[(row, 0)], m.h[(row, 1)], m.h[(row, 2)]];

    let is_pr = type_code == 0;
    let freq_band = match type_code {
        0 => 1, // PR on L1
        1 => 1, // CP on L1
        2 => 2, // CP on L2
        _ => 1,
    };

    // Find ambiguity indices from H matrix columns >= CORE_STATE_SIZE
    let mut amb_idx: Option<usize> = None;
    let mut ref_amb_idx: Option<usize> = None;
    let mut iono_sat_idx: Option<usize> = None;
    let mut iono_ref_idx: Option<usize> = None;
    let mut iono_scale: Option<f64> = None;

    if !is_pr {
        for col in CORE_STATE_SIZE..m.h.ncols() {
            let val = m.h[(row, col)];
            if val.abs() < 1e-9 {
                continue;
            }
            let key_idx = col - CORE_STATE_SIZE;
            if key_idx >= ambiguity_keys.len() {
                continue;
            }
            let (_sat, freq) = ambiguity_keys[key_idx];
            if freq == 3 {
                // Ionosphere state
                if val > 0.0 {
                    iono_sat_idx = Some(key_idx);
                    iono_scale = Some(val);
                } else {
                    iono_ref_idx = Some(key_idx);
                }
            } else if val > 0.0 {
                amb_idx = Some(key_idx);
            } else {
                ref_amb_idx = Some(key_idx);
            }
        }
    }

    // Determine reference satellite from the ambiguity pair
    let ref_sat_id = ref_amb_idx
        .and_then(|ri| ambiguity_keys.get(ri).map(|k| k.0))
        .unwrap_or(sat_id);

    let iono_pair = match (iono_sat_idx, iono_ref_idx, iono_scale) {
        (Some(si), Some(ri), Some(sc)) => Some((si, ri, sc)),
        _ => None,
    };

    HRowData {
        sat_id,
        ref_sat_id,
        z,
        h_pos,
        variance,
        is_pr,
        freq_band,
        amb_idx,
        ref_amb_idx,
        iono_pair,
    }
}

/// Map ambiguity indices from a previous epoch's `ambiguity_keys` to
/// the current combined `ambiguity_keys`. Returns a vector where
/// `map[prev_idx] = Some(current_idx)` or `None` if the key was dropped.
fn build_ambiguity_map(
    prev_keys: &[(SatelliteId, u8)],
    curr_keys: &[(SatelliteId, u8)],
) -> Vec<Option<usize>> {
    prev_keys
        .iter()
        .map(|key| curr_keys.iter().position(|k| k == key))
        .collect()
}

/// Identity ambiguity map (current epoch measurements map directly).
fn build_identity_amb_map(n_amb: usize) -> Vec<Option<usize>> {
    (0..n_amb).map(Some).collect()
}

/// Add a DD measurement factor to the optimizer.
fn add_dd_factor(
    opt: &mut FactorGraphOptimizer,
    data: &HRowData,
    pos_offset: usize,
    amb_offset: usize,
    total_dim: usize,
    amb_map: &[Option<usize>],
) {
    if data.is_pr {
        opt.add_factor(Box::new(ErrorStateDdPseudorangeFactor {
            sat_id: data.sat_id,
            ref_sat_id: data.ref_sat_id,
            z: data.z,
            h_pos: data.h_pos,
            variance: data.variance,
            index_px: pos_offset,
            index_py: pos_offset + 1,
            index_pz: pos_offset + 2,
            total_dim,
        }));
    } else {
        // Map epoch-specific ambiguity indices to factor graph indices
        let fg_sat_amb = data.amb_idx.and_then(|ai| amb_map.get(ai).copied().flatten());
        let fg_ref_amb = data.ref_amb_idx.and_then(|ri| amb_map.get(ri).copied().flatten());

        let (sat_amb, ref_amb) = match (fg_sat_amb, fg_ref_amb) {
            (Some(sa), Some(ra)) => (sa, ra),
            _ => {
                // Ambiguity not present in combined state — skip this CP
                // measurement (no ambiguity to constrain the phase).
                return;
            }
        };

        let iono_pair = data.iono_pair.and_then(|(si, ri, scale)| {
            let fg_si = amb_map.get(si).copied().flatten()?;
            let fg_ri = amb_map.get(ri).copied().flatten()?;
            Some((fg_si, fg_ri, scale))
        });

        opt.add_factor(Box::new(ErrorStateDdCarrierPhaseFactor {
            sat_id: data.sat_id,
            ref_sat_id: data.ref_sat_id,
            z: data.z,
            h_pos: data.h_pos,
            variance: data.variance,
            index_px: pos_offset,
            index_py: pos_offset + 1,
            index_pz: pos_offset + 2,
            index_amb_sat: amb_offset + sat_amb,
            index_amb_ref: amb_offset + ref_amb,
            iono_pair: iono_pair.map(|(si, ri, s)| (amb_offset + si, amb_offset + ri, s)),
            total_dim,
        }));
    }
}

/// Add a stored DD measurement factor (from previous epoch) to the optimizer.
fn add_stored_dd_factor(
    opt: &mut FactorGraphOptimizer,
    data: &crate::filter::StoredDdMeasurement,
    pos_offset: usize,
    amb_offset: usize,
    total_dim: usize,
    amb_map: &[Option<usize>],
) {
    if data.is_pr {
        opt.add_factor(Box::new(ErrorStateDdPseudorangeFactor {
            sat_id: data.sat_id,
            ref_sat_id: data.ref_sat_id,
            z: data.z,
            h_pos: data.h_pos,
            variance: data.variance,
            index_px: pos_offset,
            index_py: pos_offset + 1,
            index_pz: pos_offset + 2,
            total_dim,
        }));
    } else {
        let fg_sat_amb = data.amb_idx.and_then(|ai| amb_map.get(ai).copied().flatten());
        let fg_ref_amb = data.ref_amb_idx.and_then(|ri| amb_map.get(ri).copied().flatten());

        let (sat_amb, ref_amb) = match (fg_sat_amb, fg_ref_amb) {
            (Some(sa), Some(ra)) => (sa, ra),
            _ => return,
        };

        let iono_pair = data.iono_pair.and_then(|(si, ri, scale)| {
            let fg_si = amb_map.get(si).copied().flatten()?;
            let fg_ri = amb_map.get(ri).copied().flatten()?;
            Some((fg_si, fg_ri, scale))
        });

        opt.add_factor(Box::new(ErrorStateDdCarrierPhaseFactor {
            sat_id: data.sat_id,
            ref_sat_id: data.ref_sat_id,
            z: data.z,
            h_pos: data.h_pos,
            variance: data.variance,
            index_px: pos_offset,
            index_py: pos_offset + 1,
            index_pz: pos_offset + 2,
            index_amb_sat: amb_offset + sat_amb,
            index_amb_ref: amb_offset + ref_amb,
            iono_pair: iono_pair.map(|(si, ri, s)| (amb_offset + si, amb_offset + ri, s)),
            total_dim,
        }));
    }
}

/// Add a soft dynamics constraint between epoch k-1 and epoch k.
///
/// Models: pos_k ≈ pos_{k-1} + vel_{k-1}·dt
///         vel_k ≈ vel_{k-1}
///         clk_k ≈ clk_{k-1}
///
/// Process noise is scaled by `dt` for position and velocity.
fn add_dynamics_factor(
    opt: &mut FactorGraphOptimizer,
    pos_curr: Vector3<f64>,
    vel_curr: Vector3<f64>,
    clk_curr: f64,
    pos_prev: Vector3<f64>,
    vel_prev: Vector3<f64>,
    clk_prev: f64,
    dt: f64,
    total_dim: usize,
) {
    let dt = dt.max(0.1);

    // Process noise intensities (per √s for pos/vel, per epoch for clk).
    // sigma_pos controls how much position can deviate from the constant-velocity
    // prediction. 3.0 m/√s allows ~3m position deviation after 1s (highway
    // acceleration) while still providing weak temporal regularization.
    // sigma_vel controls velocity random walk. 1.0 m/s/√s allows mild
    // acceleration changes between epochs.
    let sigma_pos: f64 = 3.0; // m/√s — accommodates vehicle acceleration
    let sigma_vel: f64 = 1.0; // m/s/√s
    let sigma_clk: f64 = 10.0; // m — loose, DD cancels clock anyway

    // Build Q_inv (7x7 diagonal)
    let q_inv_diag = vec![
        1.0 / (sigma_pos * sigma_pos * dt), // px
        1.0 / (sigma_pos * sigma_pos * dt), // py
        1.0 / (sigma_pos * sigma_pos * dt), // pz
        1.0 / (sigma_vel * sigma_vel * dt), // vx
        1.0 / (sigma_vel * sigma_vel * dt), // vy
        1.0 / (sigma_vel * sigma_vel * dt), // vz
        1.0 / (sigma_clk * sigma_clk), // clk (white noise, no dt scaling)
    ];

    // Dynamics residual: (x_k - x_{k-1} - vel_{k-1}·dt) for position,
    // (v_k - v_{k-1}) for velocity, (clk_k - clk_{k-1}) for clock.
    // In error-state form:
    //   res = (nominal_curr + delta_curr) - A·(nominal_prev + delta_prev)
    // where A encodes the velocity integration for position.

    let nominal_prev = DVector::from_vec(vec![
        pos_prev.x, pos_prev.y, pos_prev.z,
        vel_prev.x, vel_prev.y, vel_prev.z,
        clk_prev,
    ]);

    let nominal_curr = DVector::from_vec(vec![
        pos_curr.x, pos_curr.y, pos_curr.z,
        vel_curr.x, vel_curr.y, vel_curr.z,
        clk_curr,
    ]);

    let q_inv = DMatrix::from_diagonal(&DVector::from_vec(q_inv_diag));

    opt.add_factor(Box::new(TwoEpochDynamicsFactor {
        nominal_prev,
        nominal_curr,
        dt,
        q_inv,
        total_dim,
    }));
}

/// Dynamics factor: penalises deviation from the linear motion model.
struct TwoEpochDynamicsFactor {
    nominal_prev: DVector<f64>, // 7 elements
    nominal_curr: DVector<f64>, // 7 elements
    dt: f64,
    q_inv: DMatrix<f64>, // 7×7 inverse process noise
    total_dim: usize,
}

impl Factor for TwoEpochDynamicsFactor {
    fn residual(&self, delta: &DVector<f64>) -> DVector<f64> {
        let d_prev = delta.fixed_rows::<7>(OFF_PREV_POS).clone_owned();
        let d_curr = delta.fixed_rows::<7>(OFF_CURR_POS).clone_owned();
        let full_prev = &self.nominal_prev + &d_prev;
        let full_curr = &self.nominal_curr + &d_curr;

        let mut res = DVector::zeros(7);
        // Position: pos_k - (pos_{k-1} + vel_{k-1}·dt)
        res[0] = full_curr[0] - (full_prev[0] + full_prev[3] * self.dt);
        res[1] = full_curr[1] - (full_prev[1] + full_prev[4] * self.dt);
        res[2] = full_curr[2] - (full_prev[2] + full_prev[5] * self.dt);
        // Velocity: v_k - v_{k-1}
        res[3] = full_curr[3] - full_prev[3];
        res[4] = full_curr[4] - full_prev[4];
        res[5] = full_curr[5] - full_prev[5];
        // Clock: clk_k - clk_{k-1}
        res[6] = full_curr[6] - full_prev[6];
        res
    }

    fn jacobian(&self, _delta: &DVector<f64>) -> DMatrix<f64> {
        let mut jac = DMatrix::zeros(7, self.total_dim);
        // d(res_pos)/d(dpos_{k-1}) = -1
        for i in 0..3 {
            jac[(i, OFF_PREV_POS + i)] = -1.0;
        }
        // d(res_pos)/d(dvel_{k-1}) = -dt
        for i in 0..3 {
            jac[(i, OFF_PREV_VEL + i)] = -self.dt;
        }
        // d(res_pos)/d(dpos_k) = 1
        for i in 0..3 {
            jac[(i, OFF_CURR_POS + i)] = 1.0;
        }
        // d(res_vel)/d(dvel_{k-1}) = -1
        for i in 0..3 {
            jac[(3 + i, OFF_PREV_VEL + i)] = -1.0;
        }
        // d(res_vel)/d(dvel_k) = 1
        for i in 0..3 {
            jac[(3 + i, OFF_CURR_VEL + i)] = 1.0;
        }
        // d(res_clk)/d(dclk_{k-1}) = -1
        jac[(6, OFF_PREV_CLK)] = -1.0;
        // d(res_clk)/d(dclk_k) = 1
        jac[(6, OFF_CURR_CLK)] = 1.0;
        jac
    }

    fn information(&self) -> DMatrix<f64> {
        // Return 7×7 information matrix matching the 7-dimensional residual.
        // The residual function returns a 7-element vector (pos, vel, clk),
        // so the information matrix must be 7×7, NOT total_dim×total_dim.
        self.q_inv.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::measurement::builder::EkfMeasurementMatrices;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use nalgebra::Vector3;

    fn make_sat(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn }
    }

    #[test]
    fn test_extract_h_row_pr() {
        let sat = make_sat(1);
        let h = DMatrix::from_row_slice(1, 22, &[
            0.6, -0.3, 0.8, // position
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // vel+att
            0.0, 0.0, 0.0, // accel bias
            0.0, 0.0, 0.0, // gyro bias
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // clock, ISBs, drift, ZWD
            0.0, // no ambiguity
        ]);
        let m = EkfMeasurementMatrices {
            z: DVector::from_vec(vec![2.5]),
            h,
            r: DMatrix::from_diagonal(&DVector::from_vec(vec![9.0])),
            mt: vec![(sat, 0, 0.0)],
        };
        let keys: Vec<(SatelliteId, u8)> = vec![];
        let row = extract_h_row(&m, 0, &keys);
        assert!(row.is_pr);
        assert_eq!(row.freq_band, 1);
        assert!((row.h_pos[0] - 0.6).abs() < 1e-9);
        assert!((row.h_pos[1] - (-0.3)).abs() < 1e-9);
        assert!((row.h_pos[2] - 0.8).abs() < 1e-9);
        assert!((row.z - 2.5).abs() < 1e-9);
        assert!(row.amb_idx.is_none());
    }

    #[test]
    fn test_extract_h_row_cp_with_amb() {
        let sat = make_sat(2);
        let ref_sat = make_sat(1);
        let keys = vec![(ref_sat, 1), (sat, 1)];

        let mut h = DMatrix::zeros(1, CORE_STATE_SIZE + 2);
        h[(0, 0)] = 0.6;
        h[(0, 1)] = -0.3;
        h[(0, 2)] = 0.8;
        h[(0, CORE_STATE_SIZE + 1)] = 1.0; // sat amb
        h[(0, CORE_STATE_SIZE)] = -1.0;     // ref amb

        let m = EkfMeasurementMatrices {
            z: DVector::from_vec(vec![-0.05]),
            h,
            r: DMatrix::from_diagonal(&DVector::from_vec(vec![0.01])),
            mt: vec![(sat, 1, 0.0)],
        };
        let row = extract_h_row(&m, 0, &keys);
        assert!(!row.is_pr);
        assert_eq!(row.freq_band, 1);
        assert_eq!(row.amb_idx, Some(1)); // sat is at index 1
        assert_eq!(row.ref_amb_idx, Some(0)); // ref is at index 0
    }

    #[test]
    fn test_build_ambiguity_map() {
        let prev = vec![(make_sat(1), 1), (make_sat(2), 1), (make_sat(3), 1)];
        let curr = vec![(make_sat(2), 1), (make_sat(1), 1), (make_sat(4), 1)];
        let amb_map = build_ambiguity_map(&prev, &curr);
        assert_eq!(amb_map[0], Some(1)); // G01 moves from idx 0 to 1
        assert_eq!(amb_map[1], Some(0)); // G02 moves from idx 1 to 0
        assert_eq!(amb_map[2], None);    // G03 dropped
    }

    #[test]
    fn test_dynamics_factor_residual_zero() {
        let pos = Vector3::new(1e6, 2e6, 3e6);
        let vel = Vector3::zeros();
        let dt = 1.0;

        // When both epochs have the same state, residual should be ~0
        let factor = TwoEpochDynamicsFactor {
            nominal_prev: DVector::from_vec(vec![pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, 0.0]),
            nominal_curr: DVector::from_vec(vec![pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, 0.0]),
            dt,
            q_inv: DMatrix::identity(7, 7),
            total_dim: 20,
        };

        let delta = DVector::zeros(20);
        let res = factor.residual(&delta);
        assert!(res.norm() < 1e-10, "Dynamics residual should be zero, norm={}", res.norm());
    }

    #[test]
    fn test_dynamics_factor_residual_with_motion() {
        let pos_prev = Vector3::new(1e6, 2e6, 3e6);
        let vel_prev = Vector3::new(10.0, 0.0, 0.0);
        let dt = 1.0;
        let pos_curr = pos_prev + vel_prev * dt;

        let factor = TwoEpochDynamicsFactor {
            nominal_prev: DVector::from_vec(vec![
                pos_prev.x, pos_prev.y, pos_prev.z,
                vel_prev.x, vel_prev.y, vel_prev.z, 0.0,
            ]),
            nominal_curr: DVector::from_vec(vec![
                pos_curr.x, pos_curr.y, pos_curr.z,
                vel_prev.x, vel_prev.y, vel_prev.z, 0.0,
            ]),
            dt,
            q_inv: DMatrix::identity(7, 7),
            total_dim: 20,
        };

        let delta = DVector::zeros(20);
        let res = factor.residual(&delta);
        assert!(res.norm() < 1e-10, "Dynamics residual with motion should be zero, norm={}", res.norm());
    }
}
