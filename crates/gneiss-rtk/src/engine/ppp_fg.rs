use nalgebra::{DMatrix, DVector, Vector3};
use crate::engine::processed_sat::ProcessedSat;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::engine::EngineError;

const SPEED_OF_LIGHT: f64 = 299792458.0;
const NOMINAL_SNR_DBHZ: f64 = 45.0;
const SNR_SCALE_DIVISOR: f64 = 10.0;
const PSEUDORANGE_VARIANCE_BASE: f64 = 9.0;
const CARRIER_PHASE_VARIANCE_BASE: f64 = 0.0001;
const MATRIX_REGULARIZATION: f64 = 1e-6;

/// Configuration for the Iterated EKF (Factor Graph) solver.
pub struct PppFactorGraph {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub huber_k: f64,
}

impl Default for PppFactorGraph {
    fn default() -> Self {
        Self { max_iterations: 15, convergence_threshold: 1e-3, huber_k: 3.0 }
    }
}

pub struct FgMeasurement {
    pub res: f64,
    pub h_row: DVector<f64>,
    pub weight: f64,
    pub is_phase: bool,
}

impl PppFactorGraph {
    pub fn new() -> Self { Self::default() }

    pub fn solve(&self, state: &mut RtkState, sats: &[ProcessedSat]) -> Result<(), EngineError> {
        let x_pred = extract_state_vector(state);
        let p_pred = state.covariance.clone();
        let mut x_i = x_pred.clone();

        for _ in 0..self.max_iterations {
            let meas = self.build_measurements(state, &x_i, sats);
            if meas.is_empty() { return Err(EngineError::InsufficientSatellites); }

            let (h_mat, res_vec, r_mat) = assemble_matrices(&meas, x_i.len());
            let s = &h_mat * &p_pred * h_mat.transpose() + &r_mat;
            let s_inv = invert_matrix(&s).ok_or(EngineError::StateDisappeared)?;

            let k = &p_pred * h_mat.transpose() * s_inv;
            let innov = &res_vec - &h_mat * (&x_pred - &x_i);
            
            let dx = &k * innov;
            x_i = &x_pred + &dx;

            if dx.norm() < self.convergence_threshold {
                let id = DMatrix::identity(x_i.len(), x_i.len());
                let i_kh = id - &k * h_mat;
                state.covariance = &i_kh * &p_pred * i_kh.transpose() + &k * r_mat * k.transpose();
                break;
            }
        }
        
        apply_state_vector(state, &x_i);
        Ok(())
    }

    fn build_measurements(&self, state: &RtkState, x_i: &DVector<f64>, sats: &[ProcessedSat]) -> Vec<FgMeasurement> {
        let mut meas = Vec::new();
        let rcv_pos = Vector3::new(x_i[0], x_i[1], x_i[2]);
        let rcv_clk = x_i[15];
        let ztd = if x_i.len() > 17 { x_i[17] } else { 0.0 };

        for sat in sats {
            let dist = (sat.sat_pos_rot - rcv_pos).norm();
            let los = (rcv_pos - sat.sat_pos_rot) / dist;
            let expected_p = dist + rcv_clk - (sat.sat_clock_drift * SPEED_OF_LIGHT) + sat.tropo_dry + ztd * sat.map_wet;

            if let Some(pr1) = sat.sat_obs.get_observable(1) {
                let res = pr1 - expected_p;
                let var = PSEUDORANGE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
                let w = apply_huber(res, var, self.huber_k);
                meas.push(FgMeasurement { res, h_row: build_h_row(&los, sat.map_wet, None, x_i.len()), weight: w, is_phase: false });
            }

            if let (Some(cp1), Some(cp2)) = (sat.cp1, sat.cp2) {
                if let Some(amb_idx) = find_ambiguity_index(state, sat.sat_obs.sat) {
                    let l_if = crate::combinations::iono_free(cp1 * sat.lam1, cp2 * sat.lam2, sat.f1, sat.f2);
                    let expected_l = expected_p + x_i[CORE_STATE_SIZE + amb_idx];
                    let res = l_if - expected_l;
                    
                    let var = CARRIER_PHASE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
                    let w = apply_huber(res, var, self.huber_k);
                    meas.push(FgMeasurement { res, h_row: build_h_row(&los, sat.map_wet, Some(CORE_STATE_SIZE + amb_idx), x_i.len()), weight: w, is_phase: true });
                }
            }
        }
        meas
    }
}

fn apply_huber(res: f64, var: f64, k: f64) -> f64 {
    let norm_r = res.abs() / var.sqrt();
    if norm_r > k { var * (norm_r / k) } else { var }
}

fn snr_scale(snr: i32) -> f64 {
    (10.0f64).powf((NOMINAL_SNR_DBHZ - snr as f64) / SNR_SCALE_DIVISOR)
}

fn invert_matrix(mat: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    mat.clone().cholesky().map(|c| c.inverse())
        .or_else(|| (mat.clone() + DMatrix::identity(mat.nrows(), mat.ncols()) * MATRIX_REGULARIZATION).try_inverse())
}

fn find_ambiguity_index(state: &RtkState, sat: gneiss_core::sat::SatelliteId) -> Option<usize> {
    state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 0)
}

fn build_h_row(los: &Vector3<f64>, map_wet: f64, amb_idx: Option<usize>, size: usize) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[0] = -los.x; h[1] = -los.y; h[2] = -los.z;
    h[15] = 1.0;
    if size > 17 { h[17] = map_wet; }
    if let Some(idx) = amb_idx { h[idx] = 1.0; }
    h
}

fn extract_state_vector(state: &RtkState) -> DVector<f64> {
    let mut x = DVector::zeros(CORE_STATE_SIZE + state.ambiguities.len());
    x[0] = state.position.vector.x; x[1] = state.position.vector.y; x[2] = state.position.vector.z;
    x[3] = state.velocity.x; x[4] = state.velocity.y; x[5] = state.velocity.z;
    let r = state.attitude.scaled_axis(); x[6] = r.x; x[7] = r.y; x[8] = r.z;
    x[9] = state.accel_bias.x; x[10] = state.accel_bias.y; x[11] = state.accel_bias.z;
    x[12] = state.gyro_bias.x; x[13] = state.gyro_bias.y; x[14] = state.gyro_bias.z;
    x[15] = state.rcv_clk_bias;
    if CORE_STATE_SIZE > 16 { x[16] = state.rcv_clk_drift; x[17] = state.zwd; }
    for (i, &amb) in state.ambiguities.iter().enumerate() { x[CORE_STATE_SIZE + i] = amb; }
    x
}

fn apply_state_vector(state: &mut RtkState, x: &DVector<f64>) {
    state.position.vector = Vector3::new(x[0], x[1], x[2]);
    state.velocity = Vector3::new(x[3], x[4], x[5]);
    state.attitude = nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(x[6], x[7], x[8]));
    state.accel_bias = Vector3::new(x[9], x[10], x[11]);
    state.gyro_bias = Vector3::new(x[12], x[13], x[14]);
    state.rcv_clk_bias = x[15];
    if CORE_STATE_SIZE > 16 { state.rcv_clk_drift = x[16]; state.zwd = x[17]; }
    for i in 0..state.ambiguities.len() { state.ambiguities[i] = x[CORE_STATE_SIZE + i]; }
}

fn assemble_matrices(meas: &[FgMeasurement], cols: usize) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
    let mut h = DMatrix::zeros(meas.len(), cols);
    let mut z = DVector::zeros(meas.len());
    let mut r = DMatrix::zeros(meas.len(), meas.len());
    for (i, m) in meas.iter().enumerate() {
        for j in 0..cols { h[(i, j)] = m.h_row[j]; }
        z[i] = m.res;
        r[(i, i)] = m.weight;
    }
    (h, z, r)
}
