// Shared helpers for PPP IEKF solvers (ppp_iekf.rs and ppp_ins_iekf.rs).
use crate::engine::processed_sat::ProcessedSat;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::math::inversion::invert_matrix_robust;
use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};

pub(crate) const NOMINAL_SNR_DBHZ: f64 = 45.0;
pub(crate) const SNR_SCALE_DIVISOR: f64 = 10.0;

pub(crate) fn snr_scale(snr: i32) -> f64 {
    (10.0f64).powf((NOMINAL_SNR_DBHZ - snr as f64) / SNR_SCALE_DIVISOR)
}

pub(crate) fn invert_matrix(mat: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    if mat.iter().any(|x| x.is_nan()) { return None; }
    Some(invert_matrix_robust(mat))
}

pub(crate) fn find_ambiguity_index(state: &RtkState, sat: gneiss_core::sat::SatelliteId) -> Option<usize> {
    state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 0)
}

pub(crate) fn find_amb_idx(state: &RtkState, sat: gneiss_core::sat::SatelliteId, f: u8) -> Option<usize> {
    state.ambiguity_keys.iter().position(|&x| x == (sat, f))
}

pub(crate) fn extract_state_vector(state: &RtkState) -> DVector<f64> {
    let mut x = DVector::zeros(CORE_STATE_SIZE + state.ambiguities.len());
    x[0] = state.position.vector.x;
    x[1] = state.position.vector.y;
    x[2] = state.position.vector.z;
    x[3] = state.velocity.x;
    x[4] = state.velocity.y;
    x[5] = state.velocity.z;
    let r = state.attitude.scaled_axis();
    x[6] = r.x; x[7] = r.y; x[8] = r.z;
    x[9] = state.accel_bias.x; x[10] = state.accel_bias.y; x[11] = state.accel_bias.z;
    x[12] = state.gyro_bias.x; x[13] = state.gyro_bias.y; x[14] = state.gyro_bias.z;
    x[15] = state.rcv_clk_bias;
    if CORE_STATE_SIZE > 16 {
        x[16] = state.isb_glo; x[17] = state.isb_gal; x[18] = state.isb_bds;
        x[19] = state.rcv_clk_drift; x[20] = state.zwd;
    }
    for (i, &amb) in state.ambiguities.iter().enumerate() { x[CORE_STATE_SIZE + i] = amb; }
    x
}

pub(crate) fn apply_state_vector(state: &mut RtkState, x: &DVector<f64>, p: DMatrix<f64>) {
    state.position.vector.x = x[0]; state.position.vector.y = x[1]; state.position.vector.z = x[2];
    state.velocity.x = x[3]; state.velocity.y = x[4]; state.velocity.z = x[5];
    let r_vec = Vector3::new(x[6], x[7], x[8]);
    state.attitude = if r_vec.norm() > 1e-12 { UnitQuaternion::from_scaled_axis(r_vec) } else { UnitQuaternion::identity() };
    state.accel_bias.x = x[9]; state.accel_bias.y = x[10]; state.accel_bias.z = x[11];
    state.gyro_bias.x = x[12]; state.gyro_bias.y = x[13]; state.gyro_bias.z = x[14];
    state.rcv_clk_bias = x[15];
    if CORE_STATE_SIZE > 16 {
        state.isb_glo = x[16]; state.isb_gal = x[17]; state.isb_bds = x[18];
        state.rcv_clk_drift = x[19]; state.zwd = x[20];
    }
    for i in 0..state.ambiguities.len() { state.ambiguities[i] = x[CORE_STATE_SIZE + i]; }
    state.covariance = p;
}

pub(crate) fn build_weight_matrix(meas: &[FgMeasurement], r_mat: &DMatrix<f64>) -> DMatrix<f64> {
    let mut w_mat = DMatrix::zeros(meas.len(), meas.len());
    for i in 0..meas.len() { w_mat[(i, i)] = 1.0 / r_mat[(i, i)]; }
    w_mat
}

pub(crate) fn assemble_matrices(
    meas: &[FgMeasurement], cols: usize,
) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
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

pub(crate) fn build_iono_constraint_row(size: usize, i1_idx: usize) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[i1_idx] = 1.0;
    h
}

pub struct FgMeasurement {
    pub res: f64,
    pub h_row: DVector<f64>,
    pub weight: f64,
    pub raw_var: f64,
    pub is_phase: bool,
    pub sat: Option<gneiss_core::sat::SatelliteId>,
}
