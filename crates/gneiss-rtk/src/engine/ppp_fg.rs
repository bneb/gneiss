use nalgebra::{DMatrix, DVector, Vector3};
use crate::engine::processed_sat::ProcessedSat;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::engine::EngineError;

const SPEED_OF_LIGHT: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
const NOMINAL_SNR_DBHZ: f64 = 45.0;
const SNR_SCALE_DIVISOR: f64 = 10.0;
const PSEUDORANGE_VARIANCE_BASE: f64 = 10.0;

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

        let p_inv = invert_matrix(&p_pred).ok_or(EngineError::StateDisappeared)?;

        for _iter in 0..self.max_iterations {
            let dx = self.compute_iteration_dx(state, sats, &x_i, &x_pred, &p_inv, _iter)?;
            if dx.is_none() {
                break;
            }
            let dx = dx.unwrap();
            
            x_i = &x_i + &dx;
            
            if dx.norm() < self.convergence_threshold {
                break;
            }
        }
        
        state.covariance = self.compute_final_covariance(state, sats, &x_i, &p_pred, &p_inv);

        log_ppp_convergence(state, sats, &x_i, &x_pred, &p_pred, self);
        
        apply_state_vector(state, &x_i);
        Ok(())
    }

    fn compute_iteration_dx(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, x_pred: &DVector<f64>, p_inv: &DMatrix<f64>, iter: usize) -> Result<Option<DVector<f64>>, EngineError> {
        let meas = self.build_measurements(state, sats, x_i, iter);
        if meas.is_empty() { return Err(EngineError::InsufficientSatellites); }

        let (h_mat, res_vec, r_mat) = assemble_matrices(&meas, x_i.len());
        let w_mat = build_weight_matrix(&meas, &r_mat);
        
        let h_t = h_mat.transpose();
        let htw = &h_t * &w_mat;
        let htwh = &htw * &h_mat;
        let htwr = &htw * &res_vec;
        
        let htwh_damped = &htwh + p_inv;
        let innov = &htwr + p_inv * (x_pred - x_i);
        
        match htwh_damped.cholesky() {
            Some(chol) => Ok(Some(chol.solve(&innov))),
            None => {
                tracing::warn!("Failed to solve normal equations in PPP FG!");
                Ok(None)
            }
        }
    }

    fn compute_final_covariance(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, p_pred: &DMatrix<f64>, p_inv: &DMatrix<f64>) -> DMatrix<f64> {
        let last_meas = self.build_measurements(state, sats, x_i, self.max_iterations);
        if last_meas.is_empty() { return p_pred.clone(); }
        
        let (h_mat, _, r_mat) = assemble_matrices(&last_meas, x_i.len());
        let w_mat = build_weight_matrix(&last_meas, &r_mat);
        
        let htwh = h_mat.transpose() * &w_mat * h_mat;
        let htwh_damped = htwh + p_inv;
        
        invert_matrix(&htwh_damped).unwrap_or_else(|| {
            tracing::warn!("invert_matrix(&htwh_damped) FAILED! Falling back to p_pred. htwh_damped has NaNs: {}, Infs: {}", htwh_damped.iter().any(|x| x.is_nan()), htwh_damped.iter().any(|x| x.is_infinite()));
            p_pred.clone()
        })
    }

    fn build_measurements(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, iter: usize) -> Vec<FgMeasurement> {
        let mut meas = Vec::new();
        let rcv_pos = Vector3::new(x_i[0], x_i[1], x_i[2]);
        let ztd = if x_i.len() > 17 { x_i[17] } else { 0.0 };

        for sat in sats {
            let dist = (sat.sat_pos_rot - rcv_pos).norm();
            let los = (sat.sat_pos_rot - rcv_pos) / dist;
            let expected_base = dist + x_i[15] - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet;
            let expected_pr = expected_base + sat.iono_delay;

            let res_pr = sat.p_meas - expected_pr;
            
            static mut PPP_PRINTED: bool = false;
            if unsafe { !PPP_PRINTED } {
                unsafe { PPP_PRINTED = true; }
                println!("PPP PRN{} PR res={:.3}m (meas={:.3}, expected_pr={:.3}, expected_base={:.3}, iono={:.3}, dt_sat_m={:.3}, tropo={:.3}, dist={:.3}, x_i[15]={:.3})", 
                    sat.sat_obs.sat.prn, res_pr, sat.p_meas, expected_pr, expected_base, sat.iono_delay, sat.dt_sat_m, sat.tropo_dry, dist, x_i[15]);
            }
            let var_pr = PSEUDORANGE_VARIANCE_BASE * 10.0 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
            let w_pr = apply_huber(res_pr, var_pr, self.huber_k);
            meas.push(FgMeasurement { res: res_pr, h_row: build_h_row(&los, sat.map_wet, None, x_i.len()), weight: w_pr, is_phase: false });

            if sat.doppler != 0.0 {
                let rcv_vel = Vector3::new(x_i[3], x_i[4], x_i[5]);
                let rcv_clk_drift = x_i[16];
                let meas_rr = -sat.doppler * sat.lam1;
                let expected_rr = los.dot(&sat.sat_vel) - los.dot(&rcv_vel) + rcv_clk_drift - sat.sat_clock_drift * SPEED_OF_LIGHT;
                
                let res_rr = meas_rr - expected_rr;
                if sat.sat_obs.sat.prn == 1 {
                    println!("PRN1 doppler res_rr = {:.3} m/s (meas={:.3}, expected={:.3})", res_rr, meas_rr, expected_rr);
                }
                let var_rr = 0.01; // Decreased Doppler variance (trust velocity more)
                let w_rr = apply_huber(res_rr, var_rr, 3.0);
                meas.push(FgMeasurement { res: res_rr, h_row: build_h_row_doppler(&los, x_i.len()), weight: w_rr, is_phase: false });
            }

            if let Some(cp1) = sat.cp1 {
                if cp1 == 0.0 { continue; }
                if let Some(amb_idx) = find_ambiguity_index(state, sat.sat_obs.sat) {
                    let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
                    let l_meas = if sat.is_iono_free && sat.cp2.is_some() {
                        crate::combinations::iono_free((cp1 + windup) * sat.lam1, (sat.cp2.unwrap() + windup) * sat.lam2, sat.f1, sat.f2)
                    } else {
                        (cp1 + windup) * sat.lam1
                    };
                    let expected_cp = if sat.is_iono_free && sat.cp2.is_some() {
                        expected_base + x_i[CORE_STATE_SIZE + amb_idx]
                    } else {
                        expected_base - sat.iono_delay + x_i[CORE_STATE_SIZE + amb_idx]
                    };
                    let res_cp = l_meas - expected_cp;
                    
                    if res_cp.abs() > 100.0 && iter == 0 {
                        tracing::warn!("HUGE res_cp: sat={}, l_meas={:.2}, exp={:.2}, dist={:.2}, clk={:.2}, n_amb={:.2}", 
                            sat.sat_obs.sat, l_meas, expected_cp, dist, x_i[15], x_i[CORE_STATE_SIZE + amb_idx]);
                    }
                    
                    let var_cp = 0.0001;
                    let var_cp = var_cp * snr_scale(sat.snr as i32) / libm::sin(sat.el);
                    
                    // Removed hard rejection of res_cp > 5.0m to allow convergence
                    
                    let w_cp = apply_huber(res_cp, var_cp, self.huber_k);
                    
                    meas.push(FgMeasurement { res: res_cp, h_row: build_h_row(&los, sat.map_wet, Some(CORE_STATE_SIZE + amb_idx), x_i.len()), weight: w_cp, is_phase: true });
                }
            }
        }
        meas
    }
}

fn log_ppp_convergence(state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>, x_pred: &DVector<f64>, p_pred: &DMatrix<f64>, solver: &PppFactorGraph) {
    tracing::info!("PPP Epoch: pos=[{:.2}, {:.2}, {:.2}], vel=[{:.2}, {:.2}, {:.2}], clk_d={:.2}, meas={}, P_v={:.2}, P_c={:.2}, dx_norm={:.4}", 
        x_i[0], x_i[1], x_i[2], x_i[3], x_i[4], x_i[5], x_i[16],
        sats.len(), p_pred[(3, 3)], p_pred[(16, 16)], (x_i.clone() - x_pred.clone()).norm());
    
    let pos_err = (Vector3::new(x_i[0], x_i[1], x_i[2]) - Vector3::new(x_pred[0], x_pred[1], x_pred[2])).norm();
    if pos_err > 100.0 {
        tracing::warn!("x_pred pos: {:.2}, {:.2}, {:.2}, clk: {:.2}", x_pred[0], x_pred[1], x_pred[2], x_pred[15]);
        tracing::warn!("x_i pos: {:.2}, {:.2}, {:.2}, clk: {:.2}", x_i[0], x_i[1], x_i[2], x_i[15]);
    }
    let vel_norm = Vector3::new(x_i[3], x_i[4], x_i[5]).norm();
    if vel_norm > 100.0 {
        tracing::warn!("Crazy velocity: {:.2} m/s", vel_norm);
    }
    
    let meas = solver.build_measurements(state, sats, x_i, solver.max_iterations - 1);
    let mut sum_pr = 0.0;
    let mut count_pr = 0;
    let mut sum_rr = 0.0;
    let mut count_rr = 0;
    for m in &meas {
        if m.h_row.len() == x_i.len() { // valid
            if m.is_phase { continue; } // Phase was disabled
            if m.weight > 0.05 { // It's a PR (weight = var, Huber makes weight smaller?) wait, weight is variance!
                sum_pr += m.res.abs(); count_pr += 1;
            } else {
                sum_rr += m.res.abs(); count_rr += 1;
            }
        }
    }
    if state.epoch_count.is_multiple_of(100) {
        println!("Epoch {}: Mean PR Res = {:.3} m, Mean RR Res = {:.3} m/s", state.epoch_count, sum_pr / count_pr.max(1) as f64, sum_rr / count_rr.max(1) as f64);
    }
}

fn build_weight_matrix(meas: &[FgMeasurement], r_mat: &DMatrix<f64>) -> DMatrix<f64> {
    let mut w_mat = DMatrix::zeros(meas.len(), meas.len());
    for i in 0..meas.len() {
        w_mat[(i, i)] = 1.0 / r_mat[(i, i)];
    }
    w_mat
}

fn apply_huber(res: f64, var: f64, k: f64) -> f64 {
    let norm_r = res.abs() / var.sqrt();
    if norm_r > k { var * (norm_r / k) } else { var }
}

fn snr_scale(snr: i32) -> f64 {
    (10.0f64).powf((NOMINAL_SNR_DBHZ - snr as f64) / SNR_SCALE_DIVISOR)
}

fn invert_matrix(mat: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    if mat.iter().any(|x| x.is_nan()) { return None; }
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

fn build_h_row_doppler(los: &Vector3<f64>, size: usize) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    if size > 16 {
        h[3] = -los.x; h[4] = -los.y; h[5] = -los.z;
        h[16] = 1.0;
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, DVector, Vector3};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0)),
            1.0
        )
    }

    #[test]
    fn test_solve_empty_sats() {
        let fg = PppFactorGraph::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let sats = vec![];
        let res = fg.solve(&mut state, &sats);
        assert!(matches!(res, Err(EngineError::InsufficientSatellites)));
    }

    #[test]
    fn test_ppp_factor_graph_default() {
        let fg = PppFactorGraph::default();
        assert_eq!(fg.max_iterations, 15);
        assert_eq!(fg.convergence_threshold, 1e-3);
        assert_eq!(fg.huber_k, 3.0);
        let fg2 = PppFactorGraph::new();
        assert_eq!(fg2.max_iterations, 15);
    }

    #[test]
    fn test_find_ambiguity_index() {
        let mut state = dummy_rtk_state();
        use gneiss_core::sat::{SatelliteId, Constellation};
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        state.ambiguity_keys.push((sat1, 0));
        state.ambiguity_keys.push((sat2, 1));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);
        
        assert_eq!(find_ambiguity_index(&state, sat1), Some(0));
        assert_eq!(find_ambiguity_index(&state, sat2), None);
        
        let sat3 = SatelliteId { constellation: Constellation::Gps, prn: 3 };
        assert_eq!(find_ambiguity_index(&state, sat3), None);
    }

    #[test]
    fn test_build_weight_matrix() {
        let m1 = FgMeasurement { res: 1.0, h_row: DVector::zeros(1), weight: 2.0, is_phase: false };
        let m2 = FgMeasurement { res: 2.0, h_row: DVector::zeros(1), weight: 4.0, is_phase: true };
        let meas = vec![m1, m2];
        let mut r = DMatrix::zeros(2, 2);
        r[(0, 0)] = 2.0; r[(1, 1)] = 4.0;
        let w = build_weight_matrix(&meas, &r);
        assert_eq!(w[(0, 0)], 0.5);
        assert_eq!(w[(1, 1)], 0.25);
        assert_eq!(w[(0, 1)], 0.0);
    }

    #[test]
    fn test_apply_huber() {
        let res = apply_huber(1.0, 1.0, 3.0);
        assert_eq!(res, 1.0);
        
        let res = apply_huber(-1.0, 1.0, 3.0);
        assert_eq!(res, 1.0);
        
        let res = apply_huber(4.0, 1.0, 3.0);
        assert!((res - 1.3333333333333333).abs() < 1e-10);
        
        let res = apply_huber(-6.0, 4.0, 2.0);
        assert!((res - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_snr_scale() {
        assert!((snr_scale(45) - 1.0).abs() < 1e-10);
        assert!((snr_scale(35) - 10.0).abs() < 1e-10);
        assert!((snr_scale(55) - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_invert_matrix() {
        let mut mat = DMatrix::zeros(2, 2);
        mat[(0,0)] = 2.0; mat[(1,1)] = 2.0;
        let inv = invert_matrix(&mat).unwrap();
        assert!((inv[(0,0)] - 0.5).abs() < 1e-10);
        assert!((inv[(1,1)] - 0.5).abs() < 1e-10);
        assert!((inv[(0,1)]).abs() < 1e-10);
        assert!((inv[(1,0)]).abs() < 1e-10);
        
        let mat2 = DMatrix::zeros(2, 2);
        let inv2 = invert_matrix(&mat2).unwrap();
        assert!((inv2[(0,0)] - 1e6).abs() < 1e-1);
    }

    #[test]
    fn test_build_h_row() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let h = build_h_row(&los, 4.0, Some(18), 19);
        assert_eq!(h.len(), 19);
        assert_eq!(h[0], -1.0);
        assert_eq!(h[1], -2.0);
        assert_eq!(h[2], -3.0);
        assert_eq!(h[15], 1.0);
        assert_eq!(h[17], 4.0);
        assert_eq!(h[18], 1.0);
        
        let h2 = build_h_row(&los, 4.0, None, 16);
        assert_eq!(h2.len(), 16);
        assert_eq!(h2[0], -1.0);
        assert_eq!(h2[15], 1.0);
    }

    #[test]
    fn test_build_h_row_doppler() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let h = build_h_row_doppler(&los, 17);
        assert_eq!(h.len(), 17);
        assert_eq!(h[3], -1.0);
        assert_eq!(h[4], -2.0);
        assert_eq!(h[5], -3.0);
        assert_eq!(h[16], 1.0);
        
        let h2 = build_h_row_doppler(&los, 16);
        assert_eq!(h2.len(), 16);
        assert_eq!(h2[3], 0.0);
    }

    #[test]
    fn test_assemble_matrices() {
        let m1 = FgMeasurement {
            res: 1.5,
            h_row: DVector::from_element(3, 1.0),
            weight: 2.0,
            is_phase: false,
        };
        let m2 = FgMeasurement {
            res: 2.5,
            h_row: DVector::from_element(3, 2.0),
            weight: 3.0,
            is_phase: true,
        };
        let meas = vec![m1, m2];
        let (h, z, r) = assemble_matrices(&meas, 3);
        
        assert_eq!(h.nrows(), 2);
        assert_eq!(h.ncols(), 3);
        assert_eq!(h[(0,0)], 1.0);
        assert_eq!(h[(1,2)], 2.0);
        
        assert_eq!(z.len(), 2);
        assert_eq!(z[0], 1.5);
        assert_eq!(z[1], 2.5);
        
        assert_eq!(r.nrows(), 2);
        assert_eq!(r.ncols(), 2);
        assert_eq!(r[(0,0)], 2.0);
        assert_eq!(r[(1,1)], 3.0);
        assert_eq!(r[(0,1)], 0.0);
    }

    #[test]
    fn test_extract_and_apply_state_vector() {
        let mut state = dummy_rtk_state();
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.attitude = nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(0.1, 0.2, 0.3));
        state.accel_bias = Vector3::new(10.0, 11.0, 12.0);
        state.gyro_bias = Vector3::new(13.0, 14.0, 15.0);
        state.rcv_clk_bias = 16.0;
        state.rcv_clk_drift = 17.0;
        state.zwd = 18.0;
        state.ambiguities = vec![19.0, 20.0];
        
        let x = extract_state_vector(&state);
        assert_eq!(x.len(), CORE_STATE_SIZE + 2);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[15], 16.0);
        assert_eq!(x[17], 18.0);
        assert_eq!(x[CORE_STATE_SIZE], 19.0);
        assert_eq!(x[CORE_STATE_SIZE + 1], 20.0);
        
        let mut state2 = dummy_rtk_state();
        state2.ambiguities = vec![0.0, 0.0];
        apply_state_vector(&mut state2, &x);
        
        assert_eq!(state2.position.vector, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(state2.velocity, Vector3::new(4.0, 5.0, 6.0));
        assert!((state2.attitude.scaled_axis() - Vector3::new(0.1, 0.2, 0.3)).norm() < 1e-10);
        assert_eq!(state2.accel_bias, Vector3::new(10.0, 11.0, 12.0));
        assert_eq!(state2.gyro_bias, Vector3::new(13.0, 14.0, 15.0));
        assert_eq!(state2.rcv_clk_bias, 16.0);
        assert_eq!(state2.rcv_clk_drift, 17.0);
        assert_eq!(state2.zwd, 18.0);
        assert_eq!(state2.ambiguities, vec![19.0, 20.0]);
    }
}
#[cfg(test)]
mod nan_tests {
    use super::*;
    use nalgebra::{DMatrix, DVector, Vector3};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0)),
            1.0
        )
    }

    #[test]
    fn test_solve_matrix_inversion_failure() {
        let fg = PppFactorGraph::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::from_element(CORE_STATE_SIZE, CORE_STATE_SIZE, f64::NAN);
        let sats = vec![];
        let res = fg.solve(&mut state, &sats);
        assert!(matches!(res, Err(EngineError::StateDisappeared)));
    }
}
