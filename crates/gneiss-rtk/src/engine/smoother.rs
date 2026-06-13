use crate::engine::{ProcessingEngine, EngineMode, EngineError};
use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector};

pub fn run_combined_ppk(engine: &mut ProcessingEngine) -> Result<Vec<RtkState>, EngineError> {
    let n_epochs = engine.state_history.len();
    if n_epochs == 0 { return Err(EngineError::NoObservations); }
    
    let mut smoothed_states = engine.state_history.clone();

    if matches!(engine.config.mode, EngineMode::Spp | EngineMode::Ppp | EngineMode::PppIns | EngineMode::PppInsLooselyCoupled) {
        return Ok(smoothed_states);
    }
    
    for k in (0..n_epochs - 1).rev() {
        if smoothed_states[k+1].is_reset {
            tracing::debug!("Epoch {} was reset. Breaking smoothing chain at k={}", k+1, k);
            continue;
        }
        
        let phi_k = match &smoothed_states[k+1].core_phi {
            Some(p) => p.clone(),
            None => continue,
        };
        let p_pred_k1 = match &smoothed_states[k+1].full_p_predict {
            Some(p) => p.clone(),
            None => continue,
        };
        let x_pred_k1 = match &smoothed_states[k+1].full_x_predict {
            Some(x) => x.clone(),
            None => continue,
        };

        let (left, right) = smoothed_states.split_at_mut(k + 1);
        let state_k = &mut left[k];
        let state_k1 = &right[0];

        if let Err(e) = smooth_epoch(state_k, state_k1, &phi_k, &p_pred_k1, &x_pred_k1, k) {
            tracing::debug!("RTS Smoothing skipped at k={}: {}", k, e);
            continue;
        }

        smoothed_states[k].is_fixed = false;
        smoothed_states[k].fixed_state = None;
        if !matches!(engine.config.mode, EngineMode::Spp | EngineMode::SppIns) {
            let ephemerides = &engine.ephemerides;
            let subset = engine.config.lambda_min_subset;
            let epoch_count = engine.config.ar_min_epoch_count;
            let min_lock = engine.config.ar_min_lock;
            let min_ratio = engine.config.lambda_min_ratio;
            if let Ok((fixed_state, _, _, _, _)) = smoothed_states[k].resolve_ambiguities(ephemerides, subset, epoch_count, min_lock, min_ratio) {
                smoothed_states[k].fixed_state = Some(Box::new(fixed_state));
            }
        }
    }
    Ok(smoothed_states)
}

fn smooth_epoch(
    state_k: &mut RtkState, 
    state_k1: &RtkState, 
    phi_k: &DMatrix<f64>, 
    p_pred_k1: &DMatrix<f64>, 
    x_pred_k1: &DVector<f64>,
    k_idx: usize
) -> Result<(), &'static str> {
    let core_size = crate::filter::CORE_STATE_SIZE;
    
    let (matched_k_indices, matched_k1_indices) = find_matched_ambiguities(state_k, state_k1);
    let smooth_len = core_size + matched_k_indices.len();
    
    let mut idx_k = (0..core_size).collect::<Vec<_>>();
    idx_k.extend(matched_k_indices.iter().copied());
    
    let mut idx_k1 = (0..core_size).collect::<Vec<_>>();
    idx_k1.extend(&matched_k1_indices);

    let x_k1_n = build_x_vector(state_k1, core_size, smooth_len, &matched_k1_indices);
    let p_k1_n = extract_submatrix(&state_k1.covariance, &idx_k1, &idx_k1);
    let p_k = extract_submatrix(&state_k.covariance, &idx_k, &idx_k);
    let p_pred_k1_sub = extract_submatrix(p_pred_k1, &idx_k, &idx_k);
    let phi_k_sub = build_phi_submatrix(phi_k, core_size, smooth_len);

    if p_pred_k1_sub.iter().any(|&x| !x.is_finite() || x.abs() > 1e10) || p_k.iter().any(|&x| !x.is_finite() || x.abs() > 1e10) {
        return Err("non-finite covariance");
    }

    let p_pred_inv = invert_p_pred(&p_pred_k1_sub, smooth_len)?;

    let x_pred_k1_sub = extract_subvector(x_pred_k1, &idx_k);
    let c_k = &p_k * phi_k_sub.transpose() * p_pred_inv;
    let x_k = build_x_vector(state_k, core_size, smooth_len, &matched_k_indices);
    
    let correction = &c_k * (x_k1_n - x_pred_k1_sub);
    let pos_corr_norm = correction.fixed_rows::<3>(0).norm();
    
    if pos_corr_norm > 10.0 {
        tracing::warn!("Smoother large position correction: {:.1}m at k={}", pos_corr_norm, k_idx);
    }
    
    let x_k_n = x_k + correction;
    let p_k_n = p_k + &c_k * (p_k1_n - p_pred_k1_sub) * c_k.transpose();

    update_smoothed_state(state_k, &x_k_n, &p_k_n, core_size, smooth_len, &matched_k_indices, &idx_k);
    Ok(())
}

fn find_matched_ambiguities(state_k: &RtkState, state_k1: &RtkState) -> (Vec<usize>, Vec<usize>) {
    let mut matched_k = Vec::new();
    let mut matched_k1 = Vec::new();
    for (i, key_k) in state_k.ambiguity_keys.iter().enumerate() {
        if let Some(j) = state_k1.ambiguity_keys.iter().position(|k| k == key_k) {
            if state_k.ambiguity_track_ids[i] == state_k1.ambiguity_track_ids[j] {
                let cov_idx = crate::filter::CORE_STATE_SIZE + j;
                if state_k1.covariance[(cov_idx, cov_idx)] < 10.0 {
                    matched_k.push(crate::filter::CORE_STATE_SIZE + i);
                    matched_k1.push(cov_idx);
                }
            }
        }
    }
    (matched_k, matched_k1)
}

fn build_x_vector(state: &RtkState, core_size: usize, len: usize, matched_indices: &[usize]) -> DVector<f64> {
    let mut x = DVector::zeros(len);
    x.rows_mut(0, 3).copy_from(&state.position.vector);
    x.rows_mut(3, 3).copy_from(&state.velocity);
    if core_size > 6 {
        x.rows_mut(9, 3).copy_from(&state.accel_bias);
        x.rows_mut(12, 3).copy_from(&state.gyro_bias);
    }
    if core_size > 15 {
        x[15] = state.rcv_clk_bias;
        x[16] = state.isb_glo;
        x[17] = state.isb_gal;
        x[18] = state.isb_bds;
        x[19] = state.rcv_clk_drift;
        x[20] = state.zwd;
    }
    for (i, &idx) in matched_indices.iter().enumerate() {
        x[core_size + i] = state.ambiguities[idx - crate::filter::CORE_STATE_SIZE];
    }
    x
}

fn extract_submatrix(mat: &DMatrix<f64>, rows: &[usize], cols: &[usize]) -> DMatrix<f64> {
    let mut sub = DMatrix::zeros(rows.len(), cols.len());
    for (i, &r) in rows.iter().enumerate() {
        for (j, &c) in cols.iter().enumerate() {
            sub[(i, j)] = mat[(r, c)];
        }
    }
    sub
}

fn extract_subvector(vec: &DVector<f64>, indices: &[usize]) -> DVector<f64> {
    let mut sub = DVector::zeros(indices.len());
    for (i, &idx) in indices.iter().enumerate() {
        sub[i] = vec[idx];
    }
    sub
}

fn build_phi_submatrix(phi: &DMatrix<f64>, core_size: usize, len: usize) -> DMatrix<f64> {
    let mut sub = DMatrix::zeros(len, len);
    for i in 0..core_size {
        for j in 0..core_size {
            sub[(i, j)] = phi[(i, j)];
        }
    }
    for i in core_size..len {
        sub[(i, i)] = 1.0;
    }
    sub
}

fn invert_p_pred(p_pred: &DMatrix<f64>, len: usize) -> Result<DMatrix<f64>, &'static str> {
    let active: Vec<usize> = (0..len).filter(|&i| p_pred[(i, i)] > 1e-12).collect();
    let m = active.len();
    if m == len {
        let reg = DMatrix::identity(len, len) * 1e-9;
        (p_pred + reg).try_inverse().ok_or("inversion failed")
    } else if m > 0 {
        let p_act = extract_submatrix(p_pred, &active, &active);
        let reg = DMatrix::identity(m, m) * 1e-9;
        let inv_act = (p_act + reg).try_inverse().ok_or("inversion failed")?;
        let mut inv_full = DMatrix::zeros(len, len);
        for (i, &r) in active.iter().enumerate() {
            for (j, &c) in active.iter().enumerate() {
                inv_full[(r, c)] = inv_act[(i, j)];
            }
        }
        Ok(inv_full)
    } else {
        Err("no active elements")
    }
}

fn update_smoothed_state(
    state: &mut RtkState, 
    x_k_n: &DVector<f64>, 
    p_k_n: &DMatrix<f64>, 
    core_size: usize, 
    smooth_len: usize,
    matched_indices: &[usize],
    idx_k: &[usize]
) {
    state.position.vector = x_k_n.fixed_rows::<3>(0).into_owned();
    state.velocity = x_k_n.fixed_rows::<3>(3).into_owned();
    if core_size > 6 {
        let d_theta = x_k_n.fixed_rows::<3>(6).into_owned();
        if d_theta.norm() > 1e-10 {
            let dq = nalgebra::UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_normalize(d_theta), d_theta.norm());
            state.attitude = dq * state.attitude;
            state.attitude.renormalize();
        }
        state.accel_bias = x_k_n.fixed_rows::<3>(9).into_owned();
        state.gyro_bias = x_k_n.fixed_rows::<3>(12).into_owned();
    }
    if core_size > 15 {
        state.rcv_clk_bias = x_k_n[15];
        state.isb_glo = x_k_n[16];
        state.isb_gal = x_k_n[17];
        state.isb_bds = x_k_n[18];
        state.rcv_clk_drift = x_k_n[19];
        state.zwd = x_k_n[20];
    }
    for (i, &idx) in matched_indices.iter().enumerate() {
        let amb_idx = idx - crate::filter::CORE_STATE_SIZE;
        state.ambiguities[amb_idx] = x_k_n[core_size + i];
    }
    for i in 0..smooth_len {
        for j in 0..smooth_len {
            state.covariance[(idx_k[i], idx_k[j])] = p_k_n[(i, j)];
        }
    }
}
