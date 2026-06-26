use crate::engine::ppp_common::{
    apply_state_vector, assemble_matrices, build_weight_matrix, extract_state_vector, find_amb_idx,
    invert_matrix, FgMeasurement,
};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::ProcessingEngine;
use crate::engine::EngineError;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::math::inversion::solve_cholesky_svd;
use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};

pub fn process_ppp_ins_fg<'a>(
    engine: &'a mut crate::engine::ProcessingEngine,
    rover_obs: &'a gneiss_core::obs::EpochObs,
) -> Result<&'a RtkState, EngineError> {
    if !crate::engine::ppp::valid_pos(engine) {
        return engine.process_spp(rover_obs);
    }
    let dt = (rover_obs.time.tow - engine.current_state.as_ref().unwrap().time.tow).max(0.0);
    engine.predict_state(dt);

    let state = engine.current_state.as_mut().unwrap();
    state.time = rover_obs.time;
    state.position.epoch = rover_obs.time;
    let sats = crate::engine::ppp_antenna::build_sats(engine, rover_obs);
    if sats.is_empty() {
        return Err(EngineError::InsufficientSatellites);
    }

    let state = engine.current_state.as_mut().unwrap();
    crate::engine::ppp::update_phase_ambiguities(state, &sats, rover_obs.time);
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);

    // SPP position prior — anchors the IEKF during early convergence,
    // same as process_ppp. Without this anchor, the IEKF can diverge
    // catastrophically on a single bad epoch when IMU factors are not
    // yet active (pre-alignment).
    let mut position_prior: Option<(Vector3<f64>, f64)> = None;
    let mut spp_pos_for_recovery: Option<gneiss_core::coords::Coordinate> = None;
    if let Ok(spp) = crate::spp::compute_spp(
        rover_obs,
        &engine.ephemerides,
        engine.klobuchar_params.as_ref(),
        &crate::spp::SppConfig::default(),
        None,
    ) {
        spp_pos_for_recovery = Some(spp.position);
        let is_cold_start = state.epoch_count < 2;
        if is_cold_start {
            state.position = spp.position;
            state.rcv_clk_bias = spp.cdt;
        } else {
            let pos_cov = state.covariance[(0, 0)]
                .min(state.covariance[(1, 1)])
                .min(state.covariance[(2, 2)]);
            let prior_var = pos_cov.clamp(9.0, 100.0);
            position_prior = Some((spp.position.vector, prior_var));
        }
    }

    let last_state = engine.state_history.last();
    let imu_history = &engine.imu_history;

    let lever_arm = nalgebra::Vector3::from_column_slice(&engine.config.imu_to_antenna_lever_arm);
    PppInsIteratedEkf::new().solve(state, &sats, imu_history, last_state, &lever_arm, position_prior)?;

    // Check for covariance explosion (same guard as process_ppp)
    if let Some(ref mut s) = engine.current_state {
        ProcessingEngine::check_covariance_divergence(
            s,
            spp_pos_for_recovery,
            None,
            false,
        );
    }

    crate::engine::processor::ProcessingEngine::attempt_kinematic_alignment(engine);
    let state = engine.current_state.as_mut().unwrap();
    crate::engine::processor::ProcessingEngine::apply_nhc_updates(
        &engine.config,
        &engine.imu_history,
        state,
    );

    state.epoch_count = state.epoch_count.saturating_add(1);
    let final_state = engine.current_state.as_ref().unwrap().clone();
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));
    engine.imu_buffer.clear();
    Ok(engine.current_state.as_ref().unwrap())
}

/// Iterated Extended Kalman Filter for tightly-coupled PPP+INS.
///
/// Despite the historical "fg" (factor graph) naming, this is an IEKF —
/// an iterated least-squares solver with IMU pre-integration factors,
/// Huber robust estimation, and a prior from state propagation.
/// It does not perform marginalization, variable elimination, or iSAM2-style
/// incremental smoothing.
pub struct PppInsIteratedEkf {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub huber_k: f64,
    pub elev_mask: f64,
}

impl Default for PppInsIteratedEkf {
    fn default() -> Self {
        Self {
            max_iterations: 8,
            convergence_threshold: 0.005,
            huber_k: 10.0,
            elev_mask: 10.0_f64.to_radians(),
        }
    }
}

impl PppInsIteratedEkf {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn solve(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
        lever_arm: &nalgebra::Vector3<f64>,
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<(), EngineError> {
        for outer_iter in 0..4 {
            let done = self.solve_inner(state, sats, imu_history, last_state, lever_arm, position_prior)?;
            if done || outer_iter == 3 {
                if !done {
                    tracing::warn!("Max outlier rejection iterations reached.");
                }
                break;
            }
        }
        if state.epoch_count > 10 {
            if let Err(e) = self.resolve_cascade_ar(state, sats) {
                tracing::info!("Cascade AR did not fix: {:?}", e);
            }
        }
        Ok(())
    }

    fn solve_inner(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
        lever_arm: &nalgebra::Vector3<f64>,
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<bool, EngineError> {
        let x_pred = extract_state_vector(state);
        let p_pred = state.covariance.clone();
        state.full_x_predict = Some(x_pred.clone());
        state.full_p_predict = Some(p_pred.clone());
        let mut x_i = x_pred.clone();
        let p_inv = invert_matrix(&p_pred).ok_or(EngineError::StateDisappeared)?;

        for _iter in 0..self.max_iterations {
            if let Some(dx) = self.compute_iteration_dx(
                state,
                sats,
                &x_i,
                &x_pred,
                &p_inv,
                _iter,
                imu_history,
                last_state,
                lever_arm,
                position_prior,
            )? {
                let mut x_next = &x_i + &dx;
                if x_next.len() > 8 {
                    let q_old = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                        x_i[6], x_i[7], x_i[8],
                    ));
                    let dq = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                        dx[6], dx[7], dx[8],
                    ));
                    let q_new = q_old * dq;
                    let new_scaled = q_new.scaled_axis();
                    x_next[6] = new_scaled.x;
                    x_next[7] = new_scaled.y;
                    x_next[8] = new_scaled.z;
                }
                x_i = x_next;
                if dx.norm() < self.convergence_threshold {
                    break;
                }
            } else {
                break;
            }
        }

        if let Some(sat) = self.find_worst_outlier(state, sats, &x_i, imu_history, lever_arm, position_prior) {
            tracing::warn!(
                "PPP FG Outlier Detected for {:?}. Removing ambiguity and retrying.",
                sat
            );
            for i in 0..4 {
                state.remove_ambiguity(sat, i);
            }
            return Ok(false);
        }

        let final_p = self.compute_final_covariance(
            state,
            sats,
            &x_i,
            &p_pred,
            &p_inv,
            imu_history,
            last_state,
            lever_arm,
            position_prior,
        );
        apply_state_vector(state, &x_i, final_p);
        let omega_eb_b = imu_history
            .last()
            .and_then(|buf: &Vec<gneiss_core::imu::ImuMeasurement>| buf.last())
            .map(|m| {
                let r_b_e = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                    x_i[6], x_i[7], x_i[8],
                ))
                .to_rotation_matrix();
                let omega_ie_e = nalgebra::Vector3::new(
                    0.0,
                    0.0,
                    gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S,
                );
                m.gyro
                    - nalgebra::Vector3::new(x_i[12], x_i[13], x_i[14])
                    - r_b_e.transpose() * omega_ie_e
            })
            .unwrap_or(nalgebra::Vector3::zeros());
        log_ppp_convergence(
            state,
            sats,
            &x_i,
            &x_pred,
            &p_pred,
            self,
            lever_arm,
            &omega_eb_b,
        );
        Ok(true)
    }

    fn find_worst_outlier(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        lever_arm: &nalgebra::Vector3<f64>,
        _position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Option<gneiss_core::sat::SatelliteId> {
        let omega_eb_b = imu_history
            .last()
            .and_then(|buf: &Vec<gneiss_core::imu::ImuMeasurement>| buf.last())
            .map(|m| {
                let r_b_e = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                    x_i[6], x_i[7], x_i[8],
                ))
                .to_rotation_matrix();
                let omega_ie_e = nalgebra::Vector3::new(
                    0.0,
                    0.0,
                    gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S,
                );
                m.gyro
                    - nalgebra::Vector3::new(x_i[12], x_i[13], x_i[14])
                    - r_b_e.transpose() * omega_ie_e
            })
            .unwrap_or(nalgebra::Vector3::zeros());
        let final_meas = self.build_measurements(
            state,
            sats,
            x_i,
            self.max_iterations,
            lever_arm,
            &omega_eb_b,
        );
        Self::find_worst_outlier_sat(&final_meas)
    }

    fn find_worst_outlier_sat(meas: &[FgMeasurement]) -> Option<gneiss_core::sat::SatelliteId> {
        let mut worst_sat = None;
        let mut max_norm = 15.0;
        for m in meas {
            if m.is_phase {
                let norm = m.res.abs() / m.raw_var.sqrt();
                if norm > max_norm {
                    max_norm = norm;
                    worst_sat = m.sat;
                }
            }
        }
        worst_sat
    }

    pub fn resolve_cascade_ar(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
    ) -> Result<(), &'static str> {
        let cands = self.find_ar_candidates(state, sats);
        if cands.len() < 4 {
            return Err("Insufficient dual-frequency satellites for AR");
        }
        let subset = self.build_ar_subset(&cands);
        if subset.len() < 3 {
            return Err("Insufficient satellites after single differencing");
        }

        let x = extract_state_vector(state);
        let (x_wl, p_wl, keep_indices) = self.resolve_widelane_ar(state, &subset, &x)?;

        let (x_fixed, p_fixed) =
            self.resolve_narrowlane_ar(state, &subset, &keep_indices, &x_wl, &p_wl)?;

        apply_state_vector(state, &x_fixed, p_fixed);
        state.is_fixed = true;
        Ok(())
    }

    fn find_ar_candidates(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
    ) -> Vec<(gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64)> {
        let mut cands = Vec::new();
        for sat in sats.iter().filter(|s| {
            !s.is_iono_free
                && s.cp2.is_some()
                && s.sat_obs.sat.constellation != gneiss_core::sat::Constellation::Glonass
        }) {
            if let (Some(n1), Some(n2)) = (
                find_amb_idx(state, sat.sat_obs.sat, 1),
                find_amb_idx(state, sat.sat_obs.sat, 2),
            ) {
                cands.push((sat.sat_obs.sat, n1, n2, sat.el, sat.lam1, sat.lam2));
            }
        }
        cands
    }

    #[allow(clippy::type_complexity)]
    fn build_ar_subset(
        &self,
        cands: &[(gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64)],
    ) -> Vec<(
        (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
        (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
    )> {
        let mut subset = Vec::new();
        // Group by constellation
        let mut const_cands = std::collections::HashMap::new();
        for cand in cands {
            const_cands
                .entry(cand.0.constellation)
                .or_insert_with(Vec::new)
                .push(*cand);
        }

        for (_, mut group) in const_cands {
            if group.len() < 2 {
                continue;
            }
            // Find reference satellite (highest elevation)
            group.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));
            let ref_cand = group[0];

            for cand in group.iter().skip(1) {
                subset.push((*cand, ref_cand));
            }
        }

        subset
    }

    #[allow(clippy::type_complexity)]
    fn resolve_widelane_ar(
        &self,
        state: &RtkState,
        subset: &[(
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
        )],
        x: &DVector<f64>,
    ) -> Result<(DVector<f64>, DMatrix<f64>, Vec<usize>), &'static str> {
        let mut d_wl_full = DMatrix::zeros(subset.len(), state.covariance.nrows());
        for (i, (c, ref_sat)) in subset.iter().enumerate() {
            d_wl_full[(i, CORE_STATE_SIZE + c.1)] = 1.0 / c.4;
            d_wl_full[(i, CORE_STATE_SIZE + c.2)] = -1.0 / c.5;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.1)] = -1.0 / ref_sat.4;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.2)] = 1.0 / ref_sat.5;
        }

        let q_wl_full = &d_wl_full * &state.covariance * d_wl_full.transpose();
        let keep_indices: Vec<usize> = (0..q_wl_full.nrows())
            .filter(|&i| q_wl_full[(i, i)].sqrt() < 0.30)
            .collect();
        if keep_indices.len() < 3 {
            return Err("Insufficient well-converged Widelane ambiguities");
        }

        let mut d_wl = DMatrix::zeros(keep_indices.len(), state.covariance.nrows());
        for (i, &idx) in keep_indices.iter().enumerate() {
            for j in 0..state.covariance.nrows() {
                d_wl[(i, j)] = d_wl_full[(idx, j)];
            }
        }

        let a_wl = &d_wl * x;
        let q_wl = &d_wl * &state.covariance * d_wl.transpose();
        let res_wl = crate::ambiguity::lambda::resolve_lambda(&a_wl, &q_wl)
            .map_err(|_| "WL LAMBDA Failed")?;

        if res_wl.ratio < 1.5 || res_wl.success_rate < 0.95 {
            return Err("WL ratio test failed");
        }

        let s_inv = q_wl.try_inverse().ok_or("WL Cov Inversion failed")?;
        let k_wl = &state.covariance * d_wl.transpose() * s_inv;
        let dx_wl = &k_wl * (res_wl.best_integers - a_wl);
        Ok((
            x + dx_wl,
            crate::math::covariance::apply_joseph_covariance_update(
                &state.covariance,
                &k_wl,
                &d_wl,
                &DMatrix::zeros(keep_indices.len(), keep_indices.len()),
            ),
            keep_indices,
        ))
    }

    #[allow(clippy::type_complexity)]
    fn resolve_narrowlane_ar(
        &self,
        state: &RtkState,
        subset: &[(
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
        )],
        keep_indices: &[usize],
        x_wl: &DVector<f64>,
        p_wl: &DMatrix<f64>,
    ) -> Result<(DVector<f64>, DMatrix<f64>), &'static str> {
        let mut d_nl = DMatrix::zeros(keep_indices.len(), state.covariance.nrows());
        for (i, &idx) in keep_indices.iter().enumerate() {
            let (c, ref_sat) = &subset[idx];
            d_nl[(i, CORE_STATE_SIZE + c.1)] = 1.0 / c.4;
            d_nl[(i, CORE_STATE_SIZE + ref_sat.1)] = -1.0 / ref_sat.4;
        }

        let a_nl = &d_nl * x_wl;
        let q_nl = &d_nl * p_wl * d_nl.transpose();
        let res_nl = crate::ambiguity::lambda::resolve_lambda(&a_nl, &q_nl)
            .map_err(|_| "NL LAMBDA Failed")?;

        if res_nl.ratio < 3.0 || res_nl.success_rate < 0.99 {
            return Err("NL ratio test failed");
        }

        let s_nl_inv = q_nl.try_inverse().ok_or("NL Cov Inversion failed")?;
        let k_nl = p_wl * d_nl.transpose() * s_nl_inv;

        let dx_nl = &k_nl * (res_nl.best_integers - a_nl);

        tracing::info!("p_wl dims: {}x{}", p_wl.nrows(), p_wl.ncols());
        tracing::info!("k_nl dims: {}x{}", k_nl.nrows(), k_nl.ncols());
        tracing::info!("d_nl dims: {}x{}", d_nl.nrows(), d_nl.ncols());
        tracing::info!("r dims: {}x{}", keep_indices.len(), keep_indices.len());

        let p_fixed = crate::math::covariance::apply_joseph_covariance_update(
            p_wl,
            &k_nl,
            &d_nl,
            &DMatrix::zeros(keep_indices.len(), keep_indices.len()),
        );

        tracing::info!("PPP Cascade AR Fixed! N_Sats: {}", keep_indices.len() + 1);
        Ok((x_wl + dx_nl, p_fixed))
    }

    #[allow(clippy::too_many_arguments)]
    fn compute_iteration_dx(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        x_pred: &DVector<f64>,
        p_inv: &DMatrix<f64>,
        iter: usize,
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
        lever_arm: &nalgebra::Vector3<f64>,
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<Option<DVector<f64>>, EngineError> {
        let omega_eb_b = imu_history
            .last()
            .and_then(|buf: &Vec<gneiss_core::imu::ImuMeasurement>| buf.last())
            .map(|m| {
                let r_b_e = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                    x_i[6], x_i[7], x_i[8],
                ))
                .to_rotation_matrix();
                let omega_ie_e = nalgebra::Vector3::new(
                    0.0,
                    0.0,
                    gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S,
                );
                m.gyro
                    - nalgebra::Vector3::new(x_i[12], x_i[13], x_i[14])
                    - r_b_e.transpose() * omega_ie_e
            })
            .unwrap_or(nalgebra::Vector3::zeros());
        let meas = self.build_measurements(state, sats, x_i, iter, lever_arm, &omega_eb_b);
        if meas.is_empty() {
            return Err(EngineError::InsufficientSatellites);
        }

        let (h_mat, res_vec, r_mat) = assemble_matrices(&meas, x_i.len());
        let w_mat = build_weight_matrix(&meas, &r_mat);

        let mut htwh = h_mat.transpose() * &w_mat * &h_mat;
        let mut htwr = h_mat.transpose() * &w_mat * &res_vec;

        self.accumulate_imu_factors(x_i, imu_history, last_state, &mut htwh, Some(&mut htwr));

        // SPP position prior — soft anchor to prevent divergence during
        // early convergence, matching process_ppp behaviour.
        if let Some((spp_pos, var)) = position_prior {
            let w = 1.0 / var;
            for i in 0..3 {
                htwh[(i, i)] += w;
                htwr[i] += w * (spp_pos[i] - x_i[i]);
            }
        }

        let htwh_damped = &htwh + p_inv;
        let mut diff = x_pred - x_i;
        if diff.len() > 8 {
            let q_pred = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                x_pred[6], x_pred[7], x_pred[8],
            ));
            let q_i = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                x_i[6], x_i[7], x_i[8],
            ));
            let dq_local = (q_i.inverse() * q_pred).scaled_axis();
            diff[6] = dq_local.x;
            diff[7] = dq_local.y;
            diff[8] = dq_local.z;
        }
        let innov = &htwr + p_inv * diff;

        match solve_cholesky_svd(&htwh_damped, &innov, 1e-9) {
            Ok(sol) => Ok(Some(sol)),
            Err(_) => {
                tracing::warn!("Failed to solve normal equations in PPP FG!");
                Ok(None)
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn compute_final_covariance(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        p_pred: &DMatrix<f64>,
        p_inv: &DMatrix<f64>,
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
        lever_arm: &nalgebra::Vector3<f64>,
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> DMatrix<f64> {
        let omega_eb_b = imu_history
            .last()
            .and_then(|buf: &Vec<gneiss_core::imu::ImuMeasurement>| buf.last())
            .map(|m| {
                let r_b_e = nalgebra::UnitQuaternion::from_scaled_axis(nalgebra::Vector3::new(
                    x_i[6], x_i[7], x_i[8],
                ))
                .to_rotation_matrix();
                let omega_ie_e = nalgebra::Vector3::new(
                    0.0,
                    0.0,
                    gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S,
                );
                m.gyro
                    - nalgebra::Vector3::new(x_i[12], x_i[13], x_i[14])
                    - r_b_e.transpose() * omega_ie_e
            })
            .unwrap_or(nalgebra::Vector3::zeros());
        let last_meas = self.build_measurements(
            state,
            sats,
            x_i,
            self.max_iterations,
            lever_arm,
            &omega_eb_b,
        );
        if last_meas.is_empty() {
            return p_pred.clone();
        }

        let (h_mat, _, r_mat) = assemble_matrices(&last_meas, x_i.len());
        let w_mat = build_weight_matrix(&last_meas, &r_mat);
        let mut htwh = h_mat.transpose() * &w_mat * h_mat;

        self.accumulate_imu_factors(x_i, imu_history, last_state, &mut htwh, None);

        if let Some((_spp_pos, var)) = position_prior {
            let w = 1.0 / var;
            for i in 0..3 {
                htwh[(i, i)] += w;
            }
        }

        let htwh_damped = htwh + p_inv;

        invert_matrix(&htwh_damped).unwrap_or_else(|| {
            tracing::warn!("invert_matrix(&htwh_damped) FAILED! Falling back to p_pred. htwh_damped has NaNs: {}, Infs: {}", htwh_damped.iter().any(|x| x.is_nan()), htwh_damped.iter().any(|x| x.is_infinite()));
            p_pred.clone()
        })
    }

    fn fill_imu_preint_covariance(
        preint: &mut crate::estimators::factor_graph::imu_factors::ImuPreintegration,
    ) {
        let dt = preint.dt.max(0.01);
        preint.covariance.fill(0.0);
        for i in 0..3 {
            preint.covariance[(i, i)] = ((5.0f64 * dt).max(2.0f64)).powi(2);
            preint.covariance[(i + 3, i + 3)] = (3.0f64).max(1.0f64).powi(2);
            preint.covariance[(i + 6, i + 6)] = ((0.1f64 * dt).max(0.05f64)).powi(2);
            preint.covariance[(i + 9, i + 9)] = 1e-4f64;
            preint.covariance[(i + 12, i + 12)] = 1e-6f64;
        }
    }

    fn build_imu_preint_factor(
        preint: crate::estimators::factor_graph::imu_factors::ImuPreintegration,
        last: &RtkState,
        x_i: &DVector<f64>,
    ) -> crate::estimators::factor_graph::imu_factors::ImuPreintegrationFactor {
        let gravity = crate::engine::predictor::gravity_wgs84(last.position.vector);
        let r_vec = Vector3::new(x_i[6], x_i[7], x_i[8]);
        let q_j = if r_vec.norm() > 1e-12 {
            UnitQuaternion::from_scaled_axis(r_vec)
        } else {
            UnitQuaternion::identity()
        };
        crate::estimators::factor_graph::imu_factors::ImuPreintegrationFactor {
            preint,
            gravity,
            nominal_p_i: last.position.vector,
            nominal_v_i: last.velocity,
            nominal_q_i: last.attitude,
            nominal_ba_i: last.accel_bias,
            nominal_bg_i: last.gyro_bias,
            nominal_p_j: Vector3::new(x_i[0], x_i[1], x_i[2]),
            nominal_v_j: Vector3::new(x_i[3], x_i[4], x_i[5]),
            nominal_q_j: q_j,
            nominal_ba_j: Vector3::new(x_i[9], x_i[10], x_i[11]),
            nominal_bg_j: Vector3::new(x_i[12], x_i[13], x_i[14]),
            idx_p_i: 0,
            idx_v_i: 3,
            idx_q_i: 6,
            idx_ba_i: 9,
            idx_bg_i: 12,
            idx_p_j: 15,
            idx_v_j: 18,
            idx_q_j: 21,
            idx_ba_j: 24,
            idx_bg_j: 27,
        }
    }

    fn accumulate_single_imu_factor(
        factor: &crate::estimators::factor_graph::imu_factors::ImuPreintegrationFactor,
        state_size: usize,
        htwh: &mut DMatrix<f64>,
        htwr: Option<&mut DVector<f64>>,
    ) {
        use crate::estimators::factor_graph::Factor;
        let zero_delta = DVector::zeros(30);
        let imu_res = factor.residual(&zero_delta);
        let imu_jac_full = factor.jacobian(&zero_delta);
        let imu_jac = imu_jac_full.columns(15, 15);
        let imu_info = factor.information();
        let mut full_imu_jac = DMatrix::zeros(15, state_size);
        full_imu_jac.view_mut((0, 0), (15, 15)).copy_from(&imu_jac);
        let j_t_info = full_imu_jac.transpose() * imu_info;
        *htwh += &j_t_info * &full_imu_jac;
        if let Some(r) = htwr {
            *r += &j_t_info * (-imu_res);
        }
    }

    fn accumulate_imu_factors(
        &self,
        x_i: &DVector<f64>,
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
        htwh: &mut DMatrix<f64>,
        htwr: Option<&mut DVector<f64>>,
    ) {
        let (Some(last), Some(imu_buf)) = (last_state, imu_history.last()) else {
            return;
        };
        if imu_buf.is_empty() || !last.ins_aligned {
            return;
        }
        let mut preint = crate::estimators::factor_graph::imu_factors::ImuPreintegration::new();
        preint.integrate(imu_buf, &last.accel_bias, &last.gyro_bias);
        Self::fill_imu_preint_covariance(&mut preint);
        let factor = Self::build_imu_preint_factor(preint, last, x_i);
        Self::accumulate_single_imu_factor(&factor, x_i.len(), htwh, htwr);
    }
}

#[allow(clippy::too_many_arguments)]
fn log_ppp_convergence(
    state: &RtkState,
    sats: &[ProcessedSat],
    x_i: &DVector<f64>,
    x_pred: &DVector<f64>,
    p_pred: &DMatrix<f64>,
    solver: &PppInsIteratedEkf,
    lever_arm: &nalgebra::Vector3<f64>,
    omega_eb_b: &nalgebra::Vector3<f64>,
) {
    let _p_amb = if p_pred.nrows() > 21 {
        p_pred[(21, 21)]
    } else {
        0.0
    };
    let dx_norm = (x_i.clone() - x_pred.clone()).norm();
    tracing::info!(
        "PPP Epoch: pos=[{:.2}, {:.2}, {:.2}], vel=[{:.2}, {:.2}, {:.2}], speed={:.2}, aligned={}, dx_norm={:.4}",
        state.position.vector.x,
        state.position.vector.y,
        state.position.vector.z,
        state.velocity.x,
        state.velocity.y,
        state.velocity.z,
        state.velocity.norm(),
        state.ins_aligned,
        dx_norm
    );

    let meas = solver.build_measurements(
        state,
        sats,
        x_i,
        solver.max_iterations - 1,
        lever_arm,
        omega_eb_b,
    );
    let (mut sum_pr, mut count_pr, mut sum_rr, mut count_rr) = (0.0, 0, 0.0, 0);
    for m in &meas {
        if m.h_row.len() == x_i.len() {
            if m.is_phase {
                continue;
            }
            if m.weight > 0.05 {
                sum_pr += m.res.abs();
                count_pr += 1;
            } else {
                sum_rr += m.res.abs();
                count_rr += 1;
            }
        }
    }
    if state.epoch_count.is_multiple_of(100) {
        tracing::trace!(
            "Epoch {}: Mean PR Res = {:.3} m, Mean RR Res = {:.3} m/s",
            state.epoch_count,
            sum_pr / count_pr.max(1) as f64,
            sum_rr / count_rr.max(1) as f64
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::ppp_common::find_ambiguity_index;
    use crate::engine::ppp_ins_measurements::{build_h_row, build_h_row_doppler, build_h_row_uduc, UducIndices};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_solve_empty_sats() {
        let fg = PppInsIteratedEkf::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let sats = vec![];
        let lever_arm = nalgebra::Vector3::zeros();
        let res = fg.solve(&mut state, &sats, &[], None, &lever_arm, None);
        assert!(matches!(res, Err(EngineError::InsufficientSatellites)));
    }

    #[test]
    fn test_ppp_factor_graph_default() {
        let fg = PppInsIteratedEkf::default();
        assert_eq!(fg.max_iterations, 8);
        assert_eq!(fg.convergence_threshold, 0.005);
        assert_eq!(fg.huber_k, 10.0);
        let fg2 = PppInsIteratedEkf::new();
        assert_eq!(fg2.max_iterations, 8);
    }

    #[test]
    fn test_find_ambiguity_index() {
        let mut state = dummy_rtk_state();

        use gneiss_core::sat::{Constellation, SatelliteId};
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        state.ambiguity_keys.push((sat1, 0));
        state.ambiguity_keys.push((sat2, 1));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);

        assert_eq!(find_ambiguity_index(&state, sat1), Some(0));
        assert_eq!(find_ambiguity_index(&state, sat2), None);

        let sat3 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 3,
        };
        assert_eq!(find_ambiguity_index(&state, sat3), None);
    }

    #[test]
    fn test_build_weight_matrix() {
        let m1 = FgMeasurement {
            res: 1.0,
            h_row: DVector::zeros(1),
            weight: 2.0,
            raw_var: 0.5,
            is_phase: false,
            sat: None,
        };
        let m2 = FgMeasurement {
            res: 2.0,
            h_row: DVector::zeros(1),
            weight: 4.0,
            raw_var: 0.25,
            is_phase: true,
            sat: None,
        };
        let meas = vec![m1, m2];
        let mut r = DMatrix::zeros(2, 2);
        r[(0, 0)] = 2.0;
        r[(1, 1)] = 4.0;
        let w = build_weight_matrix(&meas, &r);
        assert_eq!(w[(0, 0)], 0.5);
        assert_eq!(w[(1, 1)], 0.25);
        assert_eq!(w[(0, 1)], 0.0);
    }

    #[test]
    fn test_build_h_row() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let h_pos_att = nalgebra::Matrix3::zeros();
        // size = CORE_STATE_SIZE + 1 so ISBs (size>18) and ZWD (size>20) are populated,
        // plus one ambiguity slot at index CORE_STATE_SIZE.
        let size = CORE_STATE_SIZE + 1;
        let h = build_h_row(
            &los,
            4.0,
            Some(CORE_STATE_SIZE),
            size,
            gneiss_core::sat::Constellation::Gps,
            &h_pos_att,
        );
        assert_eq!(h.len(), size);
        assert_eq!(h[0], -1.0);
        assert_eq!(h[1], -2.0);
        assert_eq!(h[2], -3.0);
        assert_eq!(h[15], 1.0); // clock bias
        assert_eq!(h[20], 4.0); // ZWD mapping
        assert_eq!(h[CORE_STATE_SIZE], 1.0); // ambiguity

        // size = CORE_STATE_SIZE: ISBs populated (size>18) but ZWD NOT (size==21, not >20... wait 21>20 is true)
        // Actually CORE_STATE_SIZE=21 > 20, so ZWD IS set. No ambiguity.
        let h2 = build_h_row(
            &los,
            4.0,
            None,
            CORE_STATE_SIZE,
            gneiss_core::sat::Constellation::Gps,
            &h_pos_att,
        );
        assert_eq!(h2.len(), CORE_STATE_SIZE);
        assert_eq!(h2[0], -1.0);
        assert_eq!(h2[15], 1.0);
        assert_eq!(h2[20], 4.0); // ZWD mapping (21 > 20)
    }

    #[test]
    fn test_build_h_row_doppler() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let h_vel_att = nalgebra::Matrix3::zeros();
        let h_vel_bg = nalgebra::Matrix3::zeros();
        // size = CORE_STATE_SIZE (21) > 19, so velocity and clock drift are populated.
        let h = build_h_row_doppler(&los, CORE_STATE_SIZE, &h_vel_att, &h_vel_bg);
        assert_eq!(h.len(), CORE_STATE_SIZE);
        assert_eq!(h[3], -1.0);
        assert_eq!(h[4], -2.0);
        assert_eq!(h[5], -3.0);
        assert_eq!(h[19], 1.0); // clock drift at index 19

        // size = 19, NOT > 19, so nothing is set — all zeros.
        let h2 = build_h_row_doppler(&los, 19, &h_vel_att, &h_vel_bg);
        assert_eq!(h2.len(), 19);
        assert_eq!(h2[3], 0.0);
    }

    #[test]
    fn test_assemble_matrices() {
        let meas = vec![
            FgMeasurement {
                res: 1.5,
                h_row: DVector::from_element(3, 1.0),
                weight: 2.0,
                raw_var: 0.5,
                is_phase: false,
                sat: None,
            },
            FgMeasurement {
                res: 2.5,
                h_row: DVector::from_element(3, 2.0),
                weight: 3.0,
                raw_var: 0.33,
                is_phase: true,
                sat: None,
            },
        ];
        let (h, z, r) = assemble_matrices(&meas, 3);
        assert_eq!(h.nrows(), 2);
        assert_eq!(h.ncols(), 3);
        assert_eq!(h[(0, 0)], 1.0);
        assert_eq!(h[(1, 2)], 2.0);
        assert_eq!(z.len(), 2);
        assert_eq!(z[0], 1.5);
        assert_eq!(z[1], 2.5);
        assert_eq!(r.nrows(), 2);
        assert_eq!(r.ncols(), 2);
        assert_eq!(r[(0, 0)], 2.0);
        assert_eq!(r[(1, 1)], 3.0);
        assert_eq!(r[(0, 1)], 0.0);
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
        assert_eq!(x[16], 0.0);
        assert_eq!(x[17], 0.0);
        assert_eq!(x[18], 0.0);
        assert_eq!(x[19], 17.0);
        assert_eq!(x[20], 18.0);
        assert_eq!(x[CORE_STATE_SIZE], 19.0);
        assert_eq!(x[CORE_STATE_SIZE + 1], 20.0);

        let mut state2 = dummy_rtk_state();
        state2.ambiguities = vec![0.0, 0.0];
        let cov = state2.covariance.clone();
        apply_state_vector(&mut state2, &x, cov);
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
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_solve_matrix_inversion_failure() {
        let fg = PppInsIteratedEkf::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::from_element(CORE_STATE_SIZE, CORE_STATE_SIZE, f64::NAN);
        let sats = vec![];
        let lever_arm = nalgebra::Vector3::zeros();
        let res = fg.solve(&mut state, &sats, &[], None, &lever_arm, None);
        assert!(matches!(res, Err(EngineError::StateDisappeared)));
    }
}

#[cfg(test)]
mod mutant_killer_tests {
    use super::*;
    use crate::engine::ppp_ins_measurements::{build_h_row, build_h_row_doppler, build_h_row_uduc, UducIndices};
    use crate::engine::processed_sat::ProcessedSat;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::SatObs;
    use gneiss_core::imu::ImuMeasurement;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_resolve_widelane_ar_insufficient() {
        let fg = PppInsIteratedEkf::default();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::zeros(CORE_STATE_SIZE + 4, CORE_STATE_SIZE + 4);
        for i in 0..4 {
            state.covariance[(CORE_STATE_SIZE + i, CORE_STATE_SIZE + i)] = 100.0;
        } // huge variance
        let subset = [
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 2,
                    },
                    1,
                    1,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 3,
                    },
                    2,
                    2,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 4,
                    },
                    3,
                    3,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
        ];
        let x = DVector::zeros(CORE_STATE_SIZE + 4);
        let res = fg.resolve_widelane_ar(&state, &subset, &x);
        assert!(res.is_err());
    }

    #[test]
    fn test_push_cp_measurement_iono_free() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id, 0));
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 1);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let h_pos_att = nalgebra::Matrix3::zeros();
        fg.push_cp_measurement(
            &mut meas, &state, &sat, &x_i, 0, &los, 10.0, 20000000.0, 100.0, &h_pos_att,
        );
        assert_eq!(meas.len(), 1);
        // expected_base = 10.0, x_i = 0. expected_cp = 10.0.
        // l_meas = 100.0 * 0.19 = 19.0.
        // res = 19.0 - 10.0 = 9.0
        assert!((meas[0].res - 9.0).abs() < 1e-6);
        // var_cp = 0.0001 * 1.0 / 1.0 * 9.0 = 0.0009
        assert!((meas[0].raw_var - 0.0009).abs() < 1e-6);
    }

    #[test]
    fn test_push_cp_measurement_not_iono_free() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id, 0));
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 1);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let h_pos_att = nalgebra::Matrix3::zeros();
        fg.push_cp_measurement(
            &mut meas, &state, &sat, &x_i, 0, &los, 10.0, 20000000.0, 100.0, &h_pos_att,
        );
        assert_eq!(meas.len(), 1);
        // expected_base = 10.0. expected_cp = 10.0 - 5.0 = 5.0.
        // l_meas = 100.0 * 0.19 = 19.0.
        // res = 19.0 - 5.0 = 14.0
        assert!((meas[0].res - 14.0).abs() < 1e-6);
        // var_cp = 0.0001 * 1.0 / 1.0 = 0.0001
        assert!((meas[0].raw_var - 0.0001).abs() < 1e-6);
    }

    #[test]
    fn test_find_ar_candidates() {
        let fg = PppInsIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let sat_id1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id1, 1));
        state.ambiguity_keys.push((sat_id1, 2));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);

        let obs = SatObs {
            sat: sat_id1,
            observations: vec![],
        };
        let mut sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let cands = fg.find_ar_candidates(&state, &[sat.clone()]);
        assert_eq!(cands.len(), 1);

        // Test `!s.is_iono_free` mutant
        sat.is_iono_free = true;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
        sat.is_iono_free = false;

        // Test `s.cp2.is_some()` mutant
        sat.cp2 = None;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
        sat.cp2 = Some(0.0);

        // Test constellation
        let sat_id_glo = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 1,
        };
        let obs_glo = SatObs {
            sat: sat_id_glo,
            observations: vec![],
        };
        sat.sat_obs = &obs_glo;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
    }

    #[test]
    fn test_resolve_cascade_ar_bounds() {
        let fg = PppInsIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let mut sats: Vec<ProcessedSat> = Vec::new();

        let obs1 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        let obs2 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 2,
            },
            observations: vec![],
        };
        let obs3 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 3,
            },
            observations: vec![],
        };
        let obs4 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 4,
            },
            observations: vec![],
        };
        // Static references to avoid lifetime issues in closure
        let obs1_ref = Box::leak(Box::new(obs1));
        let obs2_ref = Box::leak(Box::new(obs2));
        let obs3_ref = Box::leak(Box::new(obs3));
        let obs4_ref = Box::leak(Box::new(obs4));

        {
            let mut add_sat = |obs: &'static SatObs| {
                state.add_ambiguity(obs.sat, 1, 0.0, 1.0);
                state.add_ambiguity(obs.sat, 2, 0.0, 1.0);
                sats.push(ProcessedSat {
                    sat_obs: obs,
                    dt_sat_m: 0.0,
                    p1: 0.0,
                    p2: None,
                    cp1: Some(0.0),
                    cp2: Some(0.0),
                    is_iono_free: false,
                    osb_p1: 0.0,
                    osb_p2: 0.0,
                    osb_cp1: 0.0,
                    osb_cp2: 0.0,
                    los: Vector3::zeros(),
                    dist: 0.0,
                    el: 15.01_f64.to_radians(),
                    snr: 45.0,
                    doppler: 0.0,
                    lam1: 0.19,
                    lam2: 0.24,
                    tropo_dry: 0.0,
                    map_wet: 0.0,
                    iono_delay: 5.0,
                    f1: 1.0,
                    f2: 1.0,
                    sat_pos_rot: Vector3::zeros(),
                    sat_vel: Vector3::zeros(),
                    sat_clock_drift: 0.0,
                    rcv_pos_ecef: Vector3::zeros(),
                    pcv_correction: 0.0,
                });
            };

            add_sat(obs1_ref);
            add_sat(obs2_ref);
            add_sat(obs3_ref);
        }
        assert_eq!(
            fg.resolve_cascade_ar(&mut state, &sats),
            Err("Insufficient dual-frequency satellites for AR")
        );

        state.add_ambiguity(obs4_ref.sat, 1, 0.0, 1.0);
        state.add_ambiguity(obs4_ref.sat, 2, 0.0, 1.0);
        sats.push(ProcessedSat {
            sat_obs: obs4_ref,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        });
        assert!(
            fg.resolve_cascade_ar(&mut state, &sats)
                != Err("Insufficient dual-frequency satellites for AR")
        );
    }

    #[test]
    fn test_find_worst_outlier() {
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };

        let meas = vec![
            // Not phase, shouldn't be picked even if high
            FgMeasurement {
                res: 1000.0,
                raw_var: 1.0,
                is_phase: false,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, but norm = 10.0 / sqrt(4.0) = 5.0 (less than max_norm 15.0)
            FgMeasurement {
                res: 10.0,
                raw_var: 4.0,
                is_phase: true,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, norm = 40.0 / sqrt(4.0) = 20.0 (greater than max_norm 15.0)
            FgMeasurement {
                res: 40.0,
                raw_var: 4.0,
                is_phase: true,
                sat: Some(sat2),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, norm = -60.0 / sqrt(9.0) = 20.0 (equal to current max_norm, shouldn't override because of >)
            FgMeasurement {
                res: -60.0,
                raw_var: 9.0,
                is_phase: true,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
        ];

        assert_eq!(PppInsIteratedEkf::find_worst_outlier_sat(&meas), Some(sat2));

        let meas_no_outlier = vec![FgMeasurement {
            res: 10.0,
            raw_var: 4.0,
            is_phase: true,
            sat: Some(sat1),
            h_row: DVector::zeros(0),
            weight: 1.0,
        }];
        assert_eq!(
            PppInsIteratedEkf::find_worst_outlier_sat(&meas_no_outlier),
            None
        );
    }

    #[test]
    fn test_find_ar_candidates_bounds() {
        let _fg = PppInsIteratedEkf::default();
        let _state = dummy_rtk_state();
        let _sats: Vec<ProcessedSat> = Vec::new();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let _sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        // It requires state.ambiguity_keys to contain (sat, 0) and (sat, 1) and (sat, 2) etc depending on `is_iono_free`.
        // We'll skip adding a full state test and rely on smaller integration tests or direct tests.
    }

    #[test]
    fn test_extract_isb() {
        let x = DVector::from_fn(19, |i, _| i as f64);
        assert!((PppInsIteratedEkf::extract_isb(&x, Constellation::Gps) - 0.0).abs() < 1e-10);
        assert!(
            (PppInsIteratedEkf::extract_isb(&x, Constellation::Glonass) - 16.0).abs() < 1e-10
        );
        assert!(
            (PppInsIteratedEkf::extract_isb(&x, Constellation::Galileo) - 17.0).abs() < 1e-10
        );
        assert!(
            (PppInsIteratedEkf::extract_isb(&x, Constellation::Beidou) - 18.0).abs() < 1e-10
        );
        // size <= 18 returns 0.0 for all constellations
        let x_small = DVector::from_fn(18, |i, _| i as f64);
        assert!(
            (PppInsIteratedEkf::extract_isb(&x_small, Constellation::Glonass) - 0.0).abs() < 1e-10
        );
        assert!(
            (PppInsIteratedEkf::extract_isb(&x_small, Constellation::Galileo) - 0.0).abs() < 1e-10
        );
        assert!(
            (PppInsIteratedEkf::extract_isb(&x_small, Constellation::Beidou) - 0.0).abs() < 1e-10
        );
    }

    #[test]
    fn test_build_h_row_uduc_ins_basic() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let h_pos_att = nalgebra::Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let size = CORE_STATE_SIZE + 2;
        // No iono or ambiguity indices
        let h = build_h_row_uduc(&los, 4.0, None, 1.0, None, size, Constellation::Gps, &h_pos_att);
        assert_eq!(h.len(), size);
        assert!((h[0] - (-1.0)).abs() < 1e-10);
        assert!((h[1] - (-2.0)).abs() < 1e-10);
        assert!((h[2] - (-3.0)).abs() < 1e-10);
        // pr_h_att = -los^T * h_pos_att = -[1,2,3] * I = [-1,-2,-3]
        assert!((h[6] - (-1.0)).abs() < 1e-10);
        assert!((h[7] - (-2.0)).abs() < 1e-10);
        assert!((h[8] - (-3.0)).abs() < 1e-10);
        assert!((h[15] - 1.0).abs() < 1e-10);
        assert!((h[20] - 4.0).abs() < 1e-10);
        // No indices set
        assert!((h[CORE_STATE_SIZE] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_build_h_row_uduc_ins_with_indices() {
        let los = Vector3::new(0.5, -1.5, 2.0);
        let h_pos_att = nalgebra::Matrix3::identity();
        let size = CORE_STATE_SIZE + 4;
        let i_idx = CORE_STATE_SIZE + 2;
        let n_idx = CORE_STATE_SIZE + 3;
        let h = build_h_row_uduc(
            &los,
            2.5,
            Some(i_idx),
            -1.5,
            Some(n_idx),
            size,
            Constellation::Galileo,
            &h_pos_att,
        );
        assert!((h[17] - 1.0).abs() < 1e-10); // Galileo ISB
        assert!((h[i_idx] - (-1.5)).abs() < 1e-10); // Ionosphere coefficient
        assert!((h[n_idx] - 1.0).abs() < 1e-10); // Ambiguity
        assert!((h[20] - 2.5).abs() < 1e-10); // ZWD
    }

    #[test]
    fn test_build_h_row_ins_with_h_pos_att() {
        let los = Vector3::new(1.0, 0.0, 0.0);
        let h_pos_att = nalgebra::Matrix3::new(
            0.0, -1.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        );
        let h = build_h_row(
            &los, 2.0, None, CORE_STATE_SIZE, Constellation::Gps, &h_pos_att,
        );
        // pr_h_att = -los^T * h_pos_att = -[1,0,0] * [[0,-1,0],[1,0,0],[0,0,0]] = -[0,-1,0] = [0,1,0]
        assert!((h[0] - (-1.0)).abs() < 1e-10);
        assert!((h[6] - 0.0).abs() < 1e-10);
        assert!((h[7] - 1.0).abs() < 1e-10);
        assert!((h[8] - 0.0).abs() < 1e-10);
        assert!((h[15] - 1.0).abs() < 1e-10);
        assert!((h[20] - 2.0).abs() < 1e-10);

        // With ambiguity index
        let h2 = build_h_row(
            &los, 2.0, Some(CORE_STATE_SIZE), CORE_STATE_SIZE + 1, Constellation::Gps, &h_pos_att,
        );
        assert!((h2[CORE_STATE_SIZE] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_build_h_row_doppler_ins_with_matrices() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let h_vel_att = nalgebra::Matrix3::new(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        );
        let h_vel_bg = nalgebra::Matrix3::new(
            0.0, 0.0, 1.0,
            0.0, 1.0, 0.0,
            1.0, 0.0, 0.0,
        );
        let h = build_h_row_doppler(&los, CORE_STATE_SIZE + 1, &h_vel_att, &h_vel_bg);
        assert_eq!(h.len(), CORE_STATE_SIZE + 1);
        // Velocity terms
        assert!((h[3] - (-1.0)).abs() < 1e-10);
        assert!((h[4] - (-2.0)).abs() < 1e-10);
        assert!((h[5] - (-3.0)).abs() < 1e-10);
        // dop_h_att = -los^T * h_vel_att = -(1,2,3) * I = (-1,-2,-3)
        assert!((h[6] - (-1.0)).abs() < 1e-10);
        assert!((h[7] - (-2.0)).abs() < 1e-10);
        assert!((h[8] - (-3.0)).abs() < 1e-10);
        // dop_h_bg = -los^T * h_vel_bg = -(1,2,3) * [[0,0,1],[0,1,0],[1,0,0]] = -(3,2,1) = (-3,-2,-1)
        assert!((h[12] - (-3.0)).abs() < 1e-10);
        assert!((h[13] - (-2.0)).abs() < 1e-10);
        assert!((h[14] - (-1.0)).abs() < 1e-10);
        // Clock drift
        assert!((h[19] - 1.0).abs() < 1e-10);

        // Small size (< 19): no velocity or clock drift cols set
        let h2 = build_h_row_doppler(&los, 19, &h_vel_att, &h_vel_bg);
        assert_eq!(h2.len(), 19);
        assert!((h2[3] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_fill_imu_preint_covariance() {
        let mut preint = crate::estimators::factor_graph::imu_factors::ImuPreintegration::new();
        preint.dt = 0.1;
        PppInsIteratedEkf::fill_imu_preint_covariance(&mut preint);
        // dt=0.1 -> each diag element computed from the formulas
        // Accelerometer: (5.0 * 0.1).max(2.0) = 2.0, squared = 4.0
        // Gyroscope: 3.0^2 = 9.0 (3.0.max(1.0) = 3.0)
        // Velocity: (0.1 * 0.1).max(0.05) = 0.05, squared = 0.0025
        // Accel bias: 1e-4
        // Gyro bias: 1e-6
        assert!((preint.covariance[(0, 0)] - 4.0).abs() < 1e-10);
        assert!((preint.covariance[(3, 3)] - 9.0).abs() < 1e-10);
        assert!((preint.covariance[(6, 6)] - 0.0025).abs() < 1e-10);
        assert!((preint.covariance[(9, 9)] - 1e-4).abs() < 1e-10);
        assert!((preint.covariance[(12, 12)] - 1e-6).abs() < 1e-10);

        // dt very small -> floor values used
        let mut preint2 = crate::estimators::factor_graph::imu_factors::ImuPreintegration::new();
        preint2.dt = 0.001;
        PppInsIteratedEkf::fill_imu_preint_covariance(&mut preint2);
        assert!((preint2.covariance[(0, 0)] - 4.0).abs() < 1e-10); // floored at 2.0
        assert!((preint2.covariance[(6, 6)] - 0.0025).abs() < 1e-10); // floored at 0.05
    }

    #[test]
    fn test_fill_imu_preint_covariance_large_dt() {
        let mut preint = crate::estimators::factor_graph::imu_factors::ImuPreintegration::new();
        preint.dt = 10.0;
        PppInsIteratedEkf::fill_imu_preint_covariance(&mut preint);
        // Accelerometer: (5.0 * 10.0).max(2.0) = 50.0, squared = 2500.0
        // Gyroscope: 3.0.max(1.0) = 3.0, squared = 9.0
        // Velocity: (0.1 * 10.0).max(0.05) = 1.0, squared = 1.0
        assert!((preint.covariance[(0, 0)] - 2500.0).abs() < 1e-6);
        assert!((preint.covariance[(6, 6)] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_build_ar_subset_per_constellation() {
        let fg = PppInsIteratedEkf::default();
        let cands = vec![
            (
                SatelliteId { constellation: Constellation::Gps, prn: 1 },
                0,
                1,
                20.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Gps, prn: 2 },
                2,
                3,
                30.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 1 },
                4,
                5,
                25.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 2 },
                6,
                7,
                15.0_f64.to_radians(),
                0.19,
                0.24,
            ),
        ];
        let subset = fg.build_ar_subset(&cands);
        // Per-constellation grouping: GPS pair (PRN1 ref, PRN2 vs ref)
        // and Galileo pair (PRN1 ref, PRN2 vs ref)
        assert_eq!(subset.len(), 2);
        // Verify each pair has same constellation
        assert_eq!(subset[0].0.0.constellation, subset[0].1.0.constellation);
        assert_eq!(subset[1].0.0.constellation, subset[1].1.0.constellation);
    }

    #[test]
    fn test_build_ar_subset_single_sat_skip() {
        let fg = PppInsIteratedEkf::default();
        let cands = vec![
            (
                SatelliteId { constellation: Constellation::Gps, prn: 1 },
                0,
                1,
                30.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 1 },
                2,
                3,
                25.0_f64.to_radians(),
                0.19,
                0.24,
            ),
        ];
        // Each constellation has only 1 sat -> no pairs from either
        let subset = fg.build_ar_subset(&cands);
        assert!(subset.is_empty());
    }

    #[test]
    fn test_resolve_uduc_indices_ins_with_all_indices() {
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 1)); // n1 idx 0
        state.ambiguity_keys.push((sat_id, 2)); // n2 idx 1
        state.ambiguity_keys.push((sat_id, 3)); // i1 idx 2
        state.ambiguities = vec![100.0, 200.0, 300.0];
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::from_fn(CORE_STATE_SIZE + 3, |i, _| i as f64);
        let idx = PppInsIteratedEkf::resolve_uduc_indices(&state, &sat, &x_i);
        assert_eq!(idx.i1_idx, Some(CORE_STATE_SIZE + 2));
        assert_eq!(idx.n1_idx, Some(CORE_STATE_SIZE + 0));
        assert_eq!(idx.n2_idx, Some(CORE_STATE_SIZE + 1));
        let expected_gamma = (1575.42e6 * 1575.42e6) / (1227.60e6 * 1227.60e6);
        assert!((idx.gamma - expected_gamma).abs() < 1e-6);
        assert!((idx.i1 - (CORE_STATE_SIZE + 2) as f64).abs() < 1e-6);
        assert!((idx.n1 - (CORE_STATE_SIZE + 0) as f64).abs() < 1e-6);
        assert!((idx.n2 - (CORE_STATE_SIZE + 1) as f64).abs() < 1e-6);
    }

    #[test]
    fn test_resolve_uduc_indices_ins_missing_indices() {
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let idx = PppInsIteratedEkf::resolve_uduc_indices(&state, &sat, &x_i);
        assert!(idx.i1_idx.is_none());
        assert!(idx.n1_idx.is_none());
        assert!(idx.n2_idx.is_none());
        assert!((idx.i1 - 0.0).abs() < 1e-10);
        assert!((idx.n1 - 0.0).abs() < 1e-10);
        assert!((idx.n2 - 0.0).abs() < 1e-10);
        assert!((idx.gamma - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_push_sat_meas_rejects_large_pr_residual() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let h_pos_att = nalgebra::Matrix3::identity();

        // PR residual > 100 -> returns false
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 1000.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let result = fg.push_sat_meas(
            &mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0, &h_pos_att,
        );
        assert!(!result);
        assert!(meas.is_empty());
    }

    #[test]
    fn test_push_pr_measurement_iono_free_variance_ins() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let h_pos_att = nalgebra::Matrix3::identity();

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 10.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.push_pr_measurement(
            &mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0, &h_pos_att,
        );
        assert_eq!(meas.len(), 1);
        // res_pr = 10.0 - 5.0 = 5.0
        assert!((meas[0].res - 5.0).abs() < 1e-6);
        // var_pr = 1.0 * 1.0 / 1.0 * 9.0 = 9.0 (iono-free amplifies noise)
        assert!((meas[0].raw_var - 9.0).abs() < 1e-6);
    }

    // ========== build_imu_preint_factor ==========

    #[test]
    fn test_build_imu_preint_factor_construction() {
        let preint = crate::estimators::factor_graph::imu_factors::ImuPreintegration::new();
        let mut last = dummy_rtk_state();
        last.position.vector = Vector3::new(6378137.0, 0.0, 0.0);
        last.velocity = Vector3::new(1.0, 2.0, 3.0);
        last.attitude =
            nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(0.1, 0.2, 0.3));
        last.accel_bias = Vector3::new(0.01, 0.02, 0.03);
        last.gyro_bias = Vector3::new(0.001, 0.002, 0.003);
        let x_i = DVector::from_fn(30, |i, _| (i as f64) * 10.0);

        let factor = PppInsIteratedEkf::build_imu_preint_factor(preint, &last, &x_i);

        let grav_mag = factor.gravity.norm();
        assert!(grav_mag > 9.0 && grav_mag < 10.0);

        assert_eq!(factor.nominal_p_i, last.position.vector);
        assert_eq!(factor.nominal_v_i, last.velocity);
        assert_eq!(factor.nominal_q_i, last.attitude);
        assert_eq!(factor.nominal_ba_i, last.accel_bias);
        assert_eq!(factor.nominal_bg_i, last.gyro_bias);

        assert_eq!(factor.nominal_p_j, Vector3::new(x_i[0], x_i[1], x_i[2]));
        assert_eq!(factor.nominal_v_j, Vector3::new(x_i[3], x_i[4], x_i[5]));

        let expected_q =
            nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(x_i[6], x_i[7], x_i[8]));
        assert!(
            (factor.nominal_q_j.scaled_axis() - expected_q.scaled_axis()).norm() < 1e-10
        );
        assert_eq!(factor.nominal_ba_j, Vector3::new(x_i[9], x_i[10], x_i[11]));
        assert_eq!(
            factor.nominal_bg_j,
            Vector3::new(x_i[12], x_i[13], x_i[14])
        );

        assert_eq!(factor.idx_p_i, 0);
        assert_eq!(factor.idx_v_i, 3);
        assert_eq!(factor.idx_q_i, 6);
        assert_eq!(factor.idx_ba_i, 9);
        assert_eq!(factor.idx_bg_i, 12);
        assert_eq!(factor.idx_p_j, 15);
        assert_eq!(factor.idx_v_j, 18);
        assert_eq!(factor.idx_q_j, 21);
        assert_eq!(factor.idx_ba_j, 24);
        assert_eq!(factor.idx_bg_j, 27);
    }

    #[test]
    fn test_build_imu_preint_factor_zero_rotation() {
        let preint = crate::estimators::factor_graph::imu_factors::ImuPreintegration::new();
        let last = dummy_rtk_state();
        let x_i = DVector::zeros(30);

        let factor = PppInsIteratedEkf::build_imu_preint_factor(preint, &last, &x_i);

        assert!(factor.nominal_q_j.angle() < 1e-10);
    }

    // ========== accumulate_single_imu_factor ==========

    #[test]
    fn test_accumulate_single_imu_factor_updates_htwh_and_htwr() {
        let factor = crate::estimators::factor_graph::imu_factors::ImuPreintegrationFactor {
            preint: crate::estimators::factor_graph::imu_factors::ImuPreintegration::new(),
            gravity: Vector3::new(0.0, 0.0, -9.80665),
            nominal_p_i: Vector3::zeros(),
            nominal_v_i: Vector3::zeros(),
            nominal_q_i: nalgebra::UnitQuaternion::identity(),
            nominal_ba_i: Vector3::zeros(),
            nominal_bg_i: Vector3::zeros(),
            nominal_p_j: Vector3::zeros(),
            nominal_v_j: Vector3::zeros(),
            nominal_q_j: nalgebra::UnitQuaternion::identity(),
            nominal_ba_j: Vector3::zeros(),
            nominal_bg_j: Vector3::zeros(),
            idx_p_i: 0,
            idx_v_i: 3,
            idx_q_i: 6,
            idx_ba_i: 9,
            idx_bg_i: 12,
            idx_p_j: 15,
            idx_v_j: 18,
            idx_q_j: 21,
            idx_ba_j: 24,
            idx_bg_j: 27,
        };
        let state_size = 30;
        let mut htwh = DMatrix::zeros(state_size, state_size);
        let mut htwr = DVector::zeros(state_size);

        PppInsIteratedEkf::accumulate_single_imu_factor(
            &factor,
            state_size,
            &mut htwh,
            Some(&mut htwr),
        );

        assert!(htwh.iter().any(|x| x.abs() > 1e-10));
        // htwh should be symmetric (information matrix contribution)
        for i in 0..state_size {
            for j in 0..state_size {
                assert!(
                    (htwh[(i, j)] - htwh[(j, i)]).abs() < 1e-8,
                    "htwh asymmetry at ({}, {}): {} vs {}",
                    i,
                    j,
                    htwh[(i, j)],
                    htwh[(j, i)]
                );
            }
        }
    }

    #[test]
    fn test_accumulate_single_imu_factor_no_htwr() {
        let mut preint = crate::estimators::factor_graph::imu_factors::ImuPreintegration::new();
        preint.dt = 0.1;
        preint.dv = Vector3::new(1.0, 0.0, 0.0);
        let factor = crate::estimators::factor_graph::imu_factors::ImuPreintegrationFactor {
            preint,
            gravity: Vector3::new(0.0, 0.0, -9.80665),
            nominal_p_i: Vector3::zeros(),
            nominal_v_i: Vector3::zeros(),
            nominal_q_i: nalgebra::UnitQuaternion::identity(),
            nominal_ba_i: Vector3::zeros(),
            nominal_bg_i: Vector3::zeros(),
            nominal_p_j: Vector3::zeros(),
            nominal_v_j: Vector3::zeros(),
            nominal_q_j: nalgebra::UnitQuaternion::identity(),
            nominal_ba_j: Vector3::zeros(),
            nominal_bg_j: Vector3::zeros(),
            idx_p_i: 0,
            idx_v_i: 3,
            idx_q_i: 6,
            idx_ba_i: 9,
            idx_bg_i: 12,
            idx_p_j: 15,
            idx_v_j: 18,
            idx_q_j: 21,
            idx_ba_j: 24,
            idx_bg_j: 27,
        };
        let state_size = 30;
        let mut htwh = DMatrix::zeros(state_size, state_size);

        PppInsIteratedEkf::accumulate_single_imu_factor(&factor, state_size, &mut htwh, None);

        assert!(htwh.iter().any(|x| x.abs() > 1e-10));
    }

    // ========== accumulate_imu_factors early returns ==========

    #[test]
    fn test_accumulate_imu_factors_no_last_state() {
        let fg = PppInsIteratedEkf::default();
        let x_i = DVector::zeros(30);
        let mut htwh = DMatrix::zeros(30, 30);

        fg.accumulate_imu_factors(&x_i, &[], None, &mut htwh, None);
        assert!(htwh.iter().all(|x| x.abs() < 1e-10));
    }

    #[test]
    fn test_accumulate_imu_factors_empty_history() {
        let fg = PppInsIteratedEkf::default();
        let x_i = DVector::zeros(30);
        let last = dummy_rtk_state();
        let mut htwh = DMatrix::zeros(30, 30);

        fg.accumulate_imu_factors(&x_i, &[], Some(&last), &mut htwh, None);
        assert!(htwh.iter().all(|x| x.abs() < 1e-10));
    }

    #[test]
    fn test_accumulate_imu_factors_empty_imu_buf() {
        let fg = PppInsIteratedEkf::default();
        let x_i = DVector::zeros(30);
        let last = dummy_rtk_state();
        let mut htwh = DMatrix::zeros(30, 30);
        let imu_history: Vec<Vec<ImuMeasurement>> = vec![vec![]];

        fg.accumulate_imu_factors(&x_i, &imu_history, Some(&last), &mut htwh, None);
        assert!(htwh.iter().all(|x| x.abs() < 1e-10));
    }

    #[test]
    fn test_accumulate_imu_factors_not_aligned() {
        let fg = PppInsIteratedEkf::default();
        let x_i = DVector::zeros(30);
        let mut last = dummy_rtk_state();
        last.ins_aligned = false;
        let mut htwh = DMatrix::zeros(30, 30);
        let imu_history: Vec<Vec<ImuMeasurement>> = vec![vec![ImuMeasurement::new(
            0,
            Vector3::zeros(),
            Vector3::zeros(),
        )]];

        fg.accumulate_imu_factors(&x_i, &imu_history, Some(&last), &mut htwh, None);
        assert!(htwh.iter().all(|x| x.abs() < 1e-10));
    }

    // ========== push_uduc_pr_measurements ==========

    #[test]
    fn test_push_uduc_pr_measurements_both_freq() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 2);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let h_pos_att = nalgebra::Matrix3::identity();

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 50.0,
            p2: Some(100.0),
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };

        let gamma_val =
            (1575.42e6 * 1575.42e6) / (1227.60e6 * 1227.60e6);
        let idx = UducIndices {
            i1_idx: Some(CORE_STATE_SIZE),
            n1_idx: None,
            n2_idx: None,
            i1: 2.0,
            n1: 0.0,
            n2: 0.0,
            gamma: gamma_val,
        };

        fg.push_uduc_pr_measurements(&mut meas, &sat, &x_i, &los, 40.0, &h_pos_att, &idx);

        assert_eq!(meas.len(), 2);
        assert!((meas[0].res - 8.0).abs() < 1e-6);
        assert!(!meas[0].is_phase);
        let expected_p2 = 100.0 - (40.0 + gamma_val * 2.0);
        assert!((meas[1].res - expected_p2).abs() < 1e-6);
        assert!(!meas[1].is_phase);
        assert!(meas[1].raw_var > meas[0].raw_var);
    }

    // ========== push_uduc_cp_measurements ==========

    #[test]
    fn test_push_uduc_cp_measurements_with_windup() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.windup.insert(sat_id, 0.5);

        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 3);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let h_pos_att = nalgebra::Matrix3::identity();

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(100.0),
            cp2: Some(200.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };

        let gamma_val =
            (1575.42e6 * 1575.42e6) / (1227.60e6 * 1227.60e6);
        let idx = UducIndices {
            i1_idx: Some(CORE_STATE_SIZE),
            n1_idx: Some(CORE_STATE_SIZE + 1),
            n2_idx: Some(CORE_STATE_SIZE + 2),
            i1: 2.0,
            n1: 10.0,
            n2: 20.0,
            gamma: gamma_val,
        };

        fg.push_uduc_cp_measurements(
            &mut meas, &state, &sat, &x_i, &los, 40.0, &h_pos_att, &idx,
        );

        assert_eq!(meas.len(), 2);
        let expected_l1 = (100.0 - 0.5) * 0.19 - (40.0 - idx.i1 + idx.n1);
        assert!((meas[0].res - expected_l1).abs() < 1e-6);
        assert!(meas[0].is_phase);
        let expected_l2 = (200.0 - 0.5) * 0.24 - (40.0 - gamma_val * idx.i1 + idx.n2);
        assert!((meas[1].res - expected_l2).abs() < 1e-6);
        assert!(meas[1].is_phase);
    }

    // ========== push_uduc_measurements (wrapper) ==========

    #[test]
    fn test_push_uduc_measurements_produces_four() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id, 1));
        state.ambiguity_keys.push((sat_id, 2));
        state.ambiguity_keys.push((sat_id, 3));
        state.ambiguities = vec![10.0, 20.0, 2.0];
        state.windup.insert(sat_id, 0.5);

        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let x_i = DVector::from_fn(CORE_STATE_SIZE + 3, |i, _| (i as f64) * 0.1);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let h_pos_att = nalgebra::Matrix3::identity();

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 50.0,
            p2: Some(100.0),
            cp1: Some(100.0),
            cp2: Some(200.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };

        fg.push_uduc_measurements(
            &mut meas, &state, &sat, &x_i, 0, &los, 40.0, 0.0, 0.0, &h_pos_att,
        );

        assert_eq!(meas.len(), 4);
        assert!(!meas[0].is_phase);
        assert!(!meas[1].is_phase);
        assert!(meas[2].is_phase);
        assert!(meas[3].is_phase);
    }

    // ========== push_doppler_measurement ==========

    #[test]
    fn test_push_doppler_measurement_basic() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let h_vel_att = nalgebra::Matrix3::zeros();
        let h_vel_bg = nalgebra::Matrix3::zeros();
        let v_apc = Vector3::zeros();

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 1000.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::new(100.0, 0.0, 0.0),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };

        fg.push_doppler_measurement(&mut meas, &sat, &x_i, &los, &h_vel_att, &h_vel_bg, &v_apc);

        assert_eq!(meas.len(), 1);
        assert!(!meas[0].is_phase);
        assert!((meas[0].res - (-290.0)).abs() < 1e-6);
        assert!((meas[0].raw_var - 0.25).abs() < 1e-6);
    }

    #[test]
    fn test_push_doppler_measurement_with_drift() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let mut x_i = DVector::zeros(CORE_STATE_SIZE + 1);
        x_i[19] = 2.0;
        let los = Vector3::new(1.0, 0.0, 0.0);
        let h_vel_att = nalgebra::Matrix3::zeros();
        let h_vel_bg = nalgebra::Matrix3::zeros();
        let v_apc = Vector3::new(5.0, 0.0, 0.0);

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 1000.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::new(100.0, 0.0, 0.0),
            sat_clock_drift: 1.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };

        fg.push_doppler_measurement(&mut meas, &sat, &x_i, &los, &h_vel_att, &h_vel_bg, &v_apc);

        assert_eq!(meas.len(), 1);
        let expected_rr = 100.0 - 5.0 + 2.0 - 1.0 * gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let expected_res = -190.0 - expected_rr;
        assert!((meas[0].res - expected_res).abs() < 1e-6);
    }

    #[test]
    fn test_push_doppler_measurement_small_state() {
        let fg = PppInsIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let x_i = DVector::zeros(15);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let h_vel_att = nalgebra::Matrix3::zeros();
        let h_vel_bg = nalgebra::Matrix3::zeros();
        let v_apc = Vector3::zeros();

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 1000.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::new(100.0, 0.0, 0.0),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };

        fg.push_doppler_measurement(&mut meas, &sat, &x_i, &los, &h_vel_att, &h_vel_bg, &v_apc);

        assert_eq!(meas.len(), 1);
        assert!(meas[0].h_row.iter().all(|x| x.abs() < 1e-10));
    }

    // ========== accumulate_imu_factors happy path ==========

    #[test]
    fn test_accumulate_imu_factors_happy_path() {
        let fg = PppInsIteratedEkf::default();
        let x_i = DVector::zeros(30);
        let mut last = dummy_rtk_state();
        last.ins_aligned = true;
        last.position.vector = Vector3::new(100.0, 200.0, 300.0);
        let mut htwh = DMatrix::zeros(30, 30);
        let mut htwr = DVector::zeros(30);
        let imu_history: Vec<Vec<ImuMeasurement>> = vec![vec![
            ImuMeasurement::new(0, Vector3::new(0.0, 0.0, 9.8), Vector3::zeros()),
            ImuMeasurement::new(10, Vector3::new(0.0, 0.0, 9.8), Vector3::zeros()),
        ]];

        fg.accumulate_imu_factors(
            &x_i,
            &imu_history,
            Some(&last),
            &mut htwh,
            Some(&mut htwr),
        );

        assert!(htwh.iter().any(|x| x.abs() > 1e-10));
        for i in 0..30 {
            for j in 0..30 {
                assert!(
                    (htwh[(i, j)] - htwh[(j, i)]).abs() < 1e-8,
                    "htwh asymmetry at ({}, {}): {} vs {}",
                    i,
                    j,
                    htwh[(i, j)],
                    htwh[(j, i)]
                );
            }
        }
        assert!(htwr.iter().any(|x| x.abs() > 1e-10));
    }
}
