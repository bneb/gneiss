use crate::engine::ppp_common::{
    build_iono_constraint_row, find_amb_idx, find_ambiguity_index, snr_scale, FgMeasurement,
};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::types::IonosphereModel;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::math::thresholding::apply_huber;
use nalgebra::{DMatrix, DVector, Vector3};

pub(crate) const SPEED_OF_LIGHT: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
pub(crate) const PSEUDORANGE_VARIANCE_BASE: f64 = 1.0;

impl crate::engine::ppp_iekf::PppIteratedEkf {
    pub(crate) fn build_measurements(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        iter: usize,
    ) -> Vec<FgMeasurement> {
        let mut meas = Vec::new();
        let tide_offset =
            gneiss_core::tides::solid_earth_tides_ecef(state.time, state.position.vector);
        let rcv_pos = Vector3::new(x_i[0], x_i[1], x_i[2]) + tide_offset;
        let ztd = if x_i.len() > 20 && !x_i[20].is_nan() && x_i[20] != 0.0 {
            x_i[20]
        } else {
            state.zwd
        };

        for sat in sats {
            let geometric_dist = (sat.sat_pos_rot - rcv_pos).norm();
            let dist = geometric_dist - sat.pcv_correction;
            let los = (sat.sat_pos_rot - rcv_pos) / geometric_dist;
            let isb = Self::extract_isb(x_i, sat.sat_obs.sat.constellation);
            let expected_base =
                dist + x_i[15] + isb - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet;

            if !self.push_sat_meas(
                &mut meas,
                state,
                sat,
                x_i,
                iter,
                &los,
                expected_base,
                dist,
                isb,
            ) {
                continue;
            }
            if sat.doppler != 0.0 {
                self.push_doppler_measurement(&mut meas, sat, x_i, &los);
            }
        }
        meas
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_sat_meas(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        dist: f64,
        isb: f64,
    ) -> bool {
        let expected_pr = if sat.is_iono_free {
            expected_base
        } else if sat.cp1.is_some() && sat.cp2.is_some() && sat.p2.is_some() {
            expected_base
                + find_amb_idx(state, sat.sat_obs.sat, 3)
                    .map(|i| x_i[CORE_STATE_SIZE + i])
                    .unwrap_or(sat.iono_delay)
        } else {
            expected_base + sat.iono_delay
        };

        let res_pr = sat.p1 - expected_pr;
        // Adaptive rejection threshold: during convergence (position std > 20m),
        // allow residuals up to 500m. When converged (std < 2m), tighten to 100m.
        // Formula: threshold = max(100, 5 * position_std), capped at 500m.
        // This prevents the filter death spiral where early SPP position errors
        // (common at equatorial stations) cause all measurements to be rejected,
        // leading to permanent coasting and divergence.
        let pos_var = state.covariance[(0, 0)]
            .max(state.covariance[(1, 1)])
            .max(state.covariance[(2, 2)]);
        let pos_std = pos_var.sqrt();
        let pr_threshold = (100.0_f64).max(5.0 * pos_std).min(500.0);
        if res_pr.abs() > pr_threshold {
            tracing::debug!(
                "PR rejected: sat={}, res_pr={:.1}m, threshold={:.1}m, pos_std={:.1}m",
                sat.sat_obs.sat, res_pr, pr_threshold, pos_std
            );
            return false;
        }

        if !sat.is_iono_free && sat.cp1.is_some() && sat.cp2.is_some() && sat.p2.is_some() {
            self.push_uduc_measurements(meas, state, sat, x_i, iter, los, expected_base, dist, isb);
            // Add ionospheric prior constraint: tie i1 state to Klobuchar/Ionex prediction
            if let Some(i1_idx) = find_amb_idx(state, sat.sat_obs.sat, 3) {
                let i1_est = x_i.get(CORE_STATE_SIZE + i1_idx).copied().unwrap_or(0.0);
                let res_i1 = sat.iono_delay - i1_est;
                let var_i1 = match self.iono_model {
                    IonosphereModel::Klobuchar => 9.0,    // 3m std
                    IonosphereModel::Ionex => 0.0025,     // 0.05m std (5cm)
                };
                meas.push(FgMeasurement {
                    res: res_i1,
                    h_row: build_iono_constraint_row(x_i.len(), CORE_STATE_SIZE + i1_idx),
                    weight: var_i1,
                    raw_var: var_i1,
                    is_phase: false,
                    sat: Some(sat.sat_obs.sat),
                });
            } else {
                tracing::debug!(
                    "UDUC iono prior skipped: band-3 ambiguity missing for {}",
                    sat.sat_obs.sat
                );
            }
        } else {
            self.push_pr_measurement(meas, state, sat, x_i, iter, los, expected_base, dist, isb);
            self.try_push_cp_measurement(
                meas,
                state,
                sat,
                x_i,
                iter,
                los,
                expected_base,
                dist,
            );
        }
        true
    }

    pub(crate) fn extract_isb(x_i: &DVector<f64>, constel: gneiss_core::sat::Constellation) -> f64 {
        if x_i.len() > 18 {
            match constel {
                gneiss_core::sat::Constellation::Glonass => x_i[16],
                gneiss_core::sat::Constellation::Galileo => x_i[17],
                gneiss_core::sat::Constellation::Beidou => x_i[18],
                _ => 0.0,
            }
        } else {
            0.0
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_pr_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        _dist: f64,
        _isb: f64,
    ) {
        let expected_pr = if sat.is_iono_free {
            expected_base
        } else {
            expected_base + sat.iono_delay
        };
        let res_pr = sat.p1 - expected_pr;

        if state.epoch_count == 0 && iter == 0 {
            tracing::trace!(
                "PPP {:?}{:02} PR res={:.3}m",
                sat.sat_obs.sat.constellation,
                sat.sat_obs.sat.prn,
                res_pr
            );
        }
        let mut var_pr = PSEUDORANGE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        if sat.is_iono_free {
            var_pr *= 9.0; // Iono-free combination amplifies noise
        } else {
            var_pr += 9.0; // Single frequency has ~3m Klobuchar residual iono error (3^2 = 9)
        }
        let w_pr = apply_huber(res_pr, var_pr, self.huber_k);
        meas.push(FgMeasurement {
            res: res_pr,
            h_row: build_h_row(
                los,
                sat.map_wet,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: var_pr / w_pr,
            raw_var: var_pr,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    pub(crate) fn push_doppler_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
    ) {
        let rcv_vel = Vector3::new(x_i[3], x_i[4], x_i[5]);
        let rcv_clk_drift = if x_i.len() > 19 { x_i[19] } else { 0.0 };
        let meas_rr = -sat.doppler * sat.lam1;
        let expected_rr = los.dot(&sat.sat_vel) - los.dot(&rcv_vel) + rcv_clk_drift
            - sat.sat_clock_drift * SPEED_OF_LIGHT;

        let res_rr = meas_rr - expected_rr;
        tracing::debug!("DOPPLER {}: res_rr={:.3} meas_rr={:.3} exp_rr={:.3} doppler={:.3} rcv_drift={:.3} sat_drift={:.3} los_v={:.3} sat_vel=[{:.3}, {:.3}, {:.3}]",
            sat.sat_obs.sat.to_string(), res_rr, meas_rr, expected_rr, sat.doppler, rcv_clk_drift, sat.sat_clock_drift * gneiss_core::constants::SPEED_OF_LIGHT_M_S, los.dot(&sat.sat_vel), sat.sat_vel.x, sat.sat_vel.y, sat.sat_vel.z);
        let var_rr = 0.01; // Decreased Doppler variance (trust velocity more)
        let w_rr = apply_huber(res_rr, var_rr, 3.0);
        meas.push(FgMeasurement {
            res: res_rr,
            h_row: build_h_row_doppler(los, x_i.len()),
            weight: var_rr / w_rr,
            raw_var: var_rr,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_cp_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        dist: f64,
        cp1: f64,
    ) {
        if let Some(amb_idx) = find_ambiguity_index(state, sat.sat_obs.sat) {
            let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
            let l_meas = (cp1 - windup) * sat.lam1;
            let n_amb = x_i[CORE_STATE_SIZE + amb_idx];
            let expected_cp = if sat.is_iono_free {
                expected_base + n_amb
            } else {
                expected_base - sat.iono_delay + n_amb
            };
            let res_cp = l_meas - expected_cp;
            if res_cp.abs() > 100.0 && iter == 0 {
                tracing::warn!("HUGE res_cp: sat={}, l_meas={:.2}, exp={:.2}, dist={:.2}, clk={:.2}, n_amb={:.2}", sat.sat_obs.sat, l_meas, expected_cp, dist, x_i[15], n_amb);
            }
            let mut var_cp = 0.0001 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
            if sat.is_iono_free {
                var_cp *= 9.0;
            } // Iono-free amplifies phase noise
            let w_cp = apply_huber(res_cp, var_cp, self.huber_k);
            meas.push(FgMeasurement {
                res: res_cp,
                h_row: build_h_row(
                    los,
                    sat.map_wet,
                    Some(CORE_STATE_SIZE + amb_idx),
                    x_i.len(),
                    sat.sat_obs.sat.constellation,
                ),
                weight: var_cp / w_cp,
                raw_var: var_cp,
                is_phase: true,
                sat: Some(sat.sat_obs.sat),
            });
        }
    }

    /// Push CP measurement only if sat.cp1 is Some and non-zero.
    /// Extracted as a helper to reduce nesting depth in push_sat_meas.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn try_push_cp_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        dist: f64,
    ) {
        let cp1 = match sat.cp1 {
            Some(cp) if cp != 0.0 => cp,
            _ => return,
        };
        self.push_cp_measurement(meas, state, sat, x_i, iter, los, expected_base, dist, cp1);
    }

    pub(crate) fn resolve_uduc_indices(
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
    ) -> (
        Option<usize>,
        Option<usize>,
        Option<usize>,
        f64,
        f64,
        f64,
        f64,
    ) {
        let i1_idx = find_amb_idx(state, sat.sat_obs.sat, 3).map(|idx| CORE_STATE_SIZE + idx);
        let n1_idx = find_amb_idx(state, sat.sat_obs.sat, 1).map(|idx| CORE_STATE_SIZE + idx);
        let n2_idx = find_amb_idx(state, sat.sat_obs.sat, 2).map(|idx| CORE_STATE_SIZE + idx);
        let i1 = i1_idx.map(|idx| x_i[idx]).unwrap_or(0.0);
        let n1 = n1_idx.map(|idx| x_i[idx]).unwrap_or(0.0);
        let n2 = n2_idx.map(|idx| x_i[idx]).unwrap_or(0.0);
        (
            i1_idx,
            n1_idx,
            n2_idx,
            i1,
            n1,
            n2,
            (sat.f1 * sat.f1) / (sat.f2 * sat.f2),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_uduc_pr_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
        expected_base: f64,
        i1_idx: Option<usize>,
        i1: f64,
        gamma: f64,
    ) {
        let var_p1 = PSEUDORANGE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        let res_p1 = sat.p1 - (expected_base + i1);
        meas.push(FgMeasurement {
            res: res_p1,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                1.0,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: var_p1 / apply_huber(res_p1, var_p1, self.huber_k),
            raw_var: var_p1,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
        let res_p2 = match sat.p2 {
            Some(p2) => p2 - (expected_base + gamma * i1),
            None => return,
        };
        meas.push(FgMeasurement {
            res: res_p2,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                gamma,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: (var_p1 * 1.5) / apply_huber(res_p2, var_p1 * 1.5, self.huber_k),
            raw_var: var_p1 * 1.5,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_uduc_cp_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
        expected_base: f64,
        i1_idx: Option<usize>,
        n1_idx: Option<usize>,
        n2_idx: Option<usize>,
        i1: f64,
        n1: f64,
        n2: f64,
        gamma: f64,
    ) {
        let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
        let var_l1 = 0.0001 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        let cp1_val = match sat.cp1 {
            Some(cp1) => cp1,
            None => return,
        };
        let res_l1 = (cp1_val - windup) * sat.lam1 - (expected_base - i1 + n1);
        meas.push(FgMeasurement {
            res: res_l1,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                -1.0,
                n1_idx,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: var_l1 / apply_huber(res_l1, var_l1, self.huber_k),
            raw_var: var_l1,
            is_phase: true,
            sat: Some(sat.sat_obs.sat),
        });
        let cp2_val = match sat.cp2 {
            Some(cp2) => cp2,
            None => return,
        };
        let res_l2 = (cp2_val - windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
        meas.push(FgMeasurement {
            res: res_l2,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                -gamma,
                n2_idx,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: (var_l1 * 1.5) / apply_huber(res_l2, var_l1 * 1.5, self.huber_k),
            raw_var: var_l1 * 1.5,
            is_phase: true,
            sat: Some(sat.sat_obs.sat),
        });
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_uduc_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        _iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        _dist: f64,
        _isb: f64,
    ) {
        let (i1_idx, n1_idx, n2_idx, i1, n1, n2, gamma) =
            Self::resolve_uduc_indices(state, sat, x_i);
        self.push_uduc_pr_measurements(meas, sat, x_i, los, expected_base, i1_idx, i1, gamma);
        self.push_uduc_cp_measurements(
            meas,
            state,
            sat,
            x_i,
            los,
            expected_base,
            i1_idx,
            n1_idx,
            n2_idx,
            i1,
            n1,
            n2,
            gamma,
        );
    }
}

pub(crate) fn log_ppp_convergence(
    state: &RtkState,
    sats: &[ProcessedSat],
    x_i: &DVector<f64>,
    x_pred: &DVector<f64>,
    p_pred: &DMatrix<f64>,
    solver: &crate::engine::ppp_iekf::PppIteratedEkf,
) {
    let _p_amb = if p_pred.nrows() > 21 {
        p_pred[(21, 21)]
    } else {
        0.0
    };
    tracing::info!(
        "PPP Epoch: pos=[{:.2}, {:.2}, {:.2}], dx_norm={:.4}",
        x_i[0],
        x_i[1],
        x_i[2],
        (x_i.clone() - x_pred.clone()).norm()
    );

    let meas = solver.build_measurements(state, sats, x_i, solver.max_iterations - 1);
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

pub(crate) fn build_h_row_uduc(
    los: &Vector3<f64>,
    map_wet: f64,
    i_idx: Option<usize>,
    i_coef: f64,
    n_idx: Option<usize>,
    size: usize,
    constel: gneiss_core::sat::Constellation,
) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[0] = -los.x;
    h[1] = -los.y;
    h[2] = -los.z;
    h[15] = 1.0;
    if size > 18 {
        match constel {
            gneiss_core::sat::Constellation::Glonass => h[16] = 1.0,
            gneiss_core::sat::Constellation::Galileo => h[17] = 1.0,
            gneiss_core::sat::Constellation::Beidou => h[18] = 1.0,
            _ => {}
        }
    }
    if size > 20 {
        h[20] = map_wet;
    }
    if let Some(idx) = i_idx {
        h[idx] = i_coef;
    }
    if let Some(idx) = n_idx {
        h[idx] = 1.0;
    }
    h
}

pub(crate) fn build_h_row(
    los: &Vector3<f64>,
    map_wet: f64,
    amb_idx: Option<usize>,
    size: usize,
    constel: gneiss_core::sat::Constellation,
) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[0] = -los.x;
    h[1] = -los.y;
    h[2] = -los.z;
    h[15] = 1.0;
    if size > 18 {
        match constel {
            gneiss_core::sat::Constellation::Glonass => h[16] = 1.0,
            gneiss_core::sat::Constellation::Galileo => h[17] = 1.0,
            gneiss_core::sat::Constellation::Beidou => h[18] = 1.0,
            _ => {}
        }
    }
    if size > 20 {
        h[20] = map_wet;
    }
    if let Some(idx) = amb_idx {
        h[idx] = 1.0;
    }
    h
}

pub(crate) fn build_h_row_doppler(los: &Vector3<f64>, size: usize) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    if size > 19 {
        h[3] = -los.x;
        h[4] = -los.y;
        h[5] = -los.z;
        h[19] = 1.0;
    }
    h
}
