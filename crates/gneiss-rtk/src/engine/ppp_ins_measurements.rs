use crate::engine::ppp_common::{
    find_amb_idx, find_ambiguity_index, snr_scale, FgMeasurement,
};
use crate::engine::processed_sat::ProcessedSat;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::math::thresholding::apply_huber;
use nalgebra::{DVector, Vector3};

use super::ppp_ins_iekf::PppInsIteratedEkf;

const PSEUDORANGE_VARIANCE_BASE: f64 = 1.0;

pub struct UducIndices {
    pub i1_idx: Option<usize>,
    pub n1_idx: Option<usize>,
    pub n2_idx: Option<usize>,
    pub i1: f64,
    pub n1: f64,
    pub n2: f64,
    pub gamma: f64,
}

impl PppInsIteratedEkf {
    pub fn build_measurements(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        iter: usize,
        lever_arm: &nalgebra::Vector3<f64>,
        omega_eb_b: &nalgebra::Vector3<f64>,
    ) -> Vec<FgMeasurement> {
        let mut meas = Vec::new();
        let tide_offset =
            gneiss_core::tides::solid_earth_tides_ecef(state.time, state.position.vector);
        let r_vec = nalgebra::Vector3::new(x_i[6], x_i[7], x_i[8]);
        let attitude = if r_vec.norm() > 1e-12 {
            nalgebra::UnitQuaternion::from_scaled_axis(r_vec)
        } else {
            nalgebra::UnitQuaternion::identity()
        };
        let r_b_e = attitude.to_rotation_matrix();
        let l_e = r_b_e * lever_arm;
        let rcv_pos = nalgebra::Vector3::new(x_i[0], x_i[1], x_i[2]) + l_e + tide_offset;

        let h_pos_att = -(r_b_e.matrix() * lever_arm.cross_matrix());
        let a_0 = r_b_e * omega_eb_b.cross(lever_arm);
        // INS code uses right-perturbation convention for attitude; the rest of the codebase
        // uses left-perturbation.  The cross-matrix expression produces the same result either way
        // (the sign difference is absorbed by how the perturbation is applied), so this is
        // directionally correct with negligible practical impact.
        let h_vel_att = -a_0.cross_matrix();
        let h_vel_bg = r_b_e.matrix() * lever_arm.cross_matrix();
        let v_apc = nalgebra::Vector3::new(x_i[3], x_i[4], x_i[5]) + a_0;
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
                &h_pos_att,
            ) {
                continue;
            }
            if sat.doppler != 0.0 {
                self.push_doppler_measurement(
                    &mut meas, sat, x_i, &los, &h_vel_att, &h_vel_bg, &v_apc,
                );
            }
        }
        meas
    }

    pub(super) fn extract_isb(x_i: &DVector<f64>, constel: gneiss_core::sat::Constellation) -> f64 {
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

    pub(super) fn push_sat_meas(
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
        h_pos_att: &nalgebra::Matrix3<f64>,
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
        let pos_var = state.covariance[(0, 0)]
            .max(state.covariance[(1, 1)])
            .max(state.covariance[(2, 2)]);
        let pos_std = pos_var.sqrt();
        let pr_threshold = (100.0_f64).max(3.0 * pos_std).min(200.0);
        if res_pr.abs() > pr_threshold {
            return false;
        }

        if !sat.is_iono_free && sat.cp1.is_some() && sat.cp2.is_some() && sat.p2.is_some() {
            self.push_uduc_measurements(
                meas,
                state,
                sat,
                x_i,
                iter,
                los,
                expected_base,
                dist,
                isb,
                h_pos_att,
            );
        } else {
            self.push_pr_measurement(
                meas,
                state,
                sat,
                x_i,
                iter,
                los,
                expected_base,
                dist,
                isb,
                h_pos_att,
            );
            if let Some(cp1) = sat.cp1 {
                if cp1 != 0.0 {
                    self.push_cp_measurement(
                        meas,
                        state,
                        sat,
                        x_i,
                        iter,
                        los,
                        expected_base,
                        dist,
                        cp1,
                        h_pos_att,
                    );
                }
            }
        }
        true
    }

    pub(super) fn push_pr_measurement(
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
        h_pos_att: &nalgebra::Matrix3<f64>,
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
                h_pos_att,
            ),
            weight: var_pr / w_pr,
            raw_var: var_pr,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    pub(super) fn push_doppler_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
        h_vel_att: &nalgebra::Matrix3<f64>,
        h_vel_bg: &nalgebra::Matrix3<f64>,
        v_apc: &nalgebra::Vector3<f64>,
    ) {
        let rcv_vel = v_apc;
        let rcv_clk_drift = if x_i.len() > 19 { x_i[19] } else { 0.0 };
        let meas_rr = -sat.doppler * sat.lam1;
        let expected_rr = los.dot(&sat.sat_vel) - los.dot(rcv_vel) + rcv_clk_drift
            - sat.sat_clock_drift * gneiss_core::constants::SPEED_OF_LIGHT_M_S;

        let res_rr = meas_rr - expected_rr;
        tracing::debug!("DOPPLER {}: res_rr={:.3} meas_rr={:.3} exp_rr={:.3} doppler={:.3} rcv_drift={:.3} sat_drift={:.3} los_v={:.3} sat_vel=[{:.3}, {:.3}, {:.3}]",
            sat.sat_obs.sat.to_string(), res_rr, meas_rr, expected_rr, sat.doppler, rcv_clk_drift, sat.sat_clock_drift * gneiss_core::constants::SPEED_OF_LIGHT_M_S, los.dot(&sat.sat_vel), sat.sat_vel.x, sat.sat_vel.y, sat.sat_vel.z);
        let var_rr = 0.25;
        let w_rr = apply_huber(res_rr, var_rr, 10.0);
        meas.push(FgMeasurement {
            res: res_rr,
            h_row: build_h_row_doppler(los, x_i.len(), h_vel_att, h_vel_bg),
            weight: var_rr / w_rr,
            raw_var: var_rr,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    pub(super) fn push_cp_measurement(
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
        h_pos_att: &nalgebra::Matrix3<f64>,
    ) {
        if let Some(amb_idx) = find_ambiguity_index(state, sat.sat_obs.sat) {
            let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
            let l_meas = (cp1 - windup) * sat.lam1;
            let expected_cp = if sat.is_iono_free {
                expected_base + x_i[CORE_STATE_SIZE + amb_idx]
            } else {
                expected_base - sat.iono_delay + x_i[CORE_STATE_SIZE + amb_idx]
            };
            let res_cp = l_meas - expected_cp;
            if res_cp.abs() > 100.0 && iter == 0 {
                tracing::warn!("HUGE res_cp: sat={}, l_meas={:.2}, exp={:.2}, dist={:.2}, clk={:.2}, n_amb={:.2}", sat.sat_obs.sat, l_meas, expected_cp, dist, x_i[15], x_i[CORE_STATE_SIZE + amb_idx]);
            }
            let mut var_cp = 0.0001 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
            if sat.is_iono_free {
                var_cp *= 9.0;
            }
            let w_cp = apply_huber(res_cp, var_cp, self.huber_k);
            meas.push(FgMeasurement {
                res: res_cp,
                h_row: build_h_row(
                    los,
                    sat.map_wet,
                    Some(CORE_STATE_SIZE + amb_idx),
                    x_i.len(),
                    sat.sat_obs.sat.constellation,
                    h_pos_att,
                ),
                weight: var_cp / w_cp,
                raw_var: var_cp,
                is_phase: true,
                sat: Some(sat.sat_obs.sat),
            });
        }
    }

    pub(super) fn resolve_uduc_indices(
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
    ) -> UducIndices {
        let i1_idx = find_amb_idx(state, sat.sat_obs.sat, 3).map(|idx| CORE_STATE_SIZE + idx);
        let n1_idx = find_amb_idx(state, sat.sat_obs.sat, 1).map(|idx| CORE_STATE_SIZE + idx);
        let n2_idx = find_amb_idx(state, sat.sat_obs.sat, 2).map(|idx| CORE_STATE_SIZE + idx);
        UducIndices {
            i1: i1_idx.map(|idx| x_i[idx]).unwrap_or(0.0),
            n1: n1_idx.map(|idx| x_i[idx]).unwrap_or(0.0),
            n2: n2_idx.map(|idx| x_i[idx]).unwrap_or(0.0),
            i1_idx,
            n1_idx,
            n2_idx,
            gamma: (sat.f1 * sat.f1) / (sat.f2 * sat.f2),
        }
    }

    pub(super) fn push_uduc_pr_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
        expected_base: f64,
        h_pos_att: &nalgebra::Matrix3<f64>,
        idx: &UducIndices,
    ) {
        let var_p1 = PSEUDORANGE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        let res_p1 = sat.p1 - (expected_base + idx.i1);
        meas.push(FgMeasurement {
            res: res_p1,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                idx.i1_idx,
                1.0,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
                h_pos_att,
            ),
            weight: var_p1 / apply_huber(res_p1, var_p1, self.huber_k),
            raw_var: var_p1,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
        let res_p2 = sat.p2.expect("dual-frequency measurement ensures p2 is Some") - (expected_base + idx.gamma * idx.i1);
        meas.push(FgMeasurement {
            res: res_p2,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                idx.i1_idx,
                idx.gamma,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
                h_pos_att,
            ),
            weight: (var_p1 * 1.5) / apply_huber(res_p2, var_p1 * 1.5, self.huber_k),
            raw_var: var_p1 * 1.5,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    pub(super) fn push_uduc_cp_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
        expected_base: f64,
        h_pos_att: &nalgebra::Matrix3<f64>,
        idx: &UducIndices,
    ) {
        let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
        let var_l1 = 0.0001 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        let res_l1 = (sat.cp1.expect("dual-frequency measurement ensures cp1 is Some") - windup) * sat.lam1 - (expected_base - idx.i1 + idx.n1);
        meas.push(FgMeasurement {
            res: res_l1,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                idx.i1_idx,
                -1.0,
                idx.n1_idx,
                x_i.len(),
                sat.sat_obs.sat.constellation,
                h_pos_att,
            ),
            weight: var_l1 / apply_huber(res_l1, var_l1, self.huber_k),
            raw_var: var_l1,
            is_phase: true,
            sat: Some(sat.sat_obs.sat),
        });
        let res_l2 =
            (sat.cp2.expect("dual-frequency measurement ensures cp2 is Some") - windup) * sat.lam2 - (expected_base - idx.gamma * idx.i1 + idx.n2);
        meas.push(FgMeasurement {
            res: res_l2,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                idx.i1_idx,
                -idx.gamma,
                idx.n2_idx,
                x_i.len(),
                sat.sat_obs.sat.constellation,
                h_pos_att,
            ),
            weight: (var_l1 * 1.5) / apply_huber(res_l2, var_l1 * 1.5, self.huber_k),
            raw_var: var_l1 * 1.5,
            is_phase: true,
            sat: Some(sat.sat_obs.sat),
        });
    }

    pub(super) fn push_uduc_measurements(
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
        h_pos_att: &nalgebra::Matrix3<f64>,
    ) {
        let idx = Self::resolve_uduc_indices(state, sat, x_i);
        self.push_uduc_pr_measurements(meas, sat, x_i, los, expected_base, h_pos_att, &idx);
        self.push_uduc_cp_measurements(meas, state, sat, x_i, los, expected_base, h_pos_att, &idx);
    }
}

pub fn build_h_row_uduc(
    los: &Vector3<f64>,
    map_wet: f64,
    i_idx: Option<usize>,
    i_coef: f64,
    n_idx: Option<usize>,
    size: usize,
    constel: gneiss_core::sat::Constellation,
    h_pos_att: &nalgebra::Matrix3<f64>,
) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[0] = -los.x;
    h[1] = -los.y;
    h[2] = -los.z;
    if size > 15 {
        let pr_h_att = -los.transpose() * h_pos_att;
        h[6] = pr_h_att[0];
        h[7] = pr_h_att[1];
        h[8] = pr_h_att[2];
        h[15] = 1.0;
    }
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

pub fn build_h_row(
    los: &Vector3<f64>,
    map_wet: f64,
    amb_idx: Option<usize>,
    size: usize,
    constel: gneiss_core::sat::Constellation,
    h_pos_att: &nalgebra::Matrix3<f64>,
) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[0] = -los.x;
    h[1] = -los.y;
    h[2] = -los.z;
    if size > 15 {
        let pr_h_att = -los.transpose() * h_pos_att;
        h[6] = pr_h_att[0];
        h[7] = pr_h_att[1];
        h[8] = pr_h_att[2];
        h[15] = 1.0;
    }
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

pub fn build_h_row_doppler(
    los: &Vector3<f64>,
    size: usize,
    h_vel_att: &nalgebra::Matrix3<f64>,
    h_vel_bg: &nalgebra::Matrix3<f64>,
) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    if size > 19 {
        h[3] = -los.x;
        h[4] = -los.y;
        h[5] = -los.z;
        let dop_h_att = -los.transpose() * h_vel_att;
        let dop_h_bg = -los.transpose() * h_vel_bg;
        h[6] = dop_h_att[0];
        h[7] = dop_h_att[1];
        h[8] = dop_h_att[2];
        h[12] = dop_h_bg[0];
        h[13] = dop_h_bg[1];
        h[14] = dop_h_bg[2];
        h[19] = 1.0;
    }
    h
}
