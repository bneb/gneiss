//! Tightly-Coupled Network RTK/INS Integration Pipeline.
//!
//! Fuses localized Virtual Reference Station (VRS) synthesized observables and
//! double-difference carrier-phase / pseudorange observations with a 15-state
//! Error-State Kalman Filter (ESKF) and LAMBDA ambiguity resolution.

use std::collections::HashMap;

use nalgebra::{
    DMatrix, DVector, Matrix2, Matrix3, RowVector3, SMatrix, SVector, Vector2, Vector3,
};

use gneiss_core::coords::{az_el, ecef_to_llh};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::imu::ImuMeasurement;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

use crate::ambiguity::lambda::resolve_lambda;
use crate::composite::{CompositeMode, EpochObservation, NavSolution, StationEpoch};
use crate::estimators::eskf::{
    apply_error_injection, joseph_form_update, predict, skew_symmetric,
    update_nhc, update_zupt, EngineError, EskfSmoother, EskfSnapshot, EskfState,
    Matrix15, Vector15,
};
use crate::estimators::rtk_iekf::DoubleDiffKey;
use crate::post_process::network_adj::CorsStation;
use crate::post_process::vrs::VrsSynthesizer;

const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Tracked double-difference carrier-phase ambiguity state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleDiffAmbiguity {
    pub key: DoubleDiffKey,
    pub float_cycles: f64,
    pub var_cycles2: f64,
    pub fixed_integer: Option<i32>,
    pub lock_count: u32,
}

/// Configuration parameters for Tightly-Coupled Network RTK/INS pipeline.
#[derive(Debug, Clone)]
pub struct TcRtkConfig {
    pub lever_arm: Vector3<f64>,
    pub q_diag: Vector15<f64>,
    pub sigma_dd_code: f64,
    pub sigma_dd_phase: f64,
    pub min_elevation_rad: f64,
    pub ar_ratio_threshold: f64,
    pub max_lateral_speed_nhc: f64,
    pub enable_nhc: bool,
    pub enable_zupt: bool,
    pub r_nhc: Matrix2<f64>,
    pub r_zupt: Matrix3<f64>,
    pub enable_smoother: bool,
}

impl Default for TcRtkConfig {
    fn default() -> Self {
        Self {
            lever_arm: Vector3::zeros(),
            q_diag: Vector15::from_column_slice(&[1e-4, 1e-4, 1e-4, 1e-2, 1e-2, 1e-2, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8]),
            sigma_dd_code: 0.30,
            sigma_dd_phase: 0.003,
            min_elevation_rad: 10.0_f64.to_radians(),
            ar_ratio_threshold: 2.5,
            max_lateral_speed_nhc: 1.0,
            enable_nhc: true,
            enable_zupt: true,
            r_nhc: Matrix2::from_diagonal(&Vector2::new(0.25, 0.25)),
            r_zupt: Matrix3::from_diagonal(&Vector3::new(0.001, 0.001, 0.001)),
            enable_smoother: false,
        }
    }
}

/// Tightly-Coupled Network RTK/INS navigation pipeline.
pub struct TightlyCoupledNetworkRtkIns {
    pub eskf: EskfState,
    pub vrs_synth: VrsSynthesizer,
    pub config: TcRtkConfig,
    pub ephemerides: Vec<Ephemeris>,
    pub tracked_ambiguities: HashMap<DoubleDiffKey, DoubleDiffAmbiguity>,
    pub last_time: Option<GpsTime>,
    pub smoother: Option<EskfSmoother>,
}

impl TightlyCoupledNetworkRtkIns {
    pub fn new(eskf: EskfState, vrs_synth: VrsSynthesizer, config: TcRtkConfig) -> Self {
        let smoother = if config.enable_smoother { Some(EskfSmoother::new()) } else { None };
        Self {
            eskf,
            vrs_synth,
            config,
            ephemerides: Vec::new(),
            tracked_ambiguities: HashMap::new(),
            last_time: None,
            smoother,
        }
    }

    pub fn with_ephemerides(mut self, eph: Vec<Ephemeris>) -> Self {
        self.ephemerides = eph;
        self
    }

    pub fn process_epoch(
        &mut self,
        imu_samples: &[ImuMeasurement],
        rover_obs: &EpochObservation,
        cors_obs: &[StationEpoch],
    ) -> Result<NavSolution, EngineError> {
        let phi = self.predict_inertial(imu_samples, rover_obs.time)?;
        self.apply_kinematic_constraints()?;

        let pred_state = self.eskf.clone();
        let ref_epoch = self.synthesize_vrs_reference(cors_obs)?;
        let (y, h_rows, r_diag, num_sats) = match &ref_epoch {
            Some(base) => self.formulate_dd_measurements(rover_obs, base),
            None => (Vec::new(), Vec::new(), Vec::new(), 0),
        };

        let (is_fixed, ratio) = self.attempt_lambda_ar();
        if !y.is_empty() {
            self.apply_kalman_measurement_update(&y, &h_rows, &r_diag, is_fixed)?;
        }

        if let Some(s) = &mut self.smoother {
            s.push(EskfSnapshot {
                time: rover_obs.time,
                state_pred: pred_state,
                state_post: self.eskf.clone(),
                phi: phi.unwrap_or_else(Matrix15::identity),
                is_gnss_available: num_sats >= 4,
            });
        }
        self.last_time = Some(rover_obs.time);
        Ok(self.build_solution(rover_obs.time, num_sats, is_fixed, ratio))
    }

    fn predict_inertial(
        &mut self,
        imu_samples: &[ImuMeasurement],
        epoch_time: GpsTime,
    ) -> Result<Option<Matrix15<f64>>, EngineError> {
        if imu_samples.is_empty() {
            return self.predict_dead_reckon_gap(epoch_time);
        }
        let mut last_phi = None;
        let mut last_tag = self.last_time.map_or(imu_samples[0].time_tag, |t| (t.tow * 1000.0) as u32);
        for imu in imu_samples {
            let dt_ms = (imu.time_tag.wrapping_sub(last_tag)) as f64;
            let dt = if dt_ms > 0.0 && dt_ms < 5000.0 { dt_ms * 1e-3 } else { 0.01 };
            predict(&mut self.eskf, imu, dt, &self.config.q_diag)?;
            last_tag = imu.time_tag;
            last_phi = Some(Matrix15::identity());
        }
        Ok(last_phi)
    }

    fn predict_dead_reckon_gap(
        &mut self,
        epoch_time: GpsTime,
    ) -> Result<Option<Matrix15<f64>>, EngineError> {
        let dt = self.last_time.map_or(0.1, |t| (epoch_time.tow - t.tow).clamp(0.0, 10.0));
        if dt > 1e-4 {
            let imu = ImuMeasurement::new((epoch_time.tow * 1000.0) as u32, Vector3::new(0.0, 0.0, -9.81), Vector3::zeros());
            predict(&mut self.eskf, &imu, dt, &self.config.q_diag)?;
        }
        Ok(None)
    }

    fn apply_kinematic_constraints(&mut self) -> Result<(), EngineError> {
        let r_b2e = self.eskf.attitude.to_rotation_matrix();
        let vel_b = r_b2e.transpose() * self.eskf.vel_ecef;
        if self.config.enable_zupt && self.eskf.vel_ecef.norm() < 0.05 {
            update_zupt(&mut self.eskf, &self.config.r_zupt)?;
        } else if self.config.enable_nhc && vel_b.y.abs() <= self.config.max_lateral_speed_nhc {
            update_nhc(&mut self.eskf, &self.config.lever_arm, &self.config.r_nhc)?;
        }
        Ok(())
    }

    fn synthesize_vrs_reference(
        &self,
        cors_obs: &[StationEpoch],
    ) -> Result<Option<EpochObs>, EngineError> {
        if cors_obs.is_empty() {
            return Ok(None);
        }
        let stations: Vec<CorsStation> = cors_obs.iter().map(|s| CorsStation {
            id: s.station_id.clone(),
            pos_ecef: s.station_pos,
            epochs: vec![s.obs.clone()],
        }).collect();

        let ant_pos = self.antenna_position_ecef();
        match self.vrs_synth.synthesize_vrs_stream(&stations, ant_pos, &self.ephemerides) {
            Ok(mut stream) => Ok(stream.pop()),
            Err(_) => Ok(cors_obs.first().map(|s| s.obs.clone())),
        }
    }

    fn formulate_dd_measurements(
        &mut self,
        rover_obs: &EpochObs,
        base_obs: &EpochObs,
    ) -> (Vec<f64>, Vec<SVector<f64, 15>>, Vec<f64>, usize) {
        let ant_pos = self.antenna_position_ecef();
        let l_e = self.eskf.attitude.to_rotation_matrix() * self.config.lever_arm;
        let l_skew = skew_symmetric(&l_e);

        let (mut y, mut h_rows, mut r_diag) = (Vec::new(), Vec::new(), Vec::new());
        let (ref_sat, u_ref) = match self.find_ref_sat_and_los(rover_obs, &ant_pos) {
            Some(pair) => pair,
            None => return (y, h_rows, r_diag, 0),
        };

        let mut tracked_count = 1;
        for r_sat in &rover_obs.satellites {
            if r_sat.sat == ref_sat { continue; }
            let b_sat = match base_obs.satellites.iter().find(|s| s.sat == r_sat.sat) {
                Some(s) => s,
                None => continue,
            };
            let u_i = match self.compute_sat_unit_vector(r_sat.sat, rover_obs.time, &ant_pos) {
                Some(u) => u,
                None => continue,
            };
            let (dd_geom, dd_code, dd_phase) = match self.calculate_dd_values(r_sat, b_sat, ref_sat, base_obs, rover_obs, &ant_pos) {
                Some(res) => res,
                None => continue,
            };
            tracked_count += 1;
            let (delta_u, h_dd_pos) = compute_dd_los_jacobian(&u_ref, &u_i);
            let h_dd_att = compute_dd_att_coupling_jacobian(&delta_u, &l_skew);
            let mut h_row = SVector::<f64, 15>::zeros();
            h_row.fixed_rows_mut::<3>(0).copy_from(&h_dd_pos.transpose());
            h_row.fixed_rows_mut::<3>(6).copy_from(&h_dd_att.transpose());

            if let Some(code_val) = dd_code {
                y.push(code_val - dd_geom);
                h_rows.push(h_row);
                r_diag.push(self.config.sigma_dd_code * self.config.sigma_dd_code);
            }
            if let Some(phase_val) = dd_phase {
                let key = DoubleDiffKey { constellation_id: r_sat.sat.constellation as u8, sat: r_sat.sat.prn as u16, ref_sat: ref_sat.prn as u16, freq_band: 1 };
                let amb = self.manage_dd_ambiguity(key, phase_val, dd_geom);
                y.push(compute_dd_residual(phase_val, dd_geom, amb));
                h_rows.push(h_row);
                r_diag.push(self.config.sigma_dd_phase * self.config.sigma_dd_phase);
            }
        }
        (y, h_rows, r_diag, tracked_count)
    }

    fn calculate_dd_values(
        &self,
        r_sat: &SatObs,
        b_sat: &SatObs,
        ref_sat: SatelliteId,
        base_obs: &EpochObs,
        rover_obs: &EpochObs,
        ant_pos: &Vector3<f64>,
    ) -> Option<(f64, Option<f64>, Option<f64>)> {
        let r_ref = rover_obs.satellites.iter().find(|s| s.sat == ref_sat)?;
        let b_ref = base_obs.satellites.iter().find(|s| s.sat == ref_sat)?;
        let (p_i, p_ref) = (self.get_sat_position(r_sat.sat, rover_obs.time)?, self.get_sat_position(ref_sat, rover_obs.time)?);
        let base_pos = self.vrs_synth.master_pos;

        let dd_geom = ((p_i - ant_pos).norm() - (p_i - base_pos).norm())
            - ((p_ref - ant_pos).norm() - (p_ref - base_pos).norm());
        let dd_code = match (r_sat.get_observable(1), b_sat.get_observable(1), r_ref.get_observable(1), b_ref.get_observable(1)) {
            (Some(ri), Some(bi), Some(r0), Some(b0)) => Some((ri - bi) - (r0 - b0)),
            _ => None,
        };
        let wl = get_carrier_wavelength(r_sat.sat.constellation);
        let dd_phase = match (r_sat.get_observable_phase(1), b_sat.get_observable_phase(1), r_ref.get_observable_phase(1), b_ref.get_observable_phase(1)) {
            (Some(ri), Some(bi), Some(r0), Some(b0)) => Some(((ri - bi) - (r0 - b0)) * wl),
            _ => None,
        };
        Some((dd_geom, dd_code, dd_phase))
    }

    fn manage_dd_ambiguity(&mut self, key: DoubleDiffKey, phase_m: f64, geom_m: f64) -> f64 {
        let entry = self.tracked_ambiguities.entry(key).or_insert_with(|| DoubleDiffAmbiguity {
            key,
            float_cycles: (phase_m - geom_m) / get_carrier_wavelength(Constellation::Gps),
            var_cycles2: 100.0,
            fixed_integer: None,
            lock_count: 0,
        });
        entry.lock_count = entry.lock_count.saturating_add(1);
        let wl = get_carrier_wavelength(Constellation::Gps);
        entry.fixed_integer.map_or(entry.float_cycles * wl, |n| (n as f64) * wl)
    }

    fn attempt_lambda_ar(&mut self) -> (bool, Option<f64>) {
        let n = self.tracked_ambiguities.len();
        if n < 4 {
            return (false, None);
        }
        let mut float_vec = DVector::zeros(n);
        let mut cov_mat = DMatrix::zeros(n, n);
        let mut keys = Vec::with_capacity(n);

        for (idx, (k, a)) in self.tracked_ambiguities.iter().enumerate() {
            float_vec[idx] = a.float_cycles;
            cov_mat[(idx, idx)] = a.var_cycles2.clamp(1e-4, 1.0);
            keys.push(*k);
        }
        match resolve_lambda(&float_vec, &cov_mat) {
            Ok(res) if res.ratio >= self.config.ar_ratio_threshold => {
                for (idx, k) in keys.iter().enumerate() {
                    if let Some(a) = self.tracked_ambiguities.get_mut(k) {
                        a.fixed_integer = Some(res.best_integers[idx] as i32);
                        a.var_cycles2 = 1e-4;
                    }
                }
                (true, Some(res.ratio))
            }
            Ok(res) => (false, Some(res.ratio)),
            Err(_) => (false, None),
        }
    }

    fn apply_kalman_measurement_update(
        &mut self,
        y: &[f64],
        h_rows: &[SVector<f64, 15>],
        r_diag: &[f64],
        is_fixed: bool,
    ) -> Result<(), EngineError> {
        let m = y.len().min(30);
        let mut h_mat = SMatrix::<f64, 30, 15>::zeros();
        let mut y_vec = SVector::<f64, 30>::zeros();
        let mut r_mat = SMatrix::<f64, 30, 30>::zeros();

        for i in 0..m {
            y_vec[i] = y[i];
            let r_val = if is_fixed { r_diag[i] * 0.05 } else { r_diag[i] };
            r_mat[(i, i)] = r_val;
            for j in 0..15 {
                h_mat[(i, j)] = h_rows[i][j];
            }
        }
        let h_slice = h_mat.fixed_view::<30, 15>(0, 0);
        let s = h_slice * self.eskf.cov * h_slice.transpose() + r_mat;
        let s_inv = s.try_inverse().ok_or(EngineError::InversionError)?;
        let k = self.eskf.cov * h_slice.transpose() * s_inv;

        let dx = k * y_vec;
        apply_error_injection(&mut self.eskf, &dx);
        self.eskf.cov = joseph_form_update(&self.eskf.cov, &h_slice.into_owned(), &k, &r_mat);
        Ok(())
    }

    fn find_ref_sat_and_los(&self, obs: &EpochObs, ant_pos: &Vector3<f64>) -> Option<(SatelliteId, Vector3<f64>)> {
        let llh = ecef_to_llh(*ant_pos);
        let mut best = None;
        let mut best_el = -1.0;

        for s in &obs.satellites {
            let sat_pos = match self.get_sat_position(s.sat, obs.time) {
                Some(p) => p,
                None => continue,
            };
            let (_, el) = az_el(llh, *ant_pos, sat_pos);
            if el > best_el && el >= self.config.min_elevation_rad {
                best_el = el;
                let diff = sat_pos - ant_pos;
                let u = diff / diff.norm();
                best = Some((s.sat, u));
            }
        }
        best
    }

    fn compute_sat_unit_vector(&self, sat: SatelliteId, time: GpsTime, ant_pos: &Vector3<f64>) -> Option<Vector3<f64>> {
        let sat_pos = self.get_sat_position(sat, time)?;
        let diff = sat_pos - ant_pos;
        let d = diff.norm();
        if d > 1e-3 { Some(diff / d) } else { None }
    }

    fn get_sat_position(&self, sat: SatelliteId, time: GpsTime) -> Option<Vector3<f64>> {
        let eph = self.ephemerides.iter().find(|e| e.sat() == sat)?;
        let (pos, _, _, _) = eph.position(time);
        Some(pos)
    }

    fn antenna_position_ecef(&self) -> Vector3<f64> {
        let r_b2e = self.eskf.attitude.to_rotation_matrix();
        self.eskf.pos_ecef + r_b2e * self.config.lever_arm
    }

    fn build_solution(
        &self,
        time: GpsTime,
        num_sats: usize,
        is_fixed: bool,
        ratio: Option<f64>,
    ) -> NavSolution {
        NavSolution {
            time,
            pos_ecef: self.eskf.pos_ecef,
            vel_ecef: self.eskf.vel_ecef,
            attitude: self.eskf.attitude,
            accel_bias: self.eskf.accel_bias,
            gyro_bias: self.eskf.gyro_bias,
            cov: self.eskf.cov,
            num_satellites: num_sats,
            is_fixed,
            mode: if num_sats < 4 { CompositeMode::DeadReckoning } else { CompositeMode::TightlyCoupledRtk },
            ratio,
        }
    }
}

/// Double-difference LOS difference position Jacobian: H_dd = -(u_i - u_ref)^T.
pub fn compute_dd_los_jacobian(
    u_ref: &Vector3<f64>,
    u_i: &Vector3<f64>,
) -> (Vector3<f64>, RowVector3<f64>) {
    let delta_u = u_i - u_ref;
    let h_dd = -delta_u.transpose();
    (delta_u, h_dd)
}

/// Double-difference attitude coupling Jacobian with lever arm: H_att = -delta_u^T * [l_e x].
pub fn compute_dd_att_coupling_jacobian(
    delta_u: &Vector3<f64>,
    l_skew: &Matrix3<f64>,
) -> RowVector3<f64> {
    -delta_u.transpose() * l_skew
}

/// Double-difference carrier phase residual: res = dd_meas - (dd_geom + dd_amb_m).
pub fn compute_dd_residual(dd_meas: f64, dd_geom: f64, dd_amb_m: f64) -> f64 {
    dd_meas - (dd_geom + dd_amb_m)
}

fn get_carrier_wavelength(constellation: Constellation) -> f64 {
    match constellation {
        Constellation::Gps | Constellation::Galileo => SPEED_OF_LIGHT / 1575.42e6,
        Constellation::Beidou => SPEED_OF_LIGHT / 1561.098e6,
        Constellation::Glonass => SPEED_OF_LIGHT / 1602.0e6,
        _ => SPEED_OF_LIGHT / 1575.42e6,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tc_rtk_dd_los_and_attitude_jacobians() {
        let u_s1 = Vector3::new(0.5, 0.5, 0.70);
        let u_s2 = Vector3::new(0.3, 0.7, 0.65);
        let (_delta_u, h_dd) = compute_dd_los_jacobian(&u_s1, &u_s2);
        assert!((h_dd[0] - 0.2_f64).abs() < 1e-12);
        assert!((h_dd[1] - (-0.2_f64)).abs() < 1e-12);

        let delta_u_test = Vector3::new(0.2, -0.2, 0.0);
        let lever_arm = Vector3::new(0.0, 0.0, -1.0);
        let l_skew = skew_symmetric(&lever_arm);
        let h_att = compute_dd_att_coupling_jacobian(&delta_u_test, &l_skew);
        assert!((h_att[0] - (-0.2)).abs() < 1e-12);
        assert!((h_att[1] - (-0.2)).abs() < 1e-12);
    }

    #[test]
    fn test_tc_rtk_double_difference_residual_math() {
        let dd_meas = 45.123;
        let dd_geom = 45.000;
        let dd_amb_m = 0.120;
        let res = compute_dd_residual(dd_meas, dd_geom, dd_amb_m);
        assert!((res - 0.003_f64).abs() < 1e-6);
    }

    #[test]
    fn test_tc_rtk_vrs_packet_loss_graceful_propagation() {
        let eskf = EskfState::new(Vector3::new(10.0, 20.0, 30.0), Vector3::zeros(), nalgebra::UnitQuaternion::identity());
        let synth = VrsSynthesizer::new("BASE", Vector3::zeros());
        let rtk = TightlyCoupledNetworkRtkIns::new(eskf, synth, TcRtkConfig::default());

        let res = rtk.synthesize_vrs_reference(&[]).unwrap();
        assert!(res.is_none());
    }
}
