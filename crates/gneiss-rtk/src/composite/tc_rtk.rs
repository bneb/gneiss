//! Tightly-Coupled Network RTK/INS Integration Pipeline.
use std::collections::HashMap;
use nalgebra::{Matrix2, Matrix3, RowVector3, SVector, Vector2, Vector3};

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::imu::ImuMeasurement;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

use crate::composite::tc_ambiguity::{
    find_group_ref_sat, ConstellationGroup, GroupRefSat, TcAmbiguityTracker,
};
use crate::composite::{CompositeMode, EpochObservation, NavSolution, StationEpoch};
use crate::estimators::eskf::predict::predict_with_phi;
use crate::estimators::eskf::{
    skew_symmetric, update_nhc, update_zupt, EngineError, EskfSmoother, EskfSnapshot, EskfState,
    Matrix15, Vector15,
};
use crate::estimators::rtk_iekf::ref_sat::sat_to_prn_u16;
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

struct EpochDdPair {
    key: DoubleDiffKey,
    h_x: SVector<f64, 15>,
    wavelength: f64,
    dd_code: Option<f64>,
    dd_phase: Option<f64>,
    dd_geom: f64,
    p_i: Vector3<f64>,
    p_ref: Vector3<f64>,
}

/// Tightly-Coupled Network RTK/INS navigation pipeline.
pub struct TightlyCoupledNetworkRtkIns {
    pub eskf: EskfState,
    pub vrs_synth: VrsSynthesizer,
    pub config: TcRtkConfig,
    pub ephemerides: Vec<Ephemeris>,
    pub tracker: TcAmbiguityTracker,
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
            tracker: TcAmbiguityTracker::new(),
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
        let pairs = match &ref_epoch {
            Some(base) => self.formulate_and_apply_dd(rover_obs, base),
            None => Vec::new(),
        };

        let (is_fixed, ratio) = self.attempt_ar_with_screening(&pairs);
        self.sync_tracked_ambiguities_map();

        let num_sats = pairs.len().saturating_add(if pairs.is_empty() { 0 } else { 1 });
        self.update_smoother_and_history(pred_state, phi, rover_obs.time, num_sats);
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
            let phi = predict_with_phi(&mut self.eskf, &imu.accel, &imu.gyro, dt, &self.config.q_diag)?;
            self.tracker.propagate_cross_cov(&phi);
            last_tag = imu.time_tag;
            last_phi = Some(phi);
        }
        Ok(last_phi)
    }

    fn predict_dead_reckon_gap(
        &mut self,
        epoch_time: GpsTime,
    ) -> Result<Option<Matrix15<f64>>, EngineError> {
        let dt = self.last_time.map_or(0.1, |t| (epoch_time.tow - t.tow).clamp(0.0, 10.0));
        if dt > 1e-4 {
            let accel = Vector3::new(0.0, 0.0, -9.81);
            let phi = predict_with_phi(&mut self.eskf, &accel, &Vector3::zeros(), dt, &self.config.q_diag)?;
            self.tracker.propagate_cross_cov(&phi);
            return Ok(Some(phi));
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

    fn synthesize_vrs_reference(&self, cors_obs: &[StationEpoch]) -> Result<Option<EpochObs>, EngineError> {
        if cors_obs.is_empty() { return Ok(None); }
        let stations: Vec<CorsStation> = cors_obs.iter().map(|s| CorsStation {
            id: s.station_id.clone(), pos_ecef: s.station_pos, epochs: vec![s.obs.clone()],
        }).collect();
        let ant_pos = self.antenna_position_ecef();
        match self.vrs_synth.synthesize_vrs_stream(&stations, ant_pos, &self.ephemerides) {
            Ok(mut stream) => Ok(stream.pop()),
            Err(_) => Ok(cors_obs.first().map(|s| s.obs.clone())),
        }
    }

    fn formulate_and_apply_dd(&mut self, rover_obs: &EpochObs, base_obs: &EpochObs) -> Vec<EpochDdPair> {
        let pairs = self.collect_all_dd_pairs(rover_obs, base_obs);
        if pairs.is_empty() { return pairs; }
        self.apply_dd_float_updates(&pairs);
        pairs
    }

    fn collect_all_dd_pairs(&self, rover_obs: &EpochObs, base_obs: &EpochObs) -> Vec<EpochDdPair> {
        let ant_pos = self.antenna_position_ecef();
        let l_e = self.eskf.attitude.to_rotation_matrix() * self.config.lever_arm;
        let l_skew = skew_symmetric(&l_e);
        let mut all_pairs = Vec::new();
        for group in ConstellationGroup::ALL {
            all_pairs.extend(self.collect_group_dd(group, rover_obs, base_obs, &ant_pos, &l_skew));
        }
        all_pairs
    }

    fn collect_group_dd(
        &self,
        group: ConstellationGroup,
        rover_obs: &EpochObs,
        base_obs: &EpochObs,
        ant_pos: &Vector3<f64>,
        l_skew: &Matrix3<f64>,
    ) -> Vec<EpochDdPair> {
        let Some(ref_sat) = find_group_ref_sat(
            group, rover_obs, base_obs, &self.ephemerides, ant_pos, self.config.min_elevation_rad,
        ) else { return Vec::new(); };
        let mut pairs = Vec::new();
        for r_sat in &rover_obs.satellites {
            if !group.matches(r_sat.sat) || r_sat.sat == ref_sat.sat { continue; }
            if let Some(pair) = self.build_dd_pair(r_sat, base_obs, rover_obs, ref_sat, ant_pos, l_skew) {
                pairs.push(pair);
            }
        }
        pairs
    }

    fn build_dd_pair(
        &self,
        r_sat: &SatObs,
        base_obs: &EpochObs,
        rover_obs: &EpochObs,
        ref_sat: GroupRefSat,
        ant_pos: &Vector3<f64>,
        l_skew: &Matrix3<f64>,
    ) -> Option<EpochDdPair> {
        let b_sat = base_obs.satellites.iter().find(|s| s.sat == r_sat.sat)?;
        let u_i = self.compute_sat_unit_vector(r_sat.sat, rover_obs.time, ant_pos)?;
        let (dd_geom, dd_code, dd_phase) = self.calculate_dd_values(r_sat, b_sat, ref_sat.sat, base_obs, rover_obs, ant_pos)?;
        let p_i = self.get_sat_position(r_sat.sat, rover_obs.time)?;
        let (delta_u, h_dd_pos) = compute_dd_los_jacobian(&ref_sat.u_ref, &u_i);
        let h_dd_att = compute_dd_att_coupling_jacobian(&delta_u, l_skew);
        let mut h_x = SVector::<f64, 15>::zeros();
        h_x.fixed_rows_mut::<3>(0).copy_from(&h_dd_pos.transpose());
        h_x.fixed_rows_mut::<3>(6).copy_from(&h_dd_att.transpose());
        let const_grp = ConstellationGroup::ALL.iter().find(|g| g.matches(r_sat.sat))?;
        let key = DoubleDiffKey {
            constellation_id: const_grp.constellation_id(),
            sat: sat_to_prn_u16(r_sat.sat),
            ref_sat: sat_to_prn_u16(ref_sat.sat),
            freq_band: 1,
        };
        let wavelength = get_carrier_wavelength(r_sat.sat.constellation);
        Some(EpochDdPair { key, h_x, wavelength, dd_code, dd_phase, dd_geom, p_i, p_ref: ref_sat.pos })
    }

    fn apply_dd_float_updates(&mut self, pairs: &[EpochDdPair]) {
        let active_keys: Vec<DoubleDiffKey> = pairs.iter().filter(|p| p.dd_phase.is_some()).map(|p| p.key).collect();
        let init_floats: Vec<f64> = pairs.iter().filter(|p| p.dd_phase.is_some()).map(|p| {
            let phase = p.dd_phase.unwrap_or(0.0);
            (phase - p.dd_geom) / p.wavelength
        }).collect();
        self.tracker.sync_keys(&active_keys, &init_floats);

        let var_code = self.config.sigma_dd_code * self.config.sigma_dd_code;
        for p in pairs {
            if let Some(code) = p.dd_code {
                self.tracker.update_code_float(&mut self.eskf, &p.h_x, code - p.dd_geom, var_code);
            }
        }
        let var_phase = self.config.sigma_dd_phase * self.config.sigma_dd_phase;
        for p in pairs {
            let (Some(phase), Some(idx)) = (p.dd_phase, self.tracker.key_index(&p.key)) else {
                continue;
            };
            let y = phase - p.dd_geom - self.tracker.a_float[idx] * p.wavelength;
            self.tracker.update_carrier_float(&mut self.eskf, idx, &p.h_x, p.wavelength, y, var_phase);
        }
    }

    fn attempt_ar_with_screening(&mut self, pairs: &[EpochDdPair]) -> (bool, Option<f64>) {
        let base_pos = self.vrs_synth.master_pos;
        let l_e = self.eskf.attitude.to_rotation_matrix() * self.config.lever_arm;
        let checker = |ant_pos: &Vector3<f64>, key: DoubleDiffKey, int_val: i32| -> Option<f64> {
            let p = pairs.iter().find(|pair| pair.key == key)?;
            let phase = p.dd_phase?;
            let cur_ant = *ant_pos + l_e;
            let geom = ((p.p_i - cur_ant).norm() - (p.p_i - base_pos).norm())
                - ((p.p_ref - cur_ant).norm() - (p.p_ref - base_pos).norm());
            Some(phase - geom - (int_val as f64) * p.wavelength)
        };
        self.tracker.attempt_ar_and_condition(&mut self.eskf, checker)
    }

    fn sync_tracked_ambiguities_map(&mut self) {
        for (idx, key) in self.tracker.keys.iter().enumerate() {
            let float_cycles = self.tracker.a_float[idx];
            let var_cycles2 = self.tracker.q_aa[(idx, idx)];
            let fixed_integer = self.tracker.fixed_integers.get(key).copied();
            let lock_count = self.tracker.lock_counts.get(key).copied().unwrap_or(1);
            self.tracked_ambiguities.insert(*key, DoubleDiffAmbiguity {
                key: *key, float_cycles, var_cycles2, fixed_integer, lock_count,
            });
        }
        self.tracked_ambiguities.retain(|k, _| self.tracker.keys.contains(k));
    }

    fn update_smoother_and_history(
        &mut self,
        pred_state: EskfState,
        phi: Option<Matrix15<f64>>,
        time: GpsTime,
        num_sats: usize,
    ) {
        if let Some(s) = &mut self.smoother {
            s.push(EskfSnapshot {
                time,
                state_pred: pred_state,
                state_post: self.eskf.clone(),
                phi: phi.unwrap_or_else(Matrix15::identity),
                is_gnss_available: num_sats >= 4,
            });
        }
        self.last_time = Some(time);
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
        let p_i = self.get_sat_position(r_sat.sat, rover_obs.time)?;
        let p_ref = self.get_sat_position(ref_sat, rover_obs.time)?;
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

/// Double-difference attitude coupling Jacobian with lever arm: H_att = delta_u^T * [l_e x].
pub fn compute_dd_att_coupling_jacobian(
    delta_u: &Vector3<f64>,
    l_skew: &Matrix3<f64>,
) -> RowVector3<f64> {
    delta_u.transpose() * l_skew
}

/// Double-difference carrier phase residual: res = dd_meas - (dd_geom + dd_amb_m).
pub fn compute_dd_residual(dd_meas: f64, dd_geom: f64, dd_amb_m: f64) -> f64 {
    dd_meas - (dd_geom + dd_amb_m)
}

pub fn get_carrier_wavelength(constellation: Constellation) -> f64 {
    match constellation {
        Constellation::Gps | Constellation::Qzss | Constellation::Galileo => SPEED_OF_LIGHT / 1575.42e6,
        Constellation::Beidou => SPEED_OF_LIGHT / 1561.098e6,
        Constellation::Glonass => SPEED_OF_LIGHT / 1602.0e6,
        _ => SPEED_OF_LIGHT / 1575.42e6,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::UnitQuaternion;

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
        assert!((h_att[0] - 0.2).abs() < 1e-12);
        assert!((h_att[1] - 0.2).abs() < 1e-12);
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
        let eskf = EskfState::new(Vector3::new(10.0, 20.0, 30.0), Vector3::zeros(), UnitQuaternion::identity());
        let synth = VrsSynthesizer::new("BASE", Vector3::zeros());
        let rtk = TightlyCoupledNetworkRtkIns::new(eskf, synth, TcRtkConfig::default());

        let res = rtk.synthesize_vrs_reference(&[]).unwrap();
        assert!(res.is_none());
    }

    #[test]
    fn test_tc_rtk_constellation_isolation_in_dd_formation() {
        let eskf = EskfState::new(Vector3::new(100.0, 200.0, 300.0), Vector3::zeros(), UnitQuaternion::identity());
        let synth = VrsSynthesizer::new("BASE", Vector3::new(100.0, 200.0, 300.0));
        let rtk = TightlyCoupledNetworkRtkIns::new(eskf, synth, TcRtkConfig::default());
        let sats = vec![
            SatObs { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, observations: Vec::new() },
            SatObs { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, observations: Vec::new() },
            SatObs { sat: SatelliteId { constellation: Constellation::Galileo, prn: 5 }, observations: Vec::new() },
        ];
        let rover = EpochObs { time: GpsTime::new(2200, 100.0), satellites: sats };
        assert!(rtk.collect_all_dd_pairs(&rover, &rover.clone()).is_empty());
    }
}
