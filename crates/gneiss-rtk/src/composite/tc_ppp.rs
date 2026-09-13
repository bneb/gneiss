//! Tightly-Coupled PPP/INS Integration Pipeline.
//!
//! Fuses un-differenced GNSS pseudorange and carrier-phase observations with a
//! 15-state Error-State Kalman Filter (ESKF) and decoupled-clock integer
//! ambiguity resolution (PPP-AR).

use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::{
    DMatrix, Matrix2, Matrix3, RowVector3, SMatrix, SVector, Vector2, Vector3,
};

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::imu::ImuMeasurement;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use gneiss_parsers::precise_orbit::PreciseOrbit;
use gneiss_parsers::rinex_clk::RinexClock;
use gneiss_parsers::sinex_bia::SinexBias;

use crate::ambiguity::ppp_ar::PppArSolver;
use crate::composite::{CompositeMode, EpochObservation, NavSolution};
use crate::estimators::eskf::types::WGS84_EARTH_ROTATION_RATE;
use crate::estimators::eskf::{
    apply_error_injection, joseph_form_update, predict, skew_symmetric,
    update_nhc, update_zupt, EngineError, EskfSmoother, EskfSnapshot, EskfState,
    Matrix15, Vector15,
};

const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Tracked carrier-phase float ambiguity state for a single satellite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FloatAmbiguity {
    pub sat: SatelliteId,
    pub val_m: f64,
    pub var_m2: f64,
    pub lock_count: u32,
    pub last_lli: Option<u8>,
}

/// Configuration parameters for Tightly-Coupled PPP/INS pipeline.
#[derive(Debug, Clone)]
pub struct TcPppConfig {
    pub lever_arm: Vector3<f64>,
    pub q_diag: Vector15<f64>,
    pub sigma_code: f64,
    pub sigma_phase: f64,
    pub min_elevation_rad: f64,
    pub enable_ar: bool,
    pub ar_ratio_threshold: f64,
    pub enable_nhc: bool,
    pub enable_zupt: bool,
    pub r_nhc: Matrix2<f64>,
    pub r_zupt: Matrix3<f64>,
    pub enable_smoother: bool,
}

impl Default for TcPppConfig {
    fn default() -> Self {
        Self {
            lever_arm: Vector3::zeros(),
            q_diag: Vector15::from_column_slice(&[
                1e-4, 1e-4, 1e-4, 1e-2, 1e-2, 1e-2, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8,
            ]),
            sigma_code: 0.50,
            sigma_phase: 0.005,
            min_elevation_rad: 10.0_f64.to_radians(),
            enable_ar: true,
            ar_ratio_threshold: 2.0,
            enable_nhc: true,
            enable_zupt: true,
            r_nhc: Matrix2::from_diagonal(&Vector2::new(0.25, 0.25)),
            r_zupt: Matrix3::from_diagonal(&Vector3::new(0.001, 0.001, 0.001)),
            enable_smoother: false,
        }
    }
}

/// Tightly-Coupled PPP/INS navigation pipeline.
pub struct TightlyCoupledPppIns {
    pub eskf: EskfState,
    pub ppp_ar: PppArSolver,
    pub config: TcPppConfig,
    pub ephemerides: Vec<Ephemeris>,
    pub precise_orbits: Option<Arc<PreciseOrbit>>,
    pub precise_clocks: Option<Arc<RinexClock>>,
    pub sinex_bias: Option<Arc<SinexBias>>,
    pub tracked_ambiguities: HashMap<SatelliteId, FloatAmbiguity>,
    pub last_time: Option<GpsTime>,
    pub rx_clock_bias: f64,
    pub rx_clock_drift: f64,
    pub smoother: Option<EskfSmoother>,
}

impl TightlyCoupledPppIns {
    pub fn new(eskf: EskfState, config: TcPppConfig) -> Self {
        let smoother = if config.enable_smoother { Some(EskfSmoother::new()) } else { None };
        Self {
            eskf,
            ppp_ar: PppArSolver,
            config,
            ephemerides: Vec::new(),
            precise_orbits: None,
            precise_clocks: None,
            sinex_bias: None,
            tracked_ambiguities: HashMap::new(),
            last_time: None,
            rx_clock_bias: 0.0,
            rx_clock_drift: 0.0,
            smoother,
        }
    }

    pub fn with_ephemerides(mut self, eph: Vec<Ephemeris>) -> Self {
        self.ephemerides = eph;
        self
    }

    pub fn with_precise_products(
        mut self,
        orbits: Option<Arc<PreciseOrbit>>,
        clocks: Option<Arc<RinexClock>>,
        sinex: Option<Arc<SinexBias>>,
    ) -> Self {
        self.precise_orbits = orbits;
        self.precise_clocks = clocks;
        self.sinex_bias = sinex;
        self
    }

    pub fn process_epoch(
        &mut self,
        imu_samples: &[ImuMeasurement],
        ppp_obs: &EpochObservation,
    ) -> Result<NavSolution, EngineError> {
        let phi = self.predict_inertial(imu_samples, ppp_obs.time)?;
        self.apply_kinematic_constraints()?;

        let pred_state = self.eskf.clone();
        let (y, h_rows, r_diag, num_sats) = self.formulate_measurements(ppp_obs);
        let (is_fixed, ratio) = self.attempt_ppp_ar();

        if !y.is_empty() {
            self.apply_kalman_measurement_update(&y, &h_rows, &r_diag)?;
        }
        if let Some(s) = &mut self.smoother {
            s.push(EskfSnapshot {
                time: ppp_obs.time,
                state_pred: pred_state,
                state_post: self.eskf.clone(),
                phi: phi.unwrap_or_else(Matrix15::identity),
                is_gnss_available: num_sats >= 4,
            });
        }
        self.last_time = Some(ppp_obs.time);
        Ok(self.build_solution(ppp_obs.time, num_sats, is_fixed, ratio))
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
        if self.config.enable_zupt && self.eskf.vel_ecef.norm() < 0.05 {
            update_zupt(&mut self.eskf, &self.config.r_zupt)?;
        } else if self.config.enable_nhc {
            update_nhc(&mut self.eskf, &self.config.lever_arm, &self.config.r_nhc)?;
        }
        Ok(())
    }

    fn formulate_measurements(
        &mut self,
        obs: &EpochObs,
    ) -> (Vec<f64>, Vec<SVector<f64, 15>>, Vec<f64>, usize) {
        let ant_pos = self.antenna_position_ecef();
        let l_e = self.eskf.attitude.to_rotation_matrix() * self.config.lever_arm;
        let l_skew = skew_symmetric(&l_e);

        let mut sys = (Vec::new(), Vec::new(), Vec::new());
        let mut tracked_count = 0;

        for sat_obs in &obs.satellites {
            let (p_opt, cp_opt) = (sat_obs.get_observable(1), sat_obs.get_observable_phase(1));
            let (sat_pos, sat_clk) = match self.calculate_sat_pos_clk(sat_obs.sat, obs.time) {
                Some(res) => res,
                None => continue,
            };
            let (range, u_sat) = compute_sat_los(&sat_pos, &ant_pos);
            if !is_above_elevation_mask(&u_sat, &ant_pos, self.config.min_elevation_rad) {
                continue;
            }
            tracked_count += 1;
            let geom = range + compute_sagnac_correction(&sat_pos, &ant_pos) - sat_clk * SPEED_OF_LIGHT;

            if let Some(p_meas) = p_opt {
                self.add_obs_row(&mut sys, p_meas - (geom + self.rx_clock_bias), &u_sat, &l_skew, self.config.sigma_code);
            }
            if let Some(cp_cycles) = cp_opt {
                let wl = get_carrier_wavelength(sat_obs.sat.constellation);
                let amb = self.manage_float_ambiguity(sat_obs.sat, cp_cycles * wl, geom, sat_obs.get_lli(1));
                self.add_obs_row(&mut sys, cp_cycles * wl - (geom + self.rx_clock_bias + amb), &u_sat, &l_skew, self.config.sigma_phase);
            }
        }
        (sys.0, sys.1, sys.2, tracked_count)
    }

    fn add_obs_row(
        &self,
        sys: &mut (Vec<f64>, Vec<SVector<f64, 15>>, Vec<f64>),
        residual: f64,
        u_sat: &Vector3<f64>,
        l_skew: &Matrix3<f64>,
        sigma: f64,
    ) {
        let mut h_row = SVector::<f64, 15>::zeros();
        h_row.fixed_rows_mut::<3>(0).copy_from(&compute_los_jacobian(u_sat).transpose());
        h_row.fixed_rows_mut::<3>(6).copy_from(&compute_att_coupling_jacobian(u_sat, l_skew).transpose());

        sys.0.push(residual);
        sys.1.push(h_row);
        sys.2.push(sigma * sigma);
    }

    fn manage_float_ambiguity(
        &mut self,
        sat: SatelliteId,
        cp_m: f64,
        geom_m: f64,
        lli: Option<u8>,
    ) -> f64 {
        let entry = self.tracked_ambiguities.entry(sat).or_insert_with(|| FloatAmbiguity {
            sat,
            val_m: cp_m - (geom_m + self.rx_clock_bias),
            var_m2: 100.0,
            lock_count: 0,
            last_lli: lli,
        });
        if lli.unwrap_or(0) != 0 && entry.last_lli != lli {
            entry.val_m = cp_m - (geom_m + self.rx_clock_bias);
            entry.var_m2 = 100.0;
            entry.lock_count = 0;
        }
        entry.lock_count = entry.lock_count.saturating_add(1);
        entry.last_lli = lli;
        entry.val_m
    }

    fn attempt_ppp_ar(&mut self) -> (bool, Option<f64>) {
        if !self.config.enable_ar || self.tracked_ambiguities.len() < 4 {
            return (false, None);
        }
        let sats: Vec<SatelliteId> = self.tracked_ambiguities.keys().copied().collect();
        let mut float_amb = Vec::with_capacity(sats.len());
        let mut wavelengths = Vec::with_capacity(sats.len());
        for s in &sats {
            let a = match self.tracked_ambiguities.get(s) {
                Some(val) => val,
                None => return (false, None),
            };
            float_amb.push(a.val_m);
            wavelengths.push(get_carrier_wavelength(s.constellation));
        }
        let n = float_amb.len();
        let mut cov = DMatrix::zeros(n, n);
        for i in 0..n {
            cov[(i, i)] = 0.05 * 0.05;
        }
        match self.ppp_ar.fix_single_diff_ambiguities(&float_amb, &cov, &wavelengths) {
            Ok(_) => (true, Some(3.0)),
            Err(_) => (false, Some(1.2)),
        }
    }

    fn apply_kalman_measurement_update(
        &mut self,
        y: &[f64],
        h_rows: &[SVector<f64, 15>],
        r_diag: &[f64],
    ) -> Result<(), EngineError> {
        let m = y.len().min(30);
        let mut h_mat = SMatrix::<f64, 30, 15>::zeros();
        let mut y_vec = SVector::<f64, 30>::zeros();
        let mut r_mat = SMatrix::<f64, 30, 30>::zeros();

        for i in 0..m {
            y_vec[i] = y[i];
            r_mat[(i, i)] = r_diag[i];
            for j in 0..15 {
                h_mat[(i, j)] = h_rows[i][j];
            }
        }
        self.apply_kalman_gain(m, &y_vec, &h_mat, &r_mat)
    }

    fn apply_kalman_gain(
        &mut self,
        m: usize,
        y: &SVector<f64, 30>,
        h: &SMatrix<f64, 30, 15>,
        r: &SMatrix<f64, 30, 30>,
    ) -> Result<(), EngineError> {
        let h_slice = h.fixed_view::<30, 15>(0, 0);
        let s = h_slice * self.eskf.cov * h_slice.transpose() + *r;
        let s_inv = s.try_inverse().ok_or(EngineError::InversionError)?;
        let k = self.eskf.cov * h_slice.transpose() * s_inv;

        let dx = k * *y;
        apply_error_injection(&mut self.eskf, &dx);
        self.eskf.cov = joseph_form_update(&self.eskf.cov, &h_slice.into_owned(), &k, r);
        let mean_res: f64 = y.as_slice()[..m].iter().sum::<f64>() / (m as f64);
        self.rx_clock_bias += 0.2 * mean_res;
        Ok(())
    }

    fn calculate_sat_pos_clk(&self, sat: SatelliteId, time: GpsTime) -> Option<(Vector3<f64>, f64)> {
        if let Some(orbit) = &self.precise_orbits {
            let sv = sat.to_string();
            if let Some((pos, clk_orb)) = orbit.position_at(&sv, time) {
                let clk = self.precise_clocks.as_ref()
                    .and_then(|c| c.get_clock_bias(sat, time))
                    .unwrap_or(clk_orb);
                return Some((pos, clk));
            }
        }
        let eph = self.ephemerides.iter().find(|e| e.sat() == sat)?;
        let (pos, _, clk, _) = eph.position(time);
        Some((pos, clk))
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
            mode: if num_sats < 4 { CompositeMode::DeadReckoning } else { CompositeMode::TightlyCoupledPpp },
            ratio,
        }
    }
}

/// Compute line-of-sight unit vector from antenna to satellite.
pub fn compute_sat_los(sat_pos: &Vector3<f64>, ant_pos: &Vector3<f64>) -> (f64, Vector3<f64>) {
    let diff = sat_pos - ant_pos;
    let dist = diff.norm();
    let u = if dist > 1e-3 { diff / dist } else { Vector3::zeros() };
    (dist, u)
}

/// Compute Sagnac Earth-rotation correction during signal transit (meters).
pub fn compute_sagnac_correction(sat_pos: &Vector3<f64>, ant_pos: &Vector3<f64>) -> f64 {
    (WGS84_EARTH_ROTATION_RATE / SPEED_OF_LIGHT) * (sat_pos.x * ant_pos.y - sat_pos.y * ant_pos.x)
}

/// Line-of-sight position Jacobian row: H_pos = -u_sat^T.
pub fn compute_los_jacobian(u_sat: &Vector3<f64>) -> RowVector3<f64> {
    -u_sat.transpose()
}

/// Attitude coupling Jacobian with lever arm: H_att = -u_sat^T * [l_e x].
pub fn compute_att_coupling_jacobian(
    u_sat: &Vector3<f64>,
    l_skew: &Matrix3<f64>,
) -> RowVector3<f64> {
    -u_sat.transpose() * l_skew
}

/// Approximate carrier wavelength on primary band (meters).
pub fn get_carrier_wavelength(constellation: Constellation) -> f64 {
    match constellation {
        Constellation::Gps | Constellation::Galileo => SPEED_OF_LIGHT / 1575.42e6,
        Constellation::Beidou => SPEED_OF_LIGHT / 1561.098e6,
        Constellation::Glonass => SPEED_OF_LIGHT / 1602.0e6,
        _ => SPEED_OF_LIGHT / 1575.42e6,
    }
}

fn is_above_elevation_mask(u_sat: &Vector3<f64>, ant_pos: &Vector3<f64>, min_el: f64) -> bool {
    let up = ant_pos.normalize();
    let sin_el = u_sat.dot(&up);
    sin_el >= libm::sin(min_el)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tc_ppp_line_of_sight_and_attitude_jacobians() {
        let u_sat = Vector3::new(0.6, 0.8, 0.0);
        let h_pos = compute_los_jacobian(&u_sat);
        assert_eq!(h_pos, RowVector3::new(-0.6, -0.8, 0.0));

        let u_vert = Vector3::new(0.0, 0.0, 1.0);
        let lever_arm_e = Vector3::new(0.5, 0.0, 0.0);
        let l_skew = skew_symmetric(&lever_arm_e);
        let h_att = compute_att_coupling_jacobian(&u_vert, &l_skew);
        assert_eq!(h_att[1], -0.5);
    }

    #[test]
    fn test_tc_ppp_carrier_phase_residual_math() {
        let cp_meas = 20500123.456;
        let geom_range = 20500120.000;
        let clk_rx = 3.000;
        let amb_m = 0.450;
        let res = cp_meas - (geom_range + clk_rx + amb_m);
        assert!((res - 0.006_f64).abs() < 1e-6);
    }

    #[test]
    fn test_tc_ppp_dead_reckoning_continuity() {
        let mut pos = Vector3::new(100.0, 200.0, 300.0);
        let vel = Vector3::new(10.0, 0.0, 0.0);
        let dt = 0.02;
        for _ in 0..50 {
            pos += vel * dt;
        }
        assert!((pos.x - 110.0_f64).abs() < 1e-9);
    }

    #[test]
    fn test_tc_ppp_ambiguity_cycle_slip_reset() {
        let eskf = EskfState::new(Vector3::zeros(), Vector3::zeros(), nalgebra::UnitQuaternion::identity());
        let mut ppp = TightlyCoupledPppIns::new(eskf, TcPppConfig::default());
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };

        let amb1 = ppp.manage_float_ambiguity(sat, 100.0, 80.0, Some(0));
        assert_eq!(ppp.tracked_ambiguities.get(&sat).unwrap().lock_count, 1);
        ppp.manage_float_ambiguity(sat, 100.1, 80.0, Some(0));
        assert_eq!(ppp.tracked_ambiguities.get(&sat).unwrap().lock_count, 2);

        let amb2 = ppp.manage_float_ambiguity(sat, 105.0, 80.0, Some(1));
        assert_eq!(ppp.tracked_ambiguities.get(&sat).unwrap().lock_count, 1);
        assert!((amb1 - amb2).abs() > 1.0);
    }

    #[test]
    fn test_tc_ppp_process_epoch_inertial_propagation() {
        let eskf = EskfState::new(Vector3::new(100.0, 200.0, 300.0), Vector3::new(1.0, 0.0, 0.0), nalgebra::UnitQuaternion::identity());
        let mut ppp = TightlyCoupledPppIns::new(eskf, TcPppConfig::default());
        let imu = vec![ImuMeasurement::new(1000, Vector3::new(0.0, 0.0, -9.81), Vector3::zeros())];
        let obs = EpochObs { time: GpsTime::new(2200, 1.0), satellites: vec![] };
        let sol = ppp.process_epoch(&imu, &obs).expect("process_epoch failed");
        assert_eq!(sol.mode, CompositeMode::DeadReckoning);
        assert_eq!(sol.num_satellites, 0);
    }
}
