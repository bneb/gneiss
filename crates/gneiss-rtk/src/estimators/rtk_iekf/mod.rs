//! Double-Difference Iterated Extended Kalman Filter (DD-IEKF) and RTS Smoother Engine.

pub mod ar;
pub mod iono_free;
pub mod mw;
pub mod predict;
pub mod smoother;
pub mod state;
pub mod update;
pub mod widelane;

use std::collections::HashMap;
use nalgebra::Vector3;

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::time::GpsTime;

use crate::post_process::combiner::SmoothedEpoch;
use crate::post_process::forward::FilteredEpoch;

pub use ar::{resolve_ambiguities, ArResult};
pub use smoother::{run_rts_smoother, IekfSnapshot};
pub use state::{DoubleDiffKey, RtkState};
pub use update::{iekf_update, DoubleDiffMeasurement};

/// Per-epoch double-difference measurement set (per-band + iono-free).
struct DdMeasurements {
    pub dd: Vec<DoubleDiffMeasurement>,
    pub iono_free: Vec<iono_free::IonoFreeMeasurement>,
    pub active_keys: Vec<DoubleDiffKey>,
}

/// Double-Difference IEKF and RTS Smoother engine for RTK/PPK.
pub struct GnssRtkIekf {
    pub state: RtkState,
    pub history: Vec<IekfSnapshot>,
    pub q_accel: f64,
    pub ref_sats: HashMap<u8, u16>,
    pub min_elevation_rad: f64,
    pub target_pf: f64,
    pub slip_detector: crate::post_process::screening::CycleSlipDetector,
    /// Same detector run over the BASE stream: base-side slips are invisible
    /// to rover-only checks yet shift every DD ambiguity identically once
    /// states accumulate across epochs.
    pub base_slip_detector: crate::post_process::screening::CycleSlipDetector,
    pub prev_arcs: HashMap<DoubleDiffKey, u32>,
    /// Opt-in Melbourne–Wübbena wide-lane cascade AR (long baselines).
    /// Default off; enabling changes only epochs the joint FAR/PAR left float.
    pub widelane_ar: bool,
    /// TOW of the first epoch this filter processed (session anchor for
    /// two-phase static Q timing). For the backward pass this is the END
    /// of the session, so elapsed time counts backward identically.
    pub start_tow: f64,
    /// Two-phase static Q (long-baseline mode): after this many seconds of
    /// elapsed session time, switch from `q_accel` to `static_lock_q_accel`.
    /// Loose early Q lets the float solution converge; tight late Q freezes
    /// the monument so wrong fixes surface as residual inflation instead of
    /// absorbing into position drift. Active only when `widelane_ar` is on.
    pub static_lock_after_s: Option<f64>,
    /// Post-lock acceleration process noise (m/s^2).
    pub static_lock_q_accel: f64,
    /// Rover ZWD residual (m zenith wet) scalar random-walk estimate.
    pub zwd_est_m: f64,
    pub zwd_var_m2: f64,
    prev_zwd_tow: f64,
    /// Phase-only wide-lane arc means per DD pair (long-baseline mode):
    /// pw = dd_phi1 - dd_phi2 - geom_wl, converges to N_W + satellite UPD
    /// difference with millimetre-level noise (no code term involved).
    pub pw_tracker: mw::WidelaneTracker,
    /// Rover ZWD residual (m of zenith wet delay) on top of the Saastamoinen
    /// model, tracked as a scalar random walk. Long-baseline mode only.

    pub wl_tracker: mw::WidelaneTracker,
}

impl GnssRtkIekf {
    /// Create a new RTK IEKF solver.
    pub fn new(initial_pos: Vector3<f64>, start_time: GpsTime, q_accel: f64) -> Self {
        Self {
            state: RtkState::new(initial_pos, start_time),
            history: Vec::new(),
            q_accel,
            ref_sats: HashMap::new(),
            min_elevation_rad: 0.1745, // 10 degrees
            target_pf: 0.001,
            slip_detector: crate::post_process::screening::CycleSlipDetector::new(),
            base_slip_detector: crate::post_process::screening::CycleSlipDetector::new(),
            prev_arcs: HashMap::new(),
            widelane_ar: false,
            start_tow: start_time.tow,
            static_lock_after_s: None,
            static_lock_q_accel: 1e-9,
            zwd_est_m: 0.0,
            zwd_var_m2: update::ZWD_INIT_VAR_M2,
            prev_zwd_tow: start_time.tow,
            pw_tracker: mw::WidelaneTracker::default(),

            wl_tracker: mw::WidelaneTracker::default(),
        }
    }

    /// Process a synchronous rover/base RTK epoch.
    pub fn process_epoch(
        &mut self,
        rover: &EpochObs,
        base: &EpochObs,
        base_pos: Vector3<f64>,
        ephems: &[Ephemeris],
    ) -> Result<FilteredEpoch, String> {
        self.slip_detector.check_epoch(rover);
        // Base-side slip tracking is part of the accumulating-ambiguity
        // machinery; the legacy path must stay untouched.
        if self.widelane_ar {
            self.base_slip_detector.check_epoch(base);
        }
        let dd_meas = self.build_dd_measurements(rover, base, base_pos, ephems)?;
        if dd_meas.dd.is_empty() {
            return Err("No valid double-difference measurements formed".to_string());
        }

        self.state.retain_active_ambiguities(&dd_meas.active_keys);
        if self.widelane_ar {
            self.wl_tracker.retain_active(&dd_meas.active_keys);
        }

        // Reverse-safe covariance growth on the opt-in long-baseline path:
        // the legacy signed-dt Q poisons the backward pass over long arcs.
        let q_now = match self.static_lock_after_s {
            Some(lock_s)
                if self.widelane_ar
                    && (rover.time.tow - self.start_tow).abs() > lock_s =>
            {
                self.static_lock_q_accel
            }
            _ => self.q_accel,
        };
        let f_mat = if self.widelane_ar {
            predict::predict_state_gated(&mut self.state, rover.time, q_now, true)
        } else {
            predict::predict_state(&mut self.state, rover.time, q_now)
        };
        let (x_pred, p_pred) = (self.state.to_dvector(), self.state.cov.clone());

        // Innovation gating: an undetected cycle slip (e.g. a base receiver
        // counter reset invisible to rover-side GF checks) enters the update
        // as a huge innovation against tight priors and destroys the state.
        // Re-seed offending pairs and drop their phase data for this epoch.
        let mut meas = dd_meas.dd.clone();
        if self.widelane_ar {
            let outliers = update::phase_innovation_outliers(
                &self.state, &meas, update::PHASE_INNOVATION_GATE_CYCLES,
            );
            for key in outliers {
                if let Some(m) = dd_meas.dd.iter().find(|m| m.key == key) {
                    if let Some(cp) = m.dd_cp_cycles {
                        let seed = cp - m.dd_pr_m / m.lambda;
                        self.state.reset_ambiguity(&key, seed, 100.0);
                    }
                }
                self.wl_tracker.reset_pair(&DoubleDiffKey { freq_band: 1, ..key });
                for m in meas.iter_mut() {
                    if m.key == key {
                        m.dd_cp_cycles = None;
                    }
                }
                tracing::debug!("slip-gate: tow={:.0} re-seeded sat={} band={}",
                    rover.time.tow, key.sat, key.freq_band);
            }
        }

        update::iekf_update(&mut self.state, &meas)?;
        let (x_post, p_post) = (self.state.to_dvector(), self.state.cov.clone());

        // Long-baseline mode: track the rover ZWD residual from post-update
        // phase innovations with per-epoch step saturation.
        if self.widelane_ar {
            let dt = rover.time.tow - self.prev_zwd_tow;
            let pairs = self.zwd_innovation_pairs(&dd_meas.dd);
            let (z, v) = update::update_zwd_scalar(
                self.zwd_est_m, self.zwd_var_m2,
                update::ZWD_RW_M2_PER_S, dt, &pairs,
            );
            self.zwd_est_m = z;
            self.zwd_var_m2 = v;
            self.prev_zwd_tow = rover.time.tow;
        }



        let mut ar_res = ar::resolve_ambiguities(&self.state, 3, self.target_pf);
        if self.widelane_ar {
            // A FAR fix contradicting a converged MW wide lane is a
            // confidently-wrong fix (slow iono drift dragged the per-band
            // floats past half-cycle; the ratio test cannot see it). Reject
            // it and let the iono-immune cascade try instead.
            let far_vetoed = ar_res.is_fixed
                && !widelane::far_matches_widelanes(&self.wl_tracker, &ar_res);
            if !ar_res.is_fixed || far_vetoed {
                ar_res = ar::float_result(&self.state);
                if let Some(cascade) = widelane::resolve_cascade(&self.state, &self.wl_tracker) {
                    ar_res = cascade;
                }
            }
        }

        // When fixed, re-estimate the position from iono-free phase to
        // remove the DD ionosphere bias that grows with baseline length.
        let (pos_ecef, cov_pos) = if ar_res.is_fixed {
            match iono_free::apply_fixed_iono_free(&self.state, &dd_meas.iono_free, &ar_res) {
                iono_free::IonoFreeOutcome::Solution(pos, cov) => (pos, cov),
                _ => (ar_res.position_ecef, ar_res.cov_position),
            }
        } else {
            (ar_res.position_ecef, ar_res.cov_position)
        };
        let q_flag = if ar_res.is_fixed { 1 } else { 2 };

        self.history.push(IekfSnapshot {
            time: rover.time,
            x_pred,
            p_pred,
            x_post,
            p_post,
            f_mat,
            is_fixed: ar_res.is_fixed,
            n_sats: rover.satellites.len(),
            quality: q_flag,
        });

        Ok(FilteredEpoch {
            time: rover.time,
            position_ecef: pos_ecef,
            velocity_ecef: Some(self.state.vel_ecef),
            attitude: None,
            cov_position: cov_pos,
            n_satellites: rover.satellites.len(),
            quality: q_flag,
            is_fixed: ar_res.is_fixed,
        })
    }

    /// Run RTS backward smoothing on all processed epochs.
    pub fn smooth(&self) -> Vec<SmoothedEpoch> {
        smoother::run_rts_smoother(&self.history)
    }

    /// Form double-difference measurements across all visible satellites.
    fn build_dd_measurements(
        &mut self,
        rover: &EpochObs,
        base: &EpochObs,
        base_pos: Vector3<f64>,
        ephems: &[Ephemeris],
    ) -> Result<DdMeasurements, String> {
        let sat_info = extract_sat_positions(rover, ephems, self.state.pos_ecef, self.min_elevation_rad);
        let mut meas_list = Vec::new();
        let mut if_meas = Vec::new();
        let mut active_keys = Vec::new();

        let mut constellations: Vec<u8> = sat_info.iter().map(|(s, _)| s.constellation as u8).collect();
        constellations.sort_unstable();
        constellations.dedup();

        for const_id in constellations {
            // Skip FDMA GLONASS until receiver-specific inter-channel biases are calibrated
            if const_id == gneiss_core::sat::Constellation::Glonass as u8 {
                continue;
            }

            let const_sats: Vec<(gneiss_core::sat::SatelliteId, Vector3<f64>)> = sat_info
                .iter()
                .filter(|(s, _)| s.constellation as u8 == const_id)
                .cloned()
                .collect();

            if const_sats.len() < 2 {
                continue;
            }

            let ref_sat_id = self.select_reference_satellite_hys(const_id, &const_sats);
            if ref_sat_id == 0 {
                continue;
            }

            let ref_pos = match const_sats.iter().find(|(s, _)| s.prn as u16 == ref_sat_id) {
                Some((_, p)) => *p,
                None => continue,
            };

            let ref_rov = rover.satellites.iter().find(|s| s.sat.constellation as u8 == const_id && s.sat.prn as u16 == ref_sat_id);
            let ref_bas = base.satellites.iter().find(|s| s.sat.constellation as u8 == const_id && s.sat.prn as u16 == ref_sat_id);
            let (r_rov, r_bas) = match (ref_rov, ref_bas) {
                (Some(r), Some(b)) => (r, b),
                _ => continue,
            };

            for (sat_id, sat_pos) in &const_sats {
                if sat_id.prn as u16 == ref_sat_id {
                    continue;
                }
                let rov_s = rover.satellites.iter().find(|s| s.sat == *sat_id);
                let bas_s = base.satellites.iter().find(|s| s.sat == *sat_id);
                if let (Some(rs), Some(bs)) = (rov_s, bas_s) {
                    if self.widelane_ar {
                        mw::update_tracker_from_obs(
                            &mut self.wl_tracker, *sat_id, ref_sat_id, rs, bs, r_rov, r_bas,
                        );
                    }
                    let mut pair_cp: HashMap<u8, f64> = HashMap::new();
                    for freq_band in [1, 2, 5, 7] {
                        if let Some(m) = self.build_single_dd_pair(
                            *sat_id, ref_sat_id, *sat_pos, ref_pos, base_pos, rs, bs, r_rov, r_bas, freq_band,
                        ) {
                            if let Some(cp) = m.dd_cp_cycles {
                                pair_cp.insert(freq_band, cp);
                            }
                            active_keys.push(m.key);
                            meas_list.push(m);
                        }
                    }
                    if self.widelane_ar {
                        // Phase-only wide-lane arc mean (no code term):
                        // converges to N_W + satellite UPD difference with
                        // millimetre-level noise once geometry is removed.
                        if let (Some(cp1), Some(cp2)) =
                            (pair_cp.get(&1), pair_cp.get(&2))
                        {
                            self.update_phase_wl(
                                *sat_id, *sat_pos, ref_sat_id, ref_pos,
                                base_pos, *cp1, *cp2,
                            );
                        }
                    }
                    if let Some(m) = iono_free::form_iono_free_dd(
                        *sat_id, rs, bs, r_rov, r_bas, *sat_pos, ref_pos, base_pos,
                        self.state.pos_ecef,
                        DoubleDiffKey { constellation_id: sat_id.constellation as u8, sat: sat_id.prn as u16, ref_sat: ref_sat_id, freq_band: 1 },
                    ) {
                        if_meas.push(m);
                    }
                }
            }
        }

        Ok(DdMeasurements { dd: meas_list, iono_free: if_meas, active_keys })
    }

    #[allow(clippy::too_many_arguments)]
    fn build_single_dd_pair(
        &mut self,
        sat_id: gneiss_core::sat::SatelliteId,
        ref_sat_id: u16,
        sat_pos: Vector3<f64>,
        ref_pos: Vector3<f64>,
        base_pos: Vector3<f64>,
        rov_s: &SatObs,
        bas_s: &SatObs,
        rov_ref: &SatObs,
        bas_ref: &SatObs,
        freq_band: u8,
    ) -> Option<DoubleDiffMeasurement> {
        let pr_rs = rov_s.get_observable(freq_band)?;
        let pr_rr = rov_ref.get_observable(freq_band)?;
        let pr_bs = bas_s.get_observable(freq_band)?;
        let pr_br = bas_ref.get_observable(freq_band)?;
        let dd_pr = (pr_rs - pr_rr) - (pr_bs - pr_br);

        let freq_hz = gneiss_core::signal::get_frequency(sat_id, freq_band, 0);
        let lambda = SPEED_OF_LIGHT_M_S / freq_hz;
        let (cp_rs, cp_rr, cp_bs, cp_br) = (
            rov_s.get_observable_phase(freq_band),
            rov_ref.get_observable_phase(freq_band),
            bas_s.get_observable_phase(freq_band),
            bas_ref.get_observable_phase(freq_band),
        );

        let dd_cp = match (cp_rs, cp_rr, cp_bs, cp_br) {
            (Some(s1), Some(r1), Some(s2), Some(r2)) => Some((s1 - r1) - (s2 - r2)),
            _ => None,
        };

        let key = DoubleDiffKey {
            constellation_id: sat_id.constellation as u8,
            sat: sat_id.prn as u16,
            ref_sat: ref_sat_id,
            freq_band,
        };

        let (pr_var, cp_var) = self.compute_dd_variances(sat_pos, ref_pos, lambda);
        // Wet-mapping difference at the rover: sensitivity of this DD to the
        // rover ZWD residual state (long-baseline mode).
        let dm_wet_rov = {
            let llh = gneiss_core::coords::ecef_to_llh(self.state.pos_ecef);
            let el_sat = gneiss_core::coords::az_el(llh, self.state.pos_ecef, sat_pos).1;
            let el_ref = gneiss_core::coords::az_el(llh, self.state.pos_ecef, ref_pos).1;
            let (_, w_sat) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(llh, el_sat, self.state.time);
            let (_, w_ref) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(llh, el_ref, self.state.time);
            w_sat - w_ref
        };
        let ref_sat_struct = gneiss_core::sat::SatelliteId {
            constellation: sat_id.constellation,
            prn: ref_sat_id as u8,
        };
        let mut cur_arc = self.slip_detector.get_arc(sat_id)
            + self.slip_detector.get_arc(ref_sat_struct);
        if self.widelane_ar {
            cur_arc += self.base_slip_detector.get_arc(sat_id)
                + self.base_slip_detector.get_arc(ref_sat_struct);
        }
        let arc_changed = match self.prev_arcs.get(&key) {
            Some(&prev) => prev != cur_arc,
            None => false,
        };
        self.prev_arcs.insert(key, cur_arc);

        let lli_slip = rov_s.get_lli(freq_band).is_some_and(|l| (l & 1) != 0)
            || rov_ref.get_lli(freq_band).is_some_and(|l| (l & 1) != 0)
            || arc_changed;

        self.update_dd_ambiguity(key, dd_cp, dd_pr, lambda, lli_slip);

        Some(DoubleDiffMeasurement {
            key,
            dd_pr_m: dd_pr,
            dd_cp_cycles: dd_cp,
            sat_pos,
            ref_pos,
            base_pos,
            lambda,
            pr_var_m2: pr_var,
            cp_var_cycles2: cp_var,
            dm_wet_rov,
        })
    }

    /// Phase innovations with sensitivity to the rover ZWD residual.
    fn zwd_innovation_pairs(
        &self,
        meas: &[update::DoubleDiffMeasurement],
    ) -> Vec<(f64, f64, f64)> {
        let dv = self.state.to_dvector();
        let cur = self.state.pos_ecef;
        let mut out = Vec::new();
        for m in meas {
            let Some(cp) = m.dd_cp_cycles else { continue };
            let Some(ai) = self.state.get_amb_idx(&m.key) else { continue };
            let rs = (m.sat_pos - cur).norm();
            let rr = (m.ref_pos - cur).norm();
            let base_dd =
                (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
            let tropo = update::compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur);
            let Some(zi) = self.state.zwd_idx() else { continue };
            let geom = (rs - rr) - base_dd + tropo;
            let pred = (geom + m.dm_wet_rov * dv[zi]) / m.lambda + dv[ai];
            out.push((m.dm_wet_rov / m.lambda, cp - pred, m.cp_var_cycles2.max(1e-4)));
        }
        out
    }

    /// Absorb one phase-only wide-lane observation for a DD pair.
    ///
    /// pw = dd_phi1 - dd_phi2 - geom_wl, where geom_wl uses the current
    /// float position and the Saastamoinen model. The arc mean converges to
    /// N_W + satellite UPD difference with millimetre-level noise.
    #[allow(clippy::too_many_arguments)]
    fn update_phase_wl(
        &mut self,
        sat_id: gneiss_core::sat::SatelliteId,
        sat_pos: Vector3<f64>,
        ref_sat_id: u16,
        ref_pos: Vector3<f64>,
        base_pos: Vector3<f64>,
        dd_cp1: f64,
        dd_cp2: f64,
    ) {
        let key = DoubleDiffKey {
            constellation_id: sat_id.constellation as u8,
            sat: sat_id.prn as u16,
            ref_sat: ref_sat_id,
            freq_band: 1,
        };
        let cur = self.state.pos_ecef;
        let rs = (sat_pos - cur).norm();
        let rr = (ref_pos - cur).norm();
        let base_dd =
            (sat_pos - base_pos).norm() - (ref_pos - base_pos).norm();
        let tropo = update::compute_tropo_dd(sat_pos, ref_pos, base_pos, cur);
        let f1 = gneiss_core::signal::get_frequency(sat_id, 1, 0);
        let f2 = gneiss_core::signal::get_frequency(sat_id, 2, 0);
        let lambda_wl = SPEED_OF_LIGHT_M_S / (f1 - f2);
        let pwl_cycles = dd_cp1 - dd_cp2;
        let pw = pwl_cycles - ((rs - rr) - base_dd + tropo) / lambda_wl;
        self.pw_tracker.update(key, pw, 0.0, false);
    }

    fn compute_dd_variances(&self, sat_pos: Vector3<f64>, ref_pos: Vector3<f64>, lambda: f64) -> (f64, f64) {
        let rx_llh = gneiss_core::coords::ecef_to_llh(self.state.pos_ecef);
        let (_az_s, el_s) = gneiss_core::coords::az_el(rx_llh, self.state.pos_ecef, sat_pos);
        let (_az_r, el_r) = gneiss_core::coords::az_el(rx_llh, self.state.pos_ecef, ref_pos);

        let sin_s = el_s.sin().max(0.1);
        let sin_r = el_r.sin().max(0.1);

        let pr_var = 2.0 * (0.04 / (sin_s * sin_s) + 0.04 / (sin_r * sin_r));
        let cp_base = 0.003 / lambda;
        let cp_var = 2.0 * (cp_base * cp_base / (sin_s * sin_s) + cp_base * cp_base / (sin_r * sin_r));

        (pr_var, cp_var)
    }

    fn update_dd_ambiguity(&mut self, key: DoubleDiffKey, dd_cp: Option<f64>, dd_pr: f64, lambda: f64, lli_slip: bool) {
        let init_amb = dd_cp.map_or(0.0, |cp| cp - dd_pr / lambda);
        if std::env::var("WL_TRACE").is_ok() && init_amb.abs() > 1e5 {
            eprintln!("BAD-SEED tow-file key={:?} init_amb={:.3e}", key, init_amb);
        }
        if self.state.get_amb_idx(&key).is_some() {
            if lli_slip {
                self.state.reset_ambiguity(&key, init_amb, 100.0);
            }
        } else {
            self.state.ensure_ambiguity(key, init_amb, 100.0);
        }
    }

    fn select_reference_satellite_hys(
        &mut self,
        const_id: u8,
        sats: &[(gneiss_core::sat::SatelliteId, Vector3<f64>)],
    ) -> u16 {
        if let Some(&prev) = self.ref_sats.get(&const_id) {
            if sats.iter().any(|(s, _)| s.prn as u16 == prev) {
                return prev;
            }
        }
        let chosen = sats.first().map_or(0, |(s, _)| s.prn as u16);
        if chosen > 0 {
            self.ref_sats.insert(const_id, chosen);
        }
        chosen
    }
}

fn extract_sat_positions(
    rover: &EpochObs,
    ephems: &[Ephemeris],
    rx_pos: Vector3<f64>,
    min_el: f64,
) -> Vec<(gneiss_core::sat::SatelliteId, Vector3<f64>)> {    let mut out = Vec::new();
    let rx_llh = gneiss_core::coords::ecef_to_llh(rx_pos);

    for s in &rover.satellites {
        let eph_opt = ephems
            .iter()
            .filter(|e| e.sat() == s.sat)
            .min_by(|a, b| {
                (a.toe().tow - rover.time.tow)
                    .abs()
                    .partial_cmp(&(b.toe().tow - rover.time.tow).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        if let Some(eph) = eph_opt {
            let sat_p = compute_signal_sat_pos(s, eph, rover.time);
            let (_az, el) = gneiss_core::coords::az_el(rx_llh, rx_pos, sat_p);
            if el >= min_el {
                out.push((s.sat, sat_p));
            }
        }
    }
    out
}

fn compute_signal_sat_pos(s: &SatObs, eph: &Ephemeris, time: gneiss_core::time::GpsTime) -> Vector3<f64> {
    let pr_m = s.get_observable(1).or_else(|| s.get_observable(2)).unwrap_or(20_000_000.0);
    let tau = pr_m / SPEED_OF_LIGHT_M_S;
    let t_tx = gneiss_core::time::GpsTime::new(time.week, time.tow - tau);
    let (_, _, sat_clk_err_rough, _) = eph.position(t_tx);
    let t_tx_true = gneiss_core::time::GpsTime::new(time.week, t_tx.tow - sat_clk_err_rough);
    let (sat_p, _, _, _) = eph.position(t_tx_true);

    let omega_tau = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau;
    let cos_wt = omega_tau.cos();
    let sin_wt = omega_tau.sin();
    Vector3::new(
        sat_p.x * cos_wt + sat_p.y * sin_wt,
        -sat_p.x * sin_wt + sat_p.y * cos_wt,
        sat_p.z,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::generator::{generate_simulation_dataset, SimulationConfig};

    #[test]
    fn test_gnss_rtk_iekf_runs_on_simulated_dataset() {
        let cfg = SimulationConfig {
            duration_s: 10.0,
            ..Default::default()
        };
        let sim = generate_simulation_dataset(&cfg);
        let mut engine = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);

        let mut fixed_count = 0;
        for i in 0..sim.rover_epochs.len() {
            let sol = engine.process_epoch(
                &sim.rover_epochs[i],
                &sim.base_epochs[i],
                cfg.base_ecef,
                &sim.ephemerides,
            );
            assert!(sol.is_ok());
            let s = sol.unwrap();
            if s.is_fixed {
                fixed_count += 1;
                let err = (s.position_ecef - sim.truth_positions[i].1).norm();
                assert!(err < 0.05, "Fixed epoch error should be < 5cm, got {:.4}m", err);
            }
        }

        assert!(fixed_count >= 5, "RTK engine should fix at least 5 epochs");
        let smoothed = engine.smooth();
        assert_eq!(smoothed.len(), 10);
    }
}
