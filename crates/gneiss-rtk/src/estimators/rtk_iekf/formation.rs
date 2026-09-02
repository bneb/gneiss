//! Double-difference measurement formation: turning one epoch's raw rover
//! and base observations into the [`DoubleDiffMeasurement`]s, iono-free
//! combinations, and phase wide-lane updates that `process_epoch` feeds to
//! the filter. Split out of `mod.rs` (CLAUDE.md's 500-line file standard)
//! since this is the largest single cohesive concern that lived there —
//! see `docs/PROJECT_STATUS.md` Sprint 13.

use std::collections::HashMap;

use nalgebra::Vector3;

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::time::GpsTime;
use gneiss_parsers::receiver_antenna::{compute_dd_pcv_correction, frequency_code};

use super::sat_pos::{extract_sat_positions, glo_freq_num};
use super::update::GRAD_MIN_SIN_EL;
use super::GnssRtkIekf;
use super::{iono_free, mw, update};
use super::{DoubleDiffKey, DoubleDiffMeasurement};

/// Per-epoch double-difference measurement set (per-band + iono-free).
pub(super) struct DdMeasurements {
    pub dd: Vec<DoubleDiffMeasurement>,
    pub iono_free: Vec<iono_free::IonoFreeMeasurement>,
    pub active_keys: Vec<DoubleDiffKey>,
}

impl GnssRtkIekf {
    /// Form double-difference measurements across all visible satellites.
    pub(super) fn build_dd_measurements(
        &mut self,
        rover: &EpochObs,
        base: &EpochObs,
        base_pos: Vector3<f64>,
        ephems: &[Ephemeris],
    ) -> Result<DdMeasurements, String> {
        // Divergences are a per-epoch sample: rebuilt from scratch each time.
        self.code_phase_div.clear();
        let sat_info = extract_sat_positions(rover, ephems, self.state.pos_ecef, self.min_elevation_rad, self.precise_orbits.as_ref());
        let mut meas_list = Vec::new();
        let mut if_meas = Vec::new();
        let mut active_keys = Vec::new();

        if std::env::var("GNEISS_GLO_DEBUG").is_ok() {
            let glo_in_sat_info = sat_info.iter().filter(|(s, _)| s.constellation == gneiss_core::sat::Constellation::Glonass).count();
            let selected = Self::select_constellations(sat_info.as_slice(), self.enable_glonass);
            eprintln!("GLO-DEBUG enable_glonass={} sat_info_total={} glo_in_sat_info={} selected_constellations={:?}",
                self.enable_glonass, sat_info.len(), glo_in_sat_info, selected);
        }
        for const_id in Self::select_constellations(sat_info.as_slice(), self.enable_glonass) {

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
                    // MW arcs are code-bias sensitive: FDMA inter-channel
                    // biases do not cancel between receivers, so GLONASS
                    // pairs validate but never feed the arcs.
                    let mw_eligible = sat_id.constellation
                        != gneiss_core::sat::Constellation::Glonass;
                    if self.widelane_ar && mw_eligible {
                        mw::update_tracker_from_obs(
                            &mut self.wl_tracker, *sat_id, ref_sat_id, rs, bs, r_rov, r_bas,
                            glo_freq_num(ephems, *sat_id),
                        );
                    }
                    let mut pair_cp: HashMap<u8, f64> = HashMap::new();
                    for freq_band in [1, 2, 5, 7] {
                        if let Some(m) = self.build_single_dd_pair(
                            ephems, *sat_id, ref_sat_id, *sat_pos, ref_pos, base_pos, rs, bs, r_rov, r_bas, freq_band,
                        ) {
                            if let Some(cp) = update::pcv_corrected_cp(&m) {
                                pair_cp.insert(freq_band, cp);
                            }
                            *self.pair_epochs.entry(m.key).or_insert(0) += 1;
                            active_keys.push(m.key);
                            meas_list.push(m);
                        }
                    }
                    if self.widelane_ar {
                        // Phase-only wide-lane arc mean (no code term):
                        // converges to N_W + satellite UPD difference with
                        // millimetre-level noise once geometry is removed.
                        // Secondary band: L2 for GPS/GLONASS; E5a
                        // (band 5) for Galileo exports lacking L2.
                        let b2 = match (rov_s, bas_s) {
                            (Some(r), Some(b)) if
                                r.get_observable_phase(2).is_some()
                                && b.get_observable_phase(2).is_some() => 2,
                            _ => 5,
                        };
                        if let (Some(cp1), Some(cp2)) =
                            (pair_cp.get(&1), pair_cp.get(&b2))
                        {
                            self.update_phase_wl(
                                *sat_id, *sat_pos, ref_sat_id, ref_pos,
                                base_pos, *cp1, *cp2, b2,
                                glo_freq_num(ephems, *sat_id),
                            );
                        }
                    }
                    if let Some(m) = iono_free::form_iono_free_dd(
                        *sat_id, rs, bs, r_rov, r_bas, *sat_pos, ref_pos, base_pos,
                        self.state.pos_ecef,
                        DoubleDiffKey { constellation_id: sat_id.constellation as u8, sat: sat_id.prn as u16, ref_sat: ref_sat_id, freq_band: 1 },
                        glo_freq_num(ephems, *sat_id),
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
        ephems: &[Ephemeris],
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

        let glo_k = glo_freq_num(ephems, sat_id);
        // Authoritative band->signal resolution. Policy functions
        // (secondary_signal etc.) choose WHICH band to prefer; they must
        // never override the frequency of a band actually observed.
        if std::env::var("GNEISS_FREQ_TRACE").is_ok() {
            static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if n < 12 {
                eprintln!(
                    "FREQ[{}] {:?} band={} glo_k={}",
                    n, sat_id.constellation, freq_band, glo_k
                );
            }
        }
        let freq_hz = gneiss_core::frequencies::track_c_frequency(sat_id.constellation, freq_band, glo_k);
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
        // rover ZWD residual state (long-baseline mode), plus the cot(el)
        // gradient mapping differences [north, east] for the same pair.
        // Solid Earth tide correction for this pair (DD of LOS projections).
        let tide_dd_m = {
            let tide_rov = gneiss_core::tides::solid_earth_tides_ecef(self.state.time, self.state.pos_ecef);
            let tide_base = gneiss_core::tides::solid_earth_tides_ecef(self.state.time, base_pos);
            let los_rs = (sat_pos - self.state.pos_ecef).normalize();
            let los_rr = (ref_pos - self.state.pos_ecef).normalize();
            let los_bs = (sat_pos - base_pos).normalize();
            let los_br = (ref_pos - base_pos).normalize();
            let proj = |tide: Vector3<f64>, los: Vector3<f64>| tide.dot(&los);
            (proj(tide_rov, los_rs) - proj(tide_rov, los_rr))
                - (proj(tide_base, los_bs) - proj(tide_base, los_br))
        };
        let (dm_wet_rov, dgrad_n_rov, dgrad_e_rov) = {
            let llh = gneiss_core::coords::ecef_to_llh(self.state.pos_ecef);
            let (az_sat, el_sat) = gneiss_core::coords::az_el(llh, self.state.pos_ecef, sat_pos);
            let (az_ref, el_ref) = gneiss_core::coords::az_el(llh, self.state.pos_ecef, ref_pos);
            let (_, w_sat) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(llh, el_sat, self.state.time);
            let (_, w_ref) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(llh, el_ref, self.state.time);
            let grad_term = |az: f64, el: f64| -> (f64, f64) {
                let se = el.sin().max(GRAD_MIN_SIN_EL);
                let m = el.cos().max(0.0) / se;
                (m * az.cos(), m * az.sin())
            };
            let (gn_s, ge_s) = grad_term(az_sat, el_sat);
            let (gn_r, ge_r) = grad_term(az_ref, el_ref);
            (w_sat - w_ref, gn_s - gn_r, ge_s - ge_r)
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

        // Receiver antenna PCV: the differential signature travels on the
        // measurement and is removed by every consumer through
        // pcv_corrected_cp, including the ambiguity seed below.
        // OBS-SIDE precise-clock correction: subtract c·(dt_sat − dt_ref)
        // from code and phase HERE so every downstream consumer (filter,
        // iono-free re-estimation, AR gate) sees clock-free measurements —
        // structurally consistent by construction (ledger rows 10–11).
        let dd_clk_m = self.formation_clock_corr_m(
            sat_id,
            ref_sat_id,
            self.state.pos_ecef,
            sat_pos,
            ref_pos,
        );
        let dd_pr = dd_pr - dd_clk_m;
        let dd_cp = dd_cp.map(|cp| cp - dd_clk_m / lambda);

        // Clock-datum bookkeeping must precede the ambiguity update so a
        // reference-switch transfer is in place BEFORE innovations form.
        self.apply_clock_datum(key, dd_clk_m, lambda, lli_slip);

        let dd_pcv_m = self.receiver_dd_pcv_m(sat_id, freq_band, sat_pos, ref_pos);
        let meas = DoubleDiffMeasurement {
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
            dgrad_n_rov,
            dgrad_e_rov,
            tide_dd_m,
            dd_pcv_m,
        };

        self.update_dd_ambiguity(key, update::pcv_corrected_cp(&meas), dd_pr, lambda, lli_slip, meas.pr_var_m2);

        Some(meas)
    }

    /// Differential receiver-PCV carried by one DD pair (metres).
    ///
    /// `phi_obs = rho/lam + N + phi_PCV(zenith)`, so the measurement
    /// stores its embedded antenna signature in `dd_pcv_m` and every
    /// consumer removes it via [`update::pcv_corrected_cp`] (`cp -
    /// pcv/lambda`, matching the windup convention). Elevations are
    /// computed in the rover frame and reused for both stations:
    /// baseline << orbit altitude keeps station-to-station elevation
    /// differences sub-mm in PCV terms. Returns 0.0 — a no-op correction
    /// — unless calibrations are loaded (`self.receiver_pcv.is_some()`,
    /// the caller's own opt-in signal) and the frequency is known.
    pub(super) fn receiver_dd_pcv_m(
        &self,
        sat_id: gneiss_core::sat::SatelliteId,
        freq_band: u8,
        sat_pos: Vector3<f64>,
        ref_pos: Vector3<f64>,
    ) -> f64 {
        let Some((rov, bas)) = &self.receiver_pcv else {
            return 0.0;
        };
        let Some(code) = frequency_code(sat_id.constellation, freq_band) else {
            return 0.0;
        };
        let llh = gneiss_core::coords::ecef_to_llh(self.state.pos_ecef);
        let (_az_s, el_s) = gneiss_core::coords::az_el(llh, self.state.pos_ecef, sat_pos);
        let (_az_r, el_r) = gneiss_core::coords::az_el(llh, self.state.pos_ecef, ref_pos);
        let corr = compute_dd_pcv_correction(rov, bas, &code, el_s, el_r);
        if std::env::var("GNEISS_PCV_DEBUG").is_ok() && corr.abs() > 1e-12 {
            eprintln!("PCV [{}]: {:.4} mm (el_s={:.1} el_r={:.1})", sat_id, corr*1000.0, el_s.to_degrees(), el_r.to_degrees());
        }
        corr
    }


    /// Differential precise-satellite-clock correction (metres of range),
    /// constellation-median centered with spread gating.
    ///
    /// Each satellite's clock bias is evaluated at its approximate transmit
    /// time and CENTERED by the median across its constellation mates at
    /// that instant (removes the product's arbitrary timescale datum), then
    /// the pair correction `c · (dt_sat − dt_ref)` is formed from the
    /// centered biases. If the centered inter-satellite spread exceeds 100
    /// µs the product epoch is pathological: the pair correction is
    /// suppressed (0.0) and a latched warning prints once per engine.
    /// Returns 0.0 when no product is loaded or either lookup is missing.
    pub(super) fn formation_clock_corr_m(
        &self,
        sat_id: gneiss_core::sat::SatelliteId,
        ref_sat_id: u16,
        rx_pos: Vector3<f64>,
        sat_pos: Vector3<f64>,
        ref_pos: Vector3<f64>,
    ) -> f64 {
        use gneiss_core::constants::SPEED_OF_LIGHT_M_S as C;
        let Some(clk_prod) = &self.precise_clocks else {
            return 0.0;
        };
        let tau_s = (rx_pos - sat_pos).norm() / C;
        let tau_r = (rx_pos - ref_pos).norm() / C;
        let t_s = GpsTime::new(self.state.time.week, self.state.time.tow - tau_s);
        let t_r = GpsTime::new(self.state.time.week, self.state.time.tow - tau_r);
        let ref_sv = gneiss_core::sat::SatelliteId {
            constellation: sat_id.constellation,
            prn: ref_sat_id as u8,
        };
        let cs = clk_prod.centered_clock(sat_id, t_s);
        let cr = clk_prod.centered_clock(ref_sv, t_r);
        let (corr_m, tripped) = super::centered_pair_correction(cs, cr);
        if tripped {
            self.latch_clk_gate_warning(sat_id, ref_sv);
        }
        if std::env::var("GNEISS_CLK_TRACE").is_ok() {
            super::clk_centering_trace(self.state.time.tow, sat_id, ref_sv, cs, cr, corr_m);
        }
        corr_m
    }

    /// Print the one-shot spread-gate warning (latched per engine instance
    /// via an atomic flag, so concurrent epochs cannot double-print).
    pub(super) fn latch_clk_gate_warning(
        &self,
        sat_id: gneiss_core::sat::SatelliteId,
        ref_sv: gneiss_core::sat::SatelliteId,
    ) {
        use std::sync::atomic::Ordering;
        if self
            .clk_gate_warned
            .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok()
        {
            eprintln!(
                "CLK-GATE: precise-clock DD correction disabled: centered \
                 inter-satellite spread exceeded {:.0} us (pathological \
                 clock-product epoch); all further corrections on this \
                 engine are suppressed [first trip: {:?}{:02}-{:02}, tow {:.0}]",
                gneiss_parsers::clk_centering::MAX_CENTERED_SPREAD_S * 1e6,
                sat_id.constellation,
                sat_id.prn,
                ref_sv.prn,
                self.state.time.tow
            );
        }
    }

    /// Phase innovations with sensitivity to the rover ZWD residual.
    pub(super) fn zwd_innovation_pairs(
        &self,
        meas: &[update::DoubleDiffMeasurement],
    ) -> Vec<(f64, f64, f64)> {
        let dv = self.state.to_dvector();
        let cur = self.state.pos_ecef;
        let mut out = Vec::new();
        for m in meas {
            let Some(cp) = update::pcv_corrected_cp(m) else { continue };
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
        b2: u8,
        glo_k: i8,
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
        let f1 = gneiss_core::frequencies::track_c_frequency(sat_id.constellation, 1, glo_k);
        let f2 = gneiss_core::frequencies::track_c_frequency(sat_id.constellation, b2, glo_k);
        let lambda_wl = SPEED_OF_LIGHT_M_S / (f1 - f2);
        let pwl_cycles = dd_cp1 - dd_cp2;
        let pw = pwl_cycles - ((rs - rr) - base_dd + tropo) / lambda_wl;
        self.pw_tracker.update(key, pw, 0.0, false);
    }

    pub(super) fn compute_dd_variances(&self, sat_pos: Vector3<f64>, ref_pos: Vector3<f64>, lambda: f64) -> (f64, f64) {
        let rx_llh = gneiss_core::coords::ecef_to_llh(self.state.pos_ecef);
        let (_az_s, el_s) = gneiss_core::coords::az_el(rx_llh, self.state.pos_ecef, sat_pos);
        let (_az_r, el_r) = gneiss_core::coords::az_el(rx_llh, self.state.pos_ecef, ref_pos);

        let sin_s = el_s.sin().max(0.1);
        let sin_r = el_r.sin().max(0.1);
        let pr_base = 0.20;
        let pr_var = 2.0 * (pr_base * pr_base / (sin_s * sin_s) + pr_base * pr_base / (sin_r * sin_r));
        let cp_base = 0.003 / lambda;
        let cp_var = 2.0 * (cp_base * cp_base / (sin_s * sin_s) + cp_base * cp_base / (sin_r * sin_r));

        (pr_var, cp_var)
    }
}
