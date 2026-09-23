//! Double-difference measurement formation: turning one epoch's raw rover
//! and base observations into the [`DoubleDiffMeasurement`]s, iono-free
//! combinations, and phase wide-lane updates that `process_epoch` feeds to
//! the filter. Split out of `mod.rs` (CLAUDE.md's 500-line file standard).

use std::collections::HashMap;

use nalgebra::Vector3;

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_parsers::receiver_antenna::{compute_dd_pcv_correction_2d, frequency_code};

use super::formation_cov::compute_dd_variances;
use super::ref_sat::{
    filter_reference_candidates, prn_u16_to_sat, sat_matches_id, sat_to_prn_u16,
    select_constellations, select_ref_sat_with_hysteresis,
};
use super::sat_pos::{extract_sat_positions, glo_freq_num};
use super::update::GRAD_MIN_SIN_EL;
use super::GnssRtkIekf;
use super::{iono_free, mw, update};
use super::{DoubleDiffKey, DoubleDiffMeasurement};

/// All possible RTK frequency bands across GPS, Galileo, BeiDou, and GLONASS.
pub const ALL_RTK_BANDS: [u8; 5] = [1, 2, 5, 6, 7];

/// Canonical frequency bands formed per constellation.
pub fn canonical_bands_for_constellation(c: Constellation) -> &'static [u8] {
    match c {
        Constellation::Gps | Constellation::Qzss => &[1, 2, 5],
        Constellation::Galileo => &[1, 5, 7],
        Constellation::Beidou => &[1, 6, 7],
        Constellation::Glonass => &[1, 2],
        _ => &[1, 2],
    }
}

/// Per-epoch double-difference measurement set (per-band + iono-free).
pub(super) struct DdMeasurements {
    pub dd: Vec<DoubleDiffMeasurement>,
    pub iono_free: Vec<iono_free::IonoFreeMeasurement>,
    pub active_keys: Vec<DoubleDiffKey>,
}

struct DdContext<'a> {
    rover: &'a EpochObs,
    base: &'a EpochObs,
    base_pos: Vector3<f64>,
    ephems: &'a [Ephemeris],
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
        self.code_phase_div.clear();
        let sat_info = extract_sat_positions(
            rover, ephems, self.state.pos_ecef, self.min_elevation_rad, self.precise_orbits.as_ref(),
        );
        let mut ddm = DdMeasurements {
            dd: Vec::new(),
            iono_free: Vec::new(),
            active_keys: Vec::new(),
        };
        let ctx = DdContext { rover, base, base_pos, ephems };

        for const_id in select_constellations(sat_info.as_slice(), self.enable_glonass) {
            self.form_constellation_dd(const_id, &sat_info, &ctx, &mut ddm);
        }

        Ok(ddm)
    }

    fn form_constellation_dd(
        &mut self,
        const_id: u8,
        sat_info: &[(SatelliteId, Vector3<f64>)],
        ctx: &DdContext<'_>,
        ddm: &mut DdMeasurements,
    ) {
        let const_sats = filter_constellation_sats(sat_info, const_id);
        if const_sats.len() < 2 {
            return;
        }

        let old_ref_opt = self.ref_sats.get(&const_id).copied();
        let ref_candidates = filter_reference_candidates(&const_sats, ctx.rover, ctx.base, const_id);
        let ref_sat_id = select_ref_sat_with_hysteresis(
            const_id, &ref_candidates, self.state.pos_ecef, &mut self.ref_sats,
        );
        if ref_sat_id == 0 {
            return;
        }
        if let Some(old_ref) = old_ref_opt {
            if old_ref != ref_sat_id {
                for &freq_band in &ALL_RTK_BANDS {
                    self.state.transfer_reference_satellite(const_id, freq_band, old_ref, ref_sat_id);
                }
            }
        }

        let Some((_, ref_pos)) = const_sats.iter().find(|(s, _)| sat_to_prn_u16(*s) == ref_sat_id).copied() else { return };
        let Some(r_rov) = ctx.rover.satellites.iter().find(|s| sat_matches_id(s.sat, const_id, ref_sat_id)) else { return };
        let Some(r_bas) = ctx.base.satellites.iter().find(|s| sat_matches_id(s.sat, const_id, ref_sat_id)) else { return };

        for &(sat_id, sat_pos) in &const_sats {
            if sat_to_prn_u16(sat_id) != ref_sat_id {
                self.form_pair_dd(
                    sat_id, sat_pos, ref_sat_id, ref_pos, ctx.base_pos,
                    ctx.rover, ctx.base, r_rov, r_bas, ctx.ephems, const_id, ddm,
                );
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn form_pair_dd(
        &mut self,
        sat_id: SatelliteId,
        sat_pos: Vector3<f64>,
        ref_sat_id: u16,
        ref_pos: Vector3<f64>,
        base_pos: Vector3<f64>,
        rover: &EpochObs,
        base: &EpochObs,
        r_rov: &SatObs,
        r_bas: &SatObs,
        ephems: &[Ephemeris],
        const_id: u8,
        ddm: &mut DdMeasurements,
    ) {
        let (Some(rs), Some(bs)) = (
            rover.satellites.iter().find(|s| s.sat == sat_id),
            base.satellites.iter().find(|s| s.sat == sat_id),
        ) else { return };

        let mw_ok = sat_id.constellation != Constellation::Glonass;
        if self.widelane_ar && mw_ok {
            mw::update_tracker_from_obs(
                &mut self.wl_tracker, sat_id, ref_sat_id, rs, bs, r_rov, r_bas,
                glo_freq_num(ephems, sat_id),
            );
        }

        let mut pair_cp: HashMap<u8, f64> = HashMap::new();
        for &freq_band in canonical_bands_for_constellation(sat_id.constellation) {
            if let Some(m) = self.build_single_dd_pair(
                ephems, sat_id, ref_sat_id, sat_pos, ref_pos, base_pos, rs, bs, r_rov, r_bas, freq_band,
            ) {
                if let Some(cp) = update::pcv_corrected_cp(&m) {
                    pair_cp.insert(freq_band, cp);
                }
                *self.pair_epochs.entry(m.key).or_insert(0) += 1;
                ddm.active_keys.push(m.key);
                ddm.dd.push(m);
            }
        }

        if self.widelane_ar {
            self.handle_widelane_phase_update(
                sat_id, sat_pos, ref_sat_id, ref_pos, base_pos, rs, bs, ephems, &pair_cp,
            );
        }

        let key = DoubleDiffKey { constellation_id: const_id, sat: sat_to_prn_u16(sat_id), ref_sat: ref_sat_id, freq_band: 1 };
        if let Some(m) = iono_free::form_iono_free_dd(
            sat_id, rs, bs, r_rov, r_bas, sat_pos, ref_pos, base_pos,
            self.state.pos_ecef, key, glo_freq_num(ephems, sat_id),
        ) {
            ddm.iono_free.push(m);
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn handle_widelane_phase_update(
        &mut self,
        sat_id: SatelliteId,
        sat_pos: Vector3<f64>,
        ref_sat_id: u16,
        ref_pos: Vector3<f64>,
        base_pos: Vector3<f64>,
        rov_s: &SatObs,
        bas_s: &SatObs,
        ephems: &[Ephemeris],
        pair_cp: &HashMap<u8, f64>,
    ) {
        let b2 = select_secondary_phase_band(rov_s, bas_s);
        if let (Some(&cp1), Some(&cp2)) = (pair_cp.get(&1), pair_cp.get(&b2)) {
            self.update_phase_wl(
                sat_id, sat_pos, ref_sat_id, ref_pos, base_pos, cp1, cp2, b2,
                glo_freq_num(ephems, sat_id),
            );
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn build_single_dd_pair(
        &mut self,
        ephems: &[Ephemeris],
        sat_id: SatelliteId,
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
        let (pr_rs, pr_rr) = (rov_s.get_observable(freq_band)?, rov_ref.get_observable(freq_band)?);
        let (pr_bs, pr_br) = (bas_s.get_observable(freq_band)?, bas_ref.get_observable(freq_band)?);
        let dd_pr = (pr_rs - pr_rr) - (pr_bs - pr_br);

        let glo_k = glo_freq_num(ephems, sat_id);
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

        let const_id = if sat_id.constellation == Constellation::Qzss { 0 } else { sat_id.constellation as u8 };
        let key = DoubleDiffKey { constellation_id: const_id, sat: sat_to_prn_u16(sat_id), ref_sat: ref_sat_id, freq_band };

        let snrs = (rov_s.get_snr(freq_band), rov_ref.get_snr(freq_band), bas_s.get_snr(freq_band), bas_ref.get_snr(freq_band));
        let dd_var = compute_dd_variances(self.state.pos_ecef, sat_pos, ref_pos, lambda, snrs);
        if dd_cp.is_none() && dd_var.pr_var_m2 > 100.0 {
            return None;
        }

        let tide_dd_m = compute_tide_dd(self.state.time, self.state.pos_ecef, base_pos, sat_pos, ref_pos);
        let (dm_wet_rov, dgrad_n_rov, dgrad_e_rov) = compute_tropo_gradients(self.state.pos_ecef, self.state.time, sat_pos, ref_pos);

        let ref_sat_struct = prn_u16_to_sat(const_id, ref_sat_id);
        let lli_slip = check_pair_slip(
            &self.slip_detector, &self.base_slip_detector, &mut self.prev_arcs,
            key, sat_id, ref_sat_struct, rov_s, rov_ref, self.widelane_ar,
        );

        let dd_clk_m = self.formation_clock_corr_m(sat_id, ref_sat_id, self.state.pos_ecef, sat_pos, ref_pos);
        let dd_pr = dd_pr - dd_clk_m;
        let dd_cp = dd_cp.map(|cp| cp - dd_clk_m / lambda);

        self.apply_clock_datum(key, dd_clk_m, lambda, lli_slip);
        let dd_pcv_m = self.receiver_dd_pcv_m(sat_id, freq_band, sat_pos, ref_pos);
        let meas = DoubleDiffMeasurement {
            key, dd_pr_m: dd_pr, dd_cp_cycles: dd_cp, sat_pos, ref_pos, base_pos, lambda,
            pr_var_m2: dd_var.pr_var_m2, cp_var_cycles2: dd_var.cp_var_cycles2,
            pr_ref_var_m2: dd_var.pr_ref_var_m2, cp_ref_var_cycles2: dd_var.cp_ref_var_cycles2,
            dm_wet_rov, dgrad_n_rov, dgrad_e_rov, tide_dd_m, dd_pcv_m,
        };

        self.update_dd_ambiguity(key, update::pcv_corrected_cp(&meas), dd_pr, lambda, lli_slip, meas.pr_var_m2);
        Some(meas)
    }

    /// Differential receiver-PCV carried by one DD pair (metres).
    pub(super) fn receiver_dd_pcv_m(
        &self,
        sat_id: SatelliteId,
        freq_band: u8,
        sat_pos: Vector3<f64>,
        ref_pos: Vector3<f64>,
    ) -> f64 {
        let Some((rov, bas)) = &self.receiver_pcv else { return 0.0 };
        let Some(code) = frequency_code(sat_id.constellation, freq_band) else { return 0.0 };
        let llh = gneiss_core::coords::ecef_to_llh(self.state.pos_ecef);
        let (az_s, el_s) = gneiss_core::coords::az_el(llh, self.state.pos_ecef, sat_pos);
        let (az_r, el_r) = gneiss_core::coords::az_el(llh, self.state.pos_ecef, ref_pos);
        compute_dd_pcv_correction_2d(rov, bas, &code, az_s, el_s, az_r, el_r, self.rover_heading_rad)
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
            let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
            let tropo = update::compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur);
            let Some(zi) = self.state.zwd_idx() else { continue };
            let geom = (rs - rr) - base_dd + tropo;
            let pred = (geom + m.dm_wet_rov * dv[zi]) / m.lambda + dv[ai];
            out.push((m.dm_wet_rov / m.lambda, cp - pred, m.cp_var_cycles2.max(1e-4)));
        }
        out
    }

    #[allow(clippy::too_many_arguments)]
    fn update_phase_wl(
        &mut self,
        sat_id: SatelliteId,
        sat_pos: Vector3<f64>,
        ref_sat_id: u16,
        ref_pos: Vector3<f64>,
        base_pos: Vector3<f64>,
        dd_cp1: f64,
        dd_cp2: f64,
        b2: u8,
        glo_k: i8,
    ) {
        let const_id = if sat_id.constellation == Constellation::Qzss { 0 } else { sat_id.constellation as u8 };
        let key = DoubleDiffKey { constellation_id: const_id, sat: sat_to_prn_u16(sat_id), ref_sat: ref_sat_id, freq_band: 1 };
        let cur = self.state.pos_ecef;
        let rs = (sat_pos - cur).norm();
        let rr = (ref_pos - cur).norm();
        let base_dd = (sat_pos - base_pos).norm() - (ref_pos - base_pos).norm();
        let tropo = update::compute_tropo_dd(sat_pos, ref_pos, base_pos, cur);
        let f1 = gneiss_core::frequencies::track_c_frequency(sat_id.constellation, 1, glo_k);
        let f2 = gneiss_core::frequencies::track_c_frequency(sat_id.constellation, b2, glo_k);
        let lambda_wl = SPEED_OF_LIGHT_M_S / (f1 - f2);
        let pwl_cycles = dd_cp1 - dd_cp2;
        let pw = pwl_cycles - ((rs - rr) - base_dd + tropo) / lambda_wl;
        self.pw_tracker.update(key, pw, 0.0, false);
    }
}

fn filter_constellation_sats(
    sat_info: &[(SatelliteId, Vector3<f64>)],
    const_id: u8,
) -> Vec<(SatelliteId, Vector3<f64>)> {
    sat_info.iter().filter(|(s, _)| {
        if const_id == Constellation::Gps as u8 {
            s.constellation == Constellation::Gps || s.constellation == Constellation::Qzss
        } else {
            s.constellation as u8 == const_id
        }
    }).copied().collect()
}

fn select_secondary_phase_band(rov_s: &SatObs, bas_s: &SatObs) -> u8 {
    if rov_s.get_observable_phase(2).is_some() && bas_s.get_observable_phase(2).is_some() {
        2
    } else if rov_s.get_observable_phase(7).is_some() && bas_s.get_observable_phase(7).is_some() {
        7
    } else if rov_s.get_observable_phase(6).is_some() && bas_s.get_observable_phase(6).is_some() {
        6
    } else {
        5
    }
}

fn compute_tide_dd(
    time: gneiss_core::time::GpsTime,
    rov_pos: Vector3<f64>,
    base_pos: Vector3<f64>,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
) -> f64 {
    let tide_rov = gneiss_core::tides::solid_earth_tides_ecef(time, rov_pos);
    let tide_base = gneiss_core::tides::solid_earth_tides_ecef(time, base_pos);
    let los_rs = (sat_pos - rov_pos).normalize();
    let los_rr = (ref_pos - rov_pos).normalize();
    let los_bs = (sat_pos - base_pos).normalize();
    let los_br = (ref_pos - base_pos).normalize();
    let proj = |tide: Vector3<f64>, los: Vector3<f64>| tide.dot(&los);
    (proj(tide_rov, los_rs) - proj(tide_rov, los_rr)) - (proj(tide_base, los_bs) - proj(tide_base, los_br))
}

fn compute_tropo_gradients(
    pos_ecef: Vector3<f64>,
    time: gneiss_core::time::GpsTime,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
) -> (f64, f64, f64) {
    let llh = gneiss_core::coords::ecef_to_llh(pos_ecef);
    let (az_sat, el_sat) = gneiss_core::coords::az_el(llh, pos_ecef, sat_pos);
    let (az_ref, el_ref) = gneiss_core::coords::az_el(llh, pos_ecef, ref_pos);
    let (_, w_sat) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(llh, el_sat, time);
    let (_, w_ref) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(llh, el_ref, time);
    let grad_term = |az: f64, el: f64| -> (f64, f64) {
        let se = el.sin().max(GRAD_MIN_SIN_EL);
        let m = el.cos().max(0.0) / se;
        (m * az.cos(), m * az.sin())
    };
    let (gn_s, ge_s) = grad_term(az_sat, el_sat);
    let (gn_r, ge_r) = grad_term(az_ref, el_ref);
    (w_sat - w_ref, gn_s - gn_r, ge_s - ge_r)
}

#[allow(clippy::too_many_arguments)]
fn check_pair_slip(
    slip_detector: &crate::post_process::screening::CycleSlipDetector,
    base_slip_detector: &crate::post_process::screening::CycleSlipDetector,
    prev_arcs: &mut HashMap<DoubleDiffKey, u32>,
    key: DoubleDiffKey,
    sat_id: SatelliteId,
    ref_sat: SatelliteId,
    rov_s: &SatObs,
    rov_ref: &SatObs,
    widelane_ar: bool,
) -> bool {
    let mut cur_arc = slip_detector.get_arc(sat_id) + slip_detector.get_arc(ref_sat);
    if widelane_ar {
        cur_arc += base_slip_detector.get_arc(sat_id) + base_slip_detector.get_arc(ref_sat);
    }
    let arc_changed = prev_arcs.insert(key, cur_arc).is_some_and(|prev| prev != cur_arc);
    rov_s.get_lli(key.freq_band).is_some_and(|l| (l & 1) != 0)
        || rov_ref.get_lli(key.freq_band).is_some_and(|l| (l & 1) != 0)
        || arc_changed
}
