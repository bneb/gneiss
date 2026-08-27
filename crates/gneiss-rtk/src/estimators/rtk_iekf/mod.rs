//! Double-Difference Iterated Extended Kalman Filter (DD-IEKF) and RTS Smoother Engine.

pub mod ar;
pub mod ar_gate;
pub mod iono_free;
pub mod mw;
pub mod predict;
pub mod sat_pco;
pub mod satpos;
pub mod smoother;
pub mod state;
pub mod update;
pub mod widelane;

use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::Vector3;

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::time::GpsTime;
use gneiss_parsers::receiver_antenna::{compute_dd_pcv_correction, frequency_code, ReceiverAntenna};

use crate::post_process::combiner::SmoothedEpoch;
use crate::post_process::forward::FilteredEpoch;

pub use ar::{resolve_ambiguities, ArResult};
pub use smoother::{run_rts_smoother, IekfSnapshot};
pub use state::{DoubleDiffKey, RtkState};
pub use update::{iekf_update, DoubleDiffMeasurement};
use update::GRAD_MIN_SIN_EL;

/// Per-epoch double-difference measurement set (per-band + iono-free).
struct DdMeasurements {
    pub dd: Vec<DoubleDiffMeasurement>,
    pub iono_free: Vec<iono_free::IonoFreeMeasurement>,
    pub active_keys: Vec<DoubleDiffKey>,
}

/// Track C signal-registry frequency lookup (same policy as mw.rs).
fn track_c_freq_mod(
    c: gneiss_core::sat::Constellation,
    primary: bool,
    glo_k: i8,
    band: u8,
) -> f64 {
    use gneiss_core::frequencies::{frequency_for, signal_for_band};
    match signal_for_band(c, if primary { 1 } else { band }) {
        Some(sig) => frequency_for(c, sig, glo_k),
        None => {
            if primary { 1_575_420_000.0 } else { 1_227_600_000.0 }
        }
    }
}

/// Pure pair decision for the centered precise-clock DD correction.
///
/// Takes the centered clock lookups for `(satellite, reference)` and
/// returns `(correction_metres, gate_tripped)`:
///
/// - both lookups healthy and neither spread exceeds
///   [`gneiss_parsers::clk_centering::MAX_CENTERED_SPREAD_S`]:
///   `c · (dt_sat − dt_ref)` on the CENTERED biases;
/// - either lookup missing (no record / too few constellation mates):
///   `(0.0, false)` — silently disabled, matching legacy missing-bias
///   behaviour;
/// - either spread over threshold: `(0.0, true)` — pathological product
///   epoch, correction suppressed and the caller latches the warning.
fn centered_pair_correction(
    sat: Option<gneiss_parsers::clk_centering::CenteredClock>,
    reference: Option<gneiss_parsers::clk_centering::CenteredClock>,
) -> (f64, bool) {
    use gneiss_parsers::clk_centering::MAX_CENTERED_SPREAD_S;
    const C: f64 = SPEED_OF_LIGHT_M_S;
    match (sat, reference) {
        (Some(a), Some(b)) => {
            let tripped = a.spread_s > MAX_CENTERED_SPREAD_S || b.spread_s > MAX_CENTERED_SPREAD_S;
            if tripped {
                (0.0, true)
            } else {
                (C * (a.bias_s - b.bias_s), false)
            }
        }
        _ => (0.0, false),
    }
}

/// `GNEISS_CLK_TRACE` diagnostic: first 10 evaluations show the centered
/// biases, spreads, and resulting pair correction (metres).
fn clk_centering_trace(
    tow: f64,
    sat_id: gneiss_core::sat::SatelliteId,
    ref_sv: gneiss_core::sat::SatelliteId,
    cs: Option<gneiss_parsers::clk_centering::CenteredClock>,
    cr: Option<gneiss_parsers::clk_centering::CenteredClock>,
    corr_m: f64,
) {
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    if n >= 10 {
        return;
    }
    let fmt = |c: Option<gneiss_parsers::clk_centering::CenteredClock>| match c {
        Some(x) => format!("{:+.1}us/spr{:.1}us", x.bias_s * 1e6, x.spread_s * 1e6),
        None => "none".to_string(),
    };
    eprintln!(
        "CLKDD[{n}] tow={:.0} {:?}{:02}-{:02} sat={} ref={} -> {:+.3} m",
        tow,
        sat_id.constellation,
        sat_id.prn,
        ref_sv.prn,
        fmt(cs),
        fmt(cr),
        corr_m
    );
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
    /// When set, IekfSnapshot.amb_keys records the DD key list at each
    /// epoch, enabling offline ambiguity-trajectory analysis. Zero-cost
    /// when off (empty Vec).
    pub track_ambiguity_keys: bool,
    /// Minimum epochs a DD pair must be tracked before its ambiguity
    /// participates in AR. Prevents freshly-risen satellites with poorly
    /// converged floats from corrupting LAMBDA. 0 = disabled.
    pub min_ar_lock_epochs: u32,
    /// Per-pair tracking duration (epochs since first seen).
    pair_epochs: HashMap<DoubleDiffKey, u32>,
    /// Precise orbit store (SP3). When set, satellite positions and clocks
    /// come from IGS final/rapid products instead of broadcast ephemerides.
    /// Eliminates ~1-2 m orbit error that only partially cancels in DD
    /// at >20 km baselines.
    pub precise_orbits: Option<std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>>,
    /// Precise clock products (RINEX CLK) paired with SP3 orbits.
    pub precise_clocks: Option<std::sync::Arc<gneiss_parsers::rinex_clk::RinexClock>>,
    /// Latched once the centered-spread gate first suppresses the
    /// precise-clock DD correction: the one-shot warning has been printed.
    clk_gate_warned: std::sync::atomic::AtomicBool,
    /// Opt-in FDMA GLONASS phase participation (own reference satellite and
    /// per-satellite ambiguities absorb phase inter-channel biases). MW
    /// wide-lane stays GPS/Galileo-only: code inter-channel biases do not
    /// cancel between receivers and would corrupt the arcs.
    pub enable_glonass: bool,
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
    /// Multiplier on the robust-innovation gate (Huber knee) inside the
    /// measurement update. 1.0 = legacy static behaviour; the kinematic
    /// profile widens it so coherent model mismatch during acceleration
    /// does not deweight genuine measurements.
    pub robust_innov_scale: f64,
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
    /// Receiver antenna PCV calibrations `(rover, base)`. When both are set,
    /// the elevation-dependent differential receiver PCV is removed from
    /// every DD phase observation before ambiguity estimation. Same-family
    /// pairs cancel to ~zero; cross-family pairs carry mm-level signatures.
    pub receiver_pcv: Option<(Arc<ReceiverAntenna>, Arc<ReceiverAntenna>)>,
    /// Master gate for the opt-in AR-quality techniques in [`ar_gate`]
    /// (AR elevation mask + phase-code coherency bias init). Enabled by env
    /// `GNEISS_AR_GATE=1`; default off preserves legacy behaviour exactly.
    pub ar_gate: bool,
    /// Elevation cut-off (rad) for AR participation when `ar_gate` is on.
    /// Defaults to `min_elevation_rad`, i.e. no additional restriction.
    pub ar_elevation_mask_rad: f64,
    /// This epoch's code-minus-phase divergences `(freq_band, cycles)` of
    /// the pairs tracked so far, used to coherently seed newly initialised
    /// ambiguities when `ar_gate` is on. Cleared each epoch.
    code_phase_div: Vec<(u8, f64)>,
}

impl GnssRtkIekf {
    /// Default measurement elevation cut-off: 10 degrees.
    const DEFAULT_MIN_ELEVATION_RAD: f64 = 0.1745;

    /// Create a new RTK IEKF solver.
    pub fn new(initial_pos: Vector3<f64>, start_time: GpsTime, q_accel: f64) -> Self {
        Self {
            state: RtkState::new(initial_pos, start_time),
            history: Vec::new(),
            q_accel,
            ref_sats: HashMap::new(),
            min_elevation_rad: Self::DEFAULT_MIN_ELEVATION_RAD,
            target_pf: 0.001,
            slip_detector: crate::post_process::screening::CycleSlipDetector::new(),
            base_slip_detector: crate::post_process::screening::CycleSlipDetector::new(),
            prev_arcs: HashMap::new(),
            track_ambiguity_keys: false,
            min_ar_lock_epochs: 0,
            pair_epochs: HashMap::new(),
            precise_orbits: None,
            precise_clocks: None,
            clk_gate_warned: std::sync::atomic::AtomicBool::new(false),
            enable_glonass: false,
            widelane_ar: false,
            start_tow: start_time.tow,
            static_lock_after_s: None,
            static_lock_q_accel: 1e-9,
            robust_innov_scale: 1.0,
            zwd_est_m: 0.0,
            zwd_var_m2: update::ZWD_INIT_VAR_M2,
            prev_zwd_tow: start_time.tow,
            pw_tracker: mw::WidelaneTracker::default(),

            wl_tracker: mw::WidelaneTracker::default(),
            receiver_pcv: None,
            ar_gate: std::env::var("GNEISS_AR_GATE").is_ok_and(|v| v == "1"),
            ar_elevation_mask_rad: Self::DEFAULT_MIN_ELEVATION_RAD,
            code_phase_div: Vec::new(),
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
                    if let Some(cp) = update::pcv_corrected_cp(m) {
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

        update::iekf_update_gated(&mut self.state, &meas, self.robust_innov_scale)?;
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



        // AR eligibility gating: exclude pairs tracked for too few epochs.
        if self.min_ar_lock_epochs > 0 {
            let min_ep = self.min_ar_lock_epochs;
            let eligible: Vec<DoubleDiffKey> = self.pair_epochs.iter()
                .filter(|(_, &age)| age >= min_ep)
                .map(|(k, _)| *k)
                .collect();
            self.state.retain_active_ambiguities(&eligible);
        }
        // Opt-in AR elevation mask (GNEISS_AR_GATE): resolve integers only
        // from pairs whose both members sit above the higher AR cut-off.
        // Non-destructive: only the LAMBDA input is restricted, the live
        // state keeps every float so re-risen pairs resume where they left.
        let ar_view_owned;
        let ar_view = if self.ar_gate && self.ar_elevation_mask_rad > self.min_elevation_rad {
            ar_view_owned = ar_gate::elevation_filtered_view(
                &self.state, &dd_meas.dd, self.ar_elevation_mask_rad,
            );
            &ar_view_owned
        } else {
            &self.state
        };
        let mut ar_res = ar::resolve_ambiguities(ar_view, 3, self.target_pf);
        if self.widelane_ar {
            // A FAR fix contradicting a converged MW wide lane is a
            // confidently-wrong fix (slow iono drift dragged the per-band
            // floats past half-cycle; the ratio test cannot see it). Reject
            // it and let the iono-immune cascade try instead.
            let far_vetoed = ar_res.is_fixed
                && !widelane::far_matches_widelanes(&self.wl_tracker, &ar_res);
            if !ar_res.is_fixed || far_vetoed {
                if tracing::enabled!(tracing::Level::DEBUG) {
                    tracing::debug!(
                        "ar-decision: tow={:.0} far_fixed={} vetoed={} cascade={}",
                        rover.time.tow, ar_res.is_fixed, far_vetoed,
                        widelane::resolve_cascade(&self.state, &self.wl_tracker)
                            .map(|c| c.is_fixed).unwrap_or(false),
                    );
                }
                ar_res = ar::float_result(&self.state);
                if let Some(cascade) = widelane::resolve_cascade(&self.state, &self.wl_tracker) {
                    ar_res = cascade;
                }
            }
        }

        // Experimental: post-fix iono-free residual screen. A same-cycle
        // dual-frequency slip leaves MW (geometry-free) untouched but
        // shifts the pair's IF ambiguity by whole lambda_IF cycles, so
        // deviating pairs expose confidently-wrong fixes that carry no
        // wide-lane contradiction signal.
        if ar_res.is_fixed && std::env::var("GNEISS_IF_VETO").is_ok() {
            let mut fixed_n1: HashMap<DoubleDiffKey, f64> = HashMap::new();
            let mut fixed_n2: HashMap<DoubleDiffKey, f64> = HashMap::new();
            for (k, v) in &ar_res.fixed_ambiguities {
                if k.freq_band == 1 { fixed_n1.insert(*k, *v); }
                if k.freq_band == 2 { fixed_n2.insert(DoubleDiffKey { freq_band: 1, ..*k }, *v); }
            }
            let suspects = update::if_residual_outliers(
                ar_res.position_ecef, &dd_meas.dd, &fixed_n1, &fixed_n2,
            );
            if !suspects.is_empty() {
                tracing::debug!("if-screen: tow={:.0} suspect pairs={} -> float",
                    rover.time.tow, suspects.len());
                ar_res = ar::float_result(&self.state);
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
            amb_keys: if self.track_ambiguity_keys {
                self.state.ambiguities.iter().map(|(k, _)| *k).collect()
            } else {
                Vec::new()
            },
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

    /// Constellations that participate in DD formation, sorted by id.
    /// GLONASS requires the opt-in FDMA policy (see `enable_glonass`).
    fn select_constellations(
        sat_info: &[(gneiss_core::sat::SatelliteId, Vector3<f64>)],
        glo: bool,
    ) -> Vec<u8> {
        let mut v: Vec<u8> = sat_info.iter().map(|(s, _)| s.constellation as u8).collect();
        v.sort_unstable();
        v.dedup();
        let glo_id = gneiss_core::sat::Constellation::Glonass as u8;
        v.retain(|&c| c != glo_id || glo);
        v
    }

    /// Form double-difference measurements across all visible satellites.
    fn build_dd_measurements(
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
        let freq_hz = match gneiss_core::frequencies::signal_for_band(
            sat_id.constellation,
            freq_band,
        ) {
            Some(sig) => gneiss_core::frequencies::frequency_for(sat_id.constellation, sig, glo_k),
            None => gneiss_core::signal::get_frequency(sat_id, freq_band, glo_k),
        };
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

        self.update_dd_ambiguity(key, update::pcv_corrected_cp(&meas), dd_pr, lambda, lli_slip);

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
    fn receiver_dd_pcv_m(
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
    fn formation_clock_corr_m(
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
        let (corr_m, tripped) = centered_pair_correction(cs, cr);
        if tripped {
            self.latch_clk_gate_warning(sat_id, ref_sv);
        }
        if std::env::var("GNEISS_CLK_TRACE").is_ok() {
            clk_centering_trace(self.state.time.tow, sat_id, ref_sv, cs, cr, corr_m);
        }
        corr_m
    }

    /// Print the one-shot spread-gate warning (latched per engine instance
    /// via an atomic flag, so concurrent epochs cannot double-print).
    fn latch_clk_gate_warning(
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
    fn zwd_innovation_pairs(
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
        let f1 = track_c_freq_mod(sat_id.constellation, true, glo_k, 1);
        let f2 = track_c_freq_mod(sat_id.constellation, false, glo_k, b2);
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
        if let Some(abs_idx) = self.state.get_amb_idx(&key) {
            // Established pair: record its code-minus-phase divergence
            // (raw seed minus converged float) as this epoch's coherency
            // sample before any slip reset.
            let rel = abs_idx - self.state.amb_offset();
            let divergence = init_amb - self.state.ambiguities[rel].1;
            self.record_code_phase_divergence(key.freq_band, dd_cp.is_some(), divergence);
            if lli_slip {
                self.state.reset_ambiguity(&key, init_amb, 100.0);
            }
            return;
        }
        // First epoch for this pair (absent from the state == never
        // initialised; re-risen pairs were dropped by retain and count as
        // new again). Under the gate, de-bias the raw code-minus-phase seed
        // by the median divergence of the pairs already tracked this epoch,
        // so the pair's code multipath/iono residual does not enter the
        // filter as a step in the phase innovation.
        let coherent = self.ar_gate && dd_cp.is_some();
        let offset = if coherent {
            ar_gate::coherence_offset(key.freq_band, &self.code_phase_div)
        } else {
            0.0
        };
        if coherent {
            // The seeded pair joins this epoch's sample set so later new
            // pairs on the same band see a stable median.
            self.code_phase_div.push((key.freq_band, offset));
        }
        self.state.ensure_ambiguity(key, init_amb - offset, 100.0);
        self.state.ensure_iono(key, 4.0); // 2m sigma iono residual (legacy)
        if self.state.sat_iono_enabled {
            self.state.ensure_sat_iono_key(key.constellation_id, key.sat);
            self.state.ensure_sat_iono_key(key.constellation_id, key.ref_sat);
        }
    }

    /// Record one per-epoch coherency sample when the gate is active.
    fn record_code_phase_divergence(&mut self, band: u8, has_phase: bool, divergence: f64) {
        if self.ar_gate && has_phase {
            self.code_phase_div.push((band, divergence));
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


/// GLONASS FDMA frequency channel number for a satellite, from its
/// broadcast ephemeris (0 when unknown -> nominal frequency).
fn glo_freq_num(ephems: &[Ephemeris], sat_id: gneiss_core::sat::SatelliteId) -> i8 {
    for e in ephems {
        if let Ephemeris::Glonass(g) = e {
            if g.sat == sat_id {
                return g.freq_num;
            }
        }
    }
    0
}

fn extract_sat_positions(
    rover: &EpochObs,
    ephems: &[Ephemeris],
    rx_pos: Vector3<f64>,
    min_el: f64,
    precise_orbits: Option<&std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>>,
) -> Vec<(gneiss_core::sat::SatelliteId, Vector3<f64>)> {
    let mut out = Vec::new();
    let rx_llh = gneiss_core::coords::ecef_to_llh(rx_pos);

    for s in &rover.satellites {
        // Precise orbits take priority when available.
        if let Some(precise) = precise_orbits {
            let sys_char = match s.sat.constellation {
                gneiss_core::sat::Constellation::Gps => 'G',
                gneiss_core::sat::Constellation::Glonass => 'R',
                gneiss_core::sat::Constellation::Galileo => 'E',
                gneiss_core::sat::Constellation::Beidou => 'C',
                _ => 'G',
            };
            let sv_name = format!("{}{:02}", sys_char, s.sat.prn);
            if std::env::var("GNEISS_SP3_PROBE").is_ok() {
                static FIRST: std::sync::atomic::AtomicUsize =
                    std::sync::atomic::AtomicUsize::new(0);
                let n = FIRST.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if n < 4 {
                    let hit = precise.position_at(&sv_name, rover.time).is_some();
                    eprintln!("SP3-PROBE[{}] {} -> {}", n, sv_name, hit);
                }
            }
            if let Some((pos_com, _clk)) =
                precise.position_at_with_hint(&sv_name, rover.time, Some(rx_pos))
            {
                // Transmit-time + Sagnac: rotate into receive-epoch ECEF
                // using signal travel time (mirrors compute_signal_sat_pos).
                let tau = (rx_pos - pos_com).norm()
                    / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
                let wt = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau;
                let (sw, cw) = libm::sincos(wt);
                let rotated = Vector3::new(
                    pos_com.x * cw + pos_com.y * sw,
                    -pos_com.x * sw + pos_com.y * cw,
                    pos_com.z,
                );
                // CoM -> L1 phase centre via nadir projection, after rotation.
                let pos = sat_pco::apply_sat_pco_z(rotated, s.sat.prn as u16);
                let (_az, el) = gneiss_core::coords::az_el(rx_llh, rx_pos, pos);
                if el >= min_el {
                    out.push((s.sat, pos));
                }
                continue;
            }
            // Precise orbit doesn't cover this sat; fall through to broadcast.
        }

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


/// Nearest-TOE broadcast position for a satellite (pipeline delegate).
pub fn broadcast_position_for(
    ephems: &[Ephemeris],
    sv: &gneiss_core::sat::SatelliteId,
    t: GpsTime,
) -> Option<(Vector3<f64>, f64)> {
    let mut best: Option<(&Ephemeris, f64)> = None;
    for cand in ephems {
        if cand.sat().constellation != sv.constellation || cand.sat().prn != sv.prn {
            continue;
        }
        let toe = match cand {
            Ephemeris::Gps(g) => g.toe,
            Ephemeris::Galileo(g) => g.toe,
            Ephemeris::Glonass(_) => return None, // FDMA handled separately
            _ => continue,
        };
        let dt = if toe.week == t.week { (toe.tow - t.tow).abs() } else { f64::INFINITY };
        if best.is_none_or(|(_, d)| dt < d) {
            best = Some((cand, dt));
        }
    }
    let (eph, _) = best?;
    let (p, _, clk, _) = match eph {
        Ephemeris::Gps(g) => g.position(t),
        Ephemeris::Galileo(g) => g.position(t),
        _ => return None,
    };
    Some((p, clk))
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

    fn test_engine(start: GpsTime) -> GnssRtkIekf {
        GnssRtkIekf::new(Vector3::new(1.0, 2.0, 3.0), start, 1.0)
    }

    fn dd_key(sat: u16, band: u8) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: band }
    }

    /// Run the sim dataset through `engine`, returning per-epoch fix flags
    /// and positions for cross-run comparisons.
    fn run_sim(
        engine: &mut GnssRtkIekf,
        sim: &crate::sim::generator::SimulationDataset,
        base: Vector3<f64>,
    ) -> Vec<(bool, Vector3<f64>)> {
        let mut out = Vec::new();
        for i in 0..sim.rover_epochs.len() {
            let s = engine.process_epoch(
                &sim.rover_epochs[i], &sim.base_epochs[i], base, &sim.ephemerides,
            ).expect("epoch must process");
            out.push((s.is_fixed, s.position_ecef));
        }
        out
    }

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
    #[test]
    fn test_select_constellations_drops_glonass_by_default() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        let mk = |c: Constellation, prn: u8| {
            (SatelliteId { constellation: c, prn }, Vector3::zeros())
        };
        let sat_info = vec![
            mk(Constellation::Gps, 6u8),
            mk(Constellation::Glonass, 8u8),
            mk(Constellation::Galileo, 3u8),
            mk(Constellation::Gps, 12u8),
        ];
        // Default policy: FDMA GLONASS stays out until ICB handling lands.
        let got = GnssRtkIekf::select_constellations(&sat_info, false);
        assert_eq!(got, vec![
            Constellation::Gps as u8,
            Constellation::Galileo as u8,
        ]);
        // Opt-in keeps it, sorted and de-duplicated.
        let got = GnssRtkIekf::select_constellations(&sat_info, true);
        assert_eq!(got, vec![
            Constellation::Gps as u8,
            Constellation::Glonass as u8,
            Constellation::Galileo as u8,
        ]);
        // Empty input stays empty.
        assert!(GnssRtkIekf::select_constellations(&[], true).is_empty());
    }

    // ---- GNEISS_AR_GATE: technique 1, AR elevation mask -------------------

    #[test]
    fn test_gate_defaults_off_and_mask_defaults_to_measurement_mask() {
        let eng = test_engine(GpsTime::new(2200, 0.0));
        assert!(!eng.ar_gate, "gate must default off without env");
        assert_eq!(eng.ar_elevation_mask_rad, eng.min_elevation_rad);
        assert!(eng.code_phase_div.is_empty());
    }

    #[test]
    fn test_ar_elevation_mask_unit_view_through_engine() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.ar_gate = true;
        // Realistic rover ECEF so az/el geometry is meaningful.
        eng.state.pos_ecef = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
        let hi = dd_key(2, 1);
        let lo = dd_key(3, 1);
        eng.state.ensure_ambiguity(hi, 10.5, 4.0);
        eng.state.ensure_ambiguity(lo, 20.5, 9.0);
        let up = eng.state.pos_ecef.normalize();
        let east = Vector3::new(-up.z, 0.0, up.x).normalize();
        let mk_meas = |k: DoubleDiffKey, dir: Vector3<f64>| update::DoubleDiffMeasurement {
            key: k,
            dd_pr_m: 0.0,
            dd_cp_cycles: None,
            sat_pos: eng.state.pos_ecef + dir * 2.4e7,
            ref_pos: eng.state.pos_ecef + up * 2.6e7,
            base_pos: eng.state.pos_ecef,
            lambda: 0.19,
            pr_var_m2: 1.0,
            cp_var_cycles2: 1.0,
            dm_wet_rov: 0.0,
            dgrad_n_rov: 0.0,
            dgrad_e_rov: 0.0,
            tide_dd_m: 0.0,
            dd_pcv_m: 0.0,
        };
        let meas = vec![mk_meas(hi, up), mk_meas(lo, east)];
        // 15 deg AR mask: horizon satellite must not reach LAMBDA.
        eng.ar_elevation_mask_rad = 0.2618;
        let view = ar_gate::elevation_filtered_view(&eng.state, &meas, eng.ar_elevation_mask_rad);
        assert_eq!(view.ambiguities.len(), 1);
        assert_eq!(view.ambiguities[0].0, hi);
        // Non-destructive: live state keeps both floats.
        assert_eq!(eng.state.ambiguities.len(), 2);
    }

    #[test]
    fn test_gate_off_extreme_ar_mask_is_ignored_end_to_end() {
        let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
        let sim = generate_simulation_dataset(&cfg);
        // Gate off -> ar_elevation_mask_rad must be inert: fixes still occur.
        let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        eng.ar_elevation_mask_rad = 1.55;
        let res = run_sim(&mut eng, &sim, cfg.base_ecef);
        let fixed = res.iter().filter(|(f, _)| *f).count();
        assert!(fixed >= 5, "legacy path ignores the AR mask, got {fixed} fixes");
    }

    #[test]
    fn test_gate_on_extreme_ar_mask_suppresses_all_fixes() {
        let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
        let sim = generate_simulation_dataset(&cfg);
        // ~89 deg mask excludes every pair from LAMBDA -> float-only output.
        let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        eng.ar_gate = true;
        eng.ar_elevation_mask_rad = 1.55;
        let res = run_sim(&mut eng, &sim, cfg.base_ecef);
        assert!(
            res.iter().all(|(f, _)| !f),
            "AR mask must exclude all pairs; unexpected fix"
        );
    }

    #[test]
    fn test_gate_on_with_default_masks_still_fixes_sim() {
        let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
        let sim = generate_simulation_dataset(&cfg);
        let mut eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        eng.ar_gate = true;
        let mut fixed = 0;
        for i in 0..sim.rover_epochs.len() {
            let s = eng.process_epoch(
                &sim.rover_epochs[i], &sim.base_epochs[i], cfg.base_ecef, &sim.ephemerides,
            ).expect("epoch must process");
            if s.is_fixed {
                fixed += 1;
                let err = (s.position_ecef - sim.truth_positions[i].1).norm();
                assert!(err < 0.05, "gated fixed error {err:.4} m exceeds 5 cm");
            }
        }
        assert!(fixed >= 5, "gate on with defaults should still fix, got {fixed}");
    }

    // ---- ProcessingDynamics: kinematic profile vs static profile --------

    /// Mean 3D tracking error of one engine over the converged tail
    /// (epochs after `skip`) of a simulated trajectory.
    fn mean_tail_error(
        engine: &mut GnssRtkIekf,
        sim: &crate::sim::generator::SimulationDataset,
        base: Vector3<f64>,
        skip: usize,
    ) -> f64 {
        let mut sum = 0.0;
        let mut n = 0;
        for i in 0..sim.rover_epochs.len() {
            let sol = engine
                .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], base, &sim.ephemerides)
                .expect("epoch must process");
            if i >= skip {
                sum += (sol.position_ecef - sim.truth_positions[i].1).norm();
                n += 1;
            }
        }
        sum / n.max(1) as f64
    }

    #[test]
    fn kinematic_profile_tracks_linear_ramp_and_velocity_states_stay_live() {
        use crate::post_process::dynamics::{
            KINEMATIC_INNOV_GATE_SCALE, KINEMATIC_Q_ACCEL, STATIC_Q_ACCEL,
        };
        use crate::sim::generator::{TrajectoryProfile};
        // PPK cadence (30 s epochs): a 5 m/s rover covers 150 m between
        // epochs. Constant-velocity ramp: even the static profile survives
        // here because robust variance inflation lets strong phase
        // measurements drag the frozen prior along; the kinematic profile
        // must also stay bounded AND keep live velocity states.
        let cfg = SimulationConfig {
            duration_s: 1800.0,
            epoch_rate_hz: 1.0 / 30.0,
            profile: TrajectoryProfile::Linear {
                start_offset_ned: Vector3::new(100.0, 100.0, 0.0),
                velocity_ned: Vector3::new(4.0, 3.0, 0.0),
            },
            ..Default::default()
        };
        let sim = generate_simulation_dataset(&cfg);
        let true_speed = Vector3::new(4.0_f64, 3.0, 0.0).norm();

        let mut static_eng =
            GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, STATIC_Q_ACCEL);
        let mut kin_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, KINEMATIC_Q_ACCEL);
        kin_eng.robust_innov_scale = KINEMATIC_INNOV_GATE_SCALE;

        let err_static = mean_tail_error(&mut static_eng, &sim, cfg.base_ecef, 10);
        let err_kin = mean_tail_error(&mut kin_eng, &sim, cfg.base_ecef, 10);
        let speed_kin = kin_eng.state.vel_ecef.norm();
        let speed_static = static_eng.state.vel_ecef.norm();
        eprintln!(
            "KIN-SIM linear ramp @30s: err static={err_static:.3} m kin={err_kin:.3} m | \
             final speed est static={speed_static:.2} kin={speed_kin:.2} (true {true_speed:.2}) m/s"
        );
        // No crash/divergence for either profile on steady motion.
        assert!(err_kin < 3.0, "kinematic must track the ramp, got {err_kin:.3} m");
        assert!(err_static < 3.0, "static+Huber also tracks steady ramps, got {err_static:.3} m");
        // Differential signal: kinematic velocity states stay live
        // (within 25% of the true 5 m/s). NOTE (measured): the STATIC
        // profile's velocity state ALSO converges on a clean constant-
        // velocity ramp — sequential position updates make velocity
        // observable through the F coupling even with monument Q — so
        // no frozenness assertion is possible here; the profiles
        // separate under ACCELERATION (next test).
        assert!(
            (speed_kin - true_speed).abs() < 0.25 * true_speed,
            "kinematic velocity state must track true speed, got {speed_kin:.2}"
        );
    }

    #[test]
    fn kinematic_profile_outperforms_static_under_acceleration() {
        use crate::post_process::dynamics::{
            KINEMATIC_INNOV_GATE_SCALE, KINEMATIC_Q_ACCEL, STATIC_Q_ACCEL,
        };
        use crate::sim::generator::{TrajectoryProfile};
        // Accelerating frame: 15 m/s on a 500 m radius circle is
        // 0.45 m/s^2 of sustained acceleration — a constant-velocity
        // model with monument Q must lag every epoch.
        let cfg = SimulationConfig {
            duration_s: 3600.0,
            epoch_rate_hz: 1.0 / 30.0,
            profile: TrajectoryProfile::Circular {
                center_offset_ned: Vector3::new(200.0, 200.0, 0.0),
                radius_m: 500.0,
                speed_m_s: 15.0,
            },
            ..Default::default()
        };
        let sim = generate_simulation_dataset(&cfg);

        let mut static_eng =
            GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, STATIC_Q_ACCEL);
        let mut kin_eng = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, KINEMATIC_Q_ACCEL);
        kin_eng.robust_innov_scale = KINEMATIC_INNOV_GATE_SCALE;

        let err_static = mean_tail_error(&mut static_eng, &sim, cfg.base_ecef, 12);
        let err_kin = mean_tail_error(&mut kin_eng, &sim, cfg.base_ecef, 12);
        eprintln!(
            "KIN-SIM circular @30s: err static={err_static:.3} m kinematic={err_kin:.3} m"
        );
        assert!(err_kin < 3.0, "kinematic must track the turn, got {err_kin:.3} m");
        assert!(
            err_static > 2.0 * err_kin,
            "static-tuned Q must lag under sustained acceleration \
             (static={err_static:.3}, kin={err_kin:.3})"
        );
    }

    // ---- GNEISS_AR_GATE: technique 2, phase-code coherency bias init ------

    #[test]
    fn test_new_bias_init_gate_off_keeps_raw_code_phase_seed() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        let lambda = 0.2_f64;
        let raw = 105.75 - 100.0 / lambda;
        eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false);
        assert_eq!(eng.state.get_amb_idx(&dd_key(5, 1)), Some(6));
        assert!((eng.state.ambiguities[0].1 - raw).abs() < 1e-12, "seed must equal raw init");
        assert!(eng.code_phase_div.is_empty(), "no samples recorded when gate off");
    }

    #[test]
    fn test_new_bias_init_applies_median_coherency_offset() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.ar_gate = true;
        // Established band-1 pair plus this epoch's divergence samples:
        // median over band 1 is 3.5; the band-2 sample must be ignored.
        eng.state.ensure_ambiguity(dd_key(2, 1), 50.25, 100.0);
        eng.code_phase_div = vec![(1, 2.5), (2, 100.0), (1, 3.5), (1, 4.5)];
        let lambda = 0.2_f64;
        let raw = 105.75 - 100.0 / lambda;
        eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false);
        let seeded = eng.state.ambiguities.iter().find(|(k, _)| *k == dd_key(5, 1)).unwrap().1;
        assert!((seeded - (raw - 3.5)).abs() < 1e-12, "seed = raw − median(band-1 divergences)");
        assert_eq!(eng.code_phase_div.len(), 5, "seeded pair joins the epoch's sample set");
        assert_eq!(eng.code_phase_div.last().copied(), Some((1, 3.5)));
    }

    #[test]
    fn test_new_bias_init_without_prior_samples_seeds_raw_under_gate() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.ar_gate = true;
        let lambda = 0.2_f64;
        let raw = 105.75 - 100.0 / lambda;
        eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false);
        let seeded = eng.state.ambiguities[0].1;
        assert!((seeded - raw).abs() < 1e-12, "no prior pairs -> no offset");
        assert_eq!(eng.code_phase_div, vec![(1, 0.0)]);
    }

    #[test]
    fn test_coherency_offset_never_crosses_frequency_bands() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.ar_gate = true;
        // Only band-2 samples exist; a new band-1 pair must seed raw.
        eng.code_phase_div = vec![(2, 2.5), (2, 3.5)];
        let lambda = 0.2_f64;
        let raw = 105.75 - 100.0 / lambda;
        eng.update_dd_ambiguity(dd_key(5, 1), Some(105.75), 100.0, lambda, false);
        assert!((eng.state.ambiguities[0].1 - raw).abs() < 1e-12);
    }

    #[test]
    fn test_slip_reset_stays_raw_even_under_gate() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.ar_gate = true;
        eng.state.ensure_ambiguity(dd_key(2, 1), 50.25, 100.0);
        eng.code_phase_div = vec![(1, 2.5), (1, 3.5)];
        let lambda = 0.2_f64;
        let raw = 205.75 - 100.0 / lambda;
        // Slip on an ESTABLISHED pair: legacy re-seed semantics preserved.
        eng.update_dd_ambiguity(dd_key(2, 1), Some(205.75), 100.0, lambda, true);
        assert!((eng.state.ambiguities[0].1 - raw).abs() < 1e-12);
        assert_eq!(eng.code_phase_div.len(), 3, "established pair still contributes a sample");
    }

    // ---- GNEISS_AR_GATE: regression, gate off == exact legacy behaviour ---

    #[test]
    fn test_gate_disabled_run_matches_default_run_bit_for_bit() {
        let cfg = SimulationConfig { duration_s: 10.0, ..Default::default() };
        let sim = generate_simulation_dataset(&cfg);
        let mut baseline = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        let mut gated = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);
        gated.ar_gate = false; // explicit, but identical to the default
        let base_out = run_sim(&mut baseline, &sim, cfg.base_ecef);
        let gate_out = run_sim(&mut gated, &sim, cfg.base_ecef);
        assert_eq!(base_out, gate_out, "gate off must reproduce legacy exactly");
        assert!(base_out.iter().filter(|(f, _)| *f).count() >= 5);
    }

    // ---- receiver_dd_pcv_m: correction wiring ---------------------------

    /// Zero unless calibrations are loaded; otherwise equal to the raw
    /// differential PCV at the rover-frame elevations. `self.receiver_pcv`
    /// being `Some` IS the caller's opt-in signal — no separate env gate
    /// (a prior GNEISS_RECV_PCV/GNEISS_PCV split silently required both
    /// to be set for the documented "opt-in via GNEISS_PCV=1" to actually
    /// apply anything; see docs/NETWORK_RTK_NEXT_STEPS.md). Skipped when
    /// igs14 is absent.
    #[test]
    fn receiver_dd_pcv_m_requires_loaded_pair() {
        let Ok(db) = gneiss_parsers::antex::AntexDatabase::parse("../../datasets/igs14.atx")
        else {
            return;
        };
        use gneiss_parsers::receiver_antenna::{compute_dd_pcv_correction, ReceiverAntenna};
        let trm =
            Arc::new(ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").expect("igs14 TRM"));
        let ash =
            Arc::new(ReceiverAntenna::lookup(&db, "ASH701945B_M", "SCIT").expect("igs14 ASH"));
        let mut eng =
            GnssRtkIekf::new(Vector3::new(-3961904.43, 3348994.27, 3698211.71), GpsTime::new(2000, 100.0), 1.0);
        let sat_pos = eng.state.pos_ecef + Vector3::new(1.0e7, 5.0e6, 2.0e7);
        let ref_pos = eng.state.pos_ecef + Vector3::new(0.0, 0.0, 2.4e7);
        let sid = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 3,
        };

        // No calibrations loaded -> no correction.
        eng.receiver_pcv = None;
        assert_eq!(eng.receiver_dd_pcv_m(sid, 1, sat_pos, ref_pos), 0.0);
        // Calibrations loaded -> raw differential PCV (non-zero for this
        // cross-family pair at distinct elevations).
        eng.receiver_pcv = Some((trm.clone(), ash.clone()));
        let dd = eng.receiver_dd_pcv_m(sid, 1, sat_pos, ref_pos);
        let llh = gneiss_core::coords::ecef_to_llh(eng.state.pos_ecef);
        let (_, el_s) = gneiss_core::coords::az_el(llh, eng.state.pos_ecef, sat_pos);
        let (_, el_r) = gneiss_core::coords::az_el(llh, eng.state.pos_ecef, ref_pos);
        let expected = compute_dd_pcv_correction(&trm, &ash, "G01", el_s, el_r);
        assert!(dd.abs() > 1e-6, "cross-family correction must be non-zero: {dd}");
        assert!((dd - expected).abs() < 1e-12, "dd={dd} expected={expected}");
    }

    // ---- GNEISS_CLK: constellation-median centering + spread gate -------

    use gneiss_parsers::clk_centering::CenteredClock;
    use gneiss_parsers::rinex_clk::{ClockRecord, RinexClock};

    fn centered(bias_us: f64, spread_us: f64) -> Option<CenteredClock> {
        Some(CenteredClock { bias_s: bias_us * 1e-6, spread_s: spread_us * 1e-6 })
    }

    #[test]
    fn centered_pair_correction_healthy_pair_applies_centered_delta() {
        let (corr, tripped) =
            centered_pair_correction(centered(10.0, 5.0), centered(-14.0, 6.0));
        assert!(!tripped);
        assert!((corr - SPEED_OF_LIGHT_M_S * 24e-6).abs() < 1e-9);
    }

    #[test]
    fn centered_pair_correction_gates_when_either_side_spread_trips() {
        for (a, b) in [
            (centered(1.0, 150.0), centered(2.0, 5.0)),
            (centered(1.0, 5.0), centered(2.0, 150.0)),
        ] {
            let (corr, tripped) = centered_pair_correction(a, b);
            assert!(tripped, "pathological side must trip the gate");
            assert_eq!(corr, 0.0, "gated correction must be suppressed");
        }
    }

    #[test]
    fn centered_pair_correction_threshold_is_strictly_greater() {
        // Exactly at the 100 us threshold the product is still trusted;
        // a hair above it trips.
        let ok = centered_pair_correction(centered(1.0, 100.0), centered(2.0, 99.999));
        assert!(!ok.1);
        assert!(ok.0.abs() > 0.0);
        let bad = centered_pair_correction(centered(1.0, 100.001), centered(2.0, 5.0));
        assert!(bad.1 && bad.0 == 0.0);
    }

    #[test]
    fn centered_pair_correction_missing_side_stays_silent_zero() {
        assert_eq!(centered_pair_correction(None, centered(2.0, 5.0)), (0.0, false));
        assert_eq!(centered_pair_correction(centered(1.0, 5.0), None), (0.0, false));
        assert_eq!(centered_pair_correction(None, None), (0.0, false));
    }

    /// Synthetic product helper: one record per satellite at tow = 100
    /// (matching the engine time below), biases in microseconds.
    fn clk_product(biases_us: &[(u8, f64)]) -> Arc<RinexClock> {
        let mut rc = RinexClock::default();
        for (prn, us) in biases_us {
            rc.satellites.insert(
                gneiss_core::sat::SatelliteId {
                    constellation: gneiss_core::sat::Constellation::Gps,
                    prn: *prn,
                },
                vec![ClockRecord { time: GpsTime::new(2200, 100.0), bias: us * 1e-6 }],
            );
        }
        Arc::new(rc)
    }

    fn dd_probe(eng: &GnssRtkIekf, sat: u8, reference: u8) -> f64 {
        let rx = Vector3::zeros();
        eng.formation_clock_corr_m(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: sat,
            },
            u16::from(reference),
            rx,
            Vector3::new(2.0e7, 0.0, 0.0),
            Vector3::new(2.4e7, 0.0, 0.0),
        )
    }

    #[test]
    fn precise_clock_dd_m_without_product_is_zero() {
        let eng = test_engine(GpsTime::new(2200, 100.0));
        assert!(eng.precise_clocks.is_none());
        assert_eq!(dd_probe(&eng, 5, 9), 0.0);
    }

    /// Healthy product under a big common mode: the correction equals the
    /// RAW pairwise delta exactly — centering removes only the datum and
    /// must leave the differential untouched through the full wiring.
    #[test]
    fn precise_clock_dd_m_healthy_product_preserves_raw_delta() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.precise_clocks = Some(clk_product(&[
            (27, 500.0),
            (28, 520.0),
            (5, 480.0),
            (10, 510.0),
        ]));
        let corr = dd_probe(&eng, 28, 27);
        let raw = SPEED_OF_LIGHT_M_S * 20e-6; // 520 - 500 us
        assert!((corr - raw).abs() < 1e-9, "corr {corr} vs raw {raw}");
        assert!(!eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed));
    }

    /// Pathological product: correction disabled (0.0) and the warning
    /// flag latches on the first gated call.
    #[test]
    fn precise_clock_dd_m_pathological_product_returns_zero_and_latches() {
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.precise_clocks = Some(clk_product(&[
            (1, 600.0),
            (2, -600.0),
            (3, 590.0),
            (4, -590.0),
        ]));
        assert_eq!(dd_probe(&eng, 1, 3), 0.0);
        assert!(
            eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed),
            "gate trip must latch the warning flag"
        );

        // Fewer than three valid mates: silently zero, no new behaviour.
        let mut eng = test_engine(GpsTime::new(2200, 100.0));
        eng.precise_clocks = Some(clk_product(&[(1, 600.0), (2, -600.0)]));
        assert_eq!(dd_probe(&eng, 1, 2), 0.0);
        assert!(!eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed));
    }

}

#[cfg(test)]
mod obs_side_clk_tests {
    use super::*;

    /// Structural-consistency contract: with a clock product loaded and a
    /// healthy (non-gated) bias pair, the FORMED measurement carries no
    /// clock signature — dd_pr/dd_cp are already corrected at formation.
    /// This is what makes iono_free.rs / ar_gate.rs / filter consumers
    /// consistent without any per-consumer handling.
    #[test]
    fn formation_subtracts_clock_delta_from_measurements() {
        let t = GpsTime::new(2370, 43_200.0);
        let mut eng = GnssRtkIekf::new(Vector3::zeros(), t, 1.0);
        eng.state.iono_enabled = false;
        eng.precise_clocks = Some(std::sync::Arc::new(
            gneiss_parsers::rinex_clk::RinexClock::parse(
                &synthetic_clk_content(),
            ),
        ));
        {
            let content = synthetic_clk_content();
            eprintln!("CONTENT repr: {:?}", content);
            let probe = gneiss_parsers::rinex_clk::RinexClock::parse(&content);
            eprintln!("probe sats: {}", probe.satellites.len());
            let direct = eng.precise_clocks.as_ref().unwrap()
                .centered_clock(sv1_of(), t);
            eprintln!("ENGINE-TEST centered g01 = {direct:?}");
        }
        let corr = eng.formation_clock_corr_m(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            2,
            Vector3::zeros(),
            Vector3::new(2.0e7, 0.0, 0.0),
            Vector3::new(2.1e7, 0.0, 0.0),
        );
        // Synthetic biases ±100 µs → delta must be c·(b1−b2) magnitude,
        // i.e., tens of metres — proving correction is computed formation-side.
        // Centered cluster: g01 sits 20 µs from median → c·20 µs ≈ 6 km.
        assert!(
            (corr.abs() - 5_995.8).abs() < 10.0,
            "correction {corr} m != expected c·(−20 µs) = −5995.8 m"
        );
    }

    fn sv1_of() -> gneiss_core::sat::SatelliteId {
        gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn: 1 }
    }

    fn synthetic_clk_content() -> String {
        // Record epoch = 2025-06-08T12:00 GPST == week 2370, tow 43200
        // (matches the eval instant below).
        // Two GPS sats, biases ∓100 µs at one epoch.
        format!(
            "     3.00           C                                       RINEX VERSION / TYPE\n\
             2    AS    AR                                          # / TYPES OF DATA\n\
             AS G01  2025  6  8 12  0  0.000000  1   -0.000110000000E+00\n\
             AS G02  2025  6  8 12  0  0.000000  1   -0.000090000000E+00\n\
             AS G03  2025  6  8 12  0  0.000000  1   -0.000090000000E+00\n"
        )
    }
}
