//! Double-Difference Iterated Extended Kalman Filter (DD-IEKF) and RTS Smoother Engine.

pub mod ar;
pub mod ar_gate;
pub mod clk_datum;
pub mod formation;
pub mod formation_cov;
pub mod iono_free;
pub mod mw;
pub mod predict;
pub mod ref_sat;
pub mod sat_pco;
pub mod sat_pos;
pub mod satpos;
pub mod screen;
pub mod smoother;
pub mod state;
pub mod update;
pub mod widelane;

use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::{Matrix3, Vector3};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;
use gneiss_parsers::receiver_antenna::ReceiverAntenna;

use crate::post_process::combiner::SmoothedEpoch;
use crate::post_process::iekf_pass::FilteredEpoch;

pub use ar::{resolve_ambiguities, ArResult};
pub use sat_pos::broadcast_position_for;
pub use smoother::{run_rts_smoother, IekfSnapshot};
pub use state::{DoubleDiffKey, RtkState};
pub use update::{iekf_update, DoubleDiffMeasurement};

pub(crate) use ar_gate::seed_ambiguity_variance_cycles2;
#[cfg(test)]
pub(crate) use clk_datum::centered_pair_correction;

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
    /// Per-DD-key precise-clock datum (metres): the value of
    /// `c·(dt_sat − dt_ref)` in force when that arc's ambiguity was last
    /// (re)seeded or reference-switched. Maintained by [`clk_datum`];
    /// stays empty (inert) unless a precise-clock product is loaded.
    ///
    /// [`clk_datum`]: self::clk_datum
    pub(crate) clk_datum_m: HashMap<DoubleDiffKey, f64>,
    /// Reference-switch transfers handled by THIS filter instance
    /// (diagnostic; one stderr line per event under GNEISS_CLK_TRACE).
    pub clk_ref_switches: usize,
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
    /// Kinematic profile flag for ambiguity resolution and process noise.
    pub is_kinematic: bool,
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
    /// When true, fixed ambiguities constrain the filter state (fix-and-hold).
    pub fix_and_hold: bool,
    /// Rover antenna heading offset relative to True North (radians).
    pub rover_heading_rad: f64,
    /// Number of consecutive epochs with consistent fixed ambiguities.
    pub consecutive_fixes: u32,
    last_fixed_ambs: Vec<(DoubleDiffKey, f64)>,
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
            clk_datum_m: HashMap::new(),
            clk_ref_switches: 0,
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
            is_kinematic: false,
            zwd_est_m: 0.0,
            zwd_var_m2: update::ZWD_INIT_VAR_M2,
            prev_zwd_tow: start_time.tow,
            pw_tracker: mw::WidelaneTracker::default(),

            wl_tracker: mw::WidelaneTracker::default(),
            receiver_pcv: None,
            ar_gate: std::env::var("GNEISS_AR_GATE").is_ok_and(|v| v == "1"),
            ar_elevation_mask_rad: Self::DEFAULT_MIN_ELEVATION_RAD,
            fix_and_hold: std::env::var("GNEISS_FIX_AND_HOLD").is_ok_and(|v| v == "1"),
            rover_heading_rad: 0.0,
            consecutive_fixes: 0,
            last_fixed_ambs: Vec::new(),
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
        let f_mat = predict::predict_state_gated(&mut self.state, rover.time, q_now, true);
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
                        self.state.reset_ambiguity(&key, seed, seed_ambiguity_variance_cycles2(m.pr_var_m2, m.lambda));
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

        // Hard prefit gross-error screen, evaluated at the just-predicted
        // (pre-update) position so a blunder can't hide behind an update it
        // corrupted itself -- complements iekf_update_gated's Huber-style
        // soft weighting, which has no defense against a kilometre-scale
        // blunder dragging the linearization point before it can be seen
        // as an outlier. See screen::screen_gross_pr_errors.
        let gross_rejected = screen::screen_gross_pr_errors(&mut meas, self.state.pos_ecef);
        if tracing::enabled!(tracing::Level::DEBUG) {
            for key in &gross_rejected {
                tracing::debug!("gross-error: tow={:.0} rejected sat={} band={}",
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



        let mut ar_res = self.resolve_ar_candidate(&dd_meas);
        let (pos_ecef, cov_pos) = self.finalize_fixed_position(&mut ar_res, &dd_meas);
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

    fn resolve_ar_candidate(&mut self, dd_meas: &formation::DdMeasurements) -> ar::ArResult {
        if self.min_ar_lock_epochs > 0 {
            let min_ep = self.min_ar_lock_epochs;
            let eligible: Vec<DoubleDiffKey> = self.pair_epochs.iter()
                .filter(|(_, &age)| age >= min_ep)
                .map(|(k, _)| *k)
                .collect();
            self.state.retain_active_ambiguities(&eligible);
        }
        let ar_view_owned;
        let ar_view = if self.ar_gate && self.ar_elevation_mask_rad > self.min_elevation_rad {
            ar_view_owned = ar_gate::elevation_filtered_view(
                &self.state, &dd_meas.dd, self.ar_elevation_mask_rad,
            );
            &ar_view_owned
        } else {
            &self.state
        };
        let mut ar_res = ar::resolve_ambiguities(ar_view, 3, self.target_pf, self.is_kinematic);
        if self.widelane_ar {
            let far_vetoed = ar_res.is_fixed
                && !widelane::far_matches_widelanes(&self.wl_tracker, &ar_res);
            if !ar_res.is_fixed || far_vetoed {
                ar_res = ar::float_result(&self.state);
                if let Some(cascade) = widelane::resolve_cascade(&self.state, &self.wl_tracker) {
                    ar_res = cascade;
                }
            }
        }
        self.screen_fixed_residuals(&mut ar_res, &dd_meas.dd);
        ar_res
    }

    fn screen_fixed_residuals(&self, ar_res: &mut ar::ArResult, dd: &[DoubleDiffMeasurement]) {
        if !ar_res.is_fixed {
            return;
        }
        let max_carrier_res = if self.is_kinematic { 0.08 } else { 0.05 };
        let (max_code_rms, max_code_res) = if self.is_kinematic { (6.0, 18.0) } else { (4.0, 12.0) };
        let carrier_ok = update::validate_fixed_carrier_residuals(
            ar_res.position_ecef,
            dd,
            &ar_res.fixed_ambiguities,
            max_carrier_res,
        );
        let pr_ok = update::validate_fixed_pseudorange_residuals(
            ar_res.position_ecef,
            dd,
            max_code_rms,
            max_code_res,
        );
        if !carrier_ok || !pr_ok {
            *ar_res = ar::float_result(&self.state);
        }
    }

    fn finalize_fixed_position(
        &mut self,
        ar_res: &mut ar::ArResult,
        dd_meas: &formation::DdMeasurements,
    ) -> (Vector3<f64>, Matrix3<f64>) {
        if !ar_res.is_fixed {
            ar::reset_fix_hysteresis(&mut self.consecutive_fixes, &mut self.last_fixed_ambs);
            return (ar_res.position_ecef, ar_res.cov_position);
        }
        ar::update_fix_hysteresis(&mut self.consecutive_fixes, &mut self.last_fixed_ambs, &ar_res.fixed_ambiguities);
        if self.fix_and_hold
            && ar_res.ratio >= 3.0
            && self.consecutive_fixes >= 3
            && !ar::condition_state_on_integers(&mut self.state, &ar_res.fixed_ambiguities)
        {
            *ar_res = ar::float_result(&self.state);
            ar::reset_fix_hysteresis(&mut self.consecutive_fixes, &mut self.last_fixed_ambs);
            return (ar_res.position_ecef, ar_res.cov_position);
        }
        match iono_free::apply_fixed_iono_free(&self.state, &dd_meas.iono_free, ar_res) {
            iono_free::IonoFreeOutcome::Solution(pos, cov) => (pos, cov),
            _ => (ar_res.position_ecef, ar_res.cov_position),
        }
    }

    /// Run RTS backward smoothing on all processed epochs.
    pub fn smooth(&self) -> Vec<SmoothedEpoch> {
        smoother::run_rts_smoother(&self.history)
    }

    /// Constellations that participate in DD formation, sorted by id.
    pub fn select_constellations(
        sat_info: &[(gneiss_core::sat::SatelliteId, Vector3<f64>)],
        glo: bool,
    ) -> Vec<u8> {
        ref_sat::select_constellations(sat_info, glo)
    }

    /// Seeds or resets an ambiguity using [`seed_ambiguity_variance_cycles2`]
    /// for its initial uncertainty rather than a fixed constant -- see that
    /// function's doc comment.
    fn update_dd_ambiguity(&mut self, key: DoubleDiffKey, dd_cp: Option<f64>, dd_pr: f64, lambda: f64, lli_slip: bool, pr_var_m2: f64) {
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
                self.state.reset_ambiguity(&key, init_amb, seed_ambiguity_variance_cycles2(pr_var_m2, lambda));
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
        self.state.ensure_ambiguity(key, init_amb - offset, seed_ambiguity_variance_cycles2(pr_var_m2, lambda));
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
}

#[cfg(test)]
mod tests;
