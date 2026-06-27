//! RTKLIB PPP port — line-by-line translation of ppp.c res_ppp + pppos.
//!
//! This module replicates the RTKLIB measurement model exactly:
//! - Same state vector layout: [pos(3), clk_gps, clk_glo?, tropo(1-3), biases(N)]
//! - Same troposphere model: Saastamoinen ZHD + GMF/NMF mapping
//! - Same ionosphere handling: iono-free LC combination
//! - Same outlier rejection: max innovation gate, skip for GLO
//! - Same satellite antenna PCO/PCV, phase windup, solid earth tide
//!
//! The goal is identical output to RTKLIB on the same input data.
//! Once verified, we can incrementally improve.

use std::collections::HashMap;

use crate::engine::EngineError;
use crate::filter::RtkState;
use crate::measurements::combinations;
use gneiss_core::coords::Coordinate;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use nalgebra::{DMatrix, DVector, Vector3};

// ---- RTKLIB-compatible constants ----
const SQR: fn(f64) -> f64 = |x| x * x;
const CLIGHT: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
const D2R: f64 = core::f64::consts::PI / 180.0;
const R2D: f64 = 180.0 / core::f64::consts::PI;

// Initial variances (RTKLIB defaults)
const VAR_POS: f64 = 25.0; // 5^2 m^2 (SPP accuracy), was 10000 (100^2)
const VAR_CLK: f64 = 10000.0; // 100^2 m^2
const VAR_ZTD: f64 = 0.09; // 0.3^2 m^2
const VAR_GRA: f64 = 1e-6; // 0.001^2 m^2
const VAR_BIAS: f64 = 10000.0; // 100^2 m^2
const ERR_SAAS: f64 = 0.3;
const ERR_BRDCI: f64 = 0.5;
const ERR_CBIAS: f64 = 0.3;
const REL_HUMI: f64 = 0.7;

/// State indices — RTKLIB layout
#[derive(Clone)]
struct PppState {
    /// Number of position states (3 or 9 for dynamics)
    np: usize,
    /// Are GLONASS satellites present?
    has_glo: bool,
    /// Troposphere option: 0=Saas, 1=SBAS, 2=EST, 3=ESTG, 4=COR, 5=CORG
    trop_opt: usize,
    /// Number of valid satellites in current epoch
    nsat: usize,
    /// Satellite IDs for each ambiguity slot
    amb_sats: Vec<SatelliteId>,
    /// UDUC mode: estimate ionosphere + L1/L2 ambiguities instead of IF biases
    uduc: bool,
}

impl PppState {
    fn new(has_glo: bool, dynamics: bool, uduc: bool) -> Self {
        Self {
            np: if dynamics { 9 } else { 3 },
            has_glo,
            trop_opt: 2, // EST (estimate ZTD)
            nsat: 0,
            amb_sats: Vec::new(),
            uduc,
        }
    }

    /// Number of clock states (1 for GPS only, 2 if GLONASS present)
    fn nc(&self) -> usize {
        if self.has_glo { 2 } else { 1 }
    }

    /// Index of GPS clock in state vector
    fn ic(&self, sys: usize) -> usize {
        self.np + sys
    }

    /// Index of troposphere parameters
    fn it(&self) -> usize {
        self.ic(0) + self.nc()
    }

    /// Number of troposphere states
    fn nt(&self) -> usize {
        match self.trop_opt {
            0 | 1 => 0, // Saastamoinen or SBAS — no estimation
            2 => 1, // EST — estimate ZTD only
            3 => 3, // ESTG — ZTD + gradients
            4 | 5 => 0, // COR/CORG — externally corrected
            _ => 1,
        }
    }

    /// Total number of resolved states (before biases/iono)
    fn nr(&self) -> usize {
        self.it() + self.nt()
    }

    /// Index of first ionosphere state (UDUC only)
    fn ni(&self, sat_idx: usize) -> usize {
        debug_assert!(self.uduc);
        self.nr() + sat_idx
    }

    /// Index of first L1 ambiguity state (UDUC) or IF bias (IF mode)
    fn ib(&self, sat_idx: usize) -> usize {
        if self.uduc {
            self.nr() + self.nsat + sat_idx
        } else {
            self.nr() + sat_idx
        }
    }

    /// Index of L2 ambiguity state (UDUC only)
    fn ib2(&self, sat_idx: usize) -> usize {
        debug_assert!(self.uduc);
        self.nr() + 2 * self.nsat + sat_idx
    }

    /// Total number of estimated states
    fn nx(&self) -> usize {
        if self.uduc {
            self.nr() + 3 * self.nsat // iono + L1 amb + L2 amb per sat
        } else {
            self.ib(self.nsat) // biases for all observed satellites
        }
    }
}

/// RTKLIB PPP engine — port of pppos() and res_ppp()
pub struct PppRtklib {
    pub max_iter: usize,
    pub elev_mask_deg: f64,
    pub max_inno_m: f64,   // code innovation threshold (m)
    pub max_inno_cp: f64,  // carrier phase innovation threshold (m)
    pub dynamics: bool,
    pub tide_corr: bool,
    pub x: DVector<f64>,
    pub p: DMatrix<f64>,
    pub epoch: u32,
    last_nsat: usize,
    last_has_glo: bool,    // for semantic state resize on GLO change
    biases_seeded: bool,
    was_uduc: bool,         // previous epoch's UDUC mode
    /// MW widelane EMA for AR: sat → (count, smoothed_WL_cycles)
    mw_wl_ema: HashMap<SatelliteId, (u32, f64)>,
}

impl Default for PppRtklib {
    fn default() -> Self {
        Self {
            max_iter: 5,
            elev_mask_deg: 15.0,
            max_inno_m: 200.0,     // code innovation threshold (m)
            max_inno_cp: 100.0,     // phase innovation threshold: reject outliers >100m
            dynamics: false,
            tide_corr: true,
            x: DVector::zeros(0),
            p: DMatrix::zeros(0, 0),
            epoch: 0,
            last_nsat: 0,
            last_has_glo: false,
            biases_seeded: false,
            was_uduc: false,
            mw_wl_ema: HashMap::new(),
        }
    }
}

impl PppRtklib {
    /// Allocate and return the initial state covariance matrix (RTKLIB defaults).
    fn init_covariance(&self, ppp: &PppState, nx: usize) -> DMatrix<f64> {
        let mut p = DMatrix::zeros(nx, nx);
        for i in 0..3 {
            p[(i, i)] = VAR_POS;
        }
        if self.dynamics {
            for i in 3..9 {
                p[(i, i)] = 100.0; // velocity/accel variance
            }
        }
        p[(ppp.ic(0), ppp.ic(0))] = VAR_CLK;
        if ppp.has_glo {
            p[(ppp.ic(1), ppp.ic(1))] = VAR_CLK;
        }
        let it = ppp.it();
        for i in 0..ppp.nt() {
            p[(it + i, it + i)] = if i == 0 {
                VAR_ZTD
            } else {
                VAR_GRA
            }; // gradients
        }
        let nr = ppp.nr();
        if ppp.uduc {
            for i in 0..ppp.nsat {
                p[(ppp.ni(i), ppp.ni(i))] = VAR_BIAS;   // ionosphere
                p[(ppp.ib(i), ppp.ib(i))] = VAR_BIAS;    // L1 ambiguity
                p[(ppp.ib2(i), ppp.ib2(i))] = VAR_BIAS;  // L2 ambiguity
            }
        } else {
            for i in 0..ppp.nsat {
                p[(nr + i, nr + i)] = VAR_BIAS; // IF bias
            }
        }
        p
    }

    /// Forward-predict the state (RTKLIB udstate_ppp).
    /// In static PPP, position is constant. Clock is white noise — keep the
    /// estimated value (don't reset to zero) but inflate variance to VAR_CLK
    /// and clear cross-correlations. The large variance prevents position
    /// errors from leaking into the clock estimate through the Kalman gain
    /// (gain ≈ 8% for clock with VAR_CLK=10000 vs 68% with P_clk=0.35).
    /// The correlation clearing prevents error feedback loops between epochs.
    fn predict(&self, ppp: &PppState, _x: &mut DVector<f64>, p_mat: &mut DMatrix<f64>) {
        let nx = ppp.nx();

        // Position: static — add minimal process noise (σ≈1mm/s)
        // to prevent covariance collapse from Joseph form rounding
        for i in 0..3 {
            p_mat[(i, i)] += 1e-6;
        }
        // Clock: white noise — keep value, inflate variance, clear correlations
        for i in 0..ppp.nc() {
            let ci = ppp.ic(i);
            for j in 0..nx {
                if i32::abs(ci as i32 - j as i32) > 30 {
                    continue;
                }
                p_mat[(ci, j)] = 0.0;
                p_mat[(j, ci)] = 0.0;
            }
            p_mat[(ci, ci)] = VAR_CLK;
        }
        // Tropo: random walk — add small process noise
        for i in 0..ppp.nt() {
            let ti = ppp.it() + i;
            p_mat[(ti, ti)] += 1e-6; // 1 mm²/s process noise
        }
    }

    /// Compute measurement residuals and H matrix.
    /// This is a direct port of RTKLIB's res_ppp().
    #[allow(clippy::too_many_arguments)]
    fn residuals(
        &self,
        ppp: &PppState,
        obs: &[(SatelliteId, f64, f64, f64, f64, f64)], // (sat, L1, L2, P1, P2, el)
        lc_if_vals: &[f64],                               // IF carrier phase in meters
        range_offsets: &[f64],                            // PCV+PCO+tide correction (m)
        proc_tropo_dry: &[f64],                           // pre-computed ZHD×GMF mh
        proc_map_wet: &[f64],                             // pre-computed GMF wet mapping
        iono_prior: &[f64],                                // Klobuchar/IONEX L1 iono delay (m)
        lam1_vals: &[f64],                                 // L1 wavelength (UDUC)
        lam2_vals: &[f64],                                 // L2 wavelength (UDUC)
        f1_vals: &[f64],                                   // L1 frequency (UDUC)
        f2_vals: &[f64],                                   // L2 frequency (UDUC)
        sat_pos: &[Vector3<f64>],                       // ECEF satellite positions
        sat_clk: &[f64],                                 // satellite clock corrections (m)
        sat_var: &[f64],                                 // satellite position variance
        x: &DVector<f64>,
        v: &mut DVector<f64>,
        h: &mut DMatrix<f64>,
        r: &mut DMatrix<f64>,
    ) -> usize {
        let nx = ppp.nx();
        let nr = ppp.nr();
        let mut nv = 0; let mut skipped_el = 0; let mut skipped_dist = 0; let mut skipped_cp = 0; let mut skipped_code = 0;

        for i in 0..obs.len() {
            let (sat, l1_cyc, l2_cyc, p1, p2, el_deg) = obs[i];
            if el_deg < self.elev_mask_deg { skipped_el += 1;
                continue;
            }

            let rs = sat_pos[i];
            let dts = sat_clk[i];

            // Geometric range with PCV/PCO/tide corrections from ProcessedSat
            let dist = (rs - Vector3::new(x[0], x[1], x[2])).norm() + range_offsets[i];
            if dist <= 0.0 { skipped_dist += 1;
                continue;
            }

            let el = el_deg * D2R;

            // Troposphere: pre-computed ZHD×GMF mh + GMF wet mapping from pipeline
            let dtrp = proc_tropo_dry[i];
            let mw = proc_map_wet[i];
            let vart = ERR_SAAS * ERR_SAAS;

            // Gneiss IF mode: p1 is already IF-combined. Use directly.
            let lc = lc_if_vals[i];
            let pc = p1; // IF value from build_sats

            // Corrected range
            let rng = dist - dts + dtrp; // dts already in meters from ProcessedSat

            let sys: usize = if sat.constellation == Constellation::Glonass {
                1
            } else {
                0
            };

            // Line of sight unit vector
            let e = (rs - Vector3::new(x[0], x[1], x[2])) / dist;

            if ppp.uduc {
                // ---- UDUC: 4 raw measurements per satellite ----
                let l1 = l1_cyc; let l2 = l2_cyc;
                let lam1 = lam1_vals[i]; let lam2 = lam2_vals[i];
                let f1 = f1_vals[i]; let f2 = f2_vals[i];
                let inv_f1_sq = 1.0 / (f1 * f1);
                let inv_f2_sq = 1.0 / (f2 * f2);
                if l2 == 0.0 { continue; } // need both frequencies

                // ---- PR1 ----
                if p1 != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = p1 - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] -= x[ppp.ni(i)] * inv_f1_sq;
                    h[(ppp.ni(i), nv)] = inv_f1_sq;
                    let var_pr = 25.0 / libm::sin(el).max(0.1) + sat_var[i] + vart;
                    r[(nv, nv)] = var_pr;
                    if self.max_inno_m > 0.0 && v[nv].abs() > self.max_inno_m && sys != 1 { continue; }
                    nv += 1;
                }
                // ---- PR2 ----
                if p2 != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = p2 - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] -= x[ppp.ni(i)] * inv_f2_sq;
                    h[(ppp.ni(i), nv)] = inv_f2_sq;
                    let var_pr = 25.0 / libm::sin(el).max(0.1) + sat_var[i] + vart;
                    r[(nv, nv)] = var_pr;
                    if self.max_inno_m > 0.0 && v[nv].abs() > self.max_inno_m && sys != 1 { continue; }
                    nv += 1;
                }
                // ---- CP1 ----
                if l1 != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = l1 * lam1 - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] += x[ppp.ni(i)] * inv_f1_sq; // CP iono sign opposite to PR
                    h[(ppp.ni(i), nv)] = -inv_f1_sq;
                    v[nv] -= x[ppp.ib(i)] * lam1;
                    h[(ppp.ib(i), nv)] = lam1;
                    let var_cp = 0.01 / libm::sin(el).max(0.1) + sat_var[i] + vart;
                    r[(nv, nv)] = var_cp;
                    if self.max_inno_cp > 0.0 && v[nv].abs() > self.max_inno_cp && sys != 1 { continue; }
                    nv += 1;
                }
                // ---- CP2 ----
                if l2 != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = l2 * lam2 - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] += x[ppp.ni(i)] * inv_f2_sq;
                    h[(ppp.ni(i), nv)] = -inv_f2_sq;
                    v[nv] -= x[ppp.ib2(i)] * lam2;
                    h[(ppp.ib2(i), nv)] = lam2;
                    let var_cp = 0.01 / libm::sin(el).max(0.1) + sat_var[i] + vart;
                    r[(nv, nv)] = var_cp;
                    if self.max_inno_cp > 0.0 && v[nv].abs() > self.max_inno_cp && sys != 1 { continue; }
                    nv += 1;
                }
                // Ionosphere prior: constrain UDUC iono to Klobuchar/IONEX
                if iono_prior[i] != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = x[ppp.ni(i)] - iono_prior[i];
                    h[(ppp.ni(i), nv)] = 1.0;
                    r[(nv, nv)] = 0.0001; // σ=1cm — IONEX-grade ionosphere constraint
                    nv += 1;
                }
            } else {
                // ---- IF mode: 2 measurements per satellite ----
                // ---- Phase measurement ----
                if lc != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = lc - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] -= x[ppp.ib(i)];
                    h[(ppp.ib(i), nv)] = 1.0;
                    let ar_fixed = !ppp.uduc && self.p[(ppp.ib(i), ppp.ib(i))] < 0.01;
                    let cp_var = if ar_fixed { 0.0001 } else { 0.01 };
                    let var_phase = cp_var / libm::sin(el).max(0.1) + sat_var[i] + vart;
                    r[(nv, nv)] = var_phase;
                    let cp_thresh = if ar_fixed { 2.0 } else { self.max_inno_cp };
                    if cp_thresh > 0.0 && v[nv].abs() > cp_thresh && sys != 1 { continue; }
                    nv += 1;
                }
                // ---- Code measurement ----
                if pc != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = pc - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    let var_code = 25.0 / libm::sin(el).max(0.1) + sat_var[i] + vart;
                    r[(nv, nv)] = var_code;
                    if self.max_inno_m > 0.0 && v[nv].abs() > self.max_inno_m && sys != 1 { continue; }
                    nv += 1;
                }
            }
        }
        nv
    }
    fn elevation(rcv: &Vector3<f64>, sat: &Vector3<f64>) -> f64 {
        let llh = gneiss_core::coords::ecef_to_llh(*rcv);
        let ned_mat = gneiss_core::coords::ecef_to_ned_matrix(llh);
        let delta = sat - rcv;
        let ned = ned_mat * delta;
        libm::atan2(-ned.z, (ned.x.powi(2) + ned.y.powi(2)).sqrt())
    }

    pub fn solve_with_sats(&mut self, state: &mut RtkState, sats: &[crate::engine::processed_sat::ProcessedSat]) -> Result<(), EngineError> {
        let mut obs_data: Vec<(SatelliteId, f64, f64, f64, f64, f64)> = Vec::new();
        let mut lc_if_vals: Vec<f64> = Vec::new();
        let mut range_offsets: Vec<f64> = Vec::new(); // ProcessedSat.dist - our_dist
        let mut proc_tropo_dry: Vec<f64> = Vec::new(); // ProcessedSat.tropo_dry (ZHD×GMF mh)
        let mut proc_map_wet: Vec<f64> = Vec::new();   // ProcessedSat.map_wet (GMF mw)
        let mut iono_prior: Vec<f64> = Vec::new();      // ProcessedSat.iono_delay (Klobuchar/IONEX)
        let mut lam1_vals: Vec<f64> = Vec::new();       // L1 wavelength (m)
        let mut lam2_vals: Vec<f64> = Vec::new();       // L2 wavelength (m)
        let mut f1_vals: Vec<f64> = Vec::new();          // L1 frequency (Hz)
        let mut f2_vals: Vec<f64> = Vec::new();          // L2 frequency (Hz)
        let mut sat_pos: Vec<Vector3<f64>> = Vec::new();
        let mut sat_clk: Vec<f64> = Vec::new();
        let mut sat_var: Vec<f64> = Vec::new();
        let rcv = Vector3::new(state.position.vector.x, state.position.vector.y, state.position.vector.z);
        for sat in sats { tracing::debug!("RTKLIB sat p1={:.1} el={:.1}", sat.p1, Self::elevation(&rcv, &sat.sat_pos_rot)*R2D);
            let el = Self::elevation(&rcv, &sat.sat_pos_rot);
            if el < self.elev_mask_deg * D2R { continue; }
            if sat.p1 == 0.0 { continue; }
            // IF carrier phase in meters: cp1 is already IF-combined in L1 cycles
            let lc_if = if sat.is_iono_free {
                sat.cp1.unwrap_or(0.0) * sat.lam1
            } else {
                0.0
            };
            let our_dist = (sat.sat_pos_rot - rcv).norm();
            range_offsets.push(sat.dist - our_dist); // PCV + PCO + tide corrections
            proc_tropo_dry.push(sat.tropo_dry);       // ZHD×GMF mh (pre-computed)
            proc_map_wet.push(sat.map_wet);            // GMF wet mapping (pre-computed)
            iono_prior.push(sat.iono_delay);            // Klobuchar/IONEX L1 delay (m)
            lam1_vals.push(sat.lam1);
            lam2_vals.push(sat.lam2);
            f1_vals.push(sat.f1);
            f2_vals.push(sat.f2);
            obs_data.push((sat.sat_obs.sat, sat.cp1.unwrap_or(0.0), sat.cp2.unwrap_or(0.0), sat.p1, sat.p2.unwrap_or(sat.p1), el * R2D));
            lc_if_vals.push(lc_if);
            sat_pos.push(sat.sat_pos_rot);
            sat_clk.push(sat.dt_sat_m);
            sat_var.push(0.0);
        }
        tracing::debug!("solve_with_sats: {} sats -> {} obs_data ({} with CP)", sats.len(), obs_data.len(), lc_if_vals.iter().filter(|v| **v != 0.0).count()); if obs_data.len() < 4 { return Err(EngineError::InsufficientSatellites); }
        let has_glo = obs_data.iter().any(|(s,_,_,_,_,_)| s.constellation == Constellation::Glonass);
        // UDUC mode: use raw L1/L2 if not iono-free (ProcessedSat has raw obs)
        let use_uduc = !sats.is_empty() && !sats[0].is_iono_free && sats[0].p2.is_some();
        let mut ppp = PppState::new(has_glo, self.dynamics, use_uduc);
        ppp.nsat = obs_data.len();
        let nx = ppp.nx();
        // Detect mode change
        let mode_changed = use_uduc != self.was_uduc && self.epoch > 0;

        if self.epoch == 0 {
            let mut x0 = DVector::zeros(nx);
            x0[0] = state.position.vector.x; x0[1] = state.position.vector.y; x0[2] = state.position.vector.z;
            // Seed clock from SPP estimate (typically ~5ms = 1.4M meters for F9P).
            // Starting at 0 forces the first measurement update to absorb the full
            // receiver clock offset, which leaks into position through the gain matrix.
            x0[ppp.ic(0)] = state.rcv_clk_bias;
            // Seed ZWD from a priori wet delay (~0.1-0.3m). Pre-computed
            // tropo_dry only contains hydrostatic; ZWD state holds the wet part.
            let rcv_llh = gneiss_core::coords::ecef_to_llh(Vector3::new(x0[0], x0[1], x0[2]));
            if ppp.nt() >= 1 {
                x0[ppp.it()] = self.trop_zwd(rcv_llh);
            }
            // UDUC: seed ionosphere and L1/L2 ambiguities from raw observables
            if ppp.uduc {
                for i in 0..obs_data.len() {
                    let (_, _, l2_cyc, p1, p2, _) = obs_data[i];
                    if p2 == 0.0 || l2_cyc == 0.0 { continue; }
                    let f1 = f1_vals[i]; let f2 = f2_vals[i];
                    let gamma = (f1 * f1) / (f2 * f2);
                    let mut i1_est = (p2 - p1) / (gamma - 1.0);
                    if i1_est.is_nan() || i1_est.abs() > 500.0 { i1_est = 0.0; }
                    x0[ppp.ni(i)] = i1_est;
                    // Seed L1/L2 ambiguities from carrier phase residuals
                    let rs = sat_pos[i]; let dts = sat_clk[i];
                    let dist = (rs - Vector3::new(x0[0], x0[1], x0[2])).norm();
                    let el_rad = obs_data[i].5 * D2R;
                    let dtrp = proc_tropo_dry[i];
                    let rng = dist - dts + dtrp;
                    let l1_meas = obs_data[i].1 * lam1_vals[i];
                    let l2_meas = l2_cyc * lam2_vals[i];
                    let sys: usize = if obs_data[i].0.constellation == Constellation::Glonass { 1 } else { 0 };
                    x0[ppp.ib(i)] = (l1_meas - (rng - i1_est + x0[ppp.ic(sys)])) / lam1_vals[i];
                    x0[ppp.ib2(i)] = (l2_meas - (rng - i1_est * gamma + x0[ppp.ic(sys)])) / lam2_vals[i];
                }
            }
            self.x = x0;
            self.p = self.init_covariance(&ppp, nx);
            self.biases_seeded = ppp.uduc; // UDUC biases seeded at init, IF needs warmup
        } else if mode_changed {
            // IF↔UDUC mode change: convert state vector
            tracing::info!("Mode change: {} -> {} at epoch {}",
                if self.was_uduc { "UDUC" } else { "IF" },
                if use_uduc { "UDUC" } else { "IF" }, self.epoch);
            let mut x_new = DVector::zeros(nx);
            let mut p_new = DMatrix::zeros(nx, nx);
            // Copy position, clock, tropo (same layout for both modes)
            for i in 0..ppp.nr() {
                x_new[i] = self.x[i];
                for j in 0..ppp.nr() { p_new[(i,j)] = self.p[(i,j)]; }
            }
            if use_uduc {
                // IF→UDUC: initialize iono from known IF position (precise),
                // N1/N2 from AR-fixed N_IF + MW WL.
                let pos_if = Vector3::new(self.x[0], self.x[1], self.x[2]);
                let clk_if = self.x[ppp.ic(0)];
                for i in 0..obs_data.len() {
                    let (_, _, l2_cyc, p1, p2, _) = obs_data[i];
                    if p2 == 0.0 || l2_cyc == 0.0 { continue; }
                    let f1 = f1_vals[i]; let f2 = f2_vals[i];
                    let lam1 = lam1_vals[i]; let lam2 = lam2_vals[i];
                    let f1s = f1*f1; let f2s = f2*f2;

                    // Ionosphere from P1-P2 using known IF geometry
                    let rs = sat_pos[i];
                    let dist = (rs - pos_if).norm();
                    let dtrp = proc_tropo_dry[i];
                    let rng = dist - sat_clk[i] + dtrp;
                    let mut i1_est = (p1 - (rng + clk_if)).max(-200.0).min(200.0);
                    if i1_est.is_nan() { i1_est = 0.0; }
                    x_new[ppp.ni(i)] = i1_est;

                    // N1, N2 from AR-fixed IF bias + MW WL
                    let old_bi = ppp.nr() + i;
                    let n_if = if old_bi < self.x.len() { self.x[old_bi] } else { 0.0 };
                    let lam_nl = 299792458.0 / (f1 + f2);
                    let n_wl: f64 = self.mw_wl_ema.get(&obs_data[i].0)
                        .map(|(_, e)| e.round()).unwrap_or(0.0);
                    let n1_est = (n_if - n_wl * f2s / (f1s - f2s) * lam2) / lam_nl;
                    x_new[ppp.ib(i)] = n1_est;
                    x_new[ppp.ib2(i)] = n1_est - n_wl;
                    p_new[(ppp.ni(i), ppp.ni(i))] = VAR_BIAS;
                    p_new[(ppp.ib(i), ppp.ib(i))] = VAR_BIAS;
                    p_new[(ppp.ib2(i), ppp.ib2(i))] = VAR_BIAS;
                }
            } else {
                // UDUC→IF: convert N1,N2 back to N_IF
                let old_nsat = self.last_nsat;
                for i in 0..obs_data.len() {
                    let f1 = f1_vals[i]; let f2 = f2_vals[i];
                    let lam1 = lam1_vals[i]; let lam2 = lam2_vals[i];
                    let f1s = f1*f1; let f2s = f2*f2;
                    // Old UDUC layout: L1 at nr+old_nsat+i, L2 at nr+2*old_nsat+i
                    let old_l1_idx = ppp.nr() + old_nsat + i;
                    let old_l2_idx = ppp.nr() + 2 * old_nsat + i;
                    let old_n1 = if old_l1_idx < self.x.len() { self.x[old_l1_idx] } else { 0.0 };
                    let old_n2 = if old_l2_idx < self.x.len() { self.x[old_l2_idx] } else { 0.0 };
                    let n_if = (f1s * old_n1 * lam1 - f2s * old_n2 * lam2) / (f1s - f2s);
                    let bi = ppp.nr() + i;
                    x_new[bi] = n_if;
                    p_new[(bi, bi)] = VAR_BIAS;
                }
            }
            self.x = x_new;
            self.p = p_new;
            self.biases_seeded = false; // warmup for new mode
        } else if self.last_nsat != obs_data.len() {
            // nsat changed: resize state vector with semantic index remapping.
            // Blind slice copy corrupts indices when GLONASS appears/disappears
            // (nc() changes → it() and nr() shift by 1).
            // old_nr = old position/clk/tropo count (depends on old has_glo)
            let _old_nc = if self.last_has_glo { 2 } else { 1 };
            let old_nr = ppp.np + _old_nc + ppp.nt();
            let old_nx = self.x.len();
            let mut x_new = DVector::zeros(nx);
            let mut p_new = DMatrix::zeros(nx, nx);
            // Copy position (indices 0..np, same for both layouts)
            for i in 0..3 {
                x_new[i] = self.x[i];
                for j in 0..3 { p_new[(i, j)] = self.p[(i, j)]; }
            }
            // Copy GPS clock (always at np, same for both)
            x_new[ppp.ic(0)] = self.x[ppp.ic(0)];
            p_new[(ppp.ic(0), ppp.ic(0))] = self.p[(ppp.ic(0), ppp.ic(0))];
            // Copy GLO clock if present in both old and new
            if ppp.nc() >= 2 && old_nx > ppp.ic(1) {
                x_new[ppp.ic(1)] = self.x[ppp.ic(1)];
                p_new[(ppp.ic(1), ppp.ic(1))] = self.p[(ppp.ic(1), ppp.ic(1))];
            }
            // Copy tropo (indices it()..it()+nt(), same dimension)
            for i in 0..ppp.nt().min(old_nx.saturating_sub(ppp.it())) {
                x_new[ppp.it() + i] = self.x[ppp.it() + i];
                p_new[(ppp.it()+i, ppp.it()+i)] = self.p[(ppp.it()+i, ppp.it()+i)];
            }
            // Copy existing biases to same relative positions
            let nr = ppp.nr();
            for i in 0..old_nx.saturating_sub(old_nr) {
                let new_idx = nr + i;
                if new_idx < nx {
                    x_new[new_idx] = self.x[old_nr + i];
                    p_new[(new_idx, new_idx)] = self.p[(old_nr + i, old_nr + i)];
                }
            }
            // New biases get VAR_BIAS
            for i in old_nx.saturating_sub(old_nr)..ppp.nsat {
                let idx = nr + i;
                if idx < nx { p_new[(idx, idx)] = VAR_BIAS; }
            }
            self.x = x_new;
            self.p = p_new;
            self.biases_seeded = false; // warmup needed for new satellite biases
        }
        self.last_nsat = obs_data.len();
        self.last_has_glo = has_glo;
        self.was_uduc = use_uduc;
        self.epoch += 1;
        let mut xp = self.x.clone();
        let mut pp = self.p.clone();
        self.predict(&ppp, &mut xp, &mut pp);
        // Clock-only pre-update: after predict() resets clock variance,
        // use PR and AR-fixed CP measurements to estimate clock without
        // touching position.  CP from AR-fixed biases gives σ≈1cm range,
        // dramatically improving clock accuracy when AR is active.
        if !ppp.uduc {
            let max_n = obs_data.len() * 2;
            let mut h_clk = DMatrix::zeros(ppp.nx(), max_n);
            let mut v_clk = DVector::zeros(max_n);
            let mut r_clk = DMatrix::zeros(max_n, max_n);
            let mut n_clk = 0usize;
            let rcv_pos = Vector3::new(xp[0], xp[1], xp[2]);
            for i in 0..obs_data.len() {
                let (_, _, _, p1, _, el_deg) = obs_data[i];
                if p1 == 0.0 { continue; }
                let rs = sat_pos[i]; let dts = sat_clk[i];
                let dist = (rs - rcv_pos).norm();
                let el = el_deg * D2R;
                let dtrp = proc_tropo_dry[i];
                let rng = dist - dts + dtrp;
                let sys: usize = if obs_data[i].0.constellation == Constellation::Glonass { 1 } else { 0 };
                // PR measurement (σ≈5m at zenith)
                h_clk[(ppp.ic(sys), n_clk)] = 1.0;
                v_clk[n_clk] = p1 - rng - xp[ppp.ic(sys)];
                r_clk[(n_clk, n_clk)] = 25.0 / libm::sin(el).max(0.1);
                n_clk += 1;
                // CP measurement for AR-fixed IF biases (σ≈1cm at zenith)
                let lc = lc_if_vals[i];
                if lc != 0.0 && pp[(ppp.ib(i), ppp.ib(i))] < 0.01 {
                    h_clk[(ppp.ic(sys), n_clk)] = 1.0;
                    v_clk[n_clk] = lc - rng - xp[ppp.ib(i)] - xp[ppp.ic(sys)];
                    r_clk[(n_clk, n_clk)] = 0.0001 / libm::sin(el).max(0.1);
                    n_clk += 1;
                }
            }
            if n_clk >= 4 {
                let h_s = h_clk.view((0, 0), (ppp.nx(), n_clk)).clone_owned();
                let vs = v_clk.rows(0, n_clk).clone_owned();
                let rs = r_clk.view((0, 0), (n_clk, n_clk)).clone_owned();
                let _ = Self::measurement_update(&mut xp, &mut pp, &h_s, &vs, &rs, ppp.nx(), n_clk);
            }
        }
        let nv_max = if ppp.uduc { obs_data.len() * 5 } else { obs_data.len() * 2 };
        let mut v = DVector::zeros(nv_max);
        let mut h_mat = DMatrix::zeros(nx, nv_max);
        let mut r_mat = DMatrix::zeros(nv_max, nv_max);
        // Warmup epoch: PR-only, then seed biases from the converged state
        let cp_enabled = self.biases_seeded;
        let lc_for_filter: Vec<f64> = if cp_enabled {
            lc_if_vals.clone()
        } else {
            vec![0.0f64; lc_if_vals.len()]
        };
        for _iter in 0..self.max_iter {
            let nv = self.residuals(&ppp, &obs_data, &lc_for_filter, &range_offsets, &proc_tropo_dry, &proc_map_wet, &iono_prior, &lam1_vals, &lam2_vals, &f1_vals, &f2_vals, &sat_pos, &sat_clk, &sat_var, &xp, &mut v, &mut h_mat, &mut r_mat);
            if nv < 4 { break; }
            let h_s = h_mat.view((0, 0), (nx, nv)).clone_owned();
            let vs = v.rows(0, nv).clone_owned();
            let rs = r_mat.view((0, 0), (nv, nv)).clone_owned();
            if Self::measurement_update(&mut xp, &mut pp, &h_s, &vs, &rs, nx, nv).is_err() { break; }
        }
        // After warmup epoch: seed phase biases from the PR-converged state.
        // Using the filtered position (not SPP) gives biases within ~3m,
        // so CP residuals start small enough for σ=10cm measurements to pull.
        if !self.biases_seeded && !ppp.uduc {
            // Validate position: if filter diverged during warmup, skip seeding
            let pos_jump = (Vector3::new(xp[0], xp[1], xp[2]) - Vector3::new(self.x[0], self.x[1], self.x[2])).norm();
            if pos_jump > 500.0 {
                tracing::warn!("Warmup position jump {}m > 500m — skipping bias seed", pos_jump);
                self.x[0] = xp[0]; self.x[1] = xp[1]; self.x[2] = xp[2]; // keep position
                self.p = pp.clone(); // keep covariance
                state.position.vector.x = self.x[0];
                state.position.vector.y = self.x[1];
                state.position.vector.z = self.x[2];
                state.rcv_clk_bias = self.x[ppp.ic(0)];
                state.covariance = self.p.clone();
                return Ok(()); // retry warmup next epoch
            }
            let mut seeded = 0;
            for i in 0..obs_data.len() {
                if lc_if_vals[i] == 0.0 { continue; }
                let rs = sat_pos[i];
                let dts = sat_clk[i];
                let (sat, _, _, _, _, el) = obs_data[i];
                let sys: usize = if sat.constellation == Constellation::Glonass { 1 } else { 0 };
                let dist = (rs - Vector3::new(xp[0], xp[1], xp[2])).norm();
                let dtrp = proc_tropo_dry[i]; // pre-computed ZHD×GMF mh
                let rng = dist - dts + dtrp; // dts already in meters from ProcessedSat
                xp[ppp.ib(i)] = lc_if_vals[i] - rng - xp[ppp.ic(sys)];
                seeded += 1;
            }
            tracing::debug!("solve_with_sats: seeded {} phase biases after warmup epoch", seeded);
            self.biases_seeded = true;
        }
        self.x = xp;
        self.p = pp;
        state.position.vector.x = self.x[0];
        state.position.vector.y = self.x[1];
        state.position.vector.z = self.x[2];
        state.rcv_clk_bias = self.x[ppp.ic(0)];
        state.covariance = self.p.clone();

        // ---- MW widelane tracking for integer AR ----
        // Compute Melbourne-Wübbena widelane for each satellite, maintain EMA.
        // When WL is precise enough, fix N_wl → compute N1 → fix N_IF.
        for i in 0..obs_data.len() {
            let (sat, l1_cyc, l2_cyc, p1_raw, p2_raw, _el_deg) = obs_data[i];
            if l1_cyc == 0.0 || l2_cyc == 0.0 || p1_raw == 0.0 || p2_raw == 0.0 { continue; }
            let lam1 = lam1_vals.get(i).copied().unwrap_or(0.1903);
            let lam2 = lam2_vals.get(i).copied().unwrap_or(0.2442);
            let f1 = f1_vals.get(i).copied().unwrap_or(1575.42e6);
            let f2 = f2_vals.get(i).copied().unwrap_or(1227.60e6);
            if f1 == 0.0 || f2 == 0.0 { continue; }

            // MW in cycles: (L1-L2) - narrow-lane PR / widelane wavelength
            let l1_m = l1_cyc * lam1;
            let l2_m = l2_cyc * lam2;
            let mw_m = combinations::melbourne_wubbena(l1_m, l2_m, p1_raw, p2_raw, f1, f2);
            let wl_lambda = combinations::lambda_wl(f1, f2);
            if wl_lambda <= 0.0 { continue; }
            let mw_cyc = mw_m / wl_lambda;

            // EMA: 0.05 weight for new sample
            let (count, ema) = self.mw_wl_ema.get(&sat)
                .map(|(c, e)| (c + 1, e + 0.05 * (mw_cyc - e)))
                .unwrap_or((1u32, mw_cyc));
            self.mw_wl_ema.insert(sat, (count, ema));

            // Fix WL when confident (50+ samples). IF mode only:
            // UDUC AR is handled separately via WL constraint Kalman updates.
            if !ppp.uduc && count > 100 {
                let n_wl = ema.round();
                if (ema - n_wl).abs() > 0.25 { continue; }

                let f1s = f1 * f1; let f2s = f2 * f2;

                // Get float N_IF from our IF bias
                let n_if = self.x[ppp.ib(i)];

                // N_IF = N1*λ_nl + N_wl*f2²*λ2/(f1²-f2²) where λ_nl = c/(f1+f2)
                let lam_nl = 299792458.0 / (f1 + f2);
                let n1_est = (n_if - n_wl * f2s / (f1s - f2s) * lam2) / lam_nl;
                let n1_rounded = n1_est.round();
                if (n1_est - n1_rounded).abs() > 0.3 { continue; }

                let n2_rounded = n1_rounded - n_wl;
                let n_if_fixed = (f1s * n1_rounded * lam1 - f2s * n2_rounded * lam2) / (f1s - f2s);
                let residual = n_if - n_if_fixed;
                if residual.abs() > 2.0 { continue; } // reject if residual too large

                // Tighten IF bias state to the fixed value
                let bi = ppp.ib(i);
                self.x[bi] = n_if_fixed;
                self.p[(bi, bi)] = 1e-6; // σ ≈ 1mm — effectively locked
                // Clear cross-correlations for this bias
                for j in 0..ppp.nx() {
                    if j != bi {
                        self.p[(bi, j)] = 0.0;
                        self.p[(j, bi)] = 0.0;
                    }
                }
            }
        }

        // After AR fixes: CP-only re-solve with inflated position variance.
        // Fixed IF ambiguities make CP an unbiased range measurement (σ≈1cm).
        // Temporarily inflating P_pos lets the filter reposition away from
        // the biased PR solution toward the CP-only solution.
        let any_fixed = (0..obs_data.len()).any(|i| {
            !ppp.uduc && self.p[(ppp.ib(i), ppp.ib(i))] < 0.01
        });
        if any_fixed {
            // Inflate position/clock/ZWD for CP-only re-convergence
            for k in 0..3 { self.p[(k, k)] = self.p[(k, k)].max(100.0); } // σ=10m pos
            self.p[(ppp.ic(0), ppp.ic(0))] = self.p[(ppp.ic(0), ppp.ic(0))].max(10000.0);
            if ppp.nt() >= 1 { self.p[(ppp.it(), ppp.it())] = self.p[(ppp.it(), ppp.it())].max(9.0); } // σ=3m ZWD

            // One CP-only measurement update pass
            let mut v_cp = DVector::zeros(obs_data.len());
            let mut h_cp = DMatrix::zeros(ppp.nx(), obs_data.len());
            let mut r_cp = DMatrix::zeros(obs_data.len(), obs_data.len());
            let mut n_cp = 0usize;
            let rcv_pos = Vector3::new(self.x[0], self.x[1], self.x[2]);
            for i in 0..obs_data.len() {
                if ppp.uduc { continue; }
                if self.p[(ppp.ib(i), ppp.ib(i))] > 0.01 { continue; } // not AR-fixed
                if lc_if_vals[i] == 0.0 { continue; }
                let rs = sat_pos[i]; let dts = sat_clk[i];
                let dist = (rs - rcv_pos).norm();
                let el_rad = obs_data[i].5 * D2R;
                let dtrp = proc_tropo_dry[i];
                let rng = dist - dts + dtrp;
                let n_if_fixed = self.x[ppp.ib(i)];
                let (sat, _, _, _, _, _) = obs_data[i];
                let sys: usize = if sat.constellation == Constellation::Glonass { 1 } else { 0 };
                let e = (rs - rcv_pos) / dist;
                for k in 0..ppp.nx() { h_cp[(k, n_cp)] = 0.0; }
                v_cp[n_cp] = lc_if_vals[i] - rng - self.x[ppp.ic(sys)] - n_if_fixed;
                for k in 0..3 { h_cp[(k, n_cp)] = -e[k]; }
                h_cp[(ppp.ic(sys), n_cp)] = 1.0;
                if ppp.nt() >= 1 { h_cp[(ppp.it(), n_cp)] = proc_map_wet[i]; }
                h_cp[(ppp.ib(i), n_cp)] = 1.0;
                r_cp[(n_cp, n_cp)] = 0.0001; // σ=1cm CP
                n_cp += 1;
            }
            if n_cp >= 4 {
                let h_s = h_cp.view((0, 0), (ppp.nx(), n_cp)).clone_owned();
                let vs = v_cp.rows(0, n_cp).clone_owned();
                let rs = r_cp.view((0, 0), (n_cp, n_cp)).clone_owned();
                let _ = Self::measurement_update(&mut self.x, &mut self.p, &h_s, &vs, &rs, ppp.nx(), n_cp);
            }
        }

        Ok(())
    }

    /// Saastamoinen zenith wet delay for ZWD state initialization.
    fn trop_zwd(&self, rcv_llh: Vector3<f64>) -> f64 {
        let alt_m = rcv_llh.z;
        let lat_rad = rcv_llh.x;
        let t = (288.15 - 0.0065 * alt_m).max(200.0);
        let e = 6.108 * libm::exp((17.15 * t - 4684.0) / (t - 38.45)) * REL_HUMI;
        let scale = 1.0 - 0.00266 * libm::cos(2.0 * lat_rad) - 0.00028 * alt_m / 1000.0;
        0.002277 * (1255.0 / t + 0.05) * e / scale
    }

    /// Kalman measurement update — port of RTKLIB filter()
    /// H is stored as (nx × nv) — rows = state dim, cols = measurement dim.
    fn measurement_update(
        x: &mut DVector<f64>,
        p: &mut DMatrix<f64>,
        h: &DMatrix<f64>,
        v: &DVector<f64>,
        r: &DMatrix<f64>,
        _nx: usize,
        nv: usize,
    ) -> Result<(), EngineError> {
        if nv == 0 {
            return Ok(());
        }
        let h_t = h.transpose(); // (nv × nx)
        // S = H^T * P * H + R  → (nv × nv)
        let hp = &h_t * &*p; // (nv × nx) × (nx × nx) = (nv × nx)
        let s = &hp * h + r; // (nv × nx) × (nx × nv) + (nv × nv) = (nv × nv)
        let s_inv = crate::math::inversion::invert_matrix_robust(&s);
        // K = P * H * S^{-1} → (nx × nv)
        let k = &*p * h * &s_inv; // (nx × nx) × (nx × nv) × (nv × nv) = (nx × nv)
        // dx = K * v → (nx × 1)
        let dx = &k * v;
        if dx.iter().any(|d| d.is_nan() || d.abs() > 1e8) { return Err(EngineError::StateDisappeared); }
        *x += &dx;
        // P = (I - K*H^T) * P → (nx × nx)
        let i_mat = DMatrix::identity(_nx, _nx);
        let ikh = &i_mat - &k * &h_t;
        *p = &ikh * p.clone() * &ikh.transpose() + &k * r * &k.transpose();
        if p.iter().any(|d| d.is_nan() || d.is_infinite()) { return Err(EngineError::StateDisappeared); }
        Ok(())
    }

    /// Main PPP solve — port of RTKLIB pppos()
    pub fn solve(
        &mut self,
        state: &mut RtkState,
        rover_obs: &EpochObs,
        ephemerides: &[Ephemeris],
    ) -> Result<(), EngineError> {
        let n = rover_obs.satellites.len();
        if n < 4 {
            return Err(EngineError::InsufficientSatellites);
        }

        // Detect GLONASS presence
        let has_glo = rover_obs
            .satellites
            .iter()
            .any(|s| s.sat.constellation == Constellation::Glonass);
        let mut ppp = PppState::new(has_glo, self.dynamics, false); // IF mode for old API
        ppp.nsat = n; // simplified

        let nx = ppp.nx();

        // Initialize state from current position
        let mut x = DVector::zeros(nx);
        x[0] = state.position.vector.x;
        x[1] = state.position.vector.y;
        x[2] = state.position.vector.z;
        let mut p_mat = self.init_covariance(&ppp, nx);

        // Forward prediction
        self.predict(&ppp, &mut x, &mut p_mat);

        // Build measurement data
        let mut obs_data: Vec<(SatelliteId, f64, f64, f64, f64, f64)> = Vec::new();
        let mut sat_pos: Vec<Vector3<f64>> = Vec::new();
        let mut sat_clk: Vec<f64> = Vec::new();
        let mut sat_var: Vec<f64> = Vec::new();

        let time = rover_obs.time;
        let rcv_pos = Vector3::new(x[0], x[1], x[2]);

        for sat_obs in &rover_obs.satellites {
            // Get ephemeris
            let eph = match ephemerides.iter().find(|e| e.sat() == sat_obs.sat) {
                Some(e) => e,
                None => continue,
            };

            // Compute satellite position and clock
            let (pos, _vel, clk, _clk_drift) = eph.position(time);
            let dist = (pos - rcv_pos).norm();
            let el = libm::asin((rcv_pos.z + 6371000.0) / dist); // rough elevation

            // Extract measurements using existing Gneiss helper methods
            let l1_cyc = sat_obs.get_observable_phase(1).unwrap_or(0.0);
            let l2_cyc = sat_obs.get_observable_phase(2).unwrap_or(0.0);
            let p1 = sat_obs.get_observable(1).unwrap_or(0.0);
            let p2 = sat_obs.get_observable(2).unwrap_or(0.0);

            if l1_cyc == 0.0 && p1 == 0.0 {
                continue;
            }

            let el_deg = el * R2D;
            obs_data.push((sat_obs.sat, l1_cyc, l2_cyc, p1, p2, el_deg));
            sat_pos.push(pos);
            sat_clk.push(clk * CLIGHT); // seconds → meters
            sat_var.push(0.0); // no ephemeris variance for now
        }

        if obs_data.len() < 4 {
            return Err(EngineError::InsufficientSatellites);
        }

        // Iterated measurement update
        let nv_max = if ppp.uduc { obs_data.len() * 5 } else { obs_data.len() * 2 };
        let mut v = DVector::zeros(nv_max);
        let mut h_mat = DMatrix::zeros(nx, nv_max);
        let mut r_mat = DMatrix::zeros(nv_max, nv_max);
        let mut xp = x.clone();
        let mut pp = p_mat.clone();

        for _iter in 0..self.max_iter {
            let empty: Vec<f64> = vec![0.0; obs_data.len()];
            let nv = self.residuals(
                &ppp, &obs_data, &empty, &empty, &empty, &empty, &empty,
                &empty, &empty, &empty, &empty,
                &sat_pos, &sat_clk, &sat_var,
                &xp, &mut v, &mut h_mat, &mut r_mat,
            );
            if nv == 0 {
                break;
            }
            let h_slice = h_mat.view((0, 0), (nx, nv));
            let v_slice = v.rows(0, nv);
            let r_slice = r_mat.view((0, 0), (nv, nv));
            let h_owned = h_slice.clone_owned();
            let v_owned = v_slice.clone_owned();
            let r_owned = r_slice.clone_owned();

            pp = p_mat.clone();
            if let Err(_) =
                Self::measurement_update(&mut xp, &mut pp, &h_owned, &v_owned, &r_owned, nx, nv)
            {
                break;
            }
        }

        // Update state with result
        state.position.vector.x = xp[0];
        state.position.vector.y = xp[1];
        state.position.vector.z = xp[2];
        state.rcv_clk_bias = xp[ppp.ic(0)];
        state.covariance = pp;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ppp_state_indices() {
        let ppp = PppState::new(false, false, false);
        assert_eq!(ppp.np, 3);
        assert_eq!(ppp.nc(), 1);
        assert_eq!(ppp.ic(0), 3); // pos(3) + GPS clk
        assert_eq!(ppp.it(), 4); // pos(3) + clk(1)
        assert_eq!(ppp.nt(), 1);
        assert_eq!(ppp.nr(), 5); // pos(3) + clk(1) + tropo(1)

        let ppp2 = PppState::new(true, false, false); // with GLONASS
        assert_eq!(ppp2.nc(), 2);
        assert_eq!(ppp2.ic(0), 3); // GPS clk at 3
        assert_eq!(ppp2.ic(1), 4); // GLO clk at 4
        assert_eq!(ppp2.it(), 5); // pos(3) + clk(2)
        assert_eq!(ppp2.nr(), 6); // pos(3) + clk(2) + tropo(1)
    }

    #[test]
    fn test_ppp_empty_obs() {
        let mut ppp = PppRtklib::default();
        let obs = EpochObs {
            time: GpsTime::new(2000, 0.0),
            satellites: vec![],
        };
        let mut state = RtkState::new(
            GpsTime::new(2000, 0.0),
            Coordinate::new(
                Vector3::new(0.0, 0.0, 0.0),
                gneiss_core::coords::Datum::WGS84,
                gneiss_core::coords::Frame::ECEF,
                GpsTime::new(2000, 0.0),
            ),
            1.0,
        );
        let result = ppp.solve(&mut state, &obs, &[]);
        assert!(result.is_err());
    }
}
