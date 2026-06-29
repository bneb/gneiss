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

use crate::engine::ppp_multi_epoch_batch::StaticPositionBatchSolver;
use crate::engine::EngineError;
use crate::filter::RtkState;
use crate::measurements::combinations;
use gneiss_core::coords::Coordinate;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

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

    /// Index of L2 ambiguity state (full UDUC — deprecated, use 2N hybrid)
    fn ib2(&self, _sat_idx: usize) -> usize {
        0 // not used in 2N hybrid mode
    }

    /// Total number of estimated states.
    /// UDUC hybrid: [pos, clk, tropo, iono(N), bias(N)] = nr + 2*nsat (well-determined)
    /// IF mode: [pos, clk, tropo, bias(N)] = nr + nsat
    fn nx(&self) -> usize {
        if self.uduc {
            self.nr() + 2 * self.nsat // iono + IF biases per sat
        } else {
            self.ib(self.nsat)
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
    pub enable_multi_epoch: bool, // feature flag
    pub initial_position_var: f64, // position variance for known initial pos (0 = use VAR_POS)
    pub x: DVector<f64>,
    pub p: DMatrix<f64>,
    pub epoch: u32,
    last_nsat: usize,
    last_has_glo: bool,    // for semantic state resize on GLO change
    biases_seeded: bool,
    was_uduc: bool,         // previous epoch's UDUC mode
    /// MW widelane EMA for AR: sat → (count, smoothed_WL_cycles)
    mw_wl_ema: HashMap<SatelliteId, (u32, f64)>,
    /// Multi-epoch batch position solver for static PPP.
    /// Uses CP-derived range measurements from AR-fixed satellites
    /// to jointly estimate position across epochs.
    pub batch_solver: StaticPositionBatchSolver,
    /// Position history for RTS backward smoothing (static receivers only).
    /// Stores forward-pass position and covariance at each epoch.
    position_history: Vec<PositionEntry>,
}

/// One epoch of forward-pass position data for backward smoothing.
struct PositionEntry {
    pos: Vector3<f64>,       // Updated position (after measurement update)
    cov: Matrix3<f64>,       // Updated position covariance (3×3)
    pos_pred: Vector3<f64>,  // Predicted position (before measurement update)
    cov_pred: Matrix3<f64>,  // Predicted position covariance (3×3)
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
            enable_multi_epoch: true,
            initial_position_var: 0.0, // 0 = use VAR_POS default
            mw_wl_ema: HashMap::new(),
            batch_solver: StaticPositionBatchSolver::new(),
            position_history: Vec::new(),
        }
    }
}

impl PppRtklib {
    /// Allocate and return the initial state covariance matrix (RTKLIB defaults).
    fn init_covariance(&self, ppp: &PppState, nx: usize) -> DMatrix<f64> {
        let mut p = DMatrix::zeros(nx, nx);
        let pos_var = if self.initial_position_var > 0.0 {
            self.initial_position_var
        } else {
            VAR_POS
        };
        for i in 0..3 {
            p[(i, i)] = pos_var;
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
                p[(ppp.ib(i), ppp.ib(i))] = VAR_BIAS;    // IF bias
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

            // ---- Phase measurement (IF, always) ----
            if lc != 0.0 {
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = lc - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] -= x[ppp.ib(i)];
                    h[(ppp.ib(i), nv)] = 1.0;
                    let ar_fixed = self.p[(ppp.ib(i), ppp.ib(i))] < 0.01;
                    // IF mode amplifies CP noise ~3× vs raw L1/L2 (σ=3cm vs 1cm)
                    let cp_var = if ar_fixed {
                        if ppp.uduc { 0.0001 } else { 0.001 }
                    } else {
                        0.01
                    };
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

                // Ionosphere estimation (UDUC hybrid): raw PR1/PR2 + IONEX prior
                if ppp.uduc && p2 != 0.0 && l2_cyc != 0.0 {
                    let f1 = f1_vals[i]; let f2 = f2_vals[i];
                    let inv_f1_sq = 1.0 / (f1 * f1);
                    let inv_f2_sq = 1.0 / (f2 * f2);
                    // PR1 with iono state
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = p1 - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] -= x[ppp.ni(i)] * inv_f1_sq;
                    h[(ppp.ni(i), nv)] = inv_f1_sq;
                    r[(nv, nv)] = 25.0 / libm::sin(el).max(0.1);
                    nv += 1;
                    // PR2 with iono state
                    for k in 0..nx { h[(k, nv)] = 0.0; }
                    v[nv] = p2 - rng;
                    for k in 0..3 { h[(k, nv)] = -e[k]; }
                    if sys != 1 { v[nv] -= x[ppp.ic(0)]; h[(ppp.ic(0), nv)] = 1.0; }
                    else { v[nv] -= x[ppp.ic(1)]; h[(ppp.ic(1), nv)] = 1.0; }
                    if ppp.nt() >= 1 { h[(ppp.it(), nv)] = mw; }
                    v[nv] -= x[ppp.ni(i)] * inv_f2_sq;
                    h[(ppp.ni(i), nv)] = inv_f2_sq;
                    r[(nv, nv)] = 25.0 / libm::sin(el).max(0.1);
                    nv += 1;
                    // IONEX iono prior
                    if iono_prior[i] != 0.0 {
                        for k in 0..nx { h[(k, nv)] = 0.0; }
                        v[nv] = x[ppp.ni(i)] - iono_prior[i];
                        h[(ppp.ni(i), nv)] = 1.0;
                        r[(nv, nv)] = 4.0; // σ=2m F9P-compatible
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
        let mut raw_p1_vals: Vec<f64> = Vec::new();     // Raw L1 PR (before IF override)
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
            raw_p1_vals.push(sat.p1); // Save raw L1 PR for MW WL AR (before IF override)
            lc_if_vals.push(lc_if);
            sat_pos.push(sat.sat_pos_rot);
            sat_clk.push(sat.dt_sat_m);
            sat_var.push(0.0);
        }
        tracing::debug!("solve_with_sats: {} sats -> {} obs_data ({} with CP)", sats.len(), obs_data.len(), lc_if_vals.iter().filter(|v| **v != 0.0).count()); if obs_data.len() < 4 { return Err(EngineError::InsufficientSatellites); }
        let has_glo = obs_data.iter().any(|(s,_,_,_,_,_)| s.constellation == Constellation::Glonass);

        // For static receivers with dual-freq data: compute IF combinations
        // from raw L1/L2 so that geometry-NL AR can fix IF ambiguities.
        // IF mode cancels ionosphere (1st order), reduces state dimension,
        // and enables CP re-solve with integer AR.
        // For static receivers with known position: use IF mode with NL AR.
        // The tight prior (σ=1cm) keeps position near truth while NL AR fixes
        // IF ambiguities for mm-level CP measurements.
        let force_if = self.initial_position_var > 0.0
            && !self.dynamics
            && !sats.is_empty() && !sats[0].is_iono_free && sats[0].p2.is_some();
        if force_if {
            let mut if_pr_count = 0u32;
            let mut if_cp_count = 0u32;
            for i in 0..obs_data.len() {
                let (sat, cp1_cyc, cp2_cyc, p1, p2, el) = obs_data[i];
                if p1 == 0.0 || p2 == 0.0 { continue; }
                let f1 = f1_vals[i]; let f2 = f2_vals[i];
                let f1_2 = f1 * f1; let f2_2 = f2 * f2;
                let denom = f1_2 - f2_2;
                if denom <= 0.0 { continue; }
                // IF pseudorange (meters)
                let if_pr = (f1_2 * p1 - f2_2 * p2) / denom;
                if_pr_count += 1;
                // IF carrier phase (meters)
                let if_cp_m = if cp1_cyc != 0.0 && cp2_cyc != 0.0 {
                    let cp1_m = cp1_cyc * lam1_vals[i];
                    let cp2_m = cp2_cyc * lam2_vals[i];
                    if_cp_count += 1;
                    (f1_2 * cp1_m - f2_2 * cp2_m) / denom
                } else {
                    0.0
                };
                obs_data[i] = (sat, cp1_cyc, cp2_cyc, if_pr, p2, el);
                lc_if_vals[i] = if_cp_m;
            }
            if self.epoch == 0 {
                tracing::info!(
                    "IF mode forced: {} IF PR, {} IF CP",
                    if_pr_count, if_cp_count,
                );
                // DEBUG: compare raw vs IF for first sat
                if obs_data.len() > 0 {
                    let (_, _, _, raw_p1, raw_p2, _) = obs_data[0];
                    let (_, _, _, if_p1, _, _) = obs_data[0];
                    let f1 = f1_vals[0]; let f2 = f2_vals[0];
                    let f1_2 = f1*f1; let f2_2 = f2*f2;
                    let denom = f1_2 - f2_2;
                    let if_check = (f1_2 * raw_p1 - f2_2 * raw_p2) / denom;
                    tracing::info!(
                        "IF debug: raw P1={:.1} P2={:.1} IF={:.1} stored={:.1} f1={:.3}MHz f2={:.3}MHz",
                        raw_p1, raw_p2, if_check, if_p1, f1/1e6, f2/1e6
                    );
                }
            }
        }
        // UDUC mode: use raw L1/L2 if not iono-free (ProcessedSat has raw obs).
        // Skip UDUC if we forced IF above (needed for static multi-epoch).
        let use_uduc = !force_if
            && !sats.is_empty() && !sats[0].is_iono_free && sats[0].p2.is_some();
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
            // UDUC hybrid: seed ionosphere from IONEX/Klobuchar prior (iono_prior[i])
            if ppp.uduc {
                for i in 0..obs_data.len() {
                    if iono_prior[i] != 0.0 {
                        x0[ppp.ni(i)] = iono_prior[i];
                    }
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
                // IF→UDUC hybrid: copy IF biases, seed ionosphere from IONEX
                for i in 0..obs_data.len() {
                    if iono_prior[i] != 0.0 {
                        x_new[ppp.ni(i)] = iono_prior[i];
                    }
                    p_new[(ppp.ni(i), ppp.ni(i))] = VAR_BIAS;
                    // Copy IF bias from old state
                    let old_bi = ppp.nr() + i;
                    let new_bi = ppp.ib(i);
                    if old_bi < self.x.len() {
                        x_new[new_bi] = self.x[old_bi];
                        p_new[(new_bi, new_bi)] = self.p[(old_bi, old_bi)];
                    } else {
                        p_new[(new_bi, new_bi)] = VAR_BIAS;
                    }
                }
            } else {
                // UDUC→IF: copy IF biases (skip iono states)
                for i in 0..obs_data.len() {
                    let bi = ppp.nr() + i;
                    let old_bi = ppp.ib(i);
                    if old_bi < self.x.len() {
                        x_new[bi] = self.x[old_bi];
                        p_new[(bi, bi)] = self.p[(old_bi, old_bi)].max(VAR_BIAS);
                    } else {
                        p_new[(bi, bi)] = VAR_BIAS;
                    }
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
        // Save predicted position for RTS backward smoothing
        let pos_pred = Vector3::new(xp[0], xp[1], xp[2]);
        let cov_pred = Matrix3::new(
            pp[(0,0)], pp[(0,1)], pp[(0,2)],
            pp[(1,0)], pp[(1,1)], pp[(1,2)],
            pp[(2,0)], pp[(2,1)], pp[(2,2)],
        );
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
                    let cp_var_clk = if ppp.uduc { 0.0001 } else { 0.001 };
                    r_clk[(n_clk, n_clk)] = cp_var_clk / libm::sin(el).max(0.1);
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
        // Record position for backward smoothing (static receivers only)
        if !self.dynamics {
            let pos_upd = Vector3::new(self.x[0], self.x[1], self.x[2]);
            let cov_upd = Matrix3::new(
                self.p[(0,0)], self.p[(0,1)], self.p[(0,2)],
                self.p[(1,0)], self.p[(1,1)], self.p[(1,2)],
                self.p[(2,0)], self.p[(2,1)], self.p[(2,2)],
            );
            self.record_position(pos_upd, cov_upd, pos_pred, cov_pred);
        }
        state.position.vector.x = self.x[0];
        state.position.vector.y = self.x[1];
        state.position.vector.z = self.x[2];
        state.rcv_clk_bias = self.x[ppp.ic(0)];
        state.covariance = self.p.clone();

        // ---- Multi-epoch batch solver: feed CP-derived range ----------------
        // For UDUC mode, derives geometric range from L1 carrier phase using
        // the IEKF's ionosphere and IF ambiguity estimates.  Per-satellite
        // bias elimination in the batch solver absorbs residual errors from
        // WL notching, fractional-cycle IF ambiguity, and ZWD / orbit biases.
        if self.enable_multi_epoch && !self.dynamics && ppp.uduc
            && self.biases_seeded && self.epoch > 5
        {
            let rcv_pos = Vector3::new(self.x[0], self.x[1], self.x[2]);
            let zwd = if ppp.nt() >= 1 { self.x[ppp.it()] } else { 0.0 };
            for i in 0..obs_data.len() {
                let (sat, cp1_cyc, _, _, _, el_deg) = obs_data[i];
                if cp1_cyc == 0.0 { continue; }
                let el_rad = el_deg * D2R;

                // Low-elevation satellites have stronger atmospheric and
                // multipath correlation across epochs.  Gate at 15° standard.
                if el_deg < 15.0 { continue; }

                let lam1 = lam1_vals[i];
                let l1_cp_m = cp1_cyc * lam1;
                let iono_l1 = self.x[ppp.ni(i)];        // IEKF ionosphere estimate
                let n_if = self.x[ppp.ib(i)];            // IEKF IF ambiguity estimate
                let sys: usize = if sat.constellation == Constellation::Glonass { 1 } else { 0 };
                let clock = self.x[ppp.ic(sys)];         // receiver clock
                let dts = sat_clk[i];                    // satellite clock
                let dtrp = proc_tropo_dry[i] + proc_map_wet[i] * zwd;

                // Derive geometric range from L1 carrier phase:
                //   L1_CP = ρ - dts + dtrp - I1 + clock + λ1*N1
                //   ρ = L1_CP + dts - dtrp + I1 - clock - N_IF
                //     (λ1*N1 ≈ N_IF with residual absorbed by per-sat bias)
                //
                // Pack into add_measurement:
                //   cp_if = L1_CP + dts - dtrp + I1_est
                //   n_if_fixed = N_IF_est
                //   clock = receiver_clock
                let cp_with_corrections = l1_cp_m + dts - dtrp + iono_l1;
                self.batch_solver.add_measurement(
                    sat,                  // satellite ID (per-sat bias key)
                    sat_pos[i],           // satellite ECEF position
                    cp_with_corrections,  // L1 CP + dts - dtrp + iono (meters)
                    n_if,                 // IF ambiguity estimate (meters)
                    0.0,                  // tropo (already in cp_with_corrections)
                    0.0, 0.0,             // map_wet, zwd (not used)
                    clock,                // receiver clock (subtracted)
                    el_rad,               // elevation (radians)
                );
            }

            // Periodic reset: flush early, poorly-converged measurements
            // so the solver doesn't lock onto the IEKF's initial bias.
            // 300 epochs ≈ 2.5 hours at 30s intervals — long enough for
            // geometry diversity while preventing stale bias accumulation.
            if self.epoch > 0 && self.epoch % 300 == 0 {
                self.batch_solver.clear();
            }

            // Solve for refined static position.
            if let Some((refined, _cov)) = self.batch_solver.solve(rcv_pos) {
                let delta = (refined - rcv_pos).norm();
                let iepf_sigma = (self.p[(0,0)] + self.p[(1,1)] + self.p[(2,2)]).sqrt();
                if delta < 5.0 || delta < 3.0 * iepf_sigma.max(1.0) {
                    state.position.vector.x = refined.x;
                    state.position.vector.y = refined.y;
                    state.position.vector.z = refined.z;
                }
            }
        }

        // ---- Geometry-based fast AR (epoch 30+) ----
        // Use the converged position to fix N1 directly from geometry,
        // bypassing the 100-sample MW WL accumulation. This triggers AR
        // at epoch ~40-60 instead of ~100, halving the pre-AR tail.
        if !ppp.uduc && self.epoch > 30 {
            let rcv_pos = Vector3::new(self.x[0], self.x[1], self.x[2]);
            for i in 0..obs_data.len() {
                if lc_if_vals[i] == 0.0 { continue; }
                let (sat, _, _, _, _, _) = obs_data[i];
                let sys: usize = if sat.constellation == Constellation::Glonass { 1 } else { 0 };
                let clock = self.x[ppp.ic(sys)];
                let rs = sat_pos[i]; let dts = sat_clk[i];
                let dist = (rs - rcv_pos).norm();
                if dist <= 0.0 { continue; }
                let dtrp = proc_tropo_dry[i];
                let rng = dist - dts + dtrp;
                let lam_nl = 299792458.0 / (f1_vals[i] + f2_vals[i]);
                let n1_est = (lc_if_vals[i] - rng - clock) / lam_nl;
                let pos_sigma = (self.p[(0,0)] + self.p[(1,1)] + self.p[(2,2)]).sqrt();
                // When position is well-known (σ<1m): search candidate N1 values
                // within position uncertainty range. This handles the bootstrap
                // problem where position error exceeds the NL wavelength (0.107m).
                if pos_sigma < 1.0 {
                    // Position is well-known: use wide candidate search without
                    // ratio test. The tight prior prevents wrong fixes from causing
                    // large position excursions. Any candidate within 0.2m of the
                    // float N_IF is accepted (closest wins).
                    let n1_search_radius = ((5.0 * pos_sigma / lam_nl).ceil() as i64).max(3);
                    let n1_rounded = n1_est.round();
                    let f1s = f1_vals[i]*f1_vals[i]; let f2s = f2_vals[i]*f2_vals[i];
                    let n_if_float = self.x[ppp.ib(i)];
                    let mut best_n_if: Option<f64> = None;
                    let mut best_residual = f64::MAX;
                    for dk in -n1_search_radius..=n1_search_radius {
                        let n1_cand = n1_rounded + dk as f64;
                        let n_if_cand = (f1s * n1_cand * lam1_vals[i] - f2s * n1_cand * lam2_vals[i]) / (f1s - f2s);
                        let residual = (n_if_float - n_if_cand).abs();
                        if residual < 0.2 && residual < best_residual {
                            best_n_if = Some(n_if_cand);
                            best_residual = residual;
                        }
                    }
                    if let Some(n_if_fixed) = best_n_if {
                        let bi = ppp.ib(i);
                        self.x[bi] = n_if_fixed;
                        self.p[(bi, bi)] = 0.01;
                    }
                } else {
                    // Standard NL AR: only fix when N1 is close to integer.
                    // This is more conservative but safer with loose position.
                    let n1_rounded = n1_est.round();
                    if (n1_est - n1_rounded).abs() > 0.15 { continue; }
                    let f1s = f1_vals[i]*f1_vals[i]; let f2s = f2_vals[i]*f2_vals[i];
                    let n_if_fixed = (f1s * n1_rounded * lam1_vals[i] - f2s * n1_rounded * lam2_vals[i]) / (f1s - f2s);
                    let n_if_float = self.x[ppp.ib(i)];
                    if (n_if_float - n_if_fixed).abs() > 0.2 { continue; }
                    let bi = ppp.ib(i);
                    self.x[bi] = n_if_fixed;
                    self.p[(bi, bi)] = 0.01;
                }
            }
        }

        // ---- MW widelane tracking for integer AR ----
        // Compute Melbourne-Wübbena widelane for each satellite, maintain EMA.
        // When WL is precise enough, fix N_wl → compute N1 → fix N_IF.
        for i in 0..obs_data.len() {
            let (sat, l1_cyc, l2_cyc, _if_pr, p2_raw, _el_deg) = obs_data[i];
            let p1_raw = raw_p1_vals.get(i).copied().unwrap_or(0.0); // Use raw L1 PR for MW
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
            if count > 100 {
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

                // Tighten IF bias state to the fixed value.
                // Use σ≈10cm (not 1mm hard-lock) so the filter can re-adjust
                // if position drifts. Hard-locking enables a positive feedback
                // loop where wrong fixes compound over many epochs.
                let bi = ppp.ib(i);
                self.x[bi] = n_if_fixed;
                self.p[(bi, bi)] = 0.01; // σ ≈ 10cm — fixed but adjustable
                // Keep cross-correlations: preserves position-bias coupling
                // so the filter can self-correct if the fix was wrong.
            }
        }

        // After AR fixes: CP-only re-solve.
        // Fixed IF ambiguities make CP an unbiased range measurement (σ≈1cm).
        // Do NOT inflate position variance — keeping the current IEKF covariance
        // ensures the Kalman gain distributes CP corrections between position
        // and ambiguity according to their actual uncertainties. Inflating
        // position variance (σ=10m) causes position to absorb all CP correction
        // even when ambiguities are slightly wrong, enabling a positive feedback
        // loop where wrong fixes compound over epochs.
        let any_fixed = (0..obs_data.len()).any(|i| {
            self.p[(ppp.ib(i), ppp.ib(i))] < 0.01
        });
        if any_fixed {
            // Inflate clock/ZWD for CP-only re-convergence (clock is reset each epoch)
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
                r_cp[(n_cp, n_cp)] = if ppp.uduc { 0.0001 } else { 0.001 }; // σ=1cm raw, σ=3cm IF
                n_cp += 1;

                // Feed CP measurement to static multi-epoch batch solver.
                // Only high-elevation AR-fixed satellites — low satellites have
                // stronger atmosphere and multipath correlation across epochs.
                if self.enable_multi_epoch && !self.dynamics && el_rad > 15.0_f64.to_radians() {
                    let zwd_val = if ppp.nt() >= 1 { self.x[ppp.it()] } else { 0.0 };
                    let dtrp = dtrp + proc_map_wet[i] * zwd_val;
                    self.batch_solver.add_measurement(
                        sat,                             // satellite ID for per-sat bias
                        rs,                              // satellite ECEF position
                        lc_if_vals[i],                   // IF carrier phase (meters)
                        n_if_fixed,                      // AR-fixed IF ambiguity
                        dtrp,                            // total tropo delay
                        0.0, 0.0,                        // map_wet, zwd (already in dtrp)
                        self.x[ppp.ic(sys)] - dts,       // clock - sat_clock
                        el_rad,                          // elevation (radians)
                    );
                }
            }
            if n_cp >= 4 {
                // Consistency gate: reject CP-only re-solve if it would
                // jump position by >1m. Large jumps indicate wrong AR fixes
                // that would trigger the divergence feedback loop.
                let pos_before = Vector3::new(self.x[0], self.x[1], self.x[2]);
                let h_s = h_cp.view((0, 0), (ppp.nx(), n_cp)).clone_owned();
                let vs = v_cp.rows(0, n_cp).clone_owned();
                let rs = r_cp.view((0, 0), (n_cp, n_cp)).clone_owned();
                let _ = Self::measurement_update(&mut self.x, &mut self.p, &h_s, &vs, &rs, ppp.nx(), n_cp);
                let pos_after = Vector3::new(self.x[0], self.x[1], self.x[2]);
                let pos_jump = (pos_after - pos_before).norm();
                if pos_jump > 1.0 {
                    tracing::warn!("CP-only re-solve rejected: pos jump {:.2}m > 1m ({} fixes)", pos_jump, n_cp);
                    self.x[0] = pos_before.x;
                    self.x[1] = pos_before.y;
                    self.x[2] = pos_before.z;
                }
            }

            // Periodic reset: prevent stale AR-fixed measurements from locking
            // the position to an old solution as geometry evolves.
            if self.epoch > 0 && self.epoch % 300 == 0 {
                self.batch_solver.clear();
            }

            // Solve for refined static position across all accumulated CP measurements
            if self.enable_multi_epoch && !self.dynamics {
                let n_meas = self.batch_solver.num_measurements();
                let rcv_pos = Vector3::new(self.x[0], self.x[1], self.x[2]);
                if let Some((refined, _cov)) = self.batch_solver.solve(rcv_pos) {
                    let delta = (refined - rcv_pos).norm();
                    let iepf_sigma = (self.p[(0,0)] + self.p[(1,1)] + self.p[(2,2)]).sqrt();
                    if delta < 5.0 || delta < 3.0 * iepf_sigma.max(1.0) {
                        if delta > 0.01 {
                            tracing::debug!(
                                "Batch solver: {} measurements, pos delta={:.3}m",
                                n_meas, delta
                            );
                        }
                        state.position.vector.x = refined.x;
                        state.position.vector.y = refined.y;
                        state.position.vector.z = refined.z;
                    }
                } else if n_meas > 0 && self.epoch % 60 == 0 {
                    tracing::debug!(
                        "Batch solver: {} measurements but solve failed at epoch {}",
                        n_meas, self.epoch
                    );
                }
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

    /// Store forward-pass position data for backward smoothing.
    fn record_position(&mut self, pos: Vector3<f64>, cov: Matrix3<f64>,
                        pos_pred: Vector3<f64>, cov_pred: Matrix3<f64>) {
        self.position_history.push(PositionEntry { pos, cov, pos_pred, cov_pred });
    }

    /// Run position-only RTS backward smoother over recorded history.
    /// Returns smoothed positions (ECEF, meters). Assumes static receiver
    /// (state transition = identity for position).
    pub fn smooth_positions(&self) -> Vec<Vector3<f64>> {
        let n = self.position_history.len();
        if n == 0 { return Vec::new(); }

        let mut smoothed: Vec<Vector3<f64>> = self.position_history.iter()
            .map(|e| e.pos).collect();

        for k in (0..n-1).rev() {
            let entry_k = &self.position_history[k];
            let entry_k1 = &self.position_history[k + 1];

            // C_k = P_k * inv(P_{k+1|k})
            let p_pred_inv = match crate::engine::ppp_multi_epoch_batch::try_invert_3x3(&entry_k1.cov_pred) {
                Some(inv) => inv,
                None => continue,
            };
            let c_k = entry_k.cov * p_pred_inv;

            // x_{k|N} = x_k + C_k * (x_{k+1|N} - x_{k+1|k})
            let dx = smoothed[k + 1] - entry_k1.pos_pred;
            smoothed[k] = entry_k.pos + c_k * dx;
        }

        smoothed
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
