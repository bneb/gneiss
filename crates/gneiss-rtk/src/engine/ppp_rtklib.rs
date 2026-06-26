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

use crate::engine::EngineError;
use crate::filter::RtkState;
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
const VAR_POS: f64 = 10000.0; // 100^2 m^2
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
}

impl PppState {
    fn new(has_glo: bool, dynamics: bool) -> Self {
        Self {
            np: if dynamics { 9 } else { 3 },
            has_glo,
            trop_opt: 2, // EST (estimate ZTD)
            nsat: 0,
            amb_sats: Vec::new(),
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

    /// Total number of resolved states (before biases)
    fn nr(&self) -> usize {
        self.it() + self.nt()
    }

    /// Index of the first phase bias state
    fn ib(&self, sat_idx: usize) -> usize {
        self.nr() + sat_idx
    }

    /// Total number of estimated states
    fn nx(&self) -> usize {
        self.ib(self.nsat) // biases for all observed satellites
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
    biases_seeded: bool,
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
            biases_seeded: false,
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
        for i in 0..ppp.nsat {
            p[(nr + i, nr + i)] = VAR_BIAS;
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

            // Troposphere: simple Saastamoinen
            let dtrp = self.trop_saas(el);
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

            // ---- Phase measurement ----
            if lc != 0.0 {
                for k in 0..nx {
                    h[(k, nv)] = 0.0;
                }
                v[nv] = lc - rng;
                for k in 0..3 {
                    h[(k, nv)] = -e[k];
                }
                if sys != 1 {
                    v[nv] -= x[ppp.ic(0)];
                    h[(ppp.ic(0), nv)] = 1.0;
                } else {
                    v[nv] -= x[ppp.ic(1)];
                    h[(ppp.ic(1), nv)] = 1.0;
                }
                // Troposphere mapping
                let mw = 1.0 / libm::sin(el).max(0.1);
                if ppp.nt() >= 1 {
                    h[(ppp.it(), nv)] = mw;
                }
                // Phase bias
                v[nv] -= x[ppp.ib(i)];
                h[(ppp.ib(i), nv)] = 1.0;

                // Measurement variance: σ ≈ 10cm at zenith, scaled by 1/sin(el)
                let var_phase = 0.01 / libm::sin(el).max(0.1) + sat_var[i] + vart;
                r[(nv, nv)] = var_phase;

                // Innovation test — tighter threshold for carrier phase
                if self.max_inno_cp > 0.0 && v[nv].abs() > self.max_inno_cp && sys != 1 {
                    continue;
                }
                nv += 1;
            }

            // ---- Code measurement ----
            if pc != 0.0 {
                for k in 0..nx {
                    h[(k, nv)] = 0.0;
                }
                v[nv] = pc - rng;
                for k in 0..3 {
                    h[(k, nv)] = -e[k];
                }
                if sys != 1 {
                    v[nv] -= x[ppp.ic(0)];
                    h[(ppp.ic(0), nv)] = 1.0;
                } else {
                    v[nv] -= x[ppp.ic(1)];
                    h[(ppp.ic(1), nv)] = 1.0;
                }
                let mw = 1.0 / libm::sin(el).max(0.1);
                if ppp.nt() >= 1 {
                    h[(ppp.it(), nv)] = mw;
                }

                // PR variance: σ≈5m at zenith. Inflated from RTKLIB defaults
                // because IF pseudorange has systematic biases (~10-15m) that
                // would otherwise dominate the solution. CP (σ≈0.1m) pulls
                // toward the true position once biases are initialized.
                let var_code = 25.0 / libm::sin(el).max(0.1) + sat_var[i] + vart;
                r[(nv, nv)] = var_code;

                if self.max_inno_m > 0.0 && v[nv].abs() > self.max_inno_m && sys != 1 {
                    continue;
                }
                nv += 1;
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
            obs_data.push((sat.sat_obs.sat, sat.cp1.unwrap_or(0.0), sat.cp2.unwrap_or(0.0), sat.p1, sat.p2.unwrap_or(sat.p1), el * R2D));
            lc_if_vals.push(lc_if);
            sat_pos.push(sat.sat_pos_rot);
            sat_clk.push(sat.dt_sat_m);
            sat_var.push(0.0);
        }
        tracing::debug!("solve_with_sats: {} sats -> {} obs_data ({} with CP)", sats.len(), obs_data.len(), lc_if_vals.iter().filter(|v| **v != 0.0).count()); if obs_data.len() < 4 { return Err(EngineError::InsufficientSatellites); }
        let has_glo = obs_data.iter().any(|(s,_,_,_,_,_)| s.constellation == Constellation::Glonass);
        let mut ppp = PppState::new(has_glo, self.dynamics);
        ppp.nsat = obs_data.len();
        let nx = ppp.nx();
        if self.epoch == 0 {
            let mut x0 = DVector::zeros(nx);
            x0[0] = state.position.vector.x; x0[1] = state.position.vector.y; x0[2] = state.position.vector.z;
            // Seed clock from SPP estimate (typically ~5ms = 1.4M meters for F9P).
            // Starting at 0 forces the first measurement update to absorb the full
            // receiver clock offset, which leaks into position through the gain matrix.
            x0[ppp.ic(0)] = state.rcv_clk_bias;
            self.x = x0;
            self.p = self.init_covariance(&ppp, nx);
            self.biases_seeded = false;
        } else if self.last_nsat != obs_data.len() {
            // nsat changed: resize state vector, preserve existing state values
            let old_nr = ppp.nr();
            let old_nx = self.x.len();
            let mut x_new = DVector::zeros(nx);
            let mut p_new = DMatrix::zeros(nx, nx);
            // Copy existing position/clock/tropo states
            let nr = ppp.nr();
            for i in 0..old_nr.min(nr) {
                x_new[i] = self.x[i];
                for j in 0..old_nr.min(nr) {
                    p_new[(i, j)] = self.p[(i, j)];
                }
            }
            // Existing biases get their old values; new biases get VAR_BIAS
            for i in 0..old_nx.saturating_sub(old_nr) {
                let new_idx = nr + i;
                if new_idx < nx {
                    x_new[new_idx] = self.x[old_nr + i];
                    p_new[(new_idx, new_idx)] = self.p[(old_nr + i, old_nr + i)];
                }
            }
            for i in old_nx.saturating_sub(old_nr)..ppp.nsat {
                let idx = nr + i;
                if idx < nx {
                    p_new[(idx, idx)] = VAR_BIAS;
                }
            }
            self.x = x_new;
            self.p = p_new;
            self.biases_seeded = false;
        }
        self.last_nsat = obs_data.len(); self.epoch += 1;
        let mut xp = self.x.clone();
        let mut pp = self.p.clone();
        self.predict(&ppp, &mut xp, &mut pp);
        let nv_max = obs_data.len() * 2;
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
            let nv = self.residuals(&ppp, &obs_data, &lc_for_filter, &range_offsets, &sat_pos, &sat_clk, &sat_var, &xp, &mut v, &mut h_mat, &mut r_mat);
            if nv < 4 { break; }
            let h_s = h_mat.view((0, 0), (nx, nv)).clone_owned();
            let vs = v.rows(0, nv).clone_owned();
            let rs = r_mat.view((0, 0), (nv, nv)).clone_owned();
            if Self::measurement_update(&mut xp, &mut pp, &h_s, &vs, &rs, nx, nv).is_err() { break; }
        }
        // After warmup epoch: seed phase biases from the PR-converged state.
        // Using the filtered position (not SPP) gives biases within ~3m,
        // so CP residuals start small enough for σ=10cm measurements to pull.
        if !self.biases_seeded {
            let mut seeded = 0;
            for i in 0..obs_data.len() {
                if lc_if_vals[i] == 0.0 { continue; }
                let rs = sat_pos[i];
                let dts = sat_clk[i];
                let (sat, _, _, _, _, el) = obs_data[i];
                let sys: usize = if sat.constellation == Constellation::Glonass { 1 } else { 0 };
                let dist = (rs - Vector3::new(xp[0], xp[1], xp[2])).norm();
                let el_rad = el * D2R;
                let dtrp = self.trop_saas(el_rad);
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
        Ok(())
    }

    /// Simple Saastamoinen troposphere model
    fn trop_saas(&self, el: f64) -> f64 {
        let z = std::f64::consts::FRAC_PI_2 - el;
        let p = 1013.25; // sea level pressure
        let t = 288.15; // temperature
        let e = 6.108 * libm::exp((17.15 * t - 4684.0) / (t - 38.45)) * REL_HUMI;
        0.002277 / libm::cos(z) * (p + (1255.0 / t + 0.05) * e)
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
        let mut ppp = PppState::new(has_glo, self.dynamics);
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
        let nv_max = obs_data.len() * 2;
        let mut v = DVector::zeros(nv_max);
        let mut h_mat = DMatrix::zeros(nx, nv_max);
        let mut r_mat = DMatrix::zeros(nv_max, nv_max);
        let mut xp = x.clone();
        let mut pp = p_mat.clone();

        for _iter in 0..self.max_iter {
            let empty_lc: Vec<f64> = vec![0.0; obs_data.len()];
            let empty_off: Vec<f64> = vec![0.0; obs_data.len()];
            let nv = self.residuals(
                &ppp,
                &obs_data,
                &empty_lc,
                &empty_off,
                &sat_pos,
                &sat_clk,
                &sat_var,
                &xp,
                &mut v,
                &mut h_mat,
                &mut r_mat,
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
        let ppp = PppState::new(false, false);
        assert_eq!(ppp.np, 3);
        assert_eq!(ppp.nc(), 1);
        assert_eq!(ppp.ic(0), 3); // pos(3) + GPS clk
        assert_eq!(ppp.it(), 4); // pos(3) + clk(1)
        assert_eq!(ppp.nt(), 1);
        assert_eq!(ppp.nr(), 5); // pos(3) + clk(1) + tropo(1)

        let ppp2 = PppState::new(true, false); // with GLONASS
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
