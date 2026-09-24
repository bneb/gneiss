//! Tightly-Coupled Ambiguity Tracker & Integer Conditioning Engine (TC-AR).
use std::collections::HashMap;
use nalgebra::{DMatrix, DVector, SVector, Vector3};

use gneiss_core::coords::{az_el, ecef_to_llh};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::{Constellation, SatelliteId};

use crate::ambiguity::ffrt::calculate_threshold;
use crate::ambiguity::lambda::resolve_lambda;
use crate::ambiguity::par::select_ils_subset;
use crate::estimators::eskf::condition::apply_integer_conditioning;
use crate::estimators::eskf::update::apply_error_injection;
use crate::estimators::eskf::{EskfState, Matrix15};
use crate::estimators::rtk_iekf::ref_sat::is_beidou_geo;
use crate::estimators::rtk_iekf::DoubleDiffKey;

pub const MAX_CARRIER_RESIDUAL_M: f64 = 0.05;
pub const MIN_PAR_SUBSET_SIZE: usize = 4;
pub const MIN_RATIO_FLOOR: f64 = 2.0;

/// Major GNSS constellation groups for isolated double-differencing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ConstellationGroup {
    GpsQzss,
    Galileo,
    Beidou,
    Glonass,
}

impl ConstellationGroup {
    pub const ALL: [ConstellationGroup; 4] = [
        ConstellationGroup::GpsQzss, ConstellationGroup::Galileo,
        ConstellationGroup::Beidou, ConstellationGroup::Glonass,
    ];

    #[inline]
    pub fn matches(&self, sat: SatelliteId) -> bool {
        match self {
            Self::GpsQzss => sat.constellation == Constellation::Gps || sat.constellation == Constellation::Qzss,
            Self::Galileo => sat.constellation == Constellation::Galileo,
            Self::Beidou => sat.constellation == Constellation::Beidou,
            Self::Glonass => sat.constellation == Constellation::Glonass,
        }
    }

    #[inline]
    pub fn constellation_id(&self) -> u8 {
        match self {
            ConstellationGroup::GpsQzss => Constellation::Gps as u8,
            ConstellationGroup::Galileo => Constellation::Galileo as u8,
            ConstellationGroup::Beidou => Constellation::Beidou as u8,
            ConstellationGroup::Glonass => Constellation::Glonass as u8,
        }
    }
}

/// Selected reference satellite for a constellation group.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroupRefSat {
    pub sat: SatelliteId,
    pub pos: Vector3<f64>,
    pub u_ref: Vector3<f64>,
}

fn sat_elevation(ant_pos: &Vector3<f64>, sat_pos: &Vector3<f64>) -> f64 {
    az_el(ecef_to_llh(*ant_pos), *ant_pos, *sat_pos).1
}

fn candidate_elevation(sat: SatelliteId, sat_pos: &Vector3<f64>, ant_pos: &Vector3<f64>) -> f64 {
    let el = sat_elevation(ant_pos, sat_pos);
    if is_beidou_geo(sat) { el - 0.35 } else { el }
}

fn score_candidate(
    r_sat: &SatObs,
    rover_obs: &EpochObs,
    ephems: &[Ephemeris],
    ant_pos: &Vector3<f64>,
    min_el: f64,
) -> Option<(Vector3<f64>, f64)> {
    let eph = ephems.iter().find(|e| e.sat() == r_sat.sat)?;
    let (sat_pos, _, _, _) = eph.position(rover_obs.time);
    let el = sat_elevation(ant_pos, &sat_pos);
    if el < min_el { return None; }
    Some((sat_pos, candidate_elevation(r_sat.sat, &sat_pos, ant_pos)))
}

/// Find the highest-elevation reference satellite for a constellation group.
pub fn find_group_ref_sat(
    group: ConstellationGroup,
    rover_obs: &EpochObs,
    base_obs: &EpochObs,
    ephems: &[Ephemeris],
    ant_pos: &Vector3<f64>,
    min_el: f64,
) -> Option<GroupRefSat> {
    let mut best: Option<(SatelliteId, Vector3<f64>, f64)> = None;
    let mut count = 0usize;
    for r_sat in &rover_obs.satellites {
        if !group.matches(r_sat.sat) || !base_obs.satellites.iter().any(|s| s.sat == r_sat.sat) {
            continue;
        }
        let Some((pos, score)) = score_candidate(r_sat, rover_obs, ephems, ant_pos, min_el) else { continue };
        count += 1;
        if best.as_ref().is_none_or(|b| score > b.2) {
            best = Some((r_sat.sat, pos, score));
        }
    }
    if count < 2 { return None; }
    let (sat, pos, _) = best?;
    let diff = pos - ant_pos;
    let d = diff.norm();
    if d <= 1e-3 { None } else { Some(GroupRefSat { sat, pos, u_ref: diff / d }) }
}

/// Dynamic float ambiguity and state cross-covariance tracker.
#[derive(Debug, Clone)]
pub struct TcAmbiguityTracker {
    pub keys: Vec<DoubleDiffKey>,
    pub a_float: DVector<f64>,
    pub q_aa: DMatrix<f64>,
    pub p_xa: DMatrix<f64>,
    pub fixed_integers: HashMap<DoubleDiffKey, i32>,
    pub lock_counts: HashMap<DoubleDiffKey, u32>,
}

impl TcAmbiguityTracker {
    pub fn new() -> Self {
        Self {
            keys: Vec::new(),
            a_float: DVector::zeros(0),
            q_aa: DMatrix::zeros(0, 0),
            p_xa: DMatrix::zeros(15, 0),
            fixed_integers: HashMap::new(),
            lock_counts: HashMap::new(),
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    pub fn key_index(&self, key: &DoubleDiffKey) -> Option<usize> {
        self.keys.iter().position(|k| k == key)
    }

    /// Propagate state-ambiguity cross-covariance: P_xa <- Phi_15 * P_xa.
    pub fn propagate_cross_cov(&mut self, phi: &Matrix15<f64>) {
        if self.keys.is_empty() { return; }
        for j in 0..self.keys.len() {
            let col = phi * self.p_xa.column(j);
            self.p_xa.column_mut(j).copy_from(&col);
            self.q_aa[(j, j)] += 1e-7;
        }
    }

    /// Synchronize tracked keys with active double-difference pairs.
    pub fn sync_keys(&mut self, new_keys: &[DoubleDiffKey], init_floats: &[f64]) {
        let (a, q, p) = rebuild_tracker_matrices(
            &self.keys, new_keys, &self.a_float, &self.q_aa, &self.p_xa, init_floats,
        );
        self.a_float = a;
        self.q_aa = q;
        self.p_xa = p;
        self.keys = new_keys.to_vec();
        for k in new_keys {
            let count = self.lock_counts.entry(*k).or_insert(0);
            *count = count.saturating_add(1);
        }
        self.fixed_integers.retain(|k, _| new_keys.contains(k));
        self.lock_counts.retain(|k, _| new_keys.contains(k));
    }

    /// Update ESKF state, P_xx, P_xa, and Q_aa with a float carrier observation.
    pub fn update_carrier_float(
        &mut self,
        state: &mut EskfState,
        idx: usize,
        h_x: &SVector<f64, 15>,
        wl: f64,
        y: f64,
        r_var: f64,
    ) {
        let col = self.p_xa.column(idx);
        let p_hx = state.cov * h_x;
        let h_c = h_x.dot(&col);
        let s = h_x.dot(&p_hx) + 2.0 * wl * h_c + wl * wl * self.q_aa[(idx, idx)] + r_var;
        if s <= 1e-12 { return; }
        let k_x = (p_hx + wl * col) / s;
        let q_col = self.q_aa.column(idx);
        let k_a = (self.p_xa.transpose() * h_x + wl * q_col) / s;

        let dx = k_x * y;
        apply_error_injection(state, &dx);
        self.a_float += &k_a * y;
        state.cov -= s * (k_x * k_x.transpose());
        self.p_xa -= s * (k_x * k_a.transpose());
        self.q_aa -= s * (&k_a * k_a.transpose());
        symmetrize_and_floor(&mut state.cov, &mut self.q_aa);
    }

    /// Update ESKF state, P_xx, P_xa, and Q_aa with a float pseudorange observation.
    pub fn update_code_float(
        &mut self,
        state: &mut EskfState,
        h_x: &SVector<f64, 15>,
        y: f64,
        r_var: f64,
    ) {
        let p_hx = state.cov * h_x;
        let s = h_x.dot(&p_hx) + r_var;
        if s <= 1e-12 { return; }
        let k_x = p_hx / s;
        let k_a = (self.p_xa.transpose() * h_x) / s;

        let dx = k_x * y;
        apply_error_injection(state, &dx);
        self.a_float += &k_a * y;
        state.cov -= s * (k_x * k_x.transpose());
        self.p_xa -= s * (k_x * k_a.transpose());
        self.q_aa -= s * (&k_a * k_a.transpose());
        symmetrize_and_floor(&mut state.cov, &mut self.q_aa);
    }

    /// Extract sub-matrix of P_xa for given ambiguity indices.
    pub fn submatrix_p_xa(&self, indices: &[usize]) -> DMatrix<f64> {
        let mut sub = DMatrix::zeros(15, indices.len());
        for (col_idx, &orig_idx) in indices.iter().enumerate() {
            sub.column_mut(col_idx).copy_from(&self.p_xa.column(orig_idx));
        }
        sub
    }

    /// Extract sub-matrix of Q_aa for given ambiguity indices.
    pub fn submatrix_q_aa(&self, indices: &[usize]) -> DMatrix<f64> {
        let k = indices.len();
        let mut sub = DMatrix::zeros(k, k);
        for (r, &orig_r) in indices.iter().enumerate() {
            for (c, &orig_c) in indices.iter().enumerate() {
                sub[(r, c)] = self.q_aa[(orig_r, orig_c)];
            }
        }
        sub
    }

    /// Extract sub-vector of float ambiguities for given indices.
    pub fn subvector_a_float(&self, indices: &[usize]) -> DVector<f64> {
        DVector::from_iterator(indices.len(), indices.iter().map(|&i| self.a_float[i]))
    }

    /// Attempt full-set or PAR integer AR, apply conditioning, and screen residuals.
    pub fn attempt_ar_and_condition<F>(
        &mut self,
        state: &mut EskfState,
        check_residual: F,
    ) -> (bool, Option<f64>)
    where
        F: Fn(&Vector3<f64>, DoubleDiffKey, i32) -> Option<f64>,
    {
        let n = self.keys.len();
        if n < MIN_PAR_SUBSET_SIZE { return (false, None); }
        let (indices, integers, ratio) = match self.solve_integers_full_or_par() {
            Some(res) => res,
            None => return (false, None),
        };
        let backup = state.clone();
        if !self.apply_and_screen(state, &indices, &integers, &check_residual) {
            *state = backup;
            return (false, Some(ratio));
        }
        self.record_fixed_integers(&indices, &integers);
        (true, Some(ratio))
    }

    fn solve_integers_full_or_par(&self) -> Option<(Vec<usize>, Vec<i32>, f64)> {
        let n = self.keys.len();
        let full_thresh = calculate_threshold(n, 0.001).max(MIN_RATIO_FLOOR);
        if let Ok(res) = resolve_lambda(&self.a_float, &self.q_aa) {
            if res.ratio >= full_thresh {
                let ints: Vec<i32> = res.best_integers.iter().map(|&x| x as i32).collect();
                return Some(((0..n).collect(), ints, res.ratio));
            }
        }
        let (sub_idx, sub_a, sub_q) = select_ils_subset(&self.a_float, &self.q_aa, 0.995);
        if sub_idx.len() < MIN_PAR_SUBSET_SIZE || sub_idx.len() >= n { return None; }
        let par_thresh = calculate_threshold(sub_idx.len(), 0.001).max(MIN_RATIO_FLOOR);
        let res = resolve_lambda(&sub_a, &sub_q).ok()?;
        if res.ratio >= par_thresh {
            let ints: Vec<i32> = res.best_integers.iter().map(|&x| x as i32).collect();
            Some((sub_idx, ints, res.ratio))
        } else {
            None
        }
    }

    fn apply_and_screen<F>(
        &self,
        state: &mut EskfState,
        indices: &[usize],
        integers: &[i32],
        check_residual: &F,
    ) -> bool
    where
        F: Fn(&Vector3<f64>, DoubleDiffKey, i32) -> Option<f64>,
    {
        let sub_p = self.submatrix_p_xa(indices);
        let sub_q = self.submatrix_q_aa(indices);
        let sub_af = self.subvector_a_float(indices);
        let sub_ax = DVector::from_iterator(integers.len(), integers.iter().map(|&x| x as f64));

        let cond = match apply_integer_conditioning(state, &sub_p, &sub_q, &sub_af, &sub_ax) {
            Ok(c) => c,
            Err(_) => return false,
        };
        if !cond.applied || !cond.accepted { return false; }
        let ant_pos = state.pos_ecef;
        for (i, &orig_idx) in indices.iter().enumerate() {
            let key = self.keys[orig_idx];
            let res = match check_residual(&ant_pos, key, integers[i]) {
                Some(r) => r,
                None => continue,
            };
            if res.abs() > MAX_CARRIER_RESIDUAL_M { return false; }
        }
        true
    }

    fn record_fixed_integers(&mut self, indices: &[usize], integers: &[i32]) {
        for (i, &idx) in indices.iter().enumerate() {
            let key = self.keys[idx];
            self.fixed_integers.insert(key, integers[i]);
            self.p_xa.column_mut(idx).fill(0.0);
            self.q_aa.column_mut(idx).fill(0.0);
            self.q_aa.row_mut(idx).fill(0.0);
            self.q_aa[(idx, idx)] = 1e-4;
        }
    }
}

impl Default for TcAmbiguityTracker {
    fn default() -> Self {
        Self::new()
    }
}

fn copy_cov_row(
    i: usize,
    oi: usize,
    old_keys: &[DoubleDiffKey],
    new_keys: &[DoubleDiffKey],
    old_q: &DMatrix<f64>,
    q: &mut DMatrix<f64>,
) {
    for (j, mk) in new_keys.iter().enumerate() {
        if let Some(oj) = old_keys.iter().position(|k| k == mk) {
            q[(i, j)] = old_q[(oi, oj)];
        }
    }
}

fn rebuild_tracker_matrices(
    old_keys: &[DoubleDiffKey],
    new_keys: &[DoubleDiffKey],
    old_a: &DVector<f64>,
    old_q: &DMatrix<f64>,
    old_p: &DMatrix<f64>,
    init_floats: &[f64],
) -> (DVector<f64>, DMatrix<f64>, DMatrix<f64>) {
    let m = new_keys.len();
    let mut a = DVector::zeros(m);
    let mut q = DMatrix::zeros(m, m);
    let mut p = DMatrix::zeros(15, m);

    for (i, nk) in new_keys.iter().enumerate() {
        if let Some(oi) = old_keys.iter().position(|k| k == nk) {
            a[i] = old_a[oi];
            p.column_mut(i).copy_from(&old_p.column(oi));
            copy_cov_row(i, oi, old_keys, new_keys, old_q, &mut q);
        } else {
            a[i] = init_floats[i];
            q[(i, i)] = 100.0;
        }
    }
    (a, q, p)
}

fn symmetrize_and_floor(cov15: &mut Matrix15<f64>, q_aa: &mut DMatrix<f64>) {
    *cov15 = 0.5 * (*cov15 + cov15.transpose());
    for i in 0..15 { if cov15[(i, i)] < 1e-10 { cov15[(i, i)] = 1e-10; } }
    *q_aa = 0.5 * (&*q_aa + q_aa.transpose());
    for i in 0..q_aa.nrows() { if q_aa[(i, i)] < 1e-6 { q_aa[(i, i)] = 1e-6; } }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::UnitQuaternion;

    #[test]
    fn test_constellation_group_isolation_no_cross_matching() {
        let gps1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let qzss2 = SatelliteId { constellation: Constellation::Qzss, prn: 2 };
        let gal5 = SatelliteId { constellation: Constellation::Galileo, prn: 5 };
        let bds3 = SatelliteId { constellation: Constellation::Beidou, prn: 3 };
        let glo7 = SatelliteId { constellation: Constellation::Glonass, prn: 7 };

        assert!(ConstellationGroup::GpsQzss.matches(gps1) && ConstellationGroup::GpsQzss.matches(qzss2));
        assert!(!ConstellationGroup::GpsQzss.matches(gal5));
        assert!(ConstellationGroup::Galileo.matches(gal5) && !ConstellationGroup::Galileo.matches(bds3));
        assert!(ConstellationGroup::Beidou.matches(bds3) && !ConstellationGroup::Beidou.matches(glo7));
        assert!(ConstellationGroup::Glonass.matches(glo7) && !ConstellationGroup::Glonass.matches(gps1));
    }

    #[test]
    fn test_cross_covariance_propagation_multiplication() {
        let mut tracker = TcAmbiguityTracker::new();
        let key = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        tracker.sync_keys(&[key], &[5.0]);

        tracker.p_xa[(0, 0)] = 2.0;
        tracker.p_xa[(3, 0)] = 1.0;

        let mut phi = Matrix15::identity();
        phi[(0, 3)] = 0.1; // dx = dx + 0.1 * dv -> P_xa[0] = 2.0 + 0.1 * 1.0 = 2.1
        tracker.propagate_cross_cov(&phi);

        assert!((tracker.p_xa[(0, 0)] - 2.1).abs() < 1e-12);
        assert!((tracker.p_xa[(3, 0)] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_float_carrier_and_code_update_contract_covariances() {
        let mut tracker = TcAmbiguityTracker::new();
        let key = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };
        tracker.sync_keys(&[key], &[10.0]);

        let mut state = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
        let mut h_x = SVector::<f64, 15>::zeros();
        h_x[0] = 1.0;

        let prior_var = state.cov[(0, 0)];
        tracker.update_code_float(&mut state, &h_x, 0.20, 0.09);
        assert!(state.cov[(0, 0)] < prior_var);

        let prior_amb_var = tracker.q_aa[(0, 0)];
        tracker.update_carrier_float(&mut state, 0, &h_x, 0.1903, 0.01, 1e-4);
        assert!(tracker.q_aa[(0, 0)] < prior_amb_var);
    }

    #[test]
    fn test_carrier_residual_screening_rejects_large_residual() {
        let mut tracker = TcAmbiguityTracker::new();
        let keys: Vec<DoubleDiffKey> = (2..=6).map(|prn| DoubleDiffKey {
            constellation_id: 0, sat: prn, ref_sat: 1, freq_band: 1,
        }).collect();
        tracker.sync_keys(&keys, &[1.0, 2.0, 3.0, 4.0, 5.0]);

        let mut state = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
        // Residual checker that returns an unphysical blunder (> 0.05m)
        let (fixed, _) = tracker.attempt_ar_and_condition(&mut state, |_, _, _| Some(0.12));
        assert!(!fixed);
    }

    #[test]
    fn test_tracker_post_fix_q_aa_preserves_positive_definiteness() {
        let mut tracker = TcAmbiguityTracker::new();
        let keys: Vec<DoubleDiffKey> = (2..=5).map(|sat| DoubleDiffKey {
            constellation_id: 0, sat, ref_sat: 1, freq_band: 1,
        }).collect();
        tracker.sync_keys(&keys, &[0.0, 0.0, 0.0, 0.0]);
        for i in 0..4 {
            for j in 0..4 { tracker.q_aa[(i, j)] = if i == j { 0.08 } else { 0.04 }; }
        }
        let mut state = EskfState::new(Vector3::zeros(), Vector3::zeros(), UnitQuaternion::identity());
        let (fixed, _) = tracker.attempt_ar_and_condition(&mut state, |_, _, _| Some(0.001));
        assert!(fixed);
        let eig = nalgebra::linalg::SymmetricEigen::new(tracker.q_aa.clone()).eigenvalues.min();
        assert!(eig > 0.0, "Q_aa indefinite: {eig}");
    }
}
