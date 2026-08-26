//! RTK State Vector & Covariance Management for DD-IEKF.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use gneiss_core::time::GpsTime;

/// Key uniquely identifying a double-difference satellite pair on a frequency band.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DoubleDiffKey {
    pub constellation_id: u8,
    pub sat: u16,
    pub ref_sat: u16,
    pub freq_band: u8,
}

/// Dynamic state vector for Double-Difference IEKF.
///
/// State vector layout:
///   x[0..3] — Position ECEF (m)
///   x[3..6] — Velocity ECEF (m/s)
///   x[6..]  — Double-difference carrier phase ambiguities (cycles)
///   x[last] — Rover ZWD residual (m), only when enabled (long-baseline
///     mode): non-dispersive wet-delay mis-modelling, random walk.
#[derive(Debug, Clone)]
pub struct RtkState {
    pub time: GpsTime,
    pub pos_ecef: Vector3<f64>,
    pub vel_ecef: Vector3<f64>,
    pub ambiguities: Vec<(DoubleDiffKey, f64)>,
    pub zwd_enabled: bool,
    pub zwd_m: f64,
    /// Horizontal wet-delay gradients [north, east] (metres of slant delay
    /// at m_grad(el)=cot(el)). Enabled together as a pair; captures
    /// azimuth-dependent atmospheric systematics a scalar zenith cannot.
    pub grad_enabled: bool,
    pub grad_n_m: f64,
    pub grad_e_m: f64,
    /// Per-DD-pair ionosphere residual states (metres of L1 slant delay).
    /// DEPRECATED for production use — rank-deficient against ambiguities
    /// (see NETWORK_RTK_NEXT_STEPS.md). Kept behind its legacy gate.
    pub iono_enabled: bool,
    pub ionos: Vec<(DoubleDiffKey, f64)>,
    /// Per-SATELLITE slant iono states (metres at L1). A DD pair's H row
    /// uses I_sat − I_ref so the shared reference couples all satellites:
    /// rank-deficiency broken. Replaces per-pair ionos when enabled.
    pub sat_iono_enabled: bool,
    pub sat_ionos: Vec<((u8, u16), f64)>,
    pub cov: DMatrix<f64>,
}

impl RtkState {
    /// Create a new RTK state initialized at a seed position.
    pub fn new(initial_pos: Vector3<f64>, initial_time: GpsTime) -> Self {
        let dim = 6;
        let mut cov = DMatrix::zeros(dim, dim);
        for i in 0..3 {
            cov[(i, i)] = 100.0 * 100.0; // 100m initial position uncertainty
        }
        for i in 3..6 {
            cov[(i, i)] = 10.0 * 10.0; // 10 m/s initial velocity uncertainty
        }
        Self {
            time: initial_time,
            pos_ecef: initial_pos,
            vel_ecef: Vector3::zeros(),
            ambiguities: Vec::new(),
            zwd_enabled: false,
            zwd_m: 0.0,
            grad_enabled: false,
            grad_n_m: 0.0,
            grad_e_m: 0.0,
            iono_enabled: false,
            ionos: Vec::new(),
            sat_iono_enabled: false,
            sat_ionos: Vec::new(),
            cov,
        }
    }

    /// Total dimension of the state vector.
    pub fn dim(&self) -> usize {
        6 + self.ambiguities.len()
            + self.zwd_enabled as usize
            + 2 * self.grad_enabled as usize
            + self.ionos.len()
            + self.sat_ionos.len()
    }

    /// Column offset of first per-satellite iono state.
    pub fn sat_iono_offset(&self) -> usize {
        self.iono_offset() + self.ionos.len()
    }

    /// Column index of a satellite's slant-iono state (for constellation and PRN).
    pub fn get_sat_iono_key_idx(&self, cid: u8, prn: u16) -> Option<usize> {
        self.sat_ionos.iter().position(|(k, _)| *k == (cid, prn))
            .map(|i| self.sat_iono_offset() + i)
    }

    /// Column index of a satellite's slant-iono state (defaults to primary GPS/Galileo constellation cid=0).
    pub fn get_sat_iono_idx(&self, prn: u16) -> Option<usize> {
        self.get_sat_iono_key_idx(0, prn)
    }

    /// Ensure a satellite's slant-iono state exists for specific constellation.
    pub fn ensure_sat_iono_key(&mut self, cid: u8, prn: u16) {
        if !self.sat_iono_enabled || self.get_sat_iono_key_idx(cid, prn).is_some() {
            return;
        }
        self.sat_ionos.push(((cid, prn), 0.0));
        let old_dim = self.cov.nrows();
        let new_dim = old_dim + 1;
        let mut nc = DMatrix::zeros(new_dim, new_dim);
        nc.view_range_mut(0..old_dim, 0..old_dim).copy_from(&self.cov);
        nc[(old_dim, old_dim)] = 4.0;
        self.cov = nc;
    }

    /// Ensure a satellite's slant-iono state exists.
    pub fn ensure_sat_iono(&mut self, prn: u16) {
        self.ensure_sat_iono_key(0, prn);
    }

    /// Column offset of the first iono state (after all ambiguities).
    pub fn iono_offset(&self) -> usize {
        self.amb_offset() + self.ambiguities.len()
    }

    /// Column index for a specific pair's iono state.
    pub fn get_iono_idx(&self, key: &DoubleDiffKey) -> Option<usize> {
        self.ionos.iter().position(|(k, _)| k == key)
            .map(|idx| self.iono_offset() + idx)
    }

    /// Ensure an iono state exists for this pair.
    pub fn ensure_iono(&mut self, key: DoubleDiffKey, initial_var: f64) {
        if !self.iono_enabled || self.get_iono_idx(&key).is_some() {
            return;
        }
        self.ionos.push((key, 0.0));
        let old_dim = self.cov.nrows();
        let new_dim = old_dim + 1;
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        new_cov.view_range_mut(0..old_dim, 0..old_dim).copy_from(&self.cov);
        new_cov[(old_dim, old_dim)] = initial_var.max(0.01);
        self.cov = new_cov;
    }

    /// Column offset applied to ambiguity indices.
    /// Layouts: legacy [pos,vel,ambs]; +ZWD [.., zwd@6, ambs@7..];
    /// +grad [.., zwd@6, gN@7, gE@8, ambs@9..].
    pub fn amb_offset(&self) -> usize {
        6 + self.zwd_enabled as usize + 2 * self.grad_enabled as usize
    }

    /// Column index of the ZWD state, when present.
    pub fn zwd_idx(&self) -> Option<usize> {
        self.zwd_enabled.then_some(6)
    }

    /// Column indices of the gradient states (north, east), when present.
    pub fn grad_idx(&self) -> Option<(usize, usize)> {
        if !self.grad_enabled {
            return None;
        }
        let base = 6 + self.zwd_enabled as usize;
        Some((base, base + 1))
    }

    /// Enable the rover ZWD residual state. Must be called before any
    /// ambiguities exist (i.e., at engine construction): the column is
    /// reserved at index 6 and every ambiguity shifts behind it.
    pub fn enable_zwd(&mut self, init_var: f64) {
        if self.zwd_enabled || !self.ambiguities.is_empty() {
            return;
        }
        self.zwd_enabled = true;
        let old_dim = self.cov.nrows();
        let new_dim = old_dim + 1;
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        new_cov.view_range_mut(0..old_dim, 0..old_dim).copy_from(&self.cov);
        new_cov[(old_dim, old_dim)] = init_var.max(1e-6);
        self.cov = new_cov;
    }

    /// Enable the horizontal tropospheric gradient states [north, east].
    /// Same construction-time constraint as [`enable_zwd`]: call before
    /// ambiguities exist. Reserves two columns immediately after ZWD (or at
    /// index 6 if ZWD is absent).
    pub fn enable_gradients(&mut self, init_var_each: f64) {
        if self.grad_enabled || !self.ambiguities.is_empty() {
            return;
        }
        self.grad_enabled = true;
        let old_dim = self.cov.nrows();
        let new_dim = old_dim + 2;
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        new_cov.view_range_mut(0..old_dim, 0..old_dim).copy_from(&self.cov);
        for k in old_dim..new_dim {
            new_cov[(k, k)] = init_var_each.max(1e-8);
        }
        self.cov = new_cov;
    }

    /// Find index of ambiguity in the state vector.
    pub fn get_amb_idx(&self, key: &DoubleDiffKey) -> Option<usize> {
        self.ambiguities.iter().position(|(k, _)| k == key)
            .map(|idx| self.amb_offset() + idx)
    }

    /// Pack current state into a flat DVector.
    pub fn to_dvector(&self) -> DVector<f64> {
        let mut vec = DVector::zeros(self.dim());
        vec[0] = self.pos_ecef.x;
        vec[1] = self.pos_ecef.y;
        vec[2] = self.pos_ecef.z;
        vec[3] = self.vel_ecef.x;
        vec[4] = self.vel_ecef.y;
        vec[5] = self.vel_ecef.z;
        let off = self.amb_offset();
        for (i, (_, val)) in self.ambiguities.iter().enumerate() {
            vec[off + i] = *val;
        }
        if let Some(i) = self.zwd_idx() {
            vec[i] = self.zwd_m;
        }
        if let Some((gn, ge)) = self.grad_idx() {
            vec[gn] = self.grad_n_m;
            vec[ge] = self.grad_e_m;
        }
        let io = self.iono_offset();
        for (i, (_, val)) in self.ionos.iter().enumerate() {
            vec[io + i] = *val;
        }
        let so = self.sat_iono_offset();
        for (i, (_, val)) in self.sat_ionos.iter().enumerate() {
            vec[so + i] = *val;
        }
        vec
    }

    /// Unpack a flat DVector into the structured state fields.
    pub fn update_from_dvector(&mut self, vec: &DVector<f64>) {
        self.pos_ecef = Vector3::new(vec[0], vec[1], vec[2]);
        self.vel_ecef = Vector3::new(vec[3], vec[4], vec[5]);
        let off = self.amb_offset();
        for (i, (_, val)) in self.ambiguities.iter_mut().enumerate() {
            *val = vec[off + i];
        }
        if let Some(i) = self.zwd_idx() {
            self.zwd_m = vec[i];
        }
        if let Some((gn, ge)) = self.grad_idx() {
            self.grad_n_m = vec[gn];
            self.grad_e_m = vec[ge];
        }
        let io = self.iono_offset();
        for (i, (_, val)) in self.ionos.iter_mut().enumerate() {
            *val = vec[io + i];
        }
        let so = self.sat_iono_offset();
        for (i, (_, val)) in self.sat_ionos.iter_mut().enumerate() {
            *val = vec[so + i];
        }
    }

    /// Ensure an ambiguity state exists, adding it with given seed variance if new.
    pub fn ensure_ambiguity(&mut self, key: DoubleDiffKey, initial_val: f64, initial_var: f64) {
        if self.get_amb_idx(&key).is_some() {
            return;
        }
        self.ambiguities.push((key, initial_val));
        let old_dim = self.cov.nrows();
        let new_dim = old_dim + 1;
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        new_cov.view_range_mut(0..old_dim, 0..old_dim).copy_from(&self.cov);
        new_cov[(old_dim, old_dim)] = initial_var.max(1.0);
        self.cov = new_cov;
    }

    /// Reset an ambiguity variance (e.g. after detected cycle slip).
    pub fn reset_ambiguity(&mut self, key: &DoubleDiffKey, initial_val: f64, initial_var: f64) {
        if let Some(idx) = self.get_amb_idx(key) {
            let rel = idx - self.amb_offset();
            self.ambiguities[rel].1 = initial_val;
            // Clear cross-covariances for this ambiguity
            for r in 0..self.cov.nrows() {
                self.cov[(r, idx)] = 0.0;
                self.cov[(idx, r)] = 0.0;
            }
            self.cov[(idx, idx)] = initial_var.max(1.0);
        }
    }

    /// Remove stale ambiguities not in the active set.
    pub fn retain_active_ambiguities(&mut self, active_keys: &[DoubleDiffKey]) {
        let mut keep_indices = Vec::with_capacity(6 + active_keys.len());
        for i in 0..self.amb_offset() {
            keep_indices.push(i);
        }

        let mut new_ambs = Vec::new();
        for (i, (key, val)) in self.ambiguities.iter().enumerate() {
            if active_keys.contains(key) {
                keep_indices.push(self.amb_offset() + i);
                new_ambs.push((*key, *val));
            }
        }

        // Iono states live after ambiguities; they must compact with the
        // same key set or dim() and cov desync (storm-day crash).
        let iono_base = self.iono_offset();
        for (i, (key, _)) in self.ionos.iter().enumerate() {
            if active_keys.contains(key) {
                keep_indices.push(iono_base + i);
            }
        }
        self.ionos.retain(|(k, _)| active_keys.contains(k));

        let so_base = self.sat_iono_offset();
        for (i, _) in self.sat_ionos.iter().enumerate() {
            keep_indices.push(so_base + i);
        }

        if keep_indices.len() == self.cov.nrows() {
            return;
        }

        let new_dim = keep_indices.len();
        let mut new_cov = DMatrix::zeros(new_dim, new_dim);
        for (new_r, &old_r) in keep_indices.iter().enumerate() {
            for (new_c, &old_c) in keep_indices.iter().enumerate() {
                new_cov[(new_r, new_c)] = self.cov[(old_r, old_c)];
            }
        }
        self.ambiguities = new_ambs;
        self.cov = new_cov;
    }

    /// Extract 3x3 position covariance matrix.
    pub fn extract_pos_cov(&self) -> Matrix3<f64> {
        let mut p = Matrix3::zeros();
        for r in 0..3 {
            for c in 0..3 {
                p[(r, c)] = self.cov[(r, c)];
            }
        }
        p
    }

    /// Extract float ambiguity vector and its covariance submatrix for LAMBDA.
    pub fn extract_amb_block(&self) -> (DVector<f64>, DMatrix<f64>) {
        let off = self.amb_offset();
        let n_amb = self.ambiguities.len();
        let mut a = DVector::zeros(n_amb);
        let mut q = DMatrix::zeros(n_amb, n_amb);
        for i in 0..n_amb {
            a[i] = self.ambiguities[i].1;
            for j in 0..n_amb {
                q[(i, j)] = self.cov[(off + i, off + j)];
            }
        }
        (a, q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gradient_layout_and_roundtrip() {
        let mut state = RtkState::new(Vector3::zeros(), GpsTime::new(2000, 0.0));
        assert_eq!(state.dim(), 6);
        state.enable_zwd(0.0225);
        assert_eq!(state.dim(), 7);
        assert_eq!(state.amb_offset(), 7);
        state.enable_gradients(4e-6);
        // Layout: pos(3) vel(3) zwd@6 gN@7 gE@8 ambs@9..
        assert_eq!(state.dim(), 9);
        assert_eq!(state.amb_offset(), 9);
        assert_eq!(state.grad_idx(), Some((7, 8)));

        let k = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(k, 5.5, 100.0);
        assert_eq!(state.get_amb_idx(&k), Some(9));

        state.grad_n_m = 0.002;
        state.grad_e_m = -0.001;
        let v = state.to_dvector();
        assert!((v[6] - state.zwd_m).abs() < 1e-15);
        assert!((v[7] - 0.002).abs() < 1e-15 && (v[8] + 0.001).abs() < 1e-15);
        assert!((v[9] - 5.5).abs() < 1e-15);

        state.grad_n_m = 0.0;
        state.update_from_dvector(&v);
        assert!((state.grad_n_m - 0.002).abs() < 1e-15);
        assert!((state.grad_e_m + 0.001).abs() < 1e-15);
    }

    #[test]
    fn test_gradients_disabled_keeps_legacy_layout() {
        let mut state = RtkState::new(Vector3::zeros(), GpsTime::new(2000, 0.0));
        state.enable_gradients(4e-6); // without ZWD: columns at 6,7
        assert_eq!(state.dim(), 8);
        assert_eq!(state.amb_offset(), 8);
        assert_eq!(state.grad_idx(), Some((6, 7)));
    }

    #[test]
    fn test_rtk_state_lifecycle() {
        let mut state = RtkState::new(Vector3::new(100.0, 200.0, 300.0), GpsTime::new(2000, 100.0));
        assert_eq!(state.dim(), 6);

        let k1 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        let k2 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };

        state.ensure_ambiguity(k1, 5.0, 100.0);
        state.ensure_ambiguity(k2, 10.0, 100.0);
        assert_eq!(state.dim(), 8);
        assert_eq!(state.get_amb_idx(&k1), Some(6));
        assert_eq!(state.get_amb_idx(&k2), Some(7));

        let vec = state.to_dvector();
        assert_eq!(vec[6], 5.0);
        assert_eq!(vec[7], 10.0);

        state.retain_active_ambiguities(&[k2]);
        assert_eq!(state.dim(), 7);
        assert_eq!(state.get_amb_idx(&k2), Some(6));
    }
}

#[cfg(test)]
mod iono_tests {
    use super::*;

    fn dd_key(sv: u16, band: u8) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat: sv, ref_sat: 1, freq_band: band }
    }

    #[test]
    fn test_iono_disabled_zero_impact() {
        let t = GpsTime::new(2100, 0.0);
        let st = RtkState::new(Vector3::zeros(), t);
        assert!(!st.iono_enabled);
        assert_eq!(st.ionos.len(), 0);
        // dim unchanged from legacy
        assert_eq!(st.dim(), 6);
    }

    #[test]
    fn test_iono_dim_and_offset() {
        let t = GpsTime::new(2100, 0.0);
        let mut st = RtkState::new(Vector3::zeros(), t);
        st.iono_enabled = true;
        st.ensure_ambiguity(dd_key(5, 1), 100.0, 100.0);
        st.ensure_iono(dd_key(5, 1), 4.0);
        // pos(3)+vel(3)+1 amb+1 iono = 8
        assert_eq!(st.dim(), 8);
        assert_eq!(st.iono_offset(), st.amb_offset() + 1);
        assert_eq!(st.get_iono_idx(&dd_key(5, 1)), Some(7));
    }

    #[test]
    fn test_iono_roundtrip_preserves_values() {
        let t = GpsTime::new(2100, 0.0);
        let mut st = RtkState::new(Vector3::zeros(), t);
        st.zwd_enabled = true;
        st.iono_enabled = true;
        st.ensure_ambiguity(dd_key(5, 1), 50.0, 100.0);
        st.ensure_iono(dd_key(5, 1), 4.0);
        st.ionos[0].1 = -1.234;
        let v = st.to_dvector();
        let mut st2 = st.clone();
        st2.update_from_dvector(&v);
        assert!((st2.ionos[0].1 - (-1.234)).abs() < 1e-12);
    }

    #[test]
    fn test_ensure_iono_idempotent() {
        let t = GpsTime::new(2100, 0.0);
        let mut st = RtkState::new(Vector3::zeros(), t);
        st.iono_enabled = true;
        st.ensure_iono(dd_key(5, 1), 4.0);
        st.ensure_iono(dd_key(5, 1), 4.0); // second call no-op
        assert_eq!(st.ionos.len(), 1);
        assert_eq!(st.dim(), 7);
    }
}

#[cfg(test)]
mod iono_retain_tests {
    use super::*;

    fn key(sat: u16) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: 1 }
    }

    /// Ledger row 12 candidate: retain_active_ambiguities must compact
    /// iono states alongside ambiguities. Storm-day crash reproduced here:
    /// dim() counted stale ionos while cov dropped them -> gemm mismatch.
    #[test]
    fn retain_compacts_iono_states_and_keeps_cov_consistent() {
        let t = GpsTime::new(2100, 0.0);
        let mut st = RtkState::new(Vector3::zeros(), t);
        st.iono_enabled = true;
        st.ensure_ambiguity(key(5), 10.0, 100.0);
        st.ensure_ambiguity(key(7), 20.0, 100.0);
        st.ensure_iono(key(5), 4.0);
        st.ensure_iono(key(7), 4.0);
        assert_eq!(st.dim(), st.cov.nrows(), "pre-condition");

        // Satellite 7 sets; only key(5) remains active.
        st.retain_active_ambiguities(&[key(5)]);

        assert_eq!(
            st.dim(),
            st.cov.nrows(),
            "dim/cov desync after retain — the storm-day crash"
        );
        assert_eq!(st.ambiguities.len(), 1);
        assert_eq!(st.ionos.len(), 1, "ionos must compact with ambs");
        assert_eq!(st.get_iono_idx(&key(5)), Some(st.iono_offset()));
    }
}

#[cfg(test)]
mod sat_iono_tests {
    use super::*;

    /// Per-satellite mapped-iono contract: states keyed by SATELLITE, not
    /// DD pair; a pair's H row uses I_sat − I_ref so geometry couples all
    /// satellites through the shared reference — rank deficiency broken.
    #[test]
    fn sat_iono_states_keyed_by_satellite_and_shared_across_pairs() {
        let t = GpsTime::new(2100, 0.0);
        let mut st = RtkState::new(Vector3::zeros(), t);
        st.sat_iono_enabled = true;
        st.ensure_sat_iono(5);
        st.ensure_sat_iono(7);
        st.ensure_sat_iono(9); // reference candidate
        assert_eq!(st.dim(), 6 + 3);
        // Pairs (5,9) and (7,9) both draw from sat 9's single state.
        assert_eq!(st.get_sat_iono_idx(9), Some(8));
        assert_ne!(st.get_sat_iono_idx(5), st.get_sat_iono_idx(7));
    }

    #[test]
    fn sat_iono_multi_constellation_keys() {
        let t = GpsTime::new(2100, 0.0);
        let mut st = RtkState::new(Vector3::zeros(), t);
        st.sat_iono_enabled = true;
        st.ensure_sat_iono_key(0, 10); // GPS PRN 10
        st.ensure_sat_iono_key(2, 10); // Galileo PRN 10
        assert_eq!(st.dim(), 6 + 2);
        let idx_gps = st.get_sat_iono_key_idx(0, 10);
        let idx_gal = st.get_sat_iono_key_idx(2, 10);
        assert!(idx_gps.is_some());
        assert!(idx_gal.is_some());
        assert_ne!(idx_gps, idx_gal);
    }
}
