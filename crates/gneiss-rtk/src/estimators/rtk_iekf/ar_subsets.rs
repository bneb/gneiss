//! Partial Ambiguity Resolution (PAR) candidate selection and subset partitioning.

use nalgebra::{DMatrix, DVector, Vector3};
pub use crate::ambiguity::par::{compute_cqm, AmbiguityMetadata};
use super::state::{DoubleDiffKey, RtkState};
use super::DoubleDiffMeasurement;

pub const MAX_SUBSET_SIZE: usize = 16;
pub const MIN_PAR_SUBSET_SIZE: usize = 4;
pub const MAX_ACCEPTABLE_PDOP: f64 = 15.0;

/// Stack-allocated ambiguity index subset (zero heap allocation on critical path).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AmbiguitySubset {
    pub indices: [usize; MAX_SUBSET_SIZE],
    pub len: usize,
}

impl AmbiguitySubset {
    pub fn empty() -> Self {
        Self {
            indices: [0; MAX_SUBSET_SIZE],
            len: 0,
        }
    }

    pub fn from_slice(slice: &[usize]) -> Option<Self> {
        if slice.len() > MAX_SUBSET_SIZE {
            return None;
        }
        let mut sub = Self::empty();
        sub.indices[..slice.len()].copy_from_slice(slice);
        sub.len = slice.len();
        Some(sub)
    }

    pub fn as_slice(&self) -> &[usize] {
        &self.indices[..self.len]
    }
}

/// Check if a double-difference key involves a stationary BeiDou GEO satellite.
pub fn is_beidou_geo_key(key: &DoubleDiffKey) -> bool {
    if key.constellation_id != 3 {
        return false;
    }
    let is_geo = |prn: u16| (1..=5).contains(&prn) || (59..=63).contains(&prn);
    is_geo(key.sat) || is_geo(key.ref_sat)
}

/// Check if a subset of satellite positions satisfies the DOP geometry guard.
pub fn validate_subset_geometry(
    rover_pos: Vector3<f64>,
    sat_positions: &[Vector3<f64>],
    max_pdop: f64,
) -> bool {
    if sat_positions.len() < MIN_PAR_SUBSET_SIZE {
        return false;
    }
    match gneiss_core::dop::compute_dop_from_positions(rover_pos, sat_positions) {
        Some(dop) => dop.pdop > 0.0 && dop.pdop <= max_pdop && dop.hdop <= 10.0,
        None => false,
    }
}

fn push_unique_pos(list: &mut Vec<Vector3<f64>>, pos: Vector3<f64>) {
    if !list.iter().any(|p| (p - pos).norm() < 1.0) {
        list.push(pos);
    }
}

/// Validate geometry of double-difference subset against minimum size and DOP thresholds.
pub fn validate_dd_subset_geometry(
    state: &RtkState,
    dd: &[DoubleDiffMeasurement],
    subset_idx: &[usize],
) -> bool {
    if subset_idx.len() < 3 {
        return false;
    }
    let mut sat_positions = Vec::with_capacity(subset_idx.len() + 1);
    for &idx in subset_idx {
        if idx >= state.ambiguities.len() { continue; }
        let key = state.ambiguities[idx].0;
        if let Some(m) = dd.iter().find(|m| m.key == key) {
            push_unique_pos(&mut sat_positions, m.sat_pos);
            push_unique_pos(&mut sat_positions, m.ref_pos);
        }
    }
    validate_subset_geometry(state.pos_ecef, &sat_positions, MAX_ACCEPTABLE_PDOP)
}

/// Construct ambiguity quality metadata from state and double-difference measurements.
pub fn compute_metadata_from_dd(
    state: &RtkState,
    dd: &[DoubleDiffMeasurement],
) -> Vec<AmbiguityMetadata> {
    let mut meta = Vec::with_capacity(state.ambiguities.len());
    let rover_llh = gneiss_core::coords::ecef_to_llh(state.pos_ecef);
    for (key, _) in &state.ambiguities {
        if let Some(m) = dd.iter().find(|m| m.key == *key) {
            let (_, el) = gneiss_core::coords::az_el(rover_llh, state.pos_ecef, m.sat_pos);
            let cmc_sigma = m.pr_var_m2.max(1e-4).sqrt().min(10.0);
            let snr_est = (45.0 - 10.0 * (m.cp_var_cycles2 / 0.0005).max(1e-4).log10()).clamp(20.0, 50.0);
            meta.push(AmbiguityMetadata::new(el, snr_est, 30.0, cmc_sigma));
        } else {
            meta.push(AmbiguityMetadata::default());
        }
    }
    meta
}

fn is_severe_nlos(meta: Option<&AmbiguityMetadata>) -> bool {
    let Some(m) = meta else { return false };
    m.elevation_rad < 25.0_f64.to_radians() && (m.snr_dbhz < 28.0 || m.cmc_sigma > 3.0)
}

/// Select and sort candidate ambiguity indices for PAR using CQM metadata.
fn filter_kinematic_pool(
    state: &RtkState,
    q_amb: &DMatrix<f64>,
    n_amb: usize,
    min_k: usize,
    meta: Option<&[AmbiguityMetadata]>,
) -> Vec<usize> {
    let mut c: Vec<usize> = (0..n_amb)
        .filter(|&i| {
            q_amb[(i, i)] <= 1.0
                && !is_beidou_geo_key(&state.ambiguities[i].0)
                && !is_severe_nlos(meta.and_then(|m| m.get(i)))
        })
        .collect();
    if c.len() < min_k {
        c = (0..n_amb)
            .filter(|&i| q_amb[(i, i)] <= 1.0 && !is_beidou_geo_key(&state.ambiguities[i].0))
            .collect();
    }
    if c.len() < min_k {
        c = (0..n_amb).filter(|&i| !is_beidou_geo_key(&state.ambiguities[i].0)).collect();
    }
    if c.len() < min_k {
        c = (0..n_amb).collect();
    }
    c
}

/// Select and sort candidate ambiguity indices for PAR using CQM metadata.
pub fn select_par_candidates_with_metadata(
    state: &RtkState,
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    min_ambs: usize,
    is_kinematic: bool,
    metadata: Option<&[AmbiguityMetadata]>,
) -> (Vec<usize>, usize) {
    let n_amb = a_float.len();
    if !is_kinematic {
        let mut idx: Vec<usize> = (0..n_amb).collect();
        idx.sort_by(|&i, &j| q_amb[(i, i)].total_cmp(&q_amb[(j, j)]));
        let m = (n_amb - 1).min(MAX_SUBSET_SIZE);
        return (idx, m);
    }

    let min_k = min_ambs.max(4);
    let mut c = filter_kinematic_pool(state, q_amb, n_amb, min_k, metadata);
    let score = |i: usize| {
        let m = metadata.and_then(|meta| meta.get(i));
        compute_cqm(m, q_amb[(i, i)], a_float[i])
    };
    c.sort_by(|&i, &j| score(j).total_cmp(&score(i)));
    let m = c.len().min(MAX_SUBSET_SIZE);
    (c, m)
}

/// Select and sort candidate ambiguity indices for PAR with optional DD measurements.
pub fn select_par_candidates_with_dd(
    state: &RtkState,
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    min_ambs: usize,
    is_kinematic: bool,
    dd_meas: Option<&[DoubleDiffMeasurement]>,
) -> (Vec<usize>, usize) {
    let meta = dd_meas.map(|dd| compute_metadata_from_dd(state, dd));
    select_par_candidates_with_metadata(state, a_float, q_amb, min_ambs, is_kinematic, meta.as_deref())
}

/// Select and sort candidate ambiguity indices for PAR.
pub fn select_par_candidates(
    state: &RtkState,
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    min_ambs: usize,
    is_kinematic: bool,
) -> (Vec<usize>, usize) {
    select_par_candidates_with_metadata(state, a_float, q_amb, min_ambs, is_kinematic, None)
}

fn build_cluster_subset(
    state: &RtkState,
    ranked_candidates: &[usize],
    cluster: &[u8],
) -> (usize, [usize; MAX_SUBSET_SIZE]) {
    let mut subset_indices = [0usize; MAX_SUBSET_SIZE];
    let mut count = 0;
    for &idx in ranked_candidates {
        if count >= MAX_SUBSET_SIZE {
            break;
        }
        let const_id = state.ambiguities[idx].0.constellation_id;
        if cluster.contains(&const_id) {
            subset_indices[count] = idx;
            count += 1;
        }
    }
    (count, subset_indices)
}

/// Generate constellation-partitioned subsets (e.g., GPS+Galileo, GPS+BeiDou).
pub fn partition_constellation_subsets(
    state: &RtkState,
    ranked_candidates: &[usize],
    min_k: usize,
) -> Vec<AmbiguitySubset> {
    let clusters: [&[u8]; 3] = [
        &[0, 2], // GPS (0) + Galileo (2) - isolates BeiDou multipath
        &[0, 3], // GPS (0) + BeiDou (3) - isolates Galileo multipath
        &[2, 3], // Galileo (2) + BeiDou (3) - isolates GPS canyon reflection
    ];

    let mut result = Vec::with_capacity(clusters.len());
    for &cluster in &clusters {
        let (count, subset_indices) = build_cluster_subset(state, ranked_candidates, cluster);
        if count >= min_k {
            result.push(AmbiguitySubset {
                indices: subset_indices,
                len: count,
            });
        }
    }
    result
}

/// Generate 1-omission subsets from the top candidate pool.
pub fn generate_omission_subsets(
    sorted_indices: &[usize],
    pool_len: usize,
) -> Vec<AmbiguitySubset> {
    let actual_pool = pool_len.min(sorted_indices.len()).min(MAX_SUBSET_SIZE);
    let mut subsets = Vec::with_capacity(actual_pool);
    for omit in 0..actual_pool {
        let mut indices = [0usize; MAX_SUBSET_SIZE];
        let mut count = 0;
        for (i, &idx) in sorted_indices[..actual_pool].iter().enumerate() {
            if i != omit {
                indices[count] = idx;
                count += 1;
            }
        }
        subsets.push(AmbiguitySubset {
            indices,
            len: count,
        });
    }
    subsets
}

/// Generate 2-omission subsets from the top candidate pool when 1-omission fails.
pub fn generate_two_omission_subsets(
    sorted_indices: &[usize],
    pool_len: usize,
    min_k: usize,
) -> Vec<AmbiguitySubset> {
    let actual_pool = pool_len.min(sorted_indices.len()).min(10);
    if actual_pool < min_k + 2 {
        return Vec::new();
    }
    let mut subsets = Vec::new();
    for omit1 in 0..actual_pool {
        for omit2 in (omit1 + 1)..actual_pool {
            let mut indices = [0usize; MAX_SUBSET_SIZE];
            let mut count = 0;
            for (i, &idx) in sorted_indices[..actual_pool].iter().enumerate() {
                if i != omit1 && i != omit2 {
                    indices[count] = idx;
                    count += 1;
                }
            }
            if count >= min_k {
                subsets.push(AmbiguitySubset {
                    indices,
                    len: count,
                });
            }
        }
    }
    subsets
}

/// Extract float ambiguity subvector and covariance submatrix for given indices.
pub fn extract_subset(
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    indices: &[usize],
) -> (DVector<f64>, DMatrix<f64>) {
    let k = indices.len();
    let mut sub_a = DVector::zeros(k);
    let mut sub_q = DMatrix::zeros(k, k);

    for (r, &i) in indices.iter().enumerate() {
        sub_a[r] = a_float[i];
        for (c, &j) in indices.iter().enumerate() {
            sub_q[(r, c)] = q_amb[(i, j)];
        }
    }
    (sub_a, sub_q)
}

#[cfg(test)]
mod subset_tests;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subset_empty_and_from_slice() {
        let empty = AmbiguitySubset::empty();
        assert_eq!(empty.len, 0);
        let expected_empty: &[usize] = &[];
        assert_eq!(empty.as_slice(), expected_empty);

        let slice = [1, 3, 5, 7, 9, 11];
        let sub = AmbiguitySubset::from_slice(&slice).expect("valid slice");
        assert_eq!(sub.len, 6);
        assert_eq!(sub.as_slice(), &slice);

        let over = [0; 17];
        assert!(AmbiguitySubset::from_slice(&over).is_none());
    }

    #[test]
    fn test_beidou_geo_key_detection() {
        let geo_key = DoubleDiffKey {
            constellation_id: 3,
            sat: 1, // GEO PRN 1
            ref_sat: 6, // IGSO PRN 6
            freq_band: 1,
        };
        assert!(is_beidou_geo_key(&geo_key));

        let igso_key = DoubleDiffKey {
            constellation_id: 3,
            sat: 7, // IGSO PRN 7
            ref_sat: 6, // IGSO PRN 6
            freq_band: 1,
        };
        assert!(!is_beidou_geo_key(&igso_key));

        let gps_key = DoubleDiffKey {
            constellation_id: 0,
            sat: 1,
            ref_sat: 2,
            freq_band: 1,
        };
        assert!(!is_beidou_geo_key(&gps_key));
    }

    #[test]
    fn test_generate_omission_subsets() {
        let indices = [10, 20, 30, 40, 50, 60, 70, 80];
        let subsets = generate_omission_subsets(&indices, 8);
        assert_eq!(subsets.len(), 8);
        for (omit, sub) in subsets.iter().enumerate() {
            assert_eq!(sub.len, 7);
            assert!(!sub.as_slice().contains(&indices[omit]));
        }
    }

    #[test]
    fn test_generate_two_omission_subsets() {
        let indices = [1, 2, 3, 4, 5, 6, 7, 8];
        let subsets = generate_two_omission_subsets(&indices, 8, 6);
        assert_eq!(subsets.len(), 28); // 8 choose 2 = 28
        for sub in &subsets {
            assert_eq!(sub.len, 6);
        }
    }

    #[test]
    fn test_extract_subset_exactness_and_spd() {
        let a = DVector::from_vec(vec![1.1, 2.2, 3.3, 4.4]);
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.1, 0.2, 0.3, 0.4]));
        let idx = [0, 2, 3];
        let (sub_a, sub_q) = extract_subset(&a, &q, &idx);
        assert_eq!(sub_a.len(), 3);
        assert_eq!(sub_a[0], 1.1);
        assert_eq!(sub_a[1], 3.3);
        assert_eq!(sub_a[2], 4.4);
        assert_eq!(sub_q[(0, 0)], 0.1);
        assert_eq!(sub_q[(1, 1)], 0.3);
        assert_eq!(sub_q[(2, 2)], 0.4);
        assert_eq!(sub_q[(0, 1)], 0.0);
        assert!(sub_q.cholesky().is_some());
    }

    #[test]
    fn test_cqm_candidate_prioritization_over_coincidental_float() {
        use gneiss_core::time::GpsTime;
        let mut state = RtkState::new(Vector3::new(100.0, 200.0, 300.0), GpsTime::new(2000, 100.0));
        let k1 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        let k2 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(k1, 10.25, 0.04);
        state.ensure_ambiguity(k2, 20.001, 0.04);

        let a = DVector::from_vec(vec![10.25, 20.001]);
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.04, 0.04]));
        let meta = vec![
            AmbiguityMetadata::new(1.10, 45.0, 40.0, 0.15), // Index 0: high-el clean
            AmbiguityMetadata::new(0.18, 22.0, 4.0, 3.5),   // Index 1: low-el corrupted
        ];

        let (ranked, _) = select_par_candidates_with_metadata(&state, &a, &q, 2, true, Some(&meta));
        assert_eq!(ranked[0], 0, "High elevation clean satellite must be ranked first");
        assert_eq!(ranked[1], 1, "Multipath corrupted satellite must be ranked second");
    }

    #[test]
    fn test_validate_subset_geometry_rejects_sub_four_and_collinear() {
        let rover = Vector3::new(6_378_137.0, 0.0, 0.0);
        let sats3 = vec![
            Vector3::new(20_000_000.0, 5_000_000.0, 5_000_000.0),
            Vector3::new(22_000_000.0, -5_000_000.0, 5_000_000.0),
            Vector3::new(19_000_000.0, 5_000_000.0, -5_000_000.0),
        ];
        assert!(!validate_subset_geometry(rover, &sats3, 10.0), "Subset with < 4 satellites must be rejected");

        let collinear = vec![
            Vector3::new(20_000_000.0, 0.0, 5_000_000.0),
            Vector3::new(22_000_000.0, 0.0, 10_000_000.0),
            Vector3::new(19_000_000.0, 0.0, -5_000_000.0),
            Vector3::new(21_000_000.0, 0.0, -10_000_000.0),
        ];
        assert!(!validate_subset_geometry(rover, &collinear, 10.0), "Collinear satellites must be rejected");

        let spread = vec![
            Vector3::new(26_000_000.0, 0.0, 0.0),
            Vector3::new(20_000_000.0, 0.0, 15_000_000.0),
            Vector3::new(20_000_000.0, 13_000_000.0, -10_000_000.0),
            Vector3::new(20_000_000.0, -13_000_000.0, -10_000_000.0),
        ];
        assert!(validate_subset_geometry(rover, &spread, 10.0), "Well-spread satellites must pass geometry guard");
    }
}
