//! Partial Ambiguity Resolution (PAR) candidate selection and subset partitioning.

use nalgebra::{DMatrix, DVector};
use super::state::{DoubleDiffKey, RtkState};

pub const MAX_SUBSET_SIZE: usize = 16;

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

/// Select and sort candidate ambiguity indices for PAR.
pub fn select_par_candidates(
    state: &RtkState,
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    min_ambs: usize,
    is_kinematic: bool,
) -> (Vec<usize>, usize) {
    let n_amb = a_float.len();
    if !is_kinematic {
        let mut idx: Vec<usize> = (0..n_amb).collect();
        idx.sort_by(|&i, &j| q_amb[(i, i)].total_cmp(&q_amb[(j, j)]));
        let m = (n_amb - 1).min(MAX_SUBSET_SIZE);
        return (idx, m);
    }

    let min_k = min_ambs.max(6);
    let mut c: Vec<usize> = (0..n_amb)
        .filter(|&i| {
            let key = &state.ambiguities[i].0;
            q_amb[(i, i)] <= 1.0 && !is_beidou_geo_key(key)
        })
        .collect();

    if c.len() < min_k {
        c = (0..n_amb)
            .filter(|&i| !is_beidou_geo_key(&state.ambiguities[i].0))
            .collect();
    }
    if c.len() < min_k {
        c = (0..n_amb).collect();
    }

    let score = |i: usize| {
        let frac = (a_float[i] - a_float[i].round()).abs();
        frac + q_amb[(i, i)].sqrt() * 0.5
    };
    c.sort_by(|&i, &j| score(i).total_cmp(&score(j)));
    let m = c.len().min(MAX_SUBSET_SIZE);
    (c, m)
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
}
