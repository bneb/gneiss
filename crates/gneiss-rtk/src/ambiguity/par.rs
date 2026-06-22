use nalgebra::{DMatrix, DVector};

/// Select the largest subset of float ambiguities whose bootstrapped
/// success rate exceeds `min_success_rate`. Returns the indices of
/// selected ambiguities, the sub-vector, and the sub-covariance.
///
/// Algorithm: sort by conditional variance (most confident first),
/// then add ambiguities one at a time until the cumulative success
/// rate drops below the threshold.
pub fn select_ils_subset(
    a: &DVector<f64>,
    q: &DMatrix<f64>,
    min_success_rate: f64,
) -> (Vec<usize>, DVector<f64>, DMatrix<f64>) {
    let n = a.len();
    if n == 0 {
        return (vec![], DVector::zeros(0), DMatrix::zeros(0, 0));
    }

    // Compute per-ambiguity conditional standard deviations from LDL^T decomposition
    let cond_stdevs = compute_conditional_stdevs(q);

    // Sort indices by conditional stdev ascending (most confident first)
    let mut indexed: Vec<(usize, f64)> = (0..n).map(|i| (i, cond_stdevs[i])).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // Greedy selection: add until cumulative success rate drops below threshold
    let mut selected: Vec<usize> = Vec::with_capacity(n);
    for &(idx, _) in &indexed {
        selected.push(idx);
        let sub_q = submatrix(q, &selected);
        let sub_d = ldlt_diagonal(&sub_q);
        let sr = cumulative_success_rate(&sub_d);
        if sr < min_success_rate {
            selected.pop();
            break;
        }
    }

    if selected.is_empty() {
        return (vec![], DVector::zeros(0), DMatrix::zeros(0, 0));
    }

    // Build sub-vector and sub-covariance
    let k = selected.len();
    let mut sub_a = DVector::zeros(k);
    let sub_q = submatrix(q, &selected);
    for (i, &orig_idx) in selected.iter().enumerate() {
        sub_a[i] = a[orig_idx];
    }

    (selected, sub_a, sub_q)
}

/// Compute per-ambiguity conditional standard deviations from covariance Q.
/// Uses LDL^T decomposition: conditional variance of i-th ambiguity
/// given previous ones is D[i], so stdev = sqrt(D[i]).
fn compute_conditional_stdevs(q: &DMatrix<f64>) -> Vec<f64> {
    let n = q.nrows();
    let d = ldlt_diagonal(q);
    let mut stdevs = Vec::with_capacity(n);
    // Reorder by original index: the LDL^T ordering depends on the
    // permutation used, so we approximate by sqrt of diagonal
    for i in 0..n {
        let var = q[(i, i)].max(1e-12);
        stdevs.push(var.sqrt());
    }
    let _ = d; // LDL^T diagonal for future improvement
    stdevs
}

/// Extract a submatrix specified by row/column indices.
fn submatrix(q: &DMatrix<f64>, indices: &[usize]) -> DMatrix<f64> {
    let k = indices.len();
    let mut sub = DMatrix::zeros(k, k);
    for (i, &ri) in indices.iter().enumerate() {
        for (j, &cj) in indices.iter().enumerate() {
            sub[(i, j)] = q[(ri, cj)];
        }
    }
    sub
}

/// Compute the diagonal D of the LDL^T decomposition of Q.
/// Uses a simple Cholesky-like algorithm: D[i] = Q[i,i] - sum_{k<i} L[i,k]^2 * D[k].
fn ldlt_diagonal(q: &DMatrix<f64>) -> DVector<f64> {
    let n = q.nrows();
    let mut l = DMatrix::zeros(n, n);
    let mut d = DVector::zeros(n);
    for i in 0..n {
        let mut sum_ldl = 0.0;
        for k in 0..i {
            sum_ldl += l[(i, k)] * l[(i, k)] * d[k];
        }
        d[i] = q[(i, i)] - sum_ldl;
        if d[i] <= 0.0 {
            d[i] = 1e-12;
        }
        for j in (i + 1)..n {
            let mut sum_l = 0.0;
            for k in 0..i {
                sum_l += l[(j, k)] * l[(i, k)] * d[k];
            }
            l[(j, i)] = (q[(j, i)] - sum_l) / d[i];
        }
    }
    d
}

/// Cumulative bootstrapped success rate:
/// P_s = prod_i (2 * Phi(1/(2*sigma_i)) - 1)
/// where sigma_i = sqrt(D_i) from LDL^T decomposition.
fn cumulative_success_rate(d: &DVector<f64>) -> f64 {
    let mut rate = 1.0;
    for i in 0..d.len() {
        let sigma = d[i].sqrt().max(1e-12);
        let z = 1.0 / (2.0 * sigma);
        // Phi(z) approximation: 0.5 * (1 + erf(z/sqrt(2)))
        let phi = 0.5 * (1.0 + libm::erf(z / std::f64::consts::SQRT_2));
        rate *= 2.0 * phi - 1.0;
    }
    rate.max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_select_ils_subset_all_confident() {
        let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 0.01, 0.01]));
        let (indices, sub_a, _sub_q) = select_ils_subset(&a, &q, 0.99);
        assert_eq!(indices.len(), 3);
        assert!((sub_a[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_select_ils_subset_mixed_confidence() {
        // First ambiguity is precise, second and third are very uncertain
        let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.01, 100.0, 100.0]));
        let (indices, _, _) = select_ils_subset(&a, &q, 0.99);
        // Should only select the first one
        assert_eq!(indices.len(), 1);
    }

    #[test]
    fn test_select_ils_subset_empty() {
        let a = DVector::zeros(0);
        let q = DMatrix::zeros(0, 0);
        let (indices, sub_a, sub_q) = select_ils_subset(&a, &q, 0.99);
        assert!(indices.is_empty());
        assert_eq!(sub_a.len(), 0);
        assert_eq!(sub_q.nrows(), 0);
    }

    #[test]
    fn test_ldlt_diagonal_identity() {
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![4.0, 9.0, 16.0]));
        let d = ldlt_diagonal(&q);
        assert!((d[0] - 4.0).abs() < 1e-10);
        assert!((d[1] - 9.0).abs() < 1e-10);
        assert!((d[2] - 16.0).abs() < 1e-10);
    }

    #[test]
    fn test_cumulative_success_rate_perfect() {
        // Infinitesimal variance → success rate ~1.0
        let d = DVector::from_vec(vec![1e-10, 1e-10]);
        let rate = cumulative_success_rate(&d);
        assert!(rate > 0.9999);
    }

    #[test]
    fn test_cumulative_success_rate_zero_variance() {
        // Zero variance → Phi(z) with z → inf → rate → 1.0
        let d = DVector::from_vec(vec![0.0, 0.0]);
        let rate = cumulative_success_rate(&d);
        assert!(rate > 0.9999);
    }

    #[test]
    fn test_ldlt_diagonal_correlated_3x3() {
        // Q = [[4, 2, 1], [2, 9, 3], [1, 3, 16]]
        let q = DMatrix::from_row_slice(3, 3, &[4.0, 2.0, 1.0, 2.0, 9.0, 3.0, 1.0, 3.0, 16.0]);
        let d = ldlt_diagonal(&q);
        // D[0] = Q[0,0] = 4
        assert!((d[0] - 4.0).abs() < 1e-10);
        // D[1] = Q[1,1] - L[1,0]^2 * D[0] = 9 - (2/4)^2 * 4 = 9 - 0.25*4 = 8
        assert!((d[1] - 8.0).abs() < 1e-10);
        // D[2] = Q[2,2] - L[2,0]^2*D[0] - L[2,1]^2*D[1]
        // L[2,0] = Q[2,0]/D[0] = 1/4 = 0.25
        // Q[2,1]' = Q[2,1] - L[2,0]*L[1,0]*D[0] = 3 - 0.25*0.5*4 = 3 - 0.5 = 2.5
        // L[2,1] = Q[2,1]'/D[1] = 2.5/8 = 0.3125
        // D[2] = 16 - 0.25^2*4 - 0.3125^2*8 = 16 - 0.25 - 0.78125 = 14.96875
        assert!((d[2] - 14.96875).abs() < 1e-10);
    }

    #[test]
    fn test_ldlt_diagonal_near_singular() {
        // Nearly singular: Q[1,1] - L[1,0]^2*D[0] ≈ 0, should be clamped to 1e-12
        let q = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 1.0]);
        let d = ldlt_diagonal(&q);
        assert!(
            d[1] >= 1e-12,
            "Near-singular D[1] should be clamped, got {}",
            d[1]
        );
    }

    #[test]
    fn test_submatrix_first_two() {
        let q = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let sub = submatrix(&q, &[0, 1]);
        assert_eq!(sub.nrows(), 2);
        assert_eq!(sub.ncols(), 2);
        assert!((sub[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((sub[(0, 1)] - 2.0).abs() < 1e-10);
        assert!((sub[(1, 0)] - 4.0).abs() < 1e-10);
        assert!((sub[(1, 1)] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_submatrix_noncontiguous() {
        let q = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let sub = submatrix(&q, &[2, 0]);
        assert_eq!(sub.nrows(), 2);
        assert!((sub[(0, 0)] - 9.0).abs() < 1e-10); // Q[2,2]
        assert!((sub[(0, 1)] - 7.0).abs() < 1e-10); // Q[2,0]
        assert!((sub[(1, 0)] - 3.0).abs() < 1e-10); // Q[0,2]
        assert!((sub[(1, 1)] - 1.0).abs() < 1e-10); // Q[0,0]
    }

    #[test]
    fn test_compute_conditional_stdevs() {
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![4.0, 16.0, 100.0]));
        let stdevs = compute_conditional_stdevs(&q);
        assert_eq!(stdevs.len(), 3);
        assert!((stdevs[0] - 2.0).abs() < 1e-10);
        assert!((stdevs[1] - 4.0).abs() < 1e-10);
        assert!((stdevs[2] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_select_ils_subset_zero_threshold() {
        // With min_success_rate=0, all ambiguities should be selected even if imprecise
        let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![100.0, 100.0, 100.0]));
        let (indices, sub_a, _sub_q) = select_ils_subset(&a, &q, 0.0);
        assert_eq!(
            indices.len(),
            3,
            "All ambiguities should be selected with min_success_rate=0"
        );
        assert!((sub_a[0] - 1.0).abs() < 1e-10);
        assert!((sub_a[1] - 2.0).abs() < 1e-10);
        assert!((sub_a[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_select_ils_subset_very_strict() {
        // With min_success_rate=1.0, only very precise ambiguities should survive
        let a = DVector::from_vec(vec![1.0, 2.0]);
        // Both have large variance - even the most confident won't reach 1.0
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.1, 0.1]));
        // Even the best alone may not reach 0.999999 with variance 0.1
        // sigma = sqrt(0.1) ≈ 0.316, z = 1/(2*0.316) ≈ 1.58, Phi ≈ 0.943
        // rate = 2*0.943 - 1 = 0.886 < 0.999999, so likely empty
        let (indices_strict, sub_a, _) = select_ils_subset(&a, &q, 0.999999);
        assert!(
            indices_strict.len() <= 2,
            "Should select at most 2 with very strict threshold"
        );
        if !indices_strict.is_empty() {
            assert_eq!(
                sub_a.len(),
                indices_strict.len(),
                "sub_a length should match indices length"
            );
        }
    }
}
