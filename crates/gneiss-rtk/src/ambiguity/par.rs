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
}
