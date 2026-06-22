use nalgebra::{DMatrix, DVector};

#[derive(Debug)]
pub struct LdltResult {
    pub l: DMatrix<f64>,
    pub d: DVector<f64>,
}

pub struct DecorrelateResult {
    pub z_hat: DVector<f64>,
    pub q_z: DMatrix<f64>,
    pub z_mat: DMatrix<f64>,
    pub l: DMatrix<f64>,
    pub d: DVector<f64>,
}

/// Output of the LAMBDA algorithm
#[derive(Debug, Clone, PartialEq)]
pub struct LambdaResult {
    /// The best integer ambiguity combination
    pub best_integers: DVector<f64>,
    /// The second best integer ambiguity combination (used for ratio test)
    pub second_best_integers: DVector<f64>,
    /// The ratio test value (sq_norm_second_best / sq_norm_best)
    pub ratio: f64,
    /// The bootstrapping success rate of the float solution
    pub success_rate: f64,
}

/// Resolves integer ambiguities using the LAMBDA method (Decorrelation + Search).
pub fn resolve_lambda(a: &DVector<f64>, q: &DMatrix<f64>) -> Result<LambdaResult, &'static str> {
    resolve_lambda_inner(a, q, 10000)
}

fn resolve_lambda_inner(
    a: &DVector<f64>,
    q: &DMatrix<f64>,
    max_iters: usize,
) -> Result<LambdaResult, &'static str> {
    let n = a.len();
    if n == 0 {
        return Err("Empty ambiguity vector");
    }

    let dec = decorrelate(a, q)?;
    let (best_z, second_best_z) = run_lambda_search(n, &dec, max_iters)?;

    let t_inv = dec
        .z_mat
        .transpose()
        .try_inverse()
        .ok_or("Transformation matrix inversion failed")?;
    let best_a = &t_inv * &best_z.0;
    let second_best_a = &t_inv * &second_best_z.0;

    let success_rate = bootstrapping_success_rate(&dec.d);
    let safe_best_dist = if best_z.1 < 1e-12 { 1e-12 } else { best_z.1 };

    Ok(LambdaResult {
        best_integers: best_a,
        second_best_integers: second_best_a,
        ratio: second_best_z.1 / safe_best_dist,
        success_rate,
    })
}

fn run_lambda_search(
    n: usize,
    dec: &DecorrelateResult,
    max_iters: usize,
) -> Result<((DVector<f64>, f64), (DVector<f64>, f64)), &'static str> {
    let mut best_z = DVector::zeros(n);
    let mut best_dist = f64::MAX;
    let mut second_best_z = DVector::zeros(n);
    let mut second_best_dist = f64::MAX;
    let mut current_z = DVector::zeros(n);
    let mut iter_count = 0;
    let mut y = DVector::zeros(n);

    search_recursive(
        (n - 1) as isize,
        n,
        &dec.l,
        &dec.d,
        &dec.z_hat,
        &mut y,
        &mut current_z,
        0.0,
        &mut best_z,
        &mut best_dist,
        &mut second_best_z,
        &mut second_best_dist,
        &mut iter_count,
        max_iters,
    );

    if iter_count > max_iters || best_dist >= f64::MAX || second_best_dist >= f64::MAX {
        return Err("LAMBDA search iteration limit exceeded");
    }

    Ok(((best_z, best_dist), (second_best_z, second_best_dist)))
}

/// Decorrelates the ambiguities using the LAMBDA reduction (Z-transformation).
/// Returns (z_hat, Q_z, Z_mat, L, D) where Q_z = L^T D L
fn decorrelate(a: &DVector<f64>, q: &DMatrix<f64>) -> Result<DecorrelateResult, &'static str> {
    let n = a.len();
    let mut z_mat = DMatrix::<f64>::identity(n, n);
    let mut z_hat = a.clone();

    let mut q_z = q.clone();
    for i in 0..n {
        q_z[(i, i)] += 1e-10;
    }

    let res = ldlt_lower(&q_z)?;
    let mut l = res.l;
    let mut d = res.d;

    let mut k = (n - 2) as isize;
    let mut iter = 0;
    while k >= 0 && iter < 100 {
        iter += 1;
        let k_u = k as usize;

        if apply_decorrelation_step(n, k_u, &l, &mut z_mat, &mut z_hat, &mut q_z) {
            let res = ldlt_lower(&q_z)?;
            l = res.l;
            d = res.d;
        }

        if check_swap_condition(k_u, &l, &d) {
            swap_columns(n, k_u, &mut z_mat, &mut z_hat, &mut q_z);
            let res = ldlt_lower(&q_z)?;
            l = res.l;
            d = res.d;
            k = (n - 2) as isize;
        } else {
            k -= 1;
        }
    }

    let res = ldlt_lower(&q_z)?;
    Ok(DecorrelateResult {
        z_hat,
        q_z,
        z_mat,
        l: res.l,
        d: res.d,
    })
}

fn apply_decorrelation_step(
    n: usize,
    k_u: usize,
    l: &DMatrix<f64>,
    z_mat: &mut DMatrix<f64>,
    z_hat: &mut DVector<f64>,
    q_z: &mut DMatrix<f64>,
) -> bool {
    let mut modified = false;
    for i in (k_u + 1)..n {
        let mu = l[(i, k_u)].round();
        if mu != 0.0 {
            let mut e = DMatrix::<f64>::identity(n, n);
            e[(k_u, i)] = -mu;

            *z_mat = &*z_mat * &e;
            *z_hat = e.transpose() * &*z_hat;
            *q_z = e.transpose() * &*q_z * &e;
            modified = true;
        }
    }
    modified
}

fn check_swap_condition(k_u: usize, l: &DMatrix<f64>, d: &DVector<f64>) -> bool {
    let k1 = k_u + 1;
    let delta = d[k1] + l[(k1, k_u)].powi(2) * d[k_u];
    delta < d[k_u] - 1e-6
}

fn swap_columns(
    n: usize,
    k_u: usize,
    z_mat: &mut DMatrix<f64>,
    z_hat: &mut DVector<f64>,
    q_z: &mut DMatrix<f64>,
) {
    let mut p = DMatrix::<f64>::identity(n, n);
    p.swap_columns(k_u, k_u + 1);

    *z_mat = &*z_mat * &p;
    *z_hat = p.transpose() * &*z_hat;
    *q_z = p.transpose() * &*q_z * &p;
}

fn ldlt_lower(q: &DMatrix<f64>) -> Result<LdltResult, &'static str> {
    let n = q.nrows();
    let mut l = DMatrix::<f64>::identity(n, n);
    let mut d = DVector::<f64>::zeros(n);
    let mut q_tmp = q.clone();

    for j in (0..n).rev() {
        d[j] = q_tmp[(j, j)];
        if d[j] <= 1e-18 {
            return Err("Covariance matrix is not positive definite");
        }
        for i in 0..j {
            l[(j, i)] = q_tmp[(j, i)] / d[j];
            for k in 0..=i {
                q_tmp[(i, k)] -= l[(j, i)] * q_tmp[(j, k)];
                q_tmp[(k, i)] = q_tmp[(i, k)];
            }
        }
    }
    Ok(LdltResult { l, d })
}

#[allow(clippy::too_many_arguments)]
fn search_recursive(
    k: isize,
    n: usize,
    l: &DMatrix<f64>,
    d: &DVector<f64>,
    z_hat: &DVector<f64>,
    y: &mut DVector<f64>,
    current_z: &mut DVector<f64>,
    current_dist: f64,
    best_z: &mut DVector<f64>,
    best_dist: &mut f64,
    second_best_z: &mut DVector<f64>,
    second_best_dist: &mut f64,
    iter_count: &mut usize,
    max_iters: usize,
) -> bool {
    *iter_count += 1;
    if *iter_count > max_iters {
        return true;
    }
    if current_dist >= *second_best_dist {
        return false;
    }

    if k < 0 {
        if current_dist < *best_dist {
            *second_best_dist = *best_dist;
            second_best_z.copy_from(best_z);
            *best_dist = current_dist;
            best_z.copy_from(current_z);
        } else if current_dist < *second_best_dist {
            *second_best_dist = current_dist;
            second_best_z.copy_from(current_z);
        }
        return false;
    }

    let k_u = k as usize;

    // Calculate conditional mean offset
    let mut s = 0.0;
    for j in (k_u + 1)..n {
        s += l[(j, k_u)] * y[j];
    }

    let z_cond_k = z_hat[k_u] + s; // Fixed sign: L y = z - z_hat => y_k = z_k - z_hat_k - s
    let center_z = z_cond_k.round();

    let mut offset = 0.0;
    let mut step = 1.0;
    let mut direction = if z_cond_k > center_z { 1.0 } else { -1.0 };

    loop {
        let z_test = center_z + offset;
        let y_k = z_test - z_cond_k;
        let new_dist = current_dist + (y_k * y_k) / d[k_u];

        if new_dist >= *second_best_dist {
            break;
        }

        current_z[k_u] = z_test;
        y[k_u] = y_k;

        // Removed debug print

        let aborted = search_recursive(
            k - 1,
            n,
            l,
            d,
            z_hat,
            y,
            current_z,
            new_dist,
            best_z,
            best_dist,
            second_best_z,
            second_best_dist,
            iter_count,
            max_iters,
        );
        if aborted {
            return true;
        }

        offset = step * direction;
        if direction > 0.0 {
            step += 1.0;
        }
        direction = -direction;
    }

    false
}

pub fn bootstrapping_success_rate(d: &DVector<f64>) -> f64 {
    let mut ps = 1.0;
    for &di in d.as_slice() {
        if di <= 0.0 {
            continue;
        }
        let x = 1.0 / (2.0 * f64::sqrt(di));
        ps *= libm::erf(x / std::f64::consts::SQRT_2);
    }
    ps.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambda_2d() {
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::from_row_slice(2, 2, &[6.290, 5.978, 5.978, 5.692]);

        let dec = decorrelate(&a, &q).unwrap();
        let (z_hat, z_mat, l, d) = (dec.z_hat, dec.z_mat, dec.l, dec.d);
        println!("z_hat: {:?}", z_hat);
        println!("z_mat (Z): {:?}", z_mat);
        println!("L: {:?}", l);
        println!("D: {:?}", d);

        let t_inv = z_mat.transpose().try_inverse().unwrap();
        println!("t_inv: {:?}", t_inv);

        let result = resolve_lambda(&a, &q).expect("LAMBDA should succeed");

        println!("Best a: {:?}", result.best_integers);
        println!("Second Best a: {:?}", result.second_best_integers);
        println!("Ratio: {}", result.ratio);

        assert_eq!(result.best_integers[0].fract(), 0.0);
        assert_eq!(result.best_integers[1].fract(), 0.0);

        // Expected from standard LAMBDA example: [1.0, -1.0]
        assert!(
            (result.best_integers[0] - 1.0).abs() < 1e-6,
            "Expected 1.0, got {}",
            result.best_integers[0]
        );
        assert!(
            (result.best_integers[1] - -1.0).abs() < 1e-6,
            "Expected -1.0, got {}",
            result.best_integers[1]
        );
        assert!(result.ratio > 1.0);
        assert!(result.success_rate > 0.0 && result.success_rate <= 1.0);
    }

    #[test]
    fn test_lambda_iter_limit() {
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::from_row_slice(2, 2, &[6.290, 5.978, 5.978, 5.692]);

        let result = super::resolve_lambda_inner(&a, &q, 0);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "LAMBDA search iteration limit exceeded"
        );
    }

    #[test]
    fn test_lambda_catch_gte_mutant() {
        let a = DVector::from_vec(vec![0.0, 0.0]);
        let q = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);

        // This search should take exactly 4 iterations.
        // If max_iters=4, iter_count > max_iters (4 > 4) is FALSE (succeeds).
        // If mutated to >=, (4 >= 4) is TRUE (fails).
        let result = super::resolve_lambda_inner(&a, &q, 4);
        assert!(result.is_ok());
    }

    #[test]
    fn test_lambda_empty_vector() {
        let a = DVector::zeros(0);
        let q = DMatrix::zeros(0, 0);
        let result = resolve_lambda(&a, &q);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Empty ambiguity vector");
    }

    #[test]
    fn test_ldlt_lower_known_matrix() {
        // Diagonal matrix: L = I, D = [4, 9]
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![4.0, 9.0]));
        let res = ldlt_lower(&q).unwrap();
        // L should be identity
        assert!((res.l[(0, 0)] - 1.0).abs() < 1e-15);
        assert!((res.l[(1, 0)]).abs() < 1e-15);
        assert!((res.l[(1, 1)] - 1.0).abs() < 1e-15);
        // D should be [4, 9]
        assert!((res.d[0] - 4.0).abs() < 1e-10);
        assert!((res.d[1] - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_ldlt_lower_correlated_2x2() {
        // Q = [[4, 2], [2, 9]]
        // D[1] = Q[1,1] = 9
        // L[1,0] = Q[1,0] / D[1] = 2/9
        // D[0] = Q[0,0] - L[1,0]^2 * D[1] = 4 - (4/81)*9 = 4 - 36/81 = 4 - 4/9 = 32/9 = 3.555...
        let q = DMatrix::from_row_slice(2, 2, &[4.0, 2.0, 2.0, 9.0]);
        let res = ldlt_lower(&q).unwrap();
        assert!((res.d[1] - 9.0).abs() < 1e-10);
        assert!((res.l[(1, 0)] - 2.0 / 9.0).abs() < 1e-10);
        assert!((res.d[0] - 32.0 / 9.0).abs() < 1e-10);
        // Verify Q = L^T D L
        let reconstructed = res.l.transpose() * DMatrix::from_diagonal(&res.d) * &res.l;
        for r in 0..2 {
            for c in 0..2 {
                assert!(
                    (reconstructed[(r, c)] - q[(r, c)]).abs() < 1e-10,
                    "LDL^T reconstruction failed at ({r},{c}): {} vs {}",
                    reconstructed[(r, c)],
                    q[(r, c)]
                );
            }
        }
    }

    #[test]
    fn test_ldlt_lower_not_positive_definite() {
        // Zero on diagonal
        let q = DMatrix::from_row_slice(2, 2, &[0.0, 0.0, 0.0, 1.0]);
        let res = ldlt_lower(&q);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err(), "Covariance matrix is not positive definite");
    }

    #[test]
    fn test_bootstrapping_success_rate_perfect() {
        // Very small D values -> success rate should be ~1.0
        let d = DVector::from_vec(vec![1e-10, 1e-10]);
        let sr = bootstrapping_success_rate(&d);
        assert!(
            sr > 0.9999,
            "Expected near-1.0 success rate, got {}",
            sr
        );
    }

    #[test]
    fn test_bootstrapping_success_rate_poor() {
        // Large D values -> success rate should be low
        let d = DVector::from_vec(vec![100.0, 100.0]);
        let sr = bootstrapping_success_rate(&d);
        assert!(
            sr < 0.5,
            "Expected low success rate for poor precision, got {}",
            sr
        );
        assert!(sr >= 0.0);
    }

    #[test]
    fn test_bootstrapping_success_rate_nonpositive_diagonal() {
        // Zero and negative D elements should be skipped
        let d = DVector::from_vec(vec![-1.0, 0.0, 1e-4]);
        let sr = bootstrapping_success_rate(&d);
        assert!(sr > 0.0 && sr <= 1.0);
    }

    #[test]
    fn test_decorrelate_identity() {
        // When Q = I, decorrelate should produce Z = I, L = I, D = ones
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::identity(2, 2);
        let dec = decorrelate(&a, &q).unwrap();
        // Z should be approximately I (eigenvalue refinement adds 1e-10 to diagonal, so Q_z ~ I)
        assert!((dec.z_mat[(0, 0)] - 1.0).abs() < 1e-8);
        assert!((dec.z_mat[(1, 1)] - 1.0).abs() < 1e-8);
        assert!((dec.z_mat[(0, 1)]).abs() < 1e-8);
        assert!((dec.z_mat[(1, 0)]).abs() < 1e-8);
        // z_hat should approximately equal a (z_hat is transformed a)
        // With identity covariance, z_hat should equal a since no decorrelation needed
        // But numerically there may be small differences from the LDL^T + rounding
        assert!(
            (dec.z_hat[0] - a[0]).abs() < 1e-6,
            "z_hat[0] mismatch: {} vs {}",
            dec.z_hat[0],
            a[0]
        );
    }

    #[test]
    fn test_lambda_3d() {
        // 3D test case with moderate correlation
        let a = DVector::from_vec(vec![5.45, 3.10, 7.80]);
        // Slightly correlated covariance
        let q = DMatrix::from_row_slice(
            3,
            3,
            &[6.290, 0.100, 0.050, 0.100, 5.692, 0.080, 0.050, 0.080, 4.000],
        );
        let result = resolve_lambda(&a, &q).expect("LAMBDA should succeed");
        // Best integers should be integers
        assert!(
            (result.best_integers[0] - result.best_integers[0].round()).abs() < 1e-6,
            "Best integer 0 is not integer: {}",
            result.best_integers[0]
        );
        assert!(
            (result.best_integers[1] - result.best_integers[1].round()).abs() < 1e-6,
            "Best integer 1 is not integer: {}",
            result.best_integers[1]
        );
        assert!(
            (result.best_integers[2] - result.best_integers[2].round()).abs() < 1e-6,
            "Best integer 2 is not integer: {}",
            result.best_integers[2]
        );
        assert!(result.ratio > 1.0, "Ratio should be > 1.0, got {}", result.ratio);
        assert!(
            result.success_rate > 0.0 && result.success_rate <= 1.0,
            "Success rate out of range: {}",
            result.success_rate
        );
    }

    #[test]
    fn test_lambda_highly_correlated_2d() {
        // Highly correlated ambiguities - the classic LAMBDA test case
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::from_row_slice(2, 2, &[6.290, 5.978, 5.978, 5.692]);
        let result = resolve_lambda(&a, &q).expect("LAMBDA should succeed");
        // Best integers should be [1, -1] for this well-known example
        assert!(
            (result.best_integers[0] - 1.0).abs() < 1e-6,
            "Expected best_integers[0] = 1.0, got {}",
            result.best_integers[0]
        );
        assert!(
            (result.best_integers[1] - (-1.0)).abs() < 1e-6,
            "Expected best_integers[1] = -1.0, got {}",
            result.best_integers[1]
        );
        assert!(result.ratio > 1.0, "Ratio should be > 1.0");
    }

    #[test]
    fn test_check_swap_condition_true() {
        // When delta < d[k_u] - 1e-6, swap should be triggered
        let n = 2;
        let mut l = DMatrix::identity(n, n);
        let mut d = DVector::from_vec(vec![10.0, 1.0]);
        l[(1, 0)] = 3.0;
        // delta = d[1] + l[1,0]^2 * d[0] = 1.0 + 9 * 10 = 91.0
        // 91.0 < 10.0 - 1e-6? No, that's false.
        // Let me recalculate: d[1] + l[(1,0)]^2 * d[0] = 1 + 9*10 = 91
        // 91 < 10 - 1e-6? No
        // So I need a different case where swap IS triggered:
        // delta < d[k_u] - 1e-6
        // delta = d[k1] + l[(k1,k_u)]^2 * d[k_u]
        // For this to be less than d[k_u], we need d[k1] very small and l[(k1,k_u)] near 0
        // or negative l[(k1,k_u)]^2 doesn't help since it's squared
        // Actually, d[k1] would need to be small
        let mut l2 = DMatrix::identity(n, n);
        let d2 = DVector::from_vec(vec![10.0, 0.001]);
        l2[(1, 0)] = 0.0;
        // delta = 0.001 + 0 * 10 = 0.001
        // 0.001 < 10 - 1e-6? Yes!
        assert!(check_swap_condition(0, &l2, &d2));
    }

    #[test]
    fn test_check_swap_condition_false() {
        let n = 2;
        let l = DMatrix::identity(n, n);
        let d = DVector::from_vec(vec![1.0, 10.0]);
        // delta = d[1] + l[1,0]^2 * d[0] = 10 + 0*1 = 10
        // 10 < 1 - 1e-6? No
        assert!(!check_swap_condition(0, &l, &d));
    }
}
