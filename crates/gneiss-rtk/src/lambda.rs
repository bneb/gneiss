use nalgebra::{DMatrix, DVector};

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

fn resolve_lambda_inner(a: &DVector<f64>, q: &DMatrix<f64>, max_iters: usize) -> Result<LambdaResult, &'static str> {
    let n = a.len();
    if n == 0 { return Err("Empty ambiguity vector"); }

    // 1. Decorrelation (Z-transformation)
    let dec = decorrelate(a, q)?;
    let (z_hat, z_mat, l, d) = (dec.z_hat, dec.z_mat, dec.l, dec.d);

    // 2. Search in the transformed space
    let mut best_z = DVector::zeros(n);
    let mut best_dist = f64::MAX;
    let mut second_best_z = DVector::zeros(n);
    let mut second_best_dist = f64::MAX;

    let mut current_z = DVector::zeros(n);
    let mut iter_count = 0;

    // Workspace for the search
    let mut y = DVector::zeros(n);

    search_recursive(
        (n - 1) as isize,
        n,
        &l,
        &d,
        &z_hat,
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

    if iter_count > max_iters {
        return Err("LAMBDA search iteration limit exceeded");
    }

    if best_dist == f64::MAX || second_best_dist == f64::MAX {
        return Err("LAMBDA search iteration limit exceeded");
    }

    // 3. Back-transformation to original space: a = Z^-T * z
    // z_hat = Z^T * a => a = (Z^T)^-1 * z_hat
    let t_inv = z_mat.transpose().try_inverse().ok_or("Transformation matrix inversion failed")?;
    let best_a = &t_inv * &best_z;
    let second_best_a = &t_inv * &second_best_z;

    // 4. Success Rate calculation
    let success_rate = bootstrapping_success_rate(&d);

    let safe_best_dist = if best_dist < 1e-12 { 1e-12 } else { best_dist };
    let ratio = second_best_dist / safe_best_dist;

    Ok(LambdaResult {
        best_integers: best_a,
        second_best_integers: second_best_a,
        ratio,
        success_rate,
    })
}

/// Decorrelates the ambiguities using the LAMBDA reduction (Z-transformation).
/// Returns (z_hat, Q_z, Z_mat, L, D) where Q_z = L^T D L
fn decorrelate(a: &DVector<f64>, q: &DMatrix<f64>) -> Result<DecorrelateResult, &'static str> {
    let n = a.len();
    let mut z_mat = DMatrix::<f64>::identity(n, n);
    let mut z_hat = a.clone();
    
    let mut q_z = q.clone();
    for i in 0..n { q_z[(i, i)] += 1e-10; }

    let mut l;
    let mut d;

    let mut k = (n - 2) as isize;
    let mut iter = 0;
    while k >= 0 && iter < 100 {
        iter += 1;
        let k_u = k as usize;
        let k1 = k_u + 1;

        let res = ldlt_lower(&q_z)?;
        l = res.l;
        d = res.d;

        let mut modified = false;
        for i in (k_u + 1)..n {
            let mu = l[(i, k_u)].round();
            if mu != 0.0 {
                let mut e = DMatrix::<f64>::identity(n, n);
                e[(k_u, i)] = -mu;
                
                z_mat = &z_mat * &e;
                z_hat = e.transpose() * &z_hat;
                q_z = e.transpose() * &q_z * &e;
                modified = true;
            }
        }
        
        if modified {
            let res = ldlt_lower(&q_z)?;
            l = res.l;
            d = res.d;
        }

        let delta = d[k1] + l[(k1, k_u)].powi(2) * d[k_u];
        if delta < d[k_u] - 1e-6 {
            let mut p = DMatrix::<f64>::identity(n, n);
            p.swap_columns(k_u, k1);
            
            z_mat = &z_mat * &p;
            z_hat = p.transpose() * &z_hat;
            q_z = p.transpose() * &q_z * &p;
            k = (n - 2) as isize;
        } else {
            k -= 1;
        }
    }
    
    let res = ldlt_lower(&q_z)?;
    let (l_final, d_final) = (res.l, res.d);
    Ok(DecorrelateResult { z_hat, q_z, z_mat, l: l_final, d: d_final })
}

fn ldlt_lower(q: &DMatrix<f64>) -> Result<LdltResult, &'static str> {
    let n = q.nrows();
    let mut l = DMatrix::<f64>::identity(n, n);
    let mut d = DVector::<f64>::zeros(n);
    let mut q_tmp = q.clone();

    for j in (0..n).rev() {
        d[j] = q_tmp[(j, j)];
        if d[j] <= 1e-18 { return Err("Covariance matrix is not positive definite"); }
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
    if *iter_count > max_iters { return true; }
    if current_dist >= *second_best_dist { return false; }

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

        if new_dist >= *second_best_dist { break; }

        current_z[k_u] = z_test;
        y[k_u] = y_k;
        
        // Removed debug print

        let aborted = search_recursive(
            k - 1, n, l, d, z_hat, y, current_z, new_dist,
            best_z, best_dist, second_best_z, second_best_dist, iter_count, max_iters,
        );
        if aborted { return true; }

        offset = step * direction;
        if direction > 0.0 { step += 1.0; }
        direction = -direction;
    }
    
    false
}

pub fn bootstrapping_success_rate(d: &DVector<f64>) -> f64 {
    let mut ps = 1.0;
    for &di in d.as_slice() {
        if di <= 0.0 { continue; }
        let x = 1.0 / (2.0 * f64::sqrt(di));
        ps *= libm::erf(x / std::f64::consts::SQRT_2);
    }
    ps.max(0.0).min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert!((result.best_integers[0] - 1.0).abs() < 1e-6, "Expected 1.0, got {}", result.best_integers[0]);
        assert!((result.best_integers[1] - -1.0).abs() < 1e-6, "Expected -1.0, got {}", result.best_integers[1]);
        assert!(result.ratio > 1.0);
        assert!(result.success_rate > 0.0 && result.success_rate <= 1.0);
    }

    fn test_lambda_iter_limit() {
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::from_row_slice(2, 2, &[6.290, 5.978, 5.978, 5.692]);
        
        let result = super::resolve_lambda_inner(&a, &q, 0);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "LAMBDA search iteration limit exceeded");
    }


    fn test_lambda_unreachable_candidates() {
        // If 'a' is huge, distance calculation overflows to INFINITY
        // breaking the loop immediately without hitting max iterations.
        let a = DVector::from_vec(vec![1e300, 1e300]);
        let q = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        
        let result = super::resolve_lambda_inner(&a, &q, 10000);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "LAMBDA search iteration limit exceeded");
    }


    fn test_lambda_catch_gte_mutant() {
        let a = DVector::from_vec(vec![0.0, 0.0]);
        let q = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        
        // This search should take exactly 4 iterations.
        // If max_iters=4, iter_count > max_iters (4 > 4) is FALSE (succeeds).
        // If mutated to >=, (4 >= 4) is TRUE (fails).
        let result = super::resolve_lambda_inner(&a, &q, 4);
        assert!(result.is_ok());
    }
}
