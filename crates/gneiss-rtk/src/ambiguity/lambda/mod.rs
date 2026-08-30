use nalgebra::{DMatrix, DVector, SymmetricEigen};

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
    resolve_lambda_inner(a, q, 200_000)
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

#[allow(clippy::type_complexity)]
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

/// Symmetrizes `q`, applies the same small unconditional diagonal epsilon
/// the previous approach always used (kept so the well-conditioned common
/// case stays numerically identical to the validated baseline), and, if
/// that isn't enough to make it (numerically) positive definite, shifts
/// its whole spectrum up by just enough to make the smallest eigenvalue
/// `MIN_EIGENVALUE_FLOOR`. A fixed diagonal epsilon on its own only helps
/// when the near-singularity is diagonal-aligned; ambiguity covariances
/// that have drifted numerically over hundreds of Kalman updates can be
/// near-singular along an arbitrary direction (e.g. several ambiguities
/// sharing almost all their uncertainty with the common position error),
/// which the eigenvalue-based shift fixes directly. Measured: applying
/// the eigenvalue shift unconditionally (dropping the unconditional
/// epsilon) shifts fix rates down by a few tenths of a percent across
/// both guard scripts on well-conditioned matrices that never needed the
/// shift in the first place -- keeping both terms avoids that regression.
const MIN_EIGENVALUE_FLOOR: f64 = 1e-9;
const BASE_DIAGONAL_EPSILON: f64 = 1e-10;

fn regularize_positive_definite(q: &DMatrix<f64>) -> DMatrix<f64> {
    let n = q.nrows();
    let sym = (q + q.transpose()) * 0.5 + DMatrix::identity(n, n) * BASE_DIAGONAL_EPSILON;
    let min_eig = SymmetricEigen::new(sym.clone()).eigenvalues.min();
    if min_eig >= MIN_EIGENVALUE_FLOOR {
        return sym;
    }
    sym + DMatrix::identity(n, n) * (MIN_EIGENVALUE_FLOOR - min_eig)
}

/// Decorrelates the ambiguities using the LAMBDA reduction (Z-transformation).
/// Returns (z_hat, Q_z, Z_mat, L, D) where Q_z = L^T D L
///
/// `q_z` is re-regularized ([`regularize_positive_definite`]) before every
/// `ldlt_lower` call, not just the first: the Z-transformations applied
/// each iteration are exact (unimodular, integer) only in theory — in
/// floating point they can nudge an already-barely-positive-definite
/// matrix back below zero along its worst direction, and that compounds
/// over the loop's many iterations on real (highly correlated,
/// multi-constellation) ambiguity sets.
fn decorrelate(a: &DVector<f64>, q: &DMatrix<f64>) -> Result<DecorrelateResult, &'static str> {
    let n = a.len();
    let mut z_mat = DMatrix::<f64>::identity(n, n);
    let mut z_hat = a.clone();

    let mut q_z = regularize_positive_definite(q);

    let res = ldlt_lower(&q_z)?;
    let mut l = res.l;
    let mut d = res.d;

    let mut k = (n - 2) as isize;
    let mut iter = 0;
    while k >= 0 && iter < 100 {
        iter += 1;
        let k_u = k as usize;

        if apply_decorrelation_step(n, k_u, &l, &mut z_mat, &mut z_hat, &mut q_z) {
            q_z = regularize_positive_definite(&q_z);
            let res = ldlt_lower(&q_z)?;
            l = res.l;
            d = res.d;
        }

        if check_swap_condition(k_u, &l, &d) {
            swap_columns(n, k_u, &mut z_mat, &mut z_hat, &mut q_z);
            q_z = regularize_positive_definite(&q_z);
            let res = ldlt_lower(&q_z)?;
            l = res.l;
            d = res.d;
            k = (n - 2) as isize;
        } else {
            k -= 1;
        }
    }

    q_z = regularize_positive_definite(&q_z);
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
mod tests;
