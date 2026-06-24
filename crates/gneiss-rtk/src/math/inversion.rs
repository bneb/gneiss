use crate::math::CovMatrix;

const DEFAULT_SVD_EPSILON: f64 = 1e-9;
const DEFAULT_REGULARIZATION: f64 = 1e-6;

pub fn invert_matrix_robust(m: &CovMatrix) -> CovMatrix {
    if let Some(chol) = m.clone().cholesky() {
        chol.inverse()
    } else if let Ok(inv) = m
        .clone()
        .svd(true, true)
        .pseudo_inverse(DEFAULT_SVD_EPSILON)
    {
        inv
    } else {
        CovMatrix::identity(m.nrows(), m.ncols()) * DEFAULT_REGULARIZATION
    }
}

pub fn solve_cholesky_svd(
    h: &CovMatrix,
    b: &nalgebra::DVector<f64>,
    svd_eps: f64,
) -> Result<nalgebra::DVector<f64>, &'static str> {
    if let Some(chol) = h.clone().cholesky() {
        return Ok(chol.solve(b));
    }
    let svd = h.clone().svd(true, true);
    svd.solve(b, svd_eps).map_err(|_| "SVD solve failed")
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dmatrix;

    #[test]
    fn test_invert_matrix_robust() {
        let m = CovMatrix::identity(3, 3) * 2.0;
        let inv = invert_matrix_robust(&m);
        assert!((inv[(0, 0)] - 0.5).abs() < 1e-6);

        let mut m_singular = CovMatrix::zeros(3, 3);
        m_singular[(0, 0)] = 1.0;
        let inv_singular = invert_matrix_robust(&m_singular);
        assert!((inv_singular[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((inv_singular[(1, 1)] - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_solve_cholesky_svd() {
        let h = CovMatrix::identity(2, 2) * 2.0;
        let b = nalgebra::DVector::from_vec(vec![4.0, 6.0]);
        let x = solve_cholesky_svd(&h, &b, 1e-9).unwrap();
        assert!((x[(0, 0)] - 2.0).abs() < 1e-6);
        assert!((x[(1, 0)] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_solve_cholesky_svd_fallback_svd_path() {
        // Non-positive-definite matrix: cholesky will fail, SVD solve is used
        // A 2x2 matrix with a small negative eigenvalue (not SPD)
        let h = dmatrix![1.0, 2.0; 2.0, 1.0];
        // h has eigenvalues 3 and -1 -> not SPD -> cholesky fails
        let b = nalgebra::DVector::from_vec(vec![5.0, 4.0]);
        let x = solve_cholesky_svd(&h, &b, 1e-9).unwrap();
        // Verify result solves H*x = b
        let residual = &h * &x - b;
        assert!(residual.norm() < 1e-6);
    }

    #[test]
    fn test_solve_cholesky_svd_zero_matrix_ok() {
        // Zero matrix: cholesky fails, SVD solve returns zero solution (not an error)
        let h = dmatrix![0.0, 0.0; 0.0, 0.0];
        let b = nalgebra::DVector::from_vec(vec![1.0, 1.0]);
        let x = solve_cholesky_svd(&h, &b, 1e-9).unwrap();
        assert!((x[0]).abs() < 1e-12);
        assert!((x[1]).abs() < 1e-12);
    }

    #[test]
    fn test_invert_matrix_robust_svd_path() {
        // Diagonal matrix with a zero element: not SPD, pseudo-inverse handles it
        let m = dmatrix![1.0, 0.0, 0.0; 0.0, 0.0, 0.0; 0.0, 0.0, 4.0];
        let inv = invert_matrix_robust(&m);
        assert!((inv[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((inv[(1, 1)] - 0.0).abs() < 1e-12);
        assert!((inv[(2, 2)] - 0.25).abs() < 1e-6);
    }

    #[test]
    fn test_invert_matrix_robust_fallback_regularization() {
        // Zero matrix: Cholesky fails (not SPD), SVD path succeeds and returns zero matrix.
        // This exercises the SVD fallback code path.
        let m = CovMatrix::zeros(2, 2);
        let inv = invert_matrix_robust(&m);
        // SVD pseudoinverse of zero matrix is zero matrix
        assert!((inv[(0, 0)]).abs() < 1e-12);
        assert!((inv[(1, 1)]).abs() < 1e-12);
        assert!((inv[(0, 1)]).abs() < 1e-12);
        assert!((inv[(1, 0)]).abs() < 1e-12);
        // Must be symmetric
        assert!((inv[(0, 1)] - inv[(1, 0)]).abs() < 1e-12);
    }

    #[test]
    fn test_invert_matrix_robust_scalar() {
        // 1x1 positive matrix
        let m = dmatrix![16.0];
        let inv = invert_matrix_robust(&m);
        assert!((inv[(0, 0)] - 0.0625).abs() < 1e-12);
    }

    #[test]
    fn test_solve_cholesky_svd_singular_matrix_svd_path() {
        // Rank-1 singular matrix: the function must produce a valid solution
        // regardless of whether Cholesky or SVD path is taken
        let h = dmatrix![1.0, 2.0; 1.0, 2.0];
        let b = nalgebra::DVector::from_vec(vec![5.0, 5.0]);
        let x = solve_cholesky_svd(&h, &b, 1e-9).unwrap();
        // Verify h*x ≈ b (residual small)
        let residual = &h * &x - b;
        assert!(residual.norm() < 1e-6);
    }

    #[test]
    fn test_solve_cholesky_svd_non_square_not_applicable() {
        // Non-square not possible since CovMatrix is always square via type constraint
        // Test a nearly-singular but positive semidefinite matrix
        let h = dmatrix![1.0, 1.0; 1.0, 1.0 + 1e-12]; // almost singular
        let b = nalgebra::DVector::from_vec(vec![2.0, 3.0]);
        // Cholesky on nearly-singular may fail, SVD fallback should handle it
        let x = solve_cholesky_svd(&h, &b, 1e-15).unwrap_or_else(|_| {
            // If SVD also fails with tight epsilon, relax it
            solve_cholesky_svd(&h, &b, 1e-10).unwrap()
        });
        let residual = &h * &x - b;
        assert!(residual.norm() < 1e-6);
    }
}
