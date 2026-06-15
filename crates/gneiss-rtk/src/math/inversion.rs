use crate::math::CovMatrix;
use nalgebra::DMatrix;

const DEFAULT_SVD_EPSILON: f64 = 1e-9;
const DEFAULT_REGULARIZATION: f64 = 1e-6;

pub fn invert_matrix_robust(m: &CovMatrix) -> CovMatrix {
    if let Some(chol) = m.clone().cholesky() {
        chol.inverse()
    } else if let Ok(inv) = m.clone().svd(true, true).pseudo_inverse(DEFAULT_SVD_EPSILON) {
        inv
    } else {
        CovMatrix::identity(m.nrows(), m.ncols()) * DEFAULT_REGULARIZATION
    }
}

pub fn solve_cholesky_svd(h: &CovMatrix, b: &nalgebra::DVector<f64>, svd_eps: f64) -> Result<nalgebra::DVector<f64>, &'static str> {
    if let Some(chol) = h.clone().cholesky() {
        return Ok(chol.solve(b));
    }
    let svd = h.clone().svd(true, true);
    svd.solve(b, svd_eps).map_err(|_| "SVD solve failed")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invert_matrix_robust() {
        let m = CovMatrix::identity(3, 3) * 2.0;
        let inv = invert_matrix_robust(&m);
        assert!((inv[(0, 0)] - 0.5).abs() < 1e-6);
        
        let mut m_singular = CovMatrix::zeros(3, 3);
        m_singular[(0, 0)] = 1.0;
        let inv_singular = invert_matrix_robust(&m_singular);
        assert!(inv_singular.nrows() == 3);
    }
    
    #[test]
    fn test_solve_cholesky_svd() {
        let h = CovMatrix::identity(2, 2) * 2.0;
        let b = nalgebra::DVector::from_vec(vec![4.0, 6.0]);
        let x = solve_cholesky_svd(&h, &b, 1e-9).unwrap();
        assert!((x[(0, 0)] - 2.0).abs() < 1e-6);
        assert!((x[(1, 0)] - 3.0).abs() < 1e-6);
    }
}
