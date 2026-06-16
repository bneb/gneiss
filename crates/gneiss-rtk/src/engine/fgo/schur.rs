use nalgebra::{DMatrix, DVector};

/// Computes the Schur Complement marginalization algebra.
/// 
/// We eliminate the marginalized states $x_m$ from the linear system:
/// $$ \begin{bmatrix} H_{rr} & H_{rm} \\ H_{mr} & H_{mm} \end{bmatrix} \begin{bmatrix} \Delta x_r \\ \Delta x_m \end{bmatrix} = \begin{bmatrix} b_r \\ b_m \end{bmatrix} $$
/// 
/// Solving for $\Delta x_m$ yields:
/// $$ \Delta x_m = H_{mm}^{-1} (b_m - H_{mr} \Delta x_r) $$
/// 
/// Substituting this into the first equation, we get the marginalized prior:
/// $$ H_{prior} = H_{rr} - H_{rm} H_{mm}^{-1} H_{mr} $$
/// $$ b_{prior} = b_r - H_{rm} H_{mm}^{-1} b_m $$
/// 
/// Returns `Ok((h_prior, b_prior))` on success, or an error if $H_{mm}$ is not SPD.
pub fn schur_complement(
    h_rr: &DMatrix<f64>,
    h_rm: &DMatrix<f64>,
    h_mm: &DMatrix<f64>,
    h_mr: &DMatrix<f64>,
    b_r: &DVector<f64>,
    b_m: &DVector<f64>,
) -> Result<(DMatrix<f64>, DVector<f64>), &'static str> {
    let hmm_chol = h_mm.clone().cholesky().ok_or("H_mm must be SPD")?;
    let hmm_inv = hmm_chol.inverse();
    let hmm_inv_hmr = &hmm_inv * h_mr;

    let h_prior_raw = h_rr - h_rm * &hmm_inv_hmr;
    let b_prior = b_r - h_rm * &hmm_inv * b_m;

    let h_prior = 0.5 * (&h_prior_raw + h_prior_raw.transpose());
    Ok((h_prior, b_prior))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_schur_complement_basic() {
        let h_rr = DMatrix::from_diagonal_element(2, 2, 4.0);
        let h_rm = DMatrix::from_element(2, 2, 1.0);
        let h_mm = DMatrix::from_diagonal_element(2, 2, 2.0);
        let h_mr = DMatrix::from_element(2, 2, 1.0);
        
        let b_r = DVector::from_element(2, 2.0);
        let b_m = DVector::from_element(2, 4.0);

        let (h_prior, b_prior) = schur_complement(&h_rr, &h_rm, &h_mm, &h_mr, &b_r, &b_m).unwrap();

        assert!((h_prior[(0, 0)] - 3.0).abs() < EPSILON);
        assert!((h_prior[(0, 1)] - (-1.0)).abs() < EPSILON);
        assert!((h_prior[(1, 0)] - (-1.0)).abs() < EPSILON);
        assert!((h_prior[(1, 1)] - 3.0).abs() < EPSILON);

        assert!((b_prior[0] - (-2.0)).abs() < EPSILON);
        assert!((b_prior[1] - (-2.0)).abs() < EPSILON);
    }

    #[test]
    fn test_schur_complement_symmetry() {
        let h_rr = DMatrix::from_diagonal_element(3, 3, 5.0);
        let mut h_mm = DMatrix::from_diagonal_element(2, 2, 3.0);
        h_mm[(0, 1)] = 0.5;
        h_mm[(1, 0)] = 0.5;
        let h_rm = DMatrix::from_row_slice(3, 2, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let h_mr = h_rm.transpose();
        let b_r = DVector::from_element(3, 1.0);
        let b_m = DVector::from_element(2, 1.0);

        let (h_prior, _) = schur_complement(&h_rr, &h_rm, &h_mm, &h_mr, &b_r, &b_m).unwrap();

        for i in 0..3 {
            for j in 0..3 {
                assert!((h_prior[(i, j)] - h_prior[(j, i)]).abs() < EPSILON);
            }
        }
    }

    #[test]
    fn test_schur_complement_ill_conditioned() {
        let h_rr = DMatrix::from_diagonal_element(2, 2, 4.0);
        let h_rm = DMatrix::from_element(2, 2, 1.0);
        let mut h_mm = DMatrix::from_element(2, 2, 1.0);
        h_mm[(0, 0)] = -1.0; // Definitely not SPD
        let h_mr = DMatrix::from_element(2, 2, 1.0);
        let b_r = DVector::from_element(2, 2.0);
        let b_m = DVector::from_element(2, 4.0);

        let result = schur_complement(&h_rr, &h_rm, &h_mm, &h_mr, &b_r, &b_m);
        assert_eq!(result.err(), Some("H_mm must be SPD"));
    }

    #[test]
    fn test_schur_complement_symmetrization() {
        let h_rr = DMatrix::from_diagonal_element(2, 2, 4.0);
        let mut h_rm = DMatrix::from_element(2, 2, 1.0);
        h_rm[(0, 1)] = 2.0; 
        let h_mm = DMatrix::from_diagonal_element(2, 2, 2.0);
        let mut h_mr = h_rm.transpose();
        h_mr[(1, 0)] += 1e-12; // Introduce asymmetry
        let b_r = DVector::from_element(2, 2.0);
        let b_m = DVector::from_element(2, 4.0);

        let (h_prior, _) = schur_complement(&h_rr, &h_rm, &h_mm, &h_mr, &b_r, &b_m).unwrap();
        assert_eq!(h_prior[(0, 1)], h_prior[(1, 0)], "H_prior must be strictly symmetrized");
    }
}
