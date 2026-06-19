use nalgebra::{DMatrix, Matrix3, Vector3};

pub struct PseudorangeFactor;

impl PseudorangeFactor {
    pub fn jacobian(u: &Vector3<f64>, r_b_e: &Matrix3<f64>, l_b: &Vector3<f64>) -> DMatrix<f64> {
        let mut h = DMatrix::zeros(1, 17);
        let skew_l = Matrix3::new(0.0, -l_b.z, l_b.y, l_b.z, 0.0, -l_b.x, -l_b.y, l_b.x, 0.0);
        let theta_term = -u.transpose() * r_b_e * skew_l;

        // p^e
        h.fixed_view_mut::<1, 3>(0, 0).copy_from(&u.transpose());

        // theta
        h.fixed_view_mut::<1, 3>(0, 6).copy_from(&theta_term);

        // cb_m
        h[(0, 15)] = 1.0;

        h
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, Matrix3, Vector3};

    fn skew_symmetric(v: &Vector3<f64>) -> Matrix3<f64> {
        Matrix3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
    }

    fn compute_expected_h(
        u: &Vector3<f64>,
        r_b_e: &Matrix3<f64>,
        l_b: &Vector3<f64>,
    ) -> DMatrix<f64> {
        let mut h = DMatrix::zeros(1, 17);
        let skew_l = skew_symmetric(l_b);
        let theta_term = -u.transpose() * r_b_e * skew_l;

        // p^e
        h[(0, 0)] = u.x;
        h[(0, 1)] = u.y;
        h[(0, 2)] = u.z;

        // theta
        h[(0, 6)] = theta_term[(0, 0)];
        h[(0, 7)] = theta_term[(0, 1)];
        h[(0, 8)] = theta_term[(0, 2)];

        // cb_m
        h[(0, 15)] = 1.0;

        h
    }

    #[test]
    fn test_pseudorange_jacobian() {
        let u = Vector3::new(0.6, 0.8, 0.0);
        let r_b_e = Matrix3::identity();
        let l_b = Vector3::new(1.0, 2.0, 3.0);

        let h_expected = compute_expected_h(&u, &r_b_e, &l_b);
        let h_actual = PseudorangeFactor::jacobian(&u, &r_b_e, &l_b);

        assert_eq!(h_actual.nrows(), 1);
        assert_eq!(h_actual.ncols(), 17);

        let diff = &h_actual - &h_expected;
        assert!(
            diff.norm() < 1e-6,
            "Jacobian does not match expected output"
        );
    }
}
