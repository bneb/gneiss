use nalgebra::{DMatrix, Matrix3, Vector3};

pub struct CarrierPhaseFactor;

impl CarrierPhaseFactor {
    pub fn jacobian(
        u: &Vector3<f64>,
        r_b_e: &Matrix3<f64>,
        l_b: &Vector3<f64>,
        lambda: f64,
    ) -> DMatrix<f64> {
        let mut h = DMatrix::zeros(1, 18);
        let skew_l = super::skew_symmetric(l_b);
        let theta_term = -u.transpose() * r_b_e * skew_l;

        // p^e
        h.fixed_view_mut::<1, 3>(0, 0).copy_from(&u.transpose());

        // theta
        h.fixed_view_mut::<1, 3>(0, 6).copy_from(&theta_term);

        // cb_m
        h[(0, 15)] = 1.0;

        // N^s
        h[(0, 17)] = lambda;

        h
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Matrix3, Vector3};

    #[test]
    fn test_carrier_phase_jacobian() {
        let u = Vector3::new(0.6, 0.8, 0.0);
        let r_b_e = Matrix3::identity();
        let l_b = Vector3::new(1.0, 0.0, 0.0);
        let lambda = 0.19;

        let h = CarrierPhaseFactor::jacobian(&u, &r_b_e, &l_b, lambda);

        assert_eq!(h.nrows(), 1);
        assert_eq!(h.ncols(), 18);

        // check p_e
        assert_eq!(h[(0, 0)], 0.6);
        assert_eq!(h[(0, 1)], 0.8);
        assert_eq!(h[(0, 2)], 0.0);

        // check v_e
        for i in 3..6 {
            assert_eq!(h[(0, i)], 0.0);
        }

        // check theta
        let expected_theta = -u.transpose() * r_b_e * super::super::skew_symmetric(&l_b);
        assert_eq!(h[(0, 6)], expected_theta[0]);
        assert_eq!(h[(0, 7)], expected_theta[1]);
        assert_eq!(h[(0, 8)], expected_theta[2]);

        // check ba, bg
        for i in 9..15 {
            assert_eq!(h[(0, i)], 0.0);
        }

        // check cb_m
        assert_eq!(h[(0, 15)], 1.0);

        // check cd_m
        assert_eq!(h[(0, 16)], 0.0);

        // check N^s
        assert_eq!(h[(0, 17)], lambda);
    }
}
