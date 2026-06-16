use nalgebra::{DMatrix, Matrix3, Vector3};

pub struct DopplerFactor;

impl DopplerFactor {
    pub fn jacobian(
        d: f64,
        u: &Vector3<f64>,
        v_rel: &Vector3<f64>,
        r_b_e: &Matrix3<f64>,
        l_b: &Vector3<f64>,
        v_lev_body: &Vector3<f64>,
        omega_ie_e: &Vector3<f64>,
    ) -> DMatrix<f64> {
        use super::skew_symmetric as skew;
        
        let mut h = DMatrix::zeros(1, 18);
        let skew_l_b = skew(l_b);
        let r_e_b = r_b_e.transpose();

        let h_p = (1.0 / d) * v_rel.transpose() * (Matrix3::identity() - u * u.transpose());
        let h_theta = -h_p * r_b_e * skew_l_b 
            + u.transpose() * r_b_e * (-skew(v_lev_body) + skew_l_b * skew(&(r_e_b * omega_ie_e)));
        let h_bg = u.transpose() * r_b_e * skew_l_b;

        h.fixed_view_mut::<1, 3>(0, 0).copy_from(&h_p);
        h.fixed_view_mut::<1, 3>(0, 3).copy_from(&u.transpose());
        h.fixed_view_mut::<1, 3>(0, 6).copy_from(&h_theta);
        h.fixed_view_mut::<1, 3>(0, 12).copy_from(&h_bg);
        h[(0, 16)] = 1.0;

        h
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Matrix3, RowVector3, Vector3};

    fn skew(v: &Vector3<f64>) -> Matrix3<f64> {
        Matrix3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
    }

    macro_rules! check_block {
        ($jac:expr, $col:expr, $expected:expr) => {
            for i in 0..3 {
                assert!(($jac[(0, $col + i)] - $expected[(0, i)]).abs() < 1e-10);
            }
        };
    }

    #[test]
    fn test_doppler_jacobian() {
        let u = Vector3::new(1.0, 0.0, 0.0);
        let v_rel = Vector3::new(0.0, 1.0, 2.0);
        let r_b_e = Matrix3::identity();
        let l_b = Vector3::new(0.1, 0.2, 0.3);
        let v_lev_body = Vector3::new(0.0, 0.0, 0.0);
        let omega = Vector3::new(0.0, 0.0, 7.292115e-5);
        let d = 2.0;
        
        let jac = DopplerFactor::jacobian(d, &u, &v_rel, &r_b_e, &l_b, &v_lev_body, &omega);
        assert_eq!((jac.nrows(), jac.ncols()), (1, 18));
        
        let h_p = (1.0 / d) * v_rel.transpose() * (Matrix3::identity() - u * u.transpose());
        let skew_l_b = skew(&l_b);
        let h_theta = -h_p * r_b_e * skew_l_b 
            + u.transpose() * r_b_e * (-skew(&v_lev_body) + skew_l_b * skew(&(r_b_e.transpose() * omega)));
        
        check_block!(jac, 0, h_p);
        check_block!(jac, 3, u.transpose());
        check_block!(jac, 6, h_theta);
        check_block!(jac, 9, RowVector3::<f64>::zeros());
        check_block!(jac, 12, u.transpose() * r_b_e * skew_l_b);
        
        for (i, v) in [(15, 0.0), (16, 1.0), (17, 0.0)] {
            assert!((jac[(0, i)] - v).abs() < 1e-10);
        }
    }
}
