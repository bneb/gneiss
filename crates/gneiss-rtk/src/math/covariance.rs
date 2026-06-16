use crate::math::{CovMatrix, JacMatrix};
use nalgebra::DVector;

pub fn apply_joseph_covariance_update(
    p: &CovMatrix,
    k: &CovMatrix,
    h: &JacMatrix,
    r: &CovMatrix,
) -> CovMatrix {
    let n = p.nrows();
    let i_kh = CovMatrix::identity(n, n) - k * h;
    let p_new = &i_kh * p * i_kh.transpose() + k * r * k.transpose();
    p_new.symmetric_part()
}

pub fn apply_joseph_scalar(
    p: &mut CovMatrix,
    dx: &mut DVector<f64>,
    h: &JacMatrix,
    i: usize,
    r_i: f64,
    v_i: f64,
    s_i: f64,
) {
    let mut k_i = p.clone() * h.row(i).transpose();
    k_i.unscale_mut(s_i);

    let k_h = &k_i * h.row(i);
    let n = p.nrows();
    let i_kh = CovMatrix::identity(n, n) - k_h;

    *p = &i_kh * &*p * i_kh.transpose() + &k_i * r_i * k_i.transpose();
    *p = p.symmetric_part();
    *dx += &k_i * v_i;
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{dmatrix, dvector};

    #[test]
    fn test_apply_joseph_covariance_update() {
        let p = dmatrix![2.0, 0.0; 0.0, 2.0];
        let k = dmatrix![0.5; 0.0];
        let h = dmatrix![1.0, 0.0];
        let r = dmatrix![1.0];
        
        let p_new = apply_joseph_covariance_update(&p, &k, &h, &r);
        assert!((p_new[(0,0)] - 0.75).abs() < 1e-6);
        assert!((p_new[(1,1)] - 2.0).abs() < 1e-6);
    }
}
