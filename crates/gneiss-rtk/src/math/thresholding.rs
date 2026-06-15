use crate::math::CovMatrix;
use nalgebra::DVector;
use crate::math::inversion::invert_matrix_robust;

pub fn apply_huber(res: f64, var: f64, k: f64) -> f64 {
    let abs_res = res.abs();
    let threshold = k * var.sqrt();
    if abs_res <= threshold {
        1.0
    } else {
        threshold / abs_res
    }
}

pub fn huber_scale_covariance(
    p: &CovMatrix,
    r: &CovMatrix,
    z: &DVector<f64>,
    huber_threshold_sq: f64,
) -> Result<CovMatrix, &'static str> {
    let s_raw = p + r;
    let s_raw_inv = invert_matrix_robust(&s_raw);

    let mahal_sq = (z.transpose() * &s_raw_inv * z)[(0, 0)];

    if mahal_sq <= huber_threshold_sq {
        return Ok(r.clone());
    }

    let scale = mahal_sq / huber_threshold_sq;
    Ok(r * scale)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{dmatrix, dvector};

    #[test]
    fn test_apply_huber() {
        assert!((apply_huber(1.0, 1.0, 3.0) - 1.0).abs() < 1e-6);
        assert!((apply_huber(4.0, 1.0, 3.0) - 0.75).abs() < 1e-6);
    }
    
    #[test]
    fn test_huber_scale_covariance() {
        let p = dmatrix![1.0];
        let r = dmatrix![1.0];
        let z = dvector![1.0];
        let scaled = huber_scale_covariance(&p, &r, &z, 9.0).unwrap();
        assert!((scaled[(0,0)] - 1.0).abs() < 1e-6);
        
        let z2 = dvector![4.0];
        let scaled2 = huber_scale_covariance(&p, &r, &z2, 2.0).unwrap();
        assert!((scaled2[(0,0)] - 4.0).abs() < 1e-6);
    }
}
