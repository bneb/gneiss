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
        assert!((p_new[(0, 0)] - 0.75).abs() < 1e-6);
        assert!((p_new[(1, 1)] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_apply_joseph_scalar_corrects_state_and_covariance() {
        // 2x2 diagonal covariance, 2-element state vector, single measurement
        let mut p = CovMatrix::from_diagonal(&dvector![4.0, 1.0]);
        let mut dx = dvector![0.0, 0.0];
        let h = dmatrix![1.0, 0.5];
        // measurement at index 0 with variance 1.0, innovation 2.0, innovation variance s_i
        let r_i = 1.0;
        let v_i = 2.0;
        // s_i = H * P * H^T + R = [1.0, 0.5] * diag(4, 1) * [1.0; 0.5] + 1.0
        //     = 4.0 + 0.25 + 1.0 = 5.25
        let s_i = 5.25;

        apply_joseph_scalar(&mut p, &mut dx, &h, 0, r_i, v_i, s_i);

        // Kalman gain K_i = P * H^T / s_i = [4.0; 0.5] / 5.25
        // dx += K_i * v_i = [4.0; 0.5] / 5.25 * 2.0
        let expected_dx_0 = 4.0 / 5.25 * 2.0;
        let expected_dx_1 = 0.5 / 5.25 * 2.0;
        assert!((dx[0] - expected_dx_0).abs() < 1e-10);
        assert!((dx[1] - expected_dx_1).abs() < 1e-10);

        // Covariance should be reduced after measurement update
        assert!(p[(0, 0)] < 4.0);
        assert!(p[(1, 1)] < 1.0);
        // Result must be symmetric
        assert!((p[(0, 1)] - p[(1, 0)]).abs() < 1e-12);
    }

    #[test]
    fn test_apply_joseph_scalar_zero_innovation() {
        // When v_i = 0, dx should remain unchanged
        let mut p = CovMatrix::identity(2, 2);
        let mut dx = dvector![1.0, 2.0];
        let h = dmatrix![1.0, 0.0];

        apply_joseph_scalar(&mut p, &mut dx, &h, 0, 1.0, 0.0, 2.0);

        // dx unchanged
        assert!((dx[0] - 1.0).abs() < 1e-10);
        assert!((dx[1] - 2.0).abs() < 1e-10);
        // Covariance still reduced
        assert!(p[(0, 0)] < 1.0);
    }

    #[test]
    fn test_apply_joseph_covariance_update_non_diagonal() {
        // Non-diagonal P and K, single measurement
        let p = dmatrix![4.0, 1.0; 1.0, 2.0];
        let k = dmatrix![0.8; 0.4];
        let h = dmatrix![1.0, 0.5];
        let r = dmatrix![2.0];

        let p_new = apply_joseph_covariance_update(&p, &k, &h, &r);
        // Result must be symmetric
        assert!((p_new[(0, 1)] - p_new[(1, 0)]).abs() < 1e-12);
        // Covariance should be reduced (measurement provides information)
        assert!(p_new[(0, 0)] < p[(0, 0)]);
        assert!(p_new[(1, 1)] < p[(1, 1)]);
    }

    #[test]
    fn test_apply_joseph_covariance_update_scalar() {
        // 1x1 case
        let p = dmatrix![4.0];
        let k = dmatrix![0.5];
        let h = dmatrix![1.0];
        let r = dmatrix![1.0];

        let p_new = apply_joseph_covariance_update(&p, &k, &h, &r);
        // I - K*H = 1.0 - 0.5 = 0.5
        // (I-KH)*P*(I-KH)^T = 0.5 * 4 * 0.5 = 1.0
        // K*R*K^T = 0.5 * 1 * 0.5 = 0.25
        // P_new = 1.0 + 0.25 = 1.25
        assert!((p_new[(0, 0)] - 1.25).abs() < 1e-10);
    }

    #[test]
    fn test_apply_joseph_covariance_update_identity() {
        // If K = 0, P should remain unchanged
        let p = dmatrix![2.0, 0.5; 0.5, 3.0];
        let k = dmatrix![0.0; 0.0];
        let h = dmatrix![1.0, 0.0];
        let r = dmatrix![1.0];

        let p_new = apply_joseph_covariance_update(&p, &k, &h, &r);
        assert!((p_new[(0, 0)] - p[(0, 0)]).abs() < 1e-10);
        assert!((p_new[(1, 1)] - p[(1, 1)]).abs() < 1e-10);
    }

    #[test]
    fn test_apply_joseph_scalar_larger_state() {
        // 4-element state with single measurement on the third element
        let mut p = CovMatrix::from_diagonal(&dvector![1.0, 1.0, 1.0, 1.0]);
        let mut dx = dvector![0.0, 0.0, 0.0, 0.0];
        // Measurement of third state element only
        let h = dmatrix![0.0, 0.0, 1.0, 0.0];
        let r_i = 0.5;
        let v_i = 3.0;
        // s_i = H*P*H^T + R = [0,0,1,0]*diag(1,1,1,1)*[0,0,1,0]^T + 0.5 = 1.0 + 0.5 = 1.5
        let s_i = 1.5;

        apply_joseph_scalar(&mut p, &mut dx, &h, 0, r_i, v_i, s_i);

        // Only the third element should be corrected
        assert!((dx[0] - 0.0).abs() < 1e-10);
        assert!((dx[1] - 0.0).abs() < 1e-10);
        // K_i[2] = P[2,2] / s_i = 1.0 / 1.5
        // dx[2] = K_i[2] * v_i = (1.0 / 1.5) * 3.0 = 2.0
        assert!((dx[2] - 2.0).abs() < 1e-10);
        assert!((dx[3] - 0.0).abs() < 1e-10);

        // Covariance of measured state should be reduced
        assert!(p[(2, 2)] < 1.0);
        // Off-diagonal terms should be non-zero (correlation introduced)
        assert!((p[(0, 2)] - p[(2, 0)]).abs() < 1e-12);
    }

    #[test]
    fn test_apply_joseph_scalar_high_innovation_variance() {
        // When R is very large relative to P, the gain should be near zero
        let mut p = CovMatrix::identity(2, 2);
        let mut dx = dvector![1.0, 2.0];
        let h = dmatrix![1.0, 0.0];
        // Large R_i -> s_i is dominated by R_i -> K_i is very small
        let r_i = 1e6;
        let v_i = 100.0;
        let s_i = 1.0 + r_i; // H*P*H^T = 1.0

        apply_joseph_scalar(&mut p, &mut dx, &h, 0, r_i, v_i, s_i);

        // dx should barely change (very small gain)
        assert!((dx[0] - 1.0).abs() < 1e-3);
        assert!((dx[1] - 2.0).abs() < 1e-10);

        // Covariance should remain almost unchanged
        assert!((p[(0, 0)] - 1.0).abs() < 1e-3);
        assert!((p[(1, 1)] - 1.0).abs() < 1e-10);
        assert!((p[(0, 1)] - p[(1, 0)]).abs() < 1e-12);
    }

    #[test]
    fn test_apply_joseph_scalar_perfect_measurement() {
        // r_i = 0.0 (perfect measurement): gain should fully commit to the innovation
        let mut p = CovMatrix::from_diagonal(&dvector![9.0, 4.0]);
        let mut dx = dvector![0.0, 0.0];
        let h = dmatrix![1.0, 0.0];
        // s_i = H*P*H^T + R = 9.0 + 0.0 = 9.0
        // K_i = P*H^T / s_i = [9.0; 0.0] / 9.0 = [1.0; 0.0]
        // dx[0] = K_i[0] * v_i = 1.0 * 5.0 = 5.0
        // After update:
        //   I - K*H = [[0,0],[0,1]], so P = [[0,0],[0,4]] (uncertainty on measured state goes to 0)
        apply_joseph_scalar(&mut p, &mut dx, &h, 0, 0.0, 5.0, 9.0);

        // dx reflects full innovation
        assert!((dx[0] - 5.0).abs() < 1e-10);
        assert!((dx[1] - 0.0).abs() < 1e-10);

        // Covariance of measured element goes to zero
        assert!((p[(0, 0)]).abs() < 1e-12);
        // Unmeasured element unchanged
        assert!((p[(1, 1)] - 4.0).abs() < 1e-10);
        assert!((p[(0, 1)] - p[(1, 0)]).abs() < 1e-12);
    }

    #[test]
    fn test_apply_joseph_scalar_sequential_updates() {
        // Two sequential scalar updates on a 3-state system
        let mut p = CovMatrix::from_diagonal(&dvector![4.0, 1.0, 9.0]);
        let mut dx = dvector![0.0, 0.0, 0.0];

        // Update 1: measure state 0 with v=2.0, r=1.0
        let h1 = dmatrix![1.0, 0.0, 0.0];
        let s1 = 4.0 + 1.0; // 5.0
        apply_joseph_scalar(&mut p, &mut dx, &h1, 0, 1.0, 2.0, s1);

        // After update 1: dx[0] = (4.0/5.0)*2.0 = 1.6
        assert!((dx[0] - 1.6).abs() < 1e-10);
        // Covariance should be reduced
        assert!(p[(0, 0)] < 4.0);
        assert!((p[(1, 1)] - 1.0).abs() < 1e-10);

        // Update 2: measure state 2 with v=3.0, r=4.0
        let h2 = dmatrix![0.0, 0.0, 1.0];
        // P is now non-diagonal after first update, but row 2 hasn't been touched
        // H*P*H^T = P[2,2] which is still 9.0 (first update didn't touch state 2)
        let s2 = 9.0 + 4.0; // 13.0
        apply_joseph_scalar(&mut p, &mut dx, &h2, 0, 4.0, 3.0, s2);

        // After update 2: dx[2] += (9.0/13.0)*3.0 = 2.0769...
        let expected_dx2 = 9.0 / 13.0 * 3.0;
        assert!((dx[2] - expected_dx2).abs() < 1e-10);
        // State 1 should be unchanged from first update
        assert!((dx[1] - 0.0).abs() < 1e-10);

        // Final covariance should be symmetric
        assert!((p[(0, 1)] - p[(1, 0)]).abs() < 1e-12);
        assert!((p[(0, 2)] - p[(2, 0)]).abs() < 1e-12);
        assert!((p[(1, 2)] - p[(2, 1)]).abs() < 1e-12);
    }

    #[test]
    fn test_apply_joseph_covariance_update_multiple_measurements() {
        // Batch update with a 2x2 measurement matrix (two simultaneous measurements)
        let p = CovMatrix::from_diagonal(&dvector![4.0, 1.0]);
        // K is 2x2: columns = states, rows = measurements
        let k = dmatrix![0.8, 0.1; 0.2, 0.5];
        // H is 2x2
        let h = dmatrix![1.0, 0.0; 0.0, 1.0];
        // R is 2x2 diagonal
        let r = dmatrix![1.0, 0.0; 0.0, 2.0];

        let p_new = apply_joseph_covariance_update(&p, &k, &h, &r);

        // Must be symmetric
        assert!((p_new[(0, 1)] - p_new[(1, 0)]).abs() < 1e-12);
        // Both diagonal elements should be reduced
        assert!(p_new[(0, 0)] < p[(0, 0)]);
        assert!(p_new[(1, 1)] < p[(1, 1)]);
    }
}
