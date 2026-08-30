#![allow(clippy::unwrap_used)]

use super::*;
use nalgebra::{DMatrix, DVector};

    #[test]
    fn test_lambda_2d() {
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::from_row_slice(2, 2, &[6.290, 5.978, 5.978, 5.692]);

        let dec = decorrelate(&a, &q).unwrap();
        let (z_hat, z_mat, l, d) = (dec.z_hat, dec.z_mat, dec.l, dec.d);
        println!("z_hat: {:?}", z_hat);
        println!("z_mat (Z): {:?}", z_mat);
        println!("L: {:?}", l);
        println!("D: {:?}", d);

        let t_inv = z_mat.transpose().try_inverse().unwrap();
        println!("t_inv: {:?}", t_inv);

        let result = resolve_lambda(&a, &q).expect("LAMBDA should succeed");

        println!("Best a: {:?}", result.best_integers);
        println!("Second Best a: {:?}", result.second_best_integers);
        println!("Ratio: {}", result.ratio);

        assert_eq!(result.best_integers[0].fract(), 0.0);
        assert_eq!(result.best_integers[1].fract(), 0.0);

        // Expected from standard LAMBDA example: [1.0, -1.0]
        assert!(
            (result.best_integers[0] - 1.0).abs() < 1e-6,
            "Expected 1.0, got {}",
            result.best_integers[0]
        );
        assert!(
            (result.best_integers[1] - -1.0).abs() < 1e-6,
            "Expected -1.0, got {}",
            result.best_integers[1]
        );
        assert!(result.ratio > 1.0);
        assert!(result.success_rate > 0.0 && result.success_rate <= 1.0);
    }

    #[test]
    fn test_lambda_iter_limit() {
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::from_row_slice(2, 2, &[6.290, 5.978, 5.978, 5.692]);

        let result = super::resolve_lambda_inner(&a, &q, 0);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "LAMBDA search iteration limit exceeded"
        );
    }

    #[test]
    fn test_lambda_catch_gte_mutant() {
        let a = DVector::from_vec(vec![0.0, 0.0]);
        let q = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);

        // This search should take exactly 4 iterations.
        // If max_iters=4, iter_count > max_iters (4 > 4) is FALSE (succeeds).
        // If mutated to >=, (4 >= 4) is TRUE (fails).
        let result = super::resolve_lambda_inner(&a, &q, 4);
        assert!(result.is_ok());
    }

    #[test]
    fn test_lambda_empty_vector() {
        let a = DVector::zeros(0);
        let q = DMatrix::zeros(0, 0);
        let result = resolve_lambda(&a, &q);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Empty ambiguity vector");
    }

    #[test]
    fn test_ldlt_lower_known_matrix() {
        // Diagonal matrix: L = I, D = [4, 9]
        let q = DMatrix::from_diagonal(&DVector::from_vec(vec![4.0, 9.0]));
        let res = ldlt_lower(&q).unwrap();
        // L should be identity
        assert!((res.l[(0, 0)] - 1.0).abs() < 1e-15);
        assert!((res.l[(1, 0)]).abs() < 1e-15);
        assert!((res.l[(1, 1)] - 1.0).abs() < 1e-15);
        // D should be [4, 9]
        assert!((res.d[0] - 4.0).abs() < 1e-10);
        assert!((res.d[1] - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_ldlt_lower_correlated_2x2() {
        // Q = [[4, 2], [2, 9]]
        // D[1] = Q[1,1] = 9
        // L[1,0] = Q[1,0] / D[1] = 2/9
        // D[0] = Q[0,0] - L[1,0]^2 * D[1] = 4 - (4/81)*9 = 4 - 36/81 = 4 - 4/9 = 32/9 = 3.555...
        let q = DMatrix::from_row_slice(2, 2, &[4.0, 2.0, 2.0, 9.0]);
        let res = ldlt_lower(&q).unwrap();
        assert!((res.d[1] - 9.0).abs() < 1e-10);
        assert!((res.l[(1, 0)] - 2.0 / 9.0).abs() < 1e-10);
        assert!((res.d[0] - 32.0 / 9.0).abs() < 1e-10);
        // Verify Q = L^T D L
        let reconstructed = res.l.transpose() * DMatrix::from_diagonal(&res.d) * &res.l;
        for r in 0..2 {
            for c in 0..2 {
                assert!(
                    (reconstructed[(r, c)] - q[(r, c)]).abs() < 1e-10,
                    "LDL^T reconstruction failed at ({r},{c}): {} vs {}",
                    reconstructed[(r, c)],
                    q[(r, c)]
                );
            }
        }
    }

    #[test]
    fn test_ldlt_lower_not_positive_definite() {
        // Zero on diagonal
        let q = DMatrix::from_row_slice(2, 2, &[0.0, 0.0, 0.0, 1.0]);
        let res = ldlt_lower(&q);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err(), "Covariance matrix is not positive definite");
    }

    #[test]
    fn test_bootstrapping_success_rate_perfect() {
        // Very small D values -> success rate should be ~1.0
        let d = DVector::from_vec(vec![1e-10, 1e-10]);
        let sr = bootstrapping_success_rate(&d);
        assert!(
            sr > 0.9999,
            "Expected near-1.0 success rate, got {}",
            sr
        );
    }

    #[test]
    fn test_bootstrapping_success_rate_poor() {
        // Large D values -> success rate should be low
        let d = DVector::from_vec(vec![100.0, 100.0]);
        let sr = bootstrapping_success_rate(&d);
        assert!(
            sr < 0.5,
            "Expected low success rate for poor precision, got {}",
            sr
        );
        assert!(sr >= 0.0);
    }

    #[test]
    fn test_bootstrapping_success_rate_nonpositive_diagonal() {
        // Zero and negative D elements should be skipped
        let d = DVector::from_vec(vec![-1.0, 0.0, 1e-4]);
        let sr = bootstrapping_success_rate(&d);
        assert!(sr > 0.0 && sr <= 1.0);
    }

    #[test]
    fn test_decorrelate_identity() {
        // When Q = I, decorrelate should produce Z = I, L = I, D = ones
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::identity(2, 2);
        let dec = decorrelate(&a, &q).unwrap();
        // Z should be approximately I (eigenvalue refinement adds 1e-10 to diagonal, so Q_z ~ I)
        assert!((dec.z_mat[(0, 0)] - 1.0).abs() < 1e-8);
        assert!((dec.z_mat[(1, 1)] - 1.0).abs() < 1e-8);
        assert!((dec.z_mat[(0, 1)]).abs() < 1e-8);
        assert!((dec.z_mat[(1, 0)]).abs() < 1e-8);
        // z_hat should approximately equal a (z_hat is transformed a)
        // With identity covariance, z_hat should equal a since no decorrelation needed
        // But numerically there may be small differences from the LDL^T + rounding
        assert!(
            (dec.z_hat[0] - a[0]).abs() < 1e-6,
            "z_hat[0] mismatch: {} vs {}",
            dec.z_hat[0],
            a[0]
        );
    }

    #[test]
    fn test_lambda_3d() {
        // 3D test case with moderate correlation
        let a = DVector::from_vec(vec![5.45, 3.10, 7.80]);
        // Slightly correlated covariance
        let q = DMatrix::from_row_slice(
            3,
            3,
            &[6.290, 0.100, 0.050, 0.100, 5.692, 0.080, 0.050, 0.080, 4.000],
        );
        let result = resolve_lambda(&a, &q).expect("LAMBDA should succeed");
        // Best integers should be integers
        assert!(
            (result.best_integers[0] - result.best_integers[0].round()).abs() < 1e-6,
            "Best integer 0 is not integer: {}",
            result.best_integers[0]
        );
        assert!(
            (result.best_integers[1] - result.best_integers[1].round()).abs() < 1e-6,
            "Best integer 1 is not integer: {}",
            result.best_integers[1]
        );
        assert!(
            (result.best_integers[2] - result.best_integers[2].round()).abs() < 1e-6,
            "Best integer 2 is not integer: {}",
            result.best_integers[2]
        );
        assert!(result.ratio > 1.0, "Ratio should be > 1.0, got {}", result.ratio);
        assert!(
            result.success_rate > 0.0 && result.success_rate <= 1.0,
            "Success rate out of range: {}",
            result.success_rate
        );
    }

    #[test]
    fn test_lambda_highly_correlated_2d() {
        // Highly correlated ambiguities - the classic LAMBDA test case
        let a = DVector::from_vec(vec![5.45, 3.10]);
        let q = DMatrix::from_row_slice(2, 2, &[6.290, 5.978, 5.978, 5.692]);
        let result = resolve_lambda(&a, &q).expect("LAMBDA should succeed");
        // Best integers should be [1, -1] for this well-known example
        assert!(
            (result.best_integers[0] - 1.0).abs() < 1e-6,
            "Expected best_integers[0] = 1.0, got {}",
            result.best_integers[0]
        );
        assert!(
            (result.best_integers[1] - (-1.0)).abs() < 1e-6,
            "Expected best_integers[1] = -1.0, got {}",
            result.best_integers[1]
        );
        assert!(result.ratio > 1.0, "Ratio should be > 1.0");
    }

    #[test]
    fn test_check_swap_condition_true() {
        // When delta < d[k_u] - 1e-6, swap should be triggered
        let n = 2;
        let mut l = DMatrix::identity(n, n);
        let _d = DVector::from_vec(vec![10.0, 1.0]);
        l[(1, 0)] = 3.0;
        // delta = d[1] + l[1,0]^2 * d[0] = 1.0 + 9 * 10 = 91.0
        // 91.0 < 10.0 - 1e-6? No, that's false.
        // Let me recalculate: d[1] + l[(1,0)]^2 * d[0] = 1 + 9*10 = 91
        // 91 < 10 - 1e-6? No
        // So I need a different case where swap IS triggered:
        // delta < d[k_u] - 1e-6
        // delta = d[k1] + l[(k1,k_u)]^2 * d[k_u]
        // For this to be less than d[k_u], we need d[k1] very small and l[(k1,k_u)] near 0
        // or negative l[(k1,k_u)]^2 doesn't help since it's squared
        // Actually, d[k1] would need to be small
        let mut l2 = DMatrix::identity(n, n);
        let d2 = DVector::from_vec(vec![10.0, 0.001]);
        l2[(1, 0)] = 0.0;
        // delta = 0.001 + 0 * 10 = 0.001
        // 0.001 < 10 - 1e-6? Yes!
        assert!(check_swap_condition(0, &l2, &d2));
    }

    #[test]
    fn test_check_swap_condition_false() {
        let n = 2;
        let l = DMatrix::identity(n, n);
        let d = DVector::from_vec(vec![1.0, 10.0]);
        // delta = d[1] + l[1,0]^2 * d[0] = 10 + 0*1 = 10
        // 10 < 1 - 1e-6? No
        assert!(!check_swap_condition(0, &l, &d));
    }
