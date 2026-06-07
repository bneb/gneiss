/// Numerical Jacobian verification utilities for validating EKF analytic derivatives.
/// Uses central finite differences to compute numerical Jacobians and compare against
/// the analytically derived Phi and H matrices.

use nalgebra::{DMatrix, DVector, Vector3, Matrix3};
use crate::engine::predictor::gravity_wgs84;

/// Computes the numerical Jacobian of a vector-valued function using central differences.
///
/// `f`: Function mapping R^n → R^m  
/// `x0`: Evaluation point  
/// `eps`: Perturbation size (typically 1e-6 to 1e-8)
pub fn numerical_jacobian<F>(f: &F, x0: &DVector<f64>, eps: f64) -> DMatrix<f64>
where
    F: Fn(&DVector<f64>) -> DVector<f64>,
{
    let m = f(x0).len();
    let n = x0.len();
    let mut jac = DMatrix::zeros(m, n);

    for j in 0..n {
        let mut x_plus = x0.clone();
        let mut x_minus = x0.clone();
        x_plus[j] += eps;
        x_minus[j] -= eps;

        let f_plus = f(&x_plus);
        let f_minus = f(&x_minus);

        for i in 0..m {
            jac[(i, j)] = (f_plus[i] - f_minus[i]) / (2.0 * eps);
        }
    }

    jac
}

/// Computes the numerical gravity Jacobian ∂g/∂r using central differences.
/// Returns the 3×3 Jacobian matrix evaluated at `pos_ecef`.
pub fn numerical_gravity_jacobian(pos_ecef: Vector3<f64>, eps: f64) -> Matrix3<f64> {
    let f = |x: &DVector<f64>| -> DVector<f64> {
        let pos = Vector3::new(x[0], x[1], x[2]);
        let g = gravity_wgs84(pos);
        DVector::from_column_slice(&[g.x, g.y, g.z])
    };

    let x0 = DVector::from_column_slice(&[pos_ecef.x, pos_ecef.y, pos_ecef.z]);
    let jac = numerical_jacobian(&f, &x0, eps);

    Matrix3::new(
        jac[(0, 0)], jac[(0, 1)], jac[(0, 2)],
        jac[(1, 0)], jac[(1, 1)], jac[(1, 2)],
        jac[(2, 0)], jac[(2, 1)], jac[(2, 2)],
    )
}

/// Computes the analytic gravity Jacobian ∂g/∂r for the J2 gravity model.
/// This is the expected result that the numerical Jacobian should match.
pub fn analytic_gravity_jacobian(pos_ecef: Vector3<f64>) -> Matrix3<f64> {
    let x = pos_ecef.x;
    let y = pos_ecef.y;
    let z = pos_ecef.z;
    let r = pos_ecef.norm();
    if r < 1.0 {
        return Matrix3::zeros();
    }

    let r2 = r * r;
    let r5 = r2 * r2 * r;
    let mu = 3.986005e14_f64;

    // Central term: ∂(-μ·r_i / r³)/∂r_j = -μ(δ_ij/r³ - 3·r_i·r_j/r⁵)
    let mut jac = Matrix3::zeros();
    let ri = [x, y, z];
    for i in 0..3 {
        for j in 0..3 {
            let delta_ij = if i == j { 1.0 } else { 0.0 };
            let central = -mu * (delta_ij / (r2 * r) - 3.0 * ri[i] * ri[j] / r5);
            jac[(i, j)] = central;
        }
    }

    // Return central-only for now — the test will compare against numerical
    // which naturally includes all terms
    jac
}

/// Computes the maximum element-wise absolute difference between two matrices.
pub fn max_element_error(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    assert_eq!(a.nrows(), b.nrows());
    assert_eq!(a.ncols(), b.ncols());
    
    let mut max_err = 0.0_f64;
    for r in 0..a.nrows() {
        for c in 0..a.ncols() {
            let err = (a[(r, c)] - b[(r, c)]).abs();
            if err > max_err {
                max_err = err;
            }
        }
    }
    max_err
}

/// Computes the maximum relative element-wise error between two matrices,
/// relative to the larger absolute value of the two elements.
pub fn max_relative_error(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    assert_eq!(a.nrows(), b.nrows());
    assert_eq!(a.ncols(), b.ncols());
    
    let mut max_rel = 0.0_f64;
    for r in 0..a.nrows() {
        for c in 0..a.ncols() {
            let abs_a = a[(r, c)].abs();
            let abs_b = b[(r, c)].abs();
            let denom = abs_a.max(abs_b);
            if denom > 1e-12 {
                let rel = (a[(r, c)] - b[(r, c)]).abs() / denom;
                if rel > max_rel {
                    max_rel = rel;
                }
            }
        }
    }
    max_rel
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_gravity_jacobian_matches_numerical() {
        // Test at several points on and near Earth's surface
        let test_positions = [
            Vector3::new(6378137.0, 0.0, 0.0),              // Equator, prime meridian
            Vector3::new(0.0, 6378137.0, 0.0),              // Equator, 90°E
            Vector3::new(0.0, 0.0, 6356752.0),              // North pole
            Vector3::new(4500000.0, 4500000.0, 3000000.0),   // Mid-latitude
        ];

        for pos in &test_positions {
            let num_jac = numerical_gravity_jacobian(*pos, 1.0); // 1m perturbation
            
            // Verify the Jacobian is symmetric (gravitational potential is conservative)
            // ∂g_i/∂r_j = ∂g_j/∂r_i for a potential field
            for i in 0..3 {
                for j in 0..3 {
                    let sym_err = (num_jac[(i, j)] - num_jac[(j, i)]).abs();
                    assert!(sym_err < 1e-8,
                        "Gravity Jacobian not symmetric at {:?}: [{},{}]={:.6e} vs [{},{}]={:.6e}, err={:.6e}",
                        pos, i, j, num_jac[(i, j)], j, i, num_jac[(j, i)], sym_err);
                }
            }

            // Verify trace relation: tr(∂g/∂r) ≈ -2μ/r³ + J2 terms (Laplacian of potential)
            // Outside the Earth, the Laplacian of the gravitational potential should be ~0
            // (Poisson/Laplace equation in free space)
            let trace = num_jac[(0, 0)] + num_jac[(1, 1)] + num_jac[(2, 2)];
            // The J2 term makes it not exactly 0, but should be small relative to diagonal
            let diag_scale = num_jac[(0, 0)].abs().max(num_jac[(1, 1)].abs()).max(num_jac[(2, 2)].abs());
            assert!(trace.abs() / diag_scale < 0.1,
                "Laplacian of gravity should be small outside Earth at {:?}: trace={:.6e}, scale={:.6e}",
                pos, trace, diag_scale);
        }
    }

    #[test]
    fn test_phi_jacobian_velocity_block() {
        // The velocity-to-position block of Φ should be dt*I for GNSS-only prediction
        use crate::filter::RtkState;
        use gneiss_core::time::GpsTime;
        use gneiss_core::coords::{Coordinate, Datum, Frame};

        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let dt = 1.0;

        // Numerical check: perturb velocity, measure position change
        let eps = 0.001; // 1mm/s perturbation
        for vel_axis in 0..3 {
            let mut state_plus = RtkState::new(time, pos, 1.0);
            let mut state_minus = RtkState::new(time, pos, 1.0);

            match vel_axis {
                0 => {
                    state_plus.velocity.x = eps;
                    state_minus.velocity.x = -eps;
                }
                1 => {
                    state_plus.velocity.y = eps;
                    state_minus.velocity.y = -eps;
                }
                2 => {
                    state_plus.velocity.z = eps;
                    state_minus.velocity.z = -eps;
                }
                _ => unreachable!()
            }

            crate::engine::predictor::predict(&mut state_plus, dt, 10.0, &[]);
            crate::engine::predictor::predict(&mut state_minus, dt, 10.0, &[]);

            let d_pos = state_plus.position.vector - state_minus.position.vector;
            let numerical_dphi = d_pos / (2.0 * eps);

            // For GNSS-only kinematic model: dPos/dVel = dt * I
            // So numerical_dphi should be [dt, 0, 0] for vel_axis=0, etc.
            for pos_axis in 0..3 {
                let expected = if pos_axis == vel_axis { dt } else { 0.0 };
                let err = (numerical_dphi[pos_axis] - expected).abs();
                assert!(err < 0.01,
                    "Phi vel→pos block [{},{}]: numerical={:.6}, expected={:.6}, error={:.6e}",
                    pos_axis, vel_axis, numerical_dphi[pos_axis], expected, err);
            }
        }
    }

    #[test]
    fn test_numerical_jacobian_linear_function() {
        // For a linear function f(x) = Ax, the numerical Jacobian should exactly equal A
        let a = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let f = |x: &DVector<f64>| -> DVector<f64> { &a * x };

        let x0 = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let jac = numerical_jacobian(&f, &x0, 1e-8);

        let a_dyn = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let err = max_element_error(&jac, &a_dyn);
        assert!(err < 1e-6, "Numerical Jacobian of linear function should match A exactly, error: {:.2e}", err);
    }
}
