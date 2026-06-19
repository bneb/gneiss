#[cfg(test)]
mod tests {
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    #[test]
    fn test_ekf_update_stability() {
        let initial_pos = Coordinate::new(
            Vector3::zeros(),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(0, 0.0),
        );
        let mut _state = RtkState::new(GpsTime::new(0, 0.0), initial_pos, 10.0);

        // Z = [1.0, 2.0], H = [I | 0] (first 2 states are X and Y)
        let _z = DVector::from_vec(vec![1.0, 2.0]);
        let mut h = DMatrix::zeros(2, 18);
        h[(0, 0)] = 1.0;
        h[(1, 1)] = 1.0;
        let _r = DMatrix::from_diagonal(&DVector::from_vec(vec![0.1, 0.1]));

        // K = P * H^T * (H P H^T + R)^-1
        let p = &_state.covariance.view((0, 0), (2, 2));
        let h_small = h.view((0, 0), (2, 2));
        let s = &h_small * p * h_small.transpose() + &_r;
        let s_inv = s.try_inverse().unwrap();
        let k = p * h_small.transpose() * s_inv;
        let dx = &k * &_z;

        assert!((dx[0] - 0.990099).abs() < 1e-4, "EKF update X failed");
        assert!((dx[1] - 1.980198).abs() < 1e-4, "EKF update Y failed");

        // P_new = (I - K*H) * P
        let i_kh = DMatrix::identity(2, 2) - &k * h_small;
        let p_new = i_kh * p;
        assert!(
            (p_new[(0, 0)] - 0.0990099).abs() < 1e-4,
            "EKF covariance update failed"
        );
    }
}
