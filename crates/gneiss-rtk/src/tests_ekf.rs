#[cfg(test)]
mod tests {
    use nalgebra::{DMatrix, DVector, Vector3};
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_ekf_update_stability() {
        let initial_pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let mut state = RtkState::new(GpsTime::new(0, 0.0), initial_pos, 10.0);
        
        // Z = [1.0, 2.0], H = [I | 0] (first 2 states are X and Y)
        let z = DVector::from_vec(vec![1.0, 2.0]);
        let mut h = DMatrix::zeros(2, 18); 
        h[(0, 0)] = 1.0; 
        h[(1, 1)] = 1.0;
        let r = DMatrix::from_diagonal(&DVector::from_vec(vec![0.1, 0.1]));
        
        // Small subset of states to test update logic
        // P initial = diag(10, 10)
        // K = P * H^T * (H P H^T + R)^-1
        // K = 10 * 1 * (10 + 0.1)^-1 = 10 / 10.1 = 0.990099
        // dx = K * Z = 0.99 * 1.0 = 0.990099
        
        use crate::engine::updater::update;
    }
}
