#[cfg(test)]
mod tests {
    use crate::engine::updater::update;
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    #[test]
    fn test_updater_stability() {
        let initial_pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let mut state = RtkState::new(GpsTime::new(0, 0.0), initial_pos, 10.0);
        
        let z = DVector::from_vec(vec![1.0, 2.0]);
        let mut h = DMatrix::zeros(2, 18);
        h[(0, 0)] = 1.0;
        h[(1, 1)] = 1.0;
        let r = DMatrix::from_diagonal(&DVector::from_vec(vec![0.1, 0.1]));
        
        update(&mut state, &z, &h, &r).unwrap();
        
        assert!((state.position.vector.x - 0.990099).abs() < 1e-5);
        assert!((state.position.vector.y - 1.980198).abs() < 1e-5);
    }
}
