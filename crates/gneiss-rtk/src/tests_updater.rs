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
        let _state = RtkState::new(GpsTime::new(0, 0.0), initial_pos, 10.0);
        
        let _z = DVector::from_vec(vec![1.0, 2.0]);
        let mut _h = DMatrix::zeros(2, 18);
        _h[(0, 0)] = 1.0;
        _h[(1, 1)] = 1.0;
        let _r = DMatrix::from_diagonal(&DVector::from_vec(vec![0.1, 0.1]));
        
    }
}
