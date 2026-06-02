#[cfg(test)]
mod tests {
    use crate::engine::predictor::predict;
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    #[test]
    fn test_predictor_motion() {
        let initial_pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let mut state = RtkState::new(GpsTime::new(0, 0.0), initial_pos, 1.0);
        state.velocity = Vector3::new(10.0, 0.0, 0.0);
        predict(&mut state, 1.0, 1.0, &[]);
        assert!((state.position.vector.x - 6378147.0).abs() < 1e-6);
    }
}
