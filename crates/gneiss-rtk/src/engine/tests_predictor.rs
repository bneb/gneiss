#[cfg(test)]
mod tests {
    use crate::filter::RtkState;
    use crate::engine::predictor;
    use gneiss_core::time::GpsTime;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::imu::ImuMeasurement;
    
    use nalgebra::{Vector3, UnitQuaternion, DVector, DMatrix};

    #[test]
    fn test_imu_prediction_rotation() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        
        let gyro = Vector3::new(0.0, 0.0, 90.0f64.to_radians());
        let accel = Vector3::new(0.0, 0.0, 0.0); 
        
        let imu_meas = ImuMeasurement::new(0, accel, gyro);
        predictor::predict(&mut state, 1.0, 0.1, &vec![imu_meas; 100]);
        
        let (_, _, yaw): (f64, f64, f64) = state.attitude.euler_angles();
        assert!((yaw.abs() - 1.570796).abs() < 0.1);
    }

    #[test]
    fn test_physics_stationary_gravity() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(0.0, 0.0, 6356752.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.attitude = UnitQuaternion::identity();
        
        let g = predictor::gravity_wgs84(state.position.vector);
        let accel = -g; 
        
        let imu_meas = ImuMeasurement::new(0, accel, Vector3::zeros());
        predictor::predict(&mut state, 1.0, 0.01, &vec![imu_meas; 100]);
        
        assert!(state.velocity.norm() < 1e-3, "Stationary IMU should not gain velocity, got {}", state.velocity.norm());
        assert!((state.position.vector - pos.vector).norm() < 1e-3);
    }

    #[test]
    fn test_physics_centrifugal_cancellation() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        
        let omega_ie = Vector3::new(0.0, 0.0, 7.2921151467e-5);
        let gravity = predictor::gravity_wgs84(state.position.vector);
        let centrifugal = omega_ie.cross(&(omega_ie.cross(&state.position.vector)));
        
        let f_e = -(gravity - centrifugal);
        let accel_body = state.attitude.inverse() * f_e;
        
        let imu_meas = ImuMeasurement::new(0, accel_body, Vector3::zeros());
        predictor::predict(&mut state, 1.0, 0.01, &vec![imu_meas; 100]);
        
        assert!(state.velocity.norm() < 1e-3, "Equatorial stationary IMU should be stable, got {}", state.velocity.norm());
    }

    #[test]
    fn test_coupling_lever_arm_to_attitude() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378137.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 0.01);
        
        state.attitude = UnitQuaternion::identity();
        let lever_arm = Vector3::new(10.0, 0.0, 0.0);
        
        let actual_rot = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), 1.0f64.to_radians());
        let actual_pos_apc = pos.vector + actual_rot * lever_arm;
        
        let sat_vec = actual_pos_apc + Vector3::new(0.0, 20000000.0, 0.0);
        
        let pred_pos_apc = pos.vector + state.attitude * lever_arm;
        let pred_range = (sat_vec - pred_pos_apc).norm();
        let obs_range = (sat_vec - actual_pos_apc).norm();
        
        let z = DVector::from_column_slice(&[obs_range - pred_range]);
        
        let e_sat = (sat_vec - pred_pos_apc).normalize();
        let h_pos = -e_sat.transpose(); 
        
        let lever_ecef = state.attitude * lever_arm;
        let h_att = lever_ecef.cross(&e_sat);
        
        let mut h = DMatrix::zeros(1, 18);
        h[(0, 0)] = h_pos[0]; h[(0, 1)] = h_pos[1]; h[(0, 2)] = h_pos[2];
        h[(0, 6)] = h_att.x; h[(0, 7)] = h_att.y; h[(0, 8)] = h_att.z;
        
        let r = DMatrix::from_element(1, 1, 0.0001);
        for i in 6..9 { state.covariance[(i, i)] = 1.0f64.to_radians().powi(2); }

        let h_t = h.transpose();
        let s = &h * &state.covariance * &h_t + &r;
        let s_inv = s.try_inverse().unwrap();
        let k = &state.covariance * &h_t * s_inv;
        let dx = &k * z;
        
        assert!(dx[8].abs() > 1e-6, "GNSS range should correct attitude through lever arm, got {}", dx[8]);
    }
}
