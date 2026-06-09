#[cfg(test)]
mod tests {
    use crate::engine::updater::update;
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    #[test]
    fn test_loosely_coupled_jacobian() {
        let lever_arm = Vector3::new(1.0, -0.5, 2.0);
        let omega_b = Vector3::new(0.1, -0.05, 0.2);
        
        let mut state = RtkState::new(GpsTime::new(0, 0.0), Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0)), 10.0);
        state.velocity = Vector3::new(10.0, 20.0, 30.0);
        state.attitude = nalgebra::UnitQuaternion::from_euler_angles(0.1, -0.2, 0.3);
        
        // Analytical H
        let r_b_e = state.attitude.to_rotation_matrix();
        let l_e = r_b_e * lever_arm;
        let h_pos_att = -l_e.cross_matrix();
        
        let a_b = omega_b.cross(&lever_arm);
        let a_e = r_b_e * a_b;
        let h_vel_att = -a_e.cross_matrix();
        
        let h_vel_bg = r_b_e.matrix() * lever_arm.cross_matrix();

        // Numerical H
        let epsilon = 1e-6;
        for j in 0..3 {
            // Attitude
            let mut state_pos = state.clone();
            let mut dpsi_pos = Vector3::zeros();
            dpsi_pos[j] = epsilon;
            state_pos.attitude = nalgebra::UnitQuaternion::from_scaled_axis(dpsi_pos) * state_pos.attitude;
            
            let mut state_neg = state.clone();
            let mut dpsi_neg = Vector3::zeros();
            dpsi_neg[j] = -epsilon;
            state_neg.attitude = nalgebra::UnitQuaternion::from_scaled_axis(dpsi_neg) * state_neg.attitude;
            
            let pos_apc_pos = state_pos.position.vector + state_pos.attitude.to_rotation_matrix() * lever_arm;
            let pos_apc_neg = state_neg.position.vector + state_neg.attitude.to_rotation_matrix() * lever_arm;
            let num_pos_att = (pos_apc_pos - pos_apc_neg) / (2.0 * epsilon);
            
            assert!((num_pos_att - h_pos_att.column(j)).norm() < 1e-5);
            
            let v_apc_pos = state_pos.velocity + state_pos.attitude.to_rotation_matrix() * omega_b.cross(&lever_arm);
            let v_apc_neg = state_neg.velocity + state_neg.attitude.to_rotation_matrix() * omega_b.cross(&lever_arm);
            let num_vel_att = (v_apc_pos - v_apc_neg) / (2.0 * epsilon);
            
            assert!((num_vel_att - h_vel_att.column(j)).norm() < 1e-5);
        }
        
        for j in 0..3 {
            // Gyro bias (omega_b = gyro - bg)
            let mut bg_pos = Vector3::zeros();
            bg_pos[j] = epsilon;
            let omega_b_pos = omega_b - bg_pos;
            
            let mut bg_neg = Vector3::zeros();
            bg_neg[j] = -epsilon;
            let omega_b_neg = omega_b - bg_neg;
            
            let v_apc_pos = state.velocity + state.attitude.to_rotation_matrix() * omega_b_pos.cross(&lever_arm);
            let v_apc_neg = state.velocity + state.attitude.to_rotation_matrix() * omega_b_neg.cross(&lever_arm);
            let num_vel_bg = (v_apc_pos - v_apc_neg) / (2.0 * epsilon);
            
            assert!((num_vel_bg - h_vel_bg.column(j)).norm() < 1e-5);
        }
    }
}
