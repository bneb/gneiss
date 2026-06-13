import re

with open("crates/gneiss-rtk/src/engine/tests_updater.rs", "r") as f:
    content = f.read()

# Replace test name and logic
old_test = """    #[test]
    fn test_apply_state_correction_attitude_body_frame() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        
        // Initial attitude: 90 degrees around Z
        let q_init = nalgebra::UnitQuaternion::from_axis_angle(&nalgebra::Vector3::z_axis(), core::f64::consts::FRAC_PI_2);
        state.attitude = q_init;
        
        let mut dx = DVector::zeros(state.covariance.nrows());
        // Body frame error state: roll (rotation around X axis in body frame)
        dx[6] = 0.1; 
        
        crate::engine::updater::apply_state_correction(&mut state, &dx);
        
        // The rotation should be applied in the body frame.
        // For a body-frame roll, the global rotation becomes q_init * q_roll.
        // A vector (0, 1, 0) in the body frame should be rotated by the roll around X,
        // so it becomes (0, cos(0.1), sin(0.1)) in body.
        // Then transformed by q_init (90 deg around Z), it becomes (-cos(0.1), 0, sin(0.1)).
        
        let v_b = Vector3::new(0.0, 1.0, 0.0);
        let v_e = state.attitude * v_b;
        
        assert!((v_e.x - -f64::cos(0.1)).abs() < 1e-6);
        assert!((v_e.y - 0.0).abs() < 1e-6);
        assert!((v_e.z - f64::sin(0.1)).abs() < 1e-6);
    }"""

new_test = """    #[test]
    fn test_apply_state_correction_attitude_global_frame() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        
        // Initial attitude: 90 degrees around Z
        let q_init = nalgebra::UnitQuaternion::from_axis_angle(&nalgebra::Vector3::z_axis(), core::f64::consts::FRAC_PI_2);
        state.attitude = q_init;
        
        let mut dx = DVector::zeros(state.covariance.nrows());
        // Global frame error state: rotation around X axis in ECEF frame
        dx[6] = 0.1; 
        
        crate::engine::updater::apply_state_correction(&mut state, &dx);
        
        // The rotation should be applied in the global (ECEF) frame.
        // q_new = q_roll_global * q_init
        
        let v_b = Vector3::new(0.0, 1.0, 0.0);
        // q_init rotates (0,1,0) to (-1,0,0)
        // Then roll around X by 0.1 leaves (-1,0,0) unchanged!
        let v_e = state.attitude * v_b;
        
        assert!((v_e.x - -1.0).abs() < 1e-6);
        assert!((v_e.y - 0.0).abs() < 1e-6);
        assert!((v_e.z - 0.0).abs() < 1e-6);
    }"""

if old_test in content:
    content = content.replace(old_test, new_test)
else:
    print("Could not find old test!")

with open("crates/gneiss-rtk/src/engine/tests_updater.rs", "w") as f:
    f.write(content)
