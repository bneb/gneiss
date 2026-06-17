import re

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    orig = f.read()

# Replace test_assemble_matrices
t1 = """    #[test]
    fn test_assemble_matrices() {
        let m1 = FgMeasurement {
            res: 1.5,
            h_row: DVector::from_element(3, 1.0),
            weight: 2.0,
            raw_var: 0.5,
            is_phase: false,
            sat: None,
        };
        let m2 = FgMeasurement {
            res: 2.5,
            h_row: DVector::from_element(3, 2.0),
            weight: 3.0,
            raw_var: 0.33,
            is_phase: true,
            sat: None,
        };
        let meas = vec![m1, m2];
        let (h, z, r) = assemble_matrices(&meas, 3);
        
        assert_eq!(h.nrows(), 2);
        assert_eq!(h.ncols(), 3);
        assert_eq!(h[(0,0)], 1.0);
        assert_eq!(h[(1,2)], 2.0);
        
        assert_eq!(z.len(), 2);
        assert_eq!(z[0], 1.5);
        assert_eq!(z[1], 2.5);
        
        assert_eq!(r.nrows(), 2);
        assert_eq!(r.ncols(), 2);
        assert_eq!(r[(0,0)], 2.0);
        assert_eq!(r[(1,1)], 3.0);
        assert_eq!(r[(0,1)], 0.0);
    }"""
r1 = """    #[test]
    fn test_assemble_matrices() {
        let meas = vec![
            FgMeasurement { res: 1.5, h_row: DVector::from_element(3, 1.0), weight: 2.0, raw_var: 0.5, is_phase: false, sat: None },
            FgMeasurement { res: 2.5, h_row: DVector::from_element(3, 2.0), weight: 3.0, raw_var: 0.33, is_phase: true, sat: None },
        ];
        let (h, z, r) = assemble_matrices(&meas, 3);
        assert_eq!(h.nrows(), 2); assert_eq!(h.ncols(), 3);
        assert_eq!(h[(0,0)], 1.0); assert_eq!(h[(1,2)], 2.0);
        assert_eq!(z.len(), 2); assert_eq!(z[0], 1.5); assert_eq!(z[1], 2.5);
        assert_eq!(r.nrows(), 2); assert_eq!(r.ncols(), 2);
        assert_eq!(r[(0,0)], 2.0); assert_eq!(r[(1,1)], 3.0); assert_eq!(r[(0,1)], 0.0);
    }"""
orig = orig.replace(t1, r1)

t2 = """    #[test]
    fn test_extract_and_apply_state_vector() {
        let mut state = dummy_rtk_state();
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.attitude = nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(0.1, 0.2, 0.3));
        state.accel_bias = Vector3::new(10.0, 11.0, 12.0);
        state.gyro_bias = Vector3::new(13.0, 14.0, 15.0);
        state.rcv_clk_bias = 16.0;
        state.rcv_clk_drift = 17.0;
        state.zwd = 18.0;
        state.ambiguities = vec![19.0, 20.0];
        
        let x = extract_state_vector(&state);
        assert_eq!(x.len(), CORE_STATE_SIZE + 2);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[15], 16.0);  // rcv_clk_bias
        assert_eq!(x[16], 0.0);  // isb_glo (default)
        assert_eq!(x[17], 0.0);  // isb_gal (default)
        assert_eq!(x[18], 0.0);  // isb_bds (default)
        assert_eq!(x[19], 17.0); // rcv_clk_drift
        assert_eq!(x[20], 18.0); // zwd
        assert_eq!(x[CORE_STATE_SIZE], 19.0);
        assert_eq!(x[CORE_STATE_SIZE + 1], 20.0);
        
        let mut state2 = dummy_rtk_state();
        state2.ambiguities = vec![0.0, 0.0];
        let cov = state2.covariance.clone(); apply_state_vector(&mut state2, &x, cov);
        
        assert_eq!(state2.position.vector, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(state2.velocity, Vector3::new(4.0, 5.0, 6.0));
        assert!((state2.attitude.scaled_axis() - Vector3::new(0.1, 0.2, 0.3)).norm() < 1e-10);
        assert_eq!(state2.accel_bias, Vector3::new(10.0, 11.0, 12.0));
        assert_eq!(state2.gyro_bias, Vector3::new(13.0, 14.0, 15.0));
        assert_eq!(state2.rcv_clk_bias, 16.0);
        assert_eq!(state2.rcv_clk_drift, 17.0);
        assert_eq!(state2.zwd, 18.0);
        assert_eq!(state2.ambiguities, vec![19.0, 20.0]);
    }"""
r2 = """    #[test]
    fn test_extract_and_apply_state_vector() {
        let mut state = dummy_rtk_state();
        state.position.vector = Vector3::new(1.0, 2.0, 3.0); state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.attitude = nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(0.1, 0.2, 0.3));
        state.accel_bias = Vector3::new(10.0, 11.0, 12.0); state.gyro_bias = Vector3::new(13.0, 14.0, 15.0);
        state.rcv_clk_bias = 16.0; state.rcv_clk_drift = 17.0; state.zwd = 18.0; state.ambiguities = vec![19.0, 20.0];
        let x = extract_state_vector(&state);
        assert_eq!(x.len(), CORE_STATE_SIZE + 2); assert_eq!(x[0], 1.0); assert_eq!(x[15], 16.0);
        assert_eq!(x[16], 0.0); assert_eq!(x[17], 0.0); assert_eq!(x[18], 0.0);
        assert_eq!(x[19], 17.0); assert_eq!(x[20], 18.0); assert_eq!(x[CORE_STATE_SIZE], 19.0); assert_eq!(x[CORE_STATE_SIZE + 1], 20.0);
        
        let mut state2 = dummy_rtk_state(); state2.ambiguities = vec![0.0, 0.0];
        let cov = state2.covariance.clone(); apply_state_vector(&mut state2, &x, cov);
        assert_eq!(state2.position.vector, Vector3::new(1.0, 2.0, 3.0)); assert_eq!(state2.velocity, Vector3::new(4.0, 5.0, 6.0));
        assert!((state2.attitude.scaled_axis() - Vector3::new(0.1, 0.2, 0.3)).norm() < 1e-10);
        assert_eq!(state2.accel_bias, Vector3::new(10.0, 11.0, 12.0)); assert_eq!(state2.gyro_bias, Vector3::new(13.0, 14.0, 15.0));
        assert_eq!(state2.rcv_clk_bias, 16.0); assert_eq!(state2.rcv_clk_drift, 17.0); assert_eq!(state2.zwd, 18.0);
        assert_eq!(state2.ambiguities, vec![19.0, 20.0]);
    }"""
orig = orig.replace(t2, r2)

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(orig)

