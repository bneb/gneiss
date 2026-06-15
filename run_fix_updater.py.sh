cat << 'EOF' > fix_updater.py
with open("crates/gneiss-rtk/src/engine/updater_math.rs", "r") as f:
    text = f.read()

new_test = """
    #[test]
    fn test_evaluate_post_fit_outliers() {
        use crate::engine::config::EkfTuningConfig;
        let tuning = EkfTuningConfig::default();
        
        let v = DVector::from_vec(vec![1000.0, 10.0, 10.0]);
        let s = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0, 1.0, 1.0]));
        let current_z = DVector::from_vec(vec![1000.0, 10.0, 10.0]);
        let current_valid = vec![0, 1, 2];
        
        // Mock meas_types
        // We need meas_type == 3 to bypass is_abs_outlier
        let sat = gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn: 1 };
        let meas_types = vec![(sat, 1), (sat, 3), (sat, 1)];
        
        // Test 1: meas_type != 3, so v[0] = 1000 is an abs outlier
        let (outlier, _) = evaluate_post_fit_outliers(&v, &s, &current_z, &current_valid, Some(&meas_types), 50.0, true, &tuning);
        assert_eq!(outlier, Some(0)); // Returns immediately on abs outlier
        
        // Test 2: meas_type == 3, so it's NOT an abs outlier, but it will be a ratio outlier
        let v_t2 = DVector::from_vec(vec![0.0, 1000.0, 10.0]); // v[1] corresponds to meas_types[1] which is type 3
        let current_z_t2 = DVector::from_vec(vec![0.0, 1000.0, 10.0]);
        let (outlier2, ratio) = evaluate_post_fit_outliers(&v_t2, &s, &current_z_t2, &current_valid, Some(&meas_types), 50.0, true, &tuning);
        assert_eq!(outlier2, Some(1)); // Because it has a massive ratio, but was NOT flagged as abs outlier
        assert!(ratio > 900.0);
        
        // Test 3: Multiple outliers, finds the worst ratio
        let v_t3 = DVector::from_vec(vec![0.0, 50.0, 20.0]); // None are abs outliers (assume thresh is high enough, or meas_types logic)
        let s_t3 = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0, 25.0, 1.0])); 
        // ratio 1: 50/sqrt(25) = 10
        // ratio 2: 20/sqrt(1) = 20 (worse!)
        let current_z_t3 = DVector::from_vec(vec![0.0, 50.0, 20.0]);
        let meas_types_t3 = vec![(sat, 3), (sat, 3), (sat, 3)]; // all type 3 to bypass abs
        let (outlier3, ratio3) = evaluate_post_fit_outliers(&v_t3, &s_t3, &current_z_t3, &current_valid, Some(&meas_types_t3), 1.0, true, &tuning);
        assert_eq!(outlier3, Some(2));
        assert_eq!(ratio3, 20.0);
    }
"""

text = text.replace("}\n\n", "}\n" + new_test + "\n", 1) if "test_compute_loose_coupling_innovations" in text else text + new_test

with open("crates/gneiss-rtk/src/engine/updater_math.rs", "w") as f:
    f.write(text)

EOF
python3 fix_updater.py && cargo test -p gneiss-rtk -- test_evaluate_post_fit_outliers