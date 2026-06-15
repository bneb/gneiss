cat << 'EOF' > fix_huber.py
with open("crates/gneiss-rtk/src/engine/updater_math.rs", "r") as f:
    text = f.read()

new_test = """
    #[test]
    fn test_huber_scale_covariance() {
        use crate::engine::config::EkfTuningConfig;
        let mut tuning = EkfTuningConfig::default();
        tuning.huber_threshold_loosely = 2.0; // huber_sq = 4.0
        
        let p = DMatrix::from_diagonal(&DVector::from_vec(vec![2.0, 2.0]));
        let r = DMatrix::from_diagonal(&DVector::from_vec(vec![2.0, 2.0]));
        // s_raw = p + r = diag(4.0, 4.0)
        // s_raw_inv = diag(0.25, 0.25)
        
        // Case 1: mahal_sq <= huber_sq
        // z = [2.0, 2.0] -> z^T * s_inv * z = 4*0.25 + 4*0.25 = 2.0 <= 4.0
        let z_t1 = DVector::from_vec(vec![2.0, 2.0]);
        let scaled_r = huber_scale_covariance(&p, &r, &z_t1, &tuning).unwrap();
        assert!((scaled_r[(0, 0)] - 2.0).abs() < 1e-9); // R remains unscaled
        
        // Case 2: mahal_sq > huber_sq
        // z = [4.0, 4.0] -> z^T * s_inv * z = 16*0.25 + 16*0.25 = 8.0 > 4.0
        // scale = 8.0 / 4.0 = 2.0
        // R_new = R * 2.0 = diag(4.0, 4.0)
        let z_t2 = DVector::from_vec(vec![4.0, 4.0]);
        let scaled_r_t2 = huber_scale_covariance(&p, &r, &z_t2, &tuning).unwrap();
        assert!((scaled_r_t2[(0, 0)] - 4.0).abs() < 1e-9);
        assert!((scaled_r_t2[(1, 1)] - 4.0).abs() < 1e-9);
        assert!((scaled_r_t2[(0, 1)] - 0.0).abs() < 1e-9);
    }
"""

text = text.replace("}\n\n", "}\n" + new_test + "\n", 1) if "test_evaluate_post_fit_outliers" in text else text + new_test

with open("crates/gneiss-rtk/src/engine/updater_math.rs", "w") as f:
    f.write(text)

EOF
python3 fix_huber.py && cargo test -p gneiss-rtk -- test_huber_scale_covariance