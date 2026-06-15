import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

test_add = """
    #[test]
    fn test_evaluate_post_fit_outliers_exact_ratio_1() {
        // exact ratio 1.0
        let mut nu = DVector::zeros(5);
        let mut r = DMatrix::zeros(5, 5);
        let hp = DMatrix::zeros(5, 5);
        
        nu[0] = 1.0; r[(0,0)] = 1.0; // ratio 1.0
        
        let meas_types = [(0.into(), 1), (1.into(), 1), (2.into(), 1), (3.into(), 1), (4.into(), 1)];
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, 1.0, Some(&meas_types));
        
        // strictly > 1.0 is required, so this should return None
        assert_eq!(idx, None);
        assert!((val - 1.0).abs() < 1e-9);
    }
    
    #[test]
    fn test_evaluate_post_fit_outliers_exact_valid_count_4() {
        // valid count = 4, ratio > 1.0
        let mut nu = DVector::zeros(4);
        let mut r = DMatrix::zeros(4, 4);
        let hp = DMatrix::zeros(4, 4);
        
        nu[0] = 2.0; r[(0,0)] = 1.0; // ratio 4.0
        
        let meas_types = [(0.into(), 1), (1.into(), 1), (2.into(), 1), (3.into(), 1)];
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, 1.0, Some(&meas_types));
        
        // strictly > 4 is required, so this should return None
        assert_eq!(idx, None);
        assert!((val - 4.0).abs() < 1e-9);
    }
"""

content = content.replace("mod tests {", "mod tests {" + test_add)

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
