import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

# Add a test that checks >= mutation
test_add = """
    #[test]
    fn test_evaluate_post_fit_outliers_equal_ratio() {
        let mut nu = DVector::zeros(3);
        let mut r = DMatrix::zeros(3, 3);
        let hp = DMatrix::zeros(3, 3);
        
        // item 0: ratio = 4/2 = 2.0. threshold = 1.0 -> 2.0 / 1.0 = 2.0
        nu[0] = 2.0; r[(0,0)] = 2.0;
        
        // item 1: ratio = 10/5 = 2.0. threshold = 1.0 -> 2.0 / 1.0 = 2.0
        nu[1] = 3.1622776601683795; // sqrt(10)
        r[(1,1)] = 5.0;
        
        // item 2: ratio = 1/1 = 1.0
        nu[2] = 1.0; r[(2,2)] = 1.0;
        
        let meas_types = [(0.into(), 1), (1.into(), 1), (2.into(), 1)];
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, 1.0, Some(&meas_types));
        
        // It must return the FIRST one that achieved the ratio of 2.0, so idx=0
        assert_eq!(idx, Some(0));
        assert!((val - 2.0).abs() < 1e-9);
    }
"""

content = content.replace("mod tests {", "mod tests {" + test_add)

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
