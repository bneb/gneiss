import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

bad_tests = """#[cfg(test)]
mod missed_mutant_tests {
    use super::*;
    use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};
    use gneiss_core::sat::{SatelliteId, Constellation};

    #[test]
    fn test_evaluate_post_fit_outliers_equal_ratio() {
        let mut nu = DVector::zeros(3);
        let mut r = DMatrix::zeros(3, 3);
        let hp = DMatrix::zeros(3, 3);
        
        nu[0] = 2.0; r[(0,0)] = 2.0;
        nu[1] = 3.1622776601683795; r[(1,1)] = 5.0;
        nu[2] = 1.0; r[(2,2)] = 1.0;
        
        let current_valid = vec![0, 1, 2];
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat_id, 1), (sat_id, 1), (sat_id, 1)];
        let tuning = EkfTuningConfig::default();
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, &current_valid, Some(&meas_types), 1.0, false, &tuning);
        assert_eq!(idx, Some(0));
        assert!((val - 2.0).abs() < 1e-9);
    }
    
    #[test]
    fn test_evaluate_post_fit_outliers_exact_ratio_1() {
        let mut nu = DVector::zeros(5);
        let mut r = DMatrix::zeros(5, 5);
        let hp = DMatrix::zeros(5, 5);
        
        nu[0] = 1.0; r[(0,0)] = 1.0; 
        
        let current_valid = vec![0, 1, 2, 3, 4];
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat_id, 1), (sat_id, 1), (sat_id, 1), (sat_id, 1), (sat_id, 1)];
        let tuning = EkfTuningConfig::default();
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, &current_valid, Some(&meas_types), 1.0, false, &tuning);
        assert_eq!(idx, None);
        assert!((val - 1.0).abs() < 1e-9);
    }
    
    #[test]
    fn test_evaluate_post_fit_outliers_exact_valid_count_4() {
        let mut nu = DVector::zeros(4);
        let mut r = DMatrix::zeros(4, 4);
        let hp = DMatrix::zeros(4, 4);
        
        nu[0] = 2.0; r[(0,0)] = 1.0; 
        
        let current_valid = vec![0, 1, 2, 3];
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat_id, 1), (sat_id, 1), (sat_id, 1), (sat_id, 1)];
        let tuning = EkfTuningConfig::default();
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, &current_valid, Some(&meas_types), 1.0, false, &tuning);
        assert_eq!(idx, None);
        assert!((val - 4.0).abs() < 1e-9);
    }"""

good_tests = """#[cfg(test)]
mod missed_mutant_tests {
    use super::*;
    use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};
    use gneiss_core::sat::{SatelliteId, Constellation};

    #[test]
    fn test_evaluate_post_fit_outliers_equal_ratio() {
        let mut nu = DVector::zeros(3);
        let mut r = DMatrix::zeros(3, 3);
        let hp = DVector::zeros(3);
        
        nu[0] = 2.0; r[(0,0)] = 2.0;
        nu[1] = 3.1622776601683795; r[(1,1)] = 5.0;
        nu[2] = 1.0; r[(2,2)] = 1.0;
        
        let current_valid = vec![0, 1, 2];
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat_id, 1), (sat_id, 1), (sat_id, 1)];
        let tuning = EkfTuningConfig::default();
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, &current_valid, Some(&meas_types), 1.0, false, &tuning);
        assert_eq!(idx, Some(0));
        assert!((val - 2.0).abs() < 1e-9);
    }
    
    #[test]
    fn test_evaluate_post_fit_outliers_exact_ratio_1() {
        let mut nu = DVector::zeros(5);
        let mut r = DMatrix::zeros(5, 5);
        let hp = DVector::zeros(5);
        
        nu[0] = 1.0; r[(0,0)] = 1.0; 
        
        let current_valid = vec![0, 1, 2, 3, 4];
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat_id, 1), (sat_id, 1), (sat_id, 1), (sat_id, 1), (sat_id, 1)];
        let tuning = EkfTuningConfig::default();
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, &current_valid, Some(&meas_types), 1.0, false, &tuning);
        assert_eq!(idx, None);
        assert!((val - 1.0).abs() < 1e-9);
    }
    
    #[test]
    fn test_evaluate_post_fit_outliers_exact_valid_count_4() {
        let mut nu = DVector::zeros(4);
        let mut r = DMatrix::zeros(4, 4);
        let hp = DVector::zeros(4);
        
        nu[0] = 2.0; r[(0,0)] = 1.0; 
        
        let current_valid = vec![0, 1, 2, 3];
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat_id, 1), (sat_id, 1), (sat_id, 1), (sat_id, 1)];
        let tuning = EkfTuningConfig::default();
        let (idx, val) = evaluate_post_fit_outliers(&nu, &r, &hp, &current_valid, Some(&meas_types), 1.0, false, &tuning);
        assert_eq!(idx, None);
        assert!((val - 4.0).abs() < 1e-9);
    }"""

content = content.replace(bad_tests, good_tests)
with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
