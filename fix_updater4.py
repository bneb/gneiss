import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

test_add = """
    #[test]
    fn test_populate_loosely_coupled_jacobian() {
        use nalgebra::{Vector3, UnitQuaternion};
        let mut h_mat = DMatrix::zeros(6, 15);
        let r_b_e = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
        let lever_arm = Vector3::new(1.0, 2.0, 3.0);
        let omega_b = Vector3::new(0.5, 0.6, 0.7);
        
        populate_loosely_coupled_jacobian(&mut h_mat, &r_b_e, &lever_arm, &omega_b);
        
        // Exact calculations
        let l_e = r_b_e * lever_arm;
        let h_pos_att = -l_e.cross_matrix();
        let a_e = r_b_e * omega_b.cross(&lever_arm);
        let h_vel_att = -a_e.cross_matrix();
        let h_vel_bg = r_b_e.to_rotation_matrix().matrix() * lever_arm.cross_matrix();
        
        assert_eq!(h_mat.view((0, 6), (3, 3)).clone_owned(), h_pos_att);
        assert_eq!(h_mat.view((3, 6), (3, 3)).clone_owned(), h_vel_att);
        assert_eq!(h_mat.view((3, 12), (3, 3)).clone_owned(), h_vel_bg);
        
        // Let's assert a specific value to catch the mutants directly
        assert!((h_mat[(0, 7)] - (-l_e[2])).abs() < 1e-9);
        assert!((h_mat[(1, 6)] - (l_e[2])).abs() < 1e-9);
    }
    
    #[test]
    fn test_compute_s_inverse_regularization() {
        // Singular matrix, all zeros
        let s = DMatrix::zeros(2, 2);
        
        // Should fallback to s + 1e-6 * I
        // Inverse of 1e-6 * I is 1e6 * I
        let inv = compute_s_inverse(&s).unwrap();
        
        assert!((inv[(0, 0)] - 1e6).abs() < 1e-5);
        assert!((inv[(1, 1)] - 1e6).abs() < 1e-5);
        assert!((inv[(0, 1)]).abs() < 1e-9);
        assert!((inv[(1, 0)]).abs() < 1e-9);
        
        // If + is mutated to -, inv is inverse of -1e-6 * I = -1e6 * I
        // If + is mutated to *, inv is inverse of 0 * I = singular -> error
        // If * is mutated to /, inv is inverse of 1e6 * I = 1e-6 * I
    }
"""

content = content.replace("mod tests {", "mod tests {" + test_add)

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
