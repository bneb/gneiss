with open("crates/gneiss-rtk/src/engine/updater_math.rs", "r") as f:
    text = f.read()

# Replace the incorrect test
text = text.split("#[test]\n    fn test_populate_loosely_coupled_jacobian")[0]

new_test = """
    #[test]
    fn test_populate_loosely_coupled_jacobian() {
        use nalgebra::{UnitQuaternion, Vector3, DMatrix, Rotation3};
        let mut h = DMatrix::zeros(6, 15);
        let r_b_e = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
        let lever_arm = Vector3::new(1.0, 2.0, 3.0);
        let omega_b = Vector3::new(0.1, 0.05, -0.1);
        
        populate_loosely_coupled_jacobian(&mut h, &r_b_e, &lever_arm, &omega_b);
        
        let l_e = r_b_e * lever_arm;
        let expected_h_pos_att = -l_e.cross_matrix();
        let a_e = r_b_e * omega_b.cross(&lever_arm);
        let expected_h_vel_att = -a_e.cross_matrix();
        let expected_h_vel_bg = r_b_e.to_rotation_matrix().matrix() * lever_arm.cross_matrix();
        
        assert_eq!(h.view((0, 6), (3, 3)).clone_owned(), expected_h_pos_att);
        assert_eq!(h.view((3, 6), (3, 3)).clone_owned(), expected_h_vel_att);
        assert_eq!(h.view((3, 12), (3, 3)).clone_owned(), expected_h_vel_bg);
    }
"""

text += new_test + "\n}\n"

with open("crates/gneiss-rtk/src/engine/updater_math.rs", "w") as f:
    f.write(text)

