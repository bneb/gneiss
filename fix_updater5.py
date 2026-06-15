import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

test_add = """
    #[test]
    fn test_apply_joseph_scalar() {
        let mut p = DMatrix::from_diagonal(&DVector::from_vec(vec![2.0, 3.0]));
        let mut dx = DVector::from_vec(vec![1.0, 2.0]);
        let h = DMatrix::from_row_slice(2, 2, &[
            1.0, 0.0,
            0.0, 1.0,
        ]);
        
        let i = 0;
        let r_i = 1.0;
        let v_i = 2.0;
        let s_i = 4.0;
        
        apply_joseph_scalar(&mut p, &mut dx, &h, i, r_i, v_i, s_i);
        
        // h_i = [1.0, 0.0]
        // k_i = p * h_i^T / s_i = [2.0, 0.0]^T / 4.0 = [0.5, 0.0]^T
        // dx += k_i * v_i = [1.0, 2.0] + [0.5, 0.0]*2.0 = [2.0, 2.0]
        assert!((dx[0] - 2.0).abs() < 1e-9);
        assert!((dx[1] - 2.0).abs() < 1e-9);
        
        // i_kh = I - k_i * h_i = diag(1, 1) - [0.5, 0.0]^T * [1.0, 0.0] 
        // = diag(1, 1) - [[0.5, 0.0], [0.0, 0.0]] = [[0.5, 0.0], [0.0, 1.0]]
        
        // p_new = i_kh * p * i_kh^T + k_i * r_i * k_i^T
        // i_kh * p * i_kh^T = [[0.5, 0.0], [0.0, 1.0]] * [[2.0, 0.0], [0.0, 3.0]] * [[0.5, 0.0], [0.0, 1.0]]^T
        // = [[0.5*2*0.5, 0.0], [0.0, 1.0*3*1.0]] = [[0.5, 0.0], [0.0, 3.0]]
        // k_i * r_i * k_i^T = [0.5, 0.0]^T * 1.0 * [0.5, 0.0] = [[0.25, 0.0], [0.0, 0.0]]
        // p_new = [[0.75, 0.0], [0.0, 3.0]]
        
        assert!((p[(0, 0)] - 0.75).abs() < 1e-9);
        assert!((p[(1, 1)] - 3.0).abs() < 1e-9);
        assert!((p[(0, 1)] - 0.0).abs() < 1e-9);
        assert!((p[(1, 0)] - 0.0).abs() < 1e-9);
    }
"""

content = content.replace("mod tests {", "mod tests {" + test_add)

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
