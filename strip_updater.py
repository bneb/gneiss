import re

with open("crates/gneiss-rtk/src/engine/updater.rs", "r") as f:
    text = f.read()

# We need to replace the consts at the top:
text = re.sub(r'const CP_PRE_FIT_CHI2_THRESHOLD: f64 = 100\.0;\n', '', text)
text = re.sub(r'const DOPPLER_PRE_FIT_CHI2_THRESHOLD: f64 = 50\.0;\n', '', text)
text = re.sub(r'const PR_PRE_FIT_CHI2_MULTIPLIER: f64 = 25\.0;\n', '', text)

# Now we need to remove specific functions:
def remove_func(text, func_name):
    # Find the start
    start_match = re.search(r'(pub )?fn ' + func_name + r'\(', text)
    if not start_match:
        return text
    start_idx = start_match.start()
    
    # Find the end by counting braces
    brace_count = 0
    in_func = False
    end_idx = start_idx
    for i in range(start_idx, len(text)):
        if text[i] == '{':
            brace_count += 1
            in_func = True
        elif text[i] == '}':
            brace_count -= 1
            if brace_count == 0 and in_func:
                end_idx = i + 1
                break
    return text[:start_idx] + text[end_idx:]

funcs_to_remove = [
    "apply_joseph_covariance_update",
    "enforce_symmetry",
    "get_pre_fit_threshold",
    "check_pre_fit_residual",
    "filter_pre_fit_residuals",
    "evaluate_post_fit_outliers",
    "compute_scalar_thresholds",
    "huber_scale_covariance",
    "populate_loosely_coupled_jacobian"
]

for func in funcs_to_remove:
    text = remove_func(text, func)

# For update_loosely_coupled, we need to replace the inline innovation computation
# and inline jacobian population.
# 1. z = compute_loose_coupling_innovations(...)
text = re.sub(
    r'let l_e = r_b_e \* lever_arm;\n\s*let pos_apc = state\.position\.vector \+ l_e;\n\s*let v_apc = state\.velocity \+ r_b_e \* omega_b\.cross\(&lever_arm\);\n\n\s*let mut z = DVector::zeros\(6\);\n\s*z\.rows_mut\(0, 3\)\.copy_from\(&\(gnss_state\.position\.vector - pos_apc\)\);\n\s*z\.rows_mut\(3, 3\)\.copy_from\(&\(gnss_state\.velocity - v_apc\)\);',
    r'let z = compute_loose_coupling_innovations(r_b_e.matrix(), &state.position.vector, &state.velocity, &gnss_state.position.vector, &gnss_state.velocity, &lever_arm, &omega_b);',
    text
)

# 2. populate_loosely_coupled_jacobian
text = re.sub(
    r'let h_pos_att = -l_e\.cross_matrix\(\);\n\s*let a_e = r_b_e \* omega_b\.cross\(&lever_arm\);\n\s*let h_vel_att = -a_e\.cross_matrix\(\);\n\s*let h_vel_bg = r_b_e\.matrix\(\) \* lever_arm\.cross_matrix\(\);\n\s*for i in 0\.\.3 \{\n\s*for j in 0\.\.3 \{\n\s*h_mat\[\(i, 6 \+ j\)\] = h_pos_att\[\(i, j\)\];\n\s*h_mat\[\(3 \+ i, 6 \+ j\)\] = h_vel_att\[\(i, j\)\];\n\s*h_mat\[\(3 \+ i, 12 \+ j\)\] = h_vel_bg\[\(i, j\)\];\n\s*\}\n\s*\}',
    r'populate_loosely_coupled_jacobian(&mut h_mat, &state.attitude, &lever_arm, &omega_b);',
    text
)

with open("crates/gneiss-rtk/src/engine/updater.rs", "w") as f:
    f.write(text)

