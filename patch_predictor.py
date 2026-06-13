with open("crates/gneiss-rtk/src/engine/predictor.rs", "r") as f:
    content = f.read()

# Fix compute_transition_matrix
content = content.replace("phi[(15, 16)] = dt;", "phi[(15, 19)] = dt;")

# Fix compute_process_noise
q_old = """    if crate::filter::CORE_STATE_SIZE > 15 {
        q[(15, 15)] = config.process_noise_cb * dt_abs;
        q[(16, 16)] = config.process_noise_cd * dt_abs;
        q[(17, 17)] = config.process_noise_zwd * dt_abs;
    }"""
q_new = """    if crate::filter::CORE_STATE_SIZE > 15 {
        q[(15, 15)] = config.process_noise_cb * dt_abs;
        q[(16, 16)] = 0.001 * dt_abs; // isb_glo noise
        q[(17, 17)] = 0.001 * dt_abs; // isb_gal noise
        q[(18, 18)] = 0.001 * dt_abs; // isb_bds noise
        q[(19, 19)] = config.process_noise_cd * dt_abs;
        q[(20, 20)] = config.process_noise_zwd * dt_abs;
    }"""
content = content.replace(q_old, q_new)

# Fix predict
p_old = """    if crate::filter::CORE_STATE_SIZE > 15 {
        x_pred[15] = state.rcv_clk_bias;
        x_pred[16] = state.rcv_clk_drift;
        x_pred[17] = state.zwd;
    }"""
p_new = """    if crate::filter::CORE_STATE_SIZE > 15 {
        x_pred[15] = state.rcv_clk_bias;
        x_pred[16] = state.isb_glo;
        x_pred[17] = state.isb_gal;
        x_pred[18] = state.isb_bds;
        x_pred[19] = state.rcv_clk_drift;
        x_pred[20] = state.zwd;
    }"""
content = content.replace(p_old, p_new)

with open("crates/gneiss-rtk/src/engine/predictor.rs", "w") as f:
    f.write(content)

