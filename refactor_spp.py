import re

with open("crates/gneiss-rtk/src/engine/processor.rs", "r") as f:
    content = f.read()

block_to_replace = """                    tracing::warn!("SPP EKF rejected for {} epochs. Hard resetting INS to SPP.", state.consecutive_rejections);
                    state.position = pos;
                    state.velocity = nalgebra::Vector3::zeros(); // Zero out diverged velocity
                    state.accel_bias = nalgebra::Vector3::zeros(); // Biases might be corrupted, 0 is a safer prior
                    state.gyro_bias = nalgebra::Vector3::zeros();
                    // Preserve attitude as it is far better than identity
                    if crate::filter::CORE_STATE_SIZE > 15 {
                        state.rcv_clk_bias = spp_cdt;
                        state.rcv_clk_drift = 0.0;
                    }
                    state.clear_ambiguities();
                    
                    state.covariance.fill(0.0);
                    let n = crate::filter::CORE_STATE_SIZE;
                    for i in 0..6 {
                        state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
                    }
                    let att_var = (1.0f64.to_radians()).powi(2);
                    for i in 6..9 { state.covariance[(i, i)] = att_var; }
                    for i in 9..12 { state.covariance[(i, i)] = 0.01; }
                    for i in 12..n {
                        state.covariance[(i, i)] = 1e-4;
                    }
                    if crate::filter::CORE_STATE_SIZE > 15 {
                        state.covariance[(15, 15)] = 1e6;
                    }
                    state.is_reset = true;
                    state.consecutive_rejections = 0;"""

new_block = "                    handle_spp_ins_rejection(state, pos, spp_cdt);"

content = content.replace(block_to_replace, new_block)

rtk_block_to_replace = """            if state.consecutive_rejections > 5 {
                tracing::warn!("Loose coupling rejected for {} epochs. Hard resetting INS to GNSS.", state.consecutive_rejections);
                state.position = gnss_state.position;
                state.velocity = gnss_state.velocity;
                state.accel_bias = nalgebra::Vector3::zeros();
                state.gyro_bias = nalgebra::Vector3::zeros();
                // Preserve attitude as gnss_state.attitude is likely identity
                
                state.covariance.fill(0.0);
                let n = crate::filter::CORE_STATE_SIZE;
                for i in 0..6 {
                    state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
                }
                for i in 6..n {
                    state.covariance[(i, i)] = 1e-4;
                }
                state.is_reset = true;
                state.consecutive_rejections = 0;
            }"""

rtk_new_block = """            if state.consecutive_rejections > 5 {
                handle_loose_coupling_rejection(state, gnss_state);
            }"""

loose_coupling_helper = """
fn handle_loose_coupling_rejection(state: &mut RtkState, gnss_state: &RtkState) {
    tracing::warn!("Loose coupling rejected for {} epochs. Hard resetting INS to GNSS.", state.consecutive_rejections);
    state.position = gnss_state.position;
    state.velocity = gnss_state.velocity;
    state.accel_bias = nalgebra::Vector3::zeros();
    state.gyro_bias = nalgebra::Vector3::zeros();
    
    state.covariance.fill(0.0);
    let n = crate::filter::CORE_STATE_SIZE;
    for i in 0..6 {
        state.covariance[(i, i)] = if i < 3 { 100.0 } else { 10.0 };
    }
    for i in 6..n {
        state.covariance[(i, i)] = 1e-4;
    }
    state.is_reset = true;
    state.consecutive_rejections = 0;
}
"""

content = content.replace(rtk_block_to_replace, rtk_new_block)
content = content + "\n" + loose_coupling_helper

spp_ins_loose = """            if state.consecutive_rejections > 5 {
                tracing::warn!("SPP-INS EKF rejected for {} epochs. Hard resetting INS to SPP.", state.consecutive_rejections);
                state.position = gnss_state.position;
                state.velocity = gnss_state.velocity; // Zero out diverged velocity
                state.accel_bias = nalgebra::Vector3::zeros(); // Biases might be corrupted, 0 is a safer prior
                state.gyro_bias = nalgebra::Vector3::zeros();
                // Preserve attitude as it is far better than identity
                state.covariance.fill(0.0);
                for i in 0..6 { state.covariance[(i, i)] = if i < 3 { 10.0 } else { 1.0 }; }
                let att_var = (1.0f64.to_radians()).powi(2);
                for i in 6..9 { state.covariance[(i, i)] = att_var; }
                for i in 9..12 { state.covariance[(i, i)] = 0.01; }
                for i in 12..15 { state.covariance[(i, i)] = (0.1f64.to_radians()).powi(2); }
                state.consecutive_rejections = 0;
                state.is_reset = true;
            }"""

spp_ins_loose_new = """            if state.consecutive_rejections > 5 {
                handle_spp_ins_loose_rejection(state, gnss_state);
            }"""

spp_ins_loose_helper = """
fn handle_spp_ins_loose_rejection(state: &mut RtkState, gnss_state: &RtkState) {
    tracing::warn!("SPP-INS EKF rejected for {} epochs. Hard resetting INS to SPP.", state.consecutive_rejections);
    state.position = gnss_state.position;
    state.velocity = gnss_state.velocity;
    state.accel_bias = nalgebra::Vector3::zeros();
    state.gyro_bias = nalgebra::Vector3::zeros();
    state.covariance.fill(0.0);
    for i in 0..6 { state.covariance[(i, i)] = if i < 3 { 10.0 } else { 1.0 }; }
    let att_var = (1.0f64.to_radians()).powi(2);
    for i in 6..9 { state.covariance[(i, i)] = att_var; }
    for i in 9..12 { state.covariance[(i, i)] = 0.01; }
    for i in 12..15 { state.covariance[(i, i)] = (0.1f64.to_radians()).powi(2); }
    state.consecutive_rejections = 0;
    state.is_reset = true;
}
"""

content = content.replace(spp_ins_loose, spp_ins_loose_new)
content = content + "\n" + spp_ins_loose_helper


with open("crates/gneiss-rtk/src/engine/processor.rs", "w") as f:
    f.write(content)
