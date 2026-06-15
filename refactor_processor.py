import re

with open("crates/gneiss-rtk/src/engine/processor.rs", "r") as f:
    content = f.read()

# Define helper functions to append
helpers = """
fn apply_nhc_updates(state: &mut RtkState, config: &EngineConfig, imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>]) {
    let mut is_stationary = false;
    let mut accel_var = 1.0;
    if let Some(imu_buf) = imu_history.last() {
        if imu_buf.len() > 10 {
            let mut sum_a = nalgebra::Vector3::zeros();
            let mut sum_g = nalgebra::Vector3::zeros();
            for m in imu_buf {
                sum_a += m.accel;
                sum_g += m.gyro;
            }
            let mean_a = sum_a / (imu_buf.len() as f64);
            let mean_g = sum_g / (imu_buf.len() as f64);
            
            let mut var_a = 0.0;
            let mut var_g = 0.0;
            for m in imu_buf {
                var_a += (m.accel - mean_a).norm_squared();
                var_g += (m.gyro - mean_g).norm_squared();
            }
            var_a /= imu_buf.len() as f64;
            var_g /= imu_buf.len() as f64;
            
            if var_a < 0.05 && var_g < 0.005 {
                is_stationary = true;
            }
            accel_var = var_a.max(0.001);
        }
    }
    
    if !is_stationary && state.velocity.norm() < 0.05 {
        is_stationary = true;
    }

    if is_stationary {
        let zupt_var = (accel_var * 0.1).clamp(0.001, 0.1).sqrt();
        let _ = crate::nhc::apply_zupt(state, zupt_var, &config.tuning);
    } else {
        let omega_b = if let Some(imu_buf) = imu_history.last() {
            if let Some(last_imu) = imu_buf.last() {
                last_imu.gyro - state.gyro_bias
            } else {
                nalgebra::Vector3::zeros()
            }
        } else {
            nalgebra::Vector3::zeros()
        };
        let _ = crate::nhc::apply_nhc(state, 0.1, 0.1, &config.imu_mounting_angles, &config.imu_to_nhc_lever_arm, &omega_b, &config.tuning);
    }
}

fn handle_spp_ins_rejection(state: &mut RtkState, spp_pos: Coordinate, spp_cdt: f64) {
    tracing::warn!("SPP EKF rejected for {} epochs. Hard resetting INS to SPP.", state.consecutive_rejections);
    state.position = spp_pos;
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
    state.consecutive_rejections = 0;
}
"""

content = content + "\n" + helpers

# Refactor apply_nhc logic
nhc_block = """                let mut is_stationary = false;
                let mut accel_var = 1.0;
                if let Some(imu_buf) = self.imu_history.last() {
                    if imu_buf.len() > 10 {
                        let mut sum_a = nalgebra::Vector3::zeros();
                        let mut sum_g = nalgebra::Vector3::zeros();
                        for m in imu_buf {
                            sum_a += m.accel;
                            sum_g += m.gyro;
                        }
                        let mean_a = sum_a / (imu_buf.len() as f64);
                        let mean_g = sum_g / (imu_buf.len() as f64);
                        
                        let mut var_a = 0.0;
                        let mut var_g = 0.0;
                        for m in imu_buf {
                            var_a += (m.accel - mean_a).norm_squared();
                            var_g += (m.gyro - mean_g).norm_squared();
                        }
                        var_a /= imu_buf.len() as f64;
                        var_g /= imu_buf.len() as f64;
                        
                        if var_a < 0.05 && var_g < 0.005 {
                            is_stationary = true;
                        }
                        accel_var = var_a.max(0.001);
                    }
                }
                
                if !is_stationary && state.velocity.norm() < 0.05 {
                    is_stationary = true;
                }

                if is_stationary {
                    let zupt_var = (accel_var * 0.1).clamp(0.001, 0.1).sqrt();
                    let _ = crate::nhc::apply_zupt(state, zupt_var, &self.config.tuning);
                } else {
                    let omega_b = if let Some(imu_buf) = self.imu_history.last() {
                        if let Some(last_imu) = imu_buf.last() {
                            last_imu.gyro - state.gyro_bias
                        } else {
                            nalgebra::Vector3::zeros()
                        }
                    } else {
                        nalgebra::Vector3::zeros()
                    };
                    let _ = crate::nhc::apply_nhc(state, 0.1, 0.1, &self.config.imu_mounting_angles, &self.config.imu_to_nhc_lever_arm, &omega_b, &self.config.tuning);
                }"""

content = content.replace(nhc_block, "                apply_nhc_updates(state, &self.config, &self.imu_history);")

with open("crates/gneiss-rtk/src/engine/processor.rs", "w") as f:
    f.write(content)
