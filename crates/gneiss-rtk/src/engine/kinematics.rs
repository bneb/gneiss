use crate::filter::RtkState;
use crate::engine::EngineConfig;
use gneiss_core::imu::ImuMeasurement;

pub fn apply_kinematic_constraints(
    state: &mut RtkState, 
    config: &EngineConfig, 
    imu_history: &[Vec<ImuMeasurement>]
) {
    if !config.enable_nhc {
        return;
    }

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
        let _ = crate::nhc::apply_zupt(state, zupt_var);
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
        let _ = crate::nhc::apply_nhc(state, 0.1, 0.1, &config.imu_mounting_angles, &config.imu_to_nhc_lever_arm, &omega_b);
    }
}
