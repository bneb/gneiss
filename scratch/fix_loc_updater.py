import re

with open("crates/gneiss-rtk/src/engine/updater.rs", "r") as f:
    content = f.read()

replacement = """pub fn apply_state_correction(state: &mut RtkState, dx: &DVector<f64>) {
    state.position.vector.x += dx[0];
    state.position.vector.y += dx[1];
    state.position.vector.z += dx[2];
    state.velocity.x += dx[3];
    state.velocity.y += dx[4];
    state.velocity.z += dx[5];
    if dx.len() >= crate::filter::CORE_STATE_SIZE { apply_imu_and_clock_correction(state, dx); }
    if dx.len() > crate::filter::CORE_STATE_SIZE {
        for i in 0..state.ambiguities.len() { state.ambiguities[i] += dx[crate::filter::CORE_STATE_SIZE + i]; }
    }
}

fn apply_imu_and_clock_correction(state: &mut RtkState, dx: &DVector<f64>) {
    let d_theta = Vector3::new(dx[6], dx[7], dx[8]);
    if d_theta.norm() > 1e-10 {
        let dq = UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_normalize(d_theta), d_theta.norm());
        state.attitude = state.attitude * dq;
        state.attitude.renormalize();
    }
    state.accel_bias.x += dx[9]; state.accel_bias.y += dx[10]; state.accel_bias.z += dx[11];
    state.gyro_bias.x += dx[12]; state.gyro_bias.y += dx[13]; state.gyro_bias.z += dx[14];
    if crate::filter::CORE_STATE_SIZE > 15 {
        state.rcv_clk_bias += dx[15]; state.rcv_clk_drift += dx[16];
        state.zwd = (state.zwd + dx[17]).max(0.0);
    }
}"""

start_idx = content.find("pub fn apply_state_correction(state: &mut RtkState, dx: &DVector<f64>) {")
if start_idx == -1:
    print("Could not find apply_state_correction")
    import sys; sys.exit(1)

end_idx = content.find("pub fn apply_joseph_covariance_update(", start_idx)

new_content = content[:start_idx] + replacement + "\n\n" + content[end_idx:]

with open("crates/gneiss-rtk/src/engine/updater.rs", "w") as f:
    f.write(new_content)

print("Fixed LOC for apply_state_correction in updater.rs")
