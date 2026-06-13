import re

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    content = f.read()

old_block = """        let (sat_pos, sat_vel, sat_clk, sat_drift) = eph.position(rover_obs.time);
        let dist = (sat_pos - rcv_pos_ecef).norm();"""

new_block = """        let tau_pr = pr1 / LIGHT_SPEED;
        let t_tx_nom = gneiss_core::time::GpsTime::new(rover_obs.time.week, rover_obs.time.tow - tau_pr);
        let (_, _, dt_s, _) = eph.position(t_tx_nom);
        let t_tx_true = gneiss_core::time::GpsTime::new(rover_obs.time.week, rover_obs.time.tow - tau_pr - dt_s);
        let (raw_vec, raw_vel, sat_clk, sat_drift) = eph.position(t_tx_true);
        
        let mut sat_pos = raw_vec;
        let mut sat_vel = raw_vel;
        for _ in 0..2 {
            let geometric_range = (sat_pos - rcv_pos_ecef).norm();
            let true_tau = geometric_range / LIGHT_SPEED;
            let theta = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * true_tau;
            let cos_t = libm::cos(theta);
            let sin_t = libm::sin(theta);
            sat_pos = nalgebra::Vector3::new(
                raw_vec.x * cos_t + raw_vec.y * sin_t,
                -raw_vec.x * sin_t + raw_vec.y * cos_t,
                raw_vec.z
            );
            sat_vel = nalgebra::Vector3::new(
                raw_vel.x * cos_t + raw_vel.y * sin_t,
                -raw_vel.x * sin_t + raw_vel.y * cos_t,
                raw_vel.z
            );
        }
        
        let dist = (sat_pos - rcv_pos_ecef).norm();"""

if old_block in content:
    content = content.replace(old_block, new_block)
else:
    print("Could not find old block in ppp.rs")

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(content)
