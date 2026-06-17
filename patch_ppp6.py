import os

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    content = f.read()

content = content.replace("crate::engine::ppp_math::crate::engine::ppp_math::", "crate::engine::ppp_math::")
content = content.replace("gneiss_core::ephemeris::EphemerisData", "gneiss_core::ephemeris::Ephemeris")
content = content.replace("crate::engine::ppp_math::check_clock_diff(sat, dt_s, brdc_clk);", "")
content = content.replace("crate::engine::ppp_math::check_pos_diff(sat, sp3_p, brdc_pos);", "")
content = content.replace("let mut sat_pos = brdc_pos; let mut sat_vel = brdc_vel;", "let mut sat_pos: Vector3<f64> = brdc_pos; let mut sat_vel: Vector3<f64> = brdc_vel;")
content = content.replace("crate::engine::ppp_math::detect_cycle_slip(sat.sat_obs, prev);", "crate::engine::ppp_math::detect_cycle_slip(sat.sat_obs, prev as u32);")
content = content.replace("state.locktimes.insert((sat.sat_obs.sat, 1), new_lk);", "state.locktimes.insert((sat.sat_obs.sat, 1), new_lk as u16);")

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(content)
