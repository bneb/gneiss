import sys

filename = "crates/gneiss-rtk/src/engine/mod.rs"
with open(filename, "r") as f:
    content = f.read()

# For process_spp
target1 = "if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, Some(&gneiss_core::atmosphere::KlobucharParams::default()), &crate::spp::SppConfig::default(), None) {"
replacement1 = """let prev_spp_state = self.current_state.as_ref().map(|s| {
            crate::spp::SppState {
                position: s.position,
                cdt: if crate::filter::CORE_STATE_SIZE > 15 { s.rcv_clk_bias } else { 0.0 },
                cdt_gal: 0.0,
                cdt_bds: 0.0,
                cdt_glo: 0.0,
            }
        });
        if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, Some(&gneiss_core::atmosphere::KlobucharParams::default()), &crate::spp::SppConfig::default(), prev_spp_state.as_ref()) {"""

if target1 in content:
    content = content.replace(target1, replacement1, 1)

# For process_rtk
target2 = "if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &self.ephemerides, Some(&gneiss_core::atmosphere::KlobucharParams::default()), &crate::spp::SppConfig::default(), None) {"
if target2 in content:
    content = content.replace(target2, replacement1, 1)

with open(filename, "w") as f:
    f.write(content)
