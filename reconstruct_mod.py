import re

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    text = f.read()

# 1. Add pub mod for everything!
mods = """pub mod ambiguity;
pub mod kinematics;
pub mod matcher;
pub mod measurement;
pub mod ppp;
pub mod ppp_fg;
pub mod predictor;
pub mod processed_sat;
pub mod smoother;
pub mod ssr;
pub mod tcar;
pub mod tight_fg;
pub mod updater;
pub mod jacobian_verify;
"""

# Replace existing mods
text = re.sub(r"pub mod ambiguity;\n.*pub mod ppp;\n", mods, text, flags=re.DOTALL)

# 2. Add DynamicsModel.acceleration_psd
accel_psd = """    pub fn acceleration_psd(&self) -> f64 {
        match self {
            Self::Static => 0.001,
            Self::Pedestrian => 2.5,
            Self::Marine => 2.5,
            Self::Automotive => 0.1,
            Self::Airborne => 100.0,
        }
    }
"""
if "pub fn acceleration_psd" not in text:
    text = text.replace("pub enum DynamicsModel {", accel_psd + "\npub enum DynamicsModel {")
    # Wait, it belongs in `impl DynamicsModel`
    text = text.replace(accel_psd + "\npub enum DynamicsModel {", "pub enum DynamicsModel {") # revert
    impl_dyn = """
impl DynamicsModel {
""" + accel_psd + """
}
"""
    text = text.replace("pub enum DynamicsModel {\n    Static,\n    Pedestrian,\n    Automotive,\n    Marine,\n    Airborne,\n}\n", "pub enum DynamicsModel {\n    Static,\n    Pedestrian,\n    Automotive,\n    Marine,\n    Airborne,\n}\n" + impl_dyn)


# 3. Add EngineConfig fields
if "enable_doppler: bool" not in text:
    text = text.replace("pub max_ambiguity_age_epochs: u32,", "pub max_ambiguity_age_epochs: u32,\n    pub enable_doppler: bool,\n    pub enable_tropo: bool,")
    text = text.replace("max_ambiguity_age_epochs: 10,", "max_ambiguity_age_epochs: 10,\n            enable_doppler: true,\n            enable_tropo: true,")

# 4. Add ProcessingEngine fields
if "pub klobuchar:" not in text:
    text = text.replace("pub config: EngineConfig,", "pub config: EngineConfig,\n    pub klobuchar: Option<gneiss_core::atmosphere::KlobucharParams>,\n    pub precise_orbits: Option<crate::engine::ssr::PreciseOrbits>,\n    pub precise_clocks: Option<crate::engine::ssr::PreciseClocks>,\n    pub antex: Option<gneiss_parsers::antex::Antex>,")
    text = text.replace("config,", "config,\n            klobuchar: None,\n            precise_orbits: None,\n            precise_clocks: None,\n            antex: None,")

# 5. REMOVE process_rtk, update_ekf, predict_ekf, etc. from mod.rs since they are in updater.rs, etc.
# Actually, I'll just delete them if they exist.
# process_rtk
text = re.sub(r"    pub fn process_rtk\(.*?Ok\(self\.current_state\.as_ref\(\)\.unwrap\(\)\)\n    }", "", text, flags=re.DOTALL)
# process_spp
# text = re.sub(r"    pub fn process_spp\(.*?Ok\(self\.current_state\.as_ref\(\)\.unwrap\(\)\)\n    }", "", text, flags=re.DOTALL)
# Actually, process_spp is still in mod.rs!

# Let's save and see!
with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.write(text)

print("Reconstructed mod.rs!")
