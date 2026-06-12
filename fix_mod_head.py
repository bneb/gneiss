import os

# Reset mod.rs to head
os.system("git checkout crates/gneiss-rtk/src/engine/mod.rs")

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    text = f.read()

# 1. Add missing mods at the top
missing_mods = """pub mod processed_sat;
pub mod ssr;
pub mod tcar;
pub mod tight_fg;
pub mod ppp_fg;
pub mod smoother;
"""
text = text.replace("pub mod ppp;", "pub mod ppp;\n" + missing_mods)

# 2. Add DynamicsModel.acceleration_psd
accel_psd = """    pub fn acceleration_psd(&self) -> f64 {
        match self {
            Self::Static => 0.001,
            Self::Pedestrian => 2.5,
            Self::Marine => 2.5,
            Self::Automotive => 0.1,
            Self::Airborne => 100.0,
        }
    }"""
impl_dyn = """
impl DynamicsModel {
""" + accel_psd + """
}
"""
text = text.replace("pub enum DynamicsModel {\n    Static,\n    Pedestrian,\n    Automotive,\n    Marine,\n    Airborne,\n}\n", "pub enum DynamicsModel {\n    Static,\n    Pedestrian,\n    Automotive,\n    Marine,\n    Airborne,\n}\n" + impl_dyn)

# 3. Add EngineConfig fields
text = text.replace("pub max_ambiguity_age_epochs: u32,", "pub max_ambiguity_age_epochs: u32,\n    pub enable_doppler: bool,\n    pub enable_tropo: bool,")
text = text.replace("max_ambiguity_age_epochs: 10,", "max_ambiguity_age_epochs: 10,\n            enable_doppler: true,\n            enable_tropo: true,")

# 4. Add ProcessingEngine fields
text = text.replace("pub config: EngineConfig,", "pub config: EngineConfig,\n    pub klobuchar: Option<gneiss_core::atmosphere::KlobucharParams>,\n    pub precise_orbits: Option<crate::engine::ssr::PreciseOrbits>,\n    pub precise_clocks: Option<crate::engine::ssr::PreciseClocks>,\n    pub antex: Option<gneiss_parsers::antex::AntexDatabase>,")
text = text.replace("config,", "config,\n            klobuchar: None,\n            precise_orbits: None,\n            precise_clocks: None,\n            antex: None,")

with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.write(text)

print("Fixed mod.rs from HEAD.")
