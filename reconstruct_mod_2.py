import re

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    text = f.read()

# Make sure modules are added
missing_mods = """
pub mod processed_sat;
pub mod ssr;
pub mod tcar;
pub mod tight_fg;
pub mod ppp_fg;
pub mod predictor;
pub mod smoother;
pub mod updater;
pub mod jacobian_verify;
"""
if "pub mod processed_sat;" not in text:
    text = text.replace("pub mod ppp;", "pub mod ppp;" + missing_mods)

# Add enable_doppler and enable_tropo to EngineConfig struct definition
if "pub enable_doppler: bool" not in text:
    text = text.replace("pub max_ambiguity_age_epochs: u32,", "pub max_ambiguity_age_epochs: u32,\n    pub enable_doppler: bool,\n    pub enable_tropo: bool,")

# Add to EngineConfig::default()
if "enable_doppler: true" not in text:
    text = text.replace("max_ambiguity_age_epochs: 10,", "max_ambiguity_age_epochs: 10,\n            enable_doppler: true,\n            enable_tropo: true,")

# Remove antex from ProcessingEngine struct
text = text.replace("pub antex: Option<gneiss_parsers::antex::Antex>,", "")
text = text.replace("antex: None,", "")

with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.write(text)

print("Reconstructed mod 2")
