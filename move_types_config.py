import re

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    mod_content = f.read()

# 1. Extract types
types_rs = """use serde::{Deserialize, Serialize};

/// Defines the core positioning mode
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
pub enum EngineMode {
    Spp,
    SppIns,
    SppInsLooselyCoupled,
    #[default]
    Rtk,
    RtkIns,
    RtkInsLooselyCoupled,
    Ppp,
    PppIns,
    PppInsLooselyCoupled,
}

impl EngineMode {
    pub fn is_tightly_coupled(&self) -> bool {
        matches!(self, Self::SppIns | Self::RtkIns | Self::PppIns)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum DynamicsModel {
    Static,
    Pedestrian,
    Automotive,
    Marine,
    Airborne,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EngineError {
    NoObservations,
    InitialSppFailed,
    StateDisappeared,
    InsufficientSatellites,
    MissingBasePosition,
    GeodeticMismatch(&'static str),
}

impl std::fmt::Display for EngineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EngineError::NoObservations => write!(f, "No observations available"),
            EngineError::InitialSppFailed => write!(f, "Initial SPP failed"),
            EngineError::StateDisappeared => write!(f, "EKF state disappeared mid-execution"),
            EngineError::InsufficientSatellites => write!(f, "Insufficient satellites for EKF update"),
            EngineError::MissingBasePosition => write!(f, "Base station position must be provided"),
            EngineError::GeodeticMismatch(msg) => write!(f, "Geodetic Gatekeeper Failed: {}", msg),
        }
    }
}
impl std::error::Error for EngineError {}
"""

with open("crates/gneiss-rtk/src/engine/types.rs", "w") as f:
    f.write(types_rs)

# Extract EngineConfig block and remove from mod.rs
engine_config_pattern = re.compile(r'/// Configuration for the RTK processing engine\.\n#\[derive\(Debug, Clone, serde::Deserialize, serde::Serialize\)\]\n#\[serde\(default\)\]\npub struct EngineConfig \{.*?\n\}\n\nimpl Default for EngineConfig \{.*?\n\}\n', re.DOTALL)
config_match = engine_config_pattern.search(mod_content)
if config_match:
    config_str = config_match.group(0)
    mod_content = mod_content.replace(config_str, "")
    with open("crates/gneiss-rtk/src/engine/config.rs", "a") as f:
        f.write("\n" + config_str)
else:
    print("Failed to match EngineConfig")

# Remove EngineMode, DynamicsModel, EngineError blocks from mod.rs
mod_content = re.sub(r'/// Defines the core positioning mode\n#\[derive\(Default, Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize, serde::Serialize\)\]\npub enum EngineMode \{.*?\n\}\n\nimpl EngineMode \{.*?\n\}\n', '', mod_content, flags=re.DOTALL)
mod_content = re.sub(r'#\[derive\(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize\)\]\n#\[serde\(rename_all = "lowercase"\)\]\npub enum DynamicsModel \{.*?\n\}\n', '', mod_content, flags=re.DOTALL)
mod_content = re.sub(r'#\[derive\(Debug, Clone, PartialEq\)\]\npub enum EngineError \{.*?\n\}\n\nimpl std::fmt::Display for EngineError \{.*?\n\}\nimpl std::error::Error for EngineError \{\}\n', '', mod_content, flags=re.DOTALL)

with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.write(mod_content)

