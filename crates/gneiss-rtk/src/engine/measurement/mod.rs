pub(crate) mod windup;
pub(crate) mod types;
pub(crate) mod pseudorange;
pub(crate) mod carrier_phase;
pub(crate) mod doppler;
pub(crate) mod geometry;
pub(crate) mod components;
pub(crate) mod innovations;
pub(crate) mod builder;

// Re-export measurement_math items into the measurement namespace
pub use crate::engine::measurement_math::*;

// Re-export public types and functions for the flat API (backward compatibility)
pub use self::types::{
    DdContext, DdMeasurementContext, EkfUpdates, MeasurementEnvironment, SatState,
};
pub use self::pseudorange::compute_dd_pseudorange;
pub use self::carrier_phase::compute_dd_carrier_phase;
pub use self::doppler::compute_dd_doppler;
pub use self::geometry::EkfGeometryContext;
pub use self::innovations::compute_innovations;
pub use self::components::DdComponents;
pub use self::builder::{build_measurement_model, EkfMeasurementMatrices};

// Internal re-exports for cross-submodule visibility
pub(crate) use self::windup::update_windup_state_and_obs;
pub(crate) use self::geometry::find_ephemeris;
pub(crate) use self::components::{compute_dd_components, generate_measurement_updates};

#[cfg(test)]
mod tests_phase;
#[cfg(test)]
mod tests_measurements;
#[cfg(test)]
mod tests_integration;
