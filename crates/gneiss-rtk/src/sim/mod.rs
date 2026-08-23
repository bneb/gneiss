//! High-fidelity physical simulation framework for GNSS RTK and post-processing evaluation.

pub mod generator;

pub use generator::{
    generate_simulation_dataset, generate_synthetic_ephemerides, SimulationConfig,
    SimulationDataset, TrajectoryProfile,
};
