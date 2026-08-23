//! Sliding-Window Factor Graph (SWFG) — the next-generation Gneiss estimator.
//!
//! Architecture:
//!   - Variable-based state representation (no fixed index vector)
//!   - Factor graph with LM optimization (not single-epoch EKF)
//!   - Schur complement marginalization for old epochs
//!   - Type-state configuration (invalid configs don't compile)
//!   - Unified measurement pipeline (SPP/PPP/RTK differ only in corrections)
//!
//! Sprints:
//!   1. Core graph + type-state definitions ✅ (this module)
//!   2. IMU preintegration (Forster method)
//!   3. Unified measurement pipeline
//!   4. AR-injected graph (LAMBDA)
//!   5. Schur complement marginalization
//!   6. Benchmark harness + CI

pub mod ar_integration;
pub mod benchmark;
pub mod config;
pub mod engine;
pub mod factor;
pub mod graph;
pub mod imu_preintegration;
pub mod kalman_smoother;
pub mod marginalization;
pub mod pipeline;
pub mod smoothing;
pub mod solver;
pub mod variables;

pub use config::EngineConfig;
pub use factor::Factor;
pub use graph::{EpochMetadata, EstimationGraph, MarginalPriorFactor};
pub use solver::SlidingWindowSolver;
pub use variables::{VariableDim, VariableId, VariableKind, VariableNode, VariableValues};
