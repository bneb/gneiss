//! Urban Canyon Fix Rate Expansion and Multipath Mitigation E2E Test Suite.
//!
//! Comprehensive 4-tier opaque-box integration tests covering:
//! - R1: Adaptive C/N0 (SNR) and elevation observation covariance weighting
//! - R2: Code-Minus-Carrier (CMC) multipath detection & down-weighting
//! - R3: Doppler-assisted cycle slip detection & phase continuity validation
//! - R4: C/N0- and elevation-prioritized Partial Ambiguity Resolution (PAR)
#![allow(clippy::unwrap_used)]

mod urban_canyon;
