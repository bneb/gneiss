//! Comprehensive Tier 1–4 Frame Safety & Epoch Alignment E2E Test Suite.
#![allow(clippy::unwrap_used, dead_code, unused_imports, clippy::wrong_self_convention, clippy::unnecessary_literal_unwrap)]

#[path = "test_frame_safety_e2e/types.rs"]
mod types;

#[path = "test_frame_safety_e2e/common.rs"]
mod common;

#[path = "test_frame_safety_e2e/tier1_spatial.rs"]
mod tier1_spatial;

#[path = "test_frame_safety_e2e/tier1_temporal.rs"]
mod tier1_temporal;

#[path = "test_frame_safety_e2e/tier1_geometry.rs"]
mod tier1_geometry;

#[path = "test_frame_safety_e2e/tier1_estimator.rs"]
mod tier1_estimator;

#[path = "test_frame_safety_e2e/tier2_boundaries.rs"]
mod tier2_boundaries;

#[path = "test_frame_safety_e2e/tier3_pairwise.rs"]
mod tier3_pairwise;

#[path = "test_frame_safety_e2e/tier4_scenarios.rs"]
mod tier4_scenarios;
