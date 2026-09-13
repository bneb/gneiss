//! Comprehensive Tier 1–4 Frontiers E2E Test Suite.
#![allow(clippy::unwrap_used)]

#[path = "test_frontiers_e2e/common.rs"]
mod common;

#[path = "test_frontiers_e2e/tier1_eskf.rs"]
mod tier1_eskf;

#[path = "test_frontiers_e2e/tier1_ppp.rs"]
mod tier1_ppp;

#[path = "test_frontiers_e2e/tier1_vrs.rs"]
mod tier1_vrs;

#[path = "test_frontiers_e2e/tier1_composite.rs"]
mod tier1_composite;

#[path = "test_frontiers_e2e/tier2_eskf.rs"]
mod tier2_eskf;

#[path = "test_frontiers_e2e/tier2_ppp.rs"]
mod tier2_ppp;

#[path = "test_frontiers_e2e/tier2_vrs.rs"]
mod tier2_vrs;

#[path = "test_frontiers_e2e/tier2_composite.rs"]
mod tier2_composite;

#[path = "test_frontiers_e2e/tier3_pairwise.rs"]
mod tier3_pairwise;

#[path = "test_frontiers_e2e/tier4_scenarios.rs"]
mod tier4_scenarios;
