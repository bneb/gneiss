# Dispatch: Worker R4 — Frontier R4 (Unified Composite Integration)

Working Directory: /Users/kevin/projects/gneiss/.agents/worker_r4
Original Request: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master Plan: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md

## Objective
Implement Frontier R4: Unified Composite Integration.
1. Modular interfaces composing 15-state ESKF with Integer PPP-AR:
   - `crates/gneiss-rtk/src/composite/tc_ppp.rs` (Tightly-Coupled PPP/INS: coupling PPP carrier-phase/pseudorange or PPP solution with 15-state ESKF)
2. Modular interfaces composing 15-state ESKF with Network RTK VRS:
   - `crates/gneiss-rtk/src/composite/tc_rtk.rs` (Tightly-Coupled Network RTK/INS: coupling VRS synthesized observations or double-difference RTK with 15-state ESKF)
3. Module wiring and testing:
   - `crates/gneiss-rtk/src/composite/mod.rs`
   - Wire `pub mod composite;` into `crates/gneiss-rtk/src/lib.rs`
   - Comprehensive unit and integration tests verifying composite pipelines

## Exclusive File Ownership
- `crates/gneiss-rtk/src/composite/**`
- `crates/gneiss-rtk/src/lib.rs` (only adding `pub mod composite;`)

## Mandatory Rules & Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Follow AGENTS.md standards:
- < 500 LOC per file
- < 32 LOC per function
- < 3 nesting depth
- 0 unwrap() in production code
- 0 clippy warnings (`cargo clippy -p gneiss-rtk -- -D warnings`)
- All tests pass (`cargo test -p gneiss-rtk`)

## 2026-09-12T22:25:50Z
You are Worker R4 for the Gneiss positioning engine project.
Objective: Implement Frontier R4: Unified Composite Integration.
1. Modular interfaces composing 15-state ESKF with Integer PPP-AR:
   - `crates/gneiss-rtk/src/composite/tc_ppp.rs`: Tightly-Coupled PPP/INS pipeline, coupling un-differenced carrier-phase/pseudorange innovations or PPP positioning solution with the 15-state ESKF.
2. Modular interfaces composing 15-state ESKF with Network RTK VRS:
   - `crates/gneiss-rtk/src/composite/tc_rtk.rs`: Tightly-Coupled Network RTK/INS pipeline, coupling VRS synthesized observations or double-difference RTK innovations with the 15-state ESKF.
3. Module structure and exports:
   - `crates/gneiss-rtk/src/composite/mod.rs`
   - Wire `pub mod composite;` into `crates/gneiss-rtk/src/lib.rs`.
   - Comprehensive unit and integration tests verifying composite pipelines.

Exclusive File Ownership:
- crates/gneiss-rtk/src/composite/**
- crates/gneiss-rtk/src/lib.rs (only adding `pub mod composite;`)
