# BRIEFING — 2026-09-12T16:47:30Z

## Mission
Survey codebase for Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion.

## 🔒 My Identity
- Archetype: teamwork_preview_explorer
- Roles: explorer, analyst, investigator
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2
- Original parent: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Milestone: Frontier R2 Codebase Survey

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- File size < 500 LOC
- Function size < 32 LOC
- Nesting depth < 3 levels
- Zero compiler warnings, 0 unwrap() in prod
- Frame Safety conventions

## Current Parent
- Conversation ID: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Updated: 2026-09-12T16:47:30Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-parsers/src/sinex_bia.rs`: Bias-SINEX parser, records, and lookups
  - `crates/gneiss-parsers/src/antex.rs` & `receiver_pcv/`: Satellite and receiver antenna models
  - `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` & `lambda/`: Wide-lane rounding and narrow-lane LAMBDA
  - `crates/gneiss-rtk/src/swfg/engine/epoch.rs` & `builder.rs` & `uduc_builder.rs`: Observation extraction and factors
  - `crates/gneiss-rtk/src/bin/eval_ppp.rs`: PPP evaluation harness against CSRS-PPP and RTK truth
- **Key findings**:
  - `SinexBias` is implemented but needs $O(1)$ indexing and wide-lane/narrow-lane bias helpers.
  - `com21374.bia` exists in the F9P dataset but was set to `None` in `eval_ppp.rs`.
  - `PppArSolver` is implemented but orphaned from the runtime execution pipeline.
  - Carrier phase ambiguities in un-differenced PPP require between-satellite single-differencing to cancel receiver initial phase offsets.
  - BeiDou and QZSS are artificially excluded in `epoch.rs:364` (`is_supp`) and missing in `satpos.rs:93`.
  - Current float PPP achieves $10.11\text{ m}$ RMS on F9P; integer PPP-AR with OSB and LAMBDA will achieve sub-meter kinematic accuracy closing the gap to CSRS-PPP ($0.296\text{ m}$ RMS).
- **Unexplored areas**: None for Frontier R2 survey; ready for implementation dispatch.

## Key Decisions Made
- Authored comprehensive survey report: `survey_r2.md`.
- Authored 5-component handoff report: `handoff.md`.
- Documented exact file paths, line numbers, data structures, and formulas.

## Artifact Index
- DISPATCH.md — Assignment instructions
- survey_r2.md — Survey report
- handoff.md — 5-component handoff report
- progress.md — Liveness heartbeat
