# Task Assignment — Sub-Orchestrator Milestone M2: Integer PPP-AR Engine via SINEX OSB

## Role & Mission
You are a sub-orchestrator (`teamwork_preview_orchestrator`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/sub_orch_m2/`.
Your parent orchestrator is: `1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`.

## Authoritative Documents to Read Before Starting Work
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md` (project architecture, interface contracts, and file layout)
3. `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2/survey_r2.md` (comprehensive codebase analysis, formula derivations, and benchmark investigation)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code quality standards)

## Scope & Milestone Objectives
Implement Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion:
1. Fast indexed SINEX BIA/OSB phase and code bias ingestion: Update `crates/gneiss-parsers/src/sinex_bia.rs` with indexed lookup ($O(1)$) and helper methods to compute satellite wide-lane (Melbourne-Wübbena) and narrow-lane phase biases.
2. Satellite nadir-dependent PCV interpolation (`crates/gneiss-parsers/src/antex.rs` and `crates/gneiss-rtk/src/estimators/rtk_iekf/satpos.rs`).
3. Standalone receiver PCO and PCV modeling for un-differenced rover observations in PPP mode (`crates/gneiss-parsers/src/receiver_pcv/`).
4. Multi-constellation enablement across GPS, Galileo, BeiDou, and QZSS: Update `epoch.rs` (remove restrictive `is_supp` check), `satpos.rs` (include `'J'` QZSS mapping), and frequency selection.
5. Single-differenced Wide-Lane (MW) rounding and Narrow-Lane LAMBDA integer ambiguity resolution: Wire `PppArSolver` (`crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`) into `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs` to fix single-differenced integer ambiguities without physical base stations.
6. Benchmark validation: Ingest `datasets/rtkexplorer/sample_1/f9p_ppp_1224/com21374.bia` in `crates/gneiss-rtk/src/bin/eval_ppp.rs`. Resolve integer ambiguities on the F9P kinematic drive to achieve sub-meter kinematic accuracy closing the discrepancy vs CSRS-PPP ($0.296\text{ m}$ RMS).

## Exclusive Write Ownership
You and your dispatched workers own ONLY:
- `crates/gneiss-parsers/src/sinex_bia.rs`
- `crates/gneiss-parsers/src/antex.rs`
- `crates/gneiss-parsers/src/receiver_pcv/`
- `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
- `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`
- `crates/gneiss-rtk/src/swfg/engine/epoch.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/satpos.rs`
- `crates/gneiss-rtk/src/bin/eval_ppp.rs`
Do NOT write to files owned by other milestones.

## Orchestrator Procedure & Gate Verification
Execute via the standard iteration loop: Explorer -> Worker -> Reviewer -> Challenger -> Forensic Auditor -> Gate.
- Worker dispatch prompt MUST include the mandatory integrity warning verbatim.
- Forensic Auditor verdict is a non-negotiable binary veto.
- All code must satisfy AGENTS.md (< 500 LOC/file, < 32 LOC/function, < 3 nesting, 0 unwrap in prod, 0 warnings).
- `cargo test -p gneiss-rtk --lib ambiguity::ppp_ar` and `cargo test -p gneiss-parsers --lib sinex_bia` must pass cleanly.
- `cargo run --release --bin eval_ppp` must demonstrate integer ambiguity resolution and sub-meter kinematic accuracy.

When complete, write `handoff.md` and notify parent orchestrator (`1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`).
