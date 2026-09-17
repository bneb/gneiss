## 2026-09-13T02:06:28Z

You are Worker R2 tasked with completing the Integer PPP-AR engine and achieving sub-meter kinematic accuracy on the F9P benchmark (`eval_ppp`).

Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_r2_fix/

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Authoritative user request (MUST read first): /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master project plan: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Explorer R2 Report: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2_status/report.md
Explorer R2 Handoff: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2_status/handoff.md

FILE OWNERSHIP:
You have exclusive write access to:
- crates/gneiss-rtk/src/bin/eval_ppp.rs
- crates/gneiss-rtk/src/ambiguity/ppp_ar.rs
- crates/gneiss-rtk/src/swfg/engine/ar_handler.rs
- crates/gneiss-rtk/src/swfg/engine/epoch.rs
- crates/gneiss-parsers/src/sinex_bia.rs

GOALS & ACCEPTANCE CRITERIA:
1. Fix the F9P kinematic benchmark configuration in `crates/gneiss-rtk/src/bin/eval_ppp.rs:436`:
   - Set `is_kinematic: true` for the moving F9P drive dataset so that SWFG does not collapse the trajectory into a single static point.
2. Fix arc tracking in `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs:20`:
   - Pass the actual epoch index to `PppMwTracker::update` rather than hardcoded 0.
3. Fix ambiguity constraint factor injection in SWFG:
   - Inject single-difference constraint factors between satellite pairs rather than un-differenced prior factors that lock the receiver clock.
4. Ensure receiver antenna PCO/PCV corrections from `receiver_pcv` are applied in `epoch.rs` if needed for rover-only PPP.
5. Fix Clippy warning in `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs:478:27` (`useless_vec` -> array `[12, -8, 15, 7]`).
6. Refactor `ppp_ar.rs` and touched files to strictly satisfy AGENTS.md:
   - File size < 500 LOC (currently 506 LOC, trim below 500)
   - Function size < 32 LOC (break down functions > 32 LOC into focused helper functions)
   - Nesting depth < 3 levels
   - Zero unwrap() in production code
7. Run the benchmark:
   - `PPP_ONLY=f9p cargo run --release --bin eval_ppp`
   - Verify that integer ambiguities are resolved and kinematic accuracy vs CSRS-PPP achieves sub-meter kinematic accuracy.
8. Run unit and crate tests:
   - `cargo test -p gneiss-parsers --lib sinex_bia`
   - `cargo test -p gneiss-rtk --lib ambiguity::ppp_ar`
   - `cargo clippy -p gneiss-rtk -p gneiss-parsers --all-targets -- -D warnings`
9. Write a complete handoff report to:
   /Users/kevin/projects/gneiss/.agents/worker_r2_fix/handoff.md
10. Send your completion message back to the orchestrator.

## 2026-09-13T02:20:13Z
**Context**: Checking in on Frontier R2 Integer PPP-AR implementation.
**Content**: Please report your current progress on eval_ppp benchmark, arc tracking, and AGENTS.md compliance.
**Action**: Provide a brief progress update.

## 2026-09-13T02:40:06Z
**Context**: Liveness check on Frontier R2 PPP-AR Worker.
**Content**: It has been 20 minutes since your last progress update. Please report your current status, what step you are currently running (compiling, benchmarking, or debugging), and estimated time to completion.
**Action**: Send an immediate status update.
