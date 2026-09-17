# BRIEFING — 2026-09-13T02:03:00Z

## Mission
Investigate Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion and the eval_ppp benchmark, verify status, evaluate against CSRS-PPP / RTK truth, identify gaps, check AGENTS.md compliance, and synthesize findings.

## 🔒 My Identity
- Archetype: explorer
- Roles: investigation, synthesis
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2_status
- Original parent: a6307386-3f81-4920-9a31-a6d124a2f8d6
- Milestone: Frontier R2 Survey

## 🔒 Key Constraints
- Read-only investigation — do NOT implement changes in source code (only write to own agent folder)
- Files for content delivery, Messages for coordination
- Self-contained 5-component handoff report
- AGENTS.md compliance checks (<500 LOC/file, <32 LOC/fn, <3 nesting, 0 unwrap in prod)

## Current Parent
- Conversation ID: a6307386-3f81-4920-9a31-a6d124a2f8d6
- Updated: not yet

## Investigation State
- **Explored paths**:
  - `crates/gneiss-parsers/src/sinex_bia.rs`
  - `crates/gneiss-parsers/src/antex.rs`
  - `crates/gneiss-parsers/src/receiver_pcv/`
  - `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
  - `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`
  - `crates/gneiss-rtk/src/swfg/engine/epoch.rs`
  - `crates/gneiss-rtk/src/bin/eval_ppp.rs`
- **Key findings**:
  - Acceptance criterion (sub-meter kinematic PPP-AR vs CSRS-PPP 0.296m RMS) NOT MET.
  - Current 3D RMS vs CSRS-PPP: 3.358m; vs RTK truth: 3.048m; Fix rate: 0.0%.
  - `eval_ppp.rs:436` mistakenly set `is_kinematic: false`, collapsing 600 epochs of vehicle driving to a static point.
  - `ar_handler.rs:20` hardcoded epoch to `0`, breaking continuous MW tracking.
  - `inject_fixed_priors` injected un-differenced priors that freeze the receiver clock state.
  - AGENTS.md violations: `ppp_ar.rs` is 506 LOC (> 500 limit); multiple functions > 32 LOC.
- **Unexplored areas**: None for R2 survey.

## Key Decisions Made
- Executed unit tests and release benchmark `PPP_ONLY=f9p cargo run --release --bin eval_ppp`.
- Computed exact East, North, Up component errors via `compute_enu.py`.
- Formulated complete code modification blueprint for worker implementers in `report.md` and `handoff.md`.

## Artifact Index
- DISPATCH.md — record of initial dispatch
- BRIEFING.md — persistent working memory
- progress.md — liveness heartbeat
- compute_enu.py — exact ENU statistics evaluation script
- report.md — comprehensive survey report
- handoff.md — 5-component handoff report
