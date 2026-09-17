# BRIEFING — 2026-09-12T22:25:50Z

## Mission
Complete Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion, satellite nadir PCV & standalone receiver PCO/PCV, multi-constellation support (including QZSS), single-differenced LAMBDA ambiguity resolution, and sub-meter F9P kinematic PPP benchmark evaluation.

## 🔒 My Identity
- Archetype: implementer, qa, specialist
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_r2
- Original parent: 6231be0f-4267-4418-805f-47226c64a3af
- Milestone: M2 (Frontier R2: Integer PPP-AR Engine)

## 🔒 Key Constraints
- File size < 500 LOC
- Function size < 32 LOC
- Nesting depth < 3 levels
- Exactly 0 unwrap() in production code (only allowed in #[cfg(test)])
- 0 compiler warnings and clippy warnings (`cargo clippy -p gneiss-parsers -p gneiss-rtk --bin eval_ppp -- -D warnings`)
- All tests pass (`cargo test -p gneiss-parsers`, `cargo test -p gneiss-rtk`)
- Exclusive file ownership:
  - crates/gneiss-parsers/src/sinex_bia.rs
  - crates/gneiss-parsers/src/antex.rs
  - crates/gneiss-parsers/src/receiver_pcv/**
  - crates/gneiss-rtk/src/ambiguity/ppp_ar.rs
  - crates/gneiss-rtk/src/swfg/engine/ar_handler.rs
  - crates/gneiss-rtk/src/swfg/engine/epoch.rs
  - crates/gneiss-rtk/src/bin/eval_ppp.rs
- DO NOT CHEAT: All implementations genuine, no hardcoded results, no dummy facades.

## Current Parent
- Conversation ID: 6231be0f-4267-4418-805f-47226c64a3af
- Updated: 2026-09-12T22:25:50Z

## Task Summary
- **What to build**: Fast indexed SINEX OSB parser with wide-lane/narrow-lane bias helpers; satellite nadir PCV interpolation and standalone receiver PCO/PCV evaluation; multi-constellation support including QZSS 'J'; single-differenced LAMBDA integer ambiguity resolution wired into ar_handler; eval_ppp benchmark validation on F9P drive.
- **Success criteria**: Sub-meter kinematic accuracy on F9P drive vs CSRS-PPP; all unit tests pass; clippy clean; AGENTS.md rules satisfied.
- **Interface contracts**: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
- **Code layout**: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md § Code Layout

## Key Decisions Made
- Initializing workspace and reviewing existing implementations of sinex_bia.rs, antex.rs, receiver_pcv, ppp_ar.rs, ar_handler.rs, epoch.rs, and eval_ppp.rs.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/worker_r2/progress.md — liveness and milestone checklist
- /Users/kevin/projects/gneiss/.agents/worker_r2/handoff.md — handoff report

## Change Tracker
- **Files modified**: None yet
- **Build status**: Pending initial run
- **Pending issues**: None

## Quality Status
- **Build/test result**: Pending initial run
- **Lint status**: Pending initial run
- **Tests added/modified**: None yet

## Loaded Skills
- None
