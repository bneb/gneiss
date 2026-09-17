# BRIEFING — 2026-09-13T02:06:28Z

## Mission
Complete Integer PPP-AR engine and achieve sub-meter kinematic accuracy on F9P benchmark (`eval_ppp`).

## 🔒 My Identity
- Archetype: worker
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_r2_fix
- Original parent: a6307386-3f81-4920-9a31-a6d124a2f8d6
- Milestone: Integer PPP-AR Engine and F9P Kinematic Accuracy

## 🔒 Key Constraints
- Exclusive write access to:
  - crates/gneiss-rtk/src/bin/eval_ppp.rs
  - crates/gneiss-rtk/src/ambiguity/ppp_ar.rs
  - crates/gneiss-rtk/src/swfg/engine/ar_handler.rs
  - crates/gneiss-rtk/src/swfg/engine/epoch.rs
  - crates/gneiss-parsers/src/sinex_bia.rs
- File size < 500 LOC
- Function size < 32 LOC
- Nesting depth < 3 levels
- Zero unwrap() in production code
- 0 compiler warnings, 0 clippy warnings
- Tests pass (>95% coverage, no unverified claims)
- Integrity mandate: genuine implementation, no cheating or hardcoding

## Current Parent
- Conversation ID: a6307386-3f81-4920-9a31-a6d124a2f8d6
- Updated: not yet

## Task Summary
- **What to build**:
  1. Fix F9P kinematic benchmark configuration (`is_kinematic: true` in eval_ppp.rs:436)
  2. Fix arc tracking in ar_handler.rs (pass actual epoch index to PppMwTracker::update)
  3. Fix ambiguity constraint factor injection in SWFG (single-difference constraint factors between satellite pairs)
  4. Ensure receiver antenna PCO/PCV corrections from receiver_pcv are applied in epoch.rs
  5. Fix clippy warning in ppp_ar.rs:478:27
  6. Refactor ppp_ar.rs and touched files to strictly satisfy AGENTS.md (<500 LOC, <32 LOC/fn, <3 nesting)
  7. Run eval_ppp benchmark, verify integer ambiguities resolved & sub-meter kinematic accuracy
  8. Run unit and crate tests & clippy
- **Success criteria**:
  - eval_ppp succeeds with sub-meter kinematic accuracy on F9P
  - Ambiguities resolved
  - All tests & clippy pass with 0 warnings
  - AGENTS.md standards strictly met
- **Interface contracts**: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
- **Code layout**: Gneiss workspace crates

## Key Decisions Made
- Initializing task setup.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/worker_r2_fix/DISPATCH.md
- /Users/kevin/projects/gneiss/.agents/worker_r2_fix/BRIEFING.md
- /Users/kevin/projects/gneiss/.agents/worker_r2_fix/progress.md

## Change Tracker
- **Files modified**: None yet
- **Build status**: Untested
- **Pending issues**: TBD

## Quality Status
- **Build/test result**: Untested
- **Lint status**: Untested
- **Tests added/modified**: None yet

## Loaded Skills
- None
