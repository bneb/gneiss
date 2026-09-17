# BRIEFING — 2026-09-12T16:49:00Z

## Mission
Implement Frontier R2: Integer PPP-AR Engine via SINEX OSB, PCO/PCV, multi-constellation, single-differenced LAMBDA AR, and eval_ppp benchmark.

## 🔒 My Identity
- Archetype: teamwork_preview_worker
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_m2
- Original parent: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Milestone: Milestone M2 (Frontier R2)

## 🔒 Key Constraints
- Exclusive write ownership:
  - `crates/gneiss-parsers/src/sinex_bia.rs`
  - `crates/gneiss-parsers/src/antex.rs`
  - `crates/gneiss-parsers/src/receiver_pcv/`
  - `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
  - `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`
  - `crates/gneiss-rtk/src/swfg/engine/epoch.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/satpos.rs`
  - `crates/gneiss-rtk/src/bin/eval_ppp.rs`
- Do NOT edit any files outside this set.
- All code strictly adheres to:
  - File size < 500 LOC
  - Function size < 32 LOC
  - Nesting depth < 3 levels
  - 0 unwrap() in production code
  - 0 clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)
  - Passing test coverage with unit tests in `#[cfg(test)] mod tests`

## Current Parent
- Conversation ID: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Updated: not yet

## Task Summary
- **What to build**: Integer PPP-AR Engine via SINEX OSB, PCO/PCV, multi-constellation, single-differenced LAMBDA AR, eval_ppp benchmark.
- **Success criteria**:
  - Fast indexed SINEX OSB parser with WL/NL bias helpers.
  - Satellite nadir PCV and standalone receiver PCO/PCV for PPP.
  - Multi-constellation (GPS, Galileo, BeiDou, QZSS) support.
  - Single-differenced LAMBDA AR (MW wide-lane rounding + narrow-lane LAMBDA) wired into `ar_handler.rs`.
  - Ingest `com21374.bia` into `eval_ppp.rs`, verify integer ambiguity resolution on F9P drive closing discrepancy vs CSRS-PPP to sub-meter kinematic accuracy.
  - All workspace tests and clippy pass with 0 warnings.
- **Interface contracts**: `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md`
- **Code layout**: `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md § Code Layout`

## Key Decisions Made
- Starting with codebase investigation and baseline test runs.

## Artifact Index
- `.agents/worker_m2/DISPATCH.md` — Task assignment and requirements
- `.agents/worker_m2/BRIEFING.md` — Situational awareness and working memory
- `.agents/worker_m2/progress.md` — Progress log and liveness heartbeat
- `.agents/worker_m2/handoff.md` — Final handoff report

## Change Tracker
- **Files modified**: None yet
- **Build status**: Pending baseline check
- **Pending issues**: None

## Quality Status
- **Build/test result**: Not run yet
- **Lint status**: Not run yet
- **Tests added/modified**: None yet

## Loaded Skills
- None
