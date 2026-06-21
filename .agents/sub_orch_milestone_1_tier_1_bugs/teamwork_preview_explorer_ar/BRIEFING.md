# BRIEFING — 2026-06-21T04:26:21Z

## Mission
Investigate Bug 9: Sequential AR Covariance Mismatch in `gneiss-rtk/src/engine/ppp_iekf.rs`.

## 🔒 My Identity
- Archetype: explorer
- Roles: read-only investigator, synthesis reporter
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_ar
- Original parent: f16afb25-c177-42fe-985d-6840e173046f
- Milestone: sub_orch_milestone_1_tier_1_bugs

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Network mode: CODE_ONLY (no external web access)
- Write only to own folder

## Current Parent
- Conversation ID: f16afb25-c177-42fe-985d-6840e173046f
- Updated: 2026-06-21T04:29:10Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/engine/ppp_iekf.rs`
  - `crates/gneiss-rtk/src/engine/ppp_ins_iekf.rs`
  - `crates/gneiss-rtk/src/engine/ppp_common.rs`
  - `crates/gneiss-rtk/src/math/inversion.rs`
- **Key findings**:
  - Located the `resolve_cascade_ar` method in `ppp_iekf.rs` which performs per-constellation sequential Ambiguity Resolution.
  - Found that the method mutates the RTK state `state` in-place inside the sequential loop via `apply_state_vector` so that subsequent constellations can query the updated `state.covariance`.
  - Discovered that if a later constellation fails or if the final global position validation check fails, the method returns an `Err`, but the `state` has already been permanently modified in-place, causing a state/covariance contamination and mismatch with the float solution flags.
- **Unexplored areas**: None (problem is fully identified and scope is complete).

## Key Decisions Made
- Confirmed the cause of the Sequential AR Covariance Mismatch.
- Formulated a clean refactoring plan to pass the local covariance matrix `p_current` to `resolve_widelane_ar` as a parameter instead of having it read `state.covariance`.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_ar/handoff.md — Analysis and findings handoff report
