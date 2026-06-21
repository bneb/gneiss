# BRIEFING — 2026-06-21T04:38:21Z

## Mission
Investigate Bug 2: Velocity-Attitude Transition Sign Mismatch in `gneiss-rtk/src/engine/predictor.rs`.

## 🔒 My Identity
- Archetype: explorer
- Roles: teamwork_preview_explorer, explorer
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_vel_att
- Original parent: f16afb25-c177-42fe-985d-6840e173046f
- Milestone: milestone_1_tier_1_bugs

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- CODE_ONLY network mode: no external web or service access, no curl/wget/etc.

## Current Parent
- Conversation ID: f16afb25-c177-42fe-985d-6840e173046f
- Updated: 2026-06-21T04:39:30Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/engine/predictor.rs`
  - `crates/gneiss-rtk/src/estimators/ekf/filter.rs`
  - `crates/gneiss-rtk/src/engine/updater.rs`
  - `crates/gneiss-rtk/src/engine/tests_predictor.rs`
  - `crates/gneiss-rtk/src/tests_predictor.rs`
- **Key findings**:
  - Found the sign mismatch in `predictor.rs` line 86. The coupling block from attitude error to velocity error is computed as `let vel_att = f_e_skew * dt;` but mathematically must be `let vel_att = -f_e_skew * dt;`.
  - Traced the attitude representation and update step in `updater.rs` to show that attitude error is in the ECEF frame and applied via left multiplication.
  - Formulated a proof based on physical perturbations showing the sign is inverted.
- **Unexplored areas**: None (investigation complete).

## Key Decisions Made
- Scoped fix to the sign correction in `predictor.rs`.
- Designed a unit test that isolates the transition matrix velocity-attitude coupling.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_vel_att/ORIGINAL_REQUEST.md — Original mission request
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_vel_att/BRIEFING.md — Working briefing index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_vel_att/progress.md — Progress tracker
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_vel_att/handoff.md — Final investigation handoff report
