# BRIEFING — 2026-09-25T21:10:00Z

## Mission
Survey R3 (Estimator State & Pipeline Refactoring) & R4 (Relational Coupling & Frame Safety Invariants) & Benchmarks for Gneiss Frame Safety & Epoch Alignment.

## 🔒 My Identity
- Archetype: explorer
- Roles: investigation, synthesis
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3
- Original parent: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Milestone: Survey 3 (R3 & R4 Estimators, Pipelines, Relational Coupling, and Benchmarks)

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Strictly follow AGENTS.md rules (<500 LOC, <=32 LOC functions, <3 nesting depth, 0 unwrap, 0 warnings)
- Do not modify source code; only write to working directory
- Send all results and updates via send_message to parent (db66ae0c-b21b-4e14-ac97-93509c51c4b0)

## Current Parent
- Conversation ID: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Updated: 2026-09-25T21:04:10Z

## Investigation State
- **Explored paths**: `crates/gneiss-rtk/src/estimators/eskf/{types.rs, predict.rs, update.rs, dd_update.rs, constraints.rs, alignment.rs}`, `crates/gneiss-rtk/src/swfg/engine/{mod.rs, builder.rs}`, `crates/gneiss-rtk/src/post_process/mod.rs`, `crates/gneiss-rtk/src/bin/eval_odaiba_ins/{main.rs, odaiba_helpers.rs}`, `crates/gneiss-parsers/src/receiver_pcv/mod.rs`, `crates/gneiss-parsers/src/receiver_antenna/mod.rs`, `crates/gneiss-rtk/src/estimators/rtk_iekf/{formation.rs, update/system.rs}`.
- **Key findings**:
  1. `EskfState`, `SwfgEngine`, `PostProcessOptions`, and `eval_odaiba_ins` use bare `Vector3<f64>` across coordinate frames (ECEF vs Body FRD lever arm vs NED covariance).
  2. Latent coordinate frame bug at `receiver_pcv/mod.rs:175`: `az_el(rov_llh, rov_llh, sat_pos)` passed `rov_llh` as `pos_ecef`.
  3. Relational coupling failure at `receiver_antenna/mod.rs:277-305`: rover zenith and azimuth were borrowed directly for base antenna correction.
  4. Golden baseline verified: 789 tests pass, both CI smoke scripts output `ALL CHECKS PASSED`, and `eval_odaiba_ins` achieves $p_{50} = 1.751\text{ m} \le 1.80\text{ m}$ and $\text{RMS} = 3.479\text{ m} \le 3.50\text{ m}$ with 0 false fixes.
  5. File size limits: `main.rs` (475 LOC), `formation.rs` (457 LOC), and `dd_update.rs` (454 LOC) near 500 LOC ceiling.
- **Unexplored areas**: None for R3/R4 survey scope. Ready for implementation phase handoff.

## Key Decisions Made
- Designed canonical `DoubleDiffGeometry<F>` struct ensuring geometric disagreement between rover and base is unrepresentable.
- Documented file size mitigation strategy (isolating geometry math to dedicated submodules).

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3/DISPATCH.md — Dispatch log
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3/BRIEFING.md — Persistent memory
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3/progress.md — Liveness heartbeat & task progress
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3/survey_r3_r4_estimators.md — Comprehensive survey findings
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3/handoff.md — 5-component handoff report
