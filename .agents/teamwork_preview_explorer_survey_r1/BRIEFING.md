# BRIEFING — 2026-09-12T16:46:20Z

## Mission
Survey codebase for Frontier R1 (15-State ESKF/MEKF GNSS/INS, RTS smoother, NHC/ZUPT, eval_odaiba_ins benchmark) and produce structured survey report.

## 🔒 My Identity
- Archetype: teamwork_preview_explorer
- Roles: explorer, synthesizer
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1
- Original parent: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Milestone: Frontier R1 Survey

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- AGENTS.md code standards: file size < 500 LOC, function size < 32 LOC, nesting < 3 levels, 0 unwrap() in prod, 0 warnings
- vel_att = f_e_skew * dt in predictor.rs must remain positive
- Sign conventions: predictor.rs:86-91 (attitude error state d_theta = -psi, left-multiplied global-frame)
- Phase windup: cp - windup (subtracted)

## Current Parent
- Conversation ID: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Updated: 2026-09-12T16:46:20Z

## Investigation State
- **Explored paths**: `eval_odaiba_ins.rs`, `swfg/imu_preintegration/smoother.rs`, `swfg/imu_preintegration/mod.rs`, `swfg/imu_preintegration/stationary.rs`, historical `predictor.rs` & `updater.rs`, `tests/src/inertial_outage_simulation.rs`, `scripts/check_network_benchmark.py`.
- **Key findings**:
  - Live benchmark baseline measured: GNSS raw p50=2.808m RMS=5.720m; 6-state RTS smoothed p50=2.907m RMS=5.508m. Target: p50<2.5m, RMS<5.2m.
  - Heading gap: 7.07° initial heading error is never corrected by the 6-state filter.
  - Sensor biases: ba=0 assumed; bg frozen after 350 samples.
  - Full mathematical formulation for 15-state ESKF, NHC attitude coupling, ZUPT, and 15-state RTS smoother completed in `survey_r1.md`.
- **Unexplored areas**: None for survey scope. Ready for implementation.

## Key Decisions Made
- Confirmed modular decomposition into `types.rs`, `predict.rs`, `update.rs`, `smoother.rs` to satisfy < 500 LOC constraint.
- Confirmed positive velocity-attitude coupling `vel_att = +f_e_skew * dt`.
- Output completed survey report `survey_r1.md`.

## Artifact Index
- `DISPATCH.md` — Task assignment and instructions
- `BRIEFING.md` — Persistent working memory
- `progress.md` — Progress tracker
- `survey_r1.md` — Complete Frontier R1 survey report
- `handoff.md` — 5-component handoff report
