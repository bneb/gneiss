# BRIEFING — 2026-09-12T16:47:00Z

## Mission
Survey the Gneiss codebase for Frontier R3 (Network RTK VRS Atmospheric Engine, multi-CORS adjustment, Delaunay iono/tropo models, VRS synthesis, eval_network_ppk benchmark) and Frontier R4 (Unified Composite Integration: Tightly-Coupled PPP/INS and Network RTK/INS), as well as regression guard scripts (check_network_benchmark.py and check_multignss_benchmark.py), and deliver survey_r3_r4.md and handoff.md.

## 🔒 My Identity
- Archetype: explorer
- Roles: survey, investigation, synthesis
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/
- Original parent: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Milestone: survey_r3_r4

## 🔒 Key Constraints
- Read-only investigation — do NOT implement production code
- Adhere strictly to AGENTS.md constraints: file size < 500 LOC, function size < 32 LOC, nesting depth < 3 levels, 0 warnings, 0 unwraps
- Frame safety: relational coupling structurally enforced
- Deliver findings to survey_r3_r4.md and handoff.md, then notify parent orchestrator via send_message

## Current Parent
- Conversation ID: 1bd6ce81-03bf-4c40-b8b1-3b137333b5e7
- Updated: 2026-09-12T16:41:15Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/bin/eval_network_ppk.rs`, `eval_odaiba_ins.rs`, `eval_ppp.rs`
  - `crates/gneiss-rtk/src/post_process/vrs.rs`, `network.rs`, `forward.rs`, `mod.rs`
  - `crates/gneiss-rtk/src/swfg/imu_preintegration/smoother.rs`, `mod.rs`
  - `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
  - `crates/gneiss-parsers/src/sinex_bia.rs`
  - `crates/gneiss-core/src/atmosphere/ionosphere.rs`
  - `scripts/check_network_benchmark.py`, `scripts/check_multignss_benchmark.py`
  - Datasets: `cors_short_baseline`, `cors_sf_bay_network`, `multignss_2025d160`, `urbannav/tokyo/Tokyo_Data/Odaiba`
- **Key findings**:
  - `check_network_benchmark.py --smoke` and `check_multignss_benchmark.py --smoke` both verified passing cleanly (100% green).
  - Existing VRS in `vrs.rs` has planar regression (`fit_plane_gradient`) and observation synthesis skeleton, but lacks Delaunay triangulation, multi-station network adjustment, and integration with `eval_network_ppk.rs`.
  - Existing INS in `smoother.rs` is only 6-DOF ($p, v$) with open-loop dead-reckoned attitude and static gyro bias. Needs full 15-state ESKF ($\delta \mathbf{p}^e, \delta \mathbf{v}^e, \delta \boldsymbol{\theta}, \delta \mathbf{b}_a, \delta \mathbf{b}_g$), closed-loop quaternion reset, online bias estimation, and 15-state RTS smoother.
  - To respect AGENTS.md (< 500 LOC per file), new components must be modularly isolated (`delaunay.rs`, `network_adj.rs`, `estimators/eskf/`, `composite/tc_ppp.rs`, `composite/tc_rtk.rs`).
- **Unexplored areas**: None for survey scope.

## Key Decisions Made
- Executed and validated both regression guards live.
- Designed complete modular architectural plan in `survey_r3_r4.md` compliant with all AGENTS.md constraints.
- Prepared 5-component handoff report.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/survey_r3_r4.md — Comprehensive technical survey report
- /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/handoff.md — 5-component handoff report
