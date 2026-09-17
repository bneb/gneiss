# BRIEFING — 2026-09-13T02:01:00Z

## Mission
Survey Frontier R1: 15-State Error-State Kalman Filter (ESKF/MEKF) GNSS/INS and eval_odaiba_ins benchmark results.

## 🔒 My Identity
- Archetype: explorer
- Roles: investigation, synthesis
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1_status/
- Original parent: a6307386-3f81-4920-9a31-a6d124a2f8d6
- Milestone: Frontier R1 Survey

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Verify against AGENTS.md code standards (< 500 LOC, < 32 LOC/fn, nesting < 3, 0 unwrap in prod)
- All agent metadata in .agents/ only

## Current Parent
- Conversation ID: a6307386-3f81-4920-9a31-a6d124a2f8d6
- Updated: not yet

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/estimators/eskf/` (`types.rs`, `predict.rs`, `update.rs`, `constraints.rs`, `smoother.rs`, `mod.rs`)
  - `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`
  - `ORIGINAL_REQUEST.md`, `PROJECT.md`, `worker_r1/progress.md`
- **Key findings**:
  - ESKF unit tests: 18/18 passing (0.01s).
  - Code standards: all files < 500 LOC, functions < 32 LOC, nesting < 3, 0 unwrap, 0 clippy warnings.
  - Benchmark performance: Forward filter p50=3.131m, RMS=5.774m; RTS Smoothed p50=2.761m, RMS=5.472m.
  - Acceptance criteria (p50 < 2.5m, RMS < 5.2m): NOT YET MET (gap of 0.26m on p50, 0.27m on RMS).
  - Root cause: unmodeled urban multipath with overconfident R_pos (var_p = 0.001-0.04 m^2) under elevated rail/highway in Q2 and Q4; open-sky Q1 and Q3 easily beat targets (p50 ~ 1.3m, RMS ~ 2.2m).
- **Unexplored areas**: none (investigation complete).

## Key Decisions Made
- Executed benchmark and unit tests, extracted full quartile statistics and max error values.
- Compiled concrete tuning recommendations (innovation gating, Q tuning, lever arm, velocity update disabling) for Worker R1.

## Artifact Index
- DISPATCH.md — incoming dispatch instructions
- report.md — comprehensive survey and analysis report
- handoff.md — 5-component self-contained handoff
- progress.md — activity heartbeat
