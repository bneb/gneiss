# BRIEFING — 2026-06-21T04:52:00Z

## Mission
Investigate phase wind-up sign bug in gneiss-rtk and propose fix and regression test.

## 🔒 My Identity
- Archetype: explorer
- Roles: Teamwork explorer (Read-only investigation)
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_windup_3
- Original parent: 404b157a-ed72-4236-bb0e-5434d4998c96
- Milestone: Milestone 1 Tier 1 Bugs

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Maintain BRIEFING.md and progress.md
- Produce analysis.md and handoff.md

## Current Parent
- Conversation ID: 404b157a-ed72-4236-bb0e-5434d4998c96
- Updated: 2026-06-21T04:52:00Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/engine/ppp.rs`
  - `crates/gneiss-rtk/src/engine/measurement.rs`
  - `crates/gneiss-rtk/src/engine/ppp_iekf.rs`
  - `crates/gneiss-core/src/windup.rs`
- **Key findings**:
  - Phase wind-up correction is added (`+ wup` / `+= windup`) to carrier phase observations in both RTK and PPP engines instead of being subtracted.
  - The correct physical model is `corrected_cp = raw_cp - windup`.
- **Unexplored areas**: None

## Key Decisions Made
- Analyzed codebase and identified all instances of sign discrepancy.
- Formulated fix strategy and designed unit regression tests.
- Documented findings in `analysis.md` and `handoff.md`.

## Artifact Index
- ORIGINAL_REQUEST.md — Archive of the task request
- BRIEFING.md — Persistent memory index
- progress.md — Liveness heartbeat and progress tracker
- analysis.md — Main analysis report on the bug
- handoff.md — 5-component handoff report
