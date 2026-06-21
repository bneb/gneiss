# BRIEFING — 2026-06-21T10:01:43Z

## Mission
Analyze Bug 15: Incorrect Broadcast Clock TGD Correction.

## 🔒 My Identity
- Archetype: explorer
- Roles: Teamwork explorer
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_bug_15
- Original parent: df999d12-6411-4f6a-bf1a-36a92f258c2f
- Milestone: Bug 15 Investigation

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Analyze Bug 15: Incorrect Broadcast Clock TGD Correction.
- Propose clean correct fix strategy, recommend regression test design.

## Current Parent
- Conversation ID: df999d12-6411-4f6a-bf1a-36a92f258c2f
- Updated: 2026-06-21T10:01:43Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-core/src/ephemeris.rs`
  - `crates/gneiss-rtk/src/engine/ppp.rs`
  - `crates/gneiss-rtk/src/estimators/spp.rs`
  - `crates/gneiss-rtk/src/engine/measurement_math.rs`
  - `crates/gneiss-rtk/src/engine/processed_sat.rs`
- **Key findings**:
  - `calc_keplerian` in `ephemeris.rs` always subtracts `tgd`, meaning `position()` returns a single-frequency corrected clock.
  - PPP adds `tgd` back manually to get the iono-free clock.
  - SPP dual-frequency mode uses `position()` without adding `tgd` back, causing incorrect corrections.
- **Unexplored areas**: None.

## Key Decisions Made
- Proposed introducing `position_iono_free` to avoid the manual round-trip in PPP and fix SPP dual-frequency.
- Designed a regression test `test_broadcast_clock_tgd_correct` to verify the fix across GPS, Galileo, BeiDou, and QZSS.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_bug_15/ORIGINAL_REQUEST.md — Original request containing the user prompt
- /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_bug_15/analysis.md — Full analysis report with proposed code modifications and regression test design
- /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_bug_15/handoff.md — 5-component handoff report
