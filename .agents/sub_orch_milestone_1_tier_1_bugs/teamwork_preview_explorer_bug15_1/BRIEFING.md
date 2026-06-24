# BRIEFING — 2026-06-21T20:06:04Z

## Mission
Investigate Bug 15: Incorrect Broadcast Clock TGD Correction in gneiss-core.

## 🔒 My Identity
- Archetype: Explorer
- Roles: read-only explorer
- Working directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_1
- Original parent: 68da8f76-c613-4bbf-9560-65fc9a9f87a3
- Milestone: milestone_1_tier_1_bugs

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Code-only network mode (no external network requests)

## Current Parent
- Conversation ID: 68da8f76-c613-4bbf-9560-65fc9a9f87a3
- Updated: 2026-06-21T20:10:00Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-core/src/ephemeris.rs`: analyzed `BeidouEphemeris`, `position()`, and `position_iono_free()`
  - `crates/gneiss-parsers/src/rinex.rs`: analyzed `build_beidou_ephemeris` parser mapping
  - `crates/gneiss-rtk/src/engine/ppp.rs`: checked broadcast clock usage in PPP path
  - `crates/gneiss-rtk/src/estimators/spp.rs`: checked broadcast clock usage in SPP path
  - Checked other agent folders for peer reconciliation (`teamwork_preview_explorer_bug15_2`)
- **Key findings**:
  - Beidou legacy D1/D2 broadcast clock is referenced to B3I frequency.
  - Single frequency B1I uses `tgd1`, B2I uses `tgd2`.
  - Dual frequency iono-free combination needs to apply the combined timing group delay $T_{GD\_IF} = \frac{f_1^2 T_{GD1} - f_2^2 T_{GD2}}{f_1^2 - f_2^2}$.
  - Gneiss currently passes `0.0` for Beidou `position_iono_free` and does not parse `tgd2` correctly in RINEX parser (interpreting it as `aodc`).
- **Unexplored areas**: None. The bug is fully located and the fix strategy is designed.

## Key Decisions Made
- Reconciled findings with peer agent (`teamwork_preview_explorer_bug15_2`) to confirm RINEX layout mapping.
- Formulated the exact mathematical model and verification checks for the fix.

## Artifact Index
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_1/analysis.md — detailed analysis of TGD corrections, mathematical logic, and proposed fix strategy.
- /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/teamwork_preview_explorer_bug15_1/handoff.md — handoff report detailing findings, observations, logic chain, and verification steps.
