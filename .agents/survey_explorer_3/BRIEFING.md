# BRIEFING — 2026-09-24T13:39:55Z

## Mission
Investigate and map the codebase baseline and implementation points for R3 (Doppler cycle slip validation) and R4 (C/N0- and elevation-prioritized PAR).

## 🔒 My Identity
- Archetype: explorer
- Roles: survey, read-only investigation, architectural analysis
- Working directory: /Users/kevin/projects/gneiss/.agents/survey_explorer_3
- Original parent: c1309e2d-6c95-4b14-a86d-d26a13f2a150
- Milestone: Urban Canyon Fix Rate Expansion and Multipath Mitigation

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- Obey AGENTS.md standards (file < 500 LOC, fn <= 32 LOC, nesting < 3, no unwrap in production, frame safety)
- Write only to /Users/kevin/projects/gneiss/.agents/survey_explorer_3/
- Deliver survey_r3_r4.md and handoff.md; notify parent via send_message

## Current Parent
- Conversation ID: c1309e2d-6c95-4b14-a86d-d26a13f2a150
- Updated: not yet

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/post_process/screening.rs` (Doppler cycle slip detection, threshold, bands)
  - `crates/gneiss-rtk/src/estimators/doppler.rs` (Doppler velocity estimation, range rate observables)
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/mod.rs` (epoch processing, slip propagation, AR candidate resolution)
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs` (DD formation, slip pair checking, clock datum, pair_epochs)
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs` (Melbourne-Wübbena tracker, LLI check leak)
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_gate.rs` (seed ambiguity variance)
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs` (phase innovation gate)
  - `crates/gneiss-rtk/src/ambiguity/par.rs` (select_ils_subset, LDLT diagonal, greedy failure)
  - `crates/gneiss-rtk/src/composite/tc_ambiguity.rs` (solve_integers_full_or_par, covariance zeroing)
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs` (candidate selection, omission, partitioning)
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs` (try_full_ar, try_partial_ar, condition_state_on_integers off-diagonal leak)
  - `crates/gneiss-core/src/dop.rs` (DOP computation functions)
  - `crates/gneiss-rtk/src/ambiguity/lambda/mod.rs` (spectral regularization)
- **Key findings**:
  - Doppler slip detection threshold in `screening.rs` is too coarse (`(1.0 * dt).max(1.0)` cyc), missing 0.5-cycle slips. Bands 5 & 6 omitted.
  - Slipped arcs in `formation.rs` leak into `mw.rs` and `pw_tracker` without resetting.
  - `pair_epochs` tracking duration never resets upon slip.
  - PAR candidate sorting in `ar_subsets.rs` only considers float fraction and variance, ignoring elevation, C/N0, lock time, and CMC variance.
  - `condition_state_on_integers` in `ar.rs` suffers from off-diagonal covariance leakage, making $P$ indefinite.
- **Unexplored areas**: None for survey scope.

## Key Decisions Made
- Completed detailed architectural survey for R3 and R4 in `survey_r3_r4.md`.
- Formulated composite quality metric (CQM) for PAR candidate prioritization.
- Specified DOP guard integration and strict covariance positive-definiteness enforcement.

## Artifact Index
- DISPATCH.md — Dispatch instructions and tasks
- BRIEFING.md — Working memory and investigation state
- progress.md — Liveness heartbeat
- survey_r3_r4.md — Comprehensive survey report for R3 & R4
- handoff.md — 5-component handoff report
