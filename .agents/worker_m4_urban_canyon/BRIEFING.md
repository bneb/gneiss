# BRIEFING — 2026-09-25T07:19:40Z

## Mission
Implement Milestone 4: C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR), DOP Geometry Guard, and Positive-Definite Ambiguity Conditioning (F10, F11, F12).

## 🔒 My Identity
- Archetype: teamwork_preview_worker
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon
- Original parent: 5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0
- Milestone: Milestone 4: C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR)

## 🔒 Key Constraints
- File size strictly < 500 LOC (CRITICAL: tc_ambiguity.rs 489 LOC, ar.rs 463 LOC, par.rs 284 LOC, ar_subsets.rs 285 LOC — must NOT exceed 500 LOC).
- Function size strictly <= 32 LOC.
- Nesting depth strictly < 3 levels.
- Zero `unwrap()` calls in production code (`?`, `match`, `if let`, `.expect("invariant")`).
- Zero clippy / compiler warnings under `-D warnings`.
- Exclusive write ownership:
  - `crates/gneiss-rtk/src/ambiguity/par.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`
  - `crates/gneiss-rtk/src/composite/tc_ambiguity.rs`
- Pass all tests: `cargo test -p gneiss-rtk --lib`, `cargo test --test test_urban_canyon_e2e`, CI smoke guards.

## Current Parent
- Conversation ID: 5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0
- Updated: 2026-09-25T07:11:23Z

## Task Summary
- **What to build**:
  - F10: Composite Quality Metric (CQM) candidate ranking in `ar_subsets.rs` and `par.rs` incorporating elevation, C/N0, lock time, CMC variance, and float variance.
  - F11: DOP guard and minimum subset size $\ge 4$ in `ar_subsets.rs` / `tc_ambiguity.rs` / `ar.rs`.
  - F12: Positive-definite ambiguity conditioning: zero cross-covariances and verify $\lambda_{\min}(P) \ge 10^{-6}$ in `ar.rs:condition_state_on_integers`, and preserve $Q_{aa} \succeq 0$ in `tc_ambiguity.rs`.
- **Success criteria**: Zero regressions, all tests pass, both CI smoke guards pass, all code quality standards met.
- **Interface contracts**: PROJECT.md Section 5.
- **Code layout**: `crates/gneiss-rtk/`

## Key Decisions Made
- `AmbiguityMetadata` and `compute_cqm` placed in `par.rs` with backward-compatible fallback when `metadata` is `None`.
- `ar_subsets.rs` re-exports `AmbiguityMetadata` and provides `select_par_candidates_with_dd`, deriving elevation from satellite positions, CMC noise from `pr_var_m2`, and SNR from carrier phase variance.
- DOP guard: `validate_subset_geometry` enforces minimum 4 satellites and PDOP <= 10.0 and HDOP <= 10.0 using `gneiss_core::dop::compute_dop_from_positions`.
- Covariance leakage eliminated in `ar.rs:condition_state_on_integers`: row and column cross-covariances are zeroed out for all fixed ambiguities, with $\lambda_{\min}(P) \ge 10^{-6}$ check.
- `tc_ambiguity.rs` symmetrizes $Q_{aa}$ and floors diagonals at $10^{-6}$ on integer fixing to guarantee $Q_{aa} \succeq 0$.

## Artifact Index
- `.agents/worker_m4_urban_canyon/DISPATCH.md` — Assignment instructions
- `.agents/worker_m4_urban_canyon/BRIEFING.md` — Situational awareness
- `.agents/worker_m4_urban_canyon/progress.md` — Liveness heartbeat
- `.agents/worker_m4_urban_canyon/handoff.md` — 5-component handoff report

## Change Tracker
- **Files modified**:
  - `crates/gneiss-rtk/src/ambiguity/par.rs` (377 LOC): Added `AmbiguityMetadata`, `compute_cqm`, `select_ils_subset_with_metadata`, unit tests.
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs` (423 LOC): Added geometry validation `validate_subset_geometry`, `validate_dd_subset_geometry`, `compute_metadata_from_dd`, `select_par_candidates_with_dd`, unit tests.
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs` (484 LOC): Integrated `select_par_candidates_with_dd`, geometry guard in `eval_par_subset`, zeroed cross-covariances and $\lambda_{\min} \ge 10^{-6}$ check in `condition_state_on_integers`.
  - `crates/gneiss-rtk/src/composite/tc_ambiguity.rs` (492 LOC): Added symmetrization and diagonal floor on post-fix $Q_{aa}$.
- **Build status**: Pass (all tests pass, 0 clippy warnings)
- **Pending issues**: None

## Quality Status
- **Build/test result**: Pass (449/449 lib tests, 51/51 urban canyon e2e tests, both CI smoke guards passed)
- **Lint status**: 0 compiler warnings, 0 clippy warnings
- **Tests added/modified**:
  - `test_cqm_ranking_clean_vs_corrupted` (in `par.rs`)
  - `test_select_ils_subset_with_metadata_prioritizes_high_cqm` (in `par.rs`)
  - `test_cqm_candidate_prioritization_over_coincidental_float` (in `ar_subsets.rs`)
  - `test_validate_subset_geometry_rejects_sub_four_and_collinear` (in `ar_subsets.rs`)
  - Enhanced `test_condition_state_on_integers_updates_variance_and_position` (in `ar.rs`) with cross-covariance zeroing and $\lambda_{\min} \ge 10^{-6}$ checks.

## Loaded Skills
- None
