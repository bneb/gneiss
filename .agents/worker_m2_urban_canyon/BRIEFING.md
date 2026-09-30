# BRIEFING — 2026-09-24T19:12:00Z

## Mission
Implement Milestone 2 (R2): Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting (F4, F5, F6, F7).

## 🔒 My Identity
- Archetype: teamwork_preview_worker
- Roles: implementer, qa
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon
- Original parent: 5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0
- Milestone: Milestone 2 (R2: CMC Multipath Detection & Down-Weighting)

## 🔒 Key Constraints
- Exclusive write ownership:
  - crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs
  - crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs
  - crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs
  - crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs
  - crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs (if needed)
- AGENTS.md rules:
  - File size strictly < 500 LOC
  - Function size strictly <= 32 LOC
  - Nesting depth strictly < 3 levels
  - Exactly 0 unwrap() in production code
  - 0 compiler / clippy warnings (cargo clippy --workspace --all-targets -- -D warnings)
  - All workspace tests pass
  - Both CI smoke guard scripts pass:
    python3 scripts/check_network_benchmark.py --smoke
    python3 scripts/check_multignss_benchmark.py --smoke
- No cheating, no hardcoded test values, no fake implementations.

## Current Parent
- Conversation ID: 5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0
- Updated: 2026-09-24T19:12:00Z

## Task Summary
- **What to build**:
  - F4: Decouple PR gross error screen in screen.rs so carrier phase is retained when code residual exceeds gross error threshold.
  - F5: CMC arc tracking across continuous carrier tracking arcs in formation.rs / robust.rs.
  - F6: Adaptive code de-weighting / variance inflation in formation.rs / robust.rs without inflating carrier phase variance.
  - F7: MW multipath shielding in mw.rs to prevent false cycle-slip resets during code multipath jumps.
- **Success criteria**:
  - Unit tests in screen.rs, robust.rs, mw.rs passing.
  - All 51 tests in test_urban_canyon_e2e.rs passing.
  - cargo clippy --workspace --all-targets -- -D warnings clean (0 warnings).
  - cargo test -p gneiss-rtk --lib passing (441 passed).
  - Both benchmark smoke scripts passing.
- **Interface contracts**: /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md
- **Code layout**: crates/gneiss-rtk/src/estimators/rtk_iekf/

## Key Decisions Made
- Embedded `CmcTracker` inside `WidelaneTracker` (`mw.rs`) to prevent `rtk_iekf/mod.rs` (482 LOC) from exceeding the 500 LOC limit.
- Set CMC multipath threshold at tau_mp = 2.5m with a 5-epoch warm-up running mean.
- Outlier-resistant baseline freezing: when |cmc - baseline| > 2.5m, baseline is frozen while multipath deviation is applied as code variance inflation (R_PP <- R_PP + sigma_mp^2) leaving carrier phase variance nominal.
- In `mw.rs`, wide-lane innovation jumps > 1.0 cycle do not reset the track if shielded by detected code multipath along a continuous carrier arc.
- Updated `tests/urban_canyon/tier1_features.rs:165` under parent authorization to assert that measurements with active carrier phase are retained with inflated code variance rather than removed.

## Artifact Index
- DISPATCH.md — Assignment instructions
- BRIEFING.md — Working memory
- progress.md — Liveness heartbeat
- handoff.md — Final self-contained handoff report

## Change Tracker
- **Files modified**:
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`: Decoupled PR gross error screening (187 LOC).
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`: Added CmcTrack, CmcTracker, and apply_cmc_downweighting (389 LOC).
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`: Added CMC tracking to WidelaneTracker, shielded wide-lane absorption, modularized UPD solver (485 LOC).
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`: Integrated CMC downweighting in build_single_dd_pair and active key retention (452 LOC).
  - `crates/gneiss-rtk/tests/urban_canyon/tier1_features.rs`: Aligned assertion in test_r2_screen_gross_error_preserves_carrier_phase.
- **Build status**: Pass (cargo clippy --workspace --all-targets -- -D warnings clean, cargo test -p gneiss-rtk --lib 441 passed, test_urban_canyon_e2e 51 passed)
- **Pending issues**: None

## Quality Status
- **Build/test result**: Pass (441 lib tests, 51 e2e tests, both CI benchmark smoke scripts pass)
- **Lint status**: 0 compiler warnings, 0 clippy warnings
- **Tests added/modified**:
  - `screen.rs`: `screen_preserves_carrier_phase_and_deweights_code`
  - `update/robust.rs`: `test_cmc_track_clean_arc_accumulates_without_multipath`, `test_cmc_track_detects_multipath_step_and_freezes_baseline`, `test_cmc_track_resets_on_cycle_slip`, `test_apply_cmc_downweighting_inflates_code_only`
  - `mw.rs`: `test_tracker_shields_against_multipath_jump`
  - `tier1_features.rs`: `test_r2_screen_gross_error_preserves_carrier_phase`

## Loaded Skills
None
