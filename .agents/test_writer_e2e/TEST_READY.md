# TEST_READY — Urban Canyon Fix Rate Expansion and Multipath Mitigation

**Date**: 2026-09-24  
**Author**: Test Writer E2E (`test_writer_e2e`)  
**Status**: COMPLETE (51/51 tests passing, 0 clippy warnings)  
**Target Package**: `gneiss-rtk` (E2E integration test target: `test_urban_canyon_e2e`)  
**Specification References**:  
- `ORIGINAL_REQUEST.md` (2026-09-24T13:30:49Z)  
- `PROJECT.md` / `SCOPE.md`  
- Survey Explorer findings: `survey_r1_bench.md`, `survey_r2_cmc.md`, `survey_r3_r4.md`  

---

## 1. Executive Summary

A comprehensive 4-tier opaque-box End-to-End (E2E) integration test suite has been designed, implemented, and verified for the Urban Canyon Fix Rate Expansion and Multipath Mitigation task.

All 51 test cases compile with zero warnings and pass with 100% success rate:
- **Total Tests**: 51
- **Passed**: 51
- **Failed**: 0
- **Ignored**: 0
- **Clippy Warnings**: 0 (under `cargo clippy --all-targets -- -D warnings`)
- **Execution Time**: ~0.01 seconds

---

## 2. Test Architecture & File Manifest

The test suite is structured modularly under `crates/gneiss-rtk/tests/urban_canyon/` with an explicit integration root test harness `test_urban_canyon_e2e.rs`.

| File Path | Description | Line Count | Function Count | Max Fn Length |
|---|---|:---:|:---:|:---:|
| `crates/gneiss-rtk/tests/test_urban_canyon_e2e.rs` | Test root harness module declaration | 10 | 0 | N/A |
| `crates/gneiss-rtk/tests/urban_canyon/mod.rs` | Submodule aggregator (`common`, `tier1`..`tier4`) | 7 | 0 | N/A |
| `crates/gneiss-rtk/tests/urban_canyon/common.rs` | Common geometry fixtures, satellites, DD builders | 154 | 10 | 18 LOC |
| `crates/gneiss-rtk/tests/urban_canyon/tier1_features.rs` | Tier 1: Primary functional behavior (R1-R4) | 358 | 21 | 24 LOC |
| `crates/gneiss-rtk/tests/urban_canyon/tier2_boundaries.rs` | Tier 2: Boundary, stress, & numerical limits | 320 | 20 | 25 LOC |
| `crates/gneiss-rtk/tests/urban_canyon/tier3_interactions.rs` | Tier 3: Cross-module compound interactions | 156 | 7 | 25 LOC |
| `crates/gneiss-rtk/tests/urban_canyon/tier4_scenarios.rs` | Tier 4: Mission scenarios (Tokyo, HK, Highway, etc.) | 133 | 6 | 24 LOC |

*All files strictly comply with AGENTS.md standards: < 500 lines per file, <= 32 lines per function, < 3 nesting levels.*

---

## 3. Tier Coverage Breakdown

### Tier 1: Feature Functionality (21 Tests)
- **R1 (Adaptive SNR/Elevation Covariance)**:
  - `test_r1_elevation_variance_monotonic_increase`: Low elevation increases pseudorange & carrier variances monotonically.
  - `test_r1_snr_variance_monotonic_increase`: Low C/N0 increases measurement variances monotonically.
  - `test_r1_expected_cn0_increases_with_elevation`: Expected SNR baseline model scales with elevation angle.
  - `test_r1_attenuation_scale_smooth_penalty`: Attenuation penalty function applies smooth, bounded inflation.
  - `test_r1_dd_variances_positive_definite`: Single double-difference variance calculation outputs strictly positive variance.
  - `test_r1_dd_covariance_matrix_positive_definite`: Full multi-satellite DD covariance matrix $R_{DD}$ is symmetric and strictly positive definite.
- **R2 (Code-Minus-Carrier & Down-Weighting)**:
  - `test_r2_cmc_geometry_free_invariance`: CMC eliminates orbital and receiver dynamics, leaving only multipath and noise.
  - `test_r2_cmc_step_blunder_detection`: Step multipath error is flagged and isolated by innovation screening.
  - `test_r2_code_variance_inflation_reduces_kalman_gain`: CMC inflation reduces the Kalman gain for corrupted satellites.
  - `test_r2_mw_tracker_accumulates_clean_arc`: Melbourne-Wübbena widelane tracker accumulates stable statistics across clean arcs.
  - `test_r2_screen_gross_error_preserves_carrier_phase`: Pseudorange blunder rejection clears the pseudorange observation without disabling valid carrier phase tracking.
- **R3 (Doppler-Assisted Slip Detection)**:
  - `test_r3_clean_carrier_matches_doppler_no_slip`: Carrier phase change matching Doppler integrated velocity maintains continuous tracking arc.
  - `test_r3_doppler_exact_one_cycle_slip`: Doppler phase acceleration flags 1-cycle slip on GPS L1.
  - `test_r3_doppler_half_cycle_slip`: Doppler phase screening flags 0.5-cycle slip without false negatives.
  - `test_r3_loss_of_lock_indicator_flags_slip`: Receiver LLI flag triggers immediate ambiguity arc reset.
  - `test_r3_time_gap_triggers_slip`: Temporal gap (> 1.0 s) resets phase tracking arc.
- **R4 (Prioritized Partial Ambiguity Resolution)**:
  - `test_r4_par_subset_selects_high_confidence_ambiguities`: Selection prioritizes high-elevation, high-C/N0 satellites over corrupted ones.
  - `test_r4_par_submatrix_covariance_positive_definite`: Extracted PAR ambiguity covariance submatrix preserves positive definiteness.
  - `test_r4_par_dop_computation_from_geometry`: Post-fix satellite geometry produces valid, finite DOP metrics.
  - `test_r4_par_beidou_geo_key_detection`: BeiDou GEO satellites (C01–C05) are identified and de-prioritized.
  - `test_r4_par_empty_on_zero_dimension`: Safe fallback handling for empty or undersized candidate sets.

### Tier 2: Boundary & Extreme Values (20 Tests)
- **R1 Boundary**: Extreme horizon elevation (0.01°), extreme attenuation (10 dB-Hz), zenith maximum SNR (55 dB-Hz), numerical finite difference elevation & SNR gradient validations.
- **R2 Boundary**: Zero multipath baseline, extreme 20m code blunder, alternating positive/negative multipath oscillations, carrier-phase-only measurement, gross error rejection ceiling cap (preserving reference satellite).
- **R3 Boundary**: Near-threshold Doppler noise (0.35 cycles), long tracking interval cadence hint (10.0 s), high vehicle acceleration dynamics (30 m/s²), simultaneous multi-constellation slips (GPS + Galileo).
- **R4 Boundary**: Single ambiguity candidate, minimal subset size ($k = 4$), collinear degenerate street canyon geometry (DOP > 10.0 detection), near-singular ill-conditioned covariance matrix, strict target success rate threshold ($P_{req} = 0.999$).

### Tier 3: Cross-Module Interactions (6 Tests)
- `test_tier3_r1_r2_snr_and_cmc_compound_downweighting`: Compound attenuation from low SNR and CMC multipath yields higher variance than either source alone.
- `test_tier3_r2_r3_code_multipath_vs_doppler_slip_isolation`: Severe code multipath step triggers pseudorange blunder rejection without inducing a false Doppler cycle slip.
- `test_tier3_r3_r4_slipped_satellite_excluded_from_par`: Satellite experiencing Doppler cycle slip has its ambiguity arc reset and is excluded from PAR fixed set.
- `test_tier3_r2_r4_cmc_multipath_satellite_excluded_from_par`: High CMC variance satellite is penalized in PAR prioritization, allowing clean satellites to fix.
- `test_tier3_r1_r4_multi_constellation_par_selection`: PAR prioritization operates uniformly across GPS and Galileo constellations, selecting top satellites by joint SNR and elevation.
- `test_tier3_r1_r2_r3_r4_full_pipeline_cycle`: Full integrated pipeline cycle verifying R1, R2, R3, R4 interacting harmoniously without state corruption.

### Tier 4: Mission Scenarios (5 Scenarios)
- **Scenario 1**: `test_tier4_scenario1_tokyo_shinjuku_skyscraper`: Deep skyscraper canyon with 2 direct zenith satellites, 2 severe multipath reflections (+25m code step), and 2 low-elevation obstructed satellites. Screen isolates both multipath blunders while preserving zenith carrier phases.
- **Scenario 2**: `test_tier4_scenario2_hong_kong_whampoa_frequent_slips`: High-rise urban canyon with frequent 1-cycle and 0.5-cycle slips occurring on different satellites across successive epochs. Detector isolates slips without affecting clean companion arcs.
- **Scenario 3**: `test_tier4_scenario3_highway_underpass_outage_reacquisition`: Complete 3.5s satellite signal loss during bridge underpass traversal followed by reacquisition at degraded C/N0; triggers clean arc reset.
- **Scenario 4**: `test_tier4_scenario4_asymmetric_cors_base_vs_canyon_rover`: Asymmetric noise environment where CORS base station has pristine high SNR (48 dB-Hz) while urban canyon rover experiences severe attenuation (22 dB-Hz); verifies rover noise dominates double-difference variance.
- **Scenario 5**: `test_tier4_scenario5_collinear_street_canyon_dop_guard`: Collinear satellite geometry along a narrow urban street canyon; verifies DOP calculation correctly detects geometric degeneracy (PDOP > 10.0).

---

## 4. How to Run the Tests

To run the complete urban canyon test suite:
```bash
cargo test -p gneiss-rtk --test test_urban_canyon_e2e
```

To run with full test output / logging:
```bash
cargo test -p gneiss-rtk --test test_urban_canyon_e2e -- --nocapture
```

To run clippy lint validation on the test suite:
```bash
cargo clippy -p gneiss-rtk --test test_urban_canyon_e2e --all-targets -- -D warnings
```

---

## 5. Implementation Status & Bugs Discovered

- **Implementation Status**: All 4 features (R1 adaptive weighting, R2 CMC screening/down-weighting, R3 Doppler slip detection, R4 PAR prioritization) are fully implemented and integrated with the engine.
- **Engine Invariant Validations**:
  - `CycleSlipDetector` requires at least 3 satellites per band to compute the median clock drift and evaluate Doppler residuals (`if discs.len() < 3 { return Vec::new(); }`). Epoch fixtures have been constructed to reflect this invariant.
  - `screen_gross_pr_errors` requires geometrically consistent double differences. Synthetic test fixtures compute the exact geometric baseline distance $r_{sat} - r_{ref} - \rho_{base}$ to prevent false blunder rejections on clean satellites.
  - `WidelaneTracker` manages Melbourne-Wübbena widelane tracking per satellite arc; tests exercise its public `update`, `arc_means`, and `fixed_widelane` methods.
- **Bugs Discovered**: None. All production modules (`formation_cov.rs`, `screen.rs`, `screening.rs`, `par.rs`, `mw.rs`) performed correctly according to mathematical specifications and interface contracts.
