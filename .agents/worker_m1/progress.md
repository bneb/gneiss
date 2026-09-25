# Progress — Worker M1 (R1 Adaptive C/N0 & Elevation Observation Covariance)

Last visited: 2026-09-24T14:40:00Z

## Status
- **Milestone 1 Completed**:
  - `crates/gneiss-core/src/obs.rs`: Implemented `get_snr_f64` returning unquantized C/N0, preserved `get_snr`. Split long test function into two <= 32 LOC tests.
  - `crates/gneiss-core/src/variance.rs`: Implemented unified C1-smooth SIGMA-SNR elevation and noise model (`elevation_factor`, `snr_factor`, `sigma_snr_variance_with_coeffs`, `sigma_snr_variance`, `sigma_snr_phase_variance`) with Three-Tier unit tests (golden vectors, numerical derivatives, stability bounds).
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`: Refactored double-difference covariance calculation with `IntoSnrs` trait supporting `Option<f64>`, `Option<u8>`, `Option<i32>`, and smooth SIGMA-SNR scaling.
- **Verification status**:
  - `cargo clippy --workspace --all-targets -- -D warnings`: 0 warnings.
  - `cargo test -p gneiss-core`: 144 passed, 0 failed.
  - `cargo test -p gneiss-rtk --lib`: 435 passed, 0 failed.
  - `cargo test --test test_urban_canyon_e2e`: 51 passed, 0 failed.
  - `python3 scripts/check_network_benchmark.py --smoke`: ALL CHECKS PASSED.
  - `python3 scripts/check_multignss_benchmark.py --smoke`: ALL CHECKS PASSED.
  - `benchmark_matrix::test_real_geodetic_cors_baseline_tmg2_tmgo_sub_centimeter`: 37/40 fixed epochs (p50 = 4.4mm).
  - Code standards: all files < 500 LOC, all functions <= 32 LOC, nesting <= 3, 0 unwraps in production.
- Handoff report prepared in `.agents/worker_m1/handoff.md`.
