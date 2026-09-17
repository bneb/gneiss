# Progress — Worker R2 (Frontier R2: Integer PPP-AR Engine via SINEX OSB)

Last visited: 2026-09-12T22:26:30Z

## Status
Initialized BRIEFING.md and progress.md. Ready to begin baseline check and codebase investigation.

## Milestones & Checklist
- [ ] 0. Run baseline build & tests to verify clean starting state.
- [ ] 1. Fast Indexed SINEX OSB Parser (`crates/gneiss-parsers/src/sinex_bia.rs`):
  - Indexing by `(Satellite, ObsCode)`.
  - Helpers `wide_lane_satellite_bias` & `narrow_lane_satellite_bias`.
  - Unit tests for indexed lookup & WL/NL bias computation.
- [ ] 2. Antenna PCO/PCV Corrections:
  - Satellite nadir-dependent PCV interpolation in `crates/gneiss-parsers/src/antex.rs` and projection in `satpos.rs`.
  - Standalone receiver antenna PCO/PCV evaluation for rover-only PPP in `crates/gneiss-parsers/src/receiver_pcv/`.
- [ ] 3. Multi-Constellation Support:
  - In `crates/gneiss-rtk/src/swfg/engine/epoch.rs`: Remove restrictive `is_supp` check.
  - In `satpos.rs` and `epoch.rs`: Include `'J'` QZSS constellation mapping and frequency selection.
- [ ] 4. Single-Differenced LAMBDA Integer Ambiguity Resolution:
  - Enhance `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` (Wide-Lane MW smoothing/rounding, Narrow-Lane LAMBDA).
  - Wire into `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`.
  - Back-substitute fixed integer ambiguities into the state estimator.
- [ ] 5. Benchmark Validation (`crates/gneiss-rtk/src/bin/eval_ppp.rs`):
  - Ingest `com21374.bia` (CODE MGEX OSB product) in `f9p_spec`.
  - Run `PPP_ONLY=f9p cargo run --release --bin eval_ppp`.
  - Verify sub-meter kinematic accuracy.
- [ ] 6. Comprehensive Verification:
  - `cargo build --workspace`
  - `cargo clippy --workspace --all-targets -- -D warnings`
  - Unit tests in `gneiss-parsers` and `gneiss-rtk`.
  - Handoff report generation.
