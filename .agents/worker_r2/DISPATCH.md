# Dispatch: Worker R2 — Frontier R2 (Integer PPP-AR Engine via SINEX OSB)

Working Directory: /Users/kevin/projects/gneiss/.agents/worker_r2
Original Request: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master Plan: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Previous Progress: /Users/kevin/projects/gneiss/.agents/worker_m2/progress.md

## Objective
Complete Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion.
1. Fast Indexed SINEX OSB Parser (`crates/gneiss-parsers/src/sinex_bia.rs`):
   - Fast lookup indexing by `(Satellite, ObsCode)`
   - Helpers: `wide_lane_satellite_bias` & `narrow_lane_satellite_bias`
   - Unit tests for indexed lookup and WL/NL bias computation
2. Antenna PCO/PCV Corrections:
   - Satellite nadir-dependent PCV interpolation in `antex.rs` and projection in `satpos.rs`
   - Standalone receiver antenna PCO/PCV evaluation for rover-only PPP in `crates/gneiss-parsers/src/receiver_pcv/`
3. Multi-Constellation Support:
   - In `crates/gneiss-rtk/src/swfg/engine/epoch.rs`: Remove restrictive `is_supp` check
   - In `satpos.rs` and `epoch.rs`: Include `'J'` QZSS constellation mapping and frequency selection
4. Single-Differenced LAMBDA Integer Ambiguity Resolution:
   - Enhance `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` (Wide-Lane MW smoothing/rounding, Narrow-Lane LAMBDA)
   - Wire into `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs` and back-substitute fixed integer ambiguities
5. Benchmark Validation (`crates/gneiss-rtk/src/bin/eval_ppp.rs`):
   - Ingest SINEX OSB bias product (`com21374.bia`) in `f9p_spec`
   - Run `PPP_ONLY=f9p cargo run --release --bin eval_ppp`
   - Target: Sub-meter kinematic accuracy on F9P drive vs CSRS-PPP (e.g. ~0.3 m RMS)

## Exclusive File Ownership
- `crates/gneiss-parsers/src/sinex_bia.rs`
- `crates/gneiss-parsers/src/antex.rs`
- `crates/gneiss-parsers/src/receiver_pcv/**`
- `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
- `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`
- `crates/gneiss-rtk/src/swfg/engine/epoch.rs`
- `crates/gneiss-rtk/src/bin/eval_ppp.rs`

## Mandatory Rules & Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Follow AGENTS.md standards:
- < 500 LOC per file
- < 32 LOC per function
- < 3 nesting depth
- 0 unwrap() in production code
- 0 clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)
- All tests pass (`cargo test --workspace`)

## 2026-09-12T22:25:50Z
You are Worker R2 for the Gneiss positioning engine project.
Complete Frontier R2: Integer PPP-AR Engine via SINEX OSB Ingestion.
