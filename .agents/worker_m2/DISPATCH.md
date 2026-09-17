# Task Assignment — Worker Milestone M2: Integer PPP-AR Engine via SINEX OSB

## Role & Mission
You are an implementation worker (`teamwork_preview_worker`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/worker_m2/`.
The authoritative user request is: `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`.
The project specification is: `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md`.
The survey report is: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2/survey_r2.md`.

## MANDATORY INTEGRITY WARNING
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Code Quality Standards (AGENTS.md)
- File size strictly < 500 LOC
- Function size strictly < 32 LOC
- Nesting depth strictly < 3 levels
- Exactly 0 `unwrap()` calls in production code (`match`, `if let`, `ok_or()?`, or descriptive `.expect()` only)
- Zero clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)
- Passing test coverage with comprehensive unit tests in `#[cfg(test)] mod tests`

## Exclusive Write Ownership
You own ONLY:
- `crates/gneiss-parsers/src/sinex_bia.rs`
- `crates/gneiss-parsers/src/antex.rs`
- `crates/gneiss-parsers/src/receiver_pcv/`
- `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
- `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`
- `crates/gneiss-rtk/src/swfg/engine/epoch.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/satpos.rs`
- `crates/gneiss-rtk/src/bin/eval_ppp.rs`
Do NOT edit any files outside this set.

## Implementation Requirements
Follow the exact mathematical equations from `survey_r2.md`:
1. **Fast Indexed SINEX OSB Parser (`crates/gneiss-parsers/src/sinex_bia.rs`)**:
   - Add indexing (e.g. `HashMap<(Satellite, ObsCode), Vec<BiasRecord>>` or interval tree) for $O(1)$ bias lookups.
   - Implement helper methods on `SinexBias`:
     - `wide_lane_satellite_bias(&self, sat: &Satellite, time: Epoch) -> Option<f64>`
     - `narrow_lane_satellite_bias(&self, sat: &Satellite, time: Epoch) -> Option<f64>`
2. **Antenna PCO/PCV Corrections**:
   - Add satellite nadir-dependent PCV interpolation in `crates/gneiss-parsers/src/antex.rs` and project along satellite line-of-sight in `satpos.rs`.
   - Add standalone receiver antenna PCO/PCV evaluation for rover-only PPP in `crates/gneiss-parsers/src/receiver_pcv/`.
3. **Multi-Constellation Support**:
   - In `crates/gneiss-rtk/src/swfg/engine/epoch.rs`: Remove restrictive `is_supp` check (line 364) that dropped BeiDou and QZSS.
   - In `satpos.rs` and `epoch.rs`: Include `'J'` QZSS constellation mapping and frequency selection.
4. **Single-Differenced LAMBDA Integer Ambiguity Resolution**:
   - Update `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` and wire into `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`:
     - Form single-differences between satellites to cancel the receiver phase bias.
     - Fix Wide-Lane (MW) ambiguities via integer rounding/bootstrap.
     - Fix Narrow-Lane ambiguities via LAMBDA search with ratio test validation.
     - Back-substitute fixed integer ambiguities into the state estimator.
5. **Benchmark Validation (`crates/gneiss-rtk/src/bin/eval_ppp.rs`)**:
   - Ingest `com21374.bia` (CODE MGEX OSB product) in `f9p_spec`.
   - Run `cargo run --release --bin eval_ppp`.
   - Verify integer ambiguity resolution on the F9P kinematic drive achieving sub-meter kinematic accuracy closing the discrepancy vs CSRS-PPP ($0.296\text{ m}$ RMS).

## Verification Commands
Before reporting completion, run:
```bash
cargo build --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-parsers --lib sinex_bia
cargo test -p gneiss-rtk --lib ambiguity::ppp_ar
PPP_ONLY=f9p cargo run --release --bin eval_ppp
```
Record all test and benchmark outputs in your `handoff.md`.
Send a completion message when done.

## 2026-09-12T16:48:52Z
You are Worker M2.
Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_m2/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/worker_m2/DISPATCH.md, the project scope in /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md, the survey report in /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2/survey_r2.md, and the authoritative request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md.
Implement Frontier R2 (Integer PPP-AR Engine via SINEX OSB, PCO/PCV, multi-constellation, single-differenced LAMBDA AR, eval_ppp benchmark) adhering strictly to AGENTS.md.
Run your verification commands and record all output in your handoff.md. Send a completion message when done.
