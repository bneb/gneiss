## 2026-06-22T01:20:08Z
You are teamwork_preview_worker (Worker 1 Gen 2).
Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_bug_16_gen2
Your task is to implement the fix for Bug 16: Mismatched Galileo BGD Correction.
Based on the Explorer's findings:
1. `position_e5b` needs to be defined on `Ephemeris` to delegate to `GalileoEphemeris::position_e5b` when constellation is Galileo, else `position`.
2. `get_frequency` in `crates/gneiss-core/src/signal.rs` needs to handle band 7 (returning `FREQ_GAL_E5B` for Galileo and Beidou).
3. `SppMeasurement` needs a `freq_band` field.
4. `build_single_measurement` needs to set `freq_band` based on observations (1 if p1_opt, 7 if Galileo/Beidou and p2_opt, 2 for other p2_opt, etc.).
5. Update `compute_sat_state` in `crates/gneiss-rtk/src/estimators/spp.rs` and `process_measurement` in `crates/gneiss-rtk/src/engine/spp_tight.rs` to call `m.eph.position_e5b` if `m.freq_band == 7`.
6. Add a regression unit test in `crates/gneiss-core/src/ephemeris.rs` or relevant test file. The regression test must pass now but would fail if the buggy code is restored.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Document the changes in /Users/kevin/projects/gneiss/.agents/worker_bug_16_gen2/handoff.md. Run `cargo test --workspace` using run_command to verify everything compiles and passes. Once complete, send a message to c1e1438e-1aa8-4425-a0b8-6dc0b1d37f21.
