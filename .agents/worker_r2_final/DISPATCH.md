# Dispatch: Frontier R2 (Integer PPP-AR Engine & eval_ppp Kinematic Benchmark)

## Identity & Context
- Archetype: worker (`teamwork_preview_worker`)
- Working directory: /Users/kevin/projects/gneiss/.agents/worker_r2_final/
- Parent conversation ID: cb1b403c-0ad6-4097-a034-1f63638ea37e
- Project scope document: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
- Authoritative user request: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. An independent forensic auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Objective
Complete final verification of Frontier R2:
1. Inspect the current state of Frontier R2 code in:
   - `crates/gneiss-rtk/src/bin/eval_ppp.rs`
   - `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
   - `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`
   - `crates/gneiss-rtk/src/swfg/engine/epoch.rs`
   - `crates/gneiss-parsers/src/sinex_bia.rs`
2. Note progress from previous worker:
   - AR fixing working with ratio up to 1.9B!
   - Resolved satellite Doppler transmit timing error caused by 4.775ms receiver clock bias in satpos interpolation.
   - East error is down to 0.028m.
   - Remaining: Tuning remaining North/Up offset on `eval_ppp.rs` (e.g. antenna PCO/APC offset vs monument/ARP, antenna model, or reference frame offset) to reach sub-meter kinematic accuracy vs CSRS-PPP on F9P kinematic drive.
3. Run `cargo run --bin eval_ppp --release` to see the current output and metrics.
4. Finalize the tuning/fixing so that `eval_ppp` achieves sub-meter kinematic accuracy on the F9P drive dataset vs CSRS-PPP ground truth with integer ambiguities resolved.
5. Strictly adhere to AGENTS.md standards:
   - File size < 500 LOC
   - Function size < 32 LOC
   - Nesting depth < 3 levels
   - Zero `unwrap()` in production code
   - 0 compiler and clippy warnings
6. Run `cargo test -p gneiss-rtk` and ensure all tests pass.
7. Write your handoff report to `/Users/kevin/projects/gneiss/.agents/worker_r2_final/handoff.md` and send a message back to parent when done.
