## 2026-09-13T11:57:19Z
You are worker_r2_resolution for the Gneiss positioning engine project.

Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_r2_resolution/
(Please initialize your directory, BRIEFING.md, and progress.md there).

You MUST read ORIGINAL_REQUEST.md before starting work:
Path: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Also read the project master plan:
Path: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Context from Prior Workers (worker_r2_final and worker_r2_fix):
1. WTZR 5-hour static PPP benchmark achieves sub-meter accuracy: p50 = 0.860m, RMS = 0.840m.
2. Integer PPP-AR LAMBDA fixing is verified working with ratios up to 1.74B - 1.9B!
3. Satellite Doppler transmit timing error caused by 4.775ms receiver clock bias in satpos interpolation was resolved.
4. East error is near zero (-0.06m to +0.20m).
5. On F9P kinematic vehicle drive, error vs CSRS-PPP / PPK ground truth is concentrated in ECEF Z (+2.07m), corresponding to +1.58m North and +1.28m Up at Boulder lat/lon.
6. Prefit residuals at truth reveal a North/South gradient of 3.25m (South satellites G03, G08 ~ +2.45m; North satellites G09, G26 ~ -0.80m).
7. Potential causes to investigate and resolve:
   - Antenna PCO / APC vs ARP / monument offset (does CSRS-PPP ground truth report ARP while the engine computes APC, or is receiver antenna PCO applied in the wrong direction/frame?). Check header/metadata in the RINEX / CSRS files.
   - Satellite PCO calculation/projection in `satpos.rs` (satellite nadir angle and body frame unit vectors).
   - Multi-constellation (Galileo / BeiDou) tracking in `epoch.rs` to balance satellite geometry.
   - Tropospheric mapping function or zenith delay modeling.

Your Mission:
1. Run `cargo run --release --bin eval_ppp` to inspect the exact current output and error metrics.
2. Diagnose and identify the root cause of the systematic North/Up residual on F9P kinematic drive.
3. Fix the underlying issue in `crates/gneiss-rtk/` and/or `crates/gneiss-parsers/` so that `eval_ppp` achieves sub-meter kinematic accuracy on the F9P drive dataset vs CSRS-PPP ground truth with integer ambiguities resolved.
4. Strictly satisfy AGENTS.md standards:
   - File size < 500 LOC
   - Function size < 32 LOC
   - Nesting depth < 3 levels
   - Zero `unwrap()` in production code
   - 0 compiler and clippy warnings
5. Run unit tests (`cargo test -p gneiss-rtk` and `cargo test -p gneiss-parsers`).
6. Run workspace tests and clippy:
   `cargo test --workspace`
   `cargo clippy --workspace --all-targets -- -D warnings`
7. Run regression guards:
   `python3 scripts/check_network_benchmark.py --smoke`
   `python3 scripts/check_multignss_benchmark.py --smoke`
8. Write your handoff report to `/Users/kevin/projects/gneiss/.agents/worker_r2_resolution/handoff.md` with:
   - Verbatim commands and output of `eval_ppp`.
   - Technical analysis of the root cause and fix.
   - Passing test results, clippy status, and regression guard results.
   - AGENTS.md compliance proof (< 500 LOC, < 32 LOC/fn, < 3 nesting, 0 unwrap).
9. Send a message to parent when complete using send_message.

## 2026-09-13T13:36:21Z
**Context**: Liveness status check
**Content**: Your progress.md has not been updated since 12:00:00Z. Are you running a long build/test command or stuck on eval_ppp?
**Action**: Please report your current execution status and update your progress.md.
