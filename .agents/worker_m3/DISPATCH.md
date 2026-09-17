# Task Assignment — Worker Milestone M3: Network RTK VRS Atmospheric Engine

## Role & Mission
You are an implementation worker (`teamwork_preview_worker`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/worker_m3/`.
The authoritative user request is: `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`.
The project specification is: `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md`.
The survey report is: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/survey_r3_r4.md`.

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
- `crates/gneiss-rtk/src/spatial/` (`mod.rs`, `delaunay.rs`)
- `crates/gneiss-rtk/src/post_process/network_adj.rs`
- `crates/gneiss-rtk/src/post_process/vrs.rs`
- `crates/gneiss-rtk/src/post_process/mod.rs`
- `crates/gneiss-rtk/src/bin/eval_network_ppk.rs`
Do NOT edit any files outside this set.

## Implementation Requirements
Follow the exact mathematical equations from `survey_r3_r4.md`:
1. **Spatial 2D Delaunay Triangulation (`crates/gneiss-rtk/src/spatial/delaunay.rs`)**:
   - Implement incremental Bowyer-Watson Delaunay triangulation with super-triangle initialization.
   - Point location and barycentric coordinate interpolation: $\lambda_1 p_1 + \lambda_2 p_2 + \lambda_3 p_3 = p$.
   - Comprehensive unit tests covering collinear points, duplicate points, and boundary cases.
2. **Multi-Baseline Double-Difference Network Adjustment (`crates/gneiss-rtk/src/post_process/network_adj.rs`)**:
   - Ingest 5–10 regional CORS stations simultaneously (`datasets/cors_sf_bay_network`).
   - Formulate inter-station baseline graph and double-difference integer ambiguity network adjustment.
3. **Spatial Atmospheric Modeling & VRS Observation Synthesis (`crates/gneiss-rtk/src/post_process/vrs.rs`)**:
   - Interpolate station tropospheric ZWD gradients and per-satellite ionospheric pierce points (IPP) via Delaunay triangulation.
   - Synthesize localized VRS observation data at rover approximate position, reducing effective baseline to $< 1\text{ km}$ on 15–50 km regional networks.
4. **Benchmark Validation (`crates/gneiss-rtk/src/bin/eval_network_ppk.rs`)**:
   - Integrate VRS synthesis with rover PPK.
   - Reduce baseline ppm error across CORS baselines (P181, P222, P225) toward Leica single-baseline specs ($8\text{ mm} + 1\text{ ppm}$).
   - **CRITICAL CI INVARIANT**: Both regression guard scripts (`python3 scripts/check_network_benchmark.py --smoke` and `python3 scripts/check_multignss_benchmark.py --smoke`) MUST continue to pass with `ALL CHECKS PASSED`. Do not alter stdout formatting expected by the regex matchers in those scripts!

## Verification Commands
Before reporting completion, run:
```bash
cargo build --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib spatial::delaunay
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
cargo run --release --bin eval_network_ppk
```
Record all test and benchmark outputs in your `handoff.md`.
Send a completion message when done.

## 2026-09-12T16:48:52Z
You are Worker M3.
Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_m3/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/worker_m3/DISPATCH.md, the project scope in /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md, the survey report in /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/survey_r3_r4.md, and the authoritative request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md.
Implement Frontier R3 (Network RTK VRS Atmospheric Engine, 2D Delaunay triangulation, multi-baseline network adjustment, VRS synthesis, eval_network_ppk benchmark) adhering strictly to AGENTS.md and ensuring regression guard scripts pass.
Run your verification commands and record all output in your handoff.md. Send a completion message when done.

