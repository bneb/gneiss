# Task Assignment — Sub-Orchestrator Milestone M3: Network RTK VRS Atmospheric Engine

## Role & Mission
You are a sub-orchestrator (`teamwork_preview_orchestrator`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/sub_orch_m3/`.
Your parent orchestrator is: `1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`.

## Authoritative Documents to Read Before Starting Work
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md` (project architecture, interface contracts, and file layout)
3. `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/survey_r3_r4.md` (comprehensive codebase analysis, Delaunay formulations, and benchmark investigation)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code quality standards)

## Scope & Milestone Objectives
Implement Frontier R3: Network RTK Virtual Reference Station (VRS) Atmospheric Engine:
1. Spatial 2D Delaunay Triangulation (`crates/gneiss-rtk/src/spatial/delaunay.rs`): Implement Bowyer-Watson Delaunay triangulation with point location and continuous barycentric interpolation for arbitrary 2D point sets.
2. Multi-station regional CORS base station ingestion: Ingest 5–10 regional CORS base station observation streams simultaneously from `datasets/cors_sf_bay_network`.
3. Multi-baseline double-difference network adjustment (`crates/gneiss-rtk/src/post_process/network_adj.rs`): Formulate network adjustment over inter-station baselines to solve for double-difference integer ambiguities across the network.
4. Spatial atmospheric modeling: Apply Delaunay triangulation models for ionospheric pierce points (IPP) per satellite and tropospheric zenith wet delay (ZWD) gradients across CORS stations.
5. Localized Virtual Reference Station (VRS) observation synthesis (`crates/gneiss-rtk/src/post_process/vrs.rs`): Synthesize localized VRS observation data at rover approximate position, reducing effective baseline length to $< 1\text{ km}$ on 15–50 km regional networks.
6. Benchmark validation: `crates/gneiss-rtk/src/bin/eval_network_ppk.rs` evaluating VRS-synthesized baseline against P181, P222, P225, demonstrating reduced baseline ppm error toward Leica single-baseline specifications ($8\text{ mm} + 1\text{ ppm}$).
7. CI Regression Invariant: Both regression guard scripts (`python3 scripts/check_network_benchmark.py --smoke` and `python3 scripts/check_multignss_benchmark.py --smoke`) MUST continue to pass with `ALL CHECKS PASSED`. (Do not alter stdout format expected by the guard scripts!).

## Exclusive Write Ownership
You and your dispatched workers own ONLY:
- `crates/gneiss-rtk/src/spatial/` (`mod.rs`, `delaunay.rs`)
- `crates/gneiss-rtk/src/post_process/network_adj.rs`
- `crates/gneiss-rtk/src/post_process/vrs.rs`
- `crates/gneiss-rtk/src/post_process/mod.rs`
- `crates/gneiss-rtk/src/bin/eval_network_ppk.rs`
Do NOT write to files owned by other milestones.

## Orchestrator Procedure & Gate Verification
Execute via the standard iteration loop: Explorer -> Worker -> Reviewer -> Challenger -> Forensic Auditor -> Gate.
- Worker dispatch prompt MUST include the mandatory integrity warning verbatim.
- Forensic Auditor verdict is a non-negotiable binary veto.
- All code must satisfy AGENTS.md (< 500 LOC/file, < 32 LOC/function, < 3 nesting, 0 unwrap in prod, 0 warnings).
- Regression guard scripts MUST pass cleanly.
- `cargo run --release --bin eval_network_ppk` must demonstrate reduced ppm error across CORS baselines.

When complete, write `handoff.md` and notify parent orchestrator (`1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`).
