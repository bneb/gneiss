# Task Assignment — Survey Frontier R3/R4: Network RTK VRS Engine & Unified Composite Integration

## Role & Mission
You are a read-only exploration agent (`teamwork_preview_explorer`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/`.
The authoritative user request is: `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`.

## Objective
Survey the codebase for Frontier R3: Network RTK Virtual Reference Station (VRS) Atmospheric Engine and Frontier R4: Unified Composite Integration:
1. Ingestion of 5–10 regional CORS base station observation streams simultaneously.
2. Formulation of multi-baseline double-difference network adjustment to solve for integer ambiguities across network baselines.
3. Generation of spatial 2D/3D Delaunay triangulation models for ionospheric delay pierce points and tropospheric zenith wet delay (ZWD) gradients.
4. Synthesis of localized Virtual Reference Station (VRS) observation data at rover's approximate position (< 1 km effective baseline on 15–50 km regional networks).
5. Benchmark: `crates/gneiss-rtk/src/bin/eval_network_ppk.rs` reducing baseline ppm error across CORS baselines (P181, P222, P225).
6. Frontier R4: Unified Composite Integration — modular interfaces composing 15-state ESKF with Integer PPP-AR (Tightly-Coupled PPP/INS) and with Network RTK VRS (Tightly-Coupled Network RTK/INS).
7. Regression guard scripts: `scripts/check_network_benchmark.py` and `scripts/check_multignss_benchmark.py`.

## Scope of Investigation
- Locate and examine all existing Network RTK/PPK modules (`crates/gneiss-rtk/src/network_rtk.rs` or similar, `crates/gneiss-rtk/src/bin/eval_network_ppk.rs`), multi-station CORS baseline handling (P181, P224, P225, P222, SLAC, CAPO), iono/tropo models.
- Check regression guard scripts `scripts/check_network_benchmark.py` and `scripts/check_multignss_benchmark.py`: what tests/checks are performed, what parameters they expect, and current state.
- Check composite integration architecture: how ESKF, PPP-AR, and RTK/VRS can be cleanly and modularly integrated without violating file size (< 500 LOC) and function size (< 32 LOC) limits.
- Map exact data structures, functions, interfaces, dependencies, and files that need to be created or modified.
- Verify adherence to AGENTS.md constraints.

## Deliverable
Write your complete survey report to:
`/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/survey_r3_r4.md`
and write a standard `handoff.md` in your directory.
Send a message back to the orchestrator when complete.

## 2026-09-12T16:41:15Z
You are Survey Explorer R3/R4.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/
Please read your task instructions in /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/DISPATCH.md and the authoritative request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md.
Investigate the codebase for Frontier R3 (Network RTK VRS Atmospheric Engine, multi-CORS adjustment, Delaunay iono/tropo models, VRS synthesis, eval_network_ppk benchmark) and Frontier R4 (Unified Composite Integration: Tightly-Coupled PPP/INS and Network RTK/INS), as well as regression guard scripts (check_network_benchmark.py and check_multignss_benchmark.py).
Output your survey report to /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r3_r4/survey_r3_r4.md and your handoff.md. Send a completion message to the parent orchestrator when done.

