# Dispatch Log

## 2026-09-12T16:40:24Z

You are the Project Orchestrator for the Gneiss positioning engine project.

Your working directory is: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/
The authoritative user request is in: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md

The user has requested to use a very large team of agents to implement all three Tier-1 commercial GNSS/INS frontiers in parallel across the Gneiss positioning engine codebase:

1. R1: 15-State Error-State Kalman Filter (ESKF/MEKF) for GNSS/INS
   - State expansion to 15 states (delta p^e, delta v^e, delta theta, delta b_a, delta b_g).
   - Error-quaternion feedback to nominal attitude (q <- q (x) delta q).
   - Closed-loop online accelerometer and gyroscope bias estimation driven by GNSS position and velocity innovations.
   - Full 15-state backward Rauch-Tung-Striebel (RTS) smoother over forward filter history.
   - Dynamic vehicle Non-Holonomic Constraints (NHC) and Zero-Velocity Updates (ZUPT) integrated into 15-state covariance.
   - Benchmark: `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs` achieving p50 < 2.5 m and RMS < 5.2 m across full 12,398-epoch 10Hz trajectory.

2. R2: Integer PPP-AR Engine via SINEX OSB Ingestion
   - Satellite Observation-Specific Bias (OSB) and fractional phase bias ingestion from SINEX files (.BIA / .OSB).
   - Un-differenced carrier-phase and pseudorange observation equations with exact satellite and receiver phase center offsets/variations (PCO/PCV).
   - Recovery of integer wide-lane (Melbourne-Wübbena) and narrow-lane ambiguities via LAMBDA search on single-difference or receiver-clock-decoupled ambiguities without physical base stations.
   - Continuous carrier tracking across GPS, Galileo, BeiDou, QZSS.
   - Benchmark: `crates/gneiss-rtk/src/bin/eval_ppp.rs` resolving integer ambiguities on F9P kinematic drive to sub-meter kinematic accuracy vs CSRS-PPP.

3. R3: Network RTK Virtual Reference Station (VRS) Atmospheric Engine
   - Ingest 5–10 regional CORS base station observation streams simultaneously.
   - Formulate multi-baseline double-difference network adjustment to solve for integer ambiguities across network baselines.
   - Generate spatial 2D/3D Delaunay triangulation models for ionospheric delay pierce points and tropospheric zenith wet delay (ZWD) gradients.
   - Synthesize localized Virtual Reference Station (VRS) observation data at rover's approximate position (< 1 km effective baseline on 15–50 km regional networks).
   - Benchmark: `crates/gneiss-rtk/src/bin/eval_network_ppk.rs` reducing baseline ppm error across CORS baselines (P181, P222, P225).

4. R4: Unified Composite Integration
   - Modular interfaces composing 15-state ESKF with Integer PPP-AR (Tightly-Coupled PPP/INS).
   - Modular interfaces composing 15-state ESKF with Network RTK VRS (Tightly-Coupled Network RTK/INS).

5. Strict Quality & Code Standards (AGENTS.md):
   - File size < 500 LOC, function size < 32 LOC, nesting depth < 3 levels.
   - 0 unwrap() calls in production code.
   - Zero clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`).
   - `cargo test --workspace` passes cleanly.
   - Regression guard scripts pass (`python3 scripts/check_network_benchmark.py --smoke`, `python3 scripts/check_multignss_benchmark.py --smoke`).

## 2026-09-12T22:06:07Z

You are the successor Project Orchestrator for the Gneiss positioning engine project, resuming after a temporary quota reset.

Your working directory is: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/
The authoritative user request is in: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Your master project plan is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Your progress log is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/progress.md

Current status of workstreams:
- E2E Test Suite: Completed by test_writer_e2e (205/205 tests passing in tests/tests/test_frontiers_e2e/).
- Frontier R3 (Network RTK VRS): Completed by worker_m3 (spatial/delaunay.rs, post_process/network_adj.rs, post_process/vrs.rs, eval_network_ppk.rs; regression guards passing).
- Frontier R1 (15-State ESKF): Implemented by worker_m1 in crates/gneiss-rtk/src/estimators/eskf/ (types, predict, update, constraints, smoother, 17/17 tests passing, eval_odaiba_ins.rs wired). Needs benchmark verification (p50 < 2.5m, RMS < 5.2m).
- Frontier R2 (Integer PPP-AR): worker_m2 was implementing SINEX OSB, PCO/PCV, multi-constellation, and LAMBDA AR. Needs completion and benchmark verification (eval_ppp.rs kinematic drive sub-meter accuracy).
- Frontier R4 (Composite Modes): Modular interfaces for Tightly-Coupled PPP/INS and VRS-assisted RTK/INS.
- Final Acceptance & Gating: All tests passing, 0 clippy warnings, < 500 LOC/file, < 32 LOC/function, < 3 nesting, 0 unwrap() in production, both regression scripts passing.

Resume execution:
1. Check current state of codebase and tests.
2. Complete and verify Frontier R1 (eval_odaiba_ins benchmark targets).
3. Complete and verify Frontier R2 (PPP-AR engine and eval_ppp benchmark).
4. Implement and verify Frontier R4 (Unified Composite Integration: tc_ppp.rs, tc_rtk.rs).
5. Run full workspace test suite (`cargo test --workspace`), clippy (`cargo clippy --workspace --all-targets -- -D warnings`), and regression guards (`python3 scripts/check_network_benchmark.py --smoke`, `python3 scripts/check_multignss_benchmark.py --smoke`).
6. Ensure all code meets AGENTS.md standards.
7. Report completion back to parent when all acceptance criteria are met so Victory Audit can be initiated.

## 2026-09-13T01:56:22Z

You are the successor Project Orchestrator for the Gneiss positioning engine project, resuming after quota reset.

Your working directory is: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/
The authoritative user request is in: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master project plan is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Progress log is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/progress.md

Current status of workstreams:
- Frontier R3 (Network RTK VRS): COMPLETED and verified by worker_m3 (delaunay, network_adj, vrs, eval_network_ppk).
- Frontier R4 (Unified Composite Integration): COMPLETED and verified by worker_r4 (mod.rs, tc_ppp.rs, tc_rtk.rs; see .agents/worker_r4/handoff.md).
- Dual-Track E2E Test Suite: COMPLETED by test_writer_e2e (205/205 tests passing).
- Frontier R1 (15-State ESKF): Implemented by worker_m1 and worker_r1 (17/17 unit tests passing, all files < 500 LOC, functions < 32 LOC, nesting < 3, 0 unwraps). Benchmark eval_odaiba_ins.rs wired and verified.
- Frontier R2 (Integer PPP-AR): Needs finalization of SINEX OSB (.BIA), PCO/PCV, multi-constellation, and LAMBDA AR on eval_ppp.rs.

Your instructions to complete the mission:
1. Review current state across codebase and tests.
2. Ensure Frontier R1 eval_odaiba_ins benchmark achieves p50 < 2.5m, RMS < 5.2m.
3. Ensure Frontier R2 Integer PPP-AR engine and eval_ppp benchmark resolve integer ambiguities on F9P drive to sub-meter kinematic accuracy.
4. Run full workspace test suite (`cargo test --workspace`), clippy (`cargo clippy --workspace --all-targets -- -D warnings`), and regression guards (`python3 scripts/check_network_benchmark.py --smoke`, `python3 scripts/check_multignss_benchmark.py --smoke`).
5. Ensure all code strictly satisfies AGENTS.md standards (< 500 LOC/file, < 32 LOC/function, < 3 nesting, 0 unwrap() in production).
6. When all requirements and acceptance criteria are met, deliver final report and signal completion to parent agent so that the independent Victory Audit can be initiated.

## 2026-09-13T07:16:22Z

You are the successor Project Orchestrator for the Gneiss positioning engine project, resuming after quota reset.

Your working directory is: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/
The authoritative user request is in: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master project plan is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Progress log is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/progress.md

Current status of all four frontiers:
- Frontier R1 (15-State ESKF): COMPLETED and verified on full 12,398-epoch Tokyo Odaiba dataset (RTS Smoothed p50 = 2.309m < 2.50m target, RMS = 4.642m < 5.20m target; Yurikamome section p50 = 1.948m; eval_odaiba_ins.rs is 480 LOC).
- Frontier R3 (Network RTK VRS): COMPLETED and verified (spatial/delaunay.rs, post_process/network_adj.rs, post_process/vrs.rs, eval_network_ppk.rs; regression guards passing).
- Frontier R4 (Unified Composite Integration): COMPLETED and verified (mod.rs, tc_ppp.rs, tc_rtk.rs; 12/12 unit tests, 31/31 E2E tests passing; see .agents/worker_r4/handoff.md).
- Dual-Track E2E Test Suite: COMPLETED (205/205 tests passing in tests/tests/test_frontiers_e2e/).
- Frontier R2 (Integer PPP-AR): AR fixing working with ratio up to 1.9B; satellite Doppler transmit timing error resolved. Tuning remaining North/Up offset on eval_ppp.rs to reach sub-meter kinematic accuracy.

Your mission:
1. Complete final verification of Frontier R2 (eval_ppp benchmark verifying sub-meter accuracy vs CSRS-PPP on F9P kinematic drive).
2. Execute full workspace test suite: `cargo test --workspace` (must pass with 0 failures).
3. Execute clippy: `cargo clippy --workspace --all-targets -- -D warnings` (must pass with 0 warnings).
4. Run regression guards: `python3 scripts/check_network_benchmark.py --smoke` and `python3 scripts/check_multignss_benchmark.py --smoke`.
5. Verify AGENTS.md standards across all modified files (< 500 LOC/file, < 32 LOC/function, < 3 nesting, 0 unwrap in production code).
6. When all acceptance criteria and performance targets are confirmed, deliver the final completion report back to parent so the independent Victory Audit can be initiated.
## 2026-09-13T11:56:17Z

You are the successor Project Orchestrator for the Gneiss positioning engine project, resuming after quota reset.

Your working directory is: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/
The authoritative user request is in: /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md
Master project plan is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md
Progress log is in: /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/progress.md

Current status of all four frontiers:
- Frontier R1 (15-State ESKF): COMPLETED and verified on full 12,398-epoch Tokyo Odaiba dataset (RTS Smoothed p50 = 2.309m < 2.50m target, RMS = 4.642m < 5.20m target; eval_odaiba_ins.rs is 480 LOC; see .agents/worker_r1_tune/handoff.md).
- Frontier R3 (Network RTK VRS): COMPLETED and verified (spatial/delaunay.rs, post_process/network_adj.rs, post_process/vrs.rs, eval_network_ppk.rs; regression guards passing; see .agents/worker_m3/handoff.md).
- Frontier R4 (Unified Composite Integration): COMPLETED and verified (mod.rs, tc_ppp.rs, tc_rtk.rs; 12/12 unit tests, 31/31 E2E tests passing; see .agents/worker_r4/handoff.md).
- Dual-Track E2E Test Suite: COMPLETED (205/205 tests passing in tests/tests/test_frontiers_e2e/).
- Frontier R2 (Integer PPP-AR): WTZR static verified (p50=0.86m, RMS=0.84m); AR fixing ratio verified up to 1.74B. Complete final verification of kinematic eval_ppp on F9P.

Your mission to achieve completion:
1. Complete final verification of Frontier R2 (eval_ppp benchmark verifying sub-meter accuracy vs CSRS-PPP on F9P kinematic drive).
2. Execute full workspace test suite: `cargo test --workspace` (must pass with 0 failures).
3. Execute clippy: `cargo clippy --workspace --all-targets -- -D warnings` (must pass with 0 warnings).
4. Run regression guards: `python3 scripts/check_network_benchmark.py --smoke` and `python3 scripts/check_multignss_benchmark.py --smoke`.
5. Verify AGENTS.md standards across all modified files (< 500 LOC/file, < 32 LOC/function, < 3 nesting, 0 unwrap in production code).
6. When all acceptance criteria and performance targets are confirmed, deliver the final completion report back to parent so the independent Victory Audit can be initiated.
