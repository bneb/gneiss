# Task Assignment — Worker M5: Final Integrated E2E Benchmark Validation & Code Standards Audit

## Role & Mission
You are the final verification and benchmark validation worker (`teamwork_preview_worker`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/worker_m5_final_benchmarks/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Input Documents
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (acceptance criteria)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications)
3. `/Users/kevin/projects/gneiss/.agents/worker_m1/handoff.md` (Milestone 1 handoff)
4. `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md` (Milestone 2 handoff)
5. `/Users/kevin/projects/gneiss/.agents/worker_m3_urban_canyon/handoff.md` (Milestone 3 handoff)
6. `/Users/kevin/projects/gneiss/.agents/worker_m4_urban_canyon/handoff.md` (Milestone 4 handoff)
7. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code standards)

## Scope & Verification Tasks (Milestone 5 / Acceptance Criteria)

1. **Compilation and Clippy Lints (0 warnings)**:
   ```bash
   cargo clippy --workspace --all-targets -- -D warnings
   ```
   Must pass with zero compiler warnings and zero clippy warnings.

2. **Workspace & E2E Test Suites (100% passing)**:
   ```bash
   cargo test --workspace
   cargo test --test test_urban_canyon_e2e
   ```
   All workspace tests and all 51 urban canyon E2E tests must pass.

3. **Dual CI Smoke Guard Scripts**:
   ```bash
   python3 scripts/check_network_benchmark.py --smoke
   python3 scripts/check_multignss_benchmark.py --smoke
   ```
   Both must output `ALL CHECKS PASSED` with exit code 0.

4. **Tokyo Odaiba 12,398-Epoch INS Benchmark**:
   Run:
   ```bash
   cargo run --release --bin eval_odaiba_ins
   ```
   Verify:
   - $p_{50} \le 2.134$ m
   - $\text{RMS} \le 4.156$ m

5. **UrbanNav Kinematic Rover Matrix (`eval_f9p_rover`)**:
   Run:
   ```bash
   MAX_EPOCHS=200 cargo run --release --bin eval_f9p_rover -- odaiba
   MAX_EPOCHS=200 cargo run --release --bin eval_f9p_rover -- shinjuku
   MAX_EPOCHS=200 cargo run --release --bin eval_f9p_rover -- whampoa_survey
   ```
   Verify fix rates, position accuracy, tail error suppression, and zero false fixes.

6. **AGENTS.md Audit on all Modified Production Files**:
   Audit the following files touched during the sprint:
   - `crates/gneiss-core/src/obs.rs`
   - `crates/gneiss-core/src/variance.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/mod.rs`
   - `crates/gneiss-rtk/src/post_process/screening.rs`
   - `crates/gneiss-rtk/src/ambiguity/par.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs`
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`
   - `crates/gneiss-rtk/src/composite/tc_ambiguity.rs`

   Confirm:
   - File size strictly < 500 LOC
   - Function size strictly <= 32 LOC
   - Nesting depth strictly < 3 levels
   - Zero `unwrap()` calls in production code
   - Zero dead code / warnings

Record all quantitative outputs and verdicts in `/Users/kevin/projects/gneiss/.agents/worker_m5_final_benchmarks/handoff.md`.
Then notify parent with a concise message.

## 2026-09-25T07:21:19Z
You are Worker M5 for Milestone 5: Final Integrated E2E Benchmark Validation & Code Standards Audit.
Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_m5_final_benchmarks/
Read your task assignment in /Users/kevin/projects/gneiss/.agents/worker_m5_final_benchmarks/DISPATCH.md, the authoritative user request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md, the project specifications in /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md, and all previous milestone handoffs (worker_m1, worker_m2_urban_canyon, worker_m3_urban_canyon, worker_m4_urban_canyon).

Execute the full verification and benchmark matrix:
1. cargo clippy --workspace --all-targets -- -D warnings (0 warnings)
2. cargo test --workspace and cargo test --test test_urban_canyon_e2e (all pass)
3. python3 scripts/check_network_benchmark.py --smoke (ALL CHECKS PASSED)
4. python3 scripts/check_multignss_benchmark.py --smoke (ALL CHECKS PASSED)
5. cargo run --release --bin eval_odaiba_ins (verify p50 <= 2.134m, RMS <= 4.156m)
6. MAX_EPOCHS=200 cargo run --release --bin eval_f9p_rover across odaiba, shinjuku, and whampoa_survey
7. Conduct exhaustive AGENTS.md audit across all modified files (< 500 LOC, <= 32 LOC/fn, < 3 nesting, 0 unwrap in prod).

Deliver your complete handoff report to /Users/kevin/projects/gneiss/.agents/worker_m5_final_benchmarks/handoff.md and notify your parent.

