# Dispatch: E2E Test Writer (Dual-Track Test Suite)

## Objective
Design and implement a comprehensive opaque-box E2E test suite derived strictly from user requirements in `ORIGINAL_REQUEST.md` (2026-09-24T13:30:49Z) and the project specification in `PROJECT.md`.

## Required Reading
- `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`
- `/Users/kevin/projects/gneiss/AGENTS.md`
- `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/PROJECT.md`
- `/Users/kevin/projects/gneiss/.agents/survey_explorer_1/survey_r1_bench.md`
- `/Users/kevin/projects/gneiss/.agents/survey_explorer_2/survey_r2_cmc.md`
- `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md`

## Scope & Methodology
Design test cases using the 4-tier methodology:
1. **Tier 1 - Feature Coverage (>=5 tests per feature)**: Happy-path tests verifying each of the 4 requirements (R1, R2, R3, R4) in isolation.
2. **Tier 2 - Boundary & Corner Cases (>=5 tests per feature)**: Boundary conditions (extreme low elevation < 5 deg, low SNR < 20 dB-Hz, sudden 20m code multipath steps, half-cycle slips, poor satellite geometry/DOP).
3. **Tier 3 - Cross-Feature Combinations (pairwise interactions)**: Interactions between SNR weighting, CMC multipath suppression, Doppler slip resets, and prioritized PAR.
4. **Tier 4 - Real-World Application Scenarios (>=5 scenarios)**: Realistic urban canyon vehicle tracking scenarios mimicking Tokyo Shinjuku and Hong Kong Whampoa geometry.

## Mandatory Deliverables
1. `TEST_INFRA.md` (published to `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/TEST_INFRA.md`):
   - Test architecture, feature inventory, scenario definitions, coverage thresholds.
2. Comprehensive test file in `crates/gneiss-rtk/tests/test_urban_canyon_e2e.rs`.
   - File size < 500 LOC (split into helper modules or multiple test files if needed: e.g. `tests/urban_canyon/mod.rs`, `tier1.rs`, `tier2.rs`, etc.).
   - Function size <= 32 LOC, nesting < 3 levels, 0 unwrap in non-test production code.
3. Verification:
   - Run `cargo test --test test_urban_canyon_e2e` (or workspace tests).
   - Once all tests are written and passing, publish `TEST_READY.md` to `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/TEST_READY.md`.
4. Deliver `handoff.md` and send completion message to parent.

## 2026-09-24T13:53:08Z
You are the E2E Test Writer for the Gneiss Urban Canyon Fix Rate Expansion and Multipath Mitigation task.
Your working directory is: /Users/kevin/projects/gneiss/.agents/test_writer_e2e
Your parent conversation ID is: c1309e2d-6c95-4b14-a86d-d26a13f2a150

Read your instructions in:
- /Users/kevin/projects/gneiss/.agents/test_writer_e2e/DISPATCH.md
- /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md (specifically 2026-09-24T13:30:49Z)
- /Users/kevin/projects/gneiss/AGENTS.md
- /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/PROJECT.md
- /Users/kevin/projects/gneiss/.agents/survey_explorer_1/survey_r1_bench.md
- /Users/kevin/projects/gneiss/.agents/survey_explorer_2/survey_r2_cmc.md
- /Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md

Deliverables:
1. Create TEST_INFRA.md and save to /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/TEST_INFRA.md.
2. Implement 4-tier E2E tests in crates/gneiss-rtk/tests/test_urban_canyon_e2e.rs (and helper submodules if needed to keep file size < 500 LOC).
3. Verify tests compile and run with cargo test.
4. When test suite is ready and verified, create TEST_READY.md in /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/TEST_READY.md.
5. Deliver handoff.md and send a message back to parent.

