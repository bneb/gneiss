# Task Assignment — E2E Testing Track Orchestrator

## Role & Mission
You are the E2E Testing Track Sub-Orchestrator (`teamwork_preview_orchestrator`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/sub_orch_e2e/`.
Your parent orchestrator is: `1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`.

## Authoritative Documents to Read Before Starting Work
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md` (project architecture, feature inventory, interface contracts)
3. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code quality standards)

## Scope & Track Objectives
Build a comprehensive, opaque-box, requirement-driven E2E test suite across Tiers 1–4 independently of implementation internals:
1. **Create `TEST_INFRA.md`** at `/Users/kevin/projects/gneiss/.agents/sub_orch_e2e/TEST_INFRA.md` covering test architecture, runner commands, and coverage thresholds.
2. **Design and Implement Test Cases**:
   - **Tier 1 (Feature Coverage)**: $\ge 5$ test cases per feature in `PROJECT.md § Feature Inventory` covering representative inputs in isolation.
   - **Tier 2 (Boundary & Corner Cases)**: $\ge 5$ boundary/corner cases per feature (empty inputs, max limits, zero/negative, NaN/inf guards, extreme geometry).
   - **Tier 3 (Cross-Feature Combinations)**: Pairwise interaction tests (e.g. ESKF + NHC, PPP-AR + OSB + multi-constellation, VRS + multi-CORS + Delaunay).
   - **Tier 4 (Real-World Application Scenarios)**: $\ge 5$ realistic end-to-end mission workflows (Odaiba urban canyon GNSS/INS, F9P kinematic vehicle PPP-AR, regional CORS VRS network RTK).
3. **Publish `TEST_READY.md`** at `/Users/kevin/projects/gneiss/.agents/sub_orch_e2e/TEST_READY.md` and copy/symlink to project root `/Users/kevin/projects/gneiss/TEST_READY.md` when the full test suite is implemented and ready.

## Exclusive Write Ownership
You and your dispatched test writers/workers own:
- `tests/e2e/` or `crates/gneiss-tests/tests/e2e_frontiers/` (new E2E test suites)
- `/Users/kevin/projects/gneiss/.agents/sub_orch_e2e/`
- `TEST_INFRA.md` and `TEST_READY.md`
Do NOT modify production source code in `crates/gneiss-rtk/src/` or `crates/gneiss-parsers/src/`.

## Quality Invariants
All test code must adhere to AGENTS.md: file size < 500 LOC, function size < 32 LOC, nesting < 3 levels, 0 unwrap in non-test helpers, 0 compiler/clippy warnings.

When the suite is ready, publish `TEST_READY.md`, write `handoff.md`, and notify parent orchestrator (`1bd6ce81-03bf-4c40-b8b1-3b137333b5e7`).
