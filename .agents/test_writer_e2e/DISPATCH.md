# Task Assignment — E2E Test Writer (Tiers 1–4)

## Role & Mission
You are a test writer (`teamwork_preview_test_writer`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/test_writer_e2e/`.
The authoritative user request is: `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`.
The project specification is: `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md`.

## Quality Standards (AGENTS.md)
- File size strictly < 500 LOC
- Function size strictly < 32 LOC
- Nesting depth strictly < 3 levels
- Zero clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)

## Exclusive Write Ownership
You own:
- `/Users/kevin/projects/gneiss/.agents/test_writer_e2e/`
- `crates/gneiss-tests/tests/test_frontiers_e2e.rs` (or modular test files under `crates/gneiss-tests/tests/`)
- `TEST_INFRA.md` at project root or working directory
- `TEST_READY.md` at project root or working directory
Do NOT modify production code in `crates/gneiss-rtk/src/` or `crates/gneiss-parsers/src/`.

## Scope & Deliverables
Create a comprehensive, opaque-box, requirement-driven E2E test suite across Tiers 1–4 derived directly from user requirements:
1. **`TEST_INFRA.md`**: Define test architecture, runner command, and coverage criteria following the template in `Project Pattern`.
2. **Tier 1 (Feature Coverage)**: $\ge 5$ test cases per feature in `PROJECT.md § Feature Inventory` (Features 1–19) exercising representative inputs in isolation.
3. **Tier 2 (Boundary & Corner Cases)**: $\ge 5$ test cases per feature testing limits, zeroes, empty inputs, max sizes, and domain extremes.
4. **Tier 3 (Cross-Feature Combinations)**: Pairwise interaction tests between major feature pairs (ESKF + NHC, PPP-AR + OSB + multi-constellation, VRS + multi-CORS + Delaunay).
5. **Tier 4 (Real-World Application Scenarios)**: $\ge 5$ realistic end-to-end mission workflows (Odaiba urban canyon GNSS/INS, F9P kinematic vehicle PPP-AR, regional CORS VRS network RTK).
6. **`TEST_READY.md`**: Publish `TEST_READY.md` with full coverage summary table and checklist when all tests are implemented.

## Verification Command
Run:
```bash
cargo test -p gneiss-tests --test test_frontiers_e2e
```

## 2026-09-12T16:48:52Z
You are E2E Test Writer.
Your working directory is: /Users/kevin/projects/gneiss/.agents/test_writer_e2e/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/test_writer_e2e/DISPATCH.md, the project scope in /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md, and the authoritative request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md.
Build the comprehensive opaque-box E2E test suite across Tiers 1-4 covering all features in the Feature Inventory. Publish TEST_INFRA.md and TEST_READY.md.
Run your verification commands and record all output in your handoff.md. Send a completion message when done.
