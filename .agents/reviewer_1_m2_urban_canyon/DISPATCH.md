# Reviewer 1 Task Assignment — Milestone 2: CMC Multipath Mitigation

## Role & Mission
You are Reviewer 1 (`teamwork_preview_reviewer`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/reviewer_1_m2_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Input Documents
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, Section R2)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications & interface contracts)
3. `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md` (Worker M2 handoff report)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (code standards)

## Review Scope
Review all files modified/created by Worker M2:
- `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`
- `crates/gneiss-rtk/tests/urban_canyon/tier1_features.rs`

## Review Verification Commands
Execute:
```bash
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```
Audit all files against AGENTS.md (< 500 LOC/file, <= 32 LOC/function, < 3 nesting, 0 unwrap in prod).

Write your verdict (APPROVE or REQUEST_CHANGES) and full report to `/Users/kevin/projects/gneiss/.agents/reviewer_1_m2_urban_canyon/handoff.md`.
Then send a completion message to your parent.

## 2026-09-24T19:27:47Z
<USER_REQUEST>
You are Reviewer 1 for Milestone 2: CMC Multipath Mitigation.
Your working directory is: /Users/kevin/projects/gneiss/.agents/reviewer_1_m2_urban_canyon/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/reviewer_1_m2_urban_canyon/DISPATCH.md, the worker handoff in /Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md, the user requirements in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md, and AGENTS.md.
Objectively review correctness, completeness, interface conformance, and AGENTS.md rules (< 500 LOC, <= 32 LOC/fn, < 3 nesting, 0 unwrap in prod, 0 warnings).
Run:
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke

Deliver your verdict (APPROVE or REQUEST_CHANGES) in handoff.md and notify parent.
</USER_REQUEST>
