# Forensic Auditor Task Assignment — Milestone 2: CMC Multipath Mitigation

## Role & Mission
You are the Forensic Integrity Auditor (`teamwork_preview_auditor`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/auditor_m2_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Input Documents
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, Section R2)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications & interface contracts)
3. `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md` (Worker M2 handoff report)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (code standards)

## Audit Scope
Perform exhaustive forensic integrity checks on Worker M2's implementation:
- `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`
- `crates/gneiss-rtk/tests/urban_canyon/tier1_features.rs`

## Forensic Verification Checks
1. Static analysis: Check for any hardcoded test values, magic strings tailored to specific tests, dummy or mock implementations, unwrap() calls in production code.
2. Code authenticity: Verify that mathematical formulas for CMC ($P - \Phi - 2I$), baseline accumulation, baseline freezing under multipath, pseudorange variance inflation ($R_{PP} + \sigma_{mp}^2$), and MW shielding are genuine, working physics-based algorithms.
3. Verification execution:
```bash
git diff --stat
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```

Write your verdict (CLEAN or INTEGRITY VIOLATION) and full report to `/Users/kevin/projects/gneiss/.agents/auditor_m2_urban_canyon/handoff.md`.
Then send a completion message to your parent.

## 2026-09-24T19:27:47Z
You are the Forensic Integrity Auditor for Milestone 2: CMC Multipath Mitigation.
Your working directory is: /Users/kevin/projects/gneiss/.agents/auditor_m2_urban_canyon/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/auditor_m2_urban_canyon/DISPATCH.md, the worker handoff in /Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md, and AGENTS.md.
Perform exhaustive forensic integrity checks: verify no hardcoding, no mock implementations, genuine physics-based algorithms, zero unwrap in prod, files < 500 LOC, functions <= 32 LOC.
Run:
git diff --stat
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke

Deliver your verdict (CLEAN or INTEGRITY VIOLATION) in handoff.md and notify parent.
