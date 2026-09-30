# Reviewer 2 (Adversarial) Task Assignment — Milestone 2: CMC Multipath Mitigation

## Role & Mission
You are Reviewer 2 (`teamwork_preview_reviewer`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/reviewer_2_m2_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Input Documents
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, Section R2)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications & interface contracts)
3. `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md` (Worker M2 handoff report)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (code standards)

## Adversarial Focus
Adversarially scrutinize:
- Are there edge cases where `CmcTracker` fails to reset on cycle slips?
- What happens when carrier phase is absent vs present in `screen.rs`?
- Does `retain_active` prevent unbounded memory growth in `WidelaneTracker`?
- Are code-phase cross correlations strictly avoided so Kalman gain on phase is untouched?
- Check edge values: negative elevation, zero SNR, extreme code multipath jumps.
- Verify AGENTS.md rules (< 500 LOC/file, <= 32 LOC/function, < 3 nesting, 0 unwrap in prod, 0 warnings).

## Review Verification Commands
Execute:
```bash
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```

Write your verdict (APPROVE or REQUEST_CHANGES) and full report to `/Users/kevin/projects/gneiss/.agents/reviewer_2_m2_urban_canyon/handoff.md`.
Then send a completion message to your parent.

## 2026-09-24T19:27:47Z
You are Reviewer 2 (Adversarial) for Milestone 2: CMC Multipath Mitigation.
Your working directory is: /Users/kevin/projects/gneiss/.agents/reviewer_2_m2_urban_canyon/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/reviewer_2_m2_urban_canyon/DISPATCH.md, the worker handoff in /Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md, the user requirements in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md, and AGENTS.md.
Adversarially challenge edge cases (cycle slips, missing phase, memory leaks in retain_active, code-phase covariance independence, numerical stability).
Run:
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke

Deliver your verdict (APPROVE or REQUEST_CHANGES) in handoff.md and notify parent.
