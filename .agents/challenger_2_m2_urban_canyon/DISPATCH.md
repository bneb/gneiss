# Challenger 2 Task Assignment — Milestone 2: Melbourne-Wübbena Multipath Shielding Stress Testing

## Role & Mission
You are Challenger 2 (`teamwork_preview_challenger`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/challenger_2_m2_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Input Documents
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, Section R2)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications & interface contracts)
3. `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md` (Worker M2 handoff report)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (code standards)

## Challenge Scope
Empirically stress-test the Melbourne-Wübbena shielding and tracking invariants:
1. Verify `MwTrack::absorb()` under code multipath steps:
   - In unshielded mode, an innovation jump > 1.0 cycle resets `MwTrack`.
   - In shielded mode (`shielded: true`, triggered when CMC multipath is detected on the arc), verify that a 5m–15m pseudorange jump does NOT reset `MwTrack`, does NOT corrupt `mean`, and preserves accumulated `n` epochs.
2. Verify that genuine cycle slips (`slip: true`) still reset `WidelaneTracker` and `CmcTracker` unconditionally.
3. Run:
```bash
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```

Write your verdict (APPROVE or REQUEST_CHANGES) and full report to `/Users/kevin/projects/gneiss/.agents/challenger_2_m2_urban_canyon/handoff.md`.
Then send a completion message to your parent.

## 2026-09-24T19:27:47Z
You are Challenger 2 for Milestone 2: CMC Multipath Mitigation.
Your working directory is: /Users/kevin/projects/gneiss/.agents/challenger_2_m2_urban_canyon/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/challenger_2_m2_urban_canyon/DISPATCH.md, the worker handoff in /Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md, and AGENTS.md.
Empirically stress-test Melbourne-Wübbena shielding in mw.rs (verify that code multipath steps in shielded mode do NOT reset MwTrack, verify unshielded mode resets on > 1.0 cyc, verify genuine slips reset tracker).
Run tests and smoke benchmarks.
Deliver your verdict (APPROVE or REQUEST_CHANGES) in handoff.md and notify parent.
