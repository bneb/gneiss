# Challenger 1 Task Assignment — Milestone 2: CMC Multipath Mitigation Stress Testing

## Role & Mission
You are Challenger 1 (`teamwork_preview_challenger`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/challenger_1_m2_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Input Documents
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, Section R2)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications & interface contracts)
3. `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md` (Worker M2 handoff report)
4. `/Users/kevin/projects/gneiss/AGENTS.md` (code standards)

## Challenge Scope
Empirically stress-test the CMC multipath screening and down-weighting mechanisms:
1. Verify decoupled screening in `screen.rs`: Inject +30m and -50m pseudorange gross errors on pairs with active carrier phase. Confirm carrier phase is retained with nominal variance, code variance is inflated to 1e8, and pairs without carrier phase are removed.
2. Verify CMC tracking in `robust.rs`: Feed continuous arcs with injected 5m, 10m, and 20m code multipath steps. Confirm `CmcTrack` detects the step, freezes the baseline, and flags multipath. Confirm cycle slips reset the baseline cleanly.
3. Run:
```bash
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```

Write your verdict (APPROVE or REQUEST_CHANGES) and full report to `/Users/kevin/projects/gneiss/.agents/challenger_1_m2_urban_canyon/handoff.md`.
Then send a completion message to your parent.

## 2026-09-24T19:27:47Z
You are Challenger 1 for Milestone 2: CMC Multipath Mitigation.
Your working directory is: /Users/kevin/projects/gneiss/.agents/challenger_1_m2_urban_canyon/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/challenger_1_m2_urban_canyon/DISPATCH.md, the worker handoff in /Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md, and AGENTS.md.
Empirically stress-test CMC multipath detection and decoupled gross PR screening (inject 5m, 10m, 20m code steps, verify carrier phase retention, verify baseline freezing and reset on slips).
Run tests and smoke benchmarks.
Deliver your verdict (APPROVE or REQUEST_CHANGES) in handoff.md and notify parent.

