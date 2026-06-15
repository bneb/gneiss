# plan.md — Gneiss Test Assertion Mistake Audit

This document outlines the decomposition, tracking, and aggregation strategy for auditing the test assertions across the `gneiss` workspace to identify assertion bugs or deficiencies that could hide actual bugs.

## Codebase Analysis Decomposition
We split the gneiss codebase and test suite into four main exploration scopes:

| Explorer | Target Scope | Description |
|---|---|---|
| **Explorer 1** | `crates/gneiss-core`, `tests/` | Coords, time, sat, signal, tides, dop, windup, and workspace-level integration tests |
| **Explorer 2** | `crates/gneiss-rtk` (Part A) | Ambiguity resolution, calibration, and filter utilities |
| **Explorer 3** | `crates/gneiss-rtk` (Part B) | EKF engine, updates, predictors, measurements, and jacobians |
| **Explorer 4** | `crates/gneiss-parsers`, `gneiss-fetch`, `gneiss-geodesy`, `gneiss-ntrip` | RINEX, RTCM3, SP3, UBX parsers, HTTP fetchers, geodesy, and NTRIP client |

## Verification & Audit Strategy
Explorer subagents will perform read-only static analysis and keyword searches to locate:
1. **Trivial/Always-True Assertions**: e.g., `assert!(true)`, `assert_eq!(x, x)`, `assert_ne!(x, y)` where `x == y` is statically true or comparing identical variables.
2. **Commented-out Assertions**: Assertions that were disabled or commented out, leaving tests silent.
3. **No-op / Silent Tests**: Tests that invoke logic but discard output without checking anything (unless clearly designed as a panic test).
4. **Excessive Tolerance in Approximations**: e.g., `assert_approx_eq!` or `approx` with oversized epsilon values (e.g. > 0.1 for high-precision GNSS calculations).
5. **Logic Bugs in Assertions**: e.g., using `||` instead of `&&` or vice versa, or checking `is_ok()` on errors instead of checking the expected values.
6. **Mismatched / Typo Constants**: Hardcoded values in assertions that do not correspond to the actual expected physics/math (e.g., coordinates, constants).

## Deliverables
- Individual Explorer reports saved in their respective directories: `.agents/explorer_<N>/handoff.md`.
- A consolidated synthesis document drafted in `.agents/orchestrator/suspicious_tests_report.md`.
- A final markdown file copied to `/Users/kevin/projects/gneiss/suspicious_tests_report.md` via a worker agent.

## Execution Schedule
1. **Stage 1**: Dispatch Explorer 1, 2, 3, 4.
2. **Stage 2**: Monitor and collect results.
3. **Stage 3**: Synthesize findings and write consolidated report.
4. **Stage 4**: Verify and deliver the final report.
