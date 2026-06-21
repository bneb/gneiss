# Soft Handoff: Milestone 1 Sub-Orchestrator Successor Handoff

## Milestone State
- **Bug 17: GLONASS Time Scale Discrepancy**: DONE (Fixed, verified, and audited CLEAN)
- **Bug 1: Melbourne-Wübbena Dimensional Typo**: DONE (Verified, and audited CLEAN)
- **Bug 9: Sequential AR Covariance Mismatch**: DONE (Fixed, verified, and audited CLEAN)
- **Bug 2: Velocity-Attitude Transition Sign Mismatch**: IN_PROGRESS (Final Worker `b1f32300-1121-4d28-9edf-2bee98aa0260` has successfully implemented the correct negative sign fix, updated the unit test, formatted the test file, and verified workspace tests. Reviewers and Auditor must now run on this finalized fix.)
- **Bug 18: Opposite Sign in Phase Wind-Up Correction**: PLANNED (Note: Previous worker failed; needs to be restarted after Bug 2 is finalized)
- **Bug 15: Incorrect Broadcast Clock TGD Correction**: PLANNED
- **Bug 24: Outlier Tolerance in Precise Clock Gaps**: PLANNED
- **Bug 6: GMF Legendre Unnormalized Polynomials**: PLANNED

## Active Subagents
None currently active.

## Pending Decisions
None.

## Remaining Work
The successor needs to resume execution at Bug 2's verification phase:
1. Spawn 2 Reviewer subagents to verify the finalized Bug 2 correction (`vel_att = -f_e_skew * dt;` in `predictor.rs` and the matching negative sign unit test in `tests_predictor.rs`).
2. Spawn the Auditor subagent to perform forensic integrity auditing on the final Bug 2 fix.
3. Once Bug 2 is audited CLEAN, mark Bug 2 as DONE in `progress.md` and `SCOPE.md`.
4. Proceed to Bug 18, spawning Explorer -> Worker -> Reviewers -> Auditor to fix the wind-up sign discrepancy.
5. Decompose and implement the remaining bugs (15, 24, 6) sequentially.

## Key Artifacts
- **progress.md**: `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/progress.md`
- **BRIEFING.md**: `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/BRIEFING.md`
- **SCOPE.md**: `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/SCOPE.md`
- **ORIGINAL_REQUEST.md**: `/Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs/ORIGINAL_REQUEST.md`
