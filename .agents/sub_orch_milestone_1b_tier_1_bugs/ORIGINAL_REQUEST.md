# Original User Request

## Initial Request — 2026-06-21T02:57:17-07:00

You are the Milestone 1b Sub-Orchestrator. Your role is to orchestrator the fixes for the remaining Tier 1 bugs in the gneiss navigation engine.

Working Directory: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1b_tier_1_bugs
Parent Conversation ID: 2fa793b7-d67e-47b9-8b06-31cfa02fc26b

Your scope is to resolve the following 4 bugs in priority order:
1. Bug 18: Opposite Sign in Phase Wind-Up Correction
   - File: crates/gneiss-rtk/src/engine/ppp.rs (around lines 499, etc.)
   - Bug: observed phase is corrected by adding instead of subtracting wup. wup must be subtracted.
2. Bug 15: Incorrect Broadcast Clock TGD Correction
   - File: crates/gneiss-core/src/ephemeris.rs (around lines 403, etc.)
   - Bug: tgd is subtracted for dual-frequency/iono-free clock corrections. Under these modes, TGD cancels, so it shouldn't be subtracted.
3. Bug 24: Outlier Tolerance in Precise Clock Gaps
   - File: crates/gneiss-parsers/src/rinex_clk.rs (around lines 75, etc.)
   - Bug: returning Some(r1.bias) instead of None when t - r1.time exceeds 900.0 seconds.
4. Bug 6: GMF Legendre Unnormalized Polynomials
   - File: crates/gneiss-core/src/atmosphere.rs (around lines 178, etc.)
   - Bug: standard unnormalized Legendre polynomials are used instead of fully normalized Associated Legendre Functions.

For each bug, you must:
1. Decompose the bug and plan the fix.
2. Invoke an Explorer (teamwork_preview_explorer) to analyze the relevant file, locate the bug, and recommend a strategy.
3. Invoke a Worker (teamwork_preview_worker) to implement the fix and add a regression test. 
   - MANDATORY: You must include this verbatim warning in the Worker's prompt:
     "DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected."
4. Invoke 2 Reviewers (teamwork_preview_reviewer) to independently verify the changes and regression tests.
5. Invoke a Forensic Auditor (teamwork_preview_auditor) to perform integrity verification. Note: the Forensic Auditor verdict must be CLEAN.
6. The regression test must pass with the fix, but would fail if the buggy formula were restored. It must be named to reflect the bug (e.g. test_windup_sign_correct).
7. Ensure cargo test and cargo build pass cleanly.

Reference detailed bug analyses in:
/Users/kevin/.gemini/antigravity/brain/3e07e73a-4b87-4801-b363-5d6f67bdb076/analysis_results.md

Initialize your BRIEFING.md, progress.md, and SCOPE.md in your working directory. Track your heartbeat and update progress.md regularly. Report back to your parent conversation ID with a handoff report when all 4 bugs are resolved, verified, and audited CLEAN.
