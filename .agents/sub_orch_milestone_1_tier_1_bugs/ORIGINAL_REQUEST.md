# Original User Request

## 2026-06-21T03:14:37Z

You are the Milestone 1 Sub-Orchestrator for fixing 8 Tier 1 bugs in the gneiss navigation engine.
Your working directory is: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs
Your parent's conversation ID is: e2b4cf86-7ee9-4f3c-990c-2c79b1094647

Scope:
Fix all 8 Tier 1 bugs (Ranks 1-8):
- Bug 17: GLONASS Time Scale Discrepancy (in gneiss-core/src/ephemeris.rs or relevant GLONASS orbit code)
- Bug 1: Melbourne-Wübbena Dimensional Typo (in gneiss-rtk/src/engine/ppp_math.rs)
- Bug 9: Sequential AR Covariance Mismatch (in gneiss-rtk/src/engine/ppp_iekf.rs)
- Bug 2: Velocity-Attitude Transition Sign Mismatch (in gneiss-rtk/src/engine/predictor.rs)
- Bug 18: Opposite Sign in Phase Wind-Up Correction (in gneiss-rtk/src/engine/ppp.rs)
- Bug 15: Incorrect Broadcast Clock TGD Correction (in gneiss-core/src/ephemeris.rs)
- Bug 24: Outlier Tolerance in Precise Clock Gaps (in gneiss-parsers/src/rinex_clk.rs)
- Bug 6: GMF Legendre Unnormalized Polynomials (in gneiss-core/src/atmosphere.rs)

For each bug, you must:
1. Decompose the milestone. You may group these bugs or tackle them sequentially. Sequential execution is highly recommended to avoid merge conflicts and isolate test runs.
2. Initialize BRIEFING.md, SCOPE.md, and progress.md in your working directory.
3. Start a heartbeat cron.
4. Run the iteration loop:
   a. Spawn teamwork_preview_explorer to investigate the code and propose the fix strategy.
   b. Spawn teamwork_preview_worker to implement the fix and add a regression test.
      MANDATORY: include the integrity warning in the worker prompt:
      "DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected."
      MANDATORY: The regression test must pass now but would fail if the buggy code is restored.
   c. Spawn teamwork_preview_reviewer to verify build, tests, and code layout.
   d. Spawn teamwork_preview_auditor to run integrity verification (gating).
5. Ensure `cargo test --workspace` passes clean with no failures.
6. When all 8 bugs are fixed and verified, update progress.md, write handoff.md, and send a completion message to the parent (conversation ID: e2b4cf86-7ee9-4f3c-990c-2c79b1094647).

Let's begin!

## Follow-up — 2026-06-21T04:41:41Z

Resume work at /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_1_tier_1_bugs.
Read handoff.md, BRIEFING.md, ORIGINAL_REQUEST.md, and progress.md for current state.
Your parent is e2b4cf86-7ee9-4f3c-990c-2c79b1094647 — use this ID for all escalation and status reporting (send_message).
