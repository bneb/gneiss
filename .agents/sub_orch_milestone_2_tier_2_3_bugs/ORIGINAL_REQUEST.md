# Original User Request

## 2026-06-22T01:00:13Z

You are the Milestone 2 Sub-Orchestrator for fixing remaining Tier 2 & 3 bugs in the gneiss navigation engine.
Your working directory is: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_2_tier_2_3_bugs
Your parent's conversation ID is: 2fa793b7-d67e-47b9-8b06-31cfa02fc26b

Scope:
Fix the following 8 bugs (Ranks 9-18):
- Bug 16: Mismatched Galileo BGD Correction (bgd_e1_e5b should be used for band 7 Galileo E5b observations, not bgd_e1_e5a)
- Bug 23: Klobuchar Model Evaluated at Receiver (evaluate at Ionospheric Pierce Point at ~350km altitude instead of receiver position)
- Bug 5: GMF Troposphere Longitude Omission (include longitude in the spherical harmonics annual evaluator)
- Bug 10: TOF Sat Position Clock Bias Neglect (account for receiver clock bias cdt_r in Time-of-Flight / signal transit time calculation)
- Bug 11: Sat PCV Zenith-Dependent Correction Omission (project and apply zenith-dependent PCV from ANTEX on line-of-sight)
- Bug 8: INS State AR Update Cutoff (stop zeroing out the Kalman gain for INS states in the Narrowlane constraint update)
- Bug 22: Saastamoinen Dry Delay Pressure Approx (use actual/standard surface pressure instead of a constant sea-level pressure scaled by height)
- Bug 3: L2C Phase Shift Bias (avoid hardcoding static 0.25 shift if receiver tracks in-phase, or handle properly)

Note: Bug 12 and Bug 25 are ALREADY fixed in the working tree. Verify them, but do not re-implement them.

For each bug, you must:
1. Decompose the milestone. Group these bugs or tackle them sequentially. Sequential execution is highly recommended to avoid merge conflicts and isolate test runs.
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
5. Ensure `cargo test --workspace` passes cleanly with no failures.
6. When all bugs are fixed and verified, update progress.md, write handoff.md, and send a completion message to the parent (conversation ID: 2fa793b7-d67e-47b9-8b06-31cfa02fc26b).

Let's begin!

## 2026-06-22T01:00:44Z

**Context**: Checking status of Milestone 2 Tier 2 & 3 Bugs sub-orchestration.
**Content**: Parent has resumed from a compaction/truncation checkpoint. We want to check your progress on Milestone 2.
**Action**: Please report your current status and progress. If you are active, update your progress.md and send a reply.

## 2026-06-22T01:16:25Z

**Context**: Milestone 2 Sub-Orchestration Worker monitoring.
**Content**: Worker 1 (conversation ID: b9500395-9eed-4576-83bf-8f7e99a812db) has not updated its progress.md since 01:04:40Z (approx 11 minutes).
**Action**: Please verify the liveness of Worker 1. If it has hung or is stalled, please execute the Escalation Ladder (Retry/Replace) according to the workflow instructions.

## 2026-06-22T01:20:09Z

**Context**: Milestone 2 Sub-Orchestration Worker monitoring.
**Content**: It is now 01:20:00Z (more than 3 minutes since you sent the status query to Worker 1 at 01:16:31Z). There is still no update in the worker's progress.md.
**Action**: Please run your liveness check. If Worker 1 has not responded, proceed with replacing the worker to avoid stalling the pipeline.

## 2026-06-22T06:00:45Z

You are the replacement Milestone 2 Sub-Orchestrator for fixing remaining Tier 2 & 3 bugs in the gneiss navigation engine.
Your working directory is: /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_2_tier_2_3_bugs
Your parent's conversation ID is: e2b4cf86-7ee9-4f3c-990c-2c79b1094647
Your predecessor's conversation ID was: c1e1438e-1aa8-4425-a0b8-6dc0b1d37f21

Scope:
Fix the following 8 bugs (Ranks 9-18):
- Bug 16: Mismatched Galileo BGD Correction (bgd_e1_e5b should be used for band 7 Galileo E5b observations, not bgd_e1_e5a)
- Bug 23: Klobuchar Model Evaluated at Receiver (evaluate at Ionospheric Pierce Point at ~350km altitude instead of receiver position)
- Bug 5: GMF Troposphere Longitude Omission (include longitude in the spherical harmonics annual evaluator)
- Bug 10: TOF Sat Position Clock Bias Neglect (account for receiver clock bias cdt_r in Time-of-Flight / signal transit time calculation)
- Bug 11: Sat PCV Zenith-Dependent Correction Omission (project and apply zenith-dependent PCV from ANTEX on line-of-sight)
- Bug 8: INS State AR Update Cutoff (stop zeroing out the Kalman gain for INS states in the Narrowlane constraint update)
- Bug 22: Saastamoinen Dry Delay Pressure Approx (use actual/standard surface pressure instead of a constant sea-level pressure scaled by height)
- Bug 3: L2C Phase Shift Bias (avoid hardcoding static 0.25 shift if receiver tracks in-phase, or handle properly)

Note: Bug 12 and Bug 25 are ALREADY fixed in the working tree. Verify them, but do not re-implement them.

Status:
Worker 1 Gen 2 (conversation ID: 0a8c14e1-9777-4c12-a796-d0ff576b8b8e) has completed the implementation of Bug 16 and written its handoff report to `/Users/kevin/projects/gneiss/.agents/worker_bug_16_gen2/handoff.md`.

For each bug, you must:
1. Recover state from the working directory /Users/kevin/projects/gneiss/.agents/sub_orch_milestone_2_tier_2_3_bugs. Update BRIEFING.md, SCOPE.md, and progress.md accordingly.
2. Start a fresh heartbeat cron.
3. For Bug 16, start directly at the Review and Auditor phase using Worker 1 Gen 2's output.
4. Run the iteration loop for each remaining bug:
   a. Spawn teamwork_preview_explorer to investigate the code and propose the fix strategy.
   b. Spawn teamwork_preview_worker to implement the fix and add a regression test.
      MANDATORY: include the integrity warning in the worker prompt:
      "DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected."
      MANDATORY: The regression test must pass now but would fail if the buggy code is restored.
   c. Spawn teamwork_preview_reviewer to verify build, tests, and code layout.
   d. Spawn teamwork_preview_auditor to run integrity verification (gating).
5. Ensure `cargo test --workspace` passes cleanly with no failures.
6. When all bugs are fixed and verified, update progress.md, write handoff.md, and send a completion message to the parent (conversation ID: e2b4cf86-7ee9-4f3c-990c-2c79b1094647).

Let's resume!

