# BRIEFING — 2026-09-25T21:13:15Z

## Mission
Refactor Gneiss RTK and GNSS/INS codebase to eliminate bare untyped vectors and floats, structurally enforcing compile-time frame safety, datum consistency, and temporal epoch alignment.

## 🔒 My Identity
- Archetype: teamwork_preview_orchestrator
- Roles: orchestrator, user_liaison, human_reporter, successor
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_orchestrator
- Original parent: parent
- Original parent conversation ID: 5d3514ba-c7b7-4d48-a991-062ca148f761

## 🔒 My Workflow
- **Pattern**: Project
- **Scope document**: /Users/kevin/projects/gneiss/PROJECT.md
1. **Decompose**: Decomposed into 5 Milestones (M1 Spatial Primitives, M2 Temporal Frames, M3 Relational Geometry, M4 Estimator Refactoring, M5 Verification) + E2E Testing Track.
2. **Dispatch & Execute**:
   - M1 Worker: Implementing zero-cost typestates and frame primitives in `gneiss-core::frames`.
   - M2 Worker: Implementing typed epoch systems, nanosecond time arithmetic, BeiDou alignment, and leap-second handling.
   - E2E Test Writer: Designing 4-tier requirement-driven opaque-box test suite (`TEST_INFRA.md`).
3. **On failure** (in this order):
   - Retry: nudge stuck agent or re-send task
   - Replace: spawn fresh agent with partial progress
   - Skip: proceed without (only if non-critical)
   - Redistribute: split stuck agent's remaining work
   - Redesign: re-partition decomposition
   - Escalate: report to parent (sub-orchestrators only, last resort)
4. **Succession**: at 16 spawns, write handoff.md, spawn successor
- **Work items**:
  1. Survey & Map Codebase [done]
  2. R1 Type-Safe Primitives (M1) [in-progress]
  3. R2 Temporal Frame & Epoch Safety (M2) [in-progress]
  4. R3 Estimator State & Pipeline Refactoring (M4) [pending]
  5. R4 Relational Coupling & Frame Safety (M3) [pending]
  6. E2E Dual-Track Testing & Benchmarking [in-progress]
- **Current phase**: 2
- **Current focus**: Parallel execution of M1 (Spatial), M2 (Temporal), and E2E Test Infrastructure

## 🔒 Key Constraints
- NEVER write, modify, or create source code files directly.
- NEVER run build/test commands yourself — require workers to do so.
- NEVER investigate or explore the problem at the code level — dispatch Explorers for technical investigation.
- File editing tools ONLY for metadata/state files (.md) in .agents/teamwork/ or PROJECT.md/ORIGINAL_REQUEST.md.
- Zero warnings, file size < 500 LOC, function size <= 32 LOC, nesting < 3, no unwrap() in production code.
- Pass all 789 tests + new tests + 2 CI smoke guards + eval_odaiba_ins benchmark.
- Never reuse a subagent after it has delivered its handoff — always spawn fresh.

## Current Parent
- Conversation ID: 5d3514ba-c7b7-4d48-a991-062ca148f761
- Updated: not yet

## Key Decisions Made
- Completed Survey phase with 3 parallel Explorers.
- Produced master `PROJECT.md` at root covering 21 features across 5 milestones.
- Dispatched Worker M1 (Spatial Primitives), Worker M2 (Temporal Frames), and E2E Test Writer in parallel with strict file ownership.

## Team Roster
| Agent | Type | Work Item | Status | Conv ID |
|-------|------|-----------|--------|---------|
| explorer_survey_1 | teamwork_preview_explorer | Survey R1 Spatial Frames & Datums | completed | ec60ae67-3513-47fc-b399-962f5f3707d3 |
| explorer_survey_2 | teamwork_preview_explorer | Survey R2 Temporal Frames & Epoch Safety | completed | 7448f16a-0734-4c8b-a268-c5ac6ec16581 |
| explorer_survey_3 | teamwork_preview_explorer | Survey R3 & R4 Estimators & Relational Coupling | completed | 3ded459c-7426-4c51-a0c8-dc1e1ab8341d |
| worker_m1 | teamwork_preview_worker | M1 Spatial Primitives Implementation | in-progress | aa368db6-99d7-4bcd-aada-9fb77b16d946 |
| worker_m2 | teamwork_preview_worker | M2 Temporal Frames Implementation | in-progress | 8dfe0774-dc20-4b00-9488-7405922f8bf2 |
| test_writer_e2e | teamwork_preview_test_writer | E2E Test Infra & Tiers 1-4 Test Suite | in-progress | 998fc75e-d3e6-4854-b2d3-b7aac2c291cc |

## Succession Status
- Succession required: no
- Spawn count: 6 / 16
- Pending subagents: aa368db6-99d7-4bcd-aada-9fb77b16d946, 8dfe0774-dc20-4b00-9488-7405922f8bf2, 998fc75e-d3e6-4854-b2d3-b7aac2c291cc
- Predecessor: none
- Successor: not yet spawned

## Active Timers
- Heartbeat cron: db66ae0c-b21b-4e14-ac97-93509c51c4b0/task-24
- Safety timer: managed by heartbeat cron
- On succession: kill all timers before spawning successor
- On context truncation: run manage_task(Action="list") — re-create if missing

## Artifact Index
- /Users/kevin/projects/gneiss/PROJECT.md — Master Project Specification & Roadmap
- /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md — Authoritative User Request
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_orchestrator/DISPATCH.md — Orchestrator Dispatch Record
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_orchestrator/BRIEFING.md — Persistent Working Memory
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_orchestrator/progress.md — Liveness Heartbeat & State Checkpoint
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_1/survey_r1_spatial.md — Survey 1 Full Report
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_2/survey_r2_temporal.md — Survey 2 Full Report
- /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3/survey_r3_r4_estimators.md — Survey 3 Full Report
