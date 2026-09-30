# BRIEFING — 2026-09-25T14:13:05-07:00

## Mission
Implement Milestone 1 (Features 1-5) of the Gneiss Frame Safety & Epoch Alignment Refactoring: zero-cost typestate wrappers, reference frames, primitives, relational local tangent planes, no leaky derefs, compile-fail tests, and full AGENTS.md compliance.

## 🔒 My Identity
- Archetype: Worker
- Roles: implementer, qa, specialist
- Working directory: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_worker_m1
- Original parent: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Milestone: M1 (Features 1-5)

## 🔒 Key Constraints
- File boundaries: crates/gneiss-core/src/frames/ (all files), crates/gneiss-core/src/coords/ (if needed), tests/ (e.g. tests/test_frame_safety.rs)
- Integrity mandate: genuine implementation only, no dummy/facade, no hardcoded values.
- AGENTS.md rules:
  - File size strictly < 500 LOC
  - Function size strictly <= 32 LOC
  - Nesting depth strictly < 3
  - 0 compiler and clippy warnings under cargo clippy --workspace --all-targets -- -D warnings
  - 0 unwrap() in production code
  - All workspace tests pass: cargo test --workspace

## Current Parent
- Conversation ID: db66ae0c-b21b-4e14-ac97-93509c51c4b0
- Updated: not yet

## Task Summary
- **What to build**: Coordinate frames (Ecef<R>, Ned, Enu, BodyFrd), ReferenceFrame realizations (ITRF2014, ITRF2020, WGS84, NAD83, JGD2011, PZ90 Helmert aligned), Primitives (SpatialVector, SpatialVelocity, SpatialCovariance, Point3, NedCovariance, AntennaLeverArm, Attitude), LocalTangentPlane derived from single EcefPos, disallow leaky Deref, compile-fail tests.
- **Success criteria**: All types strictly typed, no cross-frame arithmetic allowed without rotation, comprehensive unit and compile-fail tests, 0 warnings, clippy clean, AGENTS.md compliant.
- **Interface contracts**: PROJECT.md, ORIGINAL_REQUEST.md
- **Code layout**: crates/gneiss-core/src/frames/

## Key Decisions Made
- Starting investigation of existing frames/ and coords/ code, survey findings, and project specs.

## Artifact Index
- DISPATCH.md — Assignment from orchestrator
- BRIEFING.md — Situational awareness working memory

## Change Tracker
- **Files modified**: None yet
- **Build status**: Not run yet
- **Pending issues**: None

## Quality Status
- **Build/test result**: Not run yet
- **Lint status**: Not run yet
- **Tests added/modified**: None yet

## Loaded Skills
None
