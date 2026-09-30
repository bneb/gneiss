## 2026-09-25T21:13:05Z

You are Worker M1 for the Gneiss Frame Safety & Epoch Alignment Refactoring.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_worker_m1
Authoritative request file: /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md
Survey findings: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_1/survey_r1_spatial.md and /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_1/handoff.md
Project specification: /Users/kevin/projects/gneiss/PROJECT.md

MANDATORY: You MUST read /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md and /Users/kevin/projects/gneiss/AGENTS.md before starting work.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Scope & File Boundaries (M1 - Features 1–5):
You exclusively own:
- crates/gneiss-core/src/frames/ (all files, e.g. mod.rs, markers.rs, realizations.rs, primitives.rs, tangent.rs, positions.rs)
- crates/gneiss-core/src/coords/ (for local tangent plane integration if needed)
- Unit and compile-fail tests in crates/gneiss-core/tests/ (e.g. tests/test_frame_safety.rs)

Deliverables:
1. Formalize zero-cost typestate wrappers:
   - Markers: CoordinateFrame (Ecef<R>, Ned, Enu, BodyFrd)
   - Realizations: ReferenceFrame (ITRF2014, ITRF2020, WGS84, NAD83, JGD2011, and add PZ-90 Helmert aligned)
   - Primitives: SpatialVector<F>, SpatialVelocity<F>, SpatialCovariance<F>, Point3<F>, NedCovariance
   - AntennaLeverArm(SpatialVector<BodyFrd>), Attitude<From, To>
2. Disallow leaky Deref<Target=Vector3<f64>>; provide .coords() / .vector() methods.
3. Structurally prohibit cross-frame addition/subtraction without rotation at compile time.
4. Implement LocalTangentPlane<F> derived from a single EcefPos<F> to enforce relational coupling for ENU conversion.
5. Provide comprehensive unit tests and compile-fail tests verifying that invalid cross-frame arithmetic fails to compile.
6. Verify compliance with AGENTS.md:
   - File size strictly < 500 LOC
   - Function size strictly <= 32 LOC
   - Nesting depth strictly < 3
   - 0 compiler and clippy warnings under cargo clippy --workspace --all-targets -- -D warnings
   - 0 unwrap() in production code
   - All tests pass: cargo test --workspace
7. Write your handoff to handoff.md and send a completion message when done.
