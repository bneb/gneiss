## 2026-09-25T21:04:10Z

You are Survey Explorer 1 for the Gneiss Frame Safety & Epoch Alignment Refactoring.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_1
Authoritative request file: /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md

You MUST read /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md and /Users/kevin/projects/gneiss/AGENTS.md first.

Focus: R1 (Type-Safe Coordinate, Vector, and Covariance Primitives) & Body FRD lever arm.
1. Survey crates/gneiss-core and crates/gneiss-rtk for all spatial coordinates, vectors, datums, frames, and covariances.
2. Search for existing representations: Vector3<f64>, [f64; 3], Matrix3<f64>, ECEF, NED, ENU, Body FRD, and geodetic datums (WGS84, ITRF2014, ITRF2020, NAD83, JGD2011, PZ-90).
3. Investigate how antenna lever arm in Body FRD is represented and used, and where rotation C_b^e is applied.
4. Design zero-cost typestate wrappers / phantom types preventing cross-frame / cross-datum math at compile time while maintaining high performance and ergonomics.
5. Identify all modules/files that will need refactoring or new primitives.
6. Check file sizes and function sizes per AGENTS.md constraints (< 500 LOC, <= 32 LOC).
7. Save your full findings in /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_1/survey_r1_spatial.md and write a complete self-contained handoff.md. Send a completion message when finished.
