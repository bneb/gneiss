## 2026-09-25T21:13:05Z

You are the E2E Test Writer for the Gneiss Frame Safety & Epoch Alignment Refactoring.
Your working directory is: /Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_test_writer_e2e
Authoritative request file: /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md
Project specification: /Users/kevin/projects/gneiss/PROJECT.md

MANDATORY: You MUST read /Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md and /Users/kevin/projects/gneiss/AGENTS.md before starting work.

Scope: E2E Testing Track (Requirement-Driven, Opaque-Box).
1. Create /Users/kevin/projects/gneiss/TEST_INFRA.md following the template from project patterns:
   - Test Philosophy: Opaque-box, requirement-driven, derived from ORIGINAL_REQUEST.md
   - Feature Inventory: All features from PROJECT.md
   - 4-Tier Test Architecture:
     * Tier 1: Feature Coverage (>=5 tests per feature)
     * Tier 2: Boundary & Corner Cases (>=5 tests per feature)
     * Tier 3: Cross-Feature Combinations (pairwise interactions)
     * Tier 4: Real-World Application Scenarios (>=5 realistic scenarios)
2. Implement test suites in crates/gneiss-tests/tests/ (e.g. test_frame_safety_e2e.rs) covering:
   - Coordinate, vector, datum compile-time safety and boundary tests
   - Temporal epoch alignment, leap-second offsets, and week rollover tests
   - Relational coupling invariants on double-difference geometry
   - Estimator state safety invariants (lever arm rotation, attitude coupling)
3. Ensure all tests compile and pass alongside existing workspace tests:
   cargo test -p gneiss_tests
4. Create /Users/kevin/projects/gneiss/TEST_READY.md when test infrastructure and tests are ready.
5. Verify compliance with AGENTS.md:
   - File size < 500 LOC (split across multiple test files if needed)
   - Function size <= 32 LOC
   - Nesting < 3
   - 0 compiler warnings
6. Write your handoff to handoff.md and send a completion message when done.
