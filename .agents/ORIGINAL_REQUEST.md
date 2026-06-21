# Original User Request

## Initial Request — 2026-06-15T18:14:31Z

# Teamwork Project Prompt — Draft

> Status: Launched
> Goal: Craft prompt → get user approval → delegate to teamwork_preview

Federate a search through the `gneiss` workspace's test suite to find mistakes in test assertions that might be hiding actual bugs.

Working directory: /Users/kevin/projects/gneiss
Integrity mode: demo

## Requirements

### R1. Broad Workspace Scan
Scan the test suite across the entire `gneiss` workspace. Identify test assertions that contain logical errors, trivially pass when they shouldn't, shadow underlying bugs, or misuse mathematical/domain constraints.

### R2. Non-Destructive Reporting
Do not modify the test or production codebase. Focus entirely on analysis and generate a comprehensive markdown report of the uncovered assertion mistakes. 

## Acceptance Criteria

### Reporting Quality
- [ ] The team outputs a report named `suspicious_tests_report.md`.
- [ ] Every flagged test in the report explicitly cites the absolute file path and line number of the suspicious assertion.
- [ ] Every flagged test includes a brief, concrete explanation of why the assertion is logically flawed and what bug it might be obscuring.

## Follow-up — 2026-06-21T03:13:25Z

Fix 25 mathematically-identified bugs in the `gneiss` GNSS/PPP navigation engine (a Rust codebase), and guard each fix with a regression test. The bugs span EKF dynamics, ambiguity resolution, satellite geometry, atmospheric modeling, and geophysical corrections — some are sign errors, some are omissions, and some are outright wrong formulas. The goal is to close the accuracy gap versus RTKLIB.

Working directory: /Users/kevin/projects/gneiss

## Bug Reference

All 25 bugs are documented with detailed analyses, correct formulas, and file/line references in:
`/Users/kevin/.gemini/antigravity/brain/3e07e73a-4b87-4801-b363-5d6f67bdb076/analysis_results.md`

They are stack-ranked into four tiers:
- **Tier 1 (Critical)**: Ranks 1–8 — mathematically verified, high-impact errors (fix first)
- **Tier 2 (High/Conditional)**: Ranks 9–13 — large errors under specific conditions
- **Tier 3 (Medium/Systematic)**: Ranks 14–18 — centimeter-level systematic biases
- **Tier 4 (Enhancement)**: Ranks 19–25 — modeling gaps and missing features

## Requirements

### R1. Fix all Tier 1 bugs (highest priority)
Fix each of the 8 Tier 1 bugs exactly as described in the analysis document. The correct formula or fix is specified in each bug's analysis section. Tackle these first before moving to lower tiers.

### R2. Fix Tier 2 and Tier 3 bugs
After completing Tier 1, fix each of the 10 bugs in Tiers 2 and 3 in priority order (lowest rank first). These bugs involve correcting existing implementations; refactoring surrounding code is allowed and encouraged if it improves correctness or testability.

### R3. Implement Tier 4 features
After completing Tiers 1–3, implement Tier 4 items fully, even if that requires new data structures, file parsing (e.g., BLQ files for OTL), or new modules. If a full implementation is not feasible, add a well-documented `// TODO:` stub and a `#[test] #[ignore]` that documents the expected behavior.

### R4. Regression tests for every fix
For every bug fixed in R1–R3, write a Rust unit test in the same crate (in the file or its `#[cfg(test)]` module) that:
- Exercises the corrected code path with concrete inputs
- Asserts the mathematically correct result (exact or within a tight epsilon)
- **Would fail if the buggy formula were restored** — no tautological assertions
- Is named to reflect the bug it guards (e.g., `test_mw_combination_dimensionless`, `test_windup_sign_correct`, `test_glonass_orbit_time_scale`)

If a test for the bug already exists, verify its assertions against the relevant spec (IERS Conventions, ICD, or algorithm reference) and correct the test if needed.

### R5. All existing tests continue to pass
After all changes, `cargo test --workspace` must pass with zero failures. Do not remove or weaken existing test assertions.

## Acceptance Criteria

### Tier 1 Bugs Fixed (verifiable by code inspection)
- [ ] Melbourne-Wübbena: factor is `(lam2 - lam1) / (lam1 + lam2)` (dimensionless), not `(lam1 * lam2) / (lam1 + lam2)` (has units of meters)
- [ ] GLONASS orbit integration: time delta `dt` is computed consistently in one time scale (GPST or GLONASST), not mixed
- [ ] Sequential AR: Narrowlane innovation covariance computed from post-Widelane-update state covariance
- [ ] Velocity-attitude Jacobian: sign of `f_e_skew` term in the transition matrix is positive (`+f_e_skew`)
- [ ] Phase wind-up: `wup` is subtracted from the carrier phase (not added) before forming residuals
- [ ] Broadcast TGD: not subtracted when the observation is dual-frequency or ionosphere-free
- [ ] Precise clock gap: function returns `None` (not stale bias) when gap exceeds the tolerance threshold
- [ ] GMF: associated Legendre functions are fully normalized

### Regression Tests
- [ ] Every fix in R1–R3 has at least one `#[test]` that would fail with the original buggy code
- [ ] `cargo test --workspace` passes with zero failures after all changes

### Build Integrity
- [ ] `cargo build --workspace` succeeds with zero errors and zero new warnings
- [ ] No `#[allow(...)]` suppressions introduced without a code comment justifying them
