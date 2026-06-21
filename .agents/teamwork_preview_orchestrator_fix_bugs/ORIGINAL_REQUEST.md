# Original User Request

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

## 2026-06-21T03:13:58Z

You are the Project Orchestrator for fixing 25 mathematically-identified bugs in the gneiss navigation engine.
Your coordination directory is: /Users/kevin/projects/gneiss/.agents/teamwork_preview_orchestrator_fix_bugs
Please initialize your files (such as `plan.md` and `progress.md`) in that directory. Do not write source code or tests into your coordination directory; all implementation changes and tests must be in the gneiss source crates.
Your input files:
- Global user request: /Users/kevin/projects/gneiss/.agents/teamwork_preview_orchestrator_fix_bugs/ORIGINAL_REQUEST.md
- Bug Analysis Reference: /Users/kevin/.gemini/antigravity/brain/3e07e73a-4b87-4801-b363-5d6f67bdb076/analysis_results.md

Requirements:
1. Fix all 8 Tier 1 bugs first.
2. Fix all 10 Tier 2 and 3 bugs in priority order.
3. Implement all 7 Tier 4 features (if a full implementation is not feasible, add a well-documented stub and an ignored regression test).
4. Guard each fix with a regression test that passes but would have failed with the buggy code.
5. Ensure `cargo test --workspace` and `cargo build --workspace` succeed cleanly.

Please write your plan to plan.md, start spawning worker/explorer agents as needed, update progress.md continuously, and report back when all milestones are complete.

## Follow-up — 2026-06-21T09:57:18Z

**IMPORTANT STATUS UPDATE — Before starting any work, read this.**

The following bugs have ALREADY been fixed in the working tree (uncommitted) since the quota reset. Do NOT re-implement them — check the code first and verify they are present:

**Already fixed (verify, don't redo):**
- Bug 18 (Phase Wind-Up Sign): `ppp.rs` — `(cp1 - wup)` in all 4 call sites. Test `test_windup_sign_correct` added.
- Bug 6 (GMF Legendre Normalization): `atmosphere.rs` — `_legendre_norm` function added. Tests `test_legendre_normalization` and `test_gmf_longitude_variation` added.
- Bug 12 (Receiver PCV): `ppp.rs` — `compute_receiver_pcv` function added and called in `build_sats`.
- Bug 25 (Covariance Inflation on Slip): `ppp.rs` — `covariance[(i,i)] *= 4.0` after ambiguity removal. Test `test_covariance_inflated_on_slip` added.
- Galileo BGD struct: `filter.rs` — `bgd_e1_e5b: 0.0` field added.

**NOT yet fixed — please continue with these in priority order:**
1. Bug 15: TGD not applied for dual-frequency (ephemeris.rs / calc_keplerian)
2. Bug 24: Stale clock gap returns None (rinex_clk.rs)
3. Bug 16: Galileo BGD band selection (E5a vs E5b)
4. Bug 23: Klobuchar IPP at 350km altitude (atmosphere.rs)
5. Bug 5: GMF longitude term in _sh_eval_annual (atmosphere.rs)
6. Remaining Tier 3/4 bugs

**Constraint:** `vel_att = f_e_skew * dt` in predictor.rs must remain POSITIVE. Do not change it.

Run `cargo test --workspace` before starting to confirm current state. All tests should pass.
