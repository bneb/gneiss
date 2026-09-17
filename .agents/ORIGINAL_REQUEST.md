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

## Follow-up — 2026-06-21T04:52:58Z

**IMPORTANT — Bug 2 (Velocity-Attitude Jacobian Sign): Worker introduced a regression.**

The batch commit `da013e2` had already correctly fixed `predictor.rs` from `-f_e_skew * dt` → `+f_e_skew * dt`. The Bug 2 worker just re-introduced the negative sign, which is the **buggy** state.

I have manually restored `vel_att = f_e_skew * dt` (positive) in `predictor.rs`.

The correct derivation: the EKF perturbation on specific force is `δf_e = +[f_e×]ψ`, so `∂v/∂ψ = +[f_e×]` and the transition block must be `+f_e_skew * dt`.

Please ensure:
1. The auditor for Bug 2 checks that `vel_att = f_e_skew * dt` (NO leading minus sign)
2. The regression test `test_transition_matrix_velocity_attitude_coupling` asserts the **positive** coupling sign
3. Do not re-introduce the negative sign

Mark Bug 2 as fixed with the positive sign already present in the file.

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

## Follow-up — 2026-09-12T16:39:43Z

Use a very large team of agents. Implement all three Tier-1 commercial GNSS/INS frontiers in parallel across the Gneiss positioning engine codebase: (1) 15-state Error-State Kalman Filter (ESKF/MEKF) with closed-loop attitude and bias updates, (2) Integer PPP-AR engine with SINEX OSB phase bias ingestion, and (3) Network RTK Virtual Reference Station (VRS) spatial atmospheric engine, merging them into composite tightly-coupled modes.

Working directory: /Users/kevin/projects/gneiss
Integrity mode: development

## Requirements

### R1. 15-State Error-State Kalman Filter (ESKF/MEKF) for GNSS/INS
Expand the inertial filter state from 6-DOF ($p, v$) to a full 15-state Error-State Kalman Filter ($\delta \mathbf{p}^e, \delta \mathbf{v}^e, \delta \boldsymbol{\theta}, \delta \mathbf{b}_a, \delta \mathbf{b}_g$) with:
- Error-quaternion feedback to nominal attitude $\mathbf{q} \leftarrow \mathbf{q} \otimes \delta\mathbf{q}$.
- Closed-loop online accelerometer and gyroscope bias estimation driven by GNSS position and velocity innovations.
- Full 15-state backward Rauch-Tung-Striebel (RTS) smoother over forward filter history.
- Dynamic vehicle Non-Holonomic Constraints (NHC) and Zero-Velocity Updates (ZUPT) integrated into the 15-state covariance.

### R2. Integer PPP-AR Engine via SINEX OSB Ingestion
Implement an autonomous integer-fixing Precise Point Positioning (PPP-AR) engine:
- Ingest satellite Observation-Specific Bias (OSB) and fractional phase bias products from SINEX files (`.BIA` / `.OSB`).
- Formulate un-differenced carrier-phase and pseudorange observation equations with exact satellite and receiver phase center offsets/variations (PCO/PCV).
- Recover integer wide-lane (Melbourne-Wübbena) and narrow-lane ambiguities via LAMBDA search on single-difference or receiver-clock-decoupled ambiguities without physical base stations.
- Maintain continuous carrier tracking through satellite constellations (GPS, Galileo, BeiDou, QZSS).

### R3. Network RTK Virtual Reference Station (VRS) Atmospheric Engine
Build a multi-station regional CORS network atmospheric engine:
- Ingest 5–10 regional CORS base station observation streams simultaneously.
- Formulate multi-baseline double-difference network adjustment to solve for integer ambiguities across network baselines.
- Generate spatial 2D/3D Delaunay triangulation models for ionospheric delay pierce points and tropospheric zenith wet delay (ZWD) gradients.
- Synthesize localized Virtual Reference Station (VRS) observation data at the rover's approximate position, reducing effective baseline length to $< 1\text{ km}$ on 15–50 km regional networks.

### R4. Unified Composite Integration
Provide modular, clean interfaces that allow composing the 15-state ESKF with:
- Integer PPP-AR to provide Tightly-Coupled PPP/INS for base-station-free navigation.
- Network RTK VRS to provide Tightly-Coupled Network RTK/INS for metropolitan/survey navigation.

## Verification Resources

The codebase already contains benchmark evaluation harnesses and ground-truth datasets:
- `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`: Tokyo Odaiba 10Hz GNSS / 50Hz MEMS IMU urban canyon dataset with NovAtel SPAN reference truth.
- `crates/gneiss-rtk/src/bin/eval_ppp.rs`: RTK Explorer F9P vehicle dataset with Canadian Geodetic Service `rover_csrs.pos` and RTK truth.
- `crates/gneiss-rtk/src/bin/eval_network_ppk.rs`: Multi-station regional CORS baselines (P181, P224, P225, P222, SLAC, CAPO).
- `scripts/check_network_benchmark.py`: Smoke and full regression guard suite for network RTK.
- `scripts/check_multignss_benchmark.py`: Smoke and full regression guard suite for multi-constellation RTK.

## Acceptance Criteria

### Performance & Parity Targets
- [ ] **Odaiba GNSS/INS**: 15-state ESKF with RTS smoothing achieves $p_{50} < 2.5\text{ m}$ and RMS $< 5.2\text{ m}$ across the full 12,398-epoch 10Hz trajectory, improving on the 6-DOF baseline.
- [ ] **Kinematic PPP-AR**: Resolves integer ambiguities on the F9P kinematic drive, closing the discrepancy vs CSRS-PPP ($0.296\text{ m}$ RMS) to sub-meter kinematic accuracy.
- [ ] **Network RTK VRS**: Reduces baseline ppm error across 15–50 km CORS baselines (P181, P222, P225), achieving closer parity with Leica single-baseline specs ($8\text{ mm} + 1\text{ ppm})$.
- [ ] **Composite Modes**: Demonstrates successful execution of Tightly-Coupled PPP/INS and VRS-assisted RTK/INS pipelines.

### Code Quality & CI Invariants ([AGENTS.md](file:///Users/kevin/projects/gneiss/AGENTS.md))
- [ ] All code strictly adheres to file size $< 500$ LOC and function size $< 32$ LOC.
- [ ] Nesting depth strictly $< 3$ levels across all modified and new files.
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` passes with exactly 0 warnings.
- [ ] Exactly 0 `unwrap()` calls in production code (`match`, `if let`, `ok_or()?`, or descriptive `.expect()` only).
- [ ] All unit and integration tests (`cargo test --workspace`) pass with 0 failures.
- [ ] Both regression guard scripts (`check_network_benchmark.py --smoke` and `check_multignss_benchmark.py --smoke`) pass cleanly.
