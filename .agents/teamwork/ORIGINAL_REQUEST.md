# Original User Request

## 2026-09-25T21:01:47Z

Fan out sub agent teams to implement and review in parallel using a test-driven approach.

Refactor the Gneiss RTK and GNSS/INS codebase to eliminate bare, untyped vectors and floats, structurally enforcing compile-time frame safety, datum consistency, and temporal epoch alignment across all estimators and benchmark pipelines.

Working directory: /Users/kevin/projects/gneiss
Integrity mode: development

## Requirements

### R1. Type-Safe Coordinate, Vector, and Covariance Primitives
Formalize zero-cost typestate wrappers in `gneiss-core` (or `gneiss-rtk`) for spatial vectors and covariances, preventing unchecked math across coordinate frames (ECEF vs. NED vs. ENU vs. Body FRD) and geodetic datums (WGS84, ITRF2014, ITRF2020, NAD83, JGD2011, PZ-90). Disallow bare `Vector3<f64>` and raw float arrays across frame boundaries.

### R2. Temporal Frame & Epoch Alignment Safety
Enforce strictly typed epochs and time systems (`GpsTime`, `UtcTime`, `GlonassTime`, `BdtTime`). Prohibit ad-hoc float rounding or integer truncation for epoch matching. Structurally guard against leap-second offsets (e.g. GPS vs. UTC 18-second offset) between GNSS observables and IMU/base station streams.

### R3. Estimator State & Pipeline Refactoring (TDD)
Refactor `EskfState`, `SwfgEngine`, `PostProcessOptions`, and `eval_odaiba_ins` to consume typed frame primitives (e.g. `AntennaLeverArm` in Body FRD, `EcefPosition`, `EcefVelocity`, `NedCovariance`). Ensure attitude transformations $C_b^e$ are required to project body vectors to ECEF. All refactoring must follow a test-driven approach with pre-written regression and invariant tests.

### R4. Relational Coupling & Frame Safety Invariants
Structurally enforce relational coupling on shared geometry (e.g. double-difference line-of-sight unit vectors and satellite elevations derived from a single ephemeris state), ensuring geometric disagreement between rover and base is unrepresentable.

## Acceptance Criteria

### Frame Safety Invariants
- [ ] No bare `Vector3<f64>` or `[f64; 3]` used for spatial coordinates crossing frame or datum boundaries in refactored estimator modules.
- [ ] Attempting to perform vector math between mismatched frames (e.g., adding a `BodyFrameVector` directly to an `EcefPosition` without attitude rotation) fails to compile.
- [ ] Antenna lever arm is structurally typed in Body FRD and cannot be applied without explicit rotation.

### Test & Benchmark Verification
- [ ] All 789 existing workspace tests pass with 0 failures (`cargo test --workspace`).
- [ ] New unit tests in `gneiss-core` and `gneiss-rtk` specifically assert compile-time or runtime rejection of datum, frame, and epoch mismatches.
- [ ] Both CI smoke guards pass:
  - `python3 scripts/check_network_benchmark.py --smoke` outputs `ALL CHECKS PASSED`
  - `python3 scripts/check_multignss_benchmark.py --smoke` outputs `ALL CHECKS PASSED`
- [ ] `eval_odaiba_ins` reproduces the verified benchmark metrics ($p_{50} \le 1.80\text{ m}$, $\text{RMS} \le 3.50\text{ m}$, 0 false fixes).

### Code Standards Compliance (AGENTS.md)
- [ ] All modified/created files are strictly $< 500$ LOC.
- [ ] All functions are strictly $\le 32$ LOC with nesting depth $< 3$.
- [ ] Zero compiler and clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`).
- [ ] Zero `unwrap()` calls in production code.
