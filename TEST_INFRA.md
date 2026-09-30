# Test Infrastructure: Frame Safety & Epoch Alignment Refactoring

**Document Version**: 1.0.0  
**Target Architecture**: Gneiss GNSS/INS Positioning Engine (`crates/gneiss-core`, `crates/gneiss-rtk`, `tests`)  
**Scope**: 4-Tier Opaque-Box E2E Integration Test Suite for R1–R4  
**Date**: 2026-09-25  

---

## 1. Test Philosophy & Overview

Precision satellite navigation and tightly coupled inertial navigation rely on strict mathematical consistency across coordinate reference frames, geodetic datums, temporal timescales, and sensor spatial baselines. When these boundaries are represented by bare floating-point numbers or raw `Vector3<f64>` primitives, subtle frame contamination bugs (e.g. adding a Body-frame antenna lever arm directly to an ECEF position, or mixing PZ-90 GLONASS orbits with WGS-84 without Helmert alignment) compile cleanly while corrupting centimeter-level positioning.

The Gneiss Frame Safety & Epoch Alignment test architecture enforces **opaque-box, requirement-driven verification**:
- **Opaque-Box Specification Conformance**: Tests interact exclusively through public APIs and interface contracts defined in `ORIGINAL_REQUEST.md` and `PROJECT.md`.
- **Authoritative Expected Outputs**: Expected values are derived from analytical mathematical invariants, authoritative geodetic standards (IERS conventions, WGS-84 ellipsoid parameters, Helmert 14-parameter models), official GNSS ICD specifications (GPS ICD-200, BDS SIS-ICD, GLONASS ICD), and rigid physical kinematics ($C_b^e \mathbf{l}^b$, $\boldsymbol{\omega} \times \mathbf{r}$).
- **Structural Invariant Verification**: Verifies that invalid physical operations (cross-frame addition, cross-datum math without transformation, cross-timescale subtraction) fail at compile-time or are structurally impossible to represent.
- **Relational Coupling Invariant Verification**: Parameters that are distinct views of the same physical object or epoch must be derived from shared typed inputs inside the callee, making geometric disagreement unrepresentable.

---

## 2. Four-Tier Test Suite Architecture

The test suite is structured into four distinct, self-contained verification tiers:

```
┌────────────────────────────────────────────────────────────────────────┐
│                        4-TIER TEST SUITE MATRIX                        │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 1: Feature Coverage (>=5 tests per feature, R1–R4)                │
│   • R1: Spatial primitives (Point3, SpatialVector, Covariance, Ned/Enu)│
│   • R1: Datum realizations (WGS84, ITRF2014, ITRF2020, NAD83, PZ-90)   │
│   • R1: Antenna lever arm & Attitude<From, To> C_b^e transformations   │
│   • R2: Temporal typestate scales (GPS, BDT, GST, GLONASS, UTC)        │
│   • R2: Integer nanosecond arithmetic, TimeDelta, and week rollover    │
│   • R2: EpochKey alignment and tolerance-based is_within matching      │
│   • R3: EskfState typed primitives and C_b^e lever arm projection      │
│   • R4: DoubleDiffGeometry relational coupling from single ephemeris   │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 2: Boundary & Corner Cases (>=5 tests per feature)                │
│   • R1: Antimeridian, polar singularities, null lever arm, zero rot    │
│   • R1: Extreme baseline lengths (0m zero-baseline to 12,000 km GEO)   │
│   • R2: Midnight week rollover (tow = 604,800.000s boundary)           │
│   • R2: Leap second insertion epochs (18s GPS-UTC, 14s GPS-BDT)        │
│   • R2: Sub-millisecond timing jitter and tolerance window edges       │
│   • R3: Extreme angular rate lever arm centrifugal acceleration        │
│   • R4: Zenith satellites (el -> 90°), horizon cutoff (el -> 0°)       │
│   • R4: Co-located base/rover and collinear satellite geometry         │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 3: Cross-Feature Pairwise Interactions                            │
│   • Spatial + Temporal: Ephemeris propagation with typed epoch & datum │
│   • Spatial + Geometry: Double-diff geometry with ITRF2014 base/rover  │
│   • Temporal + Estimator: ESKF state propagation across week rollover  │
│   • Lever Arm + Attitude: High-rate vehicle turn with Doppler velocity │
│   • Tangent Plane + Covariance: ECEF to NED/ENU covariance projection  │
│   • Multi-Constellation: GPS + BeiDou + GLONASS joint epoch alignment  │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 4: Real-World Application Scenarios (>=5 realistic missions)      │
│   • Scenario 1: Tokyo Odaiba Kinematic INS with Offset Antenna Lever   │
│   • Scenario 2: Multi-Station CORS Baseline Across Distinct Datums     │
│   • Scenario 3: Saturday Midnight Week-Rollover Continuous Kinematic   │
│   • Scenario 4: BeiDou B1I/B3I Clock-Bias-Free Baseline Solution       │
│   • Scenario 5: High-Dynamic UAV Pitch/Roll with Doppler & DD Updates  │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Feature Inventory & Invariant Mapping

| Feature # | Feature Description | Requirement | Target Module | Mathematical Invariant / Contract |
|---|---|---|---|---|
| **F-01** | Spatial typestates | R1 | `gneiss-core::frames` | `Point3<F> + SpatialVector<F> = Point3<F>`; cross-frame fails |
| **F-02** | Reference frames & PZ-90 | R1 | `gneiss-core::frames` | PZ-90 Helmert transform to ITRF2014; $\Delta \le 3\text{ mm}$ |
| **F-03** | Antenna lever arm & Attitude | R1 | `gneiss-core::frames` | $\mathbf{l}^e = C_b^e \mathbf{l}^b$; unrotated addition unrepresentable |
| **F-04** | Non-leaky typestate deref | R1 | `gneiss-core::frames` | No `Deref<Target=Vector3<f64>>`; explicit `.coords()` access |
| **F-05** | Relational tangent plane | R1 | `gneiss-core::frames` | `LocalTangentPlane<F>` from single origin; $C_{ecef}^{ned} C_{ned}^{ecef} = I$ |
| **F-06** | Temporal scale markers | R2 | `gneiss-core::time` | `Epoch<Scale>`; cross-scale subtraction fails to compile |
| **F-07** | Integer nanosecond time | R2 | `gneiss-core::time` | `tow_nanos \in [0, 604_800 \times 10^9)`; $t_2 - t_1 = \Delta t_{\text{nanos}}$ |
| **F-08** | BeiDou broadcast scale | R2 | `gneiss-core::time` | $\text{GPST} - \text{BDT} = 14.0\text{ s}$ exact; clock bias error $\to 0$ |
| **F-09** | Explicit leap seconds | R2 | `gneiss-core::time` | $\text{GPST} - \text{UTC} = 18.0\text{ s}$; explicit parameter required |
| **F-10** | EpochKey & alignment | R2 | `gneiss-core::time` | Monotonic continuous ms key; `is_within(dt)` tolerance check |
| **F-11** | Relational DD geometry | R4 | `gneiss-core::coords` | Derived from shared station & satellite positions; disagreement 0 |
| **F-12** | Receiver PCV coordinate fix | R4 | `gneiss-parsers` | Fix `az_el(rov_llh, rov_llh, ...)` passing LLH instead of ECEF |
| **F-13** | Base antenna coupling | R4 | `gneiss-parsers` | Base PCV evaluated at base elevation, not rover elevation |
| **F-14** | Typed ESKF state | R3 | `gneiss-rtk::eskf` | `pos_ecef: EcefPos<F>`, `vel_ecef: SpatialVelocity`, biases typed |
| **F-15** | Lever arm rotation in ESKF | R3 | `gneiss-rtk::eskf` | Innovation $\mathbf{y} = \mathbf{z} - (\mathbf{p}^e + C_b^e \mathbf{l}^b)$; Jacobian $H_\theta = -[\mathbf{l}^e \times]$ |
| **F-16** | SwfgEngine refactoring | R3 | `gneiss-rtk::swfg` | Initial & previous positions typed as `EcefPos<F>` |
| **F-17** | eval_odaiba_ins refactoring| R3 | `gneiss-rtk::bin` | Odaiba INS binary fully typed while remaining $< 500$ LOC |
| **F-18** | Workspace test verification| R1–R4 | `tests` | All existing workspace tests pass with 0 failures |
| **F-19** | CI smoke benchmarks | R1–R4 | `scripts` | `check_network_benchmark` and `check_multignss_benchmark` pass |
| **F-20** | Odaiba benchmark metric | R3 | `eval_odaiba_ins` | RTS horizontal $p_{50} \le 1.80\text{ m}$, $\text{RMS} \le 3.50\text{ m}$, 0 false fixes |
| **F-21** | AGENTS.md compliance | All | Entire codebase | $< 500$ LOC/file, $\le 32$ LOC/fn, nesting $< 3$, 0 warnings |

---

## 4. Test Suite Code Layout

To comply strictly with `AGENTS.md` (all files strictly $< 500$ LOC, all functions $\le 32$ LOC, nesting depth $< 3$ levels, 0 warnings under `clippy -D warnings`), the integration test suite in `tests/tests/` (and mirror `crates/gneiss-tests/tests/`) is structured into dedicated modular files:

```
tests/tests/
├── test_frame_safety_e2e.rs          # Runner entry point (< 100 LOC)
└── test_frame_safety_e2e/
    ├── common.rs                     # Geodetic fixtures, oracles, constants (< 350 LOC)
    ├── tier1_spatial.rs              # Tier 1: Spatial & datum typestate tests (< 450 LOC)
    ├── tier1_temporal.rs             # Tier 1: Temporal typestate & scale tests (< 450 LOC)
    ├── tier1_geometry.rs             # Tier 1: Relational geometry & PCV tests (< 400 LOC)
    ├── tier1_estimator.rs            # Tier 1: ESKF & estimator state tests (< 400 LOC)
    ├── tier2_boundaries.rs           # Tier 2: Boundary, rollover & singularity tests (< 450 LOC)
    ├── tier3_pairwise.rs             # Tier 3: Cross-feature pairwise interactions (< 450 LOC)
    └── tier4_scenarios.rs            # Tier 4: Real-world mission & benchmark scenarios (< 450 LOC)
```

---

## 5. Scenario Definitions (Tier 4)

### Scenario 1: Tokyo Odaiba Kinematic INS with Offset Antenna Lever Arm
- **Operational Context**: High-speed vehicle trajectory along Tokyo Odaiba coastal highway with significant roll/pitch dynamic variation and lever arm offset $\mathbf{l}^b = [0.25, 0.10, -0.85]^T$ m.
- **Verification Criteria**:
  - Lever arm rotated via instantaneous vehicle attitude $C_b^e$ before subtraction from GNSS antenna observations.
  - Angular rate coupling $\boldsymbol{\omega}_{ib}^b \times \mathbf{l}^b$ correctly incorporated into Doppler velocity innovations.
  - Unrotated lever arm addition fails at compile time.

### Scenario 2: Multi-Station CORS Baseline Across Distinct Datums
- **Operational Context**: Baseline processing between a national geodetic monument in ITRF2014 and a regional CORS station in JGD2011 or PZ-90.
- **Verification Criteria**:
  - Coordinate mixing across datums without explicit Helmert conversion is rejected at compile time.
  - Applying Helmert transformation matches analytical IERS transformation formulas within $< 0.1$ mm.
  - Local tangent plane conversion maintains orthogonal metric distance invariance.

### Scenario 3: Saturday Midnight Week-Rollover Continuous Kinematic
- **Operational Context**: Rover continuous tracking spanning GPS week boundary (Saturday 23:59:59 to Sunday 00:00:01, `tow` wrapping from $604{,}799\text{ s}$ to $1\text{ s}$).
- **Verification Criteria**:
  - TimeDelta across week rollover evaluates to exactly $+2.000000000$ seconds ($2 \times 10^9$ ns).
  - EpochKey generation is strictly monotonic across rollover without hash collision or key inversion.
  - `is_within` tolerance matching correctly identifies synchronized rover/base epochs across the rollover.

### Scenario 4: BeiDou B1I/B3I Clock-Bias-Free Baseline Solution
- **Operational Context**: Joint GPS + BeiDou multi-constellation RTK processing where BeiDou broadcast ephemeris `toe` and `toc` must be aligned.
- **Verification Criteria**:
  - Explicit 14-second BDT-to-GPST scale conversion applied infallibly.
  - Satellite clock bias evaluated at emission epoch without 14-second offset artifact.
  - Cross-timescale subtraction between `GpsTime` and `BdtTime` rejected at compile time.

### Scenario 5: High-Dynamic UAV Pitch/Roll with Doppler & DD Updates
- **Operational Context**: Aerobatic UAV undergoing rapid $\pm 45^\circ$ bank turns and pitch maneuvers while tracking double-difference carrier phase.
- **Verification Criteria**:
  - DoubleDiffGeometry computes elevation and azimuth independently for base and rover from shared satellite positions.
  - Attitude error Jacobian $H_\theta = (\mathbf{u}^s - \mathbf{u}^{ref})^T [\mathbf{l}^e \times]$ correctly projects lever arm cross product into carrier phase innovation.
  - Covariance projection from Body to ECEF preserves positive definiteness ($\lambda_{\min} > 0$).

---

## 6. Coverage Thresholds & Quality Invariants

- **Total Test Count**: $\ge 60$ comprehensive tests across Tiers 1–4.
- **Pass Rate**: 100% (zero failures, zero ignored).
- **Compiler & Clippy Warnings**: Exactly 0 under `cargo clippy --workspace --all-targets -- -D warnings`.
- **Modularity Limits**: All files $< 500$ LOC, all functions $\le 32$ LOC, nesting depth $< 3$ levels.
- **Production `unwrap()` Policy**: 0 `unwrap()` calls in production code.

---

## 7. Execution Commands

```bash
# Run the frame safety and epoch alignment integration test suite
cargo test -p gneiss_tests --test test_frame_safety_e2e

# Run with verbose output logging
cargo test -p gneiss_tests --test test_frame_safety_e2e -- --nocapture

# Run workspace clippy verification
cargo clippy --workspace --all-targets -- -D warnings
```
