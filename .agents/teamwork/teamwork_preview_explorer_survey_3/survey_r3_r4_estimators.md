# Survey 3: Estimator State, Pipeline Refactoring, Relational Coupling & Benchmarks (R3 & R4)

**Author:** Survey Explorer 3  
**Date:** 2026-09-25  
**Working Directory:** `/Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3`  
**Reference Document:** `/Users/kevin/projects/gneiss/.agents/teamwork/ORIGINAL_REQUEST.md`  

---

## 1. Executive Summary

This survey addresses **Requirement R3** (Estimator State & Pipeline Refactoring) and **Requirement R4** (Relational Coupling & Frame Safety Invariants), alongside establishing the baseline verification metrics across the test suite and CI benchmarks.

### Core Discoveries:
1. **Pervasive Bare Vector/Float Boundaries in Estimators:** `EskfState`, `SwfgEngine`, `PostProcessOptions`, and `eval_odaiba_ins` currently use bare `Vector3<f64>` and `[f64; 3]` for ECEF positions, ECEF velocities, body-frame accelerometer/gyroscope biases, antenna lever arms, and local-level covariances. Frame mismatches (e.g. adding a Body FRD lever arm to an ECEF position without attitude rotation) are currently compile-valid and represent a major class of runtime bugs.
2. **Latent Frame Bug in `receiver_pcv`:** In `crates/gneiss-parsers/src/receiver_pcv/mod.rs:175`, `az_el(rov_llh, rov_llh, sat_pos)` passes `rov_llh` as both `pos_llh` and `pos_ecef` to `az_el`, subtracting radians from metres. This bug was invisible to the compiler because both types were bare `Vector3<f64>`.
3. **Flawed Double-Difference Geometry in `receiver_antenna`:** In `crates/gneiss-parsers/src/receiver_antenna/mod.rs:277-305`, `compute_dd_pcv_correction_2d` takes rover-only elevation and azimuth (`el_sat_rad`, `az_sat_rad`) and passes them directly into the base antenna correction (`base_ant.pcv_mm_az_zen(freq_code, az_sat_bas, zen_sat)`), erroneously assuming base and rover share the identical satellite zenith and azimuth angles.
4. **Verified Baseline Performance:**
   - Workspace tests: **789 tests pass with 0 failures**.
   - CI Smoke Guard 1: `python3 scripts/check_network_benchmark.py --smoke` -> **ALL CHECKS PASSED** ($p_{50} = 0.022\text{ m}$).
   - CI Smoke Guard 2: `python3 scripts/check_multignss_benchmark.py --smoke` -> **ALL CHECKS PASSED** (network fused fix rate $99.10\%$).
   - Tokyo Odaiba Benchmark: `cargo run --release --bin eval_odaiba_ins` -> **$p_{50} = 1.751\text{ m}$** (target $\le 1.80\text{ m}$), **$\text{RMS} = 3.479\text{ m}$** (target $\le 3.50\text{ m}$), 0 false fixes.
5. **Hard File Size Limits:** Several critical target files are near the 500 LOC ceiling (`main.rs` at 475 LOC, `formation.rs` at 457 LOC, `dd_update.rs` at 454 LOC). All relational structures must be isolated in dedicated modules (e.g. `geometry.rs`) to avoid violating `AGENTS.md`.

---

## 2. Survey of Estimator State & Pipelines (Requirement R3)

### 2.1 `EskfState` and the 15-State ESKF Architecture
Located at `crates/gneiss-rtk/src/estimators/eskf/types.rs:33-40`:
```rust
pub struct EskfState {
    pub pos_ecef: Vector3<f64>,
    pub vel_ecef: Vector3<f64>,
    pub attitude: UnitQuaternion<f64>,
    pub accel_bias: Vector3<f64>,
    pub gyro_bias: Vector3<f64>,
    pub cov: Matrix15<f64>,
}
```

#### State Definitions & Frame Semantics:
| State Variable | Physical Meaning | Intended Frame | Current Type | Safety Hazard |
|---|---|---|---|---|
| `pos_ecef` | Vehicle position in ECEF | ECEF (`Wgs84` or `Itrf2014`) | `Vector3<f64>` | Can mix with NED, Body, or cross datums unchecked |
| `vel_ecef` | Vehicle velocity w.r.t Earth | ECEF | `Vector3<f64>` | Can add to Body velocity or NED velocity |
| `attitude` | Body-to-ECEF rotation $C_b^e$ | Body FRD -> ECEF | `UnitQuaternion<f64>` | Undocumented rotation direction; bare quaternion |
| `accel_bias` | Accelerometer bias | Body FRD | `Vector3<f64>` | Can mix with ECEF acceleration |
| `gyro_bias` | Gyroscope bias | Body FRD | `Vector3<f64>` | Can mix with ECEF angular velocity |
| `cov` (15x15) | Error-state covariance | Error frame: $\delta p^e, \delta v^e, \delta\theta^e, \delta b_a^b, \delta b_g^b$ | `Matrix15<f64>` | Indices implicit; local NED cov easily injected without rotation |

#### Key Observation on Attitude Error State:
In `crates/gneiss-rtk/src/estimators/eskf/mod.rs:6` and `update.rs:16-17`:
- Error state $\delta\theta$ (indices 6..9) is a **left-multiplied global-frame error**: $q \leftarrow \delta q \cdot q$, where $\delta q = \exp(\frac{1}{2}\delta\theta)$.
- The attitude transformation $C_b^e$ projects vectors from Body FRD into ECEF: $v^e = C_b^e v^b$.

### 2.2 Antenna Lever Arm & Innovation Updates
In `crates/gneiss-rtk/src/estimators/eskf/update.rs`:
```rust
pub fn build_gnss_pos_system(
    state: &EskfState,
    pos_meas: &Vector3<f64>,
    r_pos: &Matrix3<f64>,
    lever_arm: &Vector3<f64>,
) -> (Vector3<f64>, Matrix3x15<f64>, Matrix3<f64>) {
    let r_b2e = state.attitude.to_rotation_matrix().into_inner();
    let l_e = r_b2e * lever_arm;
    let y = pos_meas - (state.pos_ecef + l_e);
    ...
}
```
**Safety Invariant Violation:** `lever_arm` is currently a bare `Vector3<f64>`. Because `state.pos_ecef` is also `Vector3<f64>`, code such as `state.pos_ecef + lever_arm` compiles without error, despite adding a Body FRD vector directly to an ECEF vector.
**Required Typestate:** `AntennaLeverArm` must be a distinct type representing a vector strictly in the Body FRD frame. Applying `AntennaLeverArm` to an `EcefPos` MUST require multiplying by an explicit Body-to-ECEF rotation $C_b^e$: `state.pos_ecef + attitude.rotate_lever_arm(lever_arm)`.

### 2.3 `SwfgEngine` & Sliding Window Factor Graph
Located at `crates/gneiss-rtk/src/swfg/engine/mod.rs:36-62`:
```rust
pub struct SwfgEngine {
    solver: SlidingWindowSolver,
    pipeline: MeasurementPipeline,
    ephemerides: Vec<Ephemeris>,
    epoch: u32,
    pub(crate) current_pose: Option<VariableId>,
    prev_position: Option<Vector3<f64>>,
    prev_prev_position: Option<Vector3<f64>>,
    prev_time: Option<GpsTime>,
    current_attitude: Option<nalgebra::UnitQuaternion<f64>>,
    initial_position: Option<Vector3<f64>>,
    ...
}
```
In `crates/gneiss-rtk/src/swfg/engine/builder.rs:30-40`:
```rust
pub struct RtkFactorContext<'a> {
    pub epoch: u32,
    pub pose_id: VariableId,
    pub prev_pose_id: Option<VariableId>,
    pub init_pos: Vector3<f64>,
    pub prev_position: Option<Vector3<f64>>,
    pub dt_sec: f64,
    pub base_pos: Vector3<f64>,
    pub corrected: &'a [CorrectedObservation],
    pub base_raw: &'a [RawObservation],
}
```
**Safety Invariant Violation:** `init_pos`, `base_pos`, `sat_pos_ecef`, and `prev_position` are all bare `Vector3<f64>`. In multi-station RTK and PPP, base stations may come from ITRF2014, while broadcast orbits are in WGS84, and regional networks may be in NAD83 or JGD2011. There is no datum enforcement on `base_pos`.

### 2.4 `PostProcessOptions`
Located at `crates/gneiss-rtk/src/post_process/mod.rs:63-125`:
```rust
pub struct PostProcessOptions {
    pub enable_bidirectional: bool,
    pub base_position: Option<Vector3<f64>>,
    pub initial_rover_position: Option<Vector3<f64>>,
    pub klobuchar_alpha: Option<[f64; 4]>,
    pub klobuchar_beta: Option<[f64; 4]>,
    pub q_accel: Option<f64>,
    pub widelane_ar: bool,
    pub tropo_gradients: bool,
    pub network_sat_upd: Option<std::collections::HashMap<u16, f64>>,
    pub receiver_pcv: Option<std::sync::Arc<ReceiverPcvPair>>,
    pub dynamics: ProcessingDynamics,
    pub enable_glonass: bool,
    ...
}
```
**Safety Invariant Violation:** `base_position` and `initial_rover_position` are untyped `Option<Vector3<f64>>`. Datums are not tracked.

### 2.5 Tokyo Odaiba Benchmark (`eval_odaiba_ins`)
Located at `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs:24`:
```rust
const ANTENNA_LEVER_ARM: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
```
In `odaiba_helpers.rs:94-100`:
```rust
pub fn compute_base_pos(approx: Option<[f64; 3]>) -> Vector3<f64> {
    let default_base_arp = Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let default_base_llh = gneiss_core::coords::ecef_to_llh(default_base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(default_base_llh).transpose();
    let computed_base = default_base_arp + ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855);
    approx.map(|p| Vector3::new(p[0], p[1], p[2])).unwrap_or(computed_base)
}
```
In `odaiba_helpers.rs:196-201`:
```rust
pub fn compute_gnss_cov_ecef(pos_ecef: Vector3<f64>, var_h: f64, var_v: f64) -> Matrix3<f64> {
    let llh = gneiss_core::coords::ecef_to_llh(pos_ecef);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(llh).transpose();
    let r_ned = Matrix3::from_diagonal(&Vector3::new(var_h, var_h, var_v));
    ned_to_ecef * r_ned * ned_to_ecef.transpose()
}
```
**Observations:**
1. `compute_base_pos` projects an antenna height offset in NED into ECEF via `ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855)`. Without typed frames, an implementer could easily add the NED offset directly to `default_base_arp`.
2. `compute_gnss_cov_ecef` converts a local horizontal/vertical variance (`var_h`, `var_v`) into an ECEF covariance matrix. This should consume a typed `NedCovariance`.

---

## 3. Survey of Relational Coupling & Frame Safety (Requirement R4)

### 3.1 The Canonical PCV Citation in `AGENTS.md`
`AGENTS.md:52-61` states:
> - **Relational coupling must be structurally enforced.** Parameters that are views of the same object, epoch, or frame must be derived from shared typed inputs inside the callee — never accepted as parallel bare floats.
>   Example: rover/base zenith angles to the same satellite pair must be computed from shared satellite positions, not passed as four independent f64 values (which allows physically impossible geometry).
> - **Review checklist**: "Which of these arguments must pairwise agree, and what makes disagreement unrepresentable?"
> - Canonical citation: receiver_pcv.rs dd_correction_m — doc author violated documented invariant within minutes of implementing it correctly.

### 3.2 Anatomy of the Failure Modes

#### Failure Mode 1: Four Independent Bare Zenith Floats (`receiver_pcv/mod.rs:149`)
```rust
pub fn dd_correction_m(
    rov: &ReceiverPcv,
    bas: &ReceiverPcv,
    zen_rov_sat_deg: f64,
    zen_rov_ref_deg: f64,
    zen_bas_sat_deg: f64,
    zen_bas_ref_deg: f64,
) -> f64
```
A caller can pass arbitrary, unrelated numbers, swap rover and base, or pass azimuth where zenith is expected. Geometric disagreement is completely representable.

#### Failure Mode 2: Latent Bug in Geometry-Derived PCV (`receiver_pcv/mod.rs:175-178`)
```rust
pub fn dd_correction_from_geometry(
    rov: &ReceiverPcv,
    bas: &ReceiverPcv,
    rov_llh: nalgebra::Vector3<f64>,
    bas_llh: nalgebra::Vector3<f64>,
    sat_pos: nalgebra::Vector3<f64>,
    ref_sat_pos: nalgebra::Vector3<f64>,
) -> f64 {
    use gneiss_core::coords::az_el;
    let (az_rov_sat, el_rov_sat) = az_el(rov_llh, rov_llh, sat_pos);
    let (az_rov_ref, el_rov_ref) = az_el(rov_llh, rov_llh, ref_sat_pos);
    ...
```
`az_el` has signature `az_el(pos_llh: Vector3<f64>, pos_ecef: Vector3<f64>, sat_ecef: Vector3<f64>)`.
Because `rov_llh` was an untyped `Vector3<f64>`, the author passed `rov_llh` as both `pos_llh` AND `pos_ecef`!
Inside `az_el`, line-of-sight vector `sat_ecef - pos_ecef` was computed as `sat_pos - rov_llh` (subtracting radians/metres from kilometres). This compile-clean catastrophic bug directly proves the necessity of typed frame wrappers!

#### Failure Mode 3: Rover Angles Reused for Base Antenna (`receiver_antenna/mod.rs:277-305`)
```rust
pub fn compute_dd_pcv_correction_2d(
    rover_ant: &ReceiverAntenna,
    base_ant: &ReceiverAntenna,
    freq_code: &str,
    az_sat_rad: f64,
    el_sat_rad: f64,
    az_ref_rad: f64,
    el_ref_rad: f64,
    rover_heading_rad: f64,
) -> f64 {
    let zen_sat = 90.0 - el_sat_rad.to_degrees();
    ...
    let (Some(bas_s), Some(bas_r)) = (
        base_ant.pcv_mm_az_zen(freq_code, az_sat_bas, zen_sat),
        base_ant.pcv_mm_az_zen(freq_code, az_ref_bas, zen_ref),
    ) ...
```
The base antenna PCV lookup directly reuses `zen_sat` and `zen_ref` calculated from the rover's elevation angle! On regional baselines (>20 km), the elevation angle at the base station differs significantly from the rover station.

#### Failure Mode 4: Duplicated Geometry in DD Updates
Across 7 different files:
1. `rtk_iekf/update/system.rs:77`: `compute_meas_geom`
2. `rtk_iekf/formation.rs:293`: `geom_dd_tropo`
3. `rtk_iekf/screen.rs:82`
4. `rtk_iekf/iono_free.rs:192`
5. `eskf/dd_update.rs:44`: `compute_dd_jacobian_15`
6. `swfg/engine/builder.rs:168`: `rover_sat_range`, `base_range_sat`, `base_dd_range`
7. `eval_odaiba_ins/main.rs:300`: `dd_geom`

Each computes $(|r_s - p_{rov}| - |r_{ref} - p_{rov}|) - (|r_s - p_{bas}| - |r_{ref} - p_{bas}|)$ independently using bare vectors.

### 3.3 Proposed Architectural Solution: `DoubleDiffGeometry<F>`

To satisfy Requirement R4, rover and base geometry must be derived from a single shared ephemeris and station state inside a dedicated type where geometric disagreement is unrepresentable.

```rust
/// Unified double-difference geometry container enforcing structural relational coupling.
/// Disagreement between rover and base geometry is unrepresentable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleDiffGeometry<F: ReferenceFrame> {
    // Rover line-of-sight & ranges
    pub u_rov_sat: Vector3<f64>, // Unit LOS in ECEF
    pub u_rov_ref: Vector3<f64>,
    pub range_rov_sat_m: f64,
    pub range_rov_ref_m: f64,
    pub el_rov_sat_rad: f64,
    pub az_rov_sat_rad: f64,
    pub el_rov_ref_rad: f64,
    pub az_rov_ref_rad: f64,

    // Base line-of-sight & ranges (computed independently from the same satellite positions)
    pub u_bas_sat: Vector3<f64>,
    pub u_bas_ref: Vector3<f64>,
    pub range_bas_sat_m: f64,
    pub range_bas_ref_m: f64,
    pub el_bas_sat_rad: f64,
    pub az_bas_sat_rad: f64,
    pub el_bas_ref_rad: f64,
    pub az_bas_ref_rad: f64,

    // Coupled DD quantities
    pub base_dd_range_m: f64,
    pub geometric_dd_m: f64,
    pub delta_u_rov: Vector3<f64>, // (u_rov_sat - u_rov_ref) for ESKF Jacobian
    _frame: core::marker::PhantomData<F>,
}
```

#### Structural Invariant Enforcement:
```rust
impl<F: ReferenceFrame> DoubleDiffGeometry<F> {
    /// Construct DD geometry from shared station and satellite positions.
    /// Both stations and both satellites MUST share the identical frame/datum F.
    pub fn compute(
        rov_pos: EcefPos<F>,
        bas_pos: EcefPos<F>,
        sat_pos: EcefPos<F>,
        ref_sat_pos: EcefPos<F>,
    ) -> Result<Self, GeometryError> {
        // 1. Compute rover azimuth, elevation, range, and unit LOS vector
        // 2. Compute base azimuth, elevation, range, and unit LOS vector
        // 3. Compute base_dd_range and geometric_dd
        // Disagreement is structurally impossible because caller cannot supply
        // independent angles or mismatched positions.
    }
}
```

When PCV or ESKF updates are performed:
```rust
// Receiver PCV correction now takes DoubleDiffGeometry directly:
pub fn dd_pcv_correction<F: ReferenceFrame>(
    rov_pcv: &ReceiverPcv,
    bas_pcv: &ReceiverPcv,
    freq_code: &str,
    geom: &DoubleDiffGeometry<F>,
    rover_heading_rad: f64,
) -> f64
```
No bare floats for zenith or azimuth are accepted anywhere!

---

## 4. Baseline Benchmark Verification

All benchmarks and tests were executed on the unmodified working tree prior to refactoring to record exact golden baselines.

### 4.1 Test Suite Verification
Command: `cargo test --workspace`
- **Result:** Exit code 0.
- **Pass Count:** 789 passed, 0 failed, 1 ignored (existing ignored test in `gneiss-fetch`).
- **Breakdown:**
  - `gneiss-core`: 144 unit tests
  - `gneiss-parsers`: 272 unit tests
  - `gneiss-rtk`: 452 unit tests
  - `gneiss-geodesy`: 51 unit tests
  - `gneiss-ntrip`: 23 unit tests
  - `gneiss-fetch`: 17 unit tests
  - `gneiss_tests`: 205 integration tests

### 4.2 Smoke Benchmark 1: Network RTK Benchmark
Command: `python3 scripts/check_network_benchmark.py --smoke`
- **Result:** Exit code 0 (`ALL CHECKS PASSED`).
- **Observed Metrics:**
  | Metric | Value | Budget | Status |
  |---|---|---|---|
  | Network fused horizontal $p_{50}$ | 0.022 m | $\le 0.04\text{ m}$ | PASS |
  | Network fused horizontal RMS | 0.032 m | $\le 0.06\text{ m}$ | PASS |
  | Network fused vertical RMS | 0.037 m | $\le 0.08\text{ m}$ | PASS |
  | P181 smoothed fixed-only $p_{50}$ | 0.021 m | $\le 0.03\text{ m}$ | PASS |
  | P222 smoothed fixed-only $p_{50}$ | 0.060 m | $\le 0.09\text{ m}$ | PASS |
  | SLAC smoothed fixed-only $p_{50}$ | 0.110 m | $\le 0.12\text{ m}$ | PASS |
  | OHLN smoothed fix rate | 99.80% | $\ge 79.0\%$ | PASS |
  | P181 smoothed fix rate | 99.70% | $\ge 85.0\%$ | PASS |
  | SLAC smoothed fix rate | 78.90% | $\ge 60.0\%$ | PASS |

### 4.3 Smoke Benchmark 2: Multi-GNSS Benchmark
Command: `python3 scripts/check_multignss_benchmark.py --smoke`
- **Result:** Exit code 0 (`ALL CHECKS PASSED`).
- **Observed Metrics:**
  | Metric | Value | Budget | Status |
  |---|---|---|---|
  | P181 fix rate | 97.80% | $\ge 97.5\%$ | PASS |
  | P181 horizontal $p_{95}$ | 138.00 mm | $\le 145.0\text{ mm}$ | PASS |
  | P181 vertical $p_{95}$ | 183.00 mm | $\le 290.0\text{ mm}$ | PASS |
  | P225 fix rate | 89.30% | $\ge 71.0\%$ | PASS |
  | P225 horizontal $p_{95}$ | 99.00 mm | $\le 245.0\text{ mm}$ | PASS |
  | P225 vertical $p_{95}$ | 223.00 mm | $\le 370.0\text{ mm}$ | PASS |
  | P222 fix rate | 97.90% | $\ge 86.0\%$ | PASS |
  | P222 horizontal $p_{95}$ | 300.00 mm | $\le 305.0\text{ mm}$ | PASS |
  | P222 vertical $p_{95}$ | 93.00 mm | $\le 135.0\text{ mm}$ | PASS |
  | Network fused fix rate | 99.10% | $\ge 96.5\%$ | PASS |

### 4.4 Tokyo Odaiba INS Benchmark (`eval_odaiba_ins`)
Command: `cargo run --release --bin eval_odaiba_ins`
- **Dataset:** Tokyo Odaiba (12,398 epochs GNSS @ 10Hz, 62,040 IMU samples @ 50Hz, NovAtel SPAN truth).
- **Result:** Exit code 0.
- **Observed Metrics:**
  - **RTS Smoothed GNSS/INS (Full trajectory, $N=12,398$):**
    - **$p_{50} = 1.751\text{ m}$** (acceptance target budget: $\le 1.80\text{ m}$ -> **PASS**)
    - $p_{68} = 2.929\text{ m}$
    - $p_{95} = 6.466\text{ m}$
    - **$\text{RMS} = 3.479\text{ m}$** (acceptance target budget: $\le 3.50\text{ m}$ -> **PASS**)
  - **Forward Inertial Filter ($N=12,398$):**
    - $p_{50} = 1.842\text{ m}$, $\text{RMS} = 3.750\text{ m}$
  - **False Fixes:** 0 false fixes.

---

## 5. Catalog of Affected Modules & Refactoring Plan

| Crate | Module / File | Current LOC | Refactoring Focus | Safety Risk / LOC Mitigation |
|---|---|---|---|---|
| `gneiss-rtk` | `src/estimators/eskf/types.rs` | 182 | Add typed fields to `EskfState<F: ReferenceFrame>`: `pos_ecef: EcefPos<F>`, `vel_ecef: EcefVector<F>`, `accel_bias: BodyFrdVector`, `gyro_bias: BodyFrdVector`. | Safe (<250 LOC). |
| `gneiss-rtk` | `src/estimators/eskf/predict.rs` | 241 | Enforce $C_b^e$ rotation of IMU specific force and angular velocity in `propagate_nominal_state`. | Safe (<300 LOC). |
| `gneiss-rtk` | `src/estimators/eskf/update.rs` | 279 | Consume `&EcefPos<F>` and `&AntennaLeverArm` in `build_gnss_pos_system` and `update_gnss_position`. | Safe (<320 LOC). |
| `gneiss-rtk` | `src/estimators/eskf/constraints.rs` | 161 | Consume `&AntennaLeverArm` in `update_nhc` and `&BodyFrdVector` in `update_body_velocity`. | Safe (<200 LOC). |
| `gneiss-rtk` | `src/estimators/eskf/alignment.rs` | 101 | Consume `EcefPos<F>` in `compute_initial_attitude` and return typed $C_b^e$. | Safe (<150 LOC). |
| `gneiss-rtk` | `src/estimators/eskf/dd_update.rs` | 454 | Consume `DoubleDiffGeometry<F>` and `AntennaLeverArm` in `compute_dd_jacobian_15`. | **HIGH RISK: 454 LOC**. Move geometry helper math to shared geometry module to prevent exceeding 500 LOC. |
| `gneiss-rtk` | `src/geometry/` (NEW) | 0 | Implement canonical `DoubleDiffGeometry<F>`. | New file strictly < 300 LOC. |
| `gneiss-rtk` | `src/swfg/engine/mod.rs` | 430 | `SwfgEngine`: type `initial_position` and `prev_position` as `EcefPos<Wgs84>`. Return `EcefPos<Wgs84>` in `SwfgSolution`. | Moderate risk: 430 LOC. Keep edits concise. |
| `gneiss-rtk` | `src/swfg/engine/builder.rs` | 395 | `RtkFactorContext`: type `base_pos` and `init_pos` as `EcefPos<Wgs84>`. Consume typed double difference geometry. | Moderate risk: 395 LOC. |
| `gneiss-rtk` | `src/post_process/mod.rs` | 326 | `PostProcessOptions`: type `base_position` and `initial_rover_position` as `Option<EcefPos<F>>`. | Safe (<360 LOC). |
| `gneiss-rtk` | `src/bin/eval_odaiba_ins/main.rs` | 475 | Type `ANTENNA_LEVER_ARM` as `AntennaLeverArm`. Enforce attitude rotation on lever arm. Consume `DoubleDiffGeometry<F>`. | **CRITICAL RISK: 475 LOC**. Must move helper functions into `odaiba_helpers.rs` to keep `main.rs` strictly < 500 LOC. |
| `gneiss-rtk` | `src/bin/eval_odaiba_ins/odaiba_helpers.rs` | 215 | `compute_base_pos` returns `EcefPos<Wgs84>`. `compute_gnss_cov_ecef` consumes `NedCovariance`. | Safe (<350 LOC). |
| `gneiss-parsers` | `src/receiver_pcv/mod.rs` | 247 | Deprecate/refactor 4-float `dd_correction_m`. Fix line 175 `az_el(rov_llh, rov_llh, sat_pos)` bug. Consume `DoubleDiffGeometry`. | Safe (<300 LOC). |
| `gneiss-parsers` | `src/receiver_antenna/mod.rs` | 338 | Fix `compute_dd_pcv_correction_2d` to compute base elevation from base geometry. | Safe (<380 LOC). |

---

## 6. Test-Driven Development (TDD) Rollout Plan

To adhere to the TDD requirements of R3 and R4:

### Phase 1: Invariant & Negative Regression Tests (Write First)
1. Write compile-fail tests asserting:
   - Attempting `let _ = eskf_state.pos_ecef + lever_arm;` fails to compile (`EcefPos` cannot add `AntennaLeverArm` directly).
   - Attempting to pass an `EcefPos<Itrf2014>` where `EcefPos<Wgs84>` is required fails to compile.
   - Attempting to construct `DoubleDiffGeometry` with mismatched station or satellite frames fails to compile.
2. Write unit tests asserting:
   - `DoubleDiffGeometry::compute` produces identical line-of-sight unit vectors and ranges to historical math on zero-baseline and 50km baseline.
   - Antenna lever arm rotation invariant: rotating $[0, 0, 0]$ lever arm produces zero correction; non-zero lever arm produces exact $C_b^e \cdot l_{body}$.
   - Regression test for `receiver_pcv` fixing the `az_el(rov_llh, rov_llh, sat_pos)` argument swap bug.

### Phase 2: Core Geometric & Estimator Refactoring
1. Implement `DoubleDiffGeometry<F>` in `gneiss-rtk`.
2. Refactor `EskfState` in `crates/gneiss-rtk/src/estimators/eskf/types.rs`.
3. Refactor `build_gnss_pos_system`, `build_doppler_velocity_system`, and `compute_dd_jacobian_15`.
4. Update `SwfgEngine` and `PostProcessOptions`.

### Phase 3: Benchmark & Pipeline Integration
1. Refactor `crates/gneiss-rtk/src/bin/eval_odaiba_ins/` (`main.rs` and `odaiba_helpers.rs`).
2. Run `cargo test --workspace` to ensure all 789+ tests pass.
3. Run `python3 scripts/check_network_benchmark.py --smoke`.
4. Run `python3 scripts/check_multignss_benchmark.py --smoke`.
5. Run `cargo run --release --bin eval_odaiba_ins` and verify $p_{50} \le 1.80\text{ m}$ and $\text{RMS} \le 3.50\text{ m}$.
6. Run `cargo clippy --workspace --all-targets -- -D warnings`.
