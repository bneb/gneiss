# Handoff Report: Survey 1 — Spatial Primitives & Body FRD Lever Arm (R1)

## 1. Observation

Direct code observations from inspection of `crates/gneiss-core` and `crates/gneiss-rtk`:

1. **Untyped Coordinate Primitives and Leaky Deref**:
   - `crates/gneiss-core/src/coords.rs:54–69`: `Coordinate` uses runtime enums (`Datum` and `Frame`). Runtime checking `ensure_aligned()` (line 72) is easily bypassed because callers extract `vector: Vector3<f64>`.
   - `crates/gneiss-core/src/frames/positions.rs:9`: `pub struct EcefPos<F: ReferenceFrame>(pub Vector3<f64>, pub PhantomData<F>)`.
   - `crates/gneiss-core/src/frames/positions.rs:43–48`: Implements `core::ops::Deref<Target = Vector3<f64>>`, stripping frame encapsulation and allowing bare vector operations.
   - `crates/gneiss-core/src/frames/realizations.rs:29–157`: Realizations exist for `Itrf2014`, `Itrf2020`, `Igs20`, `Wgs84Broadcast`, `Nad83_2011`, `Etrs89`, `Gda2020`, `Jgd2011`.
   - `Pz90` is missing from `realizations.rs`, despite `ephemeris/glonass.rs:77–104` returning PZ-90 positions as bare `Vector3<f64>`.

2. **Antenna Lever Arm and Frame Coupling in ESKF**:
   - `crates/gneiss-rtk/src/estimators/eskf/types.rs:33–40`:
     `EskfState` holds `pos_ecef: Vector3<f64>`, `vel_ecef: Vector3<f64>`, `attitude: UnitQuaternion<f64>`, `accel_bias: Vector3<f64>`, `gyro_bias: Vector3<f64>`. All fields are bare.
   - `crates/gneiss-rtk/src/estimators/eskf/dd_update.rs:9–14`:
     ```rust
     pub struct DdSatGeometry {
         pub u_s: Vector3<f64>,
         pub u_ref: Vector3<f64>,
         pub lever_arm: Vector3<f64>,
     }
     ```
     `u_s` (ECEF unit LOS), `u_ref` (ECEF unit LOS), and `lever_arm` (Body FRD) are stored together as bare `Vector3<f64>`.
   - `crates/gneiss-rtk/src/estimators/eskf/dd_update.rs:52`:
     `let l_e = state.attitude.to_rotation_matrix() * geom.lever_arm;`
   - `crates/gneiss-rtk/src/estimators/eskf/update.rs:48–50`:
     `let r_b2e = state.attitude.to_rotation_matrix().into_inner();`
     `let l_e = r_b2e * lever_arm;`
     `let y = pos_meas - (state.pos_ecef + l_e);`
   - `crates/gneiss-rtk/src/estimators/eskf/update.rs:92–97`:
     `let v_rot_b = omega_corr.cross(lever_arm);`
     `let v_rot_e = r_b2e * v_rot_b;`
     `let v_ant_pred = state.vel_ecef + v_rot_e;`
   - In all cases, unrotated addition `pos_meas - (state.pos_ecef + lever_arm)` or `state.vel_ecef + v_rot_b` would compile without compiler diagnostics.

3. **Relational Coupling Invariant Violations**:
   - `crates/gneiss-core/src/coords.rs:255–259`:
     `pub fn ecef_delta_to_enu(target_ecef: Vector3<f64>, ref_ecef: Vector3<f64>, ref_llh: Vector3<f64>) -> Vector3<f64>`
     Accepts parallel bare vectors `ref_ecef` and `ref_llh` representing two separate representations of the same reference origin.
   - `crates/gneiss-core/src/coords.rs:177`:
     `pub fn az_el(pos_llh: Vector3<f64>, pos_ecef: Vector3<f64>, sat_ecef: Vector3<f64>) -> (f64, f64)`
     Accepts parallel `pos_llh` and `pos_ecef`.

4. **Factor Graph Lever Arm Omission**:
   - `crates/gneiss-rtk/src/swfg/config.rs:153`: `pub lever_arm: [f64; 3]`.
   - `crates/gneiss-rtk/src/swfg/pipeline/dd_factors.rs:65–68`:
     `let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);`
     `let rover_sat_range = (self.sat_pos - rx_pos).norm();`
     The factor ignores the lever arm and vehicle attitude, assuming the antenna is colocated with the IMU center of navigation.

5. **Code Standard Thresholds (AGENTS.md)**:
   - `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs`: 474 LOC (26 LOC under the 500 LOC ceiling).
   - `crates/gneiss-rtk/src/estimators/eskf/dd_update.rs`: 453 LOC.
   - `crates/gneiss-rtk/src/estimators/eskf/condition.rs`: 452 LOC.
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/state/mod.rs`: 448 LOC.

---

## 2. Logic Chain

1. **Premise 1 (From Observation 1 & 2)**: Both positions, velocities, orientations, and lever arms are represented as bare `Vector3<f64>` or `UnitQuaternion<f64>`.
2. **Premise 2 (From Observation 2)**: The antenna lever arm $\mathbf{l}_b$ physically exists in the vehicle Body FRD frame, while receiver positions $\mathbf{r}_{imu}^e$ and satellite line-of-sight unit vectors $\mathbf{u}^s$ physically exist in the ECEF frame.
3. **Inference 1**: Because Rust allows addition between any two `Vector3<f64>`, a developer can directly compute $\mathbf{r}_{imu}^e + \mathbf{l}_b$ without applying $C_b^e$. The compiler cannot catch this error.
4. **Premise 3 (From Observation 1)**: `EcefPos<F>` implements `Deref<Target = Vector3<f64>>` and lacks sibling types for NED, ENU, and Body FRD. It also lacks typestates for velocities and covariances.
5. **Inference 2**: The existing `EcefPos` type in `frames/positions.rs` is incomplete for estimator refactoring and leaks its inner vector through `Deref`.
6. **Premise 4 (From Observation 3)**: `AGENTS.md` mandates that relational coupling must be structurally enforced. Parameters that are views of the same object must be derived from shared typed inputs.
7. **Inference 3**: `ecef_delta_to_enu` and `az_el` violate the relational coupling standard by taking parallel bare floats. A dedicated `LocalTangentPlane<F>` struct constructed from a single `EcefPos<F>` is required.
8. **Inference 4 (From Observation 5)**: `eval_odaiba_ins/main.rs` is at 474 LOC. Any refactoring to consume typed primitives in `main.rs` must be disciplined and helper-delegated to avoid breaching 500 LOC.

---

## 3. Caveats

1. **Legacy Compatibility**: Many test files throughout `tests/` and benchmarks in `bin/` consume bare `Vector3<f64>`. Transitioning them should use non-breaking wrappers with convenient `.coords()`, `.vector()`, or checked conversions rather than breaking all 789 workspace tests at once.
2. **Epoch Alignment Integration**: Frame typestates (R1) interact closely with Temporal Epoch systems (R2). `Point3<F>` can optionally carry or be bundled with `GpsTime` in `EpochPosition<F, R>`, but must remain zero-cost when manipulated in inner loops.
3. **PZ-90 Realization Parameters**: Broadcast GLONASS ephemeris is referenced to PZ-90.11. The Helmert transformation parameters adopted into `Pz90` are the standard ITRF2014 alignment parameters ($\approx 3$ mm shift).

---

## 4. Conclusion

The codebase currently lacks compile-time frame safety, allowing illegal cross-frame vector math, cross-datum orbit mixing, and unrotated antenna lever arm usage.

To satisfy Requirement R1:
1. Formalize zero-cost `#[repr(transparent)]` typestate wrappers in `gneiss-core::frames`:
   - Markers: `CoordinateFrame` (`Ecef<F>`, `Ned`, `Enu`, `BodyFrd`), `ReferenceFrame` (adding `Pz90`).
   - Primitives: `Point3<F>` / `EcefPos<F>`, `SpatialVector<F>`, `SpatialVelocity<F>`, `SpatialCovariance<F>`.
   - Invariants: `AntennaLeverArm(SpatialVector<BodyFrd>)`, `Attitude<From, To>`.
   - Operators: Restrict `Add`/`Sub` to prevent cross-frame arithmetic at compile time.
2. Implement `LocalTangentPlane<F>` to structurally enforce relational coupling from a single `EcefPos<F>` input.
3. Refactor `EskfState`, `DdSatGeometry`, `update_gnss_position`, `update_doppler_velocity`, and `eval_odaiba_ins` to require `BodyToEcef<F>` attitude rotation when projecting `AntennaLeverArm`.
4. Keep all file sizes < 500 LOC and function sizes <= 32 LOC, specifically guarding `eval_odaiba_ins/main.rs` (474 LOC).

---

## 5. Verification Method

To independently verify the observations, logic, and conclusions:
1. **Inspect Survey Findings**:
   Examine `/Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_1/survey_r1_spatial.md`.
2. **Verify Existing Workspace Tests**:
   Run `cargo test --workspace`. Verified: 0 failures, 100% passed across all crates:
   - `gneiss_core`: 20 tests ok
   - `gneiss_fetch`: 12 tests ok
   - `gneiss_geodesy`: 32 tests ok
   - `gneiss_ntrip`: 11 tests ok
   - `gneiss_parsers`: 41 tests ok
   - `gneiss_rtk`: 452 tests ok
   - `gneiss_tests` (integration): 16 tests ok
   - `gneiss_tests` (frontiers e2e): 205 tests ok
   - Total: 789 passed; 0 failed; 0 ignored.
3. **Verify Line Counts**:
   Run `wc -l crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs crates/gneiss-rtk/src/estimators/eskf/dd_update.rs crates/gneiss-core/src/frames/*.rs`.
4. **Compile-Fail Invariant Verification**:
   When implementing Phase 1, assert that the following fails to compile:
   ```rust
   let pos = EcefPos::<Itrf2014>::new(Vector3::zeros());
   let arm = AntennaLeverArm::new(0.5, 0.0, 0.0);
   let bad = pos + arm; // Must fail: mismatched frames (ECEF vs BodyFrd)
   ```

