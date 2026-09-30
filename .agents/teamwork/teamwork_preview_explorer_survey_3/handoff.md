# Handoff Report: Survey 3 — Estimators, Pipelines, Relational Coupling & Benchmarks (R3 & R4)

**Agent:** Survey Explorer 3  
**Working Directory:** `/Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_3`  
**Handoff Type:** Hard  

---

## 1. Observation

### 1.1 Estimator State & Pipeline Bare Vector Usage
- `crates/gneiss-rtk/src/estimators/eskf/types.rs:33-40`:
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
- `crates/gneiss-rtk/src/estimators/eskf/update.rs:42-50`:
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
  ```
- `crates/gneiss-rtk/src/estimators/eskf/dd_update.rs:10-24`:
  ```rust
  pub struct DdSatGeometry {
      pub u_s: Vector3<f64>,
      pub u_ref: Vector3<f64>,
      pub lever_arm: Vector3<f64>,
  }
  ```
- `crates/gneiss-rtk/src/swfg/engine/mod.rs:41-49`:
  `prev_position: Option<Vector3<f64>>`, `initial_position: Option<Vector3<f64>>`.
- `crates/gneiss-rtk/src/post_process/mod.rs:64-67`:
  `pub base_position: Option<Vector3<f64>>`, `pub initial_rover_position: Option<Vector3<f64>>`.
- `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs:24`:
  `const ANTENNA_LEVER_ARM: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);`

### 1.2 Latent Coordinate Frame Bug in `receiver_pcv`
- `crates/gneiss-parsers/src/receiver_pcv/mod.rs:173-178`:
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
      let (az_bas_sat, el_bas_sat) = az_el(bas_llh, bas_llh, sat_pos);
      let (az_bas_ref, el_bas_ref) = az_el(bas_llh, bas_llh, ref_sat_pos);
  ```
  `az_el` in `crates/gneiss-core/src/coords.rs:177`:
  `pub fn az_el(pos_llh: Vector3<f64>, pos_ecef: Vector3<f64>, sat_ecef: Vector3<f64>) -> (f64, f64)`
  Here, `rov_llh` (latitude/longitude in radians, height in metres) was passed as BOTH `pos_llh` AND `pos_ecef`. Inside `az_el`, line-of-sight delta was calculated as `sat_ecef - rov_llh` instead of `sat_ecef - rov_ecef`.

### 1.3 Disconnected Relational Coupling in `receiver_antenna`
- `crates/gneiss-parsers/src/receiver_antenna/mod.rs:277-305`:
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
      let zen_ref = 90.0 - el_ref_rad.to_degrees();
      ...
      let (Some(bas_s), Some(bas_r)) = (
          base_ant.pcv_mm_az_zen(freq_code, az_sat_bas, zen_sat),
          base_ant.pcv_mm_az_zen(freq_code, az_ref_bas, zen_ref),
      ) ...
  ```
  The base antenna correction directly uses `zen_sat` and `zen_ref` computed from the rover's elevation angle.

### 1.4 Golden Test & Benchmark Results
- `cargo test --workspace`: 789 passed, 0 failed, 1 ignored.
- `python3 scripts/check_network_benchmark.py --smoke`:
  Output: `ALL CHECKS PASSED`. Network fused horizontal $p_{50} = 0.022\text{ m} \le 0.04$, horizontal RMS $0.032\text{ m} \le 0.06$, vertical RMS $0.037\text{ m} \le 0.08$.
- `python3 scripts/check_multignss_benchmark.py --smoke`:
  Output: `ALL CHECKS PASSED`. Network fused fix rate $99.10\% \ge 96.5\%$, P181 fix rate $97.80\% \ge 97.5\%$.
- `cargo run --release --bin eval_odaiba_ins`:
  Output: RTS Smoothed GNSS/INS horizontal error: $p_{50} = 1.751\text{ m} \le 1.80\text{ m}$, $\text{RMS} = 3.479\text{ m} \le 3.50\text{ m}$, 0 false fixes.
- Target file line counts: `main.rs` (475 LOC), `formation.rs` (457 LOC), `dd_update.rs` (454 LOC), `swfg/engine/mod.rs` (430 LOC).

---

## 2. Logic Chain

1. **Bare Vector Insecurity (from 1.1):** `EskfState`, `build_gnss_pos_system`, and `eval_odaiba_ins` accept untyped `Vector3<f64>` across coordinate frames. Adding a Body FRD vector (`lever_arm`) directly to an ECEF vector (`pos_ecef`) compiles without error.
2. **Latent Frame Bug Consequence (from 1.2):** Because `Vector3<f64>` has no coordinate tag, `rov_llh` was passed as `pos_ecef` in `receiver_pcv::dd_correction_from_geometry`, corrupting ECEF delta calculations. Structurally typed `EcefPos<F>` and `LlhPos<F>` would have caused this call to fail at compile-time.
3. **Relational Geometry Coupling (from 1.3):** Passing 4 or 6 bare floats to `dd_correction_m` or `compute_dd_pcv_correction_2d` allows callers to provide physically contradictory angles or borrow rover zenith angles for the base station. Constructing a single `DoubleDiffGeometry<F>` struct directly from `(rov_pos, bas_pos, sat_pos, ref_sat_pos)` guarantees that rover and base view the exact same satellite positions at the same epoch in the same reference frame, making disagreement unrepresentable.
4. **Benchmark Preservation (from 1.4):** Current metrics ($p_{50} = 1.751\text{ m}$, network $p_{50} = 0.022\text{ m}$) reflect a high-performing baseline. The refactoring is purely an architectural hardening and frame safety enforcement; any regression in these benchmark metrics indicates a coordinate or math transformation defect.
5. **Code Standard Constraints (from 1.4):** Several files are within 25-45 LOC of the 500 LOC ceiling. Refactoring cannot inline large struct definitions or boilerplate into `main.rs`, `dd_update.rs`, or `formation.rs`. Shared types must reside in a separate module (`geometry.rs`).

---

## 3. Caveats

- **Network RTK Benchmark Execution Time:** Full evaluation takes ~30-40 minutes; `--smoke` evaluates 1,800 epochs in ~40 seconds. For rapid development iterations, `--smoke` must be used before running full CI sweeps.
- **Odaiba Fix Cache:** `eval_odaiba_ins` uses `target/gnss_fixes_odaiba_ar.csv` if present (taking ~7s). If `cargo clean` is run, `run_gnss_rtk_pass` executes the full RTK GNSS pass (~2 minutes) to regenerate the cache.
- **Frame Trait Integration:** Survey 1 is defining the base `EcefPos<F>` and `ReferenceFrame` traits. Survey 3 assumes `EcefPos<F>`, `AntennaLeverArm`, and `NedCovariance` will be made available in `gneiss-core`.

---

## 4. Conclusion

1. **R3 Refactoring Scope:** Refactor `EskfState`, `build_gnss_pos_system`, `build_doppler_velocity_system`, `compute_dd_jacobian_15`, `SwfgEngine`, `PostProcessOptions`, and `eval_odaiba_ins` to require:
   - `EcefPos<F>` for positions.
   - `EcefVector<F>` for ECEF velocities.
   - `BodyFrdVector` for accelerometer/gyroscope biases.
   - `AntennaLeverArm` for Body FRD offsets, requiring explicit attitude rotation $C_b^e$ before spatial addition.
   - `NedCovariance` with explicit projection to ECEF covariance.
2. **R4 Relational Coupling Scope:** Introduce `DoubleDiffGeometry<F>` constructed solely from `(rov_pos, bas_pos, sat_pos, ref_sat_pos)`:
   - Derive both rover and base azimuth, elevation, ranges, and unit LOS vectors internally.
   - Update `receiver_pcv` and `receiver_antenna` to consume `DoubleDiffGeometry<F>`.
   - Fix the argument inversion bug in `receiver_pcv::dd_correction_from_geometry`.
3. **Benchmark Budgets:** All 789 tests must remain 100% green; dual CI smoke guards must pass; `eval_odaiba_ins` must maintain $p_{50} \le 1.80\text{ m}$ and $\text{RMS} \le 3.50\text{ m}$ with 0 false fixes.

---

## 5. Verification Method

To independently verify all findings and baselines:

```bash
# 1. Run full workspace test suite (assert 789 passed, 0 failed)
cargo test --workspace

# 2. Build release binary for benchmarks
cargo build --release --bin eval_network_ppk --bin eval_odaiba_ins

# 3. Verify CI Network RTK Smoke Guard (assert "ALL CHECKS PASSED")
python3 scripts/check_network_benchmark.py --smoke

# 4. Verify CI Multi-GNSS Smoke Guard (assert "ALL CHECKS PASSED")
python3 scripts/check_multignss_benchmark.py --smoke

# 5. Verify Tokyo Odaiba INS Benchmark (assert p50 <= 1.80m, RMS <= 3.50m)
cargo run --release --bin eval_odaiba_ins

# 6. Verify Code Standard Linting (assert 0 warnings)
cargo clippy --workspace --all-targets -- -D warnings
```

**Invalidation Conditions:**
- Any compiler warning or error under `cargo clippy --all-targets -- -D warnings`.
- `eval_odaiba_ins` RTS smoothed horizontal $p_{50} > 1.80\text{ m}$ or $\text{RMS} > 3.50\text{ m}$.
- Any failure in `check_network_benchmark.py --smoke` or `check_multignss_benchmark.py --smoke`.
- Any file exceeding 500 LOC or function exceeding 32 LOC.
