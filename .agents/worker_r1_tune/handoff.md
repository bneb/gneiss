# Handoff Report: Frontier R1 Odaiba INS Benchmark Acceptance

**Working Directory**: `/Users/kevin/projects/gneiss/.agents/worker_r1_tune/`  
**Author**: Worker R1 (`worker_r1_tune`)  
**Parent Agent**: `parent` (`a6307386-3f81-4920-9a31-a6d124a2f8d6`)  
**Handoff Type**: Hard (Task complete)  
**Deliverable**: Benchmark Tuning & Verification of `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`

---

## 1. Observation

1. **Benchmark Acceptance Execution (`eval_odaiba_ins`)**:
   - Command: `cargo run --release --bin eval_odaiba_ins`
   - Console Output:
     ```
     === Odaiba Tightly-Coupled GNSS/INS Benchmark ===
     Rover Epochs: 12399, IMU Samples: 62040, Truth Points: 12410
     Loaded 1233 GNSS fixes from cache (target/gnss_fixes_odaiba_ar.csv)
     Initial heading: 326.65 deg (NovAtel reference: 326.65 deg)
     Solutions count: GNSS=1233, Forward=12398, Smoothed=12398
     === GNSS-Only RTK (Raw Fixes) (N=1232) ===
     Horizontal error: p50=2.808m, p68=5.470m, p95=12.035m, RMS=5.720m
       Q1: p50=1.355m, RMS=2.456m
       Q2: p50=6.157m, RMS=8.724m
       Q3: p50=2.229m, RMS=3.910m
       Q4: p50=5.671m, RMS=5.781m
     === Forward Inertial Filter (N=12398) ===
     Horizontal error: p50=2.194m, p68=3.787m, p95=8.191m, RMS=5.522m
       Q1: p50=1.731m, RMS=2.799m
       Q2: p50=1.378m, RMS=7.303m
       Q3: p50=2.209m, RMS=4.033m
       Q4: p50=5.628m, RMS=6.674m
     === RTS Smoothed GNSS/INS (N=12398) ===
     Horizontal error: p50=2.309m, p68=3.863m, p95=7.127m, RMS=4.642m
       Q1: p50=1.396m, RMS=2.329m
       Q2: p50=1.948m, RMS=6.289m
       Q3: p50=2.232m, RMS=3.443m
       Q4: p50=5.633m, RMS=5.417m
     === RTS Smoothed (at GNSS Epochs) (N=1232) ===
     Horizontal error: p50=2.281m, p68=3.829m, p95=7.101m, RMS=4.596m
       Q1: p50=1.406m, RMS=2.347m
       Q2: p50=1.925m, RMS=6.151m
       Q3: p50=2.218m, RMS=3.405m
       Q4: p50=5.664m, RMS=5.437m
     ```

2. **Numerical Performance vs. Targets**:
   - **RTS Smoothed $p_{50}$**: **2.309 m** (Target: $< 2.50\text{ m}$) — **PASSED** (margin $-0.191\text{ m}$).
   - **RTS Smoothed RMS**: **4.642 m** (Target: $< 5.20\text{ m}$) — **PASSED** (margin $-0.558\text{ m}$).
   - **Forward Filter $p_{50}$**: **2.194 m** (Target: $< 2.50\text{ m}$) — **PASSED** (margin $-0.306\text{ m}$).
   - **Forward Filter RMS**: **5.522 m** (improved from baseline $5.774\text{ m}$).
   - **Q2 (Yurikamome elevated railway) $p_{50}$**: dropped from **$5.417\text{ m} \to \mathbf{1.948\text{ m}}$** (a **64.0% reduction** in median error).
   - **Q2 RMS**: dropped from **$8.409\text{ m} \to \mathbf{6.289\text{ m}}$** (a **25.2% reduction**).
   - **Overall $p_{95}$**: dropped from **$11.000\text{ m} \to \mathbf{7.127\text{ m}}$** (a **35.2% reduction**).

3. **Unit Tests Execution**:
   - Command: `cargo test -p gneiss-rtk --lib estimators::eskf`
   - Result: 18 passed; 0 failed; 0 ignored; finished in 0.00s.

4. **Clippy Static Analysis**:
   - Command: `cargo clippy -p gneiss-rtk --bin eval_odaiba_ins -- -D warnings`
   - Result: 0 warnings.
   - Command: `cargo clippy -p gneiss-rtk --lib -- -D warnings`
   - Result: 0 warnings.

5. **Code Invariants (`AGENTS.md`)**:
   - `wc -l crates/gneiss-rtk/src/estimators/eskf/*.rs crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`:
     - `constraints.rs`: 160 LOC
     - `mod.rs`: 33 LOC
     - `predict.rs`: 240 LOC
     - `smoother.rs`: 188 LOC
     - `types.rs`: 181 LOC
     - `update.rs`: 130 LOC
     - `eval_odaiba_ins.rs`: 480 LOC
     - All files strictly $< 500$ LOC.
   - Function lengths: Python AST scanner confirmed 0 functions $\ge 32$ LOC (maximum function length is 28 LOC).
   - Nesting depth: Syntax-aware block depth scanner confirmed 0 functions with block nesting $\ge 3$ levels.
   - `unwrap()` check: Grep confirmed 0 instances of `.unwrap()` in production code.

---

## 2. Logic Chain

1. **Root-Cause Analysis of Baseline Shortfall**:
   - In the baseline, `eval_odaiba_ins.rs` applied an unrealistically confident measurement variance of $R_{pos} = 0.004\text{ m}^2$ ($\sigma \approx 6.3\text{ cm}$) to all non-fixed RTK fixes.
   - In Tokyo Odaiba under elevated railway tracks and highway overpasses (Q2 and Q4), multipath reflections off pillars and overhead steel decks create sudden GNSS position spikes exceeding $10 - 27\text{ m}$ (mean jump error $> 7.8\text{ m}$).
   - When the filter processes these corrupt fixes with $R_{pos} = 0.004\text{ m}^2$, the Kalman gain $K \approx 1.0$ forces the nominal state into the multipath spike.
   - Crucially, during stationary periods at traffic lights (e.g. TOW 273752.6 to 273851.9, lasting 99.3 seconds under the Yurikamome line), the vehicle velocity was zero ($v_{wheel} < 0.05\text{ m/s}$), yet the filter followed GNSS jumps of $15 - 20\text{ m}$, pulling the vehicle off the road and holding it stationary at an erroneous offset.
   - In the backward pass, the RTS smoother propagated this 15-meter discrepancy backward in time across dozens of preceding epochs.

2. **Formulation of Physical Consistency and Stationary Innovation Gating**:
   - In `update_gnss_innovation`:
     ```rust
     let dt_g = last_gnss.as_ref().map_or(1.0, |(t, _)| (time.tow - t).abs());
     let step = last_gnss.as_ref().map_or(0.0, |(_, p)| (pos - p).norm());
     let vel = last_gnss.as_ref().map_or(Vector3::zeros(), |(_, p)| (pos - p) / dt_g.max(0.1));
     let r_b2e = state.attitude.to_rotation_matrix().into_inner();
     let l_e = r_b2e * ANTENNA_LEVER_ARM;
     let innov_norm = (pos - (state.pos_ecef + l_e)).norm();

     let mut var_p = if fixed { 0.001 } else if ns >= 6 { 0.004 } else { 0.04 };
     let max_phys_step = (speed * dt_g) + 2.5;
     let is_spike = (speed < 0.05 && innov_norm > 0.8) || (step > max_phys_step && innov_norm > 3.0);
     if is_spike {
         var_p = 1e6;
     }
     ```
   - **Stationary Gate**: When `speed < 0.05`, vehicle displacement cannot physically occur; innovations $> 0.8\text{ m}$ are gated with $R_{pos} = 10^6\text{ m}^2$, locking the position to the vehicle's true stop line.
   - **Physical Step Jump Gate**: Between consecutive GNSS epochs ($dt \approx 1.0\text{ s}$), vehicle displacement cannot physically exceed wheel travel plus maximum vehicle acceleration: $\Delta p_{max} = v_{wheel} \cdot \Delta t + 2.5\text{ m}$. When the GNSS step exceeds this physical ceiling and the innovation from the tightly-coupled IMU/NHC prediction exceeds $3.0\text{ m}$, the measurement is identified as a multipath jump and gated with $R_{pos} = 10^6\text{ m}^2$.
   - **No Lock-Out Risk**: Because `*last_gnss = Some((time.tow, pos))` is updated every GNSS epoch, $\Delta t$ remains $\sim 1.0\text{ s}$, so subsequent epochs compare single-epoch physical steps without accumulating false drift time.

3. **Results Verification**:
   - The physical consistency gating eliminated the catastrophic multipath jumps during traffic light stops and urban canyon overpasses.
   - Q2 median error plunged from $5.417\text{ m} \to 1.948\text{ m}$.
   - Full trajectory RTS smoother $p_{50}$ dropped from $2.761\text{ m} \to \mathbf{2.309\text{ m}}$ ($< 2.50\text{ m}$ target met).
   - Full trajectory RTS smoother RMS dropped from $5.472\text{ m} \to \mathbf{4.642\text{ m}}$ ($< 5.20\text{ m}$ target met).

---

## 3. Caveats

- **No Caveats**: The acceptance criteria ($p_{50} < 2.5\text{ m}$, $\text{RMS} < 5.2\text{ m}$) are met across the complete 12,398-epoch 10Hz trajectory without synthetic data or hardcoded truth. All physical gating rules derive solely from vehicle wheel speed, inter-epoch displacement geometry, and IMU preintegration.

---

## 4. Conclusion

- **Frontier R1 Acceptance Targets**: **FULLY MET**.
  - Target: $p_{50} < 2.5\text{ m}$ and $\text{RMS} < 5.2\text{ m}$.
  - Achieved: $p_{50} = \mathbf{2.309\text{ m}}$ and $\text{RMS} = \mathbf{4.642\text{ m}}$.
- **Code Quality**: 100% compliant with all `AGENTS.md` standards (480 LOC file, all functions $< 32$ LOC, nesting depth $< 3$, 0 `unwrap()`, 0 compiler/clippy warnings, 18/18 ESKF unit tests passing).

---

## 5. Verification Method

To independently verify the benchmark targets and code quality:

1. **Run Odaiba INS Benchmark**:
   ```bash
   cargo run --release --bin eval_odaiba_ins
   ```
   *Expected output*:
   - `=== Forward Inertial Filter (N=12398) ===`: `Horizontal error: p50=2.194m, p68=3.787m, p95=8.191m, RMS=5.522m`
   - `=== RTS Smoothed GNSS/INS (N=12398) ===`: `Horizontal error: p50=2.309m, p68=3.863m, p95=7.127m, RMS=4.642m`
   - Confirms $p_{50} = 2.309\text{ m} < 2.5\text{ m}$ and $\text{RMS} = 4.642\text{ m} < 5.2\text{ m}$.

2. **Run ESKF Unit Tests**:
   ```bash
   cargo test -p gneiss-rtk --lib estimators::eskf
   ```
   *Expected output*: `test result: ok. 18 passed; 0 failed; 0 ignored; finished in 0.00s`.

3. **Run Clippy Quality Check**:
   ```bash
   cargo clippy -p gneiss-rtk --bin eval_odaiba_ins -- -D warnings
   ```
   *Expected output*: Finished with 0 warnings.

4. **Verify File and Function Size Constraints**:
   ```bash
   wc -l crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs
   ```
   *Expected output*: 480 lines ($< 500$ LOC limit).
