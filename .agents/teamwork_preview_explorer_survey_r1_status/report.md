# Frontier R1 Survey Report: 15-State ESKF GNSS/INS & Odaiba Benchmark

**Date**: 2026-09-13  
**Investigator**: Explorer R1 (`teamwork_preview_explorer_survey_r1_status`)  
**Target Milestone**: Frontier R1 (Milestone M1) — 15-State Error-State Kalman Filter (ESKF/MEKF) GNSS/INS  
**Target Benchmark**: `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs` (Tokyo Odaiba Urban Canyon, 12,398 epochs @ 10 Hz)  

---

## 1. Executive Summary

- **Unit Tests**: **18 / 18 passing** (`cargo test -p gneiss-rtk --lib estimators::eskf` in 0.01s).
- **Code Quality**: **100% compliant** with `AGENTS.md` (all files < 500 LOC, functions < 32 LOC, nesting < 3, 0 `unwrap()` in production, 0 Clippy warnings).
- **Benchmark Execution**: **100% complete** across the full 12,398-epoch 10Hz trajectory (62,040 IMU samples) running in **0.67s CPU time** (~1,850× faster than real-time).
- **Acceptance Criterion**: **NOT YET MET**.
  - Target: $p_{50} < 2.5\text{ m}$ and $\text{RMS} < 5.2\text{ m}$.
  - Current RTS Smoother: $p_{50} = \mathbf{2.761\text{ m}}$ (gap $+0.261\text{ m}$), $\text{RMS} = \mathbf{5.472\text{ m}}$ (gap $+0.272\text{ m}$).
  - Current Forward Filter: $p_{50} = \mathbf{3.131\text{ m}}$ (gap $+0.631\text{ m}$), $\text{RMS} = \mathbf{5.774\text{ m}}$ (gap $+0.574\text{ m}$).
- **Progress vs. 6-DOF Baseline**:
  - Forward Filter $p_{50}$ improved by **47.7%** (from $5.982\text{ m} \to 3.131\text{ m}$).
  - Forward Filter RMS improved by **40.4%** (from $9.689\text{ m} \to 5.774\text{ m}$).
  - RTS Smoothed $p_{50}$ improved from $2.907\text{ m} \to 2.761\text{ m}$.
  - RTS Smoothed RMS improved from $5.508\text{ m} \to 5.472\text{ m}$.

---

## 2. Exact Numerical Benchmark Results

The benchmark evaluates the complete 20-minute Tokyo Odaiba drive ($N=12,398$ 10Hz epochs, 62,040 IMU samples) against NovAtel SPAN ground-truth reference coordinates:

| Trajectory Stage | Epochs ($N$) | $p_{50}$ (m) | $p_{68}$ (m) | $p_{95}$ (m) | Max (m) | RMS (m) | Target Met? |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **GNSS-Only RTK (Raw Fixes)** | 1,232 | 2.808 | 5.470 | 12.035 | 27.693 | 5.720 | N/A |
| **Forward Inertial Filter (15-state ESKF)** | 12,398 | 3.131 | 5.399 | 12.936 | 29.636 | 5.774 | ❌ No ($p_{50} > 2.5$, $\text{RMS} > 5.2$) |
| **RTS Smoothed GNSS/INS (Full 10Hz)** | 12,398 | **2.761** | 5.354 | **11.000** | **21.196** | **5.472** | ❌ No ($p_{50} > 2.5$, $\text{RMS} > 5.2$) |
| **RTS Smoothed (at GNSS Epochs)** | 1,232 | **2.724** | 5.374 | 10.795 | 21.196 | 5.484 | ❌ No ($p_{50} > 2.5$, $\text{RMS} > 5.2$) |

### Quartile Breakdown (Chronological Progression)

The error distribution exhibits severe spatial non-uniformity correlated with elevated rail/expressway structures:

| Trajectory Quartile | Forward Filter $p_{50}$ | Forward Filter RMS | RTS Smoothed $p_{50}$ | RTS Smoothed RMS | Urban Environment |
|:---|:---:|:---:|:---:|:---:|:---|
| **Q1 (Epochs 1–3,099)** | **1.232 m** | **2.581 m** | **1.311 m** | **2.255 m** | Open sky / wide boulevards (**Passes target**) |
| **Q2 (Epochs 3,100–6,199)** | 6.150 m | 8.823 m | 5.417 m | 8.409 m | Under Yurikamome elevated railway (Severe NLOS multipath) |
| **Q3 (Epochs 6,200–9,299)** | **2.221 m** | **3.790 m** | **2.240 m** | **3.502 m** | Open waterfront / open sky (**Passes target**) |
| **Q4 (Epochs 9,300–12,398)** | 5.690 m | 5.874 m | 5.644 m | 5.633 m | Shuto Expressway overpass (Severe NLOS multipath) |

*Key Takeaway*: In open-sky sections (Q1 & Q3), the 15-state ESKF with RTS smoothing achieves $p_{50} \approx 1.3 - 2.2\text{ m}$ and $\text{RMS} \approx 2.2 - 3.5\text{ m}$, easily meeting the targets. The shortfall is driven entirely by the severe multipath sections under elevated structures in Q2 and Q4.

### Execution Performance

- **Execution Mode**: `cargo run --release --bin eval_odaiba_ins`
- **Elapsed Wall-Clock Time**: ~1.28 s
- **CPU User Time**: 0.56 s (0.11 s system, ~0.67 s total CPU)
- **Throughput**: ~18,500 epochs/second (~1,850× real-time continuous 10Hz/50Hz processing)

---

## 3. Code Quality & Standards Verification (`AGENTS.md`)

All production and benchmark files were inspected for compliance:

| File | LOC | Max Function LOC | Max Nesting Depth | `unwrap()` Count | Clippy Warnings |
|:---|:---:|:---:|:---:|:---:|:---:|
| `crates/gneiss-rtk/src/estimators/eskf/mod.rs` | 33 | 0 (pure exports) | 0 | 0 | 0 |
| `crates/gneiss-rtk/src/estimators/eskf/types.rs` | 181 | 16 | 2 | 0 | 0 |
| `crates/gneiss-rtk/src/estimators/eskf/predict.rs` | 240 | 18 | 2 | 0 | 0 |
| `crates/gneiss-rtk/src/estimators/eskf/update.rs` | 130 | 18 | 2 | 0 | 0 |
| `crates/gneiss-rtk/src/estimators/eskf/constraints.rs` | 160 | 23 | 2 | 0 | 0 |
| `crates/gneiss-rtk/src/estimators/eskf/smoother.rs` | 188 | 24 | 2 | 0 | 0 |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs` | 496 | 26 | 2 | 0 | 0 |

- **File size rule (< 500 LOC)**: All files comply. `eval_odaiba_ins.rs` is at 496 LOC.
- **Function size rule (< 32 LOC)**: All functions comply (max function length is 26 lines).
- **Nesting depth rule (< 3 levels)**: All functions comply (max nesting depth is 2 levels).
- **No `unwrap()` in production code**: Verified 0 instances across all production files (`unwrap_used = "deny"` in workspace lints).
- **Positive velocity-attitude coupling**: Verified strictly positive `vel_att = +f_e_skew * dt` in `predict.rs:38` as mandated by `AGENTS.md` and guarded by `test_transition_matrix_velocity_attitude_coupling_is_strictly_positive`.

---

## 4. Root-Cause Diagnosis of Performance Shortfall

### Diagnosis 1: Overconfident GNSS Measurement Covariance ($R_{\text{pos}}$)
In `eval_odaiba_ins.rs:260-262`:
```rust
let var_p = if fixed { 0.001 } else if ns >= 6 { 0.004 } else { 0.04 };
let r_pos = nalgebra::Matrix3::from_diagonal(&Vector3::new(var_p, var_p, var_p));
```
- A variance of $0.001\text{ m}^2$ corresponds to $\sigma \approx 3.1\text{ cm}$; $0.004\text{ m}^2$ is $\sigma \approx 6.3\text{ cm}$; $0.04\text{ m}^2$ is $\sigma \approx 20\text{ cm}$.
- In Q2 and Q4, raw GNSS errors exceed $12\text{ m}$ (max $27.69\text{ m}$) due to multipath reflections off elevated railway pillars.
- Because the filter assigns centimetre-level uncertainty to these corrupted GNSS fixes, the Kalman gain heavily favors the erroneous measurement over the IMU/NHC prediction, violently pulling the filter state into the multipath spike.

### Diagnosis 2: Crude Finite-Difference Velocity Updates
In `eval_odaiba_ins.rs:259`:
```rust
let vel = last_gnss.as_ref().map_or(Vector3::zeros(), |(_, p)| (pos - p) / dt_g.max(0.1));
```
- Synthesizing GNSS velocity by finite-differencing noisy positions amplifies position errors: a 10m multipath jump over 0.1s produces a 100 m/s velocity impulse.
- In `update_gnss_pos_vel`, this velocity innovation directly perturbs velocity and attitude via coupling Jacobians, destabilizing heading estimation.

### Diagnosis 3: Missing Chi-Square Innovation Gating / Outlier Rejection
In `update.rs:83-94`:
- `update_gnss_pos_vel` blindly updates the state without verifying the Mahalanobis distance $\gamma = \mathbf{y}^T \mathbf{S}^{-1} \mathbf{y}$.
- During severe multipath, fixes with innovations $> 15\text{ m}$ should be gated or downweighted. Without gating, corrupt fixes contaminate the forward filter history and cannot be fully purged by RTS smoothing.

### Diagnosis 4: Unmodeled Antenna Lever Arm
In `eval_odaiba_ins.rs:248`:
```rust
const ANTENNA_LEVER_ARM: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
```
- As documented in `TIER1_ROADMAP.md:140` and `PROJECT_STATUS.md:1251`, Tokyo Odaiba vehicle roof antenna is physically separated from the IMU by $\sim 0.70\text{ m}$.
- Every turn introduces an apparent lever-arm discrepancy of up to $1.4\text{ m}$ if unmodeled. While the 15-state ESKF has exact lever-arm Jacobian compensation implemented (`update.rs:40-53`), it is unused due to the zero lever-arm constant.

---

## 5. Actionable Tuning Recommendations

To close the remaining $0.26\text{ m}$ gap on $p_{50}$ and $0.27\text{ m}$ on RMS, the following changes are recommended for Worker R1 / M1:

1. **Robust Innovation Gating / Adaptive Measurement Covariance**:
   - Implement Chi-square ($\chi^2$) gating or Huber weighting in `update_gnss_innovation`:
     ```rust
     let innov_norm = (pos - state.pos_ecef).norm();
     let var_p = if fixed && innov_norm < 1.0 {
         0.04 // ~20cm realistic RTK fix error in canyon
     } else if ns >= 6 && innov_norm < 3.0 {
         0.25 // ~50cm float error
     } else if innov_norm < 6.0 {
         2.0  // ~1.4m degraded
     } else {
         1e6  // Reject outlier > 6m
     };
     ```
   - This allows the IMU preintegration and wheel-speed NHC to coast safely through elevated railway obstructions without being pulled by multipath spikes.

2. **Omit or Doppler-Gate GNSS Velocity Updates**:
   - Defer velocity estimation entirely to CAN wheel-speed updates (`update_body_velocity`), or use Doppler measurements if available.
   - Set `var_v = 1e6` (effectively disabled) when GNSS velocity is derived from finite differences.

3. **Tune IMU Process Noise ($Q$)**:
   - In `default_q_diag()`:
     - Reduce velocity random walk $q_a$ from $0.05$ to $0.005 - 0.01\text{ (m/s}^2)^2/\text{Hz}$ to keep the filter stiffer during GNSS multipath spikes.
     - Increase gyro bias random walk $q_{bg}$ from $10^{-9}$ to $10^{-7}\text{ (rad/s)}^2/\text{Hz}$ to accelerate gyro bias convergence.

4. **Calibrate Non-Zero Antenna Lever Arm**:
   - Configure the physical antenna lever arm (e.g. `Vector3::new(0.15, 0.0, -0.70)`) to leverage the already implemented lever-arm attitude Jacobian.
