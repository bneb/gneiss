# Survey Report: R2 Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting

**Author**: Survey Explorer 2  
**Date**: 2026-09-24  
**Working Directory**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_2`  
**Parent Task**: Gneiss Urban Canyon Fix Rate Expansion and Multipath Mitigation (2026-09-24T13:30:49Z)  
**Target Requirement**: R2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting  

---

## Executive Summary

In dense urban canyons such as Tokyo Shinjuku and Hong Kong Whampoa, Non-Line-of-Sight (NLOS) reflections and severe pseudorange multipath produce sudden 5m–20m steps in GNSS code measurements. Carrier phase observations, by contrast, have millimetre-level precision and their multipath error is mathematically bounded by $\lambda / 4 \approx 4.8\text{ cm}$ (for GPS L1, $\lambda \approx 19.0\text{ cm}$). When code multipath goes undetected:
1. It biases the float Kalman filter position and ambiguity states ($\delta N \sim 25\text{--}100\text{ cycles}$).
2. It trips the Melbourne–Wübbena slip detector (`mw.rs`), resetting wide-lane integer tracks on clean carrier arcs.
3. In existing gross-error screening (`screen.rs`), detecting a $> 15\text{m}$ pseudorange residual causes the entire `DoubleDiffMeasurement` to be dropped, inadvertently throwing away clean carrier phase measurements.
4. Float ambiguities diverge from true integers, causing LAMBDA / FFRT ratio tests to fail, collapsing urban canyon fix rates to 6%–18% (compared to > 60% in commercial Tier-1 engines).

This report presents a comprehensive investigation of the Gneiss RTK codebase baseline (`crates/gneiss-rtk/src/estimators/rtk_iekf/`), details the exact mechanics of double-difference formation, tracking arc identification, robust IEKF updates, and formulates an optimal, production-ready algorithm for Code-Minus-Carrier (CMC) multipath detection, adaptive variance inflation, and carrier phase preservation.

---

## 1. Codebase Baseline & Structural Inventory

### 1.1 Double-Difference Observation Formation (`formation.rs`)
Double-difference (DD) measurements are constructed each epoch in `GnssRtkIekf::build_dd_measurements`:
- **Constellation Selection**: Satellites are grouped by constellation (`select_constellations`). For GPS/QZSS, Galileo, BeiDou, and GLONASS, reference satellites are selected per constellation with elevation-based hysteresis (`select_ref_sat_with_hysteresis`).
- **Pair Formation**: For each non-reference satellite $s$ and reference satellite $k$, `form_pair_dd` constructs DD measurements across canonical frequency bands (`canonical_bands_for_constellation`):
  - Rover observations: $P_r^s, \Phi_r^s$ (satellite), $P_r^k, \Phi_r^k$ (reference).
  - Base observations: $P_b^s, \Phi_b^s$ (satellite), $P_b^k, \Phi_b^k$ (reference).
  - Double-difference pseudorange:
    $$\nabla\Delta P = (P_r^s - P_r^k) - (P_b^s - P_b^k)$$
  - Double-difference carrier phase (cycles):
    $$\nabla\Delta \Phi = (\Phi_r^s - \Phi_r^k) - (\Phi_b^s - \Phi_b^k)$$
- **DoubleDiffKey**: Uniquely identifies the DD pair:
  ```rust
  pub struct DoubleDiffKey {
      pub constellation_id: u8,
      pub sat: u16,
      pub ref_sat: u16,
      pub freq_band: u8,
  }
  ```
- **Measurement Struct**: Packed into `DoubleDiffMeasurement`:
  ```rust
  pub struct DoubleDiffMeasurement {
      pub key: DoubleDiffKey,
      pub dd_pr_m: f64,
      pub dd_cp_cycles: Option<f64>,
      pub sat_pos: Vector3<f64>,
      pub ref_pos: Vector3<f64>,
      pub base_pos: Vector3<f64>,
      pub lambda: f64,
      pub pr_var_m2: f64,
      pub cp_var_cycles2: f64,
      pub pr_ref_var_m2: f64,
      pub cp_ref_var_cycles2: f64,
      pub dm_wet_rov: f64,
      pub dgrad_n_rov: f64,
      pub dgrad_e_rov: f64,
      pub tide_dd_m: f64,
      pub dd_pcv_m: f64,
  }
  ```

### 1.2 Tracking Arc Identification (`formation.rs`, `screening.rs`)
Tracking arcs are tracked across epochs using a composite mechanism:
- **`CycleSlipDetector` (`crates/gneiss-rtk/src/post_process/screening.rs`)**:
  Maintains `slip_counts: HashMap<SatelliteId, u32>`. The arc ID is returned by `get_arc(sat) -> u32`.
  A cycle slip is flagged and `slip_counts` incremented if any of the following occur:
  1. **LLI bit 0**: `(lli & 1) != 0` (loss of lock indicated by receiver firmware).
  2. **Geometry-Free (GF) phase jump**: $|\Delta(\lambda_1 \Phi_1 - \lambda_2 \Phi_2)| > 0.05\text{ m}$ (5 cm jump).
  3. **Time gap**: $\Delta t > 2.0\text{ s}$ or carrier phase jump $> 10^7\text{ cycles}$.
  4. **Doppler discrepancy**: $|\Delta\Phi - (-f_D \Delta t) - \text{median}| > 1.0 \times \Delta t\text{ cycles}$.
- **Pair Slip Checking (`formation.rs:404-423`)**:
  ```rust
  fn check_pair_slip(...) -> bool {
      let mut cur_arc = slip_detector.get_arc(sat_id) + slip_detector.get_arc(ref_sat);
      if widelane_ar {
          cur_arc += base_slip_detector.get_arc(sat_id) + base_slip_detector.get_arc(ref_sat);
      }
      let arc_changed = prev_arcs.insert(key, cur_arc).is_some_and(|prev| prev != cur_arc);
      rov_s.get_lli(key.freq_band).is_some_and(|l| (l & 1) != 0)
          || rov_ref.get_lli(key.freq_band).is_some_and(|l| (l & 1) != 0)
          || arc_changed
  }
  ```
  **Key Insight**: A continuous carrier tracking arc for a satellite pair is strictly bounded between consecutive cycle slips. During a continuous arc, the integer carrier phase ambiguity $N$ is invariant:
  $$\frac{d}{dt} (\lambda N) = 0$$

### 1.3 Measurement Covariance Formulation (`formation_cov.rs`)
In `formation_cov.rs`:
- Base pseudorange and carrier variances are evaluated via `compute_dd_variances` from elevation angle and SNR:
  ```rust
  let pr_sat_var = single_diff_var(0.20, snrs.0, snrs.2, el_s);
  let pr_ref_var = single_diff_var(0.20, snrs.1, snrs.3, el_r);
  let cp_sat_var = single_diff_var(0.003 / lambda, snrs.0, snrs.2, el_s);
  let cp_ref_var = single_diff_var(0.003 / lambda, snrs.1, snrs.3, el_r);
  ```
- Nominal single-difference code $\sigma = 0.20\text{ m}$ ($0.04\text{ m}^2$).
- Nominal single-difference phase $\sigma = 0.003\text{ m} / \lambda \approx 0.015\text{ cycles}$ ($0.00025\text{ cycles}^2$).
- `single_diff_var` scales with elevation ($\sin^{-2}(\theta)$) and C/N0 attenuation (`snr_weight`, `attenuation_scale`).
- **Baseline Observation**: `formation_cov.rs` models static elevation and SNR dependencies, but has no mechanism to observe dynamic code multipath deviations or code-carrier divergence along an active tracking arc.

### 1.4 Linear System Construction & Robust Weighting (`update/system.rs`, `robust.rs`)
- **Row Separation**: In `system.rs::append_dd_meas_rows`, code and carrier phase measurements are processed into **independent rows**:
  - Code row: `append_dd_code_row` produces row in $H$, residual $y_{\text{code}}$, and diagonal variance $R_{\text{code}}$.
  - Phase row: `append_dd_phase_row` produces row in $H$, residual $y_{\text{phase}}$, and diagonal variance $R_{\text{phase}}$.
- **Existing Pseudorange Inflation (`system.rs:115-123`)**:
  ```rust
  fn effective_code_variance(m: &DoubleDiffMeasurement, pr_y: f64) -> f64 {
      let base_r = m.pr_var_m2.max(0.01);
      if m.dd_cp_cycles.is_some() && pr_y.abs() > 3.0 {
          let excess = pr_y.abs() - 3.0;
          base_r * (1.0 + excess * excess)
      } else {
          base_r
      }
  }
  ```
- **Existing Code Blunder Exclusion (`system.rs:106-113`)**:
  ```rust
  fn is_code_blunder(m: &DoubleDiffMeasurement, pr_y: f64, pr_r: f64) -> bool {
      let nis = pr_y * pr_y / pr_r;
      if m.dd_cp_cycles.is_none() {
          pr_y.abs() > 15.0 && nis > 36.0
      } else {
          pr_y.abs() > 10.0 && nis > 25.0
      }
  }
  ```
- **Huber Variance Inflation (`robust.rs:32-40`)**:
  ```rust
  pub fn robust_inflate(innovation: f64, variance: f64, gate_scale: f64) -> f64 {
      let threshold = ROBUST_INNOVATION_THRESHOLD * gate_scale; // threshold = 9.0
      let nis = innovation * innovation / variance;
      if nis > threshold {
          variance * (nis / threshold)
      } else {
          variance
      }
  }
  ```
- **Off-Diagonal Correlation Assembly (`system.rs:237-251`)**:
  ```rust
  fn fill_dd_covariances(r: &mut DMatrix<f64>, metas: &[RowMeta]) {
      for (i, mi) in metas.iter().enumerate() {
          for (j, mj) in metas.iter().enumerate().skip(i + 1) {
              let same_group = mi.kind == mj.kind
                  && mi.constellation_id == mj.constellation_id
                  && mi.ref_sat == mj.ref_sat
                  && mi.freq_band == mj.freq_band;
              if same_group {
                  let cov = mi.ref_var.min(mj.ref_var);
                  r[(i, j)] = cov;
                  r[(j, i)] = cov;
              }
          }
      }
  }
  ```
  **Critical Architectural Invariant**: Off-diagonal elements exist ONLY between rows of the same `ObsKind` (`mi.kind == mj.kind`). **Code rows NEVER cross-correlate with phase rows in $R$.**

### 1.5 Identified Defects & Vulnerabilities in Existing Pipeline

#### Critical Defect 1: Discarding Clean Phase during Code Blunders (`screen.rs:43-56`)
In `screen.rs`:
```rust
pub fn screen_gross_pr_errors(
    measurements: &mut Vec<DoubleDiffMeasurement>,
    pos_pred: Vector3<f64>,
) -> Vec<DoubleDiffKey> {
    let mut rejected = Vec::new();
    for _ in 0..MAX_GROSS_PR_REJECTIONS_PER_EPOCH {
        let Some((idx, residual)) = worst_pr_residual(measurements, pos_pred) else { break };
        if residual.abs() <= GROSS_PR_ERROR_THRESHOLD_M { // 15.0m
            break;
        }
        rejected.push(measurements.remove(idx).key); // DROPS THE ENTIRE DoubleDiffMeasurement!
    }
    rejected
}
```
**Impact**: When a 15m–20m code multipath reflection hits a satellite in an urban canyon, `screen_gross_pr_errors` calls `measurements.remove(idx)`, deleting the entire `DoubleDiffMeasurement`. This completely removes `dd_cp_cycles` from the filter update! The filter is starved of millimetre-precision carrier phase information simply because the code was reflected.

#### Critical Defect 2: Melbourne–Wübbena Slip False Alarm on Code Multipath (`mw.rs:74-80`)
In `mw.rs`, the MW observable is tracked via `WidelaneTracker`:
$$MW = (\nabla\Delta \Phi_1 - \nabla\Delta \Phi_2) - \frac{\nabla\Delta R_N}{\lambda_W}, \quad \nabla\Delta R_N = \frac{f_1 \nabla\Delta P_1 + f_2 \nabla\Delta P_2}{f_1 + f_2}, \quad \lambda_W = \frac{c}{f_1 - f_2}$$
In `MwTrack::absorb`:
```rust
fn absorb(&mut self, x: f64) {
    if self.n >= INNOVATION_ARM_EPOCHS && (x - self.mean).abs() > SLIP_INNOVATION_CYCLES {
        *self = MwTrack::new(x); // FALSE CYCLE SLIP RESET!
    } else {
        self.push(x);
    }
}
```
For GPS L1/L2, $\lambda_W \approx 0.862\text{ m}$. A 10m pseudorange multipath jump on $P_1$ causes a change in $\nabla\Delta R_N$:
$$\Delta R_N = \frac{f_1}{f_1 + f_2} \Delta P_1 \approx 0.562 \times 10\text{ m} = 5.62\text{ m}$$
In MW cycles:
$$\Delta MW = -\frac{5.62\text{ m}}{0.862\text{ m}} = -6.52\text{ cycles}$$
Since `SLIP_INNOVATION_CYCLES = 1.0`, $|-6.52| > 1.0$ triggers an immediate arc reset! The wide-lane tracker is reset, destroying wide-lane ambiguity fixing even though carrier phase tracked continuously without a slip.

#### Critical Defect 3: Ambiguity Seeding Corruption (`mod.rs:430-445`)
When a satellite pair is newly initialized or reset after an actual slip:
```rust
let init_amb = dd_cp.map_or(0.0, |cp| cp - dd_pr / lambda);
```
If `dd_pr` carries a 15m multipath blunder:
$$\Delta N_{\text{seed}} = -\frac{15.0\text{ m}}{0.1903\text{ m}} \approx -78.8\text{ cycles}$$
The ambiguity is initialized with an 79-cycle error, poisoning the float filter state.

---

## 2. Mathematical Derivations & Theoretical Foundation

### 2.1 The Physics of GNSS Multipath: Carrier Phase vs Pseudorange
Consider a line-of-sight signal and a specular reflection from a vertical building wall:
- Direct signal: $s_d(t) = A \cos(\omega t - \phi_d)$
- Reflected signal: $s_r(t) = \alpha A \cos(\omega t - \phi_d - \Delta\phi)$, where $\Delta\phi = \frac{2\pi}{\lambda} \Delta s$ and $\Delta s$ is the excess path delay.
- The composite RF carrier received at the antenna:
  $$s_c(t) = A \sqrt{1 + \alpha^2 + 2\alpha \cos\Delta\phi} \cos(\omega t - \phi_d - \theta_{\text{mp}})$$
  where the carrier phase tracking error is:
  $$\theta_{\text{mp}} = \arctan\left( \frac{\alpha \sin\Delta\phi}{1 + \alpha \cos\Delta\phi} \right)$$
  For any reflection attenuation $\alpha < 1$:
  $$|\theta_{\text{mp}}| \le \arcsin(\alpha) < \frac{\pi}{2}\text{ rad}$$
  Converting to distance:
  $$|M_\Phi| = \frac{\lambda}{2\pi} |\theta_{\text{mp}}| \le \frac{\lambda}{4} \approx 4.76\text{ cm for GPS L1}$$
  In real urban environments, carrier multipath is typically $5\text{--}20\text{ mm}$.
- By contrast, for pseudorange (code tracking via DLL correlator):
  The chipping rate (e.g. 1.023 MHz for GPS C/A, chip length $T_c \approx 293\text{ m}$) means the cross-correlation peak is distorted over the entire range $\Delta s \in [0, 1.5 T_c]$. For NLOS reflections where the direct line-of-sight is blocked, the tracking loop locks entirely onto the reflected path:
  $$M_P = \Delta s \in [5\text{m}, 50\text{m}]$$
  **Fundamental Asymmetry**: Under urban multipath, pseudorange errors jump by 5m–20m, while carrier phase errors are confined to centimeters ($\le 0.05\text{ m}$).

### 2.2 Dual-Frequency Multipath Combination ($MP$)
For dual-frequency observations on frequencies $f_1, f_2$:
The observation equations in distance units are:
$$P_1 = \rho + c(dt_r - dt^s) + T + I_1 + M_{P,1} + \epsilon_{P,1}$$
$$\Phi_1 = \rho + c(dt_r - dt^s) + T - I_1 + \lambda_1 N_1 + M_{\Phi,1} + \epsilon_{\Phi,1}$$
$$\Phi_2 = \rho + c(dt_r - dt^s) + T - I_2 + \lambda_2 N_2 + M_{\Phi,2} + \epsilon_{\Phi,2}$$
Ionospheric delay scales as $I_2 = \alpha I_1$, where $\alpha = (f_1 / f_2)^2$.
Subtracting the two carrier phases eliminates geometry, clocks, and troposphere:
$$\Phi_1 - \Phi_2 = (\alpha - 1) I_1 + (\lambda_1 N_1 - \lambda_2 N_2) + (M_{\Phi,1} - M_{\Phi,2})$$
Solving for $I_1$:
$$I_1 = \frac{\Phi_1 - \Phi_2}{\alpha - 1} - \frac{\lambda_1 N_1 - \lambda_2 N_2}{\alpha - 1}$$
Substituting $I_1$ into $P_1 - \Phi_1$:
$$P_1 - \Phi_1 = 2 I_1 + M_{P,1} - \lambda_1 N_1 + \epsilon_{P,1}$$
$$P_1 - \Phi_1 - \frac{2}{\alpha - 1}(\Phi_1 - \Phi_2) = M_{P,1} + B_1 + \epsilon_{P,1}$$
where:
$$\beta = \frac{2}{\alpha - 1}$$
$$MP_1 = P_1 - (1 + \beta) \Phi_1 + \beta \Phi_2$$
$$B_1 = -\lambda_1 N_1 - \beta (\lambda_1 N_1 - \lambda_2 N_2)$$
**Properties of $MP_1$**:
1. Geometry $\rho$, vehicle dynamics, trajectory, and acceleration: **Eliminated (0.000m)**.
2. Satellite clock $dt^s$ and receiver clock $dt_r$: **Eliminated (0.000m)**.
3. Tropospheric delay $T$: **Eliminated (0.000m)**.
4. First-order ionospheric delay $I_1$: **Eliminated (0.000m)**.
5. $B_1$ is a rigorous mathematical constant along any continuous tracking arc (where $N_1, N_2$ do not change).
6. $MP_1$ fluctuates solely due to code multipath $M_{P,1}$ and code noise $\epsilon_{P,1}$.

### 2.3 Single-Frequency Code-Minus-Carrier ($CMC$)
When only a single frequency is tracked:
$$\text{CMC}_1 = P_1 - \Phi_1 = 2 I_1 - \lambda_1 N_1 + M_{P,1} + \epsilon_{P,1}$$
The time derivative over an epoch interval $\Delta t = t_k - t_{k-1}$ is:
$$\Delta \text{CMC}_{1, k} = \text{CMC}_{1, k} - \text{CMC}_{1, k-1} = \Delta M_{P,1} + 2 \Delta I_1 + \Delta\epsilon_{P,1}$$
The temporal rate of change of slant ionospheric delay under normal to moderate conditions is:
$$|\dot{I}_1| \le 1.0\text{--}5.0\text{ mm/s}$$
For $\Delta t = 1.0\text{ s}$, $2 \Delta I_1 \le 0.010\text{ m} = 1\text{ cm}$.
Even during severe ionospheric scintillation storms, $|2 \dot{I}_1| < 0.1\text{ m/s}$.
Therefore, an instantaneous jump of $5\text{m to }20\text{m}$ in $\text{CMC}_1$ is **$500\times$ to $2000\times$ larger than any possible ionospheric rate of change**!

### 2.4 Double-Difference Epoch-to-Epoch Delta CMC ($\delta\text{DD-CMC}$)
At the double-difference level between epoch $k-1$ and epoch $k$ along a continuous carrier arc:
$$\delta\text{CMC}_k = (\nabla\Delta P_k - \nabla\Delta P_{k-1}) - \lambda (\nabla\Delta \Phi_k - \nabla\Delta \Phi_{k-1})$$
Since ambiguity $\nabla\Delta N_k = \nabla\Delta N_{k-1}$ is constant across the arc:
$$\delta\text{CMC}_k = \Delta(\nabla\Delta M_P) + \Delta\epsilon_{\nabla\Delta P}$$
Because base station multipath is near zero (stationary geodetic antenna):
$$\Delta(\nabla\Delta M_P) \approx \Delta M_{P, \text{rov}}^s - \Delta M_{P, \text{rov}}^{\text{ref}}$$
Given that the reference satellite is selected at zenith (elevation $> 60^\circ$), its multipath is negligible, so:
$$\delta\text{CMC}_k \approx \Delta M_{P, \text{rov}}^s$$
This provides a direct, calibration-free detector of sudden pseudorange multipath steps!

---

## 3. Exact Algorithmic Design for 5m–20m Multipath Step Detection

### 3.1 Two-Stage Multipath Mitigation Architecture

To completely resolve the urban canyon fix rate bottleneck, we design a two-stage architecture:
- **Stage 1 (Arc-Level CMC Tracker)**: Operates per satellite and frequency band on raw rover observations (before DD formation). Computes dual-frequency $MP$ (or single-frequency $CMC$) and tracks the running baseline using an outlier-resistant recursive filter. Detects instantaneous step jumps (5m–20m) and sustained multipath biases.
- **Stage 2 (Observation Covariance Inflation & Filter Update Protection)**:
  1. Feeds the estimated code multipath variance $\sigma_{\text{mp}}^2$ into `formation_cov.rs` to inflate `pr_var_m2`.
  2. Carrier phase variance `cp_var_cycles2` remains completely unaltered.
  3. Modifies `screen_gross_pr_errors` (`screen.rs`) to suppress the code row rather than removing the entire `DoubleDiffMeasurement`.
  4. Protects `mw.rs` by gating MW updates when code multipath is detected.

```
       Rover Observations (P, Phi, fD)
                     │
                     ▼
         ┌───────────────────────┐
         │  Cycle Slip Detector  │◄── Doppler & GF Phase Checks
         └───────────┬───────────┘
                     │ (Arc IDs / Slip Flags)
                     ▼
         ┌───────────────────────┐
         │  Arc-Level CMC / MP   │
         │  Multipath Detector   │
         └───────────┬───────────┘
                     │
         ┌───────────┴───────────────────────┐
         │                                   │
         ▼ (sigma_mp^2, step_flag)           ▼ (freeze / skip)
┌─────────────────────────────┐    ┌─────────────────────┐
│ DD Formation & Covariance   │    │ Melbourne-Wübbena   │
│ pr_var inflated, cp_var OK  │    │ Tracker (mw.rs)     │
└──────────────┬──────────────┘    └─────────────────────┘
               │
               ▼
┌─────────────────────────────┐
│ Prefit Screening (screen.rs)│
│ Exclude code only, KEEP cp  │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│ Robust IEKF Update          │
│ Down-weighted Code Row      │
│ Full-Weight Carrier Row     │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│ AR & PAR Candidate Selection│
│ Clean carriers fixed first  │
└─────────────────────────────┘
```

### 3.2 Stage 1: Robust Arc-Level CMC Tracker Algorithm

For each satellite $s$ and frequency band $b$, maintain a tracker state:
```rust
pub struct CmcArcState {
    /// Running estimate of ambiguity + slow iono baseline (meters)
    pub baseline_m: f64,
    /// Running variance of nominal code noise (meters^2)
    pub nominal_var_m2: f64,
    /// Number of epochs accumulated in current continuous arc
    pub arc_length: u32,
    /// Instantaneous multipath deviation (meters)
    pub multipath_m: f64,
    /// Timestamp of last observation (TOW seconds)
    pub last_tow: f64,
    /// Flag indicating whether current epoch is contaminated by multipath
    pub is_multipath: bool,
}
```

#### Step 1: Form Observables
For satellite $s$ at epoch $k$:
- If dual-frequency carrier phases are available (e.g. L1 and L2/E5a/B2a):
  Compute $MP_1$:
  $$MP_1 = P_1 - (1 + \beta) \Phi_1 + \beta \Phi_2, \quad \beta = \frac{2}{(f_1/f_2)^2 - 1}$$
- If single-frequency:
  Compute $CMC_1$:
  $$CMC_1 = P_1 - \lambda_1 \Phi_1$$
Let $z_k$ denote the observable ($MP_1$ or $CMC_1$).

#### Step 2: Cycle Slip & Arc Continuity Validation
Check whether satellite $s$ experienced a cycle slip:
- If `slip_detected(s)` is true (Doppler slip, GF phase slip $> 5\text{ cm}$, or LLI bit 0):
  Reset tracker:
  $$\text{baseline} \leftarrow z_k, \quad \text{arc\_length} \leftarrow 1, \quad \text{multipath\_m} \leftarrow 0.0$$
- If $\Delta t = t_k - t_{k-1} > 2.0\text{ s}$:
  Reset tracker.

#### Step 3: Outlier-Resistant Baseline Update & Multipath Estimation
When $\text{arc\_length} < 5$ (initial acquisition):
- Update baseline with simple average:
  $$\text{baseline}_k = \text{baseline}_{k-1} + \frac{1}{\text{arc\_length}} (z_k - \text{baseline}_{k-1})$$
  $$\text{multipath\_m} = 0.0$$
When $\text{arc\_length} \ge 5$ (established arc):
- Compute deviation from baseline:
  $$\delta_k = z_k - \text{baseline}_{k-1}$$
- Evaluate against multipath detection thresholds:
  - $\sigma_{\text{code}} = \sqrt{\text{nominal\_var\_m2}}$ (typically $0.20\text{--}0.40\text{ m}$).
  - Multipath threshold: $\tau_{\text{mp}} = \max(3.0 \times \sigma_{\text{code}}, 2.5\text{ m})$.
  - Gross blunder threshold: $\tau_{\text{blunder}} = 15.0\text{ m}$.

**Case A: Clean Code ($|\delta_k| \le \tau_{\text{mp}}$)**:
- Signal is clean. Update baseline using Exponential Moving Average (EMA) or Welford with slow time constant ($\tau \approx 60\text{ s}$):
  $$\alpha = \frac{1}{\min(\text{arc\_length}, 100)}$$
  $$\text{baseline}_k = \text{baseline}_{k-1} + \alpha \cdot \delta_k$$
  $$\text{multipath\_m} = 0.0, \quad \text{is\_multipath} = \text{false}$$

**Case B: Multipath Detected ($\tau_{\text{mp}} < |\delta_k| \le \tau_{\text{blunder}}$)**:
- Carrier phase is confirmed clean by Doppler & GF detectors, but code has jumped by $2.5\text{m to }15\text{m}$.
- **Do NOT update baseline with contaminated code!** (Freezing the baseline prevents baseline corruption during reflections).
  $$\text{baseline}_k = \text{baseline}_{k-1}$$
  $$\text{multipath\_m} = |\delta_k|, \quad \text{is\_multipath} = \text{true}$$

**Case C: Severe Multipath / NLOS ($|\delta_k| > \tau_{\text{blunder}}$)**:
- Code has jumped by $> 15.0\text{m}$.
- Baseline remains frozen:
  $$\text{baseline}_k = \text{baseline}_{k-1}$$
  $$\text{multipath\_m} = |\delta_k|, \quad \text{is\_multipath} = \text{true}, \quad \text{is\_blunder} = \text{true}$$

---

## 4. Integration & Interaction with IEKF Robust Update

### 4.1 Down-Weighting & Variance Inflation in `formation_cov.rs`

When forming double differences in `formation.rs`, query the CMC tracker for the multipath estimates of both satellite $s$ and reference satellite $r$:
$$\sigma_{\text{mp}, s}^2 = (\text{multipath\_m}_s)^2$$
$$\sigma_{\text{mp}, r}^2 = (\text{multipath\_m}_r)^2$$
In `formation_cov.rs`, update `compute_dd_variances`:
```rust
pub fn compute_dd_variances_with_cmc(
    rx_pos: Vector3<f64>,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
    lambda: f64,
    snrs: (Option<u8>, Option<u8>, Option<u8>, Option<u8>),
    sat_mp_m: f64,
    ref_mp_m: f64,
) -> DdVariance {
    let mut var = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, snrs);
    
    // Inflate pseudorange variance by code multipath squared
    let mp_var = sat_mp_m * sat_mp_m + ref_mp_m * ref_mp_m;
    var.pr_var_m2 += mp_var;
    var.pr_ref_var_m2 += ref_mp_m * ref_mp_m;

    // CARRIER PHASE VARIANCE REMAINS ABSOLUTELY UNTOUCHED!
    // cp_var_cycles2 and cp_ref_var_cycles2 stay at nominal mm-level
    var
}
```

### 4.2 Mathematical Proof: Preservation of Carrier Phase Weights & Positive Definiteness

Let the measurement vector for a DD pair be partitioned into code and phase:
$$\mathbf{y} = \begin{bmatrix} y_P \\ y_\Phi \end{bmatrix}, \quad H = \begin{bmatrix} H_P \\ H_\Phi \end{bmatrix}, \quad R = \begin{bmatrix} R_{PP} & \mathbf{0} \\ \mathbf{0} & R_{\Phi\Phi} \end{bmatrix}$$
where:
- $H_P = \begin{bmatrix} \mathbf{u}^T & \mathbf{0} & 1 \end{bmatrix}$ (geometry, iono)
- $H_\Phi = \begin{bmatrix} \frac{1}{\lambda}\mathbf{u}^T & 1 & -\frac{1}{\lambda} \end{bmatrix}$ (geometry, ambiguity, iono)
- $R_{PP} = \sigma_{P, \text{nom}}^2 + \sigma_{\text{mp}}^2$
- $R_{\Phi\Phi} = \sigma_\Phi^2 \approx 10^{-4}\text{ cycles}^2$

#### Kalman Gain Derivation
The innovation covariance is:
$$S = H P_0 H^T + R = \begin{bmatrix} H_P P_0 H_P^T + R_{PP} & H_P P_0 H_\Phi^T \\ H_\Phi P_0 H_P^T & H_\Phi P_0 H_\Phi^T + R_{\Phi\Phi} \end{bmatrix}$$
As code multipath $\sigma_{\text{mp}} \to \infty$ ($R_{PP} \gg 100\text{ m}^2$):
1. In $S^{-1}$, the block corresponding to code scales as $S_{PP}^{-1} \sim \frac{1}{R_{PP}} \to 0$.
2. The Kalman gain block for code:
   $$K_P = P_0 H_P^T S_{PP}^{-1} \to \mathbf{0}$$
   The corrupted pseudorange measurement has zero influence on the position and ambiguity updates!
3. The carrier phase block $H_\Phi P_0 H_\Phi^T + R_{\Phi\Phi}$ remains completely dominant. The Kalman gain for carrier phase:
   $$K_\Phi = P_0 H_\Phi^T (H_\Phi P_0 H_\Phi^T + R_{\Phi\Phi})^{-1}$$
   Carrier phase observations pull the filter states with full millimeter precision.
4. **Positive Definiteness**:
   Since $R_{PP} \ge \sigma_{P, \text{nom}}^2 > 0$ and $R_{\Phi\Phi} > 0$, and off-diagonal covariance blocks between code and phase are strictly zero, $R$ is strictly diagonally dominant:
   $$R \succ 0$$
   Therefore $S = H P_0 H^T + R \succ 0$ is strictly positive definite and its Cholesky decomposition is guaranteed to succeed.

### 4.3 Refactoring `screen_gross_pr_errors` (`screen.rs`)
To eliminate Critical Defect 1:
Instead of removing the entire `DoubleDiffMeasurement` when prefit residual $|y_P| > 15.0\text{m}$, decouple code and phase:
```rust
/// Screen gross pseudorange blunders without discarding valid carrier phase.
pub fn screen_gross_pr_errors_retaining_phase(
    measurements: &mut [DoubleDiffMeasurement],
    pos_pred: Vector3<f64>,
) -> Vec<DoubleDiffKey> {
    let mut suppressed_code = Vec::new();
    for m in measurements.iter_mut() {
        let res = pr_residual(m, pos_pred);
        if res.abs() > GROSS_PR_ERROR_THRESHOLD_M {
            // Inflate code variance to infinity (or flag to omit code row)
            // DO NOT drop the measurement: retain m.dd_cp_cycles!
            m.pr_var_m2 = 1.0e8; // Effectively drops code row in Kalman gain
            suppressed_code.push(m.key);
        }
    }
    suppressed_code
}
```
In `update/system.rs::append_dd_code_row`:
If `m.pr_var_m2 >= 1.0e6` or `is_code_blunder` triggers, the code row is simply omitted (`return;`), while `append_dd_phase_row` continues to execute normally.

### 4.4 Protecting Melbourne–Wübbena in `mw.rs`
In `mw.rs::update_tracker_from_obs`:
Before calling `tracker.update(...)`:
Check if satellite $s$ has `is_multipath == true` ($\text{multipath\_m} > 2.5\text{ m}$):
- If `is_multipath` is true:
  **Skip the MW update for this epoch.**
  Do not absorb the contaminated measurement into `MwTrack`.
  This prevents the 10m pseudorange multipath jump from injecting a 6.5-cycle blunder and tripping the false-slip reset!
- The continuous wide-lane tracking arc survives uninterrupted through the urban canyon reflection.

---

## 5. Ambiguity Resolution (AR) Protection & False Fix Prevention

### 5.1 Shielding Ambiguity Seeding (`mod.rs`)
In `GnssRtkIekf::update_dd_ambiguity`:
When an ambiguity must be freshly seeded:
```rust
let init_amb = dd_cp.map_or(0.0, |cp| cp - dd_pr / lambda);
```
If the satellite has `multipath_m > 2.5\text{ m}`:
1. De-bias the pseudorange: $dd\_pr_{\text{clean}} = dd\_pr - \text{multipath\_m} \cdot \text{sign}(y_P)$, OR
2. Inflate the initial ambiguity variance:
   ```rust
   let seed_var = seed_ambiguity_variance_cycles2(m.pr_var_m2, lambda);
   ```
   Because `m.pr_var_m2` has been inflated by $\sigma_{\text{mp}}^2$, `seed_var` increases from $70\text{ cycles}^2$ to $50,000\text{ cycles}^2$. The filter does not falsely constrain the ambiguity to the corrupted code measurement.

### 5.2 Shielding Phase-Code Coherence Bias (`ar_gate.rs`)
In `ar_gate::coherence_offset`:
```rust
pub(crate) fn coherence_offset(band: u8, divergences: &[(u8, f64)]) -> f64 {
    let same_band: Vec<f64> = divergences.iter()
        .filter(|(b, _)| *b == band)
        .map(|(_, d)| *d)
        .collect();
    median(&same_band)
}
```
Satellites currently flagged with code multipath (`is_multipath == true`) must be **excluded** from `code_phase_div` so they cannot skew the median coherence offset for freshly rising satellites.

### 5.3 SNR- and Multipath-Prioritized Partial Ambiguity Resolution (PAR)
In `ar_subsets.rs` and `par.rs`:
When full-set LAMBDA fails the FFRT ratio test ($R < \mu_{\text{req}}$):
PAR subsets must be ranked not merely by index, but by tracking quality:
$$\text{Quality}(i) = \text{LockEpochs}(i) \cdot \sin(\theta_i) \cdot \frac{1}{1.0 + \text{multipath\_m}_i}$$
Ambiguities on satellites with active code multipath are excluded from the candidate fixed subset first.
This guarantees that **only clean, uncompromised carrier phase ambiguities enter the integer search**, completely preventing false integer fixes while maximizing fix rate.

---

## 6. Implementation Architecture & File Layout

Following `AGENTS.md` and repository modularity:
- All new files $< 500$ LOC.
- All functions $\le 32$ LOC.
- Nesting depth $< 3$ levels.
- Zero `unwrap()` calls in production code.
- Zero compiler/clippy warnings.

### 6.1 Proposed New Modules & Placement
1. **`crates/gneiss-rtk/src/measurements/cmc.rs`** (~220 LOC):
   - `CmcArcTracker`: Manages per-satellite, per-band continuous arc tracking.
   - Implements `push_sample(sat, band, pr, cp, cp2_opt, lambda1, lambda2_opt, lli_slip, tow)`.
   - Computes dual-frequency $MP$ combination or single-frequency $CMC$.
   - Provides `multipath_m(sat, band) -> f64` and `is_multipath(sat, band) -> bool`.
2. **`crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`**:
   - Add `compute_dd_variances_with_cmc` (or add optional multipath variance arguments to `compute_dd_variances`).
3. **`crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`**:
   - Refactor `screen_gross_pr_errors` to flag/inflate code variance without removing the `DoubleDiffMeasurement` (retaining `dd_cp_cycles`).
4. **`crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`**:
   - Guard `update_tracker_from_obs` against code multipath contamination.

### 6.2 Three-Tier Verification Plan
1. **Tier 1 (Analytical Golden Vectors)**:
   - Synthetic constant code multipath ($5\text{m}, 10\text{m}, 20\text{m}$ step injections) with clean carrier phase.
   - Verify that `multipath_m` converges within 1–2 epochs of the step.
   - Verify that carrier phase variance remains exact nominal ($0.003\text{ m} / \lambda$).
   - Verify that a 1-cycle carrier phase slip is NOT flagged as code multipath (correctly delegated to Doppler / GF slip detector).
2. **Tier 2 (Filter Covariance Invariants)**:
   - Verify that Joseph covariance update with inflated code variance preserves $P \succ 0$ and $Q_{aa} \succ 0$.
   - Verify that Kalman gain for the multipath satellite code row drops by $\ge 99\%$, while phase Kalman gain is preserved.
3. **Tier 3 (Benchmark Regression)**:
   - Run multi-GNSS benchmark smoke script: `python3 scripts/check_multignss_benchmark.py --smoke` (must pass with `ALL CHECKS PASSED`).
   - Run network benchmark smoke script: `python3 scripts/check_network_benchmark.py --smoke` (must pass with `ALL CHECKS PASSED`).
   - Run UrbanNav benchmark (`eval_f9p_rover` on Shinjuku, Whampoa, TST1): assert fix rate expands and $p_{95}$ tail error collapses with zero false fixes.

---

## 7. Conclusion

Code-Minus-Carrier (CMC) multipath detection is the decisive missing link in Gneiss's urban canyon navigation performance. By exploiting the physical asymmetry between carrier phase multipath ($\le 4.8\text{ cm}$) and pseudorange multipath ($5\text{m--}20\text{m}$), the proposed architecture:
1. Detects code multipath steps instantaneously without false-flagging clean carrier phase.
2. Down-weights code measurements without discarding carrier phase in prefit screening.
3. Protects Melbourne–Wübbena from false slip resets.
4. Prevents float ambiguity corruption and unlocks high fix rates in skyscraper canyons.
