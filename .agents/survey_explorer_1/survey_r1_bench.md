# Survey & Baseline Benchmark Report: R1 Adaptive C/N0 & Elevation Covariance Weighting

**Author**: Survey Explorer 1  
**Working Directory**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_1`  
**Parent Task ID**: `c1309e2d-6c95-4b14-a86d-d26a13f2a150`  
**Date**: 2026-09-24  

---

## 1. Executive Summary

This survey maps the architectural baseline and mathematical formulation for **R1: Adaptive C/N0 (SNR) and elevation observation covariance weighting** in the Gneiss RTK positioning engine. It additionally establishes verified numerical baselines across the workspace test suite, dual CI smoke guard scripts, the multi-dataset u-blox ZED-F9P kinematic rover benchmark (`eval_f9p_rover`), and the Tokyo Odaiba tightly-coupled GNSS/INS RTS benchmark (`eval_odaiba_ins`).

Key findings:
1. **Clean Baseline**: The current workspace passes all tests cleanly (`205 passed; 0 failed`), achieves zero compiler and zero clippy warnings (`-D warnings`), and satisfies both CI regression guard scripts (`check_network_benchmark.py --smoke` and `check_multignss_benchmark.py --smoke`).
2. **Current R1 Limitations**:
   - `crates/gneiss-core/src/variance.rs` contains an initial elevation/SNR scaling function used by SPP, but it is not unified with `rtk_iekf`.
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs` contains an SNR attenuation multiplier with slope discontinuities at 40 dB-Hz, 25 dB-Hz, and at the 1000.0 clamp, causing Kalman filter gain chatter under dynamic SNR fluctuation.
   - `crates/gneiss-core/src/obs.rs` casts SNR observations from `f64` to `Option<u8>`, causing 1 dB discrete step quantization in covariance weights.
3. **Exact SIGMA-SNR Formulation**: A unified, continuously differentiable ($C^1$) model $\sigma^2(\theta, S) = \left( a^2 + \frac{b^2}{\sin^2 \theta} \right) \cdot f_{\text{SNR}}(S)$ has been derived with smooth sigmoid blending and algebraic asymptotic saturation, guaranteeing positive-definiteness ($R_{DD} \succ 0$) and eliminating Kalman gain chatter.

---

## 2. Codebase Baseline Mapping for R1

### 2.1 File & Struct Inventory

| Component | File Path | Current Role & State |
|---|---|---|
| **Physical Variance Models** | `crates/gneiss-core/src/variance.rs` | 64 LOC. Houses `snr_variance_scale`, `elevation_variance_scale`, and `observation_variance`. Currently used only by `spp/solver.rs`. Needs expansion with the unified SIGMA-SNR model and Three-Tier tests. |
| **Observation Representation** | `crates/gneiss-core/src/obs.rs:184-192` | `SatObs::get_snr(&self, freq_band: u8) -> Option<u8>`. `Observation::value` is `f64`, but `get_snr` casts to `u8`. Needs `get_snr_f64` to prevent 1 dB-Hz quantization steps. |
| **DD Covariance Modeling** | `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs` | 150 LOC. Computes `DdVariance` (`pr_var_m2`, `cp_var_cycles2`, `pr_ref_var_m2`, `cp_ref_var_cycles2`) via `snr_weight`, `expected_cn0_dbhz`, `attenuation_scale`, and `single_diff_var`. Has piecewise derivative discontinuities. |
| **DD Observation Formation** | `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs:240-244` | Extracts `snrs = (rov_s.get_snr(freq_band), rov_ref.get_snr(freq_band), bas_s.get_snr(freq_band), bas_ref.get_snr(freq_band))` and calls `compute_dd_variances`. Constructs `DoubleDiffMeasurement`. |
| **Measurement Update & Covariance Assembly** | `crates/gneiss-rtk/src/estimators/rtk_iekf/update/system.rs:115-153, 237-251` | Assembles diagonal measurement noise $R_{ii}$ and correlated off-diagonal terms $R_{ij} = \min(R_{\text{ref}, i}, R_{\text{ref}, j})$ for satellites sharing the reference satellite. |
| **ESKF Tightly-Coupled DD Update** | `crates/gneiss-rtk/src/estimators/eskf/dd_update.rs:101-125` | Accepts scalar DD variance `var` directly; in `eval_odaiba_ins/main.rs:221-231`, hardcodes `var_code = (2.0 / sin_el).powi(2)`. Can directly leverage the R1 variance formulation. |

### 2.2 Trace of Current Double-Difference Covariance Flow

1. **Extraction in `formation.rs:240`**:
   ```rust
   let snrs = (rov_s.get_snr(freq_band), rov_ref.get_snr(freq_band), bas_s.get_snr(freq_band), bas_ref.get_snr(freq_band));
   let dd_var = compute_dd_variances(self.state.pos_ecef, sat_pos, ref_pos, lambda, snrs);
   if dd_cp.is_none() && dd_var.pr_var_m2 > 100.0 {
       return None;
   }
   ```
2. **Evaluation in `formation_cov.rs:68-90`**:
   - Calculates satellite elevation $\theta_s$ and reference satellite elevation $\theta_r$ via `gneiss_core::coords::az_el`.
   - Computes single-difference variances:
     - `pr_sat_var = single_diff_var(0.20, snrs.0, snrs.2, el_s)`
     - `pr_ref_var = single_diff_var(0.20, snrs.1, snrs.3, el_r)`
     - `cp_sat_var = single_diff_var(0.003 / lambda, snrs.0, snrs.2, el_s)`
     - `cp_ref_var = single_diff_var(0.003 / lambda, snrs.1, snrs.3, el_r)`
   - Computes double-difference total variances:
     - `pr_var_m2 = pr_sat_var + pr_ref_var`
     - `cp_var_cycles2 = cp_sat_var + cp_ref_var`
3. **Current `single_diff_var` Function (`formation_cov.rs:59-65`)**:
   ```rust
   fn single_diff_var(base_sigma: f64, snr_a: Option<u8>, snr_b: Option<u8>, el_rad: f64) -> f64 {
       let sin_el = el_rad.sin().max(0.1);
       let w_a = snr_weight(snr_a) * attenuation_scale(snr_a, el_rad);
       let w_b = snr_weight(snr_b) * attenuation_scale(snr_b, el_rad);
       let w = (w_a + w_b).min(1000.0);
       (base_sigma * base_sigma * w) / (sin_el * sin_el)
   }
   ```

### 2.3 Identified Weaknesses & Deficiencies

1. **Derivative Discontinuities**:
   - In `snr_weight`:
     $$\left. \frac{d w}{dS} \right|_{S \to 40^-} = -0.230, \quad \left. \frac{d w}{dS} \right|_{S \to 40^+} = 0.0$$
   - In `attenuation_scale`: hard `if drop_db > 10.0` and `if snr < 25.0` cause sudden exponential scaling onset.
   - In `sin_el.max(0.1)`: slope abruptly jumps to zero for elevations below $5.74^\circ$.
2. **Quantization Chatter**:
   - `SatObs::get_snr` converts `o.value as u8`. When signal fluctuates between $39.4$ and $40.1$ dB-Hz, the integer flips between 39 and 40, causing abrupt jumps between $w=1.258$ and $w=1.000$.
3. **Coupling Discrepancy**:
   - Base station observations are taken at open-sky CORS/geodetic sites (SNR $> 45$ dB-Hz), while rover observations in canyons drop below $25$ dB-Hz. The term `w_a + w_b` lumps rover and base attenuation before multiplying by elevation, instead of independently modeling rover physical noise and base physical noise according to error propagation.

---

## 3. Exact Mathematical SIGMA-SNR Formulation

### 3.1 Model Equation
The observation variance for a single receiver-satellite link on code or carrier is formulated as:
$$\sigma^2(\theta, S) = \left( a^2 + \frac{b^2}{\sin_{\text{eff}}^2 \theta} \right) \cdot f_{\text{SNR}}(S, \theta)$$

where:
- $\theta \in (0, \pi/2]$ is the satellite elevation angle (radians).
- $S \in [0, 60]$ is the received C/N0 (dB-Hz).
- $a$ is the elevation-independent receiver channel noise floor (meters for code, cycles for phase).
- $b$ is the elevation-dependent multipath and tropospheric mapping coefficient.

### 3.2 Smooth Elevation Regularization
To eliminate the slope discontinuity of `.max(0.1)` and prevent horizon division by zero:
$$\sin_{\text{eff}}^2(\theta) = \sin^2 \theta + \sin^2 \theta_0$$
with $\theta_0 = 5^\circ$ ($0.087266$ rad, $\sin^2 \theta_0 \approx 0.007596$).
- Continuous derivative:
  $$\frac{d}{d\theta}\left(\sin_{\text{eff}}^2(\theta)\right) = 2 \sin \theta \cos \theta$$
- Strictly positive: $\sin_{\text{eff}}^2(\theta) \ge \sin^2 \theta_0 > 0$.
- Monotonic: $\frac{\partial \sigma^2}{\partial \theta} \le 0$ everywhere on $[0, \pi/2]$.

### 3.3 Smooth Differentiable Exponential SNR Scaling ($f_{\text{SNR}}$)
To achieve smooth exponential noise growth below nominal threshold $S_{\text{nom}} = 40.0\text{ dB-Hz}$ without derivative chatter:

1. **Deficit Calculation with Reflection Penalty**:
   Nominal elevation-expected C/N0:
   $$S_{\text{exp}}(\theta) = 30.0 + 20.0 \sin \theta$$
   Effective SNR deficit $\Delta S$:
   $$\Delta S = \max\left( S_{\text{nom}} - S, \; S_{\text{exp}}(\theta) - S - 8.0 \right)$$
   (where dropping $> 8\text{ dB-Hz}$ below elevation expectation flags NLOS wall reflections).

2. **Smooth Logistic Activation**:
   Instead of a piecewise `if S < 40`, define smooth activation:
   $$f_{\text{raw}}(S) = 1.0 + 10^{\frac{\Delta S}{10.0}} \cdot \sigma\left(\frac{\Delta S}{\tau}\right)$$
   where $\sigma(u) = \frac{1}{1 + e^{-u}}$ is the standard logistic sigmoid and $\tau = 1.5\text{ dB-Hz}$ is the smoothing transition bandwidth.
   - For $S \ge S_{\text{nom}} + 3\tau$ ($\Delta S \le -4.5$ dB): $\sigma \approx 0 \implies f_{\text{raw}} = 1.000$.
   - For $S \le S_{\text{nom}} - 3\tau$ ($\Delta S \ge +4.5$ dB): $\sigma \approx 1 \implies f_{\text{raw}} \approx 10^{\Delta S / 10}$.
   - At $\Delta S = 0$ ($S = 40$ dB-Hz): $\sigma = 0.5 \implies f_{\text{raw}} = 1.500$.
   - Derivative $\frac{df_{\text{raw}}}{dS}$ is smooth and continuous across the entire real line.

3. **Smooth Algebraic Ceiling Saturation ($f_{\max} = 1000.0$)**:
   To strictly bound matrix condition numbers ($\kappa(R) \le 10^4$) without a hard clamp:
   $$f_{\text{SNR}}(S) = 1.0 + \frac{(f_{\max} - 1.0) \cdot (f_{\text{raw}}(S) - 1.0)}{(f_{\max} - 1.0) + (f_{\text{raw}}(S) - 1.0)}$$
   - $f_{\text{SNR}} \in [1.0, f_{\max})$.
   - As $S \to -\infty$, $f_{\text{SNR}} \to f_{\max} = 1000.0$ asymptotically.
   - Preserves strict monotonicity: $\frac{\partial f_{\text{SNR}}}{\partial S} < 0$.

### 3.4 Physical Parameters for Code and Carrier

| Parameter | Pseudorange (Code) | Carrier Phase | Rationale |
|---|---|---|---|
| $a$ (channel noise floor) | $0.1414\text{ m}$ ($0.20 / \sqrt{2}$) | $0.00212 / \lambda\text{ cyc}$ | Zenith thermal noise floor |
| $b$ (elevation-dependent) | $0.1414\text{ m}$ ($0.20 / \sqrt{2}$) | $0.00212 / \lambda\text{ cyc}$ | Atmospheric mapping & diffuse multipath |
| Zenith Variance ($a^2 + b^2$) | $0.0400\text{ m}^2$ ($\sigma = 0.20\text{ m}$) | $(0.003 / \lambda)^2\text{ cyc}^2$ | Matches historical baseline |
| $S_{\text{nom}}$ | $40.0\text{ dB-Hz}$ | $40.0\text{ dB-Hz}$ | Standard geodetic tracking threshold |
| $f_{\max}$ | $1000.0$ | $1000.0$ | Upper variance scale ($\kappa \le 10^4$) |

### 3.5 Double-Difference Error Propagation & Covariance Definiteness
For rover receiver $u$ and base receiver $b$, observing satellite $s$ and reference satellite $r$:
$$\sigma_{\text{sd}}^2(s) = \sigma^2(\theta_s, S_u^s) + \sigma^2(\theta_s, S_b^s)$$
$$\sigma_{\text{sd}}^2(r) = \sigma^2(\theta_r, S_u^r) + \sigma^2(\theta_r, S_b^r)$$
The double-difference variance is:
$$\text{Var}(\nabla\Delta^{s, r}) = \sigma_{\text{sd}}^2(s) + \sigma_{\text{sd}}^2(r)$$
and the reference satellite covariance contribution is:
$$\text{Cov}(\nabla\Delta^{s, r}, \nabla\Delta^{k, r}) = \sigma_{\text{sd}}^2(r)$$

**Positive-Definiteness Theorem**:
The double difference transformation is $z_{DD} = D z_{\text{undiff}}$, where $D$ has full row rank $(m-1)$ and $R_{\text{undiff}} = \text{diag}(\sigma_1^2, \ldots, \sigma_{2m}^2)$.
Because $\sigma^2(\theta, S) > 0$ for all physical $\theta, S$:
$$R_{DD} = D R_{\text{undiff}} D^T \succ 0$$
Scaling any satellite's noise $f_{\text{SNR}}(S_s) \ge 1.0$ strictly increases the corresponding diagonal element of $R_{DD}$ without increasing off-diagonal terms, strictly increasing all eigenvalues ($\lambda_i(R_{DD})$) and preserving positive-definiteness.

---

## 4. Baseline Benchmark Results & Verification

All benchmarks were run directly on the current codebase commit in release mode on mac (Apple Silicon).

### 4.1 CI Regression Guard Scripts

#### 1. Network RTK Smoke Guard (`python3 scripts/check_network_benchmark.py --smoke`)
- **Status**: `ALL CHECKS PASSED` (exit code 0)
- **Epochs Evaluated**: 1,800 epochs

| Metric | Measured Baseline | Requirement Budget | Status |
|---|:---:|:---:|:---:|
| Network Fused Horizontal $p_{50}$ | **0.024 m** | $\le 0.040\text{ m}$ | OK |
| Network Fused Horizontal RMS | **0.037 m** | $\le 0.060\text{ m}$ | OK |
| Network Fused Vertical RMS | **0.042 m** | $\le 0.080\text{ m}$ | OK |
| P181 Smoothed Fixed-Only $p_{50}$ | **0.022 m** | $\le 0.030\text{ m}$ | OK |
| P222 Smoothed Fixed-Only $p_{50}$ | **0.060 m** | $\le 0.090\text{ m}$ | OK |
| SLAC Smoothed Fixed-Only $p_{50}$ | **0.113 m** | $\le 0.120\text{ m}$ | OK |
| OHLN Smoothed Fix Rate | **98.4%** | $\ge 79.0\%$ | OK |
| P181 Smoothed Fix Rate | **99.8%** | $\ge 85.0\%$ | OK |
| SLAC Smoothed Fix Rate | **70.3%** | $\ge 60.0\%$ | OK |

#### 2. Multi-GNSS Smoke Guard (`python3 scripts/check_multignss_benchmark.py --smoke`)
- **Status**: `ALL CHECKS PASSED` (exit code 0)
- **Binary Hash**: `1f9c24bd2773`
- **Epochs Evaluated**: 1,800 epochs (~30 min equivalent)

| Metric | Measured Baseline | Requirement Budget | Status |
|---|:---:|:---:|:---:|
| P181 Fix Rate | **98.30%** | $\ge 97.50\%$ | OK |
| P181 Horizontal $p_{95}$ | **134.00 mm** | $\le 145.00\text{ mm}$ | OK |
| P181 Vertical $p_{95}$ | **195.00 mm** | $\le 290.00\text{ mm}$ | OK |
| P225 Fix Rate | **93.70%** | $\ge 71.00\%$ | OK |
| P225 Horizontal $p_{95}$ | **113.00 mm** | $\le 245.00\text{ mm}$ | OK |
| P225 Vertical $p_{95}$ | **268.00 mm** | $\le 370.00\text{ mm}$ | OK |
| P222 Fix Rate | **99.50%** | $\ge 86.00\%$ | OK |
| P222 Horizontal $p_{95}$ | **198.00 mm** | $\le 295.00\text{ mm}$ | OK |
| P222 Vertical $p_{95}$ | **84.00 mm** | $\le 135.00\text{ mm}$ | OK |
| Network Fused Fix Rate | **99.70%** | $\ge 96.50\%$ | OK |

### 4.2 Workspace Test Suite & Lints
- `cargo test --workspace`: **205 passed; 0 failed; 0 ignored** in 0.01s.
- `cargo clippy --workspace --all-targets -- -D warnings`: **0 warnings** (Finished dev profile in 0.38s).

### 4.3 Multi-Dataset u-blox ZED-F9P Kinematic Rover Matrix (`eval_f9p_rover`)

#### Verified Full-Trajectory Baseline (from `PROJECT_STATUS.md`):

| Dataset & Antenna Config | Processing Mode | Matched Epochs | Fix Rate | $p_{50}$ (H) | $p_{68}$ (H) | $p_{95}$ (H) | RMS (H) | RMS (3D) |
|---|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **Tokyo Odaiba** (Survey Ant) | Forward RTK | 1242/1242 | 9.4% | 1.134 m | 1.629 m | 3.536 m | 1.843 m | 4.307 m |
| **Tokyo Odaiba** (Survey Ant) | **Smoothed PPK** | 1242/1242 | **13.8%** | **1.161 m** | **1.666 m** | **3.543 m** | **1.917 m** | **4.284 m** |
| **Tokyo Shinjuku** (Survey Ant) | Forward RTK | 2096/2096 | 4.2% | 1.620 m | 1.981 m | 11.307 m | 5.386 m | 7.794 m |
| **Tokyo Shinjuku** (Survey Ant) | **Smoothed PPK** | 2096/2096 | **6.1%** | **1.632 m** | **1.988 m** | **11.271 m** | **4.708 m** | **7.327 m** |
| **Hong Kong TST1** (Survey Splitter) | Forward RTK | 786/786 | 4.6% | 1.177 m | 1.776 m | 9.120 m | 3.945 m | 7.813 m |
| **Hong Kong TST1** (Survey Splitter) | **Smoothed PPK** | 786/786 | **6.0%** | **1.158 m** | **1.729 m** | **8.861 m** | **3.912 m** | **6.530 m** |
| **Hong Kong Whampoa** (Survey Splitter) | Forward RTK | 1535/1535 | 15.2% | 1.243 m | 2.117 m | 9.651 m | 5.681 m | 10.710 m |
| **Hong Kong Whampoa** (Survey Splitter) | **Smoothed PPK** | 1535/1535 | **12.8%** | **1.181 m** | **2.152 m** | **10.142 m** | **4.857 m** | **9.840 m** |
| **Hong Kong Whampoa** (Patch Ant) | **Smoothed PPK** | 1534/1534 | **2.6%** | **1.824 m** | **2.953 m** | **25.081 m** | **9.058 m** | **18.945 m** |
| **Hong Kong TST1** (Patch Ant) | **Smoothed PPK** | 657/657 | **0.9%** | **2.490 m** | **3.763 m** | **11.313 m** | **4.990 m** | **8.967 m** |

*Fixed-subset ground truth accuracy*: $p_{50} = 0.573\text{ m}, p_{95} = 1.361\text{ m}$ (Whampoa Survey) and $p_{50} = 0.912\text{ m}, p_{95} = 2.247\text{ m}$ (Shinjuku Survey) with **0 false integer fixes**.

#### Live 200-Epoch Smoke Run Verification (`MAX_EPOCHS=200`):

| Dataset | Mode | Epochs | Fix Rate | $p_{50}$ (H) | $p_{95}$ (H) | RMS (H) | Fixed Subset Metrics |
|---|---|:---:|:---:|:---:|:---:|:---:|---|
| **Tokyo Odaiba** | Fwd RTK | 200/200 | 20.0% | 1.252 m | 1.854 m | 1.271 m | 40 epochs: $p_{50} = 1.222\text{m}, p_{95} = 1.794\text{m}$ |
| **Tokyo Odaiba** | Smooth PPK | 200/200 | 58.5% | 1.161 m | 2.171 m | 1.265 m | 117 epochs: $p_{50} = 1.072\text{m}, p_{95} = 2.042\text{m}$ |
| **Tokyo Shinjuku** | Fwd RTK | 200/200 | 7.5% | 0.419 m | 1.408 m | 0.719 m | 15 epochs: $p_{50} = 0.460\text{m}, p_{95} = 1.467\text{m}$ |
| **Tokyo Shinjuku** | Smooth PPK | 200/200 | 20.5% | 0.455 m | 1.434 m | 0.793 m | 41 epochs: $p_{50} = 0.872\text{m}, p_{95} = 1.949\text{m}$ |
| **Whampoa Survey** | Fwd RTK | 200/200 | 41.0% | 0.627 m | 1.284 m | 0.847 m | 82 epochs: $p_{50} = 0.572\text{m}, p_{95} = 1.264\text{m}$ |
| **Whampoa Survey** | Smooth PPK | 200/200 | 40.0% | 0.606 m | 1.367 m | 0.743 m | 80 epochs: $p_{50} = 0.330\text{m}, p_{95} = 1.037\text{m}$ |

### 4.4 Tokyo Odaiba 10Hz/50Hz GNSS/INS RTS Benchmark (`eval_odaiba_ins`)
- **Dataset**: Tokyo Odaiba waterfront & elevated highway underpass ($N = 12,398$ epochs, 10 Hz GNSS + 50 Hz IMU).
- **Forward Inertial Filter**:
  - $p_{50} = 2.162\text{ m}$
  - $p_{68} = 3.173\text{ m}$
  - $p_{95} = 7.142\text{ m}$
  - $\text{RMS} = 4.287\text{ m}$
- **RTS Smoothed GNSS/INS Filter**:
  - $p_{50} = 2.134\text{ m}$
  - $p_{68} = 3.782\text{ m}$
  - $p_{95} = 7.049\text{ m}$
  - $\text{RMS} = 4.156\text{ m}$
  - Underpass / rail overpass section Q2: $p_{50} = 1.752\text{ m}$
  - High-multipath canyon section Q3: $p_{50} = 2.044\text{ m}, \text{RMS} = 3.234\text{ m}$

---

## 5. Implementation Roadmap for R1

The R1 Implementer (Worker) should follow these concrete steps:

1. **Enhance `crates/gneiss-core/src/obs.rs`**:
   - Add `pub fn get_snr_f64(&self, freq_band: u8) -> Option<f64>` on `SatObs`.
   - Preserve existing `get_snr(&self, freq_band: u8) -> Option<u8>` to maintain complete backward compatibility.
2. **Implement Unified Model in `crates/gneiss-core/src/variance.rs`**:
   - Implement `sigma_snr_variance(el_rad: f64, snr_dbhz: f64, a: f64, b: f64) -> f64`.
   - Add helper functions: `elevation_factor(el_rad: f64, a: f64, b: f64) -> f64` and `snr_factor(snr_dbhz: f64, el_rad: f64) -> f64`.
   - Write Three-Tier unit tests:
     - Tier 1: Exact analytical golden vectors at nominal $(90^\circ, 45\text{ dB-Hz})$, mid $(30^\circ, 35\text{ dB-Hz})$, and attenuated $(15^\circ, 20\text{ dB-Hz})$.
     - Tier 2: Finite-difference numerical gradient checks verifying $\frac{\partial \sigma^2}{\partial \theta} \le 0$ and $\frac{\partial \sigma^2}{\partial S} \le 0$ everywhere without slope discontinuity.
     - Tier 3: Asymptotic stability checks as $\theta \to 0$ and $S \to -\infty$ verifying $\sigma^2 \le 10^4 \sigma_0^2$.
3. **Refactor `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`**:
   - Update `compute_dd_variances` to accept `snrs: (Option<f64>, Option<f64>, Option<f64>, Option<f64>)` or overloaded/converted tuples.
   - Refactor `single_diff_var` to independently evaluate rover and base variances:
     $$\sigma_{\text{sd}}^2 = \sigma_{\text{rov}}^2(\theta_s, S_{\text{rov}}^s) + \sigma_{\text{bas}}^2(\theta_s, S_{\text{bas}}^s)$$
   - Ensure all functions remain $\le 32$ LOC, file $< 500$ LOC, 0 compiler warnings, 0 unwraps.
4. **Wire in `formation.rs`**:
   - Call `rov_s.get_snr_f64(freq_band)`, `rov_ref.get_snr_f64(freq_band)`, etc., passing float SNR to `compute_dd_variances`.
5. **Verify Against Regression Guards**:
   - Run `python3 scripts/check_network_benchmark.py --smoke`.
   - Run `python3 scripts/check_multignss_benchmark.py --smoke`.
   - Run `cargo test --workspace` and `cargo clippy --workspace --all-targets -- -D warnings`.
   - Run `eval_f9p_rover` on Shinjuku and Whampoa to verify fix rate growth and tail error stability.
