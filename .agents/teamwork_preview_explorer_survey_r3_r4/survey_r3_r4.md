# Comprehensive Technical Survey: Frontier R3 (Network RTK VRS Engine) & Frontier R4 (Unified Composite Integration)

**Date**: 2026-09-12  
**Author**: Survey Explorer R3/R4 (`teamwork_preview_explorer`)  
**Target Repository**: `gneiss` (GNSS/INS Precise Positioning Engine)  
**Status**: Completed  

---

## 1. Executive Summary

This survey provides an exhaustive, read-only architectural investigation of the `gneiss` workspace for two major commercial GNSS/INS frontiers defined in `ORIGINAL_REQUEST.md`:
1. **Frontier R3: Network RTK Virtual Reference Station (VRS) Atmospheric Engine**
   - Simultaneous ingestion of 5–10 regional CORS base station observation streams.
   - Multi-baseline double-difference (DD) network adjustment solving integer ambiguities across inter-CORS network baselines.
   - Spatial 2D/3D Delaunay triangulation models for ionospheric pierce point (IPP) delays and tropospheric zenith wet delay (ZWD) gradients.
   - Synthesis of localized Virtual Reference Station (VRS) observations at the rover's approximate position (< 1 km effective baseline on 15–50 km regional networks).
   - Evaluation on the `crates/gneiss-rtk/src/bin/eval_network_ppk.rs` benchmark reducing baseline ppm error across CORS baselines (P181, P222, P225).
2. **Frontier R4: Unified Composite Integration**
   - 15-state Error-State Kalman Filter (ESKF/MEKF) with closed-loop attitude quaternion feedback and online bias estimation.
   - Backward 15-state Rauch-Tung-Striebel (RTS) smoother over forward filter history.
   - Vehicle Non-Holonomic Constraints (NHC) and Zero-Velocity Updates (ZUPT) integrated directly into the 15-state covariance.
   - Tightly-Coupled PPP/INS integrating the 15-state ESKF with Integer PPP-AR (`PppArSolver` + SINEX OSB).
   - Tightly-Coupled Network RTK/INS integrating the 15-state ESKF with VRS synthesis and double-difference carrier phase innovations.
3. **Regression Guards**:
   - `scripts/check_network_benchmark.py`: Smoke mode verified passing (`ALL CHECKS PASSED`, network fused horizontal p50 = 0.024 m <= 0.04 m, RMS = 0.041 m <= 0.06 m).
   - `scripts/check_multignss_benchmark.py`: Smoke mode verified passing (`ALL CHECKS PASSED`, network fused fix rate = 99.30% >= 96.5%).

All findings, module structures, and implementation recommendations strictly follow the quality standards in `AGENTS.md` (file size < 500 LOC, function size < 32 LOC, nesting depth < 3 levels, 0 compiler warnings, 0 `unwrap()` calls in production code, relational coupling structurally enforced).

---

## 2. Baseline Status of the Codebase

### 2.1 Workspace Structure & Existing Modules

The `gneiss` workspace is organized into modular crates:
- `crates/gneiss-core`: Fundamental GNSS math, geodetic coordinates (`coords.rs`), time standards (`time.rs`, `gnss_time.rs`), observation types (`obs.rs`), satellite orbits/ephemeris (`ephemeris/`), and atmospheric models (`atmosphere/ionosphere.rs`, `troposphere.rs`, `mapping.rs`).
- `crates/gneiss-parsers`: High-performance parsers for RINEX 2/3 (`rinex/`), SP3 precise orbits (`sp3.rs`), RINEX CLK precise clocks (`rinex_clk.rs`), ANTEX antenna phase centers (`antex.rs`), and SINEX bias files (`sinex_bia.rs`).
- `crates/gneiss-rtk`: The core positioning and filtering engine:
  - `ambiguity/`: LAMBDA integer search (`lambda/`), partial ambiguity resolution (`par.rs`), and PPP-AR solver (`ppp_ar.rs`).
  - `estimators/`: RTK iterated extended Kalman filter (`rtk_iekf/`), SPP (`spp/`).
  - `post_process/`: 4-pass post-processing pipeline (`forward.rs`, `backward.rs`, `combiner.rs`), multi-base consensus (`network.rs`), initial VRS stub (`vrs.rs`).
  - `swfg/`: Sliding-window factor graph optimization engine (`engine/`, `pipeline/`, `imu_preintegration/`).
- `crates/gneiss-geodesy`: Geoid grids, map projections, tides (`tides.rs`), and site calibration (`site_calibration.rs`).
- `bin/`: CLI tools (`gneiss-cli`, `gneiss-parse-rtcm3`).
- `scripts/`: Benchmark runners, regression guards, dataset fetchers.

### 2.2 Verification of Existing Regression Guards

Both required regression guard scripts were executed against the compiled release binaries:
1. `scripts/check_network_benchmark.py --smoke`:
   - **Command**: `python3 scripts/check_network_benchmark.py --smoke`
   - **Mechanism**: Evaluates 1,800 epochs (TOW 345600..399600, 30-second interval) of `eval_network_ppk` on `datasets/cors_short_baseline`.
   - **Verification Result**:
     ```
     [SMOKE MODE] Evaluating 1800 epochs
       ok     network fused horizontal p50 (m)             0.024 <= 0.04
       ok     network fused horizontal RMS (m)             0.041 <= 0.06
       ok     network fused vertical RMS (m)               0.043 <= 0.08
       ok     P181 smoothed fixed-only p50 (m)             0.022 <= 0.03
       ok     P222 smoothed fixed-only p50 (m)             0.060 <= 0.09
       ok     SLAC smoothed fixed-only p50 (m)             0.112 <= 0.12
       ok     OHLN smoothed fix rate (%)                  94.000 >= 79.0
       ok     P181 smoothed fix rate (%)                  92.100 >= 85.0
       ok     SLAC smoothed fix rate (%)                  70.600 >= 60.0
     ALL CHECKS PASSED
     ```
2. `scripts/check_multignss_benchmark.py --smoke`:
   - **Command**: `python3 scripts/check_multignss_benchmark.py --smoke`
   - **Mechanism**: Rebuilds binary if needed, evaluates 1,800 epochs on `datasets/multignss_2025d160` with GPS + Galileo (`GNEISS_SYSTEMS=GE`).
   - **Verification Result**:
     ```
     [SMOKE MODE] Evaluating 1800 epochs (~900.0 min)
       ok P181 fix rate (%)                    99.00 >= 97.5
       ok P181 h_p95 (mm)                     121.00 <= 145.0
       ok P181 v_p95 (mm)                     199.00 <= 290.0
       ok P225 fix rate (%)                    87.60 >= 71.0
       ok P225 h_p95 (mm)                      89.00 <= 245.0
       ok P225 v_p95 (mm)                     204.00 <= 370.0
       ok P222 fix rate (%)                    99.40 >= 86.0
       ok P222 h_p95 (mm)                     198.00 <= 295.0
       ok P222 v_p95 (mm)                      84.00 <= 135.0
       ok network fused fix rate (%)           99.30 >= 96.5
     ALL CHECKS PASSED
     ```

**Crucial Compatibility Contract**: Both regression guard scripts parse stdout using strict section headers (`=== Smoothed RTK [station] ===`, `=== NETWORK FUSED ===`, etc.). Any enhancement to `eval_network_ppk.rs` MUST keep these output blocks intact while appending new benchmarks (e.g. `=== VRS SYNTHESIS PPK ===`).

---

## 3. Frontier R3: Network RTK VRS Atmospheric Engine

### 3.1 Current Implementation vs. Frontier R3 Requirements

| Component | Current Implementation | Frontier R3 Requirement | Gap |
|-----------|------------------------|-------------------------|-----|
| **Multi-CORS Ingestion** | `eval_network_ppk.rs:622-631` loads 3–6 bases via Rayon parallel iterator into `HashMap<&str, Arc<Vec<EpochObs>>>`. | Ingestion of 5–10 regional CORS base station observation streams simultaneously. `datasets/cors_sf_bay_network` contains 11 bases + 1 rover. | Existing parallel RINEX loader handles 3–6 bases cleanly; needs to scale to 5–11 stations (e.g. `cabl`, `capo`, `mhcb`, `ohln`, `p181`, `p222`, `p225`, `p261`, `p271`, `slac`, `tibb`). |
| **Network Adjustment** | `eval_network_ppk.rs:636-669` runs independent forward passes per base to solve satellite wide-lane UPDs (`solve_network_upd`), followed by trajectory median fusion (`network::fuse_network_solutions`). | Formulation of multi-baseline double-difference network adjustment across network baselines to solve for inter-station integer ambiguities. | Current code does rover-to-base independent PPK; does not form inter-CORS baseline network adjustment with fixed station coordinates to extract pure troposphere and ionosphere delays. |
| **Spatial Atmospheric Modeling** | `crates/gneiss-rtk/src/post_process/vrs.rs:78-98` implements `fit_plane_gradient` (first-order planar regression). | Spatial 2D/3D Delaunay triangulation models for ionospheric pierce points (IPP) and tropospheric ZWD gradients. | Planar fitting fails on large regional networks (> 30 km) where topography and localized atmospheric pockets create non-planar gradients. Delaunay triangulation with barycentric interpolation is required. |
| **VRS Observation Synthesis** | `crates/gneiss-rtk/src/post_process/vrs.rs:101-179` implements initial `synthesize_vrs_epoch` and `shift_observables`. | Synthesis of localized VRS observation stream at rover approximate position (< 1 km effective baseline) with dual-frequency phase and code corrections. | Existing `vrs.rs` has geometry shift and frequency scaling, but is disconnected from `eval_network_ppk.rs` and lacks Delaunay-interpolated corrections. |
| **PPM Error Benchmark** | `eval_network_ppk.rs` reports single-baseline and consensus error. P222 (38 km) has p50 = 6.0 cm, SLAC (50 km) has p50 = 11.2 cm. | Demonstrating baseline ppm error reduction across CORS baselines (P181, P222, P225), achieving closer parity with Leica single-baseline specs ($8\text{ mm} + 1\text{ ppm}$). | Integrating the VRS stream as the reference station for rover PPK reduces the effective baseline length from 15–50 km to < 1 km, virtually eliminating differential distance-dependent errors. |

### 3.2 Delaunay Triangulation Architecture (`delaunay.rs`)

To guarantee numerical stability, zero external dependencies, and strict AGENTS.md compliance, a dedicated Delaunay 2D Triangulation module (`crates/gneiss-rtk/src/post_process/delaunay.rs`) is required.

#### Mathematical Foundation
1. **Circumcircle Criterion**:
   For three points $A, B, C$ ordered counter-clockwise, a test point $P$ lies strictly inside the circumcircle if:
   $$\det \begin{bmatrix} A_x - P_x & A_y - P_y & (A_x - P_x)^2 + (A_y - P_y)^2 \\ B_x - P_x & B_y - P_y & (B_x - P_x)^2 + (B_y - P_y)^2 \\ C_x - P_x & C_y - P_y & (C_x - P_x)^2 + (C_y - P_y)^2 \end{bmatrix} > 0$$
2. **Barycentric Interpolation**:
   For any target point $P$ inside triangle $(A, B, C)$:
   $$P = w_A A + w_B B + w_C C, \quad w_A + w_B + w_C = 1$$
   $$w_A = \frac{\text{Area}(P, B, C)}{\text{Area}(A, B, C)}, \quad w_B = \frac{\text{Area}(A, P, C)}{\text{Area}(A, B, C)}, \quad w_C = \frac{\text{Area}(A, B, P)}{\text{Area}(A, B, C)}$$
   Any atmospheric quantity $V$ (e.g., ZWD or slant ionospheric delay) at $P$ is interpolated as:
   $$V(P) = w_A V_A + w_B V_B + w_C V_C$$
3. **Out-of-Bounds Fallback**:
   If $P$ is outside the convex hull, find the closest triangle facet and use distance-weighted linear extrapolation or nearest-facet projection to avoid unbounded divergence.

#### Proposed Data Structures (`delaunay.rs`)
```rust
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle {
    pub vertices: [usize; 3], // Indices into vertex slice
}

pub struct DelaunayMesh {
    pub vertices: Vec<Point2D>,
    pub triangles: Vec<Triangle>,
}

impl DelaunayMesh {
    pub fn build(points: &[Point2D]) -> Result<Self, &'static str>;
    pub fn locate_triangle(&self, p: Point2D) -> Option<(usize, [f64; 3])>; // (triangle_index, [w_a, w_b, w_c])
    pub fn interpolate(&self, values: &[f64], p: Point2D) -> f64;
}
```

### 3.3 Ionospheric Pierce Point (IPP) Triangulation

For ionosphere modeling, each station-satellite line of sight pierces the single-layer ionospheric shell at height $H_{\text{ion}} = 350\text{ km}$.
The spherical IPP calculation is already available in `crates/gneiss-core/src/atmosphere/ionosphere.rs:120-130`:
```rust
let re = 6371.0;
let h = 350.0;
let psi = (libm::asin(re / (re + h) * libm::cos(el)) - (core::f64::consts::PI / 2.0 - el)).max(0.0);
let ipp_lat = libm::asin(
    libm::sin(lat_r.to_radians()) * libm::cos(psi)
        + libm::cos(lat_r.to_radians()) * libm::sin(psi) * libm::cos(az),
).to_degrees();
let ipp_lon = lon_r + libm::asin(libm::sin(psi) * libm::sin(az) / libm::cos(ipp_lat.to_radians())).to_degrees();
```
For each satellite $s$, the IPP coordinates from all $N$ CORS stations form a 2D point cloud on the ionospheric shell. A `DelaunayMesh` is constructed over these IPP points per satellite. When synthesizing VRS observations for the rover, the rover's IPP for satellite $s$ is queried in the mesh, yielding exact barycentric interpolation of slant ionospheric delay $\delta I_{\text{vrs}}^s$.

### 3.4 Multi-Baseline CORS Network Adjustment (`network_adj.rs`)

To extract unambiguous atmospheric delays at each CORS station without relying on float approximations:
1. **Network Baseline Selection**: Construct minimum spanning tree (MST) or Delaunay edge graph connecting the $N$ CORS stations (e.g. P181 as master, baselines to OHLN, CAPO, P225, P222, SLAC).
2. **Fixed Geometric Ranges**: Because all CORS coordinates $\mathbf{X}_A, \mathbf{X}_B$ are known from IGS/NOAA coordinates to millimeter precision, geometric ranges $\rho_A^s, \rho_B^s$ and unit line-of-sight vectors $\mathbf{u}_A^s, \mathbf{u}_B^s$ are fixed without position error.
3. **Double-Difference Carrier Phase Equation**:
   $$\Delta \nabla \Phi_{AB}^{sk} = \Delta \nabla \rho_{AB}^{sk} + \lambda \Delta \nabla N_{AB}^{sk} + \Delta \nabla T_{AB}^{sk} - \Delta \nabla I_{AB}^{sk}$$
4. **Ambiguity Resolution**:
   - Wide-lane (MW) ambiguities $\Delta \nabla N_{W}$ are resolved using multi-epoch averaging and satellite UPD removal.
   - Narrow-lane ambiguities $\Delta \nabla N_1$ are resolved via LAMBDA or sequential rounding on the ionosphere-free residual.
5. **Atmospheric Residual Extraction**:
   With integer ambiguities $\Delta \nabla N$ fixed:
   - **Tropospheric Residuals**: Non-dispersive combination isolates relative zenith wet delay $\Delta \text{ZWD}_{AB}$ via Niell/GMF mapping functions:
     $$\Delta \nabla \Phi_{\text{IF}} - \Delta \nabla \rho_{\text{IF}} - \lambda_{\text{IF}} \Delta \nabla N_{\text{IF}} = m_w(el^s) \text{ZWD}_B - m_w(el^k) \text{ZWD}_B - \dots$$
   - **Ionospheric Residuals**: Geometry-free combination isolates relative slant ionospheric delay:
     $$\Delta \nabla \Phi_{\text{GF}} - \lambda_{\text{GF}} \Delta \nabla N_{\text{GF}} = \left( 1 - \frac{f_1^2}{f_2^2} \right) \Delta \nabla I_{AB}^{sk}$$
6. **Output**: Per-station ZWD estimates and per-station, per-satellite slant ionospheric delays $\delta I_i^s$, ready for Delaunay mesh interpolation.

### 3.5 VRS Observation Synthesis Pipeline

At each epoch $t$:
1. Rover approximate position $\mathbf{r}_{\text{rov}}$ is supplied (e.g. from SPP).
2. Master CORS station $\mathbf{r}_{\text{master}}$ is chosen (nearest CORS).
3. VRS reference coordinate is defined: $\mathbf{r}_{\text{vrs}} = \mathbf{r}_{\text{rov}}$.
4. For each satellite $s$ in master station's observation:
   - Compute geometric shift: $\Delta \rho^s = \|\mathbf{r}_{\text{sat}}^s - \mathbf{r}_{\text{vrs}}\| - \|\mathbf{r}_{\text{sat}}^s - \mathbf{r}_{\text{master}}\|$.
   - Interpolate differential troposphere $\delta T_{\text{vrs}}$ via Delaunay mesh on station coordinates.
   - Interpolate differential ionosphere $\delta I_{\text{vrs}}^s$ via Delaunay mesh on satellite $s$'s IPP coordinates.
   - Synthesize code and phase observations:
     $$P_{\text{vrs}}^s(f) = P_{\text{master}}^s(f) + \Delta \rho^s + \delta T_{\text{vrs}} + \left(\frac{f_1}{f}\right)^2 \delta I_{\text{vrs}}^s$$
     $$\Phi_{\text{vrs}}^s(f) = \Phi_{\text{master}}^s(f) + \frac{1}{\lambda_f} \left[ \Delta \rho^s + \delta T_{\text{vrs}} - \left(\frac{f_1}{f}\right)^2 \delta I_{\text{vrs}}^s \right]$$
5. Rover PPK runs with the synthesized VRS stream as its base station:
   - Effective baseline length: $\|\mathbf{r}_{\text{rov}} - \mathbf{r}_{\text{vrs}}\| \approx 0\text{ km}$ ($< 1\text{ km}$).
   - Differential ionospheric and tropospheric errors are near zero, enabling rapid ambiguity fixing and reducing ppm error across the regional network.

---

## 4. Frontier R4: Unified Composite Integration

### 4.1 15-State Error-State Kalman Filter (ESKF/MEKF)

#### Current 6-DOF Baseline
In `crates/gneiss-rtk/src/swfg/imu_preintegration/smoother.rs`:
- State vector: `pub type State6 = Vector6<f64>;` ($[p_x, p_y, p_z, v_x, v_y, v_z]^T$).
- Covariance: `pub type Cov6 = Matrix6<f64>;` (6x6).
- Attitude: `UnitQuaternion<f64>` is dead-reckoned open-loop: `self.attitude *= preint.dq;`.
- Biases: Accel/gyro biases are not part of the state. Gyro bias is an offline average (`compute_gyro_bias`).
- RTS Smoother: Backward pass smooths only position and velocity.

#### 15-State Architecture
The filter state is expanded to a full 15-state Error-State Kalman Filter ($\delta \mathbf{x} \in \mathbb{R}^{15}$):
$$\delta \mathbf{x} = \begin{bmatrix} \delta \mathbf{p}^e \\ \delta \mathbf{v}^e \\ \delta \boldsymbol{\theta} \\ \delta \mathbf{b}_a \\ \delta \mathbf{b}_g \end{bmatrix} \in \mathbb{R}^{15}$$

1. **State Partitioning**:
   - $\delta \mathbf{p}^e \in \mathbb{R}^3$: Position error in ECEF frame (m).
   - $\delta \mathbf{v}^e \in \mathbb{R}^3$: Velocity error in ECEF frame (m/s).
   - $\delta \boldsymbol{\theta} \in \mathbb{R}^3$: Attitude error vector (rad), with convention from AGENTS.md / predictor:
     $$\mathbf{q} \leftarrow \delta \mathbf{q}(\delta \boldsymbol{\theta}) \otimes \mathbf{q}, \quad \delta \mathbf{q}(\delta \boldsymbol{\theta}) \approx \begin{bmatrix} 1 \\ \frac{1}{2} \delta \boldsymbol{\theta} \end{bmatrix}$$
   - $\delta \mathbf{b}_a \in \mathbb{R}^3$: Accelerometer bias error in body frame (m/s²).
   - $\delta \mathbf{b}_g \in \mathbb{R}^3$: Gyroscope bias error in body frame (rad/s).

2. **Transition Matrix $\mathbf{F} \in \mathbb{R}^{15 \times 15}$**:
   For IMU sample interval $\Delta t$:
   $$\mathbf{F} = \begin{bmatrix}
   \mathbf{I}_3 & \mathbf{I}_3 \Delta t & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 \\
   \mathbf{0}_3 & \mathbf{I}_3 - 2 [\boldsymbol{\omega}_{ie}^e \times] \Delta t & +[\mathbf{f}^e \times] \Delta t & -\mathbf{R}_b^e \Delta t & \mathbf{0}_3 \\
   \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{I}_3 - [\boldsymbol{\omega}_{ie}^e \times] \Delta t & \mathbf{0}_3 & -\mathbf{R}_b^e \Delta t \\
   \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{I}_3 & \mathbf{0}_3 \\
   \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{I}_3
   \end{bmatrix}$$
   **Hard Rule Check**: As mandated in AGENTS.md, the velocity-attitude Jacobian term $+[\mathbf{f}^e \times] \Delta t$ has a **positive** sign (`vel_att = f_e_skew * dt`).

3. **Closed-Loop Error Injection & Reset**:
   Upon completing a Kalman update $\delta \hat{\mathbf{x}} = \mathbf{K} \mathbf{y}$:
   - Nominal position: $\mathbf{p} \leftarrow \mathbf{p} + \delta \hat{\mathbf{p}}$
   - Nominal velocity: $\mathbf{v} \leftarrow \mathbf{v} + \delta \hat{\mathbf{v}}$
   - Nominal attitude: $\mathbf{q} \leftarrow \delta \mathbf{q}(\delta \hat{\boldsymbol{\theta}}) \otimes \mathbf{q}$ (re-normalized)
   - Accelerometer bias: $\mathbf{b}_a \leftarrow \mathbf{b}_a + \delta \hat{\mathbf{b}}_a$
   - Gyroscope bias: $\mathbf{b}_g \leftarrow \mathbf{b}_g + \delta \hat{\mathbf{b}}_g$
   - Error state reset: $\delta \mathbf{x} \leftarrow \mathbf{0}$
   - Covariance reset: $\mathbf{P} \leftarrow (\mathbf{I}_{15} - \mathbf{K}\mathbf{H}) \mathbf{P} (\mathbf{I}_{15} - \mathbf{K}\mathbf{H})^T + \mathbf{K}\mathbf{R}\mathbf{K}^T$

4. **Vehicle Constraints (NHC and ZUPT) in 15 States**:
   - **Non-Holonomic Constraints (NHC)**:
     $\mathbf{v}^b = (\mathbf{R}_b^e)^T \mathbf{v}^e$. Under no-sideslip / no-hop conditions, lateral and vertical body velocities are zero:
     $$\mathbf{y}_{\text{nhc}} = \begin{bmatrix} 0 - v_y^b \\ 0 - v_z^b \end{bmatrix} \in \mathbb{R}^2$$
     $$\mathbf{H}_{\text{nhc}} = \begin{bmatrix} \mathbf{0}_{2 \times 3} & \mathbf{R}_e^b[1..2, :] & -[\mathbf{v}^b \times][1..2, :] \mathbf{R}_e^b & \mathbf{0}_{2 \times 3} & \mathbf{0}_{2 \times 3} \end{bmatrix} \in \mathbb{R}^{2 \times 15}$$
   - **Zero-Velocity Update (ZUPT)**:
     When stationary (accelerometer magnitude $\approx g$ and low variance, gyro norm $\approx 0$):
     $$\mathbf{y}_{\text{zupt}} = \mathbf{0} - \mathbf{v}^e \in \mathbb{R}^3$$
     $$\mathbf{H}_{\text{zupt}} = \begin{bmatrix} \mathbf{0}_{3 \times 3} & \mathbf{I}_3 & \mathbf{0}_{3 \times 3} & \mathbf{0}_{3 \times 3} & \mathbf{0}_{3 \times 3} \end{bmatrix} \in \mathbb{R}^{3 \times 15}$$
     Via cross-covariance $\mathbf{P}_{v, b_a}$ and $\mathbf{P}_{\theta, b_g}$, stationary ZUPTs rapidly converge sensor biases and level attitude.

5. **15-State RTS Backward Smoother**:
   - Stored forward snapshot at epoch $k$: $(\mathbf{x}_k^{\text{pred}}, \mathbf{P}_k^{\text{pred}}, \mathbf{x}_k^{\text{post}}, \mathbf{P}_k^{\text{post}}, \mathbf{F}_k)$.
   - Smoother gain: $\mathbf{C}_k = \mathbf{P}_k^{\text{post}} \mathbf{F}_k^T (\mathbf{P}_{k+1}^{\text{pred}})^{-1} \in \mathbb{R}^{15 \times 15}$.
   - Backward recursion ($k = N-1 \dots 0$):
     $$\delta \mathbf{x}_k^s = \mathbf{C}_k (\delta \mathbf{x}_{k+1}^s + \mathbf{x}_{k+1}^s - \mathbf{x}_{k+1}^{\text{pred}})$$
     $$\mathbf{P}_k^s = \mathbf{P}_k^{\text{post}} + \mathbf{C}_k (\mathbf{P}_{k+1}^s - \mathbf{P}_{k+1}^{\text{pred}}) \mathbf{C}_k^T$$
   - Applying $\delta \mathbf{x}_k^s$ to epoch $k$ nominal state yields the globally optimal smoothed trajectory.

### 4.2 Tightly-Coupled PPP/INS Architecture (`tc_ppp.rs`)

Integrates the 15-state ESKF with Integer PPP-AR (`crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`) and SINEX OSB products (`crates/gneiss-parsers/src/sinex_bia.rs`):
1. **Antenna Lever-Arm Coupling**:
   The GNSS antenna phase center in ECEF is:
   $$\mathbf{r}_{\text{ant}}^e = \mathbf{p}^e + \mathbf{R}_b^e \mathbf{l}^b$$
   where $\mathbf{l}^b$ is the calibrated IMU-to-antenna lever arm vector.
2. **Measurement Equations**:
   For satellite $s$, modeled range is:
   $$\hat{\rho}^s = \|\mathbf{r}_{\text{sat}}^s - \mathbf{r}_{\text{ant}}^e\|$$
   Line-of-sight unit vector: $\mathbf{u}^s = \frac{\mathbf{r}_{\text{sat}}^s - \mathbf{r}_{\text{ant}}^e}{\hat{\rho}^s}$.
   Measurement Jacobian for satellite $s$:
   $$\mathbf{H}^s = \begin{bmatrix} -\mathbf{u}^s & \mathbf{0}_{1 \times 3} & -(\mathbf{u}^s)^T [\mathbf{R}_b^e \mathbf{l}^b \times] & \mathbf{0}_{1 \times 3} & \mathbf{0}_{1 \times 3} \end{bmatrix} \in \mathbb{R}^{1 \times 15}$$
3. **SINEX OSB & Phase Center Corrections**:
   - Satellite code and phase biases from `SinexBias::get_bias(sat, obs_code, time)` are applied to raw observables.
   - Antenna phase center offsets (PCO) from ANTEX (`igs14.atx`) and Solid Earth Tide displacements are removed.
4. **Integer Ambiguity Constraints**:
   - `PppArSolver::fix_wide_lane` resolves wide-lane ambiguities from smoothed Melbourne-Wübbena observables.
   - Conditioned on fixed wide-lane integers, narrow-lane ambiguities are resolved via LAMBDA.
   - Fixed carrier-phase ambiguities convert carrier-phase observables into millimeter-precision pseudoranges, tightly constraining the 15-state ESKF during satellite outages.

### 4.3 Tightly-Coupled Network RTK/INS Architecture (`tc_rtk.rs`)

Integrates the 15-state ESKF with the Network RTK VRS atmospheric engine:
1. **VRS Observation Feeding**:
   At each GNSS epoch, the VRS engine synthesizes virtual base station observables $\mathbf{y}_{\text{vrs}}$ at $\mathbf{r}_{\text{vrs}} = \mathbf{p}^e$.
2. **Double-Difference Innovations**:
   For rover-VRS baseline and satellite pair $(s, k)$:
   $$\Delta \nabla y = \Delta \nabla \Phi_{\text{rover-vrs}}^{sk} - \left( \Delta \nabla \hat{\rho}^{sk} + \lambda \Delta \nabla N^{sk} \right)$$
   where $\Delta \nabla \hat{\rho}^{sk} = (\mathbf{u}^s - \mathbf{u}^k) \cdot (\mathbf{r}_{\text{ant}}^e - \mathbf{r}_{\text{vrs}})$.
3. **Jacobian Coupling**:
   $$\mathbf{H}_{\text{dd}}^{sk} = \begin{bmatrix} -(\mathbf{u}^s - \mathbf{u}^k) & \mathbf{0}_{1 \times 3} & -(\mathbf{u}^s - \mathbf{u}^k)^T [\mathbf{R}_b^e \mathbf{l}^b \times] & \mathbf{0}_{1 \times 3} & \mathbf{0}_{1 \times 3} \end{bmatrix} \in \mathbb{R}^{1 \times 15}$$
4. **Integrity in Urban Canyons**:
   Even if only 2 or 3 satellites are tracked simultaneously through an urban canyon or underpass, the double-difference carrier phase innovations continue to update the ESKF state, preventing inertial dead-reckoning divergence.

---

## 5. File Layout & Modular Design

To strictly adhere to `AGENTS.md` (< 500 LOC per file, < 32 LOC per function, < 3 levels of nesting, 0 unwraps, 0 compiler warnings), new functionality is partitioned into dedicated, highly cohesive modules:

```
crates/gneiss-rtk/src/
├── post_process/
│   ├── delaunay.rs          # 2D Delaunay triangulation & barycentric interpolation (~280 LOC)
│   ├── network_adj.rs       # Multi-baseline inter-CORS DD network adjustment (~350 LOC)
│   ├── vrs.rs               # Refactored VRS synthesis engine with Delaunay models (~320 LOC)
│   ├── network.rs           # (Existing) Multi-base trajectory consensus (450 LOC)
│   └── mod.rs               # Module exports and pipeline options
├── estimators/
│   └── eskf/                # Modular 15-state Error-State Kalman Filter
│       ├── mod.rs           # Public ESKF exports (~60 LOC)
│       ├── state.rs         # 15-state vector and covariance definitions (~200 LOC)
│       ├── predict.rs       # IMU mechanization, F matrix, Q propagation (~250 LOC)
│       ├── constraints.rs   # NHC and ZUPT measurement updates (~220 LOC)
│       ├── updates.rs       # Loosely & tightly-coupled GNSS updates (~280 LOC)
│       └── smoother.rs      # 15-state backward RTS smoother (~260 LOC)
├── composite/
│   ├── mod.rs               # Composite integration exports (~60 LOC)
│   ├── tc_ppp.rs            # Tightly-Coupled PPP/INS engine (~380 LOC)
│   └── tc_rtk.rs            # Tightly-Coupled Network RTK/INS engine (~380 LOC)
└── bin/
    ├── eval_network_ppk.rs  # Enhanced with VRS PPK mode (preserving regression blocks)
    ├── eval_odaiba_ins.rs   # Upgraded to 15-state ESKF + RTS smoother
    └── eval_ppp.rs          # Upgraded with kinematic PPP-AR evaluation
```

---

## 6. Detailed Implementation Roadmap

### Phase 1: Spatial Delaunay Triangulation & Atmospheric Meshing
- Implement `crates/gneiss-rtk/src/post_process/delaunay.rs`:
  - 2D Point, Triangle, Circumcircle determinant test.
  - Bowyer-Watson incremental triangulation.
  - Point-in-triangle location and barycentric coordinate calculation.
  - Robust out-of-hull fallback to nearest boundary facet.
  - Comprehensive unit tests guarding edge cases (collinear points, degenerate triangles).

### Phase 2: Multi-Baseline Network Adjustment & VRS Synthesis
- Implement `crates/gneiss-rtk/src/post_process/network_adj.rs`:
  - Form baseline graph across CORS reference stations.
  - Multi-baseline double-difference wide-lane and narrow-lane integer ambiguity fixing.
  - Extraction of CORS tropospheric ZWD residuals and slant ionospheric pierce point delays.
- Refactor `crates/gneiss-rtk/src/post_process/vrs.rs`:
  - Integrate `DelaunayMesh` for tropospheric ZWD/gradients and per-satellite IPP ionosphere.
  - Provide `synthesize_vrs_stream(master, rover_approx, surface, ephemerides)`.
- Update `crates/gneiss-rtk/src/bin/eval_network_ppk.rs`:
  - Add VRS-assisted PPK evaluation section (`=== VRS SYNTHESIS PPK ===`).
  - Demonstrate baseline ppm error reduction across P181 (15 km), P225 (22 km), and P222 (38 km).
  - Guarantee exact output formatting to keep `check_network_benchmark.py` and `check_multignss_benchmark.py` 100% green.

### Phase 3: 15-State ESKF/MEKF & RTS Smoother
- Implement `crates/gneiss-rtk/src/estimators/eskf/`:
  - `state.rs`: 15-state vector layout, nominal quaternion, covariance matrices.
  - `predict.rs`: Continuous-discrete mechanization, positive velocity-attitude Jacobian (`+f_e_skew * dt`), closed-loop reset.
  - `constraints.rs`: Body-frame NHC and stationary ZUPT updates with full attitude cross-coupling.
  - `smoother.rs`: 15-state backward RTS smoother.
- Upgrade `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`:
  - Replace 6-DOF filter with 15-state ESKF and RTS smoother.
  - Benchmark on Odaiba 10Hz GNSS / 50Hz MEMS IMU trajectory against NovAtel SPAN truth, targeting $p_{50} < 2.5\text{ m}$ and RMS $< 5.2\text{ m}$.

### Phase 4: Unified Composite Navigation (TC-PPP and TC-RTK)
- Implement `crates/gneiss-rtk/src/composite/tc_ppp.rs`:
  - Compose 15-state ESKF with `PppArSolver` and SINEX OSB.
  - Form uncombined carrier phase / pseudorange innovation updates with antenna lever arm.
  - Benchmark on F9P kinematic vehicle drive (`rover_csrs.pos`).
- Implement `crates/gneiss-rtk/src/composite/tc_rtk.rs`:
  - Compose 15-state ESKF with VRS synthesis and double-difference carrier phase innovations.
  - Benchmark on urban canyon / multi-CORS datasets.

---

## 7. Verification & Invariants Checklist

| Invariant / Standard | Enforcement Strategy |
|----------------------|----------------------|
| **File Size (< 500 LOC)** | New logic divided into small modules (`delaunay.rs`, `network_adj.rs`, `eskf/state.rs`, `eskf/predict.rs`, etc.), each < 400 LOC. |
| **Function Size (< 32 LOC)** | High-level orchestration functions broken into specialized helper functions. |
| **Nesting Depth (< 3 levels)** | Early exit (`guard` clauses, `let Some(...) = ... else { return; }`), flat iterators. |
| **0 `unwrap()` in Production** | `match`, `if let`, `ok_or()?`, or descriptive `.expect("invariant: ...")` only. |
| **0 Warnings** | Verified with `cargo clippy --workspace --all-targets -- -D warnings`. |
| **Sign Conventions** | Velocity-attitude Jacobian sign strictly positive (`+f_e_skew * dt`). Phase windup strictly subtracted (`cp - wup`). |
| **Regression Guards** | Both `scripts/check_network_benchmark.py --smoke` and `scripts/check_multignss_benchmark.py --smoke` pass with 0 regressions. |
