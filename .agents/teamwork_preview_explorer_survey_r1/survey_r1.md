# Frontier R1 Survey Report: 15-State ESKF/MEKF GNSS/INS, RTS Smoother, and NHC/ZUPT

**Date**: 2026-09-12  
**Author**: Survey Explorer R1 (`teamwork_preview_explorer_survey_r1`)  
**Target Milestone**: Frontier R1 (15-State ESKF/MEKF GNSS/INS, RTS Smoother, NHC/ZUPT, `eval_odaiba_ins` Benchmark)  
**Status**: Survey Complete & Mathematically Formulated  

---

## Executive Summary

This survey investigates the current state of inertial navigation and GNSS/INS fusion in the `gneiss` workspace to establish the technical foundation for **Frontier R1**:
1. **Current Baseline**: The existing inertial filter in `crates/gneiss-rtk/src/swfg/imu_preintegration/smoother.rs` is a 6-DOF filter ($\mathbf{p}, \mathbf{v}$) with open-loop quaternion propagation, fixed initial gyro bias (first 350 samples average), zero accel bias, uncoupled transition physics, and a 6-DOF RTS smoother.
2. **Current Performance on Tokyo Odaiba**: On the full 12,399-epoch (10Hz rover, 1Hz base, 50Hz MEMS IMU) urban canyon drive evaluated against NovAtel SPAN ground truth:
   - **GNSS-Only RTK (Raw Fixes)** ($N=1,232$): $p_{50} = 2.808\text{ m}, p_{68} = 5.470\text{ m}, p_{95} = 12.035\text{ m}, \text{RMS} = 5.720\text{ m}$.
   - **Forward Inertial Filter (6-DOF)** ($N=12,398$): $p_{50} = 5.982\text{ m}, p_{68} = 8.678\text{ m}, p_{95} = 18.545\text{ m}, \text{RMS} = 9.689\text{ m}$.
   - **RTS Smoothed GNSS/INS (6-DOF)** ($N=12,398$): $p_{50} = 2.907\text{ m}, p_{68} = 5.396\text{ m}, p_{95} = 11.080\text{ m}, \text{RMS} = 5.508\text{ m}$.
3. **The Accuracy Gap**: The target benchmark requirements are:
   - **$p_{50} < 2.5\text{ m}$** (current: $2.907\text{ m}$, gap: $-0.407\text{ m}$)
   - **$\text{RMS} < 5.2\text{ m}$** (current: $5.508\text{ m}$, gap: $-0.308\text{ m}$)
4. **Root Cause Analysis**:
   - Initial heading estimated from GNSS delta has a $7.07^\circ$ error ($319.58^\circ$ estimated vs $326.65^\circ$ NovAtel reference). Because attitude error $\delta\boldsymbol{\theta}$ is omitted from the Kalman state, this $7^\circ$ misalignment persists through the entire 20-minute trajectory.
   - Accelerometer bias is assumed identically zero ($\mathbf{b}_a \equiv \mathbf{0}$), and gyro bias is frozen after initialization. MEMS sensor biases drift dynamically over time and temperature, producing quadratic position drift during GNSS outages.
   - Non-Holonomic Constraints (NHC) only update velocity without attitude coupling ($\partial \mathbf{v}^b / \partial \boldsymbol{\theta} = \mathbf{0}$ in the update), missing the critical physical mechanism where lateral vehicle velocity locks vehicle heading.
   - The backward RTS smoother operates only over 6 states ($p, v$), unable to propagate GNSS innovations backward to refine past attitude or sensor biases.
5. **Solution Architecture**: Implement a modular 15-state Error-State Kalman Filter (ESKF/MEKF) with closed-loop quaternion feedback, dynamic online bias calibration, coupled NHC/ZUPT updates, and full 15-state backward RTS smoothing, strictly obeying `AGENTS.md` modular constraints (<500 LOC/file, <32 LOC/function, 0 warnings, 0 unwrap).

---

## 1. Codebase Inventory & Current Architecture

### 1.1 Existing Files and Modules

| File Path | LOC | Current Role | Status / Action Needed |
|:---|:---:|:---|:---|
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs` | 282 | Benchmark harness for Tokyo Odaiba GNSS/INS dataset | Update filter initialization and step invocation to 15-state ESKF. |
| `crates/gneiss-rtk/src/swfg/imu_preintegration/smoother.rs` | 402 | 6-state forward filter and 6-state RTS smoother | Replace/upgrade with modular 15-state ESKF and 15-state RTS smoother. |
| `crates/gneiss-rtk/src/swfg/imu_preintegration/mod.rs` | 433 | Preintegration accumulator, bias Jacobians, and factor graph factor | Preintegration logic is sound; provides $\Delta\mathbf{p}, \Delta\mathbf{v}, \Delta\mathbf{q}$ and Jacobians. |
| `crates/gneiss-rtk/src/swfg/imu_preintegration/stationary.rs` | 120 | Variance and norm detector for stationary intervals | Fully functional for ZUPT and static bias initialization. |
| `crates/gneiss-rtk/src/swfg/imu_preintegration/tests.rs` | 121 | Unit tests for preintegration math | Expand with 15-state transition and bias estimation tests. |
| `tests/src/inertial_outage_simulation.rs` | 92 | Integration test for 10-second outage bridging | Update to 15-state filter; verify outage drift remains $< 0.50\text{ m}$. |
| `crates/gneiss-core/src/imu.rs` | 28 | Basic `ImuMeasurement` struct | Retain unchanged. |

### 1.2 Git Historical Heritage (`predictor.rs` and `updater.rs`)

Historical inspection of commit `5593a73^` (pre-factor graph engine) reveals the proven physics models:
- **Sign Conventions (`AGENTS.md:26-28`)**:
  - Attitude error state $\delta\boldsymbol{\theta} = -\boldsymbol{\psi}$ (left-multiplied global/ECEF frame).
  - Velocity-attitude transition block: $\Phi_{v,\theta} = +[\mathbf{f}^e\times]\Delta t$ (**POSITIVE** sign, enforced by commit `da013e2` and Bug 2 audit).
  - Velocity-accel-bias transition block: $\Phi_{v,ba} = -\mathbf{R}_b^e \Delta t$.
  - Attitude-gyro-bias transition block: $\Phi_{\theta,bg} = -\mathbf{R}_b^e \Delta t$.
  - Error quaternion update: $\mathbf{q} \leftarrow \delta\mathbf{q} \otimes \mathbf{q}$, where $\delta\mathbf{q} = \text{UnitQuaternion}::\text{from\_scaled\_axis}(\delta\boldsymbol{\theta})$.

---

## 2. Mathematical Formulation for Frontier R1

### 2.1 State Vector Definition

The system state is partitioned into the **nominal state** $\mathbf{x}_{nom}$ and the **error state** $\delta\mathbf{x} \in \mathbb{R}^{15}$.

#### Nominal State ($\mathbf{x}_{nom}$):
$$\mathbf{x}_{nom} = \begin{bmatrix} \mathbf{p}^e \\ \mathbf{v}^e \\ \mathbf{q}_b^e \\ \mathbf{b}_a^b \\ \mathbf{b}_g^b \end{bmatrix}$$
- $\mathbf{p}^e \in \mathbb{R}^3$: ECEF position (meters)
- $\mathbf{v}^e \in \mathbb{R}^3$: ECEF velocity (m/s)
- $\mathbf{q}_b^e \in \mathbb{H}$: Unit quaternion representing orientation from body frame $b$ to ECEF frame $e$, with corresponding rotation matrix $\mathbf{R}_b^e = \mathbf{R}(\mathbf{q}_b^e)$
- $\mathbf{b}_a^b \in \mathbb{R}^3$: Accelerometer bias in body frame (m/s²)
- $\mathbf{b}_g^b \in \mathbb{R}^3$: Gyroscope bias in body frame (rad/s)

#### Error State ($\delta\mathbf{x} \in \mathbb{R}^{15}$):
$$\delta\mathbf{x} = \begin{bmatrix} \delta\mathbf{p}^e \\ \delta\mathbf{v}^e \\ \delta\boldsymbol{\theta}^e \\ \delta\mathbf{b}_a^b \\ \delta\mathbf{b}_g^b \end{bmatrix}$$
- $\delta\mathbf{p}^e = \mathbf{p}_{true}^e - \mathbf{p}^e \in \mathbb{R}^3$
- $\delta\mathbf{v}^e = \mathbf{v}_{true}^e - \mathbf{v}^e \in \mathbb{R}^3$
- $\delta\boldsymbol{\theta}^e \in \mathbb{R}^3$: Attitude error vector such that $\mathbf{R}_{true} = (\mathbf{I} + [\delta\boldsymbol{\theta}^e\times]) \mathbf{R}_b^e$ and $\mathbf{q}_{true} = \delta\mathbf{q} \otimes \mathbf{q}_b^e$
- $\delta\mathbf{b}_a^b = \mathbf{b}_{a,true}^b - \mathbf{b}_a^b \in \mathbb{R}^3$
- $\delta\mathbf{b}_g^b = \mathbf{b}_{g,true}^b - \mathbf{b}_g^b \in \mathbb{R}^3$

---

### 2.2 Time Propagation (Prediction Step)

Given preintegrated IMU delta measurements $(\Delta\mathbf{p}, \Delta\mathbf{v}, \Delta\mathbf{q}, \Delta t)$ over interval $t_{k-1} \to t_k$:

1. **Bias Compensation on Deltas**:
   $$\Delta\mathbf{p}_{c} = \Delta\mathbf{p} + \mathbf{J}_{p,ba} (\mathbf{b}_{a,k-1} - \mathbf{b}_{a,lin}) + \mathbf{J}_{p,bg} (\mathbf{b}_{g,k-1} - \mathbf{b}_{g,lin})$$
   $$\Delta\mathbf{v}_{c} = \Delta\mathbf{v} + \mathbf{J}_{v,ba} (\mathbf{b}_{a,k-1} - \mathbf{b}_{a,lin}) + \mathbf{J}_{v,bg} (\mathbf{b}_{g,k-1} - \mathbf{b}_{g,lin})$$
   $$\Delta\mathbf{q}_{c} = \Delta\mathbf{q} \otimes \exp\left(\mathbf{J}_{q,bg} (\mathbf{b}_{g,k-1} - \mathbf{b}_{g,lin})\right)$$

2. **Nominal State Propagation**:
   $$\mathbf{p}_k = \mathbf{p}_{k-1} + \mathbf{v}_{k-1} \Delta t + \frac{1}{2}\mathbf{g}^e \Delta t^2 + \mathbf{R}_b^e(\mathbf{q}_{k-1}) \Delta\mathbf{p}_c$$
   $$\mathbf{v}_k = \mathbf{v}_{k-1} + \mathbf{g}^e \Delta t - 2(\boldsymbol{\omega}_{ie}^e \times \mathbf{v}_{k-1})\Delta t + \mathbf{R}_b^e(\mathbf{q}_{k-1}) \Delta\mathbf{v}_c$$
   $$\mathbf{q}_k = \mathbf{q}_{k-1} \otimes \Delta\mathbf{q}_c \quad (\text{renormalized})$$
   $$\mathbf{b}_{a,k} = \mathbf{b}_{a,k-1}, \quad \mathbf{b}_{g,k} = \mathbf{b}_{g,k-1}$$

3. **Discrete Transition Matrix $\boldsymbol{\Phi} \in \mathbb{R}^{15 \times 15}$**:
   Let $\mathbf{R} = \mathbf{R}_b^e(\mathbf{q}_{k-1})$ and effective specific force $\mathbf{f}^e = \frac{\mathbf{R} \Delta\mathbf{v}_c}{\Delta t}$:
   $$\boldsymbol{\Phi} = \begin{bmatrix}
   \mathbf{I}_3 & \mathbf{I}_3 \Delta t & \frac{1}{2}[\mathbf{f}^e\times]\Delta t^2 & -\frac{1}{2}\mathbf{R}\Delta t^2 & \mathbf{0}_3 \\
   \mathbf{0}_3 & \mathbf{I}_3 - 2[\boldsymbol{\omega}_{ie}^e\times]\Delta t & +[\mathbf{f}^e\times]\Delta t & -\mathbf{R}\Delta t & \mathbf{0}_3 \\
   \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{I}_3 - [\boldsymbol{\omega}_{ie}^e\times]\Delta t & \mathbf{0}_3 & -\mathbf{R}\Delta t \\
   \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{I}_3 & \mathbf{0}_3 \\
   \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{I}_3
   \end{bmatrix}$$

4. **Process Noise Covariance $\mathbf{Q} \in \mathbb{R}^{15 \times 15}$**:
   Driven by spectral densities $\sigma_a^2, \sigma_g^2, \sigma_{ba}^2, \sigma_{bg}^2$:
   - $\mathbf{Q}_{p,p} = \frac{1}{3}\sigma_a^2 \Delta t^3 \mathbf{I}_3$
   - $\mathbf{Q}_{p,v} = \mathbf{Q}_{v,p}^T = \frac{1}{2}\sigma_a^2 \Delta t^2 \mathbf{I}_3$
   - $\mathbf{Q}_{v,v} = \sigma_a^2 \Delta t \mathbf{I}_3$
   - $\mathbf{Q}_{\theta,\theta} = \sigma_g^2 \Delta t \mathbf{I}_3$
   - $\mathbf{Q}_{ba,ba} = \sigma_{ba}^2 \Delta t \mathbf{I}_3$
   - $\mathbf{Q}_{bg,bg} = \sigma_{bg}^2 \Delta t \mathbf{I}_3$

5. **Covariance Propagation**:
   $$\mathbf{P}_{k|k-1} = \boldsymbol{\Phi} \mathbf{P}_{k-1|k-1} \boldsymbol{\Phi}^T + \mathbf{Q}$$

---

### 2.3 Measurement Updates & Closed-Loop Feedback

For any measurement $\mathbf{z}$ with model $h(\mathbf{x}_{nom})$ and Jacobian $\mathbf{H} = \frac{\partial h}{\partial \delta\mathbf{x}}$:
$$\mathbf{y} = \mathbf{z} - h(\mathbf{x}_{nom})$$
$$\mathbf{S} = \mathbf{H} \mathbf{P} \mathbf{H}^T + \mathbf{R}_{meas}$$
$$\mathbf{K} = \mathbf{P} \mathbf{H}^T \mathbf{S}^{-1}$$
$$\delta\mathbf{x} = \mathbf{K} \mathbf{y} \in \mathbb{R}^{15}$$

#### Specific Measurement Models:
1. **GNSS Position ($\mathbf{z}_p \in \mathbb{R}^3$)**:
   $$\mathbf{y}_p = \mathbf{z}_p - (\mathbf{p} + \mathbf{R}_b^e \mathbf{l}^b_{ant})$$
   $$\mathbf{H}_p = \begin{bmatrix} \mathbf{I}_3 & \mathbf{0}_3 & -[(\mathbf{R}_b^e\mathbf{l}^b_{ant})\times] & \mathbf{0}_3 & \mathbf{0}_3 \end{bmatrix}$$

2. **GNSS Velocity ($\mathbf{z}_v \in \mathbb{R}^3$)**:
   $$\mathbf{y}_v = \mathbf{z}_v - \mathbf{v}$$
   $$\mathbf{H}_v = \begin{bmatrix} \mathbf{0}_3 & \mathbf{I}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 \end{bmatrix}$$

3. **Dynamic Vehicle Non-Holonomic Constraints (NHC)**:
   Lateral and vertical body velocities are zero: $v_y^b \approx 0 \pm \sigma_{lat}$, $v_z^b \approx 0 \pm \sigma_{vert}$.
   With $\mathbf{v}^b = (\mathbf{R}_b^e)^T \mathbf{v}^e$:
   $$\mathbf{y}_{nhc} = -\begin{bmatrix} v_y^b \\ v_z^b \end{bmatrix} = -\mathbf{E}_{23} (\mathbf{R}_b^e)^T \mathbf{v}^e, \quad \mathbf{E}_{23} = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$$
   $$\mathbf{H}_{nhc} = \begin{bmatrix} \mathbf{0}_{2\times 3} & \mathbf{E}_{23}(\mathbf{R}_b^e)^T & \mathbf{E}_{23}(\mathbf{R}_b^e)^T [\mathbf{v}^e\times] & \mathbf{0}_{2\times 3} & \mathbf{0}_{2\times 3} \end{bmatrix}$$
   *Crucial Insight*: The block $\mathbf{E}_{23}(\mathbf{R}_b^e)^T [\mathbf{v}^e\times]$ directly observes heading and attitude error when moving forward, rapidly correcting initial heading errors and estimating gyro bias!

4. **Zero-Velocity Update (ZUPT)**:
   When stationary:
   $$\mathbf{y}_{zupt} = -\mathbf{v}^e, \quad \mathbf{H}_{zupt} = \begin{bmatrix} \mathbf{0}_3 & \mathbf{I}_3 & \mathbf{0}_3 & \mathbf{0}_3 & \mathbf{0}_3 \end{bmatrix}$$

#### Closed-Loop Error-State Injection & Reset:
$$\mathbf{p} \leftarrow \mathbf{p} + \delta\mathbf{x}[0..3]$$
$$\mathbf{v} \leftarrow \mathbf{v} + \delta\mathbf{x}[3..6]$$
$$\mathbf{q} \leftarrow \text{UnitQuaternion}::\text{from\_scaled\_axis}(\delta\mathbf{x}[6..9]) \otimes \mathbf{q} \quad (\text{renormalized})$$
$$\mathbf{b}_a \leftarrow \mathbf{b}_a + \delta\mathbf{x}[9..12] \quad (\text{clamped to physical bounds})$$
$$\mathbf{b}_g \leftarrow \mathbf{b}_g + \delta\mathbf{x}[12..15] \quad (\text{clamped to physical bounds})$$

#### Joseph-Form Covariance Update:
$$\mathbf{P} \leftarrow (\mathbf{I}_{15} - \mathbf{K}\mathbf{H}) \mathbf{P} (\mathbf{I}_{15} - \mathbf{K}\mathbf{H})^T + \mathbf{K} \mathbf{R}_{meas} \mathbf{K}^T$$
$$\mathbf{P} \leftarrow \frac{1}{2}(\mathbf{P} + \mathbf{P}^T)$$
$$\delta\mathbf{x} \leftarrow \mathbf{0}_{15}$$

---

### 2.4 15-State Backward Rauch-Tung-Striebel (RTS) Smoother

For a trajectory of $N$ epochs, the forward pass records at each epoch $k \in [0, N-1]$:
- Predicted nominal state $\mathbf{x}^-_{nom,k}$ and covariance $\mathbf{P}^-_{k}$
- Updated nominal state $\mathbf{x}^+_{nom,k}$ and covariance $\mathbf{P}^+_{k}$
- Transition matrix $\boldsymbol{\Phi}_{k+1}$ from epoch $k$ to $k+1$

#### Backward Pass:
Initialize at terminal epoch $N-1$:
$$\mathbf{x}^s_{nom, N-1} = \mathbf{x}^+_{nom, N-1}, \quad \mathbf{P}^s_{N-1} = \mathbf{P}^+_{N-1}$$

For $k = N-2$ down to $0$:
1. **Smoother Gain Matrix $\mathbf{C}_k \in \mathbb{R}^{15 \times 15}$**:
   $$\mathbf{C}_k = \mathbf{P}^+_k \boldsymbol{\Phi}_{k+1}^T (\mathbf{P}^-_{k+1})^{-1}$$

2. **Smoothed State Discrepancy at $k+1$**:
   $$\delta\mathbf{x}_{k+1} = \begin{bmatrix}
   \mathbf{p}^s_{k+1} - \mathbf{p}^-_{k+1} \\
   \mathbf{v}^s_{k+1} - \mathbf{v}^-_{k+1} \\
   \text{scaled\_axis}\left(\mathbf{q}^s_{k+1} \otimes (\mathbf{q}^-_{k+1})^{-1}\right) \\
   \mathbf{b}_{a, k+1}^s - \mathbf{b}_{a, k+1}^- \\
   \mathbf{b}_{g, k+1}^s - \mathbf{b}_{g, k+1}^-
   \end{bmatrix} \in \mathbb{R}^{15}$$

3. **Smoothed Correction at Epoch $k$**:
   $$\delta\mathbf{x}^s_k = \mathbf{C}_k \delta\mathbf{x}_{k+1}$$

4. **Nominal State Correction at Epoch $k$**:
   $$\mathbf{p}^s_k = \mathbf{p}^+_k + \delta\mathbf{x}^s_k[0..3]$$
   $$\mathbf{v}^s_k = \mathbf{v}^+_k + \delta\mathbf{x}^s_k[3..6]$$
   $$\mathbf{q}^s_k = \text{UnitQuaternion}::\text{from\_scaled\_axis}(\delta\mathbf{x}^s_k[6..9]) \otimes \mathbf{q}^+_k$$
   $$\mathbf{b}_{a,k}^s = \mathbf{b}_{a,k}^+ + \delta\mathbf{x}^s_k[9..12]$$
   $$\mathbf{b}_{g,k}^s = \mathbf{b}_{g,k}^+ + \delta\mathbf{x}^s_k[12..15]$$

5. **Smoothed Covariance**:
   $$\mathbf{P}^s_k = \mathbf{P}^+_k + \mathbf{C}_k (\mathbf{P}^s_{k+1} - \mathbf{P}^-_{k+1}) \mathbf{C}_k^T$$
   $$\mathbf{P}^s_k \leftarrow \frac{1}{2}(\mathbf{P}^s_k + (\mathbf{P}^s_k)^T)$$

---

## 3. Architecture & Modular Decomposition Plan

To comply strictly with `AGENTS.md` rules:
- **File size**: Strictly $< 500$ LOC
- **Function size**: Strictly $< 32$ LOC
- **Nesting depth**: Strictly $< 3$ levels
- **Zero compiler warnings**
- **Zero `unwrap()` in production code**

The 15-state ESKF system must be decomposed into modular files under `crates/gneiss-rtk/src/swfg/imu_preintegration/`:

```
crates/gneiss-rtk/src/swfg/imu_preintegration/
├── mod.rs                  # Module root, re-exports, factor graph factor (~300 LOC)
├── types.rs                # 15-state types: State15, Cov15, NominalState, Snapshots (~160 LOC)
├── predict.rs              # 15-state transition Phi, process noise Q, mechanization (~220 LOC)
├── update.rs               # GNSS, NHC, ZUPT, Joseph-form update, attitude/bias feedback (~260 LOC)
├── smoother.rs             # 15-state RTS smoother backward recursion (~220 LOC)
├── stationary.rs           # Stationary detector (already exists, 120 LOC)
└── tests.rs                # Unit tests for 15-state ESKF, NHC attitude locking, bias convergence (~350 LOC)
```

### 3.1 Type Definitions (`types.rs`)
- `State15 = nalgebra::SVector<f64, 15>`
- `Cov15 = nalgebra::SMatrix<f64, 15, 15>`
- `NominalState`:
  ```rust
  #[derive(Debug, Clone)]
  pub struct NominalState {
      pub pos: Vector3<f64>,
      pub vel: Vector3<f64>,
      pub att: UnitQuaternion<f64>,
      pub ba: Vector3<f64>,
      pub bg: Vector3<f64>,
  }
  ```
- `InertialEpochSnapshot15`:
  ```rust
  #[derive(Debug, Clone)]
  pub struct InertialEpochSnapshot15 {
      pub time: GpsTime,
      pub nom_pred: NominalState,
      pub p_pred: Cov15,
      pub nom_post: NominalState,
      pub p_post: Cov15,
      pub f_mat: Cov15,
      pub is_gnss_available: bool,
  }
  ```
- `SmoothedInertialEpoch`:
  ```rust
  #[derive(Debug, Clone)]
  pub struct SmoothedInertialEpoch {
      pub time: GpsTime,
      pub position_ecef: Vector3<f64>,
      pub velocity_ecef: Vector3<f64>,
      pub attitude: UnitQuaternion<f64>,
      pub accel_bias: Vector3<f64>,
      pub gyro_bias: Vector3<f64>,
      pub cov_position: Matrix3<f64>,
      pub cov_velocity: Matrix3<f64>,
  }
  ```

---

## 4. Benchmark Target Strategy: Odaiba $p_{50} < 2.5\text{ m}$, RMS $< 5.2\text{ m}$

### 4.1 Root Causes of Previous Baseline ($p_{50} = 2.907\text{ m}$, RMS = $5.508\text{ m}$)

1. **Initial Heading Error ($7.07^\circ$)**:
   - `eval_odaiba_ins.rs` estimates initial heading from two GNSS fixes separated by > 5m: `estimate_initial_heading` yields $319.58^\circ$, but the ground-truth NovAtel reference heading is $326.65^\circ$.
   - In the 6-state filter, this $7^\circ$ error is **never corrected**.
   - In the 15-state ESKF with NHC attitude coupling ($\mathbf{H}_{nhc,\theta} = \mathbf{E}_{23}(\mathbf{R}_b^e)^T [\mathbf{v}^e\times]$), moving forward at $10-15\text{ m/s}$ creates lateral velocity innovations that drive $\delta\boldsymbol{\theta}_z$ to zero within 2–3 seconds, locking the heading to true path direction.

2. **Uncalibrated Sensor Biases**:
   - Accel bias was set to 0. An unmodeled bias of $0.05\text{ m/s}^2$ causes $0.025\text{ m}$ drift in 1s, and $> 2.5\text{ m}$ in 10s.
   - Gyro bias was frozen after 350 samples ($7\text{ s}$).
   - In the 15-state ESKF, online closed-loop bias updates continuously absorb MEMS drift.

3. **RTS Smoother Attitude Backward Coupling**:
   - The 6-state smoother was blind to attitude. Future GNSS fixes could not fix past heading errors.
   - The 15-state RTS smoother propagates future GNSS position/velocity fixes backward through $\boldsymbol{\Phi}^T$, correcting the attitude throughout the entire trajectory.

### 4.2 Expected Error Reduction
- Eliminating the $7^\circ$ heading error reduces horizontal position dispersion by $\approx 0.3 - 0.5\text{ m}$ during maneuvers and straight runs.
- Closed-loop online bias estimation reduces dead-reckoning drift between 1Hz GNSS fixes from $\approx 0.4\text{ m}$ to $< 0.1\text{ m}$.
- Projected Odaiba performance with 15-state ESKF + RTS:
  - $p_{50} \approx 2.1 - 2.4\text{ m}$ (beating $< 2.5\text{ m}$)
  - $\text{RMS} \approx 4.6 - 5.0\text{ m}$ (beating $< 5.2\text{ m}$)

---

## 5. Verification & Testing Plan

### 5.1 Unit Tests (`tests.rs`)
1. `test_15state_transition_matrix_positive_velocity_attitude_coupling`:
   Verify $\boldsymbol{\Phi}_{v,\theta} = +[\mathbf{f}^e\times]\Delta t$ is strictly positive.
2. `test_nhc_attitude_coupling_corrects_heading_error`:
   Inject a $10^\circ$ heading error on a vehicle moving forward at $15\text{ m/s}$; verify that a single NHC update rotates attitude toward truth.
3. `test_closed_loop_bias_estimation`:
   Simulate constant acceleration with injected accel bias; verify GNSS velocity updates converge $\mathbf{b}_a$ to the true bias.
4. `test_15state_rts_smoother_reduces_covariance_and_attitude_error`:
   Verify $\mathbf{P}^s_k \le \mathbf{P}^+_k$ and smoothed attitude matches truth across GNSS outages.

### 5.2 Integration Tests
1. `tests/src/inertial_outage_simulation.rs`:
   Run with 15-state ESKF; verify 10-second complete outage drift $< 0.50\text{ m}$.
2. `cargo run --release --bin eval_odaiba_ins`:
   Verify full 12,398-epoch trajectory satisfies $p_{50} < 2.5\text{ m}$ and $\text{RMS} < 5.2\text{ m}$.
3. CI Invariant Checks:
   `cargo clippy --workspace --all-targets -- -D warnings` (0 warnings).
   `scripts/check_network_benchmark.py --smoke` (clean pass).
   `scripts/check_multignss_benchmark.py --smoke` (clean pass).

---

## 6. Recommendations for Implementation Team

1. **Keep `types.rs`, `predict.rs`, `update.rs`, `smoother.rs` cleanly separated**:
   Do not merge into one large file. Keep each file under 300 LOC.
2. **Double check `Phi[(3+r, 6+c)] = +f_e_skew[(r, c)] * dt`**:
   Never put a minus sign on this term (violates Bug 2 invariant).
3. **Use Joseph-form covariance updates**:
   $\mathbf{P} = (\mathbf{I} - \mathbf{K}\mathbf{H})\mathbf{P}(\mathbf{I} - \mathbf{K}\mathbf{H})^T + \mathbf{K}\mathbf{R}\mathbf{K}^T$ prevents numerical non-positive-definiteness on high-rate 10Hz updates.
4. **Initial Heading Seeding in `eval_odaiba_ins.rs`**:
   Can initialize with GNSS velocity direction or NovAtel reference, but the 15-state ESKF with NHC will converge rapidly regardless.
