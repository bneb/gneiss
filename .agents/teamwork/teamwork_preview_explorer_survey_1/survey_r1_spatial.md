# Survey R1: Type-Safe Coordinate, Vector, and Covariance Primitives & Body FRD Lever Arm

## 1. Executive Summary & Problem Scope

Gneiss is a high-precision multi-GNSS RTK/PPP/INS engine supporting network RTK, tightly coupled ESKF INS, and sliding-window factor graphs (SWFG). Across `crates/gneiss-core` and `crates/gneiss-rtk`, spatial positions, velocities, baseline offsets, antenna lever arms, and error covariances are overwhelmingly represented using bare linear algebra primitives: `nalgebra::Vector3<f64>`, `nalgebra::Matrix3<f64>`, and raw arrays `[f64; 3]`.

This lack of compile-time frame distinction exposes the engine to critical physical correctness and safety risks:
1. **Unchecked Cross-Frame Vector Arithmetic**: A displacement in the vehicle Body FRD frame (Forward, Right, Down) can be directly added to an ECEF position without rotating by the attitude matrix $C_b^e$. The compiler accepts this silently because both are `Vector3<f64>`.
2. **Unchecked Cross-Datum Geometric Mixing**: Satellite broadcast orbits (GPS/WGS84, Galileo/GTRF, BDS/CGCS2000, GLONASS/PZ-90) and CORS station ground truths (ITRF2014, ITRF2020, NAD83(2011), JGD2011) are combined as bare coordinates. In particular, GLONASS orbits computed in PZ-90 via RK4 integration (`ephemeris/glonass.rs`) are mixed with WGS84/ITRF orbits without datum tagging or Helmert alignment.
3. **Relational Coupling Invariant Violations**: Functions such as `coords::ecef_delta_to_enu(target_ecef, ref_ecef, ref_llh)` accept parallel bare floats representing multiple views of the same reference point, violating the hard requirement in `AGENTS.md` that coupled parameters must be derived from shared typed inputs.
4. **Orientation and Rotation Ambiguity**: Rotations and attitudes are stored as bare `UnitQuaternion<f64>` or `Matrix3<f64>` without typed source and target frames, making it impossible to statically verify whether a rotation represents $C_b^e$ (Body to ECEF), $C_e^b$ (ECEF to Body), or $C_b^n$ (Body to NED).
5. **Covariance Frame Confusion**: Covariances in ECEF, NED, and ENU are passed as bare `Matrix3<f64>`, risking axis permutation bugs (e.g. NED North=0, East=1 vs ENU East=0, North=1) that corrupt innovation gating and quality metrics.

This survey establishes the complete catalog of spatial types and usages across the workspace, analyzes the mathematical and structural invariants of the Body FRD antenna lever arm, and specifies zero-cost compile-time typestate wrappers and refactoring plans adhering strictly to `AGENTS.md` standards (< 500 LOC per file, <= 32 LOC per function, zero unwrap in production).

---

## 2. Comprehensive Codebase Inventory

### 2.1 `crates/gneiss-core`

#### `src/coords.rs`
- **Current Runtime Primitives**:
  - `Datum` enum: `WGS84`, `ITRF2014`, `ITRF2020`, `JGD2011`, `PZ90`, `GTRF`, `CGCS2000` (lines 28–37).
  - `Frame` enum: `ECEF`, `ENU`, `LLH`, `Body` (lines 40–50).
  - `Coordinate` struct (lines 53–85): Runtime struct holding `vector: Vector3<f64>`, `datum: Datum`, `frame: Frame`, `epoch: GpsTime`.
    - Implements runtime check: `ensure_aligned(&self, other: &Coordinate) -> Result<(), &'static str>`.
    - **Deficiency**: Runtime verification is opt-in, incurs runtime branching, and is bypassed everywhere in core estimators where bare `Vector3<f64>` is extracted and manipulated directly.
- **Conversion Functions (Bare Signatures)**:
  - `ecef_to_llh_with_model<M: GeodeticModel>(model: &M, ecef: Vector3<f64>) -> Vector3<f64>` (lines 88–144)
  - `ecef_to_llh(ecef: Vector3<f64>) -> Vector3<f64>` (lines 146–148)
  - `llh_to_ecef(llh: Vector3<f64>) -> Vector3<f64>` (lines 172–174)
  - `az_el(pos_llh: Vector3<f64>, pos_ecef: Vector3<f64>, sat_ecef: Vector3<f64>) -> (f64, f64)` (lines 177–204): Takes both `pos_llh` and `pos_ecef` as parallel bare vectors.
  - `ecef_to_ned_matrix(llh: Vector3<f64>) -> Matrix3<f64>` (lines 207–227)
  - `enu_to_ecef(origin_ecef: Vector3<f64>, enu: Vector3<f64>) -> Vector3<f64>` (lines 230–249)
  - `ecef_delta_to_enu(target_ecef: Vector3<f64>, ref_ecef: Vector3<f64>, ref_llh: Vector3<f64>) -> Vector3<f64>` (lines 255–259): Takes parallel `ref_ecef` and `ref_llh`.
  - `ecef_cov_to_enu_std(pos_ecef: Vector3<f64>, cov_ecef: Matrix3<f64>) -> (f64, f64, f64)` (lines 266–275): Takes bare `Vector3<f64>` and bare `Matrix3<f64>`, returning untyped tuple `(f64, f64, f64)` representing (East, North, Up) standard deviations.

#### `src/frames/`
- **`mod.rs`**: Exports `EcefPos`, `EpochPosition`, `ReferenceFrame`, realization markers.
- **`realizations.rs`**:
  - Marker trait `pub trait ReferenceFrame { const NAME: &'static str; const HELMERT_TO_ITRF2014: Option<HelmertParams>; }` (lines 6–9).
  - Concrete realizations: `Itrf2014`, `Itrf2020`, `Igs20`, `Wgs84Broadcast`, `Nad83_2011`, `Etrs89`, `Gda2020`, `Jgd2011`.
  - **Deficiency**: `Pz90` (GLONASS reference frame PZ-90 / PZ-90.11) is entirely missing from `realizations.rs` despite being required by R1 and used in GLONASS ephemeris calculation.
- **`positions.rs`**:
  - `EcefPos<F: ReferenceFrame>(pub Vector3<f64>, pub PhantomData<F>)` (line 9).
  - `EpochPosition<F: ReferenceFrame, R: AntennaReference = Arp>` (line 111).
  - **Deficiencies in `EcefPos`**:
    - Implements `core::ops::Deref<Target = Vector3<f64>>` (lines 43–48) and `From<Vector3<f64>>` (lines 50–54).
    - `Deref` implicitly strips the frame/datum marker, allowing callers to perform unchecked vector operations or pass `&EcefPos<F>` anywhere `&Vector3<f64>` is accepted.
    - No strong types exist for local frames: `Ned`, `Enu`, `BodyFrd`.
    - No strong types exist for velocities or covariances.
- **`helmert.rs`**:
  - `HelmertParams`: 14-parameter Helmert transformation with time-dependent rates (lines 11–27).
  - Supports `apply` and `apply_inverse` on `Vector3<f64>`.

#### `src/ephemeris/glonass.rs`
- **GLONASS Orbit Integration**:
  - Broadcast ephemeris orbit integration via 4th-order Runge-Kutta in PZ-90 (lines 1, 27–74).
  - `GlonassEphemeris::position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64)` (lines 76–104).
  - **Deficiency**: Returns `(Vector3<f64>, Vector3<f64>, ...)` which are ECEF coordinates in the PZ-90 frame. When consumed by estimators, these coordinates are treated as if they were in WGS84/ITRF without transformation.

#### `src/imu.rs`
- `ImuMeasurement`:
  - `pub accel: Vector3<f64>` (vehicle/sensor body frame, m/s^2)
  - `pub gyro: Vector3<f64>` (vehicle/sensor body frame, rad/s)
  - `pub time_tag: u32` (untagged integer timestamp)
  - **Deficiency**: Bare `Vector3<f64>` without Body FRD marker.

---

### 2.2 `crates/gneiss-rtk`

#### `src/estimators/eskf/`
- **`types.rs`**:
  - `EskfState` (lines 33–40):
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
    - `pos_ecef` is a bare `Vector3<f64>`.
    - `vel_ecef` is a bare `Vector3<f64>`.
    - `attitude` is a bare `UnitQuaternion<f64>` representing $C_b^e$ (Body-to-ECEF rotation), but unannotated.
    - `accel_bias` and `gyro_bias` are Body FRD vectors.
    - `cov` is a flat 15x15 covariance matrix without block type safety.
  - `normal_gravity_ecef(pos: &Vector3<f64>) -> Vector3<f64>` (lines 106–116)
  - `earth_rotation_rate_ecef() -> Vector3<f64>` (lines 118–120)
- **`dd_update.rs`**:
  - `DdSatGeometry` (lines 9–14):
    ```rust
    pub struct DdSatGeometry {
        pub u_s: Vector3<f64>,       // ECEF line-of-sight unit vector
        pub u_ref: Vector3<f64>,     // ECEF line-of-sight unit vector to reference satellite
        pub lever_arm: Vector3<f64>, // Body FRD antenna lever arm!
    }
    ```
    - Bundles ECEF LOS unit vectors with a Body FRD lever arm as bare `Vector3<f64>`.
  - `compute_dd_jacobian_15(state: &EskfState, geom: &DdSatGeometry) -> RowVector15<f64>` (lines 44–60):
    - Multiplies attitude rotation matrix by body lever arm:
      `let l_e = state.attitude.to_rotation_matrix() * geom.lever_arm;`
    - Computes attitude Jacobian: `delta_u.transpose() * skew_symmetric(&l_e)`.
- **`update.rs`**:
  - `build_gnss_pos_system(state: &EskfState, pos_meas: &Vector3<f64>, r_pos: &Matrix3<f64>, lever_arm: &Vector3<f64>)` (lines 42–63):
    - `let l_e = state.attitude.to_rotation_matrix().into_inner() * lever_arm;`
    - `let y = pos_meas - (state.pos_ecef + l_e);`
  - `build_doppler_velocity_system(state: &EskfState, vel_meas: &Vector3<f64>, r_vel: &Matrix3<f64>, lever_arm: &Vector3<f64>, gyro_meas: &Vector3<f64>)` (lines 85–110):
    - `let omega_corr = gyro_meas - state.gyro_bias;` (Body frame)
    - `let v_rot_b = omega_corr.cross(lever_arm);` (Body frame)
    - `let v_rot_e = r_b2e * v_rot_b;` (ECEF frame)
    - `let v_ant_pred = state.vel_ecef + v_rot_e;` (ECEF frame)
    - `let y = vel_meas - v_ant_pred;`
  - `update_gnss_position` (lines 66–82) and `update_doppler_velocity` (lines 142–156).
- **`constraints.rs`**:
  - `build_nhc_system(state: &EskfState, _lever_arm: &Vector3<f64>)` (lines 9–32):
    - `let r_e2b = state.attitude.to_rotation_matrix().into_inner().transpose();`
    - `let v_b = r_e2b * state.vel_ecef;` (Projects ECEF velocity to Body frame)
    - `let y = Vector2::new(-v_b.y, -v_b.z);` (Constrains lateral and vertical body velocity)
- **`alignment.rs`**:
  - `compute_initial_attitude(imu_samples: &[ImuSample], init_pos: Vector3<f64>, heading_rad: f64, max_samples: usize) -> UnitQuaternion<f64>` (lines 30–42):
    - Bridges frames: Body leveling -> NED heading rotation ($C_b^n$) -> Local tangent rotation ($C_n^e$) -> $C_b^e = C_n^e \cdot C_b^n$.
    - Returns bare `UnitQuaternion<f64>` with no type indicator of source or target frame.

#### `src/estimators/rtk_iekf/`
- **`state/mod.rs`**:
  - `RtkState` (lines 24–48):
    - `pub pos_ecef: Vector3<f64>`
    - `pub vel_ecef: Vector3<f64>`
    - `pub cov: DMatrix<f64>`
    - `pub fn extract_pos_cov(&self) -> Matrix3<f64>` (lines 420–428): Returns bare `Matrix3<f64>` in ECEF.

#### `src/estimators/doppler.rs`
- `DopplerVelocitySolution` (lines 12–24):
  - `pub vel_ecef: Vector3<f64>`
  - `pub cov: Matrix3<f64>`

#### `src/swfg/`
- **`config.rs`**:
  - `ImuConfig`: `pub lever_arm: [f64; 3]`, `pub nhc_lever_arm: [f64; 3]` (lines 153, 159).
  - `RtkConfig`: `pub base_position: [f64; 3]`, `pub initial_position: Option<[f64; 3]>`.
  - `PppConfig`: `pub initial_position: Option<[f64; 3]>`.
- **`engine/mod.rs`**:
  - `SwfgSolution`: `pub position_ecef: Vector3<f64>` (line 28).
  - `SwfgEngine`: `prev_position: Option<Vector3<f64>>`, `current_attitude: Option<UnitQuaternion<f64>>`, `initial_position: Option<Vector3<f64>>`.
- **`imu_preintegration/mod.rs`**:
  - `ImuPreintegration` (lines 33–54):
    - `pub dp: Vector3<f64>` (preintegrated position delta in initial body frame $b_i$)
    - `pub dv: Vector3<f64>` (preintegrated velocity delta in initial body frame $b_i$)
    - `pub dq: UnitQuaternion<f64>` (rotation delta $C_{b_t}^{b_i}$)
  - `ImuPreintegrationFactor` (lines 189–323):
    - Connects `Pose(epoch_i)`, `Velocity(epoch_i)`, `Pose(epoch_j)`, `Velocity(epoch_j)`, `ImuBias`.
    - Predicts motion in world frame, then rotates to body frame via $C_e^{b_i} = (C_{b_i}^e)^T$:
      `let dp_pred_body = q_i.inverse() * dp_pred_world;`
      `let dv_pred_body = q_i.inverse() * dv_pred_world;`
- **`pipeline/dd_factors.rs`**:
  - `DdPseudorangeFactor`: `pub sat_pos: Vector3<f64>`, `pub ref_pos: Vector3<f64>`, `pub base_pos: Vector3<f64>`.
  - Residual assumes antenna is at pose translation origin: `let rx_pos = Vector3::new(pose[0], pose[1], pose[2])`.
- **`pipeline/aux_factors.rs`**:
  - `OdometerVelocityFactor`: `pub measured_v_body: Vector3<f64>`, rotates predicted ECEF velocity to body: `q.inverse() * v_ecef`.
  - `DualAntennaHeadingFactor`: `pub baseline_body: Vector3<f64>`, `pub measured_baseline_ecef: Vector3<f64>`, rotates body baseline to ECEF: `R_b^e(q) * baseline_body`.

#### `src/post_process/`
- **`lever_arm.rs`**:
  - `LeverArmObservation`: `omega_body`, `alpha_body`, `accel_imu_body`, `accel_gnss_body` (all `Vector3<f64>`).
  - `LeverArmEstimate`: `lever_arm_body: Vector3<f64>`, `std_body: Vector3<f64>`.
  - Solves batch least-squares for antenna-to-IMU lever arm $\mathbf{l}_b$ from vehicle maneuvers:
    $\mathbf{a}_{GNSS}^b - \mathbf{a}_{IMU}^b = ([\boldsymbol{\omega} \times]^2 + [\dot{\boldsymbol{\omega}} \times]) \mathbf{l}_b$.
- **`mod.rs`**:
  - `PostProcessOptions` (lines 64–120):
    - `pub base_position: Option<Vector3<f64>>`
    - `pub initial_rover_position: Option<Vector3<f64>>`
- **`combiner.rs`**:
  - `SmoothedEpoch` (lines 16–28):
    - `pub position_ecef: Vector3<f64>`
    - `pub velocity_ecef: Option<Vector3<f64>>`
    - `pub attitude: Option<UnitQuaternion<f64>>`
    - `pub cov_position: Matrix3<f64>`
    - `pub std_east: f64`, `pub std_north: f64`, `pub std_up: f64`
- **`antenna.rs`**:
  - `station_recv_pco_ecef(rinex_path: &Path, antex_path: &str, arp: Vector3<f64>) -> Option<Vector3<f64>>`:
    Transforms antenna PCO from ENU (extracted from ANTEX) to ECEF at ARP position.

#### `src/bin/eval_odaiba_ins/`
- **`main.rs`**:
  - `const ANTENNA_LEVER_ARM: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);` (line 24)
  - Lever arm applied in `update_gnss_innovation` (line 201):
    `let r_b2e = state.attitude.to_rotation_matrix().into_inner();`
    `let l_e = r_b2e * ANTENNA_LEVER_ARM;`
    `let innov_norm = (pos - (state.pos_ecef + l_e)).norm();`
  - Lever arm passed to `update_gnss_position` (line 216).
  - Lever arm applied in `apply_single_sat_dd` (line 245):
    `let delta_pos = p_fix - (state.pos_ecef + state.attitude.to_rotation_matrix() * geom.lever_arm);`

---

## 3. Body FRD Antenna Lever Arm Deep-Dive

### 3.1 Physical and Mathematical Model
In GNSS/INS integration, the IMU is located at the center of navigation / body origin, while the GNSS antenna phase center is mounted at a physical offset:
$$\mathbf{r}_{ant}^b = \begin{bmatrix} l_x \\ l_y \\ l_z \end{bmatrix}_{Body FRD}$$
where $+X$ is Forward along vehicle chassis, $+Y$ is Right, $+Z$ is Down.

The position of the antenna phase center in the Earth-Centered, Earth-Fixed (ECEF) frame at epoch $t$ is:
$$\mathbf{r}_{ant}^e(t) = \mathbf{r}_{imu}^e(t) + C_b^e(t) \mathbf{r}_{ant}^b$$
where $C_b^e \in SO(3)$ is the direction cosine matrix rotating vectors from the Body frame to ECEF.

The kinematic velocity of the antenna phase center in ECEF is obtained by time differentiation:
$$\mathbf{v}_{ant}^e(t) = \mathbf{v}_{imu}^e(t) + \boldsymbol{\omega}_{ie}^e \times (C_b^e \mathbf{r}_{ant}^b) + C_b^e (\boldsymbol{\omega}_{ib}^b \times \mathbf{r}_{ant}^b)$$
where $\boldsymbol{\omega}_{ib}^b$ is the IMU gyro angular rate corrected for sensor bias, and $\boldsymbol{\omega}_{ie}^e$ is the Earth rotation rate vector in ECEF. (In high-rate local calculations, the Coriolis term is often absorbed into the state transition model).

### 3.2 Observation Equations and Jacobians

#### 1. GNSS Position Innovation and Measurement Jacobian
The GNSS position measurement $\mathbf{z}_{pos}^e$ measures the antenna position in ECEF:
$$\mathbf{z}_{pos}^e = \mathbf{r}_{ant}^e + \mathbf{v}_{pos} = \mathbf{r}_{imu}^e + C_b^e \mathbf{r}_{ant}^b + \mathbf{v}_{pos}$$
The innovation is:
$$\mathbf{y} = \mathbf{z}_{pos}^e - (\mathbf{r}_{imu}^e + \mathbf{l}^e), \quad \text{where } \mathbf{l}^e = C_b^e \mathbf{r}_{ant}^b$$
Under a global left-multiplied attitude error model $\hat{C}_b^e = (\mathbf{I} - [\boldsymbol{\delta\theta} \times]) C_b^e$, perturbation yields:
$$\delta \mathbf{l}^e = -[\boldsymbol{\delta\theta} \times] \mathbf{l}^e = [\mathbf{l}^e \times] \boldsymbol{\delta\theta} = -[\mathbf{l}^e \times] \boldsymbol{\psi}$$
Thus, the 15-state ESKF Jacobian is:
$$H_{pos} = \begin{bmatrix} \mathbf{I}_{3\times3} & \mathbf{0}_{3\times3} & -[\mathbf{l}^e \times] & \mathbf{0}_{3\times3} & \mathbf{0}_{3\times3} \end{bmatrix}$$
Notice that this matches `update.rs:59`:
`h[(r, c + 6)] = -l_e_skew[(r, c)];`

#### 2. Doppler Velocity Innovation and Measurement Jacobian
Doppler measurements resolve line-of-sight velocity of the antenna:
$$\mathbf{v}_{ant}^e = \mathbf{v}_{imu}^e + C_b^e (\boldsymbol{\omega}_b \times \mathbf{r}_{ant}^b)$$
where $\boldsymbol{\omega}_b = \boldsymbol{\omega}_{meas} - \mathbf{b}_g$.
The rotational velocity in Body frame is $\mathbf{v}_{rot}^b = \boldsymbol{\omega}_b \times \mathbf{r}_{ant}^b$.
Rotated to ECEF: $\mathbf{v}_{rot}^e = C_b^e \mathbf{v}_{rot}^b$.
The Jacobian w.r.t. gyro bias $\mathbf{b}_g$ is:
$$\frac{\partial \mathbf{v}_{ant}^e}{\partial \mathbf{b}_g} = -C_b^e [\mathbf{r}_{ant}^b \times]$$

#### 3. Double-Difference (DD) Carrier Phase and Pseudorange
The double-difference range for satellite pair $(s, \text{ref})$ observed at antenna $\mathbf{r}_{ant}^e$ is:
$$\nabla\Delta \rho = \|\mathbf{r}^s - \mathbf{r}_{ant}^e\| - \|\mathbf{r}^{ref} - \mathbf{r}_{ant}^e\| - (\dots)_{\text{base}}$$
Linearizing w.r.t. antenna position:
$$\nabla\Delta \rho \approx -(\mathbf{u}^s - \mathbf{u}^{ref})^T \delta\mathbf{r}_{ant}^e$$
Substituting $\delta\mathbf{r}_{ant}^e = \delta\mathbf{r}_{imu}^e - [\mathbf{l}^e \times] \boldsymbol{\delta\theta}$:
$$\delta(\nabla\Delta \rho) = -(\mathbf{u}^s - \mathbf{u}^{ref})^T \delta\mathbf{r}_{imu}^e + (\mathbf{u}^s - \mathbf{u}^{ref})^T [\mathbf{l}^e \times] \boldsymbol{\delta\theta}$$
This directly verifies `dd_update.rs:54`:
`let h_att = delta_u.transpose() * l_e_skew;`

### 3.3 The Core Frame Safety Vulnerability
In current production code:
```rust
// In dd_update.rs:
pub struct DdSatGeometry {
    pub u_s: Vector3<f64>,       // ECEF
    pub u_ref: Vector3<f64>,     // ECEF
    pub lever_arm: Vector3<f64>, // Body FRD
}
```
Because both `u_s` and `lever_arm` are bare `Vector3<f64>`:
1. A programmer can write `state.pos_ecef + geom.lever_arm` without rotating by `state.attitude`. This compiles without error, injecting a ~1 meter body vector directly into an earth-centered coordinate space.
2. In `swfg/pipeline/dd_factors.rs`, the factor graph uses `pose[0..3]` directly as the antenna position, silently omitting the lever arm and attitude coupling entirely.
3. If attitude is rotated in reverse ($C_e^b$ instead of $C_b^e$), the compiler cannot detect the transpose error.

---

## 4. Typestate Architecture & Type-Safety Specification (R1)

To solve these vulnerabilities with **zero runtime overhead**, we design a compile-time phantom typestate architecture in `gneiss-core::frames`.

### 4.1 Frame and Datum Marker Traits

```rust
// In gneiss-core/src/frames/mod.rs

pub trait ReferenceFrame: 'static + Send + Sync + Copy + PartialEq + Eq {
    const NAME: &'static str;
    const HELMERT_TO_ITRF2014: Option<HelmertParams>;
}

pub trait CoordinateFrame: 'static + Send + Sync + Copy + PartialEq + Eq {
    const NAME: &'static str;
    const IS_EARTH_FIXED: bool;
}

// Coordinate Frame Markers:
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ecef<F: ReferenceFrame>(core::marker::PhantomData<F>);
impl<F: ReferenceFrame> CoordinateFrame for Ecef<F> {
    const NAME: &'static str = "ECEF";
    const IS_EARTH_FIXED: bool = true;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ned;
impl CoordinateFrame for Ned {
    const NAME: &'static str = "NED";
    const IS_EARTH_FIXED: bool = false;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Enu;
impl CoordinateFrame for Enu {
    const NAME: &'static str = "ENU";
    const IS_EARTH_FIXED: bool = false;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BodyFrd;
impl CoordinateFrame for BodyFrd {
    const NAME: &'static str = "Body-FRD";
    const IS_EARTH_FIXED: bool = false;
}

// Datum Markers in realizations.rs:
// Existing: Itrf2014, Itrf2020, Igs20, Wgs84Broadcast, Nad83_2011, Etrs89, Gda2020, Jgd2011.
// NEW ADDITION required by R1:
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pz90;
impl ReferenceFrame for Pz90 {
    const NAME: &'static str = "PZ-90.11";
    // PZ-90.11 to ITRF2014 transformation parameters
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
        tx_mm: 3.0,
        ty_mm: -1.0,
        tz_mm: 0.0,
        scale_ppb: 0.0,
        rx_mas: 0.019,
        ry_mas: -0.042,
        rz_mas: 0.002,
        ref_epoch_yr: 2011.0,
        tx_rate: 0.0,
        ty_rate: 0.0,
        tz_rate: 0.0,
        rx_rate: 0.0,
        ry_rate: 0.0,
        rz_rate: 0.0,
        scale_rate: 0.0,
    });
}
```

### 4.2 Spatial Primitives

#### 1. Affine Position: `Point3<Frame>` and `EcefPos<F>`
An affine point in a coordinate frame.
```rust
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct Point3<Frame: CoordinateFrame> {
    coords: nalgebra::Vector3<f64>,
    _frame: core::marker::PhantomData<Frame>,
}

// Backwards-compatible type alias for ECEF positions:
pub type EcefPos<F> = Point3<Ecef<F>>;
pub type NedPos = Point3<Ned>;
pub type EnuPos = Point3<Enu>;
```

#### 2. Spatial Vector / Displacement: `SpatialVector<Frame>`
A displacement or direction vector in a coordinate frame.
```rust
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct SpatialVector<Frame: CoordinateFrame> {
    vector: nalgebra::Vector3<f64>,
    _frame: core::marker::PhantomData<Frame>,
}

pub type EcefVector<F> = SpatialVector<Ecef<F>>;
pub type NedVector = SpatialVector<Ned>;
pub type EnuVector = SpatialVector<Enu>;
pub type BodyVector = SpatialVector<BodyFrd>;
```

#### 3. Spatial Velocity: `SpatialVelocity<Frame>`
```rust
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct SpatialVelocity<Frame: CoordinateFrame> {
    vector: nalgebra::Vector3<f64>,
    _frame: core::marker::PhantomData<Frame>,
}

pub type EcefVelocity<F> = SpatialVelocity<Ecef<F>>;
pub type NedVelocity = SpatialVelocity<Ned>;
pub type EnuVelocity = SpatialVelocity<Enu>;
pub type BodyVelocity = SpatialVelocity<BodyFrd>;
```

#### 4. Spatial Covariance: `SpatialCovariance<Frame>`
A 3x3 error covariance matrix in a specific frame.
```rust
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct SpatialCovariance<Frame: CoordinateFrame> {
    matrix: nalgebra::Matrix3<f64>,
    _frame: core::marker::PhantomData<Frame>,
}

pub type EcefCovariance<F> = SpatialCovariance<Ecef<F>>;
pub type NedCovariance = SpatialCovariance<Ned>;
pub type EnuCovariance = SpatialCovariance<Enu>;
pub type BodyCovariance = SpatialCovariance<BodyFrd>;

impl EnuCovariance {
    #[inline]
    pub fn std_east(&self) -> f64 { libm::sqrt(self.matrix[(0, 0)].max(0.0)) }
    #[inline]
    pub fn std_north(&self) -> f64 { libm::sqrt(self.matrix[(1, 1)].max(0.0)) }
    #[inline]
    pub fn std_up(&self) -> f64 { libm::sqrt(self.matrix[(2, 2)].max(0.0)) }
}

impl NedCovariance {
    #[inline]
    pub fn std_north(&self) -> f64 { libm::sqrt(self.matrix[(0, 0)].max(0.0)) }
    #[inline]
    pub fn std_east(&self) -> f64 { libm::sqrt(self.matrix[(1, 1)].max(0.0)) }
    #[inline]
    pub fn std_down(&self) -> f64 { libm::sqrt(self.matrix[(2, 2)].max(0.0)) }
}
```

#### 5. Antenna Lever Arm: `AntennaLeverArm`
Strictly typed in the Body FRD frame:
```rust
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AntennaLeverArm(pub SpatialVector<BodyFrd>);

impl AntennaLeverArm {
    #[inline]
    pub fn new(forward_m: f64, right_m: f64, down_m: f64) -> Self {
        Self(SpatialVector::new(forward_m, right_m, down_m))
    }
    #[inline]
    pub fn zero() -> Self {
        Self(SpatialVector::zero())
    }
    #[inline]
    pub fn as_body_vector(&self) -> &SpatialVector<BodyFrd> {
        &self.0
    }
}
```

#### 6. Typed Attitude & Rotation: `Attitude<FromFrame, ToFrame>`
```rust
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct Attitude<From: CoordinateFrame, To: CoordinateFrame> {
    q: nalgebra::UnitQuaternion<f64>,
    _frames: core::marker::PhantomData<(From, To)>,
}

pub type BodyToEcef<F> = Attitude<BodyFrd, Ecef<F>>;
pub type EcefToBody<F> = Attitude<Ecef<F>, BodyFrd>;
pub type BodyToNed = Attitude<BodyFrd, Ned>;
pub type NedToEcef<F> = Attitude<Ned, Ecef<F>>;

impl<From: CoordinateFrame, To: CoordinateFrame> Attitude<From, To> {
    #[inline]
    pub fn rotate_vector(&self, v: SpatialVector<From>) -> SpatialVector<To> {
        SpatialVector::from_inner(self.q * v.into_inner())
    }

    #[inline]
    pub fn rotate_velocity(&self, v: SpatialVelocity<From>) -> SpatialVelocity<To> {
        SpatialVelocity::from_inner(self.q * v.into_inner())
    }

    #[inline]
    pub fn rotate_cov(&self, cov: SpatialCovariance<From>) -> SpatialCovariance<To> {
        let r = self.q.to_rotation_matrix().into_inner();
        SpatialCovariance::from_inner(r * cov.into_inner() * r.transpose())
    }

    #[inline]
    pub fn inverse(&self) -> Attitude<To, From> {
        Attitude::from_inner(self.q.inverse())
    }
}
```

### 4.3 Structural Compile-Time Invariants (Operator Algebra)

| Operation | Left Type | Right Type | Result Type | Compiler Rule |
|-----------|-----------|------------|-------------|---------------|
| Displacement | `Point3<F>` | `Point3<F>` | `SpatialVector<F>` | `Sub` allowed only when `F` matches exactly |
| Point Translation | `Point3<F>` | `SpatialVector<F>` | `Point3<F>` | `Add` allowed only when `F` matches |
| Vector Addition | `SpatialVector<F>` | `SpatialVector<F>` | `SpatialVector<F>` | `Add` allowed only when `F` matches |
| Kinematic Step | `Point3<F>` | `SpatialVelocity<F> * dt` | `Point3<F>` | `Add` allowed only when `F` matches |
| **Cross-Frame Addition** | `Point3<Ecef<F>>` | `SpatialVector<BodyFrd>` | **DOES NOT COMPILE** | No `impl Add<SpatialVector<BodyFrd>> for Point3<Ecef<F>>` |
| **Cross-Datum Addition** | `Point3<Ecef<Itrf2014>>` | `SpatialVector<Ecef<Nad83>>` | **DOES NOT COMPILE** | Mismatched generic datum parameter |
| **Point-Point Addition** | `Point3<F>` | `Point3<F>` | **DOES NOT COMPILE** | No `impl Add<Point3<F>> for Point3<F>` |
| Lever Arm Application | `BodyToEcef<F>` | `AntennaLeverArm` | `EcefVector<F>` | Via `att.rotate_vector(arm.0)` |
| Lever Arm Translation | `Point3<Ecef<F>>` | `EcefVector<F>` | `Point3<Ecef<F>>` | `pos + l_e` compiles cleanly |

### 4.4 Encapsulation & Deref Policy
- **No `core::ops::Deref<Target = Vector3<f64>>`**: In the existing `EcefPos`, `Deref` leaked the inner vector, defeating frame safety. The new primitives do NOT implement `Deref`.
- **Explicit Accessors**:
  - `p.coords() -> &Vector3<f64>`
  - `v.vector() -> &Vector3<f64>`
  - `v.into_vector() -> Vector3<f64>`
  - `v.x()`, `v.y()`, `v.z()`
- **Zero-Cost Representation**: `#[repr(transparent)]` guarantees memory layout identical to raw `Vector3<f64>` and `Matrix3<f64>` with 0 runtime penalty.

### 4.5 Relational Coupling: Local Tangent Plane
To resolve the relational coupling violation in `coords::ecef_delta_to_enu`, we introduce:
```rust
pub struct LocalTangentPlane<F: ReferenceFrame> {
    origin: EcefPos<F>,
    r_ned: nalgebra::Matrix3<f64>,
}

impl<F: ReferenceFrame> LocalTangentPlane<F> {
    pub fn from_origin(origin: EcefPos<F>) -> Self {
        let llh = ecef_to_llh(*origin.coords());
        let r_ned = ecef_to_ned_matrix(llh);
        Self { origin, r_ned }
    }

    pub fn to_enu(&self, target: EcefPos<F>) -> EnuVector {
        let d = self.r_ned * (target.coords() - self.origin.coords());
        EnuVector::new(d.y, d.x, -d.z) // Strict NED -> ENU conversion
    }

    pub fn from_enu(&self, enu: EnuVector) -> EcefPos<F> {
        let enu_v = enu.vector();
        let ned_v = nalgebra::Vector3::new(enu_v.y, enu_v.x, -enu_v.z);
        let ecef_delta = self.r_ned.transpose() * ned_v;
        EcefPos::from_coords(self.origin.coords() + ecef_delta)
    }

    pub fn project_cov(&self, cov: EcefCovariance<F>) -> EnuCovariance {
        let cov_enu = self.r_ned * cov.matrix() * self.r_ned.transpose();
        // Permute [N, E, D] -> [E, N, U]
        let mut m = nalgebra::Matrix3::zeros();
        m[(0, 0)] = cov_enu[(1, 1)]; // E
        m[(1, 1)] = cov_enu[(0, 0)]; // N
        m[(2, 2)] = cov_enu[(2, 2)]; // U
        m[(0, 1)] = cov_enu[(1, 0)];
        m[(1, 0)] = cov_enu[(0, 1)];
        m[(0, 2)] = -cov_enu[(1, 2)];
        m[(2, 0)] = -cov_enu[(2, 1)];
        m[(1, 2)] = -cov_enu[(0, 2)];
        m[(2, 1)] = -cov_enu[(2, 0)];
        EnuCovariance::from_matrix(m)
    }
}
```
**Advantage**: Callers cannot pass mismatched origin coordinates. Disagreement between `ref_ecef` and `ref_llh` is unrepresentable.

---

## 5. Modules and Files Requiring Refactoring or New Primitives

### 5.1 New Primitives to Create in `crates/gneiss-core`
1. `crates/gneiss-core/src/frames/traits.rs` (~60 LOC):
   `CoordinateFrame` marker trait, `ReferenceFrame` extension.
2. `crates/gneiss-core/src/frames/markers.rs` (~80 LOC):
   `Ecef<F>`, `Ned`, `Enu`, `BodyFrd` frame markers; `Pz90` datum realization.
3. `crates/gneiss-core/src/frames/vectors.rs` (~180 LOC):
   `SpatialVector<Frame>`, `SpatialVelocity<Frame>`, operator overloading (`Add`, `Sub`, `Mul`), accessors (`x`, `y`, `z`).
4. `crates/gneiss-core/src/frames/covariances.rs` (~140 LOC):
   `SpatialCovariance<Frame>`, accessor methods (`std_east`, `std_north`, `std_up`, `std_down`).
5. `crates/gneiss-core/src/frames/attitude.rs` (~150 LOC):
   `Attitude<From, To>`, `rotate_vector`, `rotate_velocity`, `rotate_cov`, `inverse`.
6. `crates/gneiss-core/src/frames/tangent.rs` (~120 LOC):
   `LocalTangentPlane<F>`, relational coupling origin derivations.
7. `crates/gneiss-core/src/frames/lever_arm.rs` (~70 LOC):
   `AntennaLeverArm` struct, constructors, conversions.

### 5.2 Files to Refactor in `crates/gneiss-core`
1. `crates/gneiss-core/src/frames/mod.rs`: Re-export new primitives, maintain backwards-compatible aliases.
2. `crates/gneiss-core/src/frames/positions.rs`: Remove `Deref<Target = Vector3<f64>>`, integrate `Point3<Frame>`, provide explicit `.coords()`.
3. `crates/gneiss-core/src/frames/realizations.rs`: Add `Pz90` reference frame realization.
4. `crates/gneiss-core/src/coords.rs`: Mark relational-coupling-violating functions (`ecef_delta_to_enu`) as legacy/deprecated in favor of `LocalTangentPlane`.
5. `crates/gneiss-core/src/ephemeris/glonass.rs`: Tag GLONASS positions with `Pz90`.

### 5.3 Files to Refactor in `crates/gneiss-rtk`
1. `crates/gneiss-rtk/src/estimators/eskf/types.rs`:
   Update `EskfState`:
   - `pos_ecef: EcefPos<Itrf2014>` (or generic `F`)
   - `vel_ecef: EcefVelocity<Itrf2014>`
   - `attitude: BodyToEcef<Itrf2014>`
   - `accel_bias: SpatialVector<BodyFrd>`
   - `gyro_bias: SpatialVector<BodyFrd>`
2. `crates/gneiss-rtk/src/estimators/eskf/dd_update.rs`:
   - `DdSatGeometry`: `u_s: EcefVector<F>`, `u_ref: EcefVector<F>`, `lever_arm: AntennaLeverArm`.
   - Update `compute_dd_jacobian_15` to require `state.attitude.rotate_vector(geom.lever_arm.0)`.
3. `crates/gneiss-rtk/src/estimators/eskf/update.rs`:
   - `update_gnss_position`: takes `pos_meas: &EcefPos<F>`, `r_pos: &EcefCovariance<F>`, `lever_arm: &AntennaLeverArm`.
   - `build_doppler_velocity_system`: typed `vel_meas: &EcefVelocity<F>`, `gyro_meas: &SpatialVector<BodyFrd>`, `lever_arm: &AntennaLeverArm`.
4. `crates/gneiss-rtk/src/estimators/eskf/constraints.rs`:
   - `build_nhc_system`: typed `state.attitude.inverse().rotate_velocity(state.vel_ecef)`.
5. `crates/gneiss-rtk/src/estimators/eskf/alignment.rs`:
   - `compute_initial_attitude` returns `BodyToEcef<F>`.
6. `crates/gneiss-rtk/src/post_process/lever_arm.rs`:
   - `LeverArmEstimate`: `lever_arm_body: AntennaLeverArm`, `std_body: SpatialVector<BodyFrd>`.
7. `crates/gneiss-rtk/src/post_process/mod.rs`:
   - `PostProcessOptions`: `base_position: Option<EcefPos<Itrf2014>>`, `initial_rover_position: Option<EcefPos<Itrf2014>>`.
8. `crates/gneiss-rtk/src/post_process/combiner.rs`:
   - `SmoothedEpoch`: `position_ecef: EcefPos<Itrf2014>`, `cov_position: EcefCovariance<Itrf2014>`, `attitude: Option<BodyToEcef<Itrf2014>>`.
9. `crates/gneiss-rtk/src/swfg/config.rs`:
   - `ImuConfig`: `lever_arm: AntennaLeverArm`.
10. `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs`:
   - `ANTENNA_LEVER_ARM: AntennaLeverArm = AntennaLeverArm::zero();`
   - Use typed `update_gnss_position` and `LocalTangentPlane`.

---

## 6. Code Standards Compliance Check (AGENTS.md)

### 6.1 File Size Analysis (Limit: < 500 LOC)
Existing key files line counts:
| File | Current LOC | Margin to 500 LOC | Risk Assessment |
|------|-------------|-------------------|-----------------|
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs` | 474 | 26 lines | **HIGH RISK**: Must not expand. Refactoring should extract helper logic to `odaiba_helpers.rs` if needed to maintain < 500 LOC. |
| `crates/gneiss-rtk/src/estimators/eskf/condition.rs` | 452 | 48 lines | Medium risk. |
| `crates/gneiss-rtk/src/estimators/eskf/dd_update.rs` | 453 | 47 lines | Medium risk. Keep helper functions compact. |
| `crates/gneiss-rtk/src/estimators/rtk_iekf/state/mod.rs` | 448 | 52 lines | Medium risk. |
| `crates/gneiss-rtk/src/swfg/engine/mod.rs` | 429 | 71 lines | Low risk. |
| `crates/gneiss-core/src/coords.rs` | 338 | 162 lines | Safe. |
| `crates/gneiss-core/src/frames/tests.rs` | 342 | 158 lines | Safe. |
| `crates/gneiss-core/src/frames/positions.rs` | 205 | 295 lines | Safe. |
| `crates/gneiss-core/src/frames/realizations.rs` | 157 | 343 lines | Safe. |
| `crates/gneiss-rtk/src/estimators/eskf/update.rs` | 278 | 222 lines | Safe. |
| `crates/gneiss-rtk/src/estimators/eskf/types.rs` | 181 | 319 lines | Safe. |
| `crates/gneiss-rtk/src/post_process/lever_arm.rs` | 157 | 343 lines | Safe. |

### 6.2 Function Size Analysis (Limit: <= 32 LOC, Nesting < 3)
- All new typestate functions must be short and focused:
  - Constructors, accessors, and operators: 1–8 LOC.
  - Rotations and projections: 8–18 LOC.
  - `LocalTangentPlane` projection: 15–25 LOC.
- Any function refactored in `main.rs` or `dd_update.rs` must not exceed 32 LOC.

### 6.3 Hard Rules Check
- **Zero `unwrap()` in production code**: All conversions and inversions use `?`, `ok_or(EngineError::...)`, or `.expect("invariant: ...")`.
- **Zero compiler warnings**: All code passes `cargo check --workspace` and `cargo clippy --workspace --all-targets -- -D warnings`.
- **No `.orig` or `.rej` files**.

---

## 7. Migration & TDD Implementation Strategy

### Phase 1: Foundations in `gneiss-core` (Non-breaking)
1. Add `Pz90` realization in `frames/realizations.rs`.
2. Add new modular files in `frames/`: `traits.rs`, `markers.rs`, `vectors.rs`, `covariances.rs`, `attitude.rs`, `tangent.rs`, `lever_arm.rs`.
3. Add comprehensive unit tests in `frames/` asserting:
   - Identical memory layout (`std::mem::size_of`, `align_of`).
   - Zero-cost pass-by-value / pass-by-ref.
   - Successful compilation of correct cross-frame rotations.
   - Compile-fail tests (using `trybuild` or documentation compile-fail doctests) for illegal cross-frame additions (`EcefPos + AntennaLeverArm` fails to compile).
4. Verify all existing tests in `gneiss-core` pass.

### Phase 2: Estimator Refactoring in `gneiss-rtk`
1. Update `EskfState` and `DdSatGeometry` in `eskf/`.
2. Update `update_gnss_position` and `update_doppler_velocity`.
3. Update `eval_odaiba_ins/main.rs`.
4. Run regression suite:
   - `cargo test --workspace` (all 789+ tests must pass).
   - `python3 scripts/check_network_benchmark.py --smoke` (ALL CHECKS PASSED).
   - `python3 scripts/check_multignss_benchmark.py --smoke` (ALL CHECKS PASSED).
   - Verify `eval_odaiba_ins` metrics ($p_{50} \le 1.80$ m, $\text{RMS} \le 3.50$ m, 0 false fixes).
