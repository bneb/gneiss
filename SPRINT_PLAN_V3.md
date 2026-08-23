# Sprint Plan v11 — Merged Architecture

**Sources**: Dellaert/Teunissen/Humphreys/Molteno/Dampf/Bisnath conference + type-state/Forster/AR refinement
**Status**: Sprint 1 executing

---

## Architecture Decisions (Final)

1. **Pose is a single variable** (position + attitude quaternion). This is the right granularity — IMU preintegration connects Pose(tₖ₋₁) → Pose(tₖ), and GNSS factors observe Pose through their geometric models.

2. **VariableId is a newtype** (`struct VariableId(u64)`). Compact, hashable, zero-cost. Variables are created once and persist for the window duration.

3. **Type-state config**: `EngineMode` is an enum where each variant carries its own required config. `RtkInsConfig` embeds `ImuConfig` — you can't construct it without IMU parameters. Compile-time safety.

4. **Forster IMU preintegration**: Δp, Δv, Δq accumulated between GNSS epochs. Bias Jacobians maintained for first-order correction during LM iterations (no reintegration needed).

5. **Two-step AR**: Float LM solve → extract ambiguity marginal covariance → LAMBDA → inject `FixedAmbiguityPriorFactor` → partial LM re-optimize. The fixed prior has near-infinite information, snapping pose to the correct integer.

6. **Correction stack as middleware**: `CorrectionPass` trait. `PreciseOrbits`, `SaastamoinenTropo`, `IonexCorrection` are composable passes. RTK vs PPP differ only in which passes are in the stack.

7. **No SPP reset after epoch 1**. Variables carry forward. Covariance naturally shrinks. The only exception is detected divergence (cov trace > 10,000 m²).

8. **Schur complement marginalization** when epoch slides out of window. Old variables condensed into a `MarginalPriorFactor` on the remaining variables.

---

## Sprint Execution Order

| # | Sprint | Sessions | LOC |
|---|--------|----------|-----|
| 1 | Core graph + type-state definitions | 2 | +2000 |
| 2 | IMU preintegration (Forster) | 2 | +800 |
| 3 | Unified measurement pipeline | 2 | +1500 |
| 4 | AR-injected graph (LAMBDA) | 2 | +800 |
| 5 | Schur complement marginalization | 1 | +500 |
| 6 | Benchmark harness + CI | 1 | +1000 |
| 7 | Hardening + dead code removal | 1 | -8000 |
| **Total** | | **11** | **+6600, -8000** |

Net: 1400 fewer lines than current, dramatically more capable.

---

## Sprint 1: Core Graph & Type-State Definitions

### 1.1 Variable System

```rust
use nalgebra::{Vector3, Vector6, UnitQuaternion, DMatrix, DVector, Matrix3, Matrix6};
use std::collections::{HashMap, VecDeque, BTreeMap};

/// Unique identifier for a variable in the factor graph.
/// Created once per variable, persists for the window duration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VariableId(u64);

impl VariableId {
    pub fn new(id: u64) -> Self { Self(id) }
    pub fn as_u64(self) -> u64 { self.0 }
}

/// Dimension of a variable: 6 for Pose, 3 for Velocity, 6 for ImuBias, etc.
#[derive(Debug, Clone, Copy)]
pub enum VariableDim {
    Scalar,     // 1
    Vector3,    // 3
    Vector6,    // 6
    Matrix3x3,  // 9
}

impl VariableDim {
    pub fn size(&self) -> usize {
        match self {
            Self::Scalar => 1,
            Self::Vector3 => 3,
            Self::Vector6 => 6,
            Self::Matrix3x3 => 9,
        }
    }
}

/// A variable in the factor graph.  Each variable has a current estimate
/// (the linearization point) and a dimension.
#[derive(Debug, Clone)]
pub struct VariableNode {
    pub id: VariableId,
    pub dim: VariableDim,
    /// Current estimate at which Jacobians are evaluated.  Updated after
    /// each LM iteration or after marginalization.
    pub value: DVector<f64>,
}

/// The key variable types in the graph.  Ambiguities and IMU biases
/// persist across ALL epochs.  Pose, Velocity, and ClockBias are
/// per-epoch (but the graph connects them across epochs via factors).
#[derive(Debug, Clone)]
pub enum VariableKind {
    /// 6-DOF pose at a specific epoch: [x, y, z, qx, qy, qz] in ECEF.
    /// Quaternion scalar qw is recovered from unit-norm constraint.
    Pose { epoch: u32 },
    /// 3-DOF velocity at a specific epoch (ECEF, m/s).
    Velocity { epoch: u32 },
    /// IMU bias (6-DOF: accel + gyro).  One per session, shared across epochs.
    ImuBias,
    /// Receiver clock bias per constellation per epoch.
    /// 2-DOF: [bias_m, drift_m_s].
    ClockBias { epoch: u32, constellation: Constellation },
    /// Troposphere zenith wet delay (1-DOF, meters).  Per epoch.
    TropoZwd { epoch: u32 },
    /// GLONASS IFB slope (1-DOF, m/freq_num).  One per session.
    IfbGlonass,
    /// Carrier-phase ambiguity (1-DOF, cycles).  One per satellite-frequency.
    /// Persists across ALL epochs — this is the key difference from EKF.
    Ambiguity { satellite: u16, frequency: u8 },
}

impl VariableKind {
    pub fn dim(&self) -> VariableDim {
        match self {
            Self::Pose { .. } => VariableDim::Vector6,
            Self::Velocity { .. } => VariableDim::Vector3,
            Self::ImuBias => VariableDim::Vector6,
            Self::ClockBias { .. } => VariableDim::Vector3, // [bias, drift, isb]
            Self::TropoZwd { .. } => VariableDim::Scalar,
            Self::IfbGlonass => VariableDim::Scalar,
            Self::Ambiguity { .. } => VariableDim::Scalar,
        }
    }
}
```

### 1.2 Factor Trait

```rust
/// A factor is a cost term in the nonlinear least-squares problem:
///   E(x) = Σ ||r_i(x)||²_{W_i}
/// where r_i is the residual and W_i is the information matrix.
pub trait Factor: std::fmt::Debug {
    /// The variables this factor connects to, in order.
    fn variables(&self) -> &[VariableId];

    /// Residual r(x) = h(x) - z.  Dimension: measurement_dim × 1.
    fn residual(&self, values: &VariableValues) -> DVector<f64>;

    /// Jacobian J = ∂r/∂x at the current linearization point.
    /// Dimension: measurement_dim × total_dim.
    /// `total_dim` is the sum of dimensions of all connected variables.
    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64>;

    /// Information matrix (inverse measurement covariance).
    /// Dimension: measurement_dim × measurement_dim.
    fn information(&self) -> DMatrix<f64>;

    /// Huber/Cauchy loss threshold.  None = no robust loss.
    fn robust_threshold(&self) -> Option<f64> { None }
}

/// Efficient lookup of variable values during factor evaluation.
pub struct VariableValues {
    /// Maps VariableId → (start_index_in_state_vector, dimension, current_value)
    index: HashMap<VariableId, (usize, usize)>,
    state: DVector<f64>,
}

impl VariableValues {
    pub fn new(variables: &BTreeMap<VariableId, VariableNode>) -> Self {
        let mut index = HashMap::new();
        let mut offset = 0;
        for (id, node) in variables.iter() {
            let dim = node.dim.size();
            index.insert(*id, (offset, dim));
            offset += dim;
        }
        let state = DVector::zeros(offset);
        // Copy current values into the state vector
        for (id, node) in variables.iter() {
            let (start, dim) = index[id];
            state.rows_mut(start, dim).copy_from(&node.value);
        }
        Self { index, state }
    }

    pub fn get(&self, id: VariableId) -> Option<nalgebra::DVectorView<f64>> {
        let (start, dim) = self.index.get(&id)?;
        Some(self.state.rows(*start, *dim))
    }

    pub fn total_dim(&self) -> usize { self.state.len() }
}
```

### 1.3 Estimation Graph

```rust
pub struct EstimationGraph {
    /// All active variables, ordered by VariableId for deterministic iteration.
    pub variables: BTreeMap<VariableId, VariableNode>,
    /// Factors connecting variables.
    pub factors: Vec<Box<dyn Factor>>,
    /// Dense prior from marginalized epochs.  Applied as an additional
    /// quadratic cost term: 1/2 ||J_prior * x - r_prior||².
    pub marginal_prior: Option<MarginalPriorFactor>,
    /// Window metadata: which epochs are currently in the window.
    pub window: VecDeque<EpochMetadata>,
    /// Next VariableId to assign.
    next_id: u64,
}

#[derive(Debug, Clone)]
pub struct EpochMetadata {
    pub time: GpsTime,
    /// Which variables were created at this epoch (for marginalization).
    pub variables: Vec<VariableId>,
    /// Number of satellite observations in this epoch.
    pub n_satellites: usize,
}

#[derive(Debug, Clone)]
pub struct MarginalPriorFactor {
    /// Variables that the prior connects to (subset of active variables).
    pub variables: Vec<VariableId>,
    /// Information matrix (prior_dim × prior_dim).
    pub hessian: DMatrix<f64>,
    /// Information-weighted residual (prior_dim × 1).  The prior cost is:
    ///   1/2 * x^T H x - r^T x + const
    pub gradient: DVector<f64>,
}

impl EstimationGraph {
    pub fn new() -> Self {
        Self {
            variables: BTreeMap::new(),
            factors: Vec::new(),
            marginal_prior: None,
            window: VecDeque::new(),
            next_id: 0,
        }
    }

    pub fn add_variable(&mut self, kind: VariableKind) -> VariableId {
        let id = VariableId::new(self.next_id);
        self.next_id += 1;
        let dim = kind.dim().size();
        self.variables.insert(id, VariableNode {
            id,
            dim: kind.dim(),
            value: DVector::zeros(dim),
        });
        id
    }

    pub fn add_factor(&mut self, factor: Box<dyn Factor>) {
        self.factors.push(factor);
    }
}
```

### 1.4 Type-State Configuration

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "mode")]
pub enum EngineConfig {
    Spp(SppConfig),
    Ppp(PppConfig),
    Rtk(RtkConfig),
    RtkIns(RtkInsConfig),
    PppIns(PppInsConfig),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtkConfig {
    pub base_position: [f64; 3],
    pub base_source: BaseSource,
    #[serde(default = "default_max_base_age")]
    pub max_base_age_s: f64,
    #[serde(default)]
    pub ar: Option<ArConfig>,
    pub initial_position: Option<[f64; 3]>,
    #[serde(default)]
    pub window_size: usize,
    #[serde(default)]
    pub ionosphere: IonosphereConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtkInsConfig {
    #[serde(flatten)]
    pub rtk: RtkConfig,
    pub imu: ImuConfig, // REQUIRED — deserialization fails without it
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImuConfig {
    pub lever_arm: [f64; 3],
    #[serde(default)]
    pub mounting_angles: [f64; 3],
    #[serde(default)]
    pub enable_nhc: bool,
    #[serde(default)]
    pub nhc_lever_arm: [f64; 3],
    #[serde(default)]
    pub tuning: ImuTuning,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArConfig {
    #[serde(default = "default_lambda_ratio")]
    pub min_ratio: f64,
    #[serde(default = "default_min_subsets")]
    pub min_subsets: usize,
    #[serde(default = "default_ffrt_prob")]
    pub ffrt_prob: f64,
}

// ImuTuning, IonosphereConfig, BaseSource, etc. follow the same pattern.
// Every sub-config has sensible Defaults.
// The key invariant: if an RtkInsConfig deserializes successfully,
// it is guaranteed to be valid.  No runtime validation needed.
```

### 1.5 Solver

```rust
pub struct SlidingWindowSolver {
    pub graph: EstimationGraph,
    pub window_size: usize,
    pub max_lm_iterations: usize,
    pub convergence_tol: f64,
    current_epoch: u32,
}

impl SlidingWindowSolver {
    pub fn new(config: &EngineConfig) -> Self { /* ... */ }

    /// Add GNSS + IMU factors for a new epoch, solve the graph,
    /// marginalize the oldest epoch if window is full.
    pub fn process_epoch(
        &mut self,
        rover: &EpochObs,
        base: Option<&EpochObs>,
        imu: Option<&[ImuMeasurement]>,
    ) -> Result<Solution, EngineError> {
        // 1. Create variables for the new epoch
        // 2. Push IMU measurements through preintegrator
        // 3. Build GNSS factors via measurement pipeline
        // 4. Add IMU preintegration factor connecting epochs
        // 5. Run LM solver
        // 6. If window is full, Schur-complement marginalize oldest epoch
        // 7. Extract and return Solution
        todo!("Sprint 1: graph construction and solve")
    }

    /// Levenberg-Marquardt optimizer over the active variable set.
    fn lm_solve(&mut self) -> Result<DVector<f64>, SolveError> {
        let values = VariableValues::new(&self.graph.variables);
        let total_dim = values.total_dim();

        for iteration in 0..self.max_lm_iterations {
            let (jtj, jtr, total_error) = self.build_normal_equations(&values);
            let lambda = self.lambda * jtj.diagonal().max();

            // Damped normal equations: (J^T W J + λI) Δx = -J^T W r
            let mut damped = jtj.clone();
            for i in 0..total_dim {
                damped[(i, i)] += lambda;
            }

            let delta = solve_linear_system(&damped, &(-jtr))?;
            self.apply_delta(&delta);
            let new_error = self.compute_total_error(&VariableValues::new(&self.graph.variables));

            if new_error < total_error {
                self.lambda *= 0.5; // Accept step, decrease damping
                if delta.norm() < self.convergence_tol { break; }
            } else {
                self.lambda *= 2.0; // Reject step, increase damping
            }
        }

        Ok(self.extract_state_vector())
    }
}
```

## Sprint 1 Success Criteria

After Sprint 1, the following must compile and pass tests:

1. `VariableId`, `VariableNode`, `VariableValues`, `VariableKind` — all defined, documented, tested
2. `Factor` trait — defined with residual, jacobian, information methods
3. `EstimationGraph` — can add variables and factors, iterate in order
4. `EngineConfig` with type-state variants — `RtkConfig`, `RtkInsConfig`, etc. — deserializes from JSON
5. `SlidingWindowSolver` — skeleton with process_epoch and lm_solve stubs
6. `VariableValues` — correct offset computation for heterogeneous variable dimensions
7. Unit tests: VariableValues offset correctness, EngineConfig serde round-trip, EstimationGraph add/remove
