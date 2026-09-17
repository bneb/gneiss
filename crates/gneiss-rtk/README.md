# gneiss-rtk

`gneiss-rtk` provides the state estimation, factor graph optimization, and integer ambiguity resolution engines for high-precision GNSS and GNSS/INS positioning.

## Architecture

The crate is structured into four primary estimation and processing subsystems:

### 1. Sliding Window Factor Graph (SWFG)
- Formulates multi-epoch state estimation as a non-linear least squares factor graph.
- Jointly optimizes kinematic receiver position, velocity, receiver clock biases, atmospheric delays (tropospheric ZWD), and real-valued carrier phase ambiguities.
- Supports bidirectional forward/backward smoothing via covariance intersection.

### 2. 15-State Error-State Kalman Filter (ESKF / MEKF)
- Full 15-state error formulation: position error ($\delta\mathbf{p}^e$), velocity error ($\delta\mathbf{v}^e$), attitude error vector ($\delta\boldsymbol{\theta}$), accelerometer bias ($\delta\mathbf{b}_a$), and gyroscope bias ($\delta\mathbf{b}_g$).
- Closed-loop error quaternion feedback $\mathbf{q} \leftarrow \mathbf{q} \otimes \delta\mathbf{q}$ with error-state resetting.
- Vehicle motion constraints: Non-Holonomic Constraints (NHC) and Zero-Velocity Updates (ZUPT).
- Full Rauch-Tung-Striebel (RTS) backward smoothing across forward filter covariance history.

### 3. Ambiguity Resolution (LAMBDA & PPP-AR)
- Least-squares AMBiguity Decorrelation Adjustment (LAMBDA) for integer search.
- Partial Ambiguity Resolution (PAR) with dimension-dependent ratio tests ($R_{\text{req}} = f(k, P_{\text{fail}})$).
- Post-fix carrier phase residual screening and autonomous float fallback.
- Multi-frequency single-differenced PPP-AR using SINEX OSB satellite phase bias corrections and Melbourne-Wübbena wide-lane/narrow-lane cascades.

### 4. Network RTK / Virtual Reference Station (VRS)
- 2D Delaunay network triangulation across regional CORS reference stations.
- Multi-baseline double-difference network adjustment isolating tropospheric ZWD and ionospheric gradients.
- Localized virtual reference station synthesis with satellite transmit-time iteration and Sagnac rotation.

## Testing & Quality Standards

- Adheres to `AGENTS.md`: files $< 500$ LOC, functions $\le 32$ LOC, nesting depth $< 3$.
- Zero `unwrap()` calls in production code.
- Continuous regression gates enforced via `scripts/check_network_benchmark.py` and `scripts/check_multignss_benchmark.py`.
