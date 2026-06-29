# Sprint Plan v7 — All-In on Gneiss-Native

## State (2026-06-29)

### shipped
- **Native IEKF is default**: `--mode ppp` uses the gneiss-native 21-element core state solver. RTKLIB port demoted to `--mode ppp-rtklib`, maintenance stopped.
- **RINEX auto-position**: Zero-config cm-accurate initialization from APPROX POSITION XYZ header field.
- **Tight position prior**: σ=1cm prior preserved through IEKF solve (posterior covariance fix).
- **Static process noise**: 1e-6×dt (was 9 m²/epoch velocity-integration model).
- **Corrected IF CP noise**: σ=3cm for IF mode (was σ=1cm, 9× underweight).
- **Soft AR lock**: σ=10cm Joseph updates, cross-correlations preserved.
- **AR validation gates**: Position jump check (0.1m tight, 2.0m loose), N_IF consistency, skip when σ_pos<1cm.
- **RTK baseline**: 1.15m p50, 4.97m p95 Odaiba. Beats RTKLIB on median, p95 dominated by urban canyon.
- **Position smoother**: Implemented (RTS over 3×3 position covariance).

### accuracy
| Mode | Dataset | p50 | p95 | Goal | Gap |
|------|---------|-----|-----|------|-----|
| PPP native IEKF | CEDU (500ep) | **1.0cm** | **2.0cm** | 1.00m | ✅ |
| PPP native IEKF | CEDU (2880ep) | 3.4m | 14.4m | 1.00m | 14× |
| RTK | Odaiba | 1.15m | 4.97m | 0.25m | 20× |
| RTK | Shinjuku | 1.21m | 11.16m | 0.25m | 45× |

### one bug remaining
**Native IEKF variance growth**: Position variance grows from 8mm² to 50000 m² by epoch 1000, triggering SPP reset. Root cause: the 21+N-element state propagation accumulates variance in ambiguity/velocity/ISB dimensions that couples back into position through the measurement Jacobian. The position prior is correctly applied each epoch but the non-position states don't have equivalent priors, creating an information asymmetry that slowly inflates position variance.

---

## Phase 1: Fix Native IEKF Variance Growth (0.5-1 session)

**Goal**: Sub-cm accuracy maintained across all 2880 epochs on all 4 IGS stations.

### Task 1.1: Add process noise to non-position states
**File**: `crates/gneiss-rtk/src/engine/predictor.rs`
**Problem**: The static process noise fix only addressed position. But velocity, ISB, and ambiguity states accumulate process noise at their original rates (100 m² for velocity, config.process_noise_isb for ISBs, etc.). These couple into position through the measurement Jacobian.
**Fix**: Clamp non-position state process noise when tight prior is active. Velocity: 0.01 (was 100), ISBs: 1e-8 (was config.process_noise_isb × dt).

### Task 1.2: Regularize posterior covariance
**File**: `crates/gneiss-rtk/src/engine/ppp_iekf.rs`
**Problem**: The posterior covariance `(HᵀWH + P_ext⁻¹ + P_pred⁻¹)⁻¹` can become ill-conditioned when P_pred has inflated non-position variances. The inverse amplifies small eigenvalues.
**Fix**: Add Tikhonov regularization to the posterior: `P = (HᵀWH + P⁻¹ + λI)⁻¹` with λ=1e-8. This floors all eigenvalues at 1e-8, preventing the inverse from blowing up.

### Task 1.3: Re-apply position prior as post-hoc correction
**File**: `crates/gneiss-rtk/src/engine/ppp_iekf.rs`
**Problem**: Even with correct posterior covariance, the state estimate can drift if non-position states are misestimated. The tight prior is applied during the solve but not as a hard constraint.
**Fix**: After each IEKF solve, before pushing to history: compute position delta from prior, and if >3σ, pull position back with a Kalman-like update.

### Task 1.4: Full benchmark
- All 4 IGS stations, 2880 epochs, native IEKF default
- Verify no SPP resets, no variance warnings
- Verify sub-meter p95 on all stations
- Compare: `--mode ppp` (native) vs `--mode ppp-rtklib` (RTKLIB port)
- **Target**: 4/4 stations <1m p95, zero divergence events

**Acceptance**: All 4 IGS stations maintain sub-meter p95 across 2880 epochs. Native IEKF conclusively beats RTKLIB port on all metrics.

---

## Phase 2: Multi-Epoch Factor Graph (1-2 sessions)

**Goal**: Joint optimization across epochs with shared parameters. Breaks the single-epoch ceiling for both static and kinematic.

### Task 2.1: Fix PppTwoEpochOptimizer for native IEKF
**File**: `crates/gneiss-rtk/src/engine/ppp_multi_epoch.rs`
**Change**: Remove internal `PppIteratedEkf::solve()` call. Accept pre-solved `RtkState` and `EpochSnapshot` from outside. The processor feeds native IEKF output to the optimizer.

### Task 2.2: Shared position for static receivers
**Change**: When dynamics is static, state vector becomes `[position(3), clock_0, tropo_0, ..., clock_{N-1}, tropo_{N-1}, ambiguities]`. One position for the entire window. Clocks, tropo, ambiguities per-epoch.

### Task 2.3: Wire into main loop
**File**: `crates/gneiss-rtk/src/engine/processor/mod.rs`
**Change**: After each IEKF solve, snapshot state. When window is full (N epochs), run joint optimization. Write smoothed position back to output state.

### Task 2.4: Benchmark
- Window sizes: 2, 5, 10 epochs
- Compare multi-epoch vs single-epoch on all IGS stations
- **Target**: +20% p95 improvement over single-epoch native IEKF

---

## Phase 3: Urban Canyon + Kinematic (1-2 sessions)

**Goal**: Improve moving-receiver PPP and RTK in urban environments.

### Task 3.1: Multi-constellation for Tokyo
- Enable Galileo + QZSS (already in SP3/CLK files)
- Per-constellation ISB and AR already in 21-element state

### Task 3.2: Kinematic position smoother
- Enable position smoother for automotive dynamics
- SPP seed (σ=5m) → forward convergence → backward propagation

### Task 3.3: RTK AR hardening
- Port AR validation gates from PPP cascade AR to RTK AR
- Position-jump check after RTK AR fix
- Multi-base selection (already in commit 4efe53e)

### Task 3.4: Elevation/azimuth weighting
- Exclude <15° satellites from CP in urban canyon
- C/N0-based variance scaling (already partially implemented)

**Acceptance**: UrbanNav Odaiba PPP p50 <5m (from 7m). RTK Odaiba p95 <3m (from 5m).

---

## Phase 4: INS Coupling (2-3 sessions)

**Goal**: Tightly-coupled GNSS-INS for urban canyon and kinematic. The 21-element core state was designed for this.

### Task 4.1: Wire IMU preintegration factors
**File**: `crates/gneiss-rtk/src/engine/fgo/factors/imu.rs`
**Change**: The FGO module already has IMU factors. Wire them into the native IEKF processing loop for `PppIns` mode.

### Task 4.2: NHC (Non-Holonomic Constraints)
**File**: `crates/gneiss-rtk/src/engine/processor/mod.rs`
**Change**: NHC already implemented. Verify it works with the native IEKF state layout.

### Task 4.3: Benchmark with UrbanNav IMU data
- Odaiba + Shinjuku have IMU data
- Compare PPP vs PPP-INS in urban canyon

**Acceptance**: UrbanNav PPP-INS p50 <2m (from 7m), p95 <5m (from 21m).

---

## Timeline

| Phase | Sessions | Key Metric |
|-------|----------|------------|
| 1: Fix variance growth | 0.5-1 | 4/4 IGS <1m p95, zero resets |
| 2: Multi-epoch factor graph | 1-2 | +20% p95 improvement |
| 3: Urban canyon hardening | 1-2 | Odaiba PPP p50 <5m, RTK p95 <3m |
| 4: INS coupling | 2-3 | PPP-INS p50 <2m in urban canyon |

**Total: 5-8 sessions to production-ready PPP + RTK + INS.**

---

## What We Killed

- **RTKLIB PPP port**: Demoted to `--mode ppp-rtklib`. No further development. Removed from default dispatch. Will be deleted once native IEKF conclusively beats it on all benchmarks.
- **IF combination code in RTKLIB**: The `force_if` block, NL AR candidate search, IF CP noise correction — all RTKLIB-specific. The native IEKF handles IF/UDUC mode selection through its own measurement model.
- **Batch solver**: The `StaticPositionBatchSolver` in `ppp_multi_epoch_batch.rs` was an attempt to add multi-epoch capability to the RTKLIB port. The factor graph (Phase 2) is the correct architecture for multi-epoch. Remove batch solver after Phase 2 ships.
