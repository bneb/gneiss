# Sprint Plan v8

## State (2026-06-29)

Native IEKF is default `--mode ppp`. RTKLIB port is `--mode ppp-rtklib` (regression only, no active development).

### Accuracy
| Mode | Dataset | Key Metric | Goal | Gap |
|------|---------|-----------|------|-----|
| PPP native IEKF | CEDU (1000ep) | **1.7cm p95** | 1.00m | ✅ |
| PPP native IEKF | CEDU (2880ep) | crashes at ep1500 | 1.00m | 1 bug |
| RTK | Odaiba | 1.15m p50, 4.97m p95 | 0.25m p95 | 20× |

### One bug remaining

Constellation rotation at epoch 1000-1500 causes the 21+N-element state to resize (satellites rise/set). The state resize produces ill-conditioned normal equations with the tight position prior (weight 10000 on position vs ~0.01 on ambiguities). The SVD solve produces large dx, the per-iteration clamp fires, and the solve returns `StateDisappeared`. Root cause: ambiguity initialization during resize uses the tight position variance for new ambiguities, which is correct for the first N sats but too tight when combined with existing converged states — the condition number of the information matrix spikes.

---

## Phase 1: Fix Constellation Rotation (0.5 session)

### Solution: regularize per-ambiguity, not globally

The SVD regularization (1e-4) is applied globally, discarding valid ambiguity information. Instead: add a per-ambiguity regularization term that scales with the ambiguity variance. New ambiguities get large regularization (they're unknown), converged ambiguities get small regularization.

**Implementation**: In `compute_iteration_dx` and `compute_final_covariance`, add `λI` to the normal equations where `λ` is proportional to `1/σ²_amb` for each ambiguity. This is Tikhonov regularization with a diagonal matrix instead of a scalar.

**Alternative**: Simpler — before the SVD solve, check the condition number of `htwh_damped`. If >1e8, add incremental regularization until condition number drops below threshold. This is Levenberg-Marquardt style adaptive damping.

**Verification**: Run all 4 IGS stations, 2880 epochs. Zero crashes. Target p95 <0.1m on all stations with known position.

---

## Phase 2: Factor Graph (1-2 sessions)

With the native IEKF stable across all epochs, wire the multi-epoch optimizer.

### 2.1: Fix PppTwoEpochOptimizer
Remove internal `PppIteratedEkf::solve()`. Accept pre-solved `RtkState`. The processor feeds native IEKF output to the optimizer each epoch.

### 2.2: Shared position for static
State vector: `[position(3), clock_0, tropo_0, ..., clock_{N-1}, tropo_{N-1}, ambiguities]`. One position shared across window.

### 2.3: Benchmark
Target +20% p95 over single-epoch native IEKF.

---

## Phase 3: Urban Canyon (1-2 sessions)

### 3.1: Multi-constellation
Galileo + QZSS for Tokyo datasets. Already in SP3/CLK files, 21-element state has ISB slots.

### 3.2: Kinematic smoother
Enable position smoother for automotive dynamics. SPP seed (σ=5m) → forward convergence → backward propagation of converged info.

### 3.3: RTK AR hardening
Port PPP cascade AR validation to RTK. Multi-base selection (already in commit 4efe53e).

Target: Odaiba PPP p50 <5m (from 7m). RTK p95 <3m (from 5m).

---

## Phase 4: INS Coupling (2-3 sessions)

Wire IMU preintegration factors from FGO module into native IEKF loop. NHC already implemented. The 21-element state was designed for this.

Target: UrbanNav PPP-INS p50 <2m in urban canyon (from 7m).

---

## Timeline

| Phase | Sessions | Key Metric |
|-------|----------|------------|
| 1: Constellation rotation | 0.5 | 4/4 IGS stable 2880ep, p95 <0.1m |
| 2: Factor graph | 1-2 | +20% p95 over single-epoch |
| 3: Urban canyon | 1-2 | PPP p50 <5m, RTK p95 <3m |
| 4: INS coupling | 2-3 | PPP-INS p50 <2m |

**Total: 5-8 sessions to production-ready PPP + RTK + INS.**
