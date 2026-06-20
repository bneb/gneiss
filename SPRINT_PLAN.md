# Sprint Plan v4 — Multi-Epoch Factor Graph PPP

## 🔴 Red Team Review

### Why This Sprint Exists

Nine hypotheses eliminated. The bottleneck is architectural: a **single-epoch IEKF cannot accumulate carrier-phase information** across epochs because position is re-anchored to SPP each epoch. The SPP prior (Sprint B) helped the tail but couldn't break the median floor. The only remaining path to <4m is multi-epoch optimization.

### What We're Building

A sliding-window factor graph that jointly optimizes position, clock, tropo, and ambiguities across N epochs. Between-epoch dynamics constraints allow carrier phase to contribute relative position information at mm precision — no integer AR needed.

### The Existing Infrastructure (Reusable)

| Component | File | Lines | What it does |
|:----------|:-----|:------|:-------------|
| `Factor` trait | `mod.rs` | 7-25 | residual(), jacobian(), information() |
| `PriorFactor` | `mod.rs` | 28-46 | State prior via error-state formulation |
| `FactorGraphOptimizer` | `mod.rs` | 49-247 | LM optimization with Huber loss |
| GNSS factors | `gnss_factors.rs` | 615 | PR/CP/Doppler error-state factors for tightly-coupled GNSS-INS |
| IMU factors | `imu_factors.rs` | 307 | IMU preintegration between epochs |
| Schur complement | (in mod.rs) | — | For marginalization |

The IMU factors demonstrate exactly the pattern we need: a between-epoch constraint that relates state at time k to state at time k+1. We replace IMU with a simpler dynamics constraint.

### Three-Phase Plan

#### Phase 1: 2-Epoch Joint Optimization (1 session)

**Goal:** Prove the architecture works with minimal scope.

Build a `PppTwoEpochOptimizer`:
- Accumulate state for epochs k-1 and k
- State vector: [x_{k-1} (21), x_k (21), ambiguities (shared)]
- Factors for each epoch: PR, CP, Doppler (reuse existing GNSS factor code)
- Dynamics factor between epochs: x_k ≈ Phi * x_{k-1}
- SPP prior on current epoch position
- LM optimization → extract smoothed state at k

**Benchmark target:** Odaiba Hz 50th < 5.0m (from 5.3m)

**Risk:** May not improve if dynamics constraint too loose or too tight.
**Mitigation:** Sweep process noise (0.1, 1.0, 10.0 m²/s) to find optimal.

#### Phase 2: N-Epoch Sliding Window (1-2 sessions)

Extend to 5-10 epochs:
- Window management with VecDeque
- Schur complement marginalization of oldest epoch
- Re-linearize each new epoch within the window
- Computational cost: O(N³) dense → O(N·M³) with marginalization where M is state size

**Benchmark target:** Odaiba Hz 50th < 4.0m

**Risk:** Numerical instability from repeated marginalization.
**Mitigation:** Add diagonal regularization in Schur complement.

#### Phase 3: Shared Ambiguities Across Window (1 session)

Currently, ambiguities are per-epoch state elements. Move them to shared parameters:
- One ambiguity parameter per satellite per frequency for the entire window
- Ambiguity factor: a_k = a_{k-1} (constant unless cycle slip)
- This allows CP to constrain position across ALL epochs in the window, not just adjacent pairs

**Benchmark target:** Odaiba Hz 50th < 3.5m

### What Could Kill This Sprint

1. **The factor graph doesn't improve accuracy.** Like GMF and UDUC before it, the architecture change may not matter if the dominant error is code multipath (~2-5m), not single-epoch limitations. *Probability: 30%*

2. **Numerical issues.** The existing factor graph was tested on tightly-coupled GNSS-INS with IMU, not on GNSS-only PPP. The state may be poorly conditioned without IMU constraints. *Probability: 20%*

3. **Implementation bugs.** The factor graph is 1,172 lines of existing code. Integrating PPP-specific factors correctly requires deep understanding of the error-state formulation and the coordinate frames. *Probability: 40%*

4. **Run time explodes.** Even with marginalization, 10 epochs × 21 states = 210-dimensional LM optimization. Each iteration requires building normal equations from scratch. May be too slow for 6,200 epochs. *Probability: 25%*

### What NOT to Do

- ❌ Don't build a full iSAM2 incremental solver — overkill for a 10-epoch window
- ❌ Don't rewrite the GNSS factors — reuse the existing ones
- ❌ Don't touch the IEKF — it remains the per-epoch initial guess provider
- ❌ Don't optimize for speed yet — correctness first

### Success Criteria

| Phase | Odaiba Hz 50th | Shinjuku Hz 50th | Gate |
|:------|:---------------|:-----------------|:-----|
| Current (SPP prior) | 6.4m | 7.5m | — |
| Phase 1 (2-epoch) | < 5.0m | < 7.0m | Proceed to Phase 2 |
| Phase 2 (N-epoch) | < 4.0m | < 6.0m | Proceed to Phase 3 |
| Phase 3 (shared amb) | < 3.5m | < 5.0m | Tied or beating RTKLIB |

If Phase 1 doesn't improve Odaiba by at least 0.3m, abort the sprint — the factor graph path is a dead end and we accept the 5.3m ceiling for broadcast-ephemeris PPP.
