# Sprint Plan v9

## State (2026-06-29, post-commit 4272587)

Single-epoch RTK architecture is at its fundamental limit. AR is fully functional with
87% fix rate and 40-sat full-constellation resolution. The solution is perfectly stable
(0.000m position jumps after first fix, p50-p95 spread = 1.4cm) but systematically biased
by 0.9m from wrong first-fix NL integers.

### Accuracy
| Mode | Dataset | Key Metric | Goal | Gap |
|------|---------|-----------|------|-----|
| RTK single-epoch | Odaiba 4km | **0.91m p50, 0.94m p95** | 0.25m p95 | **3.8×** |
| RTK single-epoch | F9P 8km | 1.94m p50, 3.83m p95 | 0.25m p95 | 15× |
| PPP native IEKF | CEDU (1000ep) | 1.7cm p95 | 1.00m | ✅ |
| PPP native IEKF | CEDU (2880ep) | crashes at ep1500 | 1.00m | 1 bug |

### Root cause of 0.9m RTK bias
Code multipath (σ≈1m) biases the float solution by ~1m. The first AR fix picks the
LAMBDA integer set closest to the biased float estimates. MW and NL validations both pass
because code multipath biases their EMAs identically. PR residuals can't distinguish
correct from wrong fixes (ratio only 1.18×). **The bias source contaminates all
single-epoch validation pathways equally.**

### Why multi-epoch code averaging breaks the deadlock
Averaging DD pseudorange over N epochs reduces code multipath noise by √N.
With N=100 epochs (10s at 10Hz), PR noise drops from σ=1m to σ=0.1m. This is
precise enough to resolve NL integers (±5.4cm half-cycle at 2σ confidence).

Critically, the averaged PR provides an INDEPENDENT validation pathway — unlike MW/NL
EMAs which are derived from the same biased code measurements, the position-constrained
PR average uses satellite geometry to break the ambiguity-code correlation.

---

## Phase 5: Multi-Epoch Code Averaging (1-2 sessions)

**Goal**: p95 horizontal ≤ 25cm on Odaiba 4km by resolving the first-fix NL integer bias.

### 5.1: Sliding-window PR accumulator (0.5 session)

Add a per-satellite-pair ring buffer to `RtkState` that accumulates DD pseudorange
over the last N epochs. At each epoch, push the current DD PR and pop the oldest.
Maintain the running mean and variance.

```
state.pr_dd_accum: HashMap<(SatId, SatId), RingBuffer<f64>>
state.pr_dd_mean: HashMap<(SatId, SatId), f64>
```

The PR DD for satellite pair (rov, ref) is already computed in the measurement model.
Store it after each successful EKF update.

### 5.2: Position-constrained PR validation (0.5 session)

At AR time, use the time-averaged PR to validate the LAMBDA NL integers through
geometry rather than through the code-minus-phase NL combination:

1. Compute the expected DD PR at the fixed position using satellite geometry
2. Compare against the time-averaged DD PR from the accumulator
3. If the residual exceeds the expected noise (σ/√N), reject the fix

The validation function:
```
fn validate_geometry_pr(
    fixed_state: &RtkState,
    pr_accum: &PrAccumulator,
    ephemerides: &[Ephemeris],
) -> bool {
    for each (sat, ref) pair with N >= 10 accumulated epochs:
        expected_pr_dd = geometric_dd(fixed_pos, sat_pos, ref_pos, base_pos)
        residual = expected_pr_dd - pr_accum.mean[(sat, ref)]
        if residual > 3.0 * pr_accum.std[(sat, ref)] / sqrt(N):
            return false
    true
}
```

**Expected outcome**: Wrong first-fix NL integers (off by 1-2 cycles = 10-21cm at NL)
produce PR residuals of 10-21cm, which exceed the 3σ threshold of 3×0.1m=0.3m with
N=100. Correct fixes produce residuals < 0.3m. The first wrong fix is rejected,
the float solution continues to converge, and a subsequent correct fix is accepted.

### 5.3: Shared-position batch estimator (optional, 0.5 session)

If 5.2 alone doesn't close the gap, implement a lightweight batch estimator that
jointly solves for position across the sliding window using accumulated PR:

```
minimize Σ ||PR_dd_observed(t) - PR_dd_predicted(pos, t)||²
```

This replaces the single-epoch float position with a multi-epoch smoothed position.
The improved float position reduces the LAMBDA search space, making the first fix
more likely correct even before the geometry validation in 5.2.

### 5.4: Benchmark (0.5 session)

Run full Odaiba + F9P datasets. Target p95 ≤ 0.25m horizontal on Odaiba 4km.
If achieved, test on highway dataset (UrbanNav Odaiba open-sky segments or GSDC
after fixing the 60m vertical issue).

---

## Phase 6: Highway + Suburban Validation (1 session)

With multi-epoch averaging delivering 25cm on Odaiba:

### 6.1: Fix GSDC 60m vertical
Root cause identified: smartphone L1-only, min_elevation_deg silently ignored.
Config aliases already fixed. Need to verify on full dataset.

### 6.2: Odaiba open-sky segments
Identify highway sections in UrbanNav Odaiba (Rainbow Bridge crossing, Bayshore Route).
Evaluate p95 on these segments separately from urban canyon.

### 6.3: RTKLIB comparison
Run RTKLIB on same datasets with equivalent settings. Target: Gneiss p95 ≤ RTKLIB p95.

---

## Phase 7: INS Coupling (2-3 sessions)

Wire IMU preintegration factors into the multi-epoch estimator. NHC already implemented.
The 21-element state was designed for this.

Target: RTK-INS p50 < 0.15m in suburban, < 0.5m in urban canyon.

---

## Timeline

| Phase | Sessions | Key Metric |
|-------|----------|------------|
| **5: Multi-epoch code averaging** | **1-2** | **RTK p95 ≤ 0.25m on Odaiba 4km** |
| 6: Highway + suburban validation | 1 | RTK p95 ≤ 0.25m on highway/suburban |
| 7: INS coupling | 2-3 | RTK-INS p50 < 0.15m suburban |

**Next session: Phase 5.1-5.2 — sliding-window PR accumulator + geometry-based validation.**
