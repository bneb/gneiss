# Sprint Plan v3 — Closing the PPP Accuracy Gap

## 🔴 Red Team Review

### The Gap We're Chasing

| Dataset | Gneiss PPP-FG | RTKLIB PPP | Gap | Notes |
|:--------|:--------------|:-----------|:----|:------|
| Odaiba | 5.3m Hz 50th | 3.4m Hz 50th | **~2m** | June 19 baseline, different code version |
| Shinjuku | 7.5m Hz 50th | 5.0m Hz 50th | **~2.5m** | Same caveat |

### Eight Hypotheses Tested and Eliminated

| # | Hypothesis | Test | Result |
|:--|:-----------|:-----|:-------|
| 1 | NL ratio test blocking valid fixes | Removed NL gate | NL fixes are perfect, don't move solution |
| 2 | WL ratio threshold too strict | Lowered 1.3→1.1 | +32% fix rate, zero accuracy change |
| 3 | Inter-constellation ISB destroys WL | Per-constellation AR | Works for 4+ sats/constellation, otherwise worse |
| 4 | Smoother bugs cause 10⁶m blowups | Fixed 3 bugs | Blowups eliminated, median unchanged |
| 5 | Precise products improve accuracy | CODE MGEX SP3+CLK | +20% on f9p, zero on Odaiba |
| 6 | UDUC AR strategy matters | Enabled uduc_ar | Zero impact without precise products |
| 7 | NMF mapping function is the bottleneck | Implemented GMF | 0.02m change — mapping not dominant error |
| 8 | SPP anchoring limits convergence | Removed per-epoch reset | Degradation to 12.7m — anchor is load-bearing |

### What These Results Prove

**The SPP anchor is a double-edged sword.** It prevents divergence but sets a per-epoch accuracy floor. Gneiss computes a single-epoch IEKF solution from an SPP seed — it cannot accumulate carrier-phase information across epochs because the position is reset each time. RTKLIB's EKF propagates position across epochs, accumulating information from every carrier-phase measurement ever observed. That's the architectural difference.

**The factor graph path is correct but underutilized.** The codebase has `estimators/factor_graph/` with IMU factors, GNSS factors, Schur complement, and LM optimization. The PPP IEKF (`PppIteratedEkf`) is labeled "Despite the historical 'fg' naming, this is an IEKF". The infrastructure for true multi-epoch optimization exists — it just hasn't been connected to PPP.

### Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|:-----|:-----------|:-------|:-----------|
| Path B doesn't converge | Medium | Wasted 2-3 sessions | Timebox to 1 session, fall back to Path C |
| Path C has hidden complexity | Medium | 5+ sessions | Start with 2-epoch window, scale up |
| Gap is in eval methodology | Low | Path A catches it | Do Path A first |
| Gap is in RTKLIB config, not code | Low | False premise | Path A verifies |
| Current code has regressed since baseline | Medium | We're chasing a moving target | Re-benchmark everything with current code |

### What NOT to Do

- ❌ Don't tune more AR thresholds — AR doesn't move the solution
- ❌ Don't chase better mapping functions — GMF proved mapping isn't dominant
- ❌ Don't try to remove SPP anchoring — it's necessary for stability
- ❌ Don't download more precise products — they don't help Odaiba

---

## 🎯 Three-Sprint Plan

### Sprint A: Verify the Gap (1 session)

**Goal:** Confirm the Gneiss-vs-RTKLIB gap exists with current code, same data, same eval.

**Tasks:**
- [ ] A.1 Re-run RTKLIB PPP on Odaiba dataset with broadcast ephemeris
- [ ] A.2 Re-run RTKLIB PPP on Shinjuku dataset with broadcast ephemeris
- [ ] A.3 Evaluate both with the same `gneiss-cli eval` tool
- [ ] A.4 Re-run Gneiss PPP-FG on both with current code
- [ ] A.5 Produce a single table: Gneiss vs RTKLIB on identical data/eval
- [ ] A.6 Update COMPARISON.md with verified numbers

**Success:** Quantified gap (or discovered it doesn't exist with current code).

**Goal command:**
```
/goal "Verify the PPP accuracy gap between Gneiss and RTKLIB on identical data and evaluation. Re-run RTKLIB PPP (broadcast ephemeris, kinematic, forward-only) on Odaiba and Shinjuku datasets. Re-run Gneiss PPP-FG on both with current code. Evaluate both engines with the same gneiss-cli eval tool against ground truth. Produce a single comparison table. If the gap is < 1m, the problem is solved — update COMPARISON.md. If > 1m, quantify it precisely for Sprints B and C."
```

### Sprint B: SPP Prior Instead of SPP Reset (2-3 sessions)

**Goal:** Replace the hard SPP position reset with a loose SPP prior, enabling multi-epoch carrier-phase convergence.

**Architecture:**
- Instead of `state.position = spp.position` (hard reset)
- Add SPP position as a measurement in the IEKF with variance ~25 m² (5m std)
- The dynamics model predicts position across epochs
- The IEKF refines with both code and carrier phase measurements
- The SPP prior keeps the solution loosely anchored without resetting

**Tasks:**
- [ ] B.1 Add `spp_prior_variance` field to EngineConfig (default 25.0)
- [ ] B.2 In `process_ppp`, instead of `state.position = spp.position`, push SPP position as a prior measurement in the IEKF measurement vector
- [ ] B.3 The SPP prior contribution: `h_row` = identity for position, `res` = state_pos - spp_pos, `weight` = 1/spp_prior_variance
- [ ] B.4 Reduce SPP prior variance gradually across epochs (from 25 → 1 m²) as the IEKF converges
- [ ] B.5 Keep the SPP hard reset as fallback for epoch 0 and after detected divergence
- [ ] B.6 Benchmark Odaiba: target Hz 50th < 4m
- [ ] B.7 Benchmark Shinjuku: target Hz 50th < 6m

**Success:** Hz 50th improves by ≥1m on either dataset.

**Goal command:**
```
/goal "Replace hard SPP position reset with a loose SPP prior to enable multi-epoch convergence. In process_ppp(), instead of state.position = spp.position, inject the SPP position as a measurement in the IEKF with initial variance ~25 m² (5m std). The dynamics model predicts position across epochs; the IEKF refines with carrier phase while the SPP prior prevents divergence. Reduce prior variance gradually (25→1 m²) as the filter converges. Keep hard reset for epoch 0 and divergence recovery. Run Odaiba: target Hz 50th < 4m."
```

### Sprint C: True Sliding-Window Factor Graph (3-5 sessions)

**Goal:** Replace the single-epoch IEKF with a multi-epoch sliding-window factor graph that jointly optimizes position, clock, tropo, and ambiguities across a window of epochs.

**Architecture:**
- Window size: 10-30 epochs
- Factors: pseudorange, carrier phase, Doppler, iono constraint, SPP prior
- Variables: position/velocity/attitude per epoch, clock per epoch, ZWD shared, ambiguities shared
- Marginalization: Schur complement on oldest epoch when window slides
- Reuse existing `estimators/factor_graph/` infrastructure (LM optimizer, Schur complement, factor types)

**Tasks:**
- [ ] C.1 Design the factor graph variable structure (per-epoch states + shared parameters)
- [ ] C.2 Implement `PppFactorGraph` struct with window management
- [ ] C.3 Port PR/CP/Doppler factors from IEKF to factor graph form
- [ ] C.4 Implement SPP prior factor (from Sprint B)
- [ ] C.5 Implement Schur complement marginalization for sliding window
- [ ] C.6 Benchmark Odaiba with window=10: target Hz 50th < 3m
- [ ] C.7 Benchmark Shinjuku with window=10: target Hz 50th < 5m
- [ ] C.8 Sweep window sizes (5, 10, 20, 30) to find optimal

**Success:** Hz 50th ≤ RTKLIB on both datasets.

**Goal command:**
```
/goal "Implement sliding-window factor graph PPP. Replace the single-epoch IEKF with a multi-epoch factor graph using the existing estimators/factor_graph/ infrastructure. Window of 10 epochs with PR/CP/Doppler/iono/SPP-prior factors. Joint optimization via LM, Schur complement marginalization on oldest epoch. Reuse existing factor types where possible. Benchmark Odaiba window=10: target Hz 50th < 3m."
```

---

## 📊 Success Metrics

| Metric | Current | Sprint A | Sprint B | Sprint C |
|:-------|:--------|:---------|:---------|:---------|
| Verified gap Odaiba | ~2m (stale) | **Measured** | — | — |
| Odaiba Hz 50th | 5.3m | 5.3m | < 4m | < 3m |
| Shinjuku Hz 50th | 7.5m | 7.5m | < 6m | < 5m |
| Odaiba vs RTKLIB | −1.9m | **Verified** | −1m | tied |

## 🔄 Loop Commands

**Sprint A loop:**
```
/loop /goal "Verify the PPP accuracy gap between Gneiss and RTKLIB on identical data and evaluation. Re-run RTKLIB PPP (broadcast ephemeris, kinematic, forward-only) on Odaiba and Shinjuku datasets. Re-run Gneiss PPP-FG on both with current code. Evaluate both engines with the same gneiss-cli eval tool against ground truth. Produce a single comparison table. If the gap is < 1m, the problem is solved — update COMPARISON.md. If > 1m, quantify it precisely for Sprints B and C."
```

**Sprint B loop:**
```
/loop /goal "Replace hard SPP position reset with an SPP prior measurement in the IEKF. Instead of state.position = spp.position, inject SPP position as a measurement with initial variance 25 m², decreasing to 1 m² as filter converges. Keep hard reset for epoch 0. Run Odaiba: target Hz 50th < 4m."
```

**Sprint C loop:**
```
/loop /goal "Implement sliding-window factor graph PPP using existing factor_graph/ infrastructure. Window of 10 epochs with PR/CP/Doppler/iono/SPP-prior factors. Schur complement marginalization. Benchmark Odaiba window=10: target Hz 50th < 3m."
```
