# Sprint Plan — 2026-06-20

## 📊 Progress (2026-06-20 Session)

| Sprint | Status | Key Result |
|:-------|:-------|:-----------|
| 1: Fix Smoother | 🟡 Partial | Catastrophic blowups (10⁶m) eliminated. AR re-resolution removed, symmetry enforced, NaN guards added. Median vertical degraded (5.8→17.5m). |
| 2: Precise Products | 🔴 Blocked | IGS BKG doesn't archive 2018 (Week 2032). Needs CDDIS authentication. |
| 3: Per-Constellation AR | 🟢 Done | Implemented with inter-const fallback. Galileo fixes consistently (4 sats → 3 pairs). Accuracy unchanged — NL barely moves float solution. |
| 4: PPP-EKF Investigation | ⬜ Todo | 10-15× gap vs RTKLIB remains. |
| 5: Code Quality | ⬜ Todo | Quick wins: fix comments, rename types, MW filter. |

**Key finding:** AR is NOT the PPP accuracy bottleneck. NL fixes have sub-cycle residuals and move the position <0.05m. Float solution quality — limited by broadcast ephemeris (~2-5m orbit+clock error) — determines accuracy. The 2m gap vs RTKLIB PPP-FG is from the float solution, not AR.

## 🔴 Red Team Assessment (Original)

### What's Actually Broken

1. **PPP-EKF is 10-15× worse than RTKLIB.** Shinjuku 30m vs 2m, Odaiba 16m vs 4m. This is the single biggest performance gap in the entire benchmark matrix. PPP-FG (IEKF) is better at 5m / 3.4m but still ~2m behind RTKLIB's float PPP.

2. **Backward smoothing produces 10⁶m errors** on PPP-FG. The RTS smoother (`smoother.rs`) has at least two bugs: (a) it forces AR re-resolution on smoothed states which can pick different integers than the forward pass, creating state inconsistencies; (b) the smoother covariance update `p_k + C*(p_k1 - p_pred)*C^T` can produce negative eigenvalues when the smoothed covariance increases relative to the predicted (mathematically invalid for an RTS update).

3. **Inter-constellation WL AR has critically low ratio.** The WL LAMBDA search mixes GPS, Galileo, and QZSS satellites with uncalibrated ISB offsets of 10-60 cycles. WL ratio is 1.0-1.1 for most epochs, making LAMBDA unable to discriminate between integer sets. Shinjuku WL pass rate is only 5%.

4. **No precise products for UrbanNav datasets.** PPP falls back to broadcast ephemeris (~1-2m orbit error, ~1.5m clock error), directly limiting accuracy floor to ~3-5m. RTKLIB's 2-4m PPP results likely use IGS rapid/ultra-rapid SP3+CLK.

5. **GPS PRN 10 in Odaiba has corrupt MW EMA.** MW value of -115 cycles with only 60-70 samples (vs 2600+ for other sats) — this outlier corrupts the WL LAMBDA search for all pairs. The MW threshold of 10 samples is too low to reject it.

### What's Working

- **NL fix quality is excellent.** All residuals sub-cycle (<0.8σ), zero position jump rejections, position corrections mean 0.03m. The NL ratio test removal was correct.
- **SPP and RTK already beat RTKLIB.** Gneiss wins SPP (1.8m vs 2.2m) and RTK (1.3m vs 2.9m) on UrbanNav datasets.
- **Unit tests are solid.** 317 tests pass, 0 fail.
- **Float PPP solution converges.** Even without AR, the IEKF converges to 5m (Odaiba) / 7m (Shinjuku) from a cold SPP start of 191m.

### Wrong Assumptions to Discard

- ❌ "AR is the key to closing the PPP accuracy gap" → AR barely moves the position (0.03m). Float solution quality is the bottleneck.
- ❌ "Inter-constellation AR gives more pairs → better fixes" → ISB mixing destroys WL ratio. Per-constellation AR is cleaner.
- ❌ "Lower WL threshold increases fix rate without risk" → True for Odaiba (37%→69%), but Shinjuku stays at 5% — the real problem is ISB, not the threshold.
- ❌ "Backward smoothing improves accuracy" → Currently BROKEN for PPP-FG. Must fix before assuming it helps.

### Known Code Issues

| Issue | Location | Severity |
|:------|:---------|:---------|
| Position jump comment says 5m, code uses 20m | `ppp_iekf.rs:157,161` | Low |
| "PPP-FG" naming is misleading (it's IEKF, not Factor Graph) | `ppp_iekf.rs:14` comment | Low |
| MW minimum sample threshold (10) is too low for reliable WL | `ppp_iekf.rs:252-253` | Medium |
| Smoother re-fixes AR on smoothed states | `smoother.rs:57-64` | High |
| Smoother covariance update can go negative-definite | `smoother.rs:137` | High |

---

## 🎯 Sprint Plan

### Sprint 1: Fix Backward Smoothing (est. 1-2 sessions)

**Why first:** Recovers ~2m accuracy for free. Enables proper comparison with RTKLIB baseline. Unblocks all subsequent accuracy work.

**Tasks:**
- [ ] 1.1 Remove AR re-resolution in smoother (`smoother.rs:57-64`) — use forward-pass AR fixes as-is
- [ ] 1.2 Fix RTS covariance update to use Joseph form: `P_smooth = P_k - C*(P_pred - P_k1)*C^T`
- [ ] 1.3 Add covariance symmetry enforcement after smooth
- [ ] 1.4 Run Odaiba PPP-FG with smoothing → verify Hz 50th < 4m
- [ ] 1.5 Run Shinjuku PPP-FG with smoothing → verify Hz 50th < 6m
- [ ] 1.6 Update COMPARISON.md with fresh smoothed numbers

**Goal command:**
```
/goal "Fix backward smoothing for PPP-FG. Remove AR re-resolution in smoother.rs, fix RTS covariance update to use correct Joseph form P_smooth = P_k - C*(P_pred - P_k1)*C^T instead of P_k + C*(P_k1 - P_pred)*C^T, add symmetry enforcement. Run Odaiba PPP-FG with --enable-backward-smoothing and verify Hz 50th < 4m (baseline 3.4m). Update COMPARISON.md."
```

### Sprint 2: Download Precise Products (est. 1 session)

**Why:** Directly addresses the broadcast-ephemeris accuracy floor. IGS rapid SP3+CLK are free.

**Tasks:**
- [ ] 2.1 Determine GPS week/tow for Shinjuku and Odaiba datasets
- [ ] 2.2 Download IGS rapid SP3 (.sp3) and CLK (.clk) files from CDDIS/NASA
- [ ] 2.3 Download CODE or CNES phase bias (.bia) files for FCB-capable PPP-AR
- [ ] 2.4 Run both datasets with --sp3, --clk, --bia flags
- [ ] 2.5 Compare accuracy vs broadcast-only runs
- [ ] 2.6 Update COMPARISON.md

**Goal command:**
```
/goal "Download IGS rapid precise products (SP3+CLK) for UrbanNav Shinjuku and Odaiba datasets. Determine GPS week/TOW from nav files, fetch from CDDIS, run PPP-FG benchmarks with --sp3 --clk flags, compare accuracy vs broadcast-only baseline."
```

### Sprint 3: Per-Constellation WL AR (est. 1-2 sessions)

**Why:** Eliminates ISB mixing in WL LAMBDA. Should increase Shinjuku WL pass rate from 5% to >>50%.

**Tasks:**
- [ ] 3.1 Modify `build_ar_subset` to support per-constellation mode via config flag
- [ ] 3.2 Modify `resolve_cascade_ar` to loop over per-constellation groups
- [ ] 3.3 Add `ar_inter_constellation: bool` to EngineConfig (default true for back-compat)
- [ ] 3.4 Benchmark Odaiba + Shinjuku with inter=false → measure WL pass rate
- [ ] 3.5 If WL pass rate > 50%, make per-constellation the default
- [ ] 3.6 Update COMPARISON.md

**Goal command:**
```
/goal "Implement per-constellation WL AR as an alternative to inter-constellation. Add ar_inter_constellation flag to EngineConfig, modify build_ar_subset and resolve_cascade_ar to support per-constellation groups. Run Odaiba and Shinjuku benchmarks with per-constellation mode. If WL pass rate improves (target >50%), make it the default."
```

### Sprint 4: PPP-EKF Investigation (est. 1-2 sessions)

**Why:** 10-15× gap vs RTKLIB is unacceptable. The IEKF variant (PPP-FG) is much better, suggesting the EKF variant has a specific bug or missing feature.

**Tasks:**
- [ ] 4.1 Compare EKF vs IEKF code paths — identify what IEKF has that EKF doesn't
- [ ] 4.2 Check if EKF uses iono-free combination (IEKF does)
- [ ] 4.3 Check if EKF has proper outlier rejection (IEKF does)
- [ ] 4.4 Check EKF convergence diagnostics vs IEKF
- [ ] 4.5 Port the most impactful IEKF features to EKF
- [ ] 4.6 Benchmark both datasets EKF vs IEKF

**Goal command:**
```
/goal "Diagnose why PPP-EKF is 10-15× worse than RTKLIB (30m vs 2m Shinjuku). Compare EKF and IEKF code paths in ppp_iekf.rs vs the EKF estimator. Identify specific missing features in EKF: iono-free combination, outlier rejection, iteration, convergence criteria. Port critical features. Benchmark before/after."
```

### Sprint 5: Code Quality Cleanup (est. 0.5 session)

**Why:** Low effort, prevents future confusion.

**Tasks:**
- [ ] 5.1 Fix position jump comment: 5m → 20m (or change threshold to 5m)
- [ ] 5.2 Rename PppIteratedEkf → PppIekf (or update doc comment)
- [ ] 5.3 Increase MW minimum sample threshold from 10 to 30
- [ ] 5.4 Fix any remaining compiler warnings

**Goal command:**
```
/goal "Code quality cleanup: fix position jump threshold comment (5m→20m), rename PppIteratedEkf→PppIekf, increase MW minimum samples from 10→30, fix any warnings."
```

---

## 📊 Success Metrics

| Metric | Current | Sprint 1 Target | Sprint 2 Target | Sprint 3 Target | Final Target |
|:-------|:--------|:----------------|:----------------|:----------------|:-------------|
| Odaiba PPP-FG Hz 50th | 5.3m | ≤4.0m | ≤2.5m | ≤2.5m | ≤2.5m |
| Shinjuku PPP-FG Hz 50th | 7.5m | ≤6.0m | ≤3.5m | ≤3.5m | ≤3.5m |
| Odaiba PPP-EKF Hz 50th | 16.5m | — | — | — | ≤5.0m |
| Shinjuku PPP-EKF Hz 50th | 30.3m | — | — | — | ≤5.0m |
| WL pass rate (Odaiba) | 69% | — | — | ≥80% | ≥80% |
| WL pass rate (Shinjuku) | 5% | — | — | ≥50% | ≥50% |
| Tests passing | 317 | 317 | 317 | 317 | 317 |

## 🔄 Loop Commands

Run each sprint sequentially. After each sprint completes, commit results and update this document.

**Sprint 1 loop:**
```
/loop /goal "Fix backward smoothing for PPP-FG. Remove AR re-resolution in smoother.rs, fix RTS covariance update to use correct Joseph form P_smooth = P_k - C*(P_pred - P_k1)*C^T instead of P_k + C*(P_k1 - P_pred)*C^T, add symmetry enforcement. Run Odaiba PPP-FG with --enable-backward-smoothing and verify Hz 50th < 4m (baseline 3.4m). Update COMPARISON.md."
```

**Sprint 2 loop:**
```
/loop /goal "Download IGS rapid precise products (SP3+CLK) for UrbanNav Shinjuku and Odaiba datasets. Determine GPS week/TOW from nav files, fetch from CDDIS, run PPP-FG benchmarks with --sp3 --clk flags, compare accuracy vs broadcast-only baseline."
```

**Sprint 3 loop:**
```
/loop /goal "Implement per-constellation WL AR as an alternative to inter-constellation. Add ar_inter_constellation flag to EngineConfig, modify build_ar_subset and resolve_cascade_ar to support per-constellation groups. Run Odaiba and Shinjuku benchmarks with per-constellation mode. If WL pass rate improves (target >50%), make it the default."
```

**Sprint 5 loop (can run anytime):**
```
/loop /goal "Code quality cleanup: fix position jump threshold comment (5m→20m), rename PppIteratedEkf→PppIekf, increase MW minimum samples from 10→30, fix any warnings."
```
