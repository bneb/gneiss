# Sprint Plan v2 — 2026-06-20

## 🔴 Red Team Pass (Second Look)

### What We've Established

1. **AR is not the accuracy bottleneck.** NL fixes have sub-cycle residuals and move the position mean 0.03m. WL fixes gate NL but don't change the position either. The float solution IS the accuracy-determining step.

2. **The PPP-EKF path no longer exists.** Both `--mode ppp` and `--mode ppp-fg` dispatch to the same `process_ppp()` → `PppIteratedEkf::solve()` and produce identical results (verified to 4dp). The 15× gap vs RTKLIB was a historical artifact from the old EKF code.

3. **The smoother is partially fixed.** Catastrophic blowups (10⁶m) eliminated. But median vertical degrades (5.8m → 17.5m Odaiba). White-noise states (clock bias, ISB, ZWD) are incorrectly propagated through the phi matrix.

4. **f9p dataset has precise products on disk.** SP3/CLK/BIA files from COD, ESA, GRG, SHAO sit unused in `datasets/rtkexplorer/sample_1/f9p_ppp_1224/`. The f9p PPP-FG result (3.5m) might already use them, or might not — we need to measure the precise product impact.

5. **UrbanNav (Shinjuku/Odaiba) lacks precise products.** GPS Week 2032 (Dec 2018) archives are not on BKG's public server. CDDIS registration required.

### Fresh Assessment

| Gap | Current | Root Cause | Fix |
|:----|:--------|:-----------|:----|
| Float PPP accuracy | 5.2m Odaiba, 7.5m Shinjuku | Broadcast ephemeris (1-2m orbit + 1.5m clock error) | Precise SP3+CLK |
| Smoothing vertical | 17.5m Vt 50th | White-noise states propagated through phi | Zero clk/ISB/ZWD in smoother phi |
| WL pass rate (Shinjuku) | 5% | Inter-constellation ISB, few sats/constellation | Acceptable — NL doesn't change solution anyway |
| AR fix rate (Odaiba) | 69% WL, 100% NL | Galileo+GPS+QZSS pooling works OK | Acceptable |

### What Matters and What Doesn't

**Matters for accuracy:**
- Precise orbit/clock products → ~2-3m gain (broadcast → precise)
- Float solution convergence quality
- Smoother correctness (for post-processing applications)

**Doesn't matter for accuracy (confirmed by benchmarks):**
- WL ratio threshold tuning
- NL ratio test presence/absence
- Per-constellation vs inter-constellation AR
- AR fix rate (NL barely moves the solution)

---

## 🎯 Revised Sprint Plan

### Sprint 1: Measure Precise Product Impact (f9p) 🎯

**Dataset:** `datasets/rtkexplorer/sample_1/f9p_ppp_1224/` (GPS Week 2137, Dec 2020)
**Has:** SP3, CLK, BIA files from COD/ESA/GRG/SHAO

**Tasks:**
- [ ] 1.1 Run f9p PPP-FG WITHOUT SP3/CLK (broadcast-only baseline)
- [ ] 1.2 Run f9p PPP-FG WITH SP3+CLK+BIA (precise products)
- [ ] 1.3 Compare accuracy → quantify precise product gain
- [ ] 1.4 If gain >2m, prioritize CDDIS registration for UrbanNav
- [ ] 1.5 Update COMPARISON.md with fresh f9p benchmarks

**Goal:**
```
/goal "Measure the accuracy impact of precise products using the f9p_ppp dataset. Run f9p PPP-FG twice: (1) broadcast-only baseline, (2) with --sp3 --clk --bia using the COD products already in datasets/rtkexplorer/sample_1/f9p_ppp_1224/. Compare Hz 50th and Vt 50th. Quantify the precise product gain to decide whether CDDIS registration is worth pursuing for UrbanNav. Update COMPARISON.md."
```

### Sprint 2: Fix Smoother Vertical Degradation

**Problem:** The RTS smoother propagates white-noise states (clock bias, ISB, ZWD) through the phi matrix. These should be zeroed in the state transition — they don't follow dynamics.

**Tasks:**
- [ ] 2.1 In `smooth_epoch`, set phi elements for white-noise states to 0
- [ ] 2.2 Run Odaiba PPP-FG with smoothing → verify Vt 50th < 8m
- [ ] 2.3 If successful, run Shinjuku with smoothing
- [ ] 2.4 Update COMPARISON.md

**Goal:**
```
/goal "Fix smoother vertical degradation. In smooth_epoch(), zero out the phi matrix rows/cols for white-noise states: rcv_clk_bias (idx 15), rcv_clk_drift (idx 19), isb_glo/gal/bds (idx 16-18), and zwd (idx 20). These states don't follow dynamics — propagating them through phi introduces bias. Run Odaiba PPP-FG with --enable-backward-smoothing, verify Vt 50th drops from 17.5m to < 8m."
```

### Sprint 3: Code Quality & Documentation

**Tasks:**
- [ ] 3.1 Remove all [deprecated] and [stale] rows from COMPARISON.md
- [ ] 3.2 Add MW value outlier filter (reject sats with MW EMA > 4σ from constellation median)
- [ ] 3.3 Rename `PppIteratedEkf` to `PppIekf` (or just fix the doc comment)
- [ ] 3.4 Fix the dead `process.rs` file (references non-existent `EngineMode::PppFg`)
- [ ] 3.5 Run `cargo clippy` and fix warnings

**Goal:**
```
/goal "Code quality cleanup: (1) remove all [deprecated] and [stale] rows from COMPARISON.md, (2) add MW value outlier filter to reject satellites with MW EMA > 4σ from constellation median in resolve_widelane_ar, (3) fix dead code in bin/gneiss-cli/src/process.rs (references non-existent EngineMode variants), (4) run cargo clippy and fix warnings."
```

---

## 📊 Success Metrics

| Metric | Current | Sprint 1 Target | Sprint 2 Target | Sprint 3 Target |
|:-------|:--------|:----------------|:----------------|:----------------|
| f9p PPP-FG (broadcast) | ? | Measure | — | — |
| f9p PPP-FG (precise) | 3.5m (stale) | ≤2.5m | — | — |
| Odaiba smoothed Vt 50th | 17.5m | — | <8m | — |
| Compiler warnings | 1 | — | — | 0 |

## 🔄 Loop Commands

**Sprint 1 loop:**
```
/loop /goal "Measure the accuracy impact of precise products using the f9p_ppp dataset. Run f9p PPP-FG twice: (1) broadcast-only baseline, (2) with --sp3 --clk --bia using the COD products already in datasets/rtkexplorer/sample_1/f9p_ppp_1224/. Compare Hz 50th and Vt 50th. Quantify the precise product gain to decide whether CDDIS registration is worth pursuing for UrbanNav. Update COMPARISON.md."
```

**Sprint 2 loop:**
```
/loop /goal "Fix smoother vertical degradation. In smooth_epoch(), zero out the phi matrix for white-noise states: rcv_clk_bias (15), isb_glo/gal/bds (16-18), rcv_clk_drift (19), zwd (20). Run Odaiba PPP-FG with --enable-backward-smoothing, verify Vt 50th drops from 17.5m to < 8m."
```

**Sprint 3 loop:**
```
/loop /goal "Code quality cleanup: remove stale COMPARISON.md rows, add MW outlier filter (4σ from constellation median), fix dead process.rs references, run clippy and fix warnings."
```
