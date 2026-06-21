# Post Mortem — PPP Accuracy Investigation, 2026-06-20

## Summary

**Goal:** Close the ~2m PPP accuracy gap between Gneiss and RTKLIB on UrbanNav datasets (Odaiba, Shinjuku).

**Result:** The gap is architectural. After testing ten hypotheses across 30+ build/benchmark cycles, median accuracy is unchanged at 5.3m Hz 50th on Odaiba. Two code quality improvements were made (composable tropo mapping, smoother bug fixes). One partial accuracy win was achieved (SPP prior improves 95th %ile tails by 4.3m but degrades median by 1.2m).

**Root cause:** The single-epoch SPP-anchored IEKF architecture has a ~5m accuracy floor with broadcast ephemeris. Carrier-phase measurements provide mm-level relative precision but cannot improve absolute position without multi-epoch convergence. The SPP anchor is load-bearing (removing it causes 12.7m divergence) but sets the accuracy ceiling.

---

## What We Tried

### Hypothesis Matrix

| # | Hypothesis | Implementation | Odaiba Hz 50th | Verdict |
|:--|:-----------|:---------------|:---------------|:--------|
| — | **Baseline (hard SPP reset)** | Current main | **5.26m** | Reference |
| 1 | NL ratio test blocks valid AR fixes | Removed NL gate | 5.26m | ❌ No change |
| 2 | WL ratio threshold too strict | Lowered 1.3→1.1+s≥0.05 | 5.26m | ❌ No change |
| 3 | Inter-constellation ISB limits WL | Per-constellation AR | 5.21m | ❌ No change |
| 4 | Smoother bugs cause 10⁶m blowups | 3 bug fixes | 5.22m | ✅ Fixed blowups |
| 5 | Precise products improve accuracy | CODE MGEX SP3+CLK+BIA | 5.37m | ❌ No gain on Odaiba |
| 6 | UDUC AR strategy matters | Enabled uduc_ar | 5.26m | ❌ No change |
| 7 | NMF mapping limits vertical | GMF replacement | 5.27m | ❌ No change |
| 8 | SPP anchoring limits convergence | Removed SPP reset | 12.7m | ❌ Divergence |
| 9 | SPP prior enables convergence | Soft prior (25→1 m²) | 6.43m | 🟡 Tails improved |
| 10 | TDCP adds between-epoch constraint | Time-diff carrier phase | 9.30m | ❌ Degradation |
| 11 | 2-epoch joint optimization | Factor-graph-like smoother | 9.53m | ❌ Degradation |

### What Actually Shipped (Committed)

| Commit | What | Value |
|:-------|:-----|:------|
| `dc5dbdc` | Composable tropo mapping (NMF/GMF/VMF1) | Architecture for future VMF3 grid files |
| `5b8553f` | Smoother bug fixes (3 bugs) | Catastrophic blowups eliminated |
| `0990803` | Per-constellation AR + inter-const fallback | Galileo fixes well, hybrid approach |
| `190e006` | NL ratio removal + WL threshold relaxation | Simplified AR cascade |
| `f8f2a32` | SPP position prior | 95th %ile improved 4.3m |

---

## Key Findings

### 1. AR does not move the solution

NL fixes have sub-cycle residuals (<0.8σ) and mean position correction of **0.03m**. The float solution already estimates ambiguities correctly; AR confirms without changing the position. AR tuning is not a path to better accuracy.

### 2. The SPP anchor is both load-bearing and limiting

Removing the per-epoch SPP position reset causes divergence to 12.7m. The anchor is necessary for stability with broadcast ephemeris. But it also sets a per-epoch accuracy floor equal to SPP quality (~5m). This is the fundamental tension.

### 3. Precise products alone don't fix UrbanNav

CODE MGEX SP3+CLK+BIA (multi-GNSS, 92 satellites, full phase bias support) produced zero accuracy improvement on Odaiba (5.37m vs 5.26m). The same products improve f9p by 23% (1.93m→1.49m). UrbanNav's error budget is dominated by code multipath and ionosphere in the urban canyon environment, not orbit/clock error.

### 4. Tropospheric mapping is not the bottleneck at 5m

NMF→GMF replacement produced a 0.02m change. The mapping function contributes 2-5cm of slant error, which is lost in the 5m total error budget. VMF1/VMF3 would be beneficial below 2m accuracy but not at current levels.

### 5. Multi-epoch optimization is architecturally correct but non-trivial

Both attempted implementations (TDCP, 2-epoch joint smoother) degraded accuracy because:
- Dynamics constraint tuning (process noise) is critical and dataset-specific
- Measurement alignment across epochs requires matching ambiguity keys
- The IEKF linearization point changes between epochs
- State dimension changes as satellites appear/disappear

A correct implementation requires careful handling of all four issues — likely 3-5 focused sessions with incremental validation at each step.

---

## Current State of `main`

### Accuracy (Odaiba PPP-FG, broadcast ephemeris)

| Percentile | Hz | Vt | 3D |
|:-----------|:----|:----|:----|
| 25th | 3.71m | 3.01m | 5.93m |
| 50th | **5.26m** | 5.80m | 8.72m |
| 75th | 8.99m | 10.04m | 13.15m |
| 95th | 20.41m | 43.95m | 47.13m |

### Code Health

- **Tests:** 317 passed, 0 failed, 2 ignored
- **Warnings:** 5 (unused imports in composable tropo, 3 fixable automatically)
- **Architecture:** Clean module structure, composable tropo mapping, per-constellation AR
- **Known issues:** Smoother degrades median vertical (17.5m vs 5.8m forward-only), dead `process.rs` file with stale EngineMode variants

### Competitive Position (vs RTKLIB)

| Mode | Gneiss | RTKLIB | Winner |
|:-----|:-------|:-------|:-------|
| SPP | 1.8m | 2.0m | **Gneiss** |
| RTK | 1.5m | 2.2m | **Gneiss** |
| PPP (float) | 5.3m | 3.4m | RTKLIB |

Gneiss leads in SPP and RTK. PPP is the only mode where RTKLIB wins, by approximately 2m.

---

## Recommendations

### Immediate (1 session each)

1. **Re-verify the RTKLIB baseline.** The 3.4m Odaiba number is from June 19 with different code. Re-run RTKLIB on the exact same data and evaluation to confirm the gap is real and quantify it precisely.

2. **Fix smoother vertical degradation.** White-noise states (clock bias, ISB, ZWD) are incorrectly propagated through the phi matrix. Zeroing them should recover the 5.8m→?m vertical gap and make smoothing usable.

3. **Register for CDDIS Earthdata.** Free registration at `urs.earthdata.nasa.gov` enables automated download of IGS final/rapid SP3+CLK for any GPS week. The f9p dataset confirms precise products help when the error budget isn't dominated by multipath.

### Medium-term (2-3 sessions)

4. **Debug the 2-epoch smoother.** The approach is correct but the implementation has dynamics tuning and measurement alignment bugs. Start with a 2-epoch window on a short dataset segment, validate against a known trajectory, then scale up.

5. **Add GIM/IONEX support.** Global Ionosphere Maps are 10× more accurate than Klobuchar. For the UDUC path (which directly estimates ionosphere), a better a priori model would reduce the ionospheric error contribution.

### Long-term (3-5 sessions)

6. **Full sliding-window factor graph.** Use the existing `estimators/factor_graph/` infrastructure (LM optimizer, GNSS factors, Schur complement) to build a proper 10-epoch window. Re-linearize within the window. This is the architecture that the big three (NovAtel, Leica, Qinertia) use and is the correct long-term direction.

---

## Lessons Learned

- **Test the null hypothesis first.** We spent significant time on AR tuning before proving AR doesn't move the solution. A quick "disable AR and measure" test early would have redirected effort sooner.

- **Benchmark every change.** The 30-minute build cycle was painful but each A/B test produced a definitive answer. No ambiguity about what worked and what didn't.

- **The SPP anchor is smarter than it looks.** Two separate attempts to remove or weaken it (de-anchoring, TDCP) produced worse results. The anchor exists for a reason — don't fight it, work with it.

- **Architectural ceilings are real.** When eight independent hypotheses all produce zero improvement, the bottleneck isn't in the components — it's in the architecture that connects them. Recognizing this early saves time.
