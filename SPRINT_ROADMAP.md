# Sprint Roadmap — Industry-Leading PPP Accuracy

**Updated: 2026-06-24**  
**Status: Two critical bugs fixed. CP now works. AR achieves 0.84m Hz / 2.2m Up. File refactoring 25% done.**

---

## Current Competitive Position

| Mode | Gneiss | RTKLIB | Status |
|------|--------|--------|--------|
| SPP (urban) | 1.8–2.1m | 2.8m | ✅ We win |
| RTK (urban) | 1.3–1.5m | 2.2–2.9m | ✅ We win |
| PPP float (IGS) | 1.2m Hz / 7.5m Up | — | 🟡 Float only |
| **PPP AR (IGS)** | **0.84m Hz / 2.2m Up** | — | **🟢 Competitive** |
| PPP (urban) | not yet benchmarked | 2.0–4.0m | ⏳ Unknown |
| **PPP AR (3 IGS stations)** | **0.73–1.94m Hz / 2.2–9.2m Up** | — | **🟢 Working on 3/7** |

### What Changed This Session

Two critical bugs were root-caused and fixed:

1. **CP silently dropped (since inception).** `push_cp_measurement()` called `find_ambiguity_index()` (band-0 only), but `update_phase_ambiguities()` created UDUC bands 1/2/3 whenever raw L1/L2 existed. The lookup returned `None`, the `if-let` guard silently skipped, and **the IEKF was pseudorange-only**. Fix: check `is_iono_free` first, create band-0 before UDUC band intercept.

2. **NaN in covariance (cycle slip inflation).** Bug 25 inflates position/velocity covariance ×4 per slip. After many slips, values approach f64 overflow. `Phi*P*Phi^T + Q` produces NaN entries (~90 per event, 11 events per 2880-epoch run). Each caused a `StateDisappeared` error → 30m position jump. Fix: pre-clamp covariance at ±1e10 before the predict step.

---

## Accuracy Progression (ALIC, GPS-only PPP, 2880 epochs)

| Stage | Hz 50% | Hz RMS | Up RMS | 3D RMS | Errors | Converge |
|-------|--------|--------|--------|--------|--------|----------|
| PR-only (bug) | 1.61m | 1.83m | 5.36m | 5.66m | 43 | 1.03 ✗ |
| + CP fix | 1.45m | 1.60m | 5.17m | 5.41m | 7 | 0.86 |
| + NaN fix + pre-clamp | 1.22m | 2.48m | 7.46m | 7.86m | 0 | 0.84 |
| **+ AR enabled** | **0.84m** | **1.54m** | **2.18m** | **2.67m** | **0** | **0.28** |

AR converged scatter: East σ=1.09m (50%=0.49m), North σ=1.09m (50%=0.45m), Up σ=2.18m (50%=1.26m).

---

## Phase 1: Productionize AR (ACTIVE)

### 1.1 — AR Stability & Robustness
**File:** `crates/gneiss-rtk/src/engine/ppp_ar.rs`

AR has been tested on one station (ALIC). Need multi-station validation.

- [ ] Run AR on all 7 IGS stations (ALIC, CEDU, HOB2, NKLG, PARK, PERT, YARR)
- [ ] Fix WTZR divergence (3,940km — likely antenna or P2/C2 mapping)
- [ ] Add AR fix-count and ratio monitoring in CLI output
- [ ] **Gate:** AR fixes on ≥5 stations with Hz50 < 1.0m

### 1.2 — IONEX + AR
**Files:** `crates/gneiss-core/src/atmosphere.rs`, `crates/gneiss-parsers/src/ionex.rs`

IONEX provides 5cm iono prior (vs 3m Klobuchar) in UDUC mode. With AR forcing UDUC, IONEX becomes directly applicable.

- [ ] Benchmark IONEX + AR on ALIC (compare vs Klobuchar + AR)
- [ ] Fix IONEX interpolation performance (cache temporal window, pre-compute IPP)
- [ ] **Gate:** IONEX + AR Hz50 < 0.8m

### 1.3 — Multi-GNSS AR
- [ ] Test GPS+GLONASS AR (ALIC has 2 GLONASS sats)
- [ ] Test GPS+Galileo AR if data available
- [ ] **Gate:** Multi-GNSS AR Hz50 < 0.6m

---

## Phase 2: Code Quality Standards (CLAUDE.md compliance)

Current state vs targets:

| Standard | Current | Target | Gap |
|----------|---------|--------|-----|
| Compiler warnings | 0 | 0 | ✅ |
| Test failures | 0 | 0 | ✅ |
| Clippy warnings | 182 | 0 | 🔴 |
| Line coverage | 89.8% | >95% | 🟡 |
| Mutation survivors | ? | 0 | ⚪ |
| File size (ppp_iekf.rs) | 3,264 | <500 | 🔴 |
| File size (ppp.rs) | 2,931 | <500 | 🔴 |

### 2.1 — File Size Reduction (IN PROGRESS)

ppp_iekf.rs: 4353 → 3264 (-25%), ppp.rs: 3460 → 2931 (-15%)
Extracted: ppp_measurements.rs (580), ppp_ar.rs (535), ppp_antenna.rs (553)

- [ ] Move tests from ppp_iekf.rs to ppp_iekf_tests.rs (~1500 lines)
- [ ] Move tests from ppp.rs to ppp_tests.rs (~1200 lines)
- [ ] Extract ppp_spp_anchor.rs from ppp.rs (SPP prior logic, ~200 lines)
- [ ] **Gate:** All files <500 LOC except test files

### 2.2 — Clippy Cleanup
- [ ] Fix 182 clippy warnings (mostly `unwrap_used`, `too_many_arguments`)
- [ ] Add `#![deny(clippy::all)]` to lib.rs
- [ ] **Gate:** 0 clippy warnings

### 2.3 — Coverage
- [ ] Write ~30 targeted tests for uncovered branches (identified by workflow)
- [ ] Focus on ppp_iekf.rs (252 uncovered) and ppp.rs (193 uncovered)
- [ ] **Gate:** >95% line coverage

### 2.4 — Mutation Testing
- [ ] Run `cargo mutants` on high-coverage modules
- [ ] Kill all survivors or document equivalent mutants
- [ ] **Gate:** 0 mutation survivors

---

## Phase 3: Urban Benchmark

### 3.1 — Odaiba F9P Dataset
- [ ] Run SPP baseline on Odaiba (compare against 5.3m baseline)
- [ ] Run PPP float on Odaiba
- [ ] Run PPP AR on Odaiba
- [ ] **Gate:** PPP AR Hz50 < 2.0m on Odaiba (ties RTKLIB)

### 3.2 — Competitive Matrix
- [ ] Run all modes across all working stations
- [ ] Compare vs RTKLIB where ground truth available
- [ ] Publish accuracy matrix

---

## Deferred: Multi-Epoch Factor Graph

The original Phase A (2-epoch sliding window) was deprioritized after discovering:
1. The CP bug meant all prior testing was on PR-only — the "5m architectural floor" was actually a software bug
2. AR provides larger accuracy gains (71% vertical improvement) with less complexity
3. The full-state factor graph requires shared ambiguities to function, which is a significant refactor

**Revisit when:** AR is productionized and the accuracy limit of single-epoch AR is understood.

---

## Dataset Inventory

| Station | File | Epochs | Systems | Status |
|---------|------|--------|---------|--------|
| ALIC | `alic3350.19o` | 2,880 | GPS+GLO | ✅ AR working, 0.84m Hz50 |
| CEDU | `cedu3350.19o` | 2,880 | GPS | ✅ AR working, 0.73m Hz50 |
| YARR | `yarr3350.19o` | 2,880 | GPS | ✅ AR working, 1.94m Hz50 |
| HOB2 | `hob23350.19o` | 2,880 | GPS | 🟡 AR fixes 3×, 3.99m |
| SUTH | `suth3350.19o` | 2,880 | GPS+GLO | 🟡 AR 66% fix, 0.79m |
| NKLG | `nklg3350.19o` | 2,880 | Multi | ❌ Diverges (0% AR fix) |
| PARK | `park3350.19o` | 2,880 | GPS | ❌ Diverges (55m drift) |
| PERT | `pert3350.19o` | 2,880 | GPS | ❌ Diverges (167m drift) |
| WTZR | `wtzr3350.19o` | 2,880 | GPS+GLO+GAL | ❌ Not yet tested |

**Products for day 335, 2019:**
SP3 `cod20820.sp3`, CLK `gfz20820.clk`, NAV `brdc3350.19n`, ANTEX `igs14.atx`, IONEX `codg3350.19i`

---

## Known Bugs (All Fixed)

| Bug | Root Cause | Fix | Session |
|-----|-----------|-----|---------|
| IEKF was PR-only | `find_ambiguity_index` band-0 mismatch | Check `is_iono_free` first | 2026-06-24 |
| NaN in covariance | Cycle slip inflation → overflow | Pre-clamp at ±1e10 | 2026-06-24 |
| `auto_detect_dynamics` | Never implemented, defaulted to Automotive | Default to Static | 2026-06-24 |
| AR WL covariance | From EKF state P (277k cyc²) | Build Q_WL from MW statistics | Prior session |
| AR too early | Unconverged MW at epoch 12 | MW threshold 10→50 | Prior session |
| AR re-fixing | Fix every epoch after success | Add `is_fixed && !has_new_sats` | Prior session |
| ANTEX O(n) lookup | Linear scan per satellite | HashMap indexing | Prior session |
| SP3 clock fallback | Missing clock → satellite dropped | Fall back to broadcast clock | Prior session |
