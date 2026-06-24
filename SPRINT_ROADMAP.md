# Sprint Roadmap — Industry-Leading PPP Accuracy

**Updated: 2026-06-24**  
**Status: Phase A1 feature-flag created, IONEX implemented, AR working, dataset acquired**

---

## Current Competitive Position

| Mode | Gneiss | RTKLIB | Status |
|------|--------|--------|--------|
| SPP (urban) | 1.8–2.1m | 2.8m | ✅ We win |
| RTK (urban) | 1.3–1.5m | 2.2–2.9m | ✅ We win |
| PPP (urban) | 5.0–5.3m | 2.0–4.0m | ❌ Lose by 1–3m |
| PPP (IGS station) | 0.95m | ? | 🟡 No precise truth |

### The Fundamental Problem

The single-epoch IEKF has a **~5m architectural accuracy floor** with broadcast ephemeris. The POST_MORTEM proved this conclusively — 10 hypotheses tested, none broke through. The SPP position anchor is simultaneously load-bearing (removing causes divergence to 12.7m) and limiting (re-anchoring each epoch to SPP quality sets the ceiling).

**The fix:** Multi-epoch sliding-window factor graph (SPRINT_PLAN.md). Jointly optimize across N epochs so carrier phase contributes relative position at mm precision.

---

## What Was Accomplished (Last 3 Sessions)

| Achievement | Why It Matters |
|-------------|---------------|
| SP3/CLK fallback bugs (×5) | Eliminated 76.9km PPP error on f9p |
| ANTEX O(1) HashMap indexing | PPP went from >5min to 1.9s per 200 epochs |
| AR bugs (×4) root-caused & fixed | AR now finds candidates, WL LAMBDA converges, fixes correct |
| IONEX parser + TEC interpolation | Infrastructure for cm-level iono prior (replaces 3m Klobuchar) |
| ALIC 24h IGS dataset | Proper benchmark data (2,880 epochs GPS+GLO+GAL) |
| 7 IGS stations downloaded | ALIC, CEDU, HOB2, NKLG, PARK, PERT, YARR |
| Multi-epoch feature flag | `EngineMode::PppMultiEpoch` + CLI `--mode ppp-me` |

---

## Phase A: Break the 5m Floor (CRITICAL PATH)

### A1 — 2-Epoch Joint Optimization ← BLOCKED (2026-06-24)
**File:** `crates/gneiss-rtk/src/engine/ppp_multi_epoch.rs` (production code, but diverges)

~~State vector `[x_{k-1}, x_k, ambiguities]`. Factors: PR/CP/Doppler per epoch, dynamics constraint between epochs, SPP prior on current epoch. LM optimization. Reuses existing GNSS factor code.~~

**Status: Architecturally blocked.** Full-state 2-epoch optimization with per-epoch ambiguities diverges after ~60 epochs. Root cause: without shared ambiguities across epochs, the phase-derived inter-epoch constraint is too weak — the prior from predicted covariance cannot prevent drift. The Schur complement covariance fix and LM damping slow but don't stop the drift.

**What was tried:**
- Full-state implementation with LM iteration, Schur complement covariance, priors on all state elements
- Solver is stable for first ~60 epochs (0.2-1.8m corrections), then diverges exponentially
- 95.6% of epochs diverged in full run

**What's needed:** Shared ambiguities across epochs (Phase A3) must come BEFORE the full-state optimization. The current per-epoch ambiguity model loses the mm-level carrier phase constraint between epochs — which is the entire purpose of multi-epoch optimization.

**New plan:** Skip Phase A1 full-state. Instead:
1. Implement shared ambiguities (Phase A3)
2. Then return to 2-epoch optimization with shared ambiguities

**Meanwhile:** The IEKF baseline on ALIC (IGS station) achieves 2.2m Hz50, 4.1m Hz95 — already below the 5m Odaiba threshold. The Odaiba urban dataset should be benchmarked with current IEKF to see if it already passes.

### A2 — N-Epoch Sliding Window
**Depends on:** A1 succeeding

Window of 5–10 epochs. VecDeque management, Schur complement marginalization of oldest epoch, re-linearization.

- **Gate:** Odaiba Hz50 < 4.0m

### A3 — Shared Ambiguities Across Window
**Depends on:** A2

Move ambiguities from per-epoch to shared. One ambiguity per satellite per frequency for entire window. `a_k = a_{k-1}` unless cycle slip.

- **Gate:** Odaiba Hz50 < 3.5m (ties RTKLIB)

---

## Phase B: Complete IONEX Pipeline

**Status:** Parser done, interpolation done, CLI flag done. **NOT YET BENCHMARKED** (full run was too slow, still processing).

### Remaining work:
1. **Fix IONEX performance** — binary search + bilinear interp is ~6ms per sat-epoch. Cache temporal window between epochs, pre-compute IPP per satellite.
2. **Set IONEX iono prior variance** — currently hardcoded 9.0 (3m std for Klobuchar) in UDUC measurements. Should be 0.0025 (0.05m std for IONEX) when iono_model=Ionex.
3. **Auto-download IONEX** — script to fetch matching IONEX file for benchmark day from AIUB FTP.

- **Key file:** `crates/gneiss-core/src/atmosphere.rs` → `iono_ionex()`
- **Parser:** `crates/gneiss-parsers/src/ionex.rs`
- **IONEX file:** `datasets/igs/codg3350.19i` (downloaded from ftp.aiub.unibe.ch)

---

## Phase C: Decoupled-Clock AR

**Depends on:** Phase A1 + Phase B

Keep IF combination for positioning (eliminates iono, 0.95m accuracy). Run parallel UDUC ambiguity observer that tracks L1/L2/i1 states solely for AR. After WL/NL AR fixes, apply integer constraints to IF solution via MW relationship.

- **Gate:** IF positioning accuracy preserved while AR converges and holds
- **Why not first:** AR doesn't help until float solution < ~2m. Phase A gets us there.

---

## Phase D: Fix Remaining Bugs

### D1 — WTZR Divergence (3,940km)
WTZR (Wettzell, Germany) diverges catastrophically while SUTH and ALIC work. GPS-only also diverges, so it's not multi-constellation. Suspect: LEIAR25.R3 antenna PCO/PCV, or P2 vs C2 observation type mapping quirk.

### D2 — RINEX Epoch Flag `&`
ALIC and CEDU files use `&` as epoch flag (external clock indicator in Hatanaka format). After CRX2RNX decompression, this should be fixed, but verify.

### D3 — Precise IGS Coordinates
ALIC and other stations need precise IGS14 SINEX coordinates (not RINEX APPROX POS) for accurate benchmarking.

---

## Phase E: Competitive Benchmark Sweep

Run all modes (SPP, PPP IF, PPP UDUC+AR, PPP IONEX+AR, PPP Multi-Epoch) across all working IGS stations vs RTKLIB. Publish matrix.

---

## Key Files Modified This Session

| File | Change |
|------|--------|
| `crates/gneiss-rtk/src/engine/ppp.rs` | IONEX dispatch + enable_ar UDUC gate |
| `crates/gneiss-rtk/src/engine/ppp_iekf.rs` | WL covariance from MW, MW threshold 50, is_fixed guard |
| `crates/gneiss-rtk/src/engine/types.rs` | PppMultiEpoch variant, IonosphereModel enum |
| `crates/gneiss-rtk/src/engine/config.rs` | enable_ar, iono_model, enable_tropo_gradients |
| `crates/gneiss-rtk/src/engine/processor/mod.rs` | ionex_grid, ionex_maps, ppp_multi_epoch_opt fields |
| `crates/gneiss-rtk/src/engine/mod.rs` | ppp_multi_epoch module registration |
| `crates/gneiss-rtk/src/engine/ppp_multi_epoch.rs` | **NEW** — stub solver (needs nalgebra fix) |
| `crates/gneiss-core/src/atmosphere.rs` | iono_ionex() — TEC bilinear+temporal interpolation |
| `crates/gneiss-parsers/src/ionex.rs` | **NEW** — IONEX v1.0 parser |
| `crates/gneiss-parsers/src/lib.rs` | ionex module registration |
| `bin/gneiss-cli/src/main.rs` | --ionex flag, --mode ppp-me, PppMultiEpoch dispatch |

---

## Immediate Next Step

Two parallel tracks:

**Track 1 — Shared Ambiguities (unblock multi-epoch):**
```
1. Move ambiguities from per-epoch RtkState to a shared window structure
2. One ambiguity per satellite per frequency for entire window
3. a_k = a_{k-1} unless cycle slip
4. Re-test 2-epoch joint optimization with shared ambiguities
```

**Track 2 — IONEX Pipeline (Phase B):**
```
1. Profile IONEX TEC interpolation performance
2. Set IONEX iono prior variance to 0.0025 (0.05m std)
3. Benchmark IONEX + IEKF vs Klobuchar + IEKF on ALIC
```

---

## Dataset Inventory

| Station | File | Epochs | Systems | Status |
|---------|------|--------|---------|--------|
| SUTH | `suth3350.19o` | 2,880 | GPS+GLO | ⚠️ Only ~200 have ≥4 DF sats |
| WTZR | `wtzr3350.19o` | 2,880 | GPS+GLO+GAL+SBAS | ❌ Diverges to 3,940km |
| ALIC | `alic3350.19o` | 2,880 | GPS+GLO+GAL | ✅ Works (CRX2RNX decompressed) |
| CEDU | `cedu3350.19o` | 2,880 | GPS | ✅ Downloaded (CRX2RNX) |
| HOB2 | `hob23350.19o` | 2,880 | GPS | ✅ Downloaded (CRX2RNX) |
| NKLG | `nklg3350.19o` | 2,880 | Multi | ✅ Downloaded (CRX2RNX) |
| PARK | `park3350.19o` | 2,880 | GPS | ✅ Downloaded (CRX2RNX) |
| PERT | `pert3350.19o` | 2,880 | GPS | ✅ Downloaded (CRX2RNX) |
| YARR | `yarr3350.19o` | 2,880 | GPS | ✅ Downloaded (CRX2RNX) |

**Products for day 335, 2019:**
- SP3: `cod20820.sp3` (CODE)
- CLK: `gfz20820.clk` (GFZ)
- NAV: `brdc3350.19n` (broadcast)
- ANTEX: `igs14.atx`
- IONEX: `codg3350.19i` (CODE, from ftp.aiub.unibe.ch)

---

## AR Status (All 4 Bugs Fixed)

| Bug | Fix |
|-----|-----|
| `find_ar_candidates` filtered `!s.is_iono_free` | `enable_ar` forces UDUC mode in `get_obs_and_corrections` |
| WL LAMBDA covariance from EKF state P (277k cyc²) | Build Q_WL from MW statistics (0.18/N cyc²) |
| AR at epoch 12 with unconverged MW | MW threshold increased 10→50 samples |
| Re-fixing every epoch after successful AR | Added `is_fixed && !has_new_sats` guard |
