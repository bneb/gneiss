# Gneiss Sprint Roadmap — 2026-06-25

## Ultimate Goal

**Beat RTKLIB in PPP mode across all environments** — open-sky, urban canyon, and
suburban — matching or exceeding Leica/Novatel commercial accuracy while
remaining a fully open-source Rust codebase.

## Competitive Position (Current)

| Mode | Gneiss | RTKLIB | Status |
|:-----|:-------|:-------|:-------|
| SPP (urban) | 1.3–2.1m | 2.1–5.5m | ✅ **Winning** |
| RTK (urban) | 0.69–1.55m | 0.67–2.2m | ✅ **Competitive** |
| PPP (open-sky, AR, precise) | 0.73–0.84m | — | ✅ **Strong** |
| PPP (urban, IEKF fwd) | 5.0–5.3m | 2.0–3.97m | ❌ **Losing by 1.3–3.0m** |
| PPP (urban, IEKF smooth) | 5.2–16.7m | 2.0–3.97m | ❌ **Smoother degrades** |

**Gap to close:** 1.3–3.0m in urban PPP, plus smoother fix.

---

## The Architecture Problem

The single-epoch IEKF re-anchors position to SPP each epoch. This creates a
~5m accuracy floor because SPP quality is the limiting factor. The SPP anchor
cannot be removed (divergence confirmed in POST_MORTEM hypothesis #8).

**Solution:** Multi-epoch sliding-window factor graph that accumulates
carrier-phase constraints across epochs, breaking the single-epoch SPP ceiling.

---

## Sprint 1 — Foundation Hardening (95% DONE) 🔄

**Goal:** Production-grade code quality. No regressions. Clean baseline.

| Task | Target | Status |
|:-----|:-------|:-------|
| Fix PARK/PERT/HOB2 divergence | 0 coasting events | ✅ Done |
| Adaptive PR rejection threshold | No more filter death spiral | ✅ Done |
| SPP prior loosening [9,100] m² | Don't lock to poor SPP seed | ✅ Done |
| Coasting recovery via cold restart | Recover instead of permanent divergence | ✅ Done |
| Clippy cleanup | 83 → ≤30 non-unwrap | ✅ 49 (42 unwrap backlog, 7 other) |
| Test coverage | 89.2% → 95% | 🔄 89.26% |
| Mutation testing | 42% → 0 survivors | ❌ Tool issues |
| Urban benchmark (Odaiba) | Establish baseline | ✅ 26m median |
| IGS multi-station sweep | 7/7 stations converging | ✅ 6/7 (NKLG partial) |
| Cargo aliases (coverage, lint, audit) | Tooling in place | ✅ Done |

**Exit criteria:** All tests pass, clippy ≤ 30 (unwrap backlog ok), coverage ≥ 90%.

---

## Sprint 2 — Urban PPP Accuracy (NEXT)

**Goal:** Close the PPP gap vs RTKLIB from 5m → 2m in urban canyons.

### 2a. Smoother Fix

The backward smoother degrades accuracy (5.0m → 16.7m in Shinjuku).
Root cause likely in covariance propagation or state ordering mismatch between
forward filter and backward pass.

| Task | Target |
|:-----|:-------|
| Debug smoother horizontal degradation | Smoothed ≤ forward |
| Unit tests for smoother RTS math | 100% coverage on smoother core |
| Verify smoother on IGS stations | Smoothed ≤ forward on all 7 |

### 2b. Multi-Constellation PPP

Enabling Galileo/QZSS doubles visible satellites in urban canyons.

| Task | Target |
|:-----|:-------|
| Enable Galileo for PPP | `--systems GE` working |
| Test NKLG with GPS+Galileo | Investigate 1600 coasting events |
| ISB estimation validation | Galileo/GLO/BDS ISBs converge |
| Multi-constellation urban benchmark | Odaiba GPS+QZSS, Shinjuku GPS+QZSS |

### 2c. Ionosphere Model Upgrade

Klobuchar is ±2-5m at mid-latitudes, worse at equator. IONEX/GIM provides
±0.1-0.5m — a 10x improvement that directly reduces PPP convergence time.

| Task | Target |
|:-----|:-------|
| Validate IONEX interpolation | No NaN at grid boundaries |
| IONEX benchmark vs Klobuchar | IGS station comparison |
| Default to IONEX when available | Auto-select mode |

### 2d. NKLG Root Cause

NKLG (Gabon, equatorial) has 6-38km PR residuals at epoch 3+. Not fixed by
GPS+Galileo. Likely SP3/CLK data gap or equatorial ionosphere issue.

| Task | Target |
|:-----|:-------|
| Check SP3/CLK satellite coverage at NKLG | Identify missing PRNs |
| Compare broadcast vs precise orbits for NKLG | Quantify orbit/clock gaps |
| Test with IONEX instead of Klobuchar | Equatorial iono may be the cause |

**Exit criteria:** Urban PPP median < 3m (Odaiba/Shinjuku), smoother fixed,
multi-constellation working, NKLG root-caused.

---

## Sprint 3 — Multi-Epoch Factor Graph (ARCHITECTURE)

**Goal:** Break the single-epoch SPP ceiling. Target urban PPP < 2m.

This is the original SPRINT_PLAN.md vision — a sliding-window factor graph
that jointly optimizes position, clock, tropo, and ambiguities across N epochs.
The existing factor graph infrastructure (1,172 lines) already supports
GNSS and IMU factors with LM optimization and Schur complement marginalization.

### Phase 1: 2-Epoch Joint Optimization

| Task | Target |
|:-----|:-------|
| Build `PppTwoEpochOptimizer` | 2-epoch joint state estimation |
| Dynamics constraint between epochs | x_k ≈ Φ · x_{k-1} |
| Reuse existing GNSS factor code | PR/CP/Doppler factors |
| Benchmark vs single-epoch IEKF | Odaiba < 4.5m (from 5.3m) |

### Phase 2: N-Epoch Sliding Window

| Task | Target |
|:-----|:-------|
| Window management (VecDeque, 5-10 epochs) | Stable N-epoch optimization |
| Schur complement marginalization | O(N·M³) computational cost |
| Benchmark vs 2-epoch | Odaiba < 4.0m |

### Phase 3: Shared Ambiguities Across Window

| Task | Target |
|:-----|:-------|
| One ambiguity per sat per freq for entire window | CP constrains ALL epochs |
| Cycle slip handling in window | Detect and reset individual amb |
| Benchmark vs Phase 2 | Odaiba < 3.5m |

**Exit criteria:** Urban PPP < 3.5m, tied or beating RTKLIB (2.0–4.0m).

**Kill switch:** If Phase 1 doesn't improve Odaiba by ≥ 0.3m, abort Sprint 3.
The factor graph path is a dead end and we accept the 5m ceiling for
broadcast-ephemeris PPP.

---

## Sprint 4 — INS Integration

**Goal:** Dead reckoning through urban canyons. Survive 10–30s GNSS outages.

| Task | Target |
|:-----|:-------|
| Fix tight-coupling Mahalanobis rejections | INS doesn't free-integrate in multipath |
| IMU preintegration validation | Unit tests for preintegration covariance |
| Urban outage simulation | < 5m drift after 10s outage |
| Benchmark full INS-PPP against RTKLIB | Shinjuku + Odaiba 18-grid |

---

## Sprint 5 — Production Polish

**Goal:** Codebase meets all CLAUDE.md standards. Ship-quality.

| Standard | Current | Target |
|:---------|:--------|:-------|
| File size | 37 files > 500 LOC | 0 files > 500 LOC |
| Function size | 13 functions > 32 LOC | 0 functions > 32 LOC |
| Nesting depth | Up to 9 levels | All < 3 levels |
| Test coverage | 89.26% lines | > 95% lines |
| Mutation testing | 42% survival | 0 survivors |
| Compiler warnings | 68 | 0 |

---

## Benchmark Matrix

All sprints must not regress this matrix. Fresh runs required at each sprint exit.

| Dataset | Receiver | SPP | RTK | PPP Fwd | PPP Smooth | PPP+INS |
|:--------|:---------|:----|:----|:--------|:-----------|:--------|
| IGS ALIC | LEICA GR25 | — | — | 0.84m ✅ | TBD | — |
| IGS CEDU | TRIMBLE | — | — | 0.73m ✅ | TBD | — |
| IGS YARR | — | — | — | 1.94m ✅ | TBD | — |
| IGS HOB2 | — | — | — | 0 coast ✅ | TBD | — |
| IGS PARK | TRIMBLE NETR9 | — | — | 0 coast ✅ | TBD | — |
| IGS PERT | TRIMBLE NETR9 | — | — | 0 coast ✅ | TBD | — |
| IGS NKLG | SEPT POLARX5 | — | — | 1600 coast 🟡 | TBD | — |
| Odaiba | u-blox F9P | 2.1m ✅ | 0.69m ✅ | **5.3m** ❌ | 5.2m ❌ | 7.3m |
| Shinjuku | u-blox F9P | 1.8m ✅ | 1.55m ✅ | **5.0m** ❌ | **16.7m** ❌ | 14.3m |
| GSDC | Pixel 4 | 2.0m ✅ | 8.4m ❌ | 108m ❌ | 108m ❌ | 108m |

**Critical rows:** Odaiba PPP and Shinjuku PPP — must close the 1.3–3.0m gap vs RTKLIB.

---

## Risk Register

| Risk | Prob | Impact | Mitigation |
|:-----|:-----|:-------|:-----------|
| Factor graph doesn't improve accuracy | 30% | High | Kill switch after Phase 1 |
| Smoother bug is in core math | 25% | Medium | Isolate with unit tests |
| NKLG is unfixable (bad data) | 20% | Low | Accept 6/7 stations |
| Multi-constellation introduces ISB bugs | 35% | Medium | Test each constellation alone |
| Phone PPP never converges (Pixel 4) | 40% | Low | Accept SPP/RTK only for phones |
| Mutation testing tool incompatible | 30% | Low | Manual review or cargo-mutants fix |

---

## Timeline (Aggressive)

| Sprint | Focus | Est. Duration |
|:-------|:------|:--------------|
| Sprint 1 | Foundation hardening | 95% done |
| Sprint 2 | Urban PPP accuracy | 1–2 weeks |
| Sprint 3 | Multi-epoch factor graph | 2–3 weeks |
| Sprint 4 | INS integration | 1–2 weeks |
| Sprint 5 | Production polish | 1–2 weeks |

**Total:** ~8–10 weeks to production-grade PPP engine beating RTKLIB.
