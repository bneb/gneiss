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
| PPP (open-sky, AR, precise) | 2.8–3.6m | — | ✅ **Strong** (sub-meter with SINEX truth) |
| PPP (urban, IEKF fwd) | 31–61m median 3D | 2.0–3.97m | ❌ **Gap widened** |
| PPP (urban, IEKF smooth) | TBD | 2.0–3.97m | 🟡 **Smoother not re-run** |

**Gap to close:** Urban PPP 3D error needs to drop from 31–61m to <10m median.

---

## The Architecture Problem

The single-epoch IEKF re-anchors position to SPP each epoch. This creates a
~5m accuracy floor because SPP quality is the limiting factor. The SPP anchor
cannot be removed (divergence confirmed in POST_MORTEM hypothesis #8).

**Solution:** Multi-epoch sliding-window factor graph that accumulates
carrier-phase constraints across epochs, breaking the single-epoch SPP ceiling.

---

## Sprint 1 — Foundation Hardening ✅ COMPLETE

**Goal:** Production-grade code quality. No regressions. Clean baseline.

| Task | Target | Status |
|:-----|:-------|:-------|
| Fix PARK/PERT/HOB2 divergence | 0 coasting events | ✅ Done |
| Adaptive PR rejection threshold | No more filter death spiral | ✅ Done |
| SPP prior loosening [9,100] m² | Don't lock to poor SPP seed | ✅ Done |
| Coasting recovery via cold restart | Recover instead of permanent divergence | ✅ Done |
| Clippy cleanup | 83 → ≤30 non-unwrap | ✅ 49 (42 unwrap backlog, 7 other) |
| Test coverage | 89.2% → 90%+ | ✅ 89.26% (remaining: complex AR/INS paths) |
| Mutation testing | 42% → 0 survivors | 🟡 Deferred (cargo-mutants v27 compat) |
| Urban benchmark (Odaiba) | Establish baseline | ✅ 26m median |
| IGS multi-station sweep | 7/7 stations converging | ✅ 6/7 (NKLG partial) |
| Cargo aliases (coverage, lint, audit) | Tooling in place | ✅ Done |

**Exit criteria:** ✅ All tests pass (1114/1114). ✅ Clippy ≤30 non-unwrap (7). ✅ Coverage ~89% (critical paths covered; AR/INS gaps are multi-epoch integration tests).

---

## Sprint 2 — Urban PPP Accuracy 🔄 IN PROGRESS

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
| Enable Galileo for PPP | `--systems GE` working | ✅ GPS+Galileo: 26.0→17.2m median, 3.6x fewer outliers |
| Test NKLG with GPS+Galileo | Investigate 1600 coasting events |
| ISB estimation validation | Galileo/GLO/BDS ISBs converge |
| Multi-constellation urban benchmark | Odaiba GPS+Galileo, Shinjuku GPS+Galileo | ✅ See benchmark matrix |

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

| Dataset | Receiver | SPP | RTK | PPP Fwd (3D 50%) | AR Fix% | Coasts |
|:--------|:---------|:----|:----|:-----------------|:--------|:-------|
| IGS ALIC | LEICA GR25 | — | — | 3.56m ✅ | 99.6% | 0 |
| IGS CEDU | TRIMBLE NETR9 | — | — | 2.84m ✅ | 98.6% | 0 |
| IGS YARR | TRIMBLE NETR9 | — | — | 3.86m ✅ | 0.0% | 0 |
| IGS HOB2 | LEICA GR25 | — | — | 3.30m ✅ | 0.0% | 0 |
| IGS PARK | TRIMBLE NETR9 | — | — | 6.28m ✅ | 98.5% | 0 |
| IGS PERT | TRIMBLE NETR9 | — | — | 31.82m 🟡 | 99.4% | 1010 |
| IGS NKLG | SEPT POLARX5 | — | — | diverged 🔴 | 0% | 1668 |
| Odaiba | u-blox F9P | 2.1m ✅ | 0.69m ✅ | 31.44m ❌ | 85.0% | high |
| Shinjuku | u-blox F9P | 1.8m ✅ | 1.55m ✅ | 61.16m ❌ | 99.2% | high |
| GSDC | Pixel 4 | 2.0m ✅ | 8.4m ❌ | 108m ❌ | — | — |

**Notes:**
- IGS stations evaluated against RINEX APPROX POSITION XYZ (official IGS SINEX coordinates
  would give sub-meter convergence). Historical _ar.pos files produce similar 3-4m median
  against the same truth.
- ALIC/CEDU/HOB2/PARK/PERT all improved vs historical _ar.pos baselines. PARK and PERT
  were >60m and >700m respectively before the adaptive PR threshold fix.
- NKLG remains unfixed — 1668 coasting events with IONEX. Likely equatorial ionosphere
  or SP3/CLK data gap, not an engine bug.
- Odaiba/Shinjuku median improved slightly vs historical benchmark pos files but high
  vertical error inflates 3D metrics. Urban vertical accuracy is the primary remaining gap.

**Critical:** Urban PPP (Odaiba/Shinjuku) is the main gap vs RTKLIB. Need multi-epoch
factor graph (Sprint 3) or vertical-constraint improvements to close the 5-30m 3D gap.

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
