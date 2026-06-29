# Sprint Plan v6 — Ship Gneiss-Native PPP, Hit 1m p95 on All Stations

## State at Start (2026-06-29, end of session 3)

### shipped
- **RINEX auto-position**: Every RINEX APPROX POSITION XYZ parsed automatically, seeds PPP with cm-accurate coordinates. Zero config.
- **RTKLIB port**: 0.67m CEDU, 0.98m HOB2 p95. Meets 1m PPP goal on 2/4 stations.
- **Static process noise fix**: Was 9 m²/epoch (velocity-integration dt³ model). Now 3e-5 m²/epoch (constant 1e-6×dt), matching RTKLIB.
- **IF CP noise correction**: IF-mode CP σ=3cm (was σ=1cm), matching actual IF combination noise amplification. Prevents Kalman gain from over-weighting CP 9×.
- **NL AR candidate search** (RTKLIB port): Searches ±N candidates within position uncertainty. With tight prior (σ=1cm), search radius is ±3 candidates minimum.
- **Soft AR lock** (both solvers): σ=10cm instead of 1mm hard-lock. Preserves cross-correlations.
- **Position RTS smoother**: Implemented, negligible with tight prior, useful for SPP-seeded case.
- **Cascade AR validation** (native IEKF): Position jump gate (0.1m tight, 2.0m loose), per-satellite N_IF consistency check, post-AR re-anchoring.

### known bugs blocking native IEKF
- **IEKF posterior covariance inflation**: Position prior (σ=1cm, info=10000) is applied as normal-equation damping, but posterior covariance grows to 0.27 m² (σ=0.52m) after solve. Root cause: `solve_inner` in `ppp_iekf.rs` computes `P_new = (H^T W H + P_prior^{-1})^{-1}` but the prior information may be overwritten or the covariance computation has a bug. Needs line-by-line trace.
- **State history double-push**: `process_ppp` and main loop both push to `state_history`, corrupting first epoch (180km error) and producing wrong epoch counts.
- **Cascade AR still makes wrong fixes**: AR validation catches most but a single wrong fix at σ=0.5m position variance slips through (jump=0.24m < 2.0m threshold) and cascades.

### what we learned
1. The tight position prior is THE unlock. Everything else (AR, noise models, process noise) is secondary.
2. The RTKLIB port preserves the prior because its predict step adds only 1e-6 to position diagonal. The native IEKF's predict added 9.0 m² via the velocity-integration model. Fixed.
3. The native IEKF's cascade AR is more sophisticated than RTKLIB's NL AR (LAMBDA vs geometry-only), but lacks position-constrained search. LAMBDA searches the full integer space regardless of position certainty.
4. IF mode + NL AR works on stations with good antenna calibrations (SEPT/Trimble). Leica antennas (ALIC) may need antenna-specific PCO/PCV corrections.
5. The smoother and factor graph are designed for the 21-element core state — they CANNOT work with RTKLIB layout. Must ship native IEKF to unlock them.

---

## Phase 1: Fix Native IEKF Posterior Covariance (0.5 session)

**Goal**: Native IEKF preserves σ=1cm position prior across all epochs. AR validation uses correct tight threshold (0.1m max jump). No divergence.

### Task 1.1: Trace posterior covariance computation
**File**: `crates/gneiss-rtk/src/engine/ppp_iekf.rs`, `solve_inner` method (~line 170-210)
**Problem**: Position prior is added to normal equations as `htwh[(i,i)] += 1.0/var`, but the posterior covariance `P_new = (H^T W H + P_prior^{-1})^{-1}` may be computed from a different H^T W H that doesn't include the prior term. Or the prior information is overwritten by a subsequent computation.

**Debug approach**:
- Add log: prior variance, H^T W H diagonal, posterior position variance
- Compare first epoch (cold start, no prior) vs second epoch (prior active)
- Expected: posterior σ_pos ≈ 1/sqrt(10000 + 9×0.25) ≈ 1cm

### Task 1.2: Verify predict preserves prior
**File**: `crates/gneiss-rtk/src/engine/predictor.rs`
**Already fixed**: Static process noise changed from q_acc×dt³/3 to 1e-6×dt.
**Verify**: With 30s epochs, q_pos = 3e-5 m². After 100 epochs, accumulated σ_pos = sqrt(0.0001 + 100×3e-5) = sqrt(0.0031) = 5.6cm. Still below 0.01 threshold for tight AR validation until ~30 epochs. Acceptable — AR doesn't run until epoch 10+.

### Task 1.3: Fix state history double-push
**File**: `crates/gneiss-rtk/src/engine/ppp.rs:247`, `processor/mod.rs:427`
**Fix**: Remove `engine.state_history.push(final_state)` from `process_ppp`. The main loop handles history. Or gate on a flag.

### Task 1.4: Validate on all 4 IGS stations
- Run 2880 epochs on CEDU, ALIC, YARR, HOB2
- Verify no divergence, position variance stays <0.01
- Compare p95 vs RTKLIB port baseline
- **Target**: All stations <1m p95

**Acceptance**: 4/4 IGS stations <1m p95 with native IEKF. No epoch has position variance >0.01 after epoch 10. RTKLIB port deprecated.

---

## Phase 2: Position-Constrained Cascade AR (0.5 session)

**Goal**: Cascade AR uses the tight position prior to constrain LAMBDA search. Eliminates remaining wrong fixes.

### Task 2.1: Add position constraint to LAMBDA ambiguity covariance
**File**: `crates/gneiss-rtk/src/engine/ppp_ar.rs`, `resolve_narrowlane_ar`
**Problem**: LAMBDA searches the full integer space defined by the ambiguity covariance. With tight position prior, the ambiguity covariance should be inflated by the position uncertainty projected into ambiguity space.

**Fix**: Before NL LAMBDA, compute the position→ambiguity projection:
```
P_amb_constrained = P_amb + D * P_pos * D^T
```
where D is the position→ambiguity Jacobian (rows = sats, cols = 3). This inflates the ambiguity covariance to reflect position uncertainty. LAMBDA then only considers integers consistent with the known position.

If position is known to 1cm and NL wavelength is 10.7cm, position uncertainty contributes <0.1 cycles to NL — LAMBDA should find a single unambiguous integer set.

### Task 2.2: Tighten AR validation with known position
- When σ_pos < 0.1m: max_jump = 0.05m (was 0.1m)
- N_IF consistency: reject if ANY sat fails (was 50% threshold)
- Add WL consistency: reject if WL integer changes from MW EMA by >0.5 cycles

### Task 2.3: Fallback to WL-only AR when NL fails
**Problem**: If NL AR is rejected, the current code returns Err and keeps float solution. Float solution drifts.
**Fix**: When NL fails but WL succeeds, apply WL-only constraint. WL fixes reduce the NL search space from ±140 candidates to ±20. The float solution with WL constraints is more stable than pure float.

**Acceptance**: Zero AR-induced divergence events across all 4 stations, 2880 epochs. AR fix rate >50%, fix rejection rate <10%.

---

## Phase 3: Ship Native IEKF as Default (0.5 session)

**Goal**: Native IEKF replaces RTKLIB port as the default PPP solver.

### Task 3.1: Route EngineMode::Ppp to native IEKF
**File**: `crates/gneiss-rtk/src/engine/processor/mod.rs`
**Change**: `EngineMode::Ppp` dispatches to `ppp::process_ppp` (native IEKF). Keep `EngineMode::PppRtklib` for regression testing.

### Task 3.2: Remove dead code
- Delete `force_if` IF combination code from RTKLIB port (no longer needed)
- Remove RTKLIB-specific AR code paths that are duplicated in native IEKF
- Clean up debug logging

### Task 3.3: Enable backward smoother by default for static PPP
**File**: `bin/gneiss-cli/src/process.rs`
**Change**: When mode is PPP and dynamics is static, enable backward smoothing automatically. The position smoother propagates well-converged epoch info backward through the position history. With tight prior, the effect is small but positive.

### Task 3.4: Full regression test
- All 4 IGS stations, 2880 epochs
- UrbanNav Odaiba/Shinjuku PPP kinematic
- Verify RTKLIB port produces identical results to before (no regression)
- Run full test suite (`cargo test --workspace`)

**Acceptance**: Native IEKF is default. All existing tests pass. IGS accuracy unchanged or improved. No regressions.

---

## Phase 4: Multi-Epoch Factor Graph (1-2 sessions)

**Goal**: Joint estimation of position across multiple epochs with shared ambiguities. Breaks the single-epoch accuracy ceiling for kinematic PPP.

### Task 4.1: Wire PppTwoEpochOptimizer into native IEKF loop
**File**: `crates/gneiss-rtk/src/engine/ppp_multi_epoch.rs`
**Change**: Remove internal `iekf.solve()` call. Accept pre-solved `RtkState` from outside. The optimizer snapshots state after each native IEKF solve and periodically runs joint optimization.

### Task 4.2: Shared position parameters for static receivers
**Change**: When dynamics is static, replace per-epoch position states with a single shared position parameter. State vector: `[position(3), clock_0, tropo_0, ..., clock_{N-1}, tropo_{N-1}, ambiguities]`.

### Task 4.3: Benchmark
- Compare multi-epoch vs single-epoch on all IGS stations
- Tune window size and process noise
- **Target**: +20% p95 improvement over single-epoch native IEKF

**Acceptance**: Multi-epoch mode runs without errors. Accuracy improvement or neutral. No regressions.

---

## Phase 5: Urban Canyon PPP + RTK (1-2 sessions)

**Goal**: Improve kinematic PPP and RTK in urban canyons.

### Task 5.1: Multi-constellation for Tokyo datasets
- Enable Galileo + QZSS (already in SP3 files)
- Per-constellation ISB and AR already implemented

### Task 5.2: Kinematic position smoother
- Enable backward smoothing for automotive dynamics
- With SPP seed (σ=5m), the forward filter converges from 5m→2m over time
- Backward pass propagates converged information to early epochs

### Task 5.3: RTK AR hardening
- Port AR validation gates from PPP to RTK
- Multi-base selection (already in commit 4efe53e)

**Acceptance**: UrbanNav Odaiba PPP p50 <5m (from 7m). RTK Odaiba p95 <3m (from 5m).

---

## Timeline

| Phase | Sessions | Key Metric |
|-------|----------|------------|
| 1: Fix IEKF covariance | 0.5 | pos_var <0.01 after epoch 10 |
| 2: Constrained cascade AR | 0.5 | Zero AR divergence events |
| 3: Ship native IEKF | 0.5 | 4/4 IGS <1m p95, native default |
| 4: Multi-epoch factor graph | 1-2 | +20% p95 improvement |
| 5: Urban canyon hardening | 1-2 | Odaiba PPP p50 <5m, RTK p95 <3m |

**Total: 4-5 sessions to ship native IEKF + factor graph + urban improvements.**

---

## Risk Register

| Risk | Prob | Impact | Mitigation |
|------|------|--------|------------|
| IEKF covariance bug deeper than identified | 30% | Phase 1 delayed | Fall back: keep RTKLIB port, fix only its AR divergence |
| LAMBDA position constraint doesn't prevent wrong fixes | 20% | Phase 2 partial | Accept partial AR (fix only high-confidence sats) |
| Factor graph no improvement over single-epoch | 25% | Phase 4 wasted | Validate on 2-epoch first, abort if no gain |
| Leica antenna bias (ALIC) not fixable in software | 40% | ALIC stays >1m | Document as known limitation, require ANTEX update |
| RTK fundamentally limited by urban canyon multipath | 60% | Phase 5 partial | Accept 1-3m p95, document INS coupling as requirement for 25cm |
