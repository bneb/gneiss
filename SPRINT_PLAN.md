# Sprint Plan v5 — Phase Out RTKLIB, Ship Gneiss-Native PPP

## State at Start (2026-06-29)

### What works
- **RINEX auto-position**: Every IGS RINEX header parsed automatically, seeds PPP with cm-accurate coordinates. Zero config needed.
- **RTKLIB port**: 0.62m p95 CEDU, 0.97m HOB2. Meets 1m PPP goal on 2/4 IGS stations.
- **Native IEKF (first 500 epochs)**: 4mm median, 1cm p95 on CEDU. 100× better than RTKLIB.
- **Position smoother**: Implemented, negligible with tight prior, will be useful for kinematic.
- **RTK base stations**: Already use RINEX header position (verified matching).

### What's broken
- **Native IEKF diverges after epoch 500+**: Cascade AR feedback loop (same root cause as original RTKLIB IF mode divergence).
- **RTKLIB port state layout**: Incompatible with smoother, factor graph, FGO. Architectural dead end.
- **ALIC/YARR**: Don't meet 1m PPP goal (1.86m, 1.21m with RTKLIB; untested on native IEKF).
- **PPP kinematic**: 21m p95 Odaiba. Urban canyon multipath.
- **RTK**: 5m p95 Odaiba, 11m Shinjuku. 20-50× from 25cm goal.

### Key architectural insight
The gneiss-native 21-element core state `[pos(3), vel(3), att(3), acc_bias(3), gyro_bias(3), clk, isb_glo, isb_gal, isb_bds, clk_drift, zwd, ambiguities]` is the correct design. The RTKLIB layout `[pos(3), clk, tropo, biases]` was scaffolding. Every gneiss subsystem (smoother, factor graph, FGO, INS coupling) expects the 21-element state. We must complete the migration.

---

## Phase 1: Ship Native IEKF (1 session)

**Goal**: Native IEKF matches or exceeds RTKLIB port on all 4 IGS stations. RTKLIB port deprecated.

### Task 1.1: NL AR candidate search in cascade AR
**File**: `crates/gneiss-rtk/src/engine/ppp_ar.rs`
**Problem**: The cascade AR (LAMBDA WL + NL) doesn't constrain the integer search using position uncertainty. With tight position prior (σ=1cm), the search space is small — but LAMBDA searches the full integer space anyway. When a wrong integer passes the ratio test at epoch 500+ (constellation rotation), the soft lock (σ=10cm) limits but doesn't prevent drift.

**Fix**: Add a geometry-based NL validation step after LAMBDA:
- For each NL-fixed satellite, compute the implied N_IF from the fixed N1/N2
- Check that N_IF is within 0.2m of the float N_IF estimate
- If >50% of satellites fail this check, reject the entire AR fix
- This is equivalent to the RTKLIB port's check at ppp_rtklib.rs:817

### Task 1.2: Fix state history double-push
**File**: `crates/gneiss-rtk/src/engine/ppp.rs:247`, `processor/mod.rs:427`
**Problem**: `process_ppp` pushes to `state_history` internally AND the main loop pushes again. This causes duplicate epochs in output.

**Fix**: Remove the internal push from `process_ppp`. The main loop handles history management. Or gate: only push in `process_ppp` when called from the multi-epoch path.

### Task 1.3: Run full IGS benchmark suite
- All 4 IGS stations, 2880 epochs, native IEKF
- Compare vs RTKLIB port baseline
- Verify no divergence on any station
- Target: All stations <1m p95 (meeting PPP goal)

### Task 1.4: Wire native IEKF as default
- Change `EngineMode::Ppp` to route to `ppp::process_ppp` (native IEKF)
- Keep `EngineMode::PppRtklib` for regression testing
- Remove `force_if` and IF-mode-specific code from RTKLIB port (simplify)

**Acceptance**: All 4 IGS stations <1m p95 with native IEKF. No divergence. RTKLIB port still available for comparison.

---

## Phase 2: Multi-Epoch Factor Graph (2 sessions)

**Goal**: Native IEKF + multi-epoch joint optimization. Position shared across epochs, ambiguities shared across epochs. This directly addresses the remaining error sources (code multipath bias, cascade AR wrong fixes).

### Task 2.1: Fix PppTwoEpochOptimizer to use native IEKF output
**File**: `crates/gneiss-rtk/src/engine/ppp_multi_epoch.rs`
**Problem**: Currently calls `PppIteratedEkf` internally. Should receive state from outside.

**Fix**: 
- Remove the internal `iekf.solve()` call from `PppTwoEpochOptimizer::solve`
- Change signature to accept an already-solved `RtkState`
- The processor calls native IEKF first, then feeds state to optimizer

### Task 2.2: Wire optimizer into main processing loop
**File**: `crates/gneiss-rtk/src/engine/processor/mod.rs`
**Problem**: Multi-epoch optimizer exists but isn't called.
**Fix**: After each epoch's IEKF solve, snapshot state and feed to optimizer. After N epochs accumulate, run joint optimization, write smoothed position back.

### Task 2.3: Shared position parameters
**File**: `crates/gneiss-rtk/src/engine/ppp_multi_epoch.rs`
**Problem**: Currently uses dynamics constraints between epochs (x_k ≈ x_{k-1}). For static receivers, position should be a single shared parameter.
**Fix**: Modify state vector from `[epoch_0_core, ..., epoch_{N-1}_core, ambiguities]` to `[position(3), clock_0, tropo_0, ..., clock_{N-1}, tropo_{N-1}, ambiguities]`. Position is shared. Clocks, tropo, ambiguities are per-epoch (with ambiguity constancy constraints).

### Task 2.4: Benchmark and tune
- Compare multi-epoch vs single-epoch native IEKF on all 4 IGS stations
- Tune window size (2, 5, 10 epochs)
- Tune process noise (currently pos=0.1, too tight for convergence)

**Acceptance**: Multi-epoch improves p95 by ≥20% over single-epoch native IEKF on at least 3/4 stations.

---

## Phase 3: Integer AR Hardening (1 session)

**Goal**: Robust integer ambiguity resolution that doesn't diverge.

### Task 3.1: AR validation gate
**Problem**: Wrong AR fixes are the primary failure mode. Current ratio test (LAMBDA best/second-best) is the only validation.
**Fix**: Add multi-epoch consistency check:
- After AR fix, compute position from CP-only measurements (using fixed ambiguities)  
- If position jumps >0.5m from float position, reject fix
- Track per-satellite fix history: if a satellite's fixed ambiguity changes by >1 cycle within 10 epochs, mark as unreliable

### Task 3.2: Partial AR
**Problem**: When geometry is poor, fixing ALL ambiguities can force wrong integers on some sats.
**Fix**: Implement partial AR — fix only the subset of ambiguities that:
- Have WL convergence >100 samples
- Have NL LAMBDA ratio >3.0
- Are consistent with the tight position prior
Fix the reliable subset, leave unreliable sats as float.

### Task 3.3: WL-only mode for bootstrap
**Problem**: NL AR needs position to ~5cm, but position may be 2m off (without tight prior).
**Fix**: WL AR is geometry-free — works regardless of position error. After WL fixing:
- WL constraint reduces the NL search space from ±140 candidates to ±20
- NL can then be fixed with lower-confidence position
- This provides a bootstrap path when RINEX position isn't available (e.g., smartphone, kinematic)

**Acceptance**: Zero AR-induced divergence events across all 4 IGS stations, 2880 epochs. AR fix rate >50%.

---

## Phase 4: PPP Kinematic + Multi-Constellation (1-2 sessions)

**Goal**: Improve UrbanNav PPP from 21m p95 toward 1m goal.

### Task 4.1: Enable RINEX position for kinematic rover
**Problem**: Kinematic rover doesn't have a known position. But the FIRST epoch might have a good SPP fix.
**Fix**: Use SPP position with appropriate variance (σ=5m) as initial prior. The position smoother becomes important here — backward pass propagates converged information to early epochs.

### Task 4.2: Multi-constellation (GPS + Galileo + QZSS)
**Problem**: UrbanNav Tokyo has QZSS visible. Currently only GPS is used.
**Fix**: Enable Galileo and QZSS constellations in the config. Requires:
- Multi-constellation ISB handling (already in 21-element state)
- Per-constellation AR (already implemented in cascade AR)
- QZSS ephemeris loading (already in SP3 parser)

### Task 4.3: Elevation-dependent CP weighting
**Problem**: Low-elevation satellites have worse multipath in urban canyons.
**Fix**: Already partially done (C/N0-based variance scaling). Add:
- Elevation mask: exclude <15° satellites from CP measurements entirely
- Azimuth-dependent weighting: downweight satellites in directions with known obstructions (from ephemeris)

**Acceptance**: UrbanNav Odaiba PPP p50 <5m (from 7m), p95 <15m (from 21m).

---

## Phase 5: RTK Urban Canyon (1-2 sessions)

**Goal**: Improve RTK p95 from 5m toward 1m. 25cm goal requires INS coupling.

### Task 5.1: Multi-base selection
**Problem**: Single base station at 25km has ionospheric residuals. Multiple base stations could provide redundancy.
**Fix**: The multi-base pipeline (commit 4efe53e) already exists. Wire it into the RTK processor. Select the nearest base station or combine measurements from multiple bases.

### Task 5.2: Cycle slip detection hardening
**Problem**: GF cycle slip detection was disabled by default (now enabled). Still may miss slips in urban canyon.
**Fix**: Add Doppler-phase consistency check (already implemented, threshold was 5.0 cycles → tightened to 2.0). Add time-difference CP check: if ΔCP between consecutive epochs exceeds 0.5m after geometry removal, flag as slip.

### Task 5.3: AR validation for RTK
**Problem**: RTK AR ratio test (1.5-1.6) is permissive. Wrong fixes cause 5m+ errors.
**Fix**: Same validation as PPP Phase 3: position-jump check after AR, multi-epoch consistency.

**Acceptance**: RTK Odaiba p95 <3m (from 5m), Shinjuku p95 <8m (from 11m).

---

## Timeline

| Phase | Sessions | Cumulative | Key Metric |
|-------|----------|------------|------------|
| 1: Ship native IEKF | 1 | 1 | 4/4 IGS <1m p95 |
| 2: Factor graph | 2 | 3 | +20% p95 improvement |
| 3: AR hardening | 1 | 4 | Zero divergence events |
| 4: Kinematic PPP | 1-2 | 5-6 | Odaiba p50 <5m |
| 5: RTK hardening | 1-2 | 6-8 | Odaiba p95 <3m |

## Risk Register

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Native IEKF cascade AR unfixable | 20% | Must keep RTKLIB port for PPP static | Try partial AR first; fall back to WL-only AR |
| Factor graph doesn't improve over single-epoch | 30% | Phase 2 wasted | Validate on 2-epoch case first; abort if no gain |
| RTK fundamentally limited by multipath | 60% | Phase 5 wasted | Accept 1-3m p95 for urban canyon; document INS coupling as required for 25cm |
| RINEX position unavailable for kinematic | 100% | Phase 4 harder | SPP seed + WL bootstrap (Phase 3.3) provides alternative path |
