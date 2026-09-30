# Survey Report: R3 Doppler-Assisted Cycle Slip Detection & R4 Prioritized Partial Ambiguity Resolution (PAR)

**Author**: Survey Explorer 3  
**Date**: 2026-09-24  
**Working Directory**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_3`  
**Parent Task**: Urban Canyon Fix Rate Expansion and Multipath Mitigation (`ORIGINAL_REQUEST.md` @ 2026-09-24T13:30:49Z)  
**Target File**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md`  

---

## Executive Summary

This survey maps the existing codebase baseline and defines the concrete mathematical, algorithmic, and architectural implementation points for:
- **Requirement R3**: Doppler-assisted cycle slip detection & phase continuity validation across epochs.
- **Requirement R4**: C/N0-, elevation-, lock-duration-, and CMC-prioritized Partial Ambiguity Resolution (PAR) with geometry (DOP) guards and strict covariance positive-definiteness ($Q_{aa} \succeq 0, P \succ 0$).

### Key Findings & Existing Baseline Deficiencies

1. **Doppler Slip Detection Exists Only in Pre-Processing Screening, Not in Real-Time DD Formation**:
   - `crates/gneiss-rtk/src/post_process/screening.rs` contains an initial implementation of `CycleSlipDetector::check_band_doppler_slips`, but its detection threshold is set to `(1.0 * dt).max(1.0)` cycles. This threshold **completely misses** half-cycle slips ($0.5$ cycles) and frequently fails to catch 1-cycle slips under noisy urban conditions.
   - It only checks frequency bands 1, 2, and 7, ignoring band 5 (Galileo E5a) and band 6 (BeiDou B3I).
   - In GLONASS, it subtracts a single scalar median across all satellites, ignoring FDMA inter-channel frequency offsets ($f_i = 1602 + k_i \times 0.5625\text{ MHz}$) which cause receiver clock drift in cycles to vary per satellite.
   - It is inactive if fewer than 3 satellites are tracked on a band.

2. **Severe Slip Flag Propagation Leaks in `rtk_iekf`**:
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs:416`, base-side slip detection is **only checked when `widelane_ar` is active**. In standard short/medium baseline RTK, base slips are never evaluated!
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs:304`, the Melbourne-Wübbena tracker `WidelaneTracker::update_tracker_from_obs` checks **only raw observation LLI bit 0**. It is never informed of Doppler-detected slips or GF slips! The MW arc continues accumulating corrupted phase for up to 5 epochs before an innovation gate can trigger.
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs:336`, the phase-only wide-lane tracker `pw_tracker` is hardcoded to `slip: false`.
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs:161`, tracking arc duration `pair_epochs` is **never reset upon cycle slip**. Slipped satellites continue to be treated as converged by `min_ar_lock_epochs`.
   - The phase innovation gate `PHASE_INNOVATION_GATE_CYCLES` in `robust.rs` is set to $500.0$ cycles, functioning only as a catastrophic failure guard.

3. **Existing PAR Implementation Ignores Physical Signal Quality**:
   - `crates/gneiss-rtk/src/ambiguity/par.rs:select_ils_subset` sorts candidates purely by diagonal float variance $Q_{ii}$ and uses a naive greedy loop that breaks upon the first candidate failure. It lacks DOP geometry guards and minimum subset size guards.
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs:select_par_candidates` sorts candidates by $score = |a_i - [a_i]| + 0.5\sqrt{Q_{ii}}$. It ignores elevation, C/N0, tracking lock duration, and CMC multipath residual variance.
   - In urban canyons, a multipath-corrupted satellite can have near-zero fractional float error ($|a_i - [a_i]| \approx 0$) by coincidence, while a clean high-elevation satellite with high C/N0 has slight bias ($0.15$ cyc) due to common filter pull. The existing sorter prioritizes the corrupted satellite, leading to ratio test rejection or false integer fixes.

4. **Off-Diagonal Covariance Leakage & Loss of Positive Definiteness**:
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs:condition_state_on_integers`, when ambiguities are fixed, their diagonal variances are clamped to `1e-4`, but their off-diagonal cross-covariances with other states are **not zeroed**. This violates the Cauchy-Schwarz inequality ($|P_{ij}| \le \sqrt{P_{ii} P_{jj}}$) and produces an indefinite state covariance matrix ($P \not\succ 0$).
   - In contrast, `crates/gneiss-rtk/src/composite/tc_ambiguity.rs:340` properly zeroes out row and column cross-covariances.
   - `crates/gneiss-core/src/dop.rs` provides `compute_dop_from_positions` which must be integrated as a hard geometry guard before subset acceptance.

---

## Part I: R3 — Doppler-Assisted Cycle Slip Detection & Phase Continuity Validation

### 1. Doppler Measurement Representation & Physical Invariants

In GNSS signal processing, Doppler frequency shift $f_D$ (in Hz) is the difference between received carrier frequency $f_{\text{rx}}$ and nominal carrier frequency $f_0$:
$$f_D = f_{\text{rx}} - f_0 = -\frac{1}{\lambda} \dot{\rho}_{\text{geom}} + f_0 (\dot{\delta t}_{\text{rx}} - \dot{\delta t}_{\text{sat}})$$

- **Sign Invariant**: When the receiver and satellite approach each other ($\dot{\rho} < 0$), the received frequency increases ($f_D > 0$).
- **Carrier Phase Observable**: In RINEX and u-blox raw observables, carrier phase $\Phi$ is expressed in cycles, where $\Phi = \frac{1}{\lambda} \rho + f_0 (\delta t_{\text{rx}} - \delta t_{\text{sat}}) + N - I + T$.
- **Carrier Phase Rate**:
  $$\dot{\Phi} = \frac{d\Phi}{dt} = -f_D \quad [\text{cycles/s}]$$
  $$\dot{\rho} = -\lambda f_D \quad [\text{m/s}]$$

#### Codebase Storage & Retrieval Baseline:
- `crates/gneiss-core/src/obs.rs`:
  - `ObsType::Doppler` ('D') in `ObsCode`.
  - `SatObs::get_doppler(freq_band: u8) -> Option<f64>` retrieves Doppler in Hz.
  - `SatObs::get_observable_phase(freq_band: u8) -> Option<f64>` retrieves carrier phase in cycles.
  - `SatObs::get_lli(freq_band: u8) -> Option<u8>` retrieves Loss of Lock Indicator.
- `crates/gneiss-parsers/src/ubx/mod.rs:217`:
  - Parses `UBX-RXM-RAWX` payload `meas.do_mes` as `ObsType::Doppler`.
- `crates/gneiss-parsers/src/rinex/obs/mod.rs:73`:
  - Parses standard RINEX 'D' observable records.
- `crates/gneiss-rtk/src/estimators/doppler.rs:107-109`:
  - Computes $\dot{\rho}_{\text{meas}} = -\lambda \cdot f_D$ and compares against satellite line-of-sight velocity:
    `let rho_meas = -lambda * doppler_hz;`
    `let rho_sat = e_los.dot(&sat_vel) - c_sat_drift;`

---

### 2. Phase Increment vs. Integrated Doppler Range Rate Comparison

#### Discrete Integration Mathematics:
Between epochs $k-1$ and $k$ separated by interval $\Delta t = t_k - t_{k-1}$ (typically $0.1\text{ s}$ to $1.0\text{ s}$):
- Observed phase increment:
  $$\Delta \Phi_k = \Phi_k - \Phi_{k-1} \quad [\text{cycles}]$$
- Trapezoidal integration of Doppler range rate:
  $$\Delta \Phi_{\text{dopp}, k} = \int_{t_{k-1}}^{t_k} \dot{\Phi}(t) dt \approx -\frac{f_{D, k} + f_{D, k-1}}{2} \Delta t \quad [\text{cycles}]$$
- Discrepancy per satellite $i$:
  $$D_i = \Delta \Phi_{i, k} - \Delta \Phi_{\text{dopp}, i, k} = \Delta N_i + f_{0, i} \Delta \delta t_{\text{rx}} + \epsilon_{\text{int}}$$

#### Receiver Clock Drift Isolation:
In low-cost receivers (e.g. u-blox F9P TCXO), receiver clock drift rate $\dot{\delta t}_{\text{rx}}$ ranges between $\pm 1\text{ to } 5\text{ ppm}$ ($\sim 1.5\text{ to } 7.8\text{ kHz}$ at L1), contributing several thousand cycles per second common-mode shift.
- To cancel $\Delta \delta t_{\text{rx}}$ across CDMA satellites on the same frequency band:
  $$r_i = D_i - \text{median}_{j \in \mathcal{S}_{\text{band}}}(D_j)$$
  where $\mathcal{S}_{\text{band}}$ is the set of satellites tracked on that band.
- For GLONASS FDMA satellites, where $f_{0, i} = f_0 + k_i \Delta f$:
  The receiver clock drift in seconds is first estimated:
  $$\Delta \delta t_{\text{rx}} = \text{median}_{j \in \text{CDMA}}\left( \frac{D_j}{f_{0, j}} \right)$$
  Then for GLONASS satellite $i$:
  $$r_i = D_i - f_{0, i} \cdot \Delta \delta t_{\text{rx}}$$

#### Current Implementation in `screening.rs`:
```rust
// crates/gneiss-rtk/src/post_process/screening.rs:87-91
let d_phi = cp - prev_cp;
let dop_avg = 0.5 * (dop + prev_dop);
let pred_d_phi = -dop_avg * dt;
discs.push((sat, d_phi - pred_d_phi, dt));

// crates/gneiss-rtk/src/post_process/screening.rs:150-154
let thresh = (1.0 * dt).max(1.0);
if residual > thresh {
    slipped.push(sat);
}
```

#### Deficiencies of Baseline:
1. **Threshold Over-Conservatism**: `thresh = (1.0 * dt).max(1.0)`:
   - Half-cycle slips ($0.5$ cycles) produce $r_i \approx 0.5$. With threshold $\ge 1.0$, they are **100% missed**.
   - 1-cycle slips with slight measurement noise (residual $\approx 0.88\text{--}0.95$) are missed.
2. **Missing Constellation Bands**: Only checks `[1, 2, 7]`. Bands 5 (Galileo E5a) and 6 (BeiDou B3I) are ignored.
3. **Double-Difference Doppler Bypass**: In RTK double-difference processing (`formation.rs`), double-difference Doppler:
   $$\nabla\Delta \dot{\Phi} = -\nabla\Delta f_D$$
   is completely unformed and unused. In double-difference, receiver and satellite clock drift **cancel out analytically without needing median estimation**.

#### Proposed High-Precision Doppler Detector:
1. **Single-Receiver Screening Threshold**:
   $$\text{Thresh}_{\text{slip}}(\Delta t) = \max\left(0.28, \; 0.20 + 0.10 \cdot \Delta t\right) \quad [\text{cycles}]$$
   For $\Delta t = 1.0\text{ s}$, $\text{Thresh} = 0.30$ cycles. This provides a $6\sigma$ clearance above Doppler noise ($\sigma \approx 0.035$ cyc) while reliably catching half-cycle ($0.50$ cyc) and 1-cycle ($1.0$ cyc) slips.
2. **Double-Difference Phase-Rate Audit**:
   Form double-difference Doppler in `formation.rs`:
   $$\nabla\Delta D = (\Phi_{\text{rov}}^s - \Phi_{\text{rov}}^r) - (\Phi_{\text{bas}}^s - \Phi_{\text{bas}}^r)$$
   $$\nabla\Delta \dot{\Phi} = -(f_{D, \text{rov}}^s - f_{D, \text{rov}}^r) + (f_{D, \text{bas}}^s - f_{D, \text{bas}}^r)$$
   Residual: $|\Delta \nabla\Delta \Phi_k - \overline{\nabla\Delta \dot{\Phi}} \Delta t| > 0.35\text{ cycles} \implies \text{Slip}$.

---

### 3. Slip Flag Propagation & Covariance Re-seeding Pipeline

```
[Raw EpochObs (Rover & Base)]
            │
            ▼
┌──────────────────────────────────────────────┐
│ CycleSlipDetector::check_epoch               │
│ - LLI Bit 0 check                            │
│ - Geometry-Free (GF) phase step (> 5 cm)     │
│ - Time gap check (> 2 s)                     │
│ - Enhanced Doppler residual check (> 0.30 cyc)│
└──────────────────────┬───────────────────────┘
                       │
       Increment sat_arc counter
                       │
                       ▼
┌──────────────────────────────────────────────┐
│ rtk_iekf::formation::check_pair_slip         │
│ - Compare prev_arcs with cur_arc             │
│ - MUST check both rover AND base slip counts │  <-- FIX: Remove `if widelane_ar` gate
└──────────────────────┬───────────────────────┘
                       │
                  `lli_slip` flag
                       │
       ┌───────────────┴───────────────┐
       ▼                               ▼
┌──────────────────────────────┐ ┌─────────────────────────────────────────┐
│ rtk_iekf::update_dd_ambiguity│ │ mw::update_tracker_from_obs / pw_tracker│
│ - Reset state.ambiguities    │ │ - MUST reset MwTrack when `lli_slip`!   │ <-- FIX: Pass `lli_slip`
│ - Re-seed cov:               │ │ - Drop accumulated 20-epoch average     │
│   P[(i,i)] = (4*sigma_pr/lam)│ └─────────────────────────────────────────┘
│ - Zero cross-covariances     │
│ - Reset pair_epochs = 0      │  <-- FIX: Reset tracking duration to 0
└──────────────────────────────┘
```

#### Ambiguity Covariance Re-seeding Formulation:
When a slip is detected on pair `key` with wavelength $\lambda$ and double-difference pseudorange variance $\sigma^2_{\nabla\Delta P}$:
- Initial float ambiguity seed:
  $$\hat{a}_{\text{seed}} = \nabla\Delta \Phi - \frac{\nabla\Delta P}{\lambda}$$
- Prior variance calculation (`ar_gate.rs:56-59`):
  $$\sigma_{\text{seed}} = \frac{\sqrt{\sigma^2_{\nabla\Delta P}}}{\lambda}$$
  $$P_{aa, \text{seed}} = \left( k_{\text{margin}} \cdot \sigma_{\text{seed}} \right)^2 \quad \text{with } k_{\text{margin}} = 4.0$$
- In `RtkState::reset_ambiguity` (`state/mod.rs:361-372`):
  - Row and column cross-covariances are set to 0.0:
    `self.cov[(r, idx)] = 0.0; self.cov[(idx, r)] = 0.0;`
  - Diagonal variance is set to $P_{aa, \text{seed}}.\max(1.0)$.
- **Missing Action to Implement**:
  - `self.pair_epochs.insert(key, 0);` to immediately reset the lock duration counter so `min_ar_lock_epochs` will properly prevent premature AR fixing on the freshly seeded arc.

---

## Part II: R4 — C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR)

### 1. Existing Ambiguity Machinery Baseline

| Component | File Path | Current Mechanism | Critical Weakness |
|---|---|---|---|
| **PAR ILS Subset** | `crates/gneiss-rtk/src/ambiguity/par.rs` | `select_ils_subset`: sorts by diagonal variance $Q_{ii}$; greedy break-on-first-failure | No DOP/geometry checks; no min subset size guard; breaks on first drop |
| **TC-PAR** | `crates/gneiss-rtk/src/composite/tc_ambiguity.rs` | `solve_integers_full_or_par`: calls `select_ils_subset(0.995)`; checks `MIN_PAR_SUBSET_SIZE` (4) | Dependent on flawed `par.rs` sorting |
| **IEKF-PAR Candidates** | `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs` | `select_par_candidates`: $score = |a - [a]| + 0.5\sqrt{Q_{ii}}$ | Ignores elevation, SNR, lock duration, CMC multipath |
| **IEKF-PAR Search** | `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs` | `try_partial_ar`: prefix subsets, constellation clusters, 1-omission, 2-omission | Tests subsets without checking subset geometry (DOP) |
| **Fix-and-Hold Update** | `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs` | `condition_state_on_integers`: $P_{\text{new}} = P - P_{xa} Q_{aa}^{-1} P_{ax}$; clamps diagonal to $10^{-4}$ | **Off-diagonal leakage**: does not zero cross-covariances; $P$ becomes indefinite |

---

### 2. Multi-Metric Candidate Prioritization (Elevation, C/N0, Lock Duration, CMC Variance)

In urban canyons, float ambiguity fractional offsets can be artificially small on corrupted satellites. The candidate pool must be sorted by a physically grounded **Composite Quality Metric (CQM)**.

For each candidate ambiguity $i \in \{1, \dots, n\}$ associated with double-difference satellite pair $(s_i, r)$:

#### Factor 1: Satellite Elevation Weight ($w_{\text{el}}$)
Low elevation angles suffer longer atmospheric paths and higher building obstruction/multipath probability:
$$w_{\text{el}}(\theta_i) = \sin^2\left(\max(\theta_i, \; 10^\circ)\right)$$

#### Factor 2: Receiver Carrier-to-Noise Ratio ($w_{\text{snr}}$)
Direct LOS signals in urban environments typically have $S_i \ge 38\text{--}45\text{ dB-Hz}$, while NLOS/reflected signals drop below $32\text{ dB-Hz}$:
$$w_{\text{snr}}(S_i) = \frac{1}{1 + \exp\left(-0.25 \cdot (S_i - 35.0)\right)}$$
Values: $S=45 \implies 0.92$, $S=35 \implies 0.50$, $S=25 \implies 0.076$.

#### Factor 3: Continuous Tracking Lock Duration ($w_{\text{lock}}$)
Ambiguities with long continuous tracking have converged float filters and negligible filter initialization bias:
$$w_{\text{lock}}(N_{\text{epochs}}) = \min\left(1.0, \; \frac{N_{\text{epochs}}}{30.0}\right)$$
Reaches 1.0 after 30 epochs (e.g. 6 seconds at 5 Hz, or 30 seconds at 1 Hz).

#### Factor 4: Code-Minus-Carrier (CMC) Multipath Variance ($w_{\text{cmc}}$)
$$\text{CMC}_k = P_k - \Phi_k \lambda - 2 I_k$$
On clean satellites, sample variance $\sigma^2_{\text{cmc}} \le 0.25\text{ m}^2$. Under severe multipath, $\sigma^2_{\text{cmc}} \ge 10\text{--}50\text{ m}^2$:
$$w_{\text{cmc}}(\sigma^2_{\text{cmc}, i}) = \frac{\sigma_0^2}{\sigma_0^2 + \sigma^2_{\text{cmc}, i}} \quad \text{with } \sigma_0^2 = 0.50\text{ m}^2$$

#### Factor 5: Mathematical Float Quality ($w_{\text{math}}$)
Combines diagonal float variance $Q_{ii}$ and fractional deviation from the nearest integer:
$$w_{\text{math}}(a_i, Q_{ii}) = \frac{1}{\sqrt{Q_{ii}} \cdot \left(1.0 + 3.0 \cdot |a_i - \text{round}(a_i)|\right)}$$

#### Composite Quality Score:
$$\mathcal{Q}_i = w_{\text{el}}(\theta_i) \cdot w_{\text{snr}}(S_i) \cdot w_{\text{lock}}(N_i) \cdot w_{\text{cmc}}(\sigma^2_{\text{cmc}, i}) \cdot w_{\text{math}}(a_i, Q_{ii})$$

**Ordering Rule**: Candidates are ranked in descending order of $\mathcal{Q}_i$. The highest-quality, physically verified satellites are prioritized for integer fixing.

---

### 3. Geometry and DOP Guard

A subset of ambiguities, even if highly confident, must never be fixed if the remaining satellites form a degenerate geometry (e.g. all satellites along one collinear azimuth line).

1. **Minimum Subset Size Guard**:
   $$k_{\text{subset}} \ge 4 \quad \text{(minimum 4 satellites / 3 DD pairs for 3D fix)}$$
2. **DOP Computation**:
   Using `gneiss_core::dop::compute_dop_from_positions(rover_pos, &subset_sat_positions)`:
   - Compute PDOP of the candidate subset.
   - **DOP Guard Rule**:
     $$\text{PDOP}_{\text{subset}} \le \text{PDOP}_{\max} \quad (\text{default: } 6.0)$$
     $$\text{PDOP}_{\text{subset}} \le 2.0 \cdot \text{PDOP}_{\text{full}}$$
   If the subset violates the DOP guard, the candidate subset is rejected, preventing geometric explosion of the conditional position solution.

---

### 4. Preserving Covariance Definiteness ($Q_{aa} \succeq 0, P \succ 0$) and Eliminating Off-Diagonal Leakage

When a subset of $k$ ambiguities is fixed to integers $\check{a}_{\text{sub}}$:

#### 1. Ambiguity Sub-Covariance Definiteness ($Q_{aa} \succeq 0$):
When extracting submatrix $Q_{\text{sub}} \in \mathbb{R}^{k \times k}$:
- Symmetrize: $Q_{\text{sub}} = \frac{1}{2}(Q_{\text{sub}} + Q_{\text{sub}}^T)$.
- Apply regularized spectral floor via `crates/gneiss-rtk/src/ambiguity/lambda/mod.rs:regularize_positive_definite`:
  $$\lambda_{\min}(Q_{\text{sub}}) \ge 10^{-9}$$
  This prevents singular LDL^T decomposition and ill-conditioned LAMBDA search volume distortions.

#### 2. Position Projection Definiteness ($P_{xx} \succ 0$):
In `ar.rs:project_subset_fixed`:
$$P_{xx|\text{fix}} = P_{xx} - P_{xa} Q_{\text{sub}}^{-1} P_{ax}$$
- Must check eigenvalue floor:
  $$\lambda_{\min}(P_{xx|\text{fix}}) \ge 10^{-6} \text{ m}^2$$

#### 3. State Covariance Conditioning Without Off-Diagonal Leakage:
In `ar.rs:condition_state_on_integers`:
Let the full covariance be updated by the Kalman reduction:
$$\Delta P = P_{*, \text{sub}} \cdot Q_{\text{sub}}^{-1} \cdot P_{\text{sub}, *}$$
$$P_{\text{post}} = P - \Delta P$$
To avoid off-diagonal leakage and preserve $P_{\text{post}} \succ 0$:
```rust
// For every fixed ambiguity index `idx`:
for r in 0..state.dim() {
    new_cov[(r, idx)] = 0.0;
    new_cov[(idx, r)] = 0.0;
}
new_cov[(idx, idx)] = 1e-4; // Small pseudo-variance
```
- Validate minimum eigenvalue across the non-fixed subspace:
  $$\lambda_{\min}(P_{\text{non-fixed}}) \ge 10^{-6}$$
  If $\lambda_{\min} < 10^{-6}$, reject the conditioning update and retain the float solution.

---

## Part III: Codebase Implementation Mapping

### Detailed Modification Map

| Module / File | Changes Required | AGENTS.md Standards Compliance |
|---|---|---|
| `crates/gneiss-core/src/obs.rs` | Add helper `SatObs::get_primary_doppler() -> Option<f64>` and multi-band query helper. | LOC < 400, no unwrap |
| `crates/gneiss-rtk/src/post_process/screening.rs` | Update `check_band_doppler_slips` to cover bands 1, 2, 5, 6, 7. Lower threshold to $0.30$ cycles with adaptive scaling. Fix GLONASS FDMA clock drift subtraction. | Split helper functions to keep fn $\le 32$ LOC |
| `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs` | 1. In `check_pair_slip`, evaluate `base_slip_detector` unconditionally (remove `if widelane_ar`).<br>2. In `update_dd_ambiguity`, reset `self.pair_epochs.insert(key, 0)` on slip.<br>3. In `handle_widelane_phase_update`, pass `lli_slip` to `pw_tracker.update`. | Keep functions strictly $\le 32$ LOC |
| `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs` | Modify `update_tracker_from_obs` to accept `slip: bool` parameter from `formation.rs` instead of only checking LLI. | Function size $\le 32$ LOC |
| `crates/gneiss-rtk/src/ambiguity/par.rs` | 1. Replace naive variance-only sort with multi-metric ordering.<br>2. Replace break-on-first greedy loop with omission / prefix search.<br>3. Add minimum subset size ($\ge 4$) guard.<br>4. Integrate `compute_dop_from_positions`. | New file size $\approx 350$ LOC, fn $\le 32$ LOC |
| `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_subsets.rs` | Enhance `select_par_candidates` to take elevation, SNR, lock duration, and CMC variance from measurement bundle. Implement composite quality metric scoring. | Function sizes $\le 32$ LOC, nesting depth $< 3$ |
| `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs` | 1. Integrate DOP guard check in `eval_par_subset`.<br>2. In `condition_state_on_integers`, zero out fixed ambiguity row/column cross-covariances to eliminate off-diagonal leakage.<br>3. Check $\lambda_{\min}(P) \ge 10^{-6}$ before committing update. | Existing file 463 LOC; extract conditioning helper if needed to stay $< 500$ LOC |

---

## Part IV: Verification & Testing Strategy

To ensure zero regressions and satisfy the Three-Tier Verification Standard:

1. **Tier 1 (Analytical Golden Vectors)**:
   - Unit test Doppler cycle slip detector with synthetic carrier phase increments having exact $0.5$ cycle, $1.0$ cycle, and $2.0$ cycle slips.
   - Verify that clean signals with Doppler noise $\le 0.1$ Hz pass without false slip alarms.
   - Unit test PAR candidate sorting on synthesized multi-satellite scenarios (e.g. high elevation + high SNR vs low elevation + low SNR + multipath).
   - Unit test covariance conditioning: verify all cross-covariances are 0.0 and $\lambda_{\min}(P) \ge 10^{-6}$.

2. **Tier 2 (Numerical Finite Differences & Invariants)**:
   - Verify positive semi-definiteness ($Q_{aa} \succeq 0$) on all extracted PAR subsets using symmetric eigenvalue decomposition.
   - Verify that DOP guard correctly rejects ill-conditioned satellite geometries (coplanar / near-collinear satellites).

3. **Tier 3 (Closed-Loop Benchmarks & Guard Scripts)**:
   - Verify that both CI smoke guard scripts pass with zero warnings:
     - `python3 scripts/check_network_benchmark.py --smoke`
     - `python3 scripts/check_multignss_benchmark.py --smoke`
   - Run `eval_f9p_rover` across UrbanNav datasets (Shinjuku, Whampoa, TST1) and verify fix rate expansion and $p_{95}$ tail error reduction without introducing false fixes.
