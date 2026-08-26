# Kinematic Processing Mode — Landing Report

Branch `exp/kinematic-mode`, worktree `.worktrees/kinematic`.
Commits: `ab22cc6` (feature) + `1a7a47c` (measured calibration fix).
Base: `63b0b29`. Binary hashes: baseline `8607dcdd265d`,
post-feature `4159a20dc6b4`, final (calibrated) `0ba61aa63deb`.

## 1. Mechanism

New module `crates/gneiss-rtk/src/post_process/dynamics.rs`:
`ProcessingDynamics { Static, Kinematic }`, selected by
`GNEISS_DYNAMICS=kinematic` (exact match; unset/garbage ⇒ Static,
byte-identical legacy). Threaded through `PostProcessOptions` into
forward pass, backward pass, combiner, and the network continuity
gate. The UPD pre-pass stays Static by design (satellite wide-lane
biases are motion-independent; static Q maximises arc convergence).

| Knob | Static (legacy, bitwise) | Kinematic |
|---|---|---|
| q_accel | 1e-6 (binaries); `None` fallback stays legacy 1.0 | 1.0 (~55 m/epoch re-randomization @30 s; ≥100× static; live velocity Q = dt·q) |
| Monument lock (two-phase 900 s → 1e-8) | active | suppressed (no monument to lock), fwd + bwd |
| Robust innovation knee | NIS > 9 (×1.0 bitwise) | ×3 ⇒ NIS > 27 before Huber inflation |
| Combiner disagreement limit | fixed 0.50 m | clamp(6σ_epoch, **0.50**, 10.0) m |
| Combiner one-sided fusion window | fixed 0.50 m | clamp(6σ, 0.50, 10.0) m |
| Combiner both-fixed fusion window | fixed 0.20 m | clamp(2σ, 0.20, 2.0) m |
| Continuity gate jump allowance | fixed 0.20 m / 90 s | max(0.20, 1.5·\|v\|·dt + 6σ) m |

σ_epoch = √(trace(Σ_pos)/3) of the passes' reported formal covariance.
Floors equal the audited static constants: the kinematic rule can only
**widen** validated tolerances, never tighten below them.

## 2. Tests (red→green, all green; workspace 719 passed / 0 failed)

| Requirement | Test | Red evidence |
|---|---|---|
| env unset ⇒ Static legacy Q | `env_unset_profile_selects_static_legacy_q_bitwise` (predict.rs): profile-resolved Q bitwise == hard-coded 1e-6 Q | new API: compile-fail red, then green |
| Kinematic Q ≥100×, velocity Q live | `kinematic_q_position_at_least_100x_static_and_velocity_nonzero`; exact dt³/3·q magnitude check | same |
| Gate scale semantics | update.rs: `robust_gate_scale_one_is_bitwise_legacy`, `kinematic_gate_keeps_nominal_weight_between_legacy_and_widened_knee` (NIS=18 nominal under ×3, inflated under ×1), gated-vs-legacy state equality | assertion-fail red at 1 ULP first, then green |
| Combiner σ scaling (synthetic sigmas) | combiner.rs: monotone growth, floor==audited bound, cap, never-stricter-than-static property, strict-off no-op | first version with 5 cm floor was **empirically falsified** by benchmark (see §3) — test constants updated with the calibration, mechanics tests unchanged |
| Continuity gate motion consistency | network.rs: 150 m travel at claimed 5 m/s survives; 10× overshoot downgraded + anchored; idle floor preserves static rule; Static delegates bit-identically | same pattern |
| Sim differential @30 s epochs | mod.rs: linear ramp — both profiles bounded (static 0.283 m, kin 1.110 m mean tail err; **finding**: static velocity states converge on clean CV ramps too); circular 15 m/s r=500 m — **static Q diverges to 26 236 m mean error, kinematic tracks at 0.738 m** | initial "static velocity frozen" hypothesis asserted, measured false, removed honestly |

## 3. Benchmark protocol (dataset A replay as pseudo-kinematic)

```
cargo build --release --bin eval_network_ppk && shasum -a 256 …   # 0ba61aa63deb
GNEISS_DATASET=multi2025 GNEISS_SYSTEMS=GE            … > static
GNEISS_DATASET=multi2025 GNEISS_SYSTEMS=GE GNEISS_DYNAMICS=kinematic … > kin
```

**Static-path unchanged proof:** every headline statistic AND every
per-base epoch-level fix count is identical to the pre-change run from
binary `8607dcdd265d` built at 63b0b29 (e.g. P181 sm 2839/2880, fused
2808/2880). Both guards (`check_network_benchmark.py`,
`check_multignss_benchmark.py`) report ALL CHECKS PASSED on the final
binary. Per-epoch quality flags unchanged ⇒ not merely aggregate parity.

**Paired table, static truth (P224 monument), horizontal metres:**

| Block | fix% static → kin | h p50 | h p95 | RMS |
|---|---|---|---|---|
| P181 15 km fwd | 97.3 → 97.2 | 0.101 → 0.101 | 0.142 → 0.144 | 0.106 → 0.108 |
| P181 15 km sm | 98.6 → 98.3 | 0.106 → 0.105 | 0.129 → 0.129 | 0.107 → 0.108 |
| P225 21.9 km fwd | 83.2 → 82.1 | 0.058 → 0.058 | 0.240 → **0.341** | 0.111 → **0.168** |
| P225 21.9 km sm | 73.4 → 73.5 | 0.056 → 0.055 | 0.214 → 0.238 | 0.135 → 0.125 |
| P222 38 km fwd | 88.0 → 89.1 | 0.161 → 0.165 | 0.341 → 0.355 | 0.195 → 0.204 |
| P222 38 km sm | 88.4 → 87.6 | 0.125 → 0.126 | 0.267 → 0.266 | 0.157 → 0.163 |
| NETWORK FUSED | 97.5 → **99.0** | 0.083 → **0.081** | 0.168 → 0.169 | 0.096 → 0.096 |

Reading: on a dataset whose rover never moves, a mobile prior costs a
little forward-pass float accuracy at the longest baselines (P225 fwd
p95) and buys mobility headroom everywhere else; smoothed products and
the fused product hold or improve.

**Measured calibration event (honest negative):** the first cut used a
5 cm floor for σ-scaled limits. WL_DUMP separation analysis on P181
showed honest fwd/bwd separations p50 9.6 cm / p90 22.6 cm / p99 40.6 cm
— i.e. formal sigmas are >10× optimistic — and smoothed fix rates
collapsed (P181 34.7%, P225 15.8%, **P222 1.9%**) while accuracy of
surviving fixes improved (selection effect). Fix: floors pinned to the
audited 0.50/0.20 m constants; scaling only widens. This is exactly the
combiner-threshold failure mode flagged in §5, caught by the benchmark
before landing.

## 4. Validation limits (stated plainly)

- **No ground-truth moving data exists offline in this repo.** The
  multi2025 kinematic run proves non-divergence and preserved accuracy
  when a STATIC dataset is processed through the mobile prior; it does
  NOT measure kinematic accuracy. The sim differential tests provide
  the only moving-truth evidence (linear ramp + sustained acceleration).
- AR behaviour under acceleration is exercised by simulation only;
  integer-fix integrity on real dynamic data (vibration, lever-arm
  swing, real cycle-slip statistics) is unvalidated.
- k=6σ steepnesses and continuity factors (1.5×|v|dt) are principled
  defaults, not tuned against measured kinematic disagreement
  distributions — blocked on moving-truth data.
- The Huber ×3 gate scale trades outlier sensitivity for dynamics
  tolerance; gross-blunder rejection now leans harder on the cycle-slip
  machinery (slip gate 500 cycles ≫ any plausible blunder that survives).

## 5. SELF RED TEAM

**AR under acceleration?**
- Kinematic Q inflates float ambiguity variance between 30 s epochs
  (position prior ~95 m/epoch), so LAMBDA ratio tests see noisier
  floats during dynamics: expect fix-rate dips in turns vs static
  processing of the same data. Mitigations already present: MW cascade
  veto is geometry-free (motion-independent) and still rejects FAR
  fixes contradicting converged wide-lanes; PAR falls back to subsets.
- Unresolved risk kept out of scope deliberately: `min_ar_lock_epochs`
  is disabled engine-wide; freshly-risen pairs entering AR immediately
  is riskier under fast geometry change. Follow-up: enable per-profile.
- The iono-free post-fix re-estimation assumes the fixed integers are
  right; under acceleration a confidently-wrong joint fix would bias
  it. The widelane veto covers wide-lane-contradictory cases only.

**Combiner threshold scaling failure modes?**
- *Over-optimistic sigmas* (realized): threshold collapses toward the
  floor; honest epochs mass-downgraded. Now bounded by pinning floors
  to audited constants — cost is that scaling is inert until 6σ exceeds
  0.50 m, i.e. honesty tightening below the static rule is impossible
  by construction. Acceptable until covariances are recalibrated.
- *Inflated sigmas*: cap 10 m bounds excuse-making, but 10 m is
  arbitrary; a wrong-but-diffuse pass pair within 10 m stays "fixed".
- *Both-fixed averaging*: capped at 2 m so diverged integer sets can't
  be averaged into a phantom position; between 0.20 m and 2 m the fuse
  decision trusts 2σ — unvalidated for real kinematic error tails.
- *Continuity gate self-anchoring*: allowance consumes reported
  velocity, so a consistently-wrong-but-smooth trajectory can anchor
  itself (same weakness class as the static gate, not worsened here).

## 6. Verdict: **SHIP** (behind flag)

Inert by default (bitwise-identical static path proven end-to-end),
structurally complete mechanism (Q, lock suppression, gates, combiner,
continuity), differential sim proof of necessity (26 km → 0.74 m under
sustained acceleration), guards green, honest failure-mode ledger above.
NOT production-ready for revenue kinematic workflows until validated
against true moving-truth data; the flag boundary is the product line.
