# Sprint 61 — session handoff

Written at the request to compact context and resume. Everything a fresh session
needs to continue without re-deriving it.

## Repo state

| | |
|---|---|
| branch | `main` |
| last commit | `8cdad51` — "test: coverage sprint 81.9% -> 90.2%; six production bugs found and fixed" |
| unpushed | no (main is in sync with `origin/main` at `8cdad51`) |
| tags | `checkpoint/sprint-61` (pushed) |
| review branches | `experiments/f9p-broken-constellation-screen`, `experiments/f9p-glonass-ambiguity-exclusion` (both pushed) |

**IN-FLIGHT, UNCOMMITTED:** ~26 modified + 8 untracked files from three running
subagents. The tree does **not** build at the moment of writing (2 issues) because
the agents are mid-edit, and there is 1 staged deletion. **Do not commit or
revert this work blindly** — wait for the agents to finish.

## Where the objective stands

**Not met.** Urban canyon p95 is **4.8–24.7 m** against Tier-1 targets of
**0.15–3.50 m**. CORS/PPP/datum-tie rows do meet or exceed Tier-1 spec; the six
canyon rows are the entire remaining gap.

## What has actually been achieved

Two benchmark-moving changes:
- `31dee7e` per-epoch median-relative NLOS code gate — mean p95 −10.5% across six
  UrbanNav datasets; TST1 Survey p95 −40% (13.739 → 8.185 m)
- `e985330` GLONASS broadcast ephemerides were propagated as Kepler elements — a
  real correctness bug; F9P went from resolving **0%** of ambiguities to **69%**

Plus: Wilcoxon + Weibull + paired bootstrap CDF band (`gneiss_core::stats`),
`eval_compare`, two repaired CI guards, a stale-binary defect in a guard, and
`docs/TIER1_ROADMAP.md` §1 corrected to quote all-epoch CDF instead of
fixed-subset numbers against commercial specs.

## The single most important finding

**Several zero-false-fix and frame-safety claims rest on tests that cannot fail.**

- `swfg/engine/setup_tests.rs:236` **names a test that was never written**
  (`imu_residual_rotates_only_the_prediction`) and diagnoses the exact defect in
  prose, while the helper above it forces `UnitQuaternion::identity()` so the bug
  cancels. All 8 IMU-residual tests route through it.
- `tests/src/inertial_outage_simulation.rs:31` bakes in the same bug for
  `smoother.rs` by pre-rotating the delta into the opposite frame.
- `gneiss_core::frames` — a 2,620-LOC typed frame module — is imported by
  **zero of 157** library files in `gneiss-rtk`. `FRAME_SAFETY_PLAN.md` explicitly
  scoped out "Attitude/body frames", and `PROJECT_STATUS.md:138` logged it done
  anyway. Its success criterion "compiler prevents mixing incompatible frames" is
  **false**.
- `AGENTS.md` cites `predictor.rs:86-91` for sign conventions; **that file does
  not exist**, and its `d_theta = -psi` claim appears nowhere in the code.

**Coverage is not the gap — input diversity is.** Every preintegration test used
identity attitude, so body ≡ ECEF and any frame error is arithmetically invisible.

## Open defect (committed as a RED test)

`crates/gneiss-rtk/tests/frame_red_test.rs` is `#[ignore]`d with a re-enable
condition. IMU preintegration mixes frames: free-fall position residual is 17.68 m
where 0 is required. **Not fixable without an owner decision** — `dp`/`dv` doc
says ECEF, the arithmetic produces body-0, and three consumers independently
chose different frames. Fixing it changes published Odaiba INS numbers.
See `docs/P95_TAIL_ROOT_CAUSE.md` and section 11.1 of that file.

## Known-but-unfixed (each found during the sprint, deliberately not fixed)

1. `post_process/backward.rs` — reversing the IMU sample list makes
   `time_diff_us` read the GPS-week wrap and clamp `dt` to 0.1 s per step, so a
   1 s epoch integrates as 10 s. `dp`/`dv` are ~10× too large.
2. `eskf/alignment.rs` — `clamp(1, max_samples)` returns 1 for an empty slice,
   then `&imu_samples[..1]` panics.
3. `EngineConfig::PppIns` never sets `is_kinematic`.
4. `swfg/imu_preintegration/smoother.rs` — **zero production callers**, dead code
   (violates AGENTS.md). Deleting it removes 2 of 3 frame-defect sites free.
5. `DdCarrierPhaseFactor::information` floors at 1e-4 while the zenith DD carrier
   variance is 7.2e-5 — the floor binds at the best geometry, de-weighting phase
   by 28% exactly where the model is most trustworthy.
6. `resolve_dd_ambiguities` omits the DD troposphere; exact only while |T| < 53.5 mm.

## The F9P finding (do not re-litigate)

With a paired bootstrap CDF band (2000 reps), F9P AR-on vs float:
signs **cross between p68 and p75**. Every level is significant:
AR better at p10–p68, float better at p75–p99. p95 regression is real
(+0.129 m, CI [+0.112, +0.142]), not noise.
Per-epoch fallback was tried and is **strictly dominated** (best selector 57.6%
accurate; "float everywhere" has zero error and beats it).
**Do not use fix rate as a gate** — it reported a `1.4e-208` "win" for a change
that costs 81 mm at p95.

## Rules learned here at cost

1. **Validate on the raw all-epoch CDF first**, with intervals. Fix rate alone
   will ship a regression.
2. **Coverage is not safety.** Ask whether a test can fail at a non-degenerate
   input.
3. **When an agent corrects the brief, believe it.** Three agents did; all were
   right (UDC basis det is −1 not +1; `R_z(pi)` maps X→−X not Y).
4. **A `sed -i '' Ns/.../...` that matches nothing exits 0 and prints nothing.**
   This cost three wrong conclusions (rounds 8, 9, 12). Locate by content, assert
   the value changed.
5. **Never assert a configuration default as an invariant** — AR enable/disable
   is still an open question.

## Recommended next steps, in order

1. Wait for the three running subagents (core/time+signal, rtcm3+lnav parsers,
   composite+streaming) to finish. Commit their work. Re-verify:
   `cargo test --workspace --no-fail-fast`, `cargo clippy --workspace --all-targets`,
   and the three guard scripts.
2. **Audit the new tests the same way as the last batch** — 300+ agent-written
   tests is exactly where slop hides. Re-derive a sample of numeric claims
   independently.
3. Fix the frame convention (needs an owner decision), then green
   `frame_red_test.rs`.
4. Delete the dead `imu_preintegration/smoother.rs`.
5. Fix the backward-pass `dt` inflation — it is ~10× on real preintegration.
6. Close the calibration-convergence gap (round 5 finding, still open) and the
   Odaiba SWFG routing question (`forward.rs` only takes the IEKF path when
   `imu_samples.is_none()`, so Odaiba is routed to SWFG).

## Verification commands

```bash
cargo test --workspace --no-fail-fast
cargo clippy --workspace --all-targets          # must be 0 warnings
cargo llvm-cov --workspace --lib --summary-only # was 90.17% at 8cdad51
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
python3 scripts/check_f9p_benchmark.py --smoke   # FAILS by design: fix rate floor
cargo run --release --bin eval_f9p_rover -- all # ~10 min, 6 UrbanNav datasets
```

---

# ADDENDUM — impact-audit agents (rounds 16-17)

Three agents audited high-impact code NOT for coverage but for bugs. **16 further
production defects found.** All uncommitted; nothing integrated yet.

## gneiss-parsers — the correction path was badly broken

| # | Defect | Impact |
|---|---|---|
| P1 | **RTCM 1019 discarded the 10-bit GPS week** (`rtcm3/ephemeris.rs:39`), stamping `GpsTime::new(0,..)` — week 0 = **1980-01-06**. No 10-bit rollover resolution existed anywhere. | Every broadcast ephemeris dated ~46 years in the past. Now `adjust_gps_week(raw, reference)` + `parse_1019_at`. |
| P2 | **MSM epochs all stamped week 0** (`msm/decoder.rs:12`) | Same failure for MSM corrections. |
| P3 | **MSM GLONASS epoch used the wrong field layout** (RTKLIB: DF393 day-of-week + DF392 time-of-day + `adjday_glot`); **BeiDou missing its 14 s BDT offset** | GLONASS/BeiDou correction epochs wrong. |
| P4 | **SSR accepted every message type** — per-satellite record width is constellation-dependent, so a Galileo SSR orbit was decoded as GPS | Silent corruption. Now rejects non-1057/8/9. |
| P5 | **LNAV handover read out of the parity half** — parity is bits 0-5, info is 6-29 | Subframe IDs garbage. Fixed via an info-bit accessor that makes reading parity **unrepresentable**. |
| P6 | **LNAV: 18 of 23 ephemeris fields read parity bits** | Clock + orbit corrupt. **Deliberately NOT fixed** — subframe-1 offsets under-determined without IS-GPS-200 Table 20-I. `build_ephemeris_from_sfrbx` has no production caller; must not be wired up first. |
| P7 | **10 RINEX defects** — Navic `'I'` dropped (8994 of 122460 records on the real WTZR file); RINEX-2 nav hardcoded GPS; two-digit-year `80` → 2080 vs 1980; epoch flag never read; blank sat-id columns fabricated GPS PRN-0 satellites and desynchronised the file; byte-index panics on `MARKER NAME` | Silent data loss/desync. |
| P8 | **Bernese DCB validity stamps read from wrong columns** → both fell back to a hard-coded 2000-2050 window; measured `get_bias(G01,C1C,2020-12-20)` returning a value for a product valid only from 2020-12-24 | Wrong biases for 50 years of dates. |
| P9 | **`receiver_pcv` passed geodetic LLH as the ECEF origin to `az_el`** — every line of sight taken from (0.78,0,0) | 20 000 km due-north lookup: zenith 72.30° instead of 90°. |

`rtcm3/ephemeris.rs` went **0% → 100%** and had **no production caller** — dead code.
Every RTCM layout claim checked against RTKLIB `rtcm3.c`; `crc24q` differentially
validated against RTKLIB's table (6/6 golden vectors, bit-for-bit).

## gneiss-rtk composite/streaming

| # | Defect | Impact |
|---|---|---|
| C1 | **`tc_ppp` cannot process ANY epoch with <30 observation rows** (`tc_ppp.rs:306-345`) — fixed 30x30 matrices, only `m` rows written, so `S` has zero rows and `try_inverse()` returns `None`. Invisible because the existing test passes `satellites: vec![]`. | **CRITICAL.** A 6-satellite epoch is 12 rows and fails outright. |
| C2 | **`tc_ppp` hands the RTS smoother an identity transition matrix** (`tc_ppp.rs:177`) — `predict()` discards Phi, then fabricates `I`. `tc_rtk.rs:153` uses `predict_with_phi` correctly. | Every backward gain computed from Phi=I. |
| C3 | **`streaming.rs` discards the GPS week from the base-epoch buffer key** (`:86-87`) | Buffer matches across a 604 800 s separation. |
| C4 | **A single non-finite observation permanently NaNs the DD-IEKF** (`update/kalman.rs:11-28`) — neither Cholesky nor the LU fallback rejects NaN; `is_zero()` is false for NaN | Currently masked by `total_cmp` sorting NaN last in `screen_gross_pr_errors` — **one comparison away from poisoning the filter.** |

**Retraction convention resolved:** it is **LEFT / global-frame everywhere**, and
all components agree. Four tests now pin it (algebraic, physical, golden vector,
filter/smoother agreement). `AGENTS.md`'s `predictor.rs:86-91` citation is stale on
BOTH counts — the file does not exist, and `d_theta = -psi` has no counterpart
because this ESKF state is ECEF, not local-level, so there is no yaw error to negate.

## Agent discipline worth noting

- Two agents reported being **wrong** and said so unprompted (an MSM header-width
  guess, a GLONASS epoch expectation, a truncated-record premise that would have
  dropped 13146 records to 2193, and a `sqrt(26560)=162.97 km` error conflating km
  with km^½).
- 10 bug-demonstration tests ship as `#[ignore]`d **failing** tests with re-enable
  conditions — `cargo test -p gneiss-rtk --lib -- --ignored` fails all 10 today.
- Weak existing tests found: `test_real_world_rtcm3_parsing` **passes having parsed
  nothing** (its dataset does not exist); RINEX nav fixtures used space-padded
  `"G 1"` sat IDs no producer emits.

## Integration checklist (NOT DONE)

1. Wait for the BeiDou implementer.
2. `cargo test --workspace --no-fail-fast`, `cargo clippy --workspace --all-targets`.
3. Re-measure coverage.
4. **Audit the new tests the way the last batch was audited** — re-derive a sample
   of numeric claims independently. These agents are strong, but 300+ more
   agent-written tests is exactly where slop hides.
5. Decide what to do with the ~20 `#[ignore]`d failing bug-demonstration tests:
   they are valuable evidence but must not accumulate as permanent decoration.
