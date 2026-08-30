# Roadmap to Tier-1 PPK Parity

**Purpose of this document**: Sprint 16 (`docs/PROJECT_STATUS.md`) covers *accuracy*
in depth. This document is the cross-cutting companion: performance, UX,
features, and code quality, synthesized into one prioritized plan alongside
accuracy. Where Sprint 16 already has the detail, this points there instead
of repeating it.

**Comparison target, made explicit.** The standing goal names Leica, NovAtel,
and Qinertia. Those three make both real-time RTK receivers *and* PPK
post-processing software, and they are different product categories with
different requirements. Gneiss is a post-processing engine — it reads
complete RINEX sessions and produces a trajectory, the same shape of product
as **NovAtel Waypoint (GrafNav / Inertial Explorer)**, **Leica Infinity**,
and **SBG Qinertia**. This document benchmarks against *those specifically*,
not against RTK receiver firmware. That reframing matters: it lowers the
urgency of real-time streaming (Sprint 12f — still valuable, still real, but
not a blocking requirement for this product category the way it is for a
receiver) and raises the urgency of things PPK tools are specifically judged
on: multi-format export, batch/multi-mission workflows, QC reporting, and
exposing the engine's full capability through the actual shipped tool.

Sources for the tier-1 feature claims below: [NovAtel Waypoint overview](https://novatel.com/products/waypoint-post-processing-software), [Inertial Explorer product page](https://novatel.com/products/waypoint-post-processing-software/inertial-explorer), [GrafNav product page](https://novatel.com/products/waypoint-post-processing-software/grafnav), [Qinertia product page](https://www.sbg-systems.com/software/qinertia-gnss-ins-ppk/), [Qinertia CLI examples repo](https://github.com/SBG-Systems/qinertiaExamples). Leica Infinity claims are lower-confidence (marketing pages only, not user-manual-verified) — flagged inline where used.

---

## 1. Accuracy — see Sprint 16 for the full analysis

Summary only: gneiss runs 1.5-2.6x worse on horizontal RMS and 1.9-3.2x
worse on vertical RMS than Leica's *receiver* datasheet spec at 15-50km
baselines (the closest measured comparison point available; no PPK-software
accuracy spec has been found published in comparable form — vendors publish
receiver specs, not post-processing-software specs, because the post-
processing step is assumed to reach the receiver's noise floor). Sprint
16's own diagnostic (run this session) found the gap is a **steady-state
noise floor problem**, not a wrong-fix problem: ordinary, correctly-fixed
epochs already run 10-55mm horizontal at P181's 15km baseline, straddling
Leica's 23mm spec. Next step per Sprint 16: investigate receiver/antenna
modeling depth and noise-weighting tuning, not more ambiguity-resolution
work. Not repeated further here — see Sprint 16 for evidence and detail.

---

## 2. Performance

**Measured this session** (release build, single core, `gneiss-cli process`,
full 24h/2880-epoch single-baseline RTK session, P224 rover + P181 base):
**20.2 seconds wall-clock**, 99% CPU (effectively single-threaded). That's
~143 epochs/second, a ~4300x real-time processing speedup for one baseline.
This is *fast* in absolute terms — a full day's single-baseline survey
processes in under half a minute. Tier-1 PPK tools are also fast (batch
post-processing at many-times-real-time is table stakes for the category,
not a differentiator), so raw single-baseline throughput is not believed to
be a competitive gap. Two real gaps exist underneath that headline number:

- **No parallelism, despite an embarrassingly-parallel structure.**
  `eval_network_ppk.rs`'s per-base forward and backward passes are
  independent of each other until the network-fusion step combines them
  (`for base in bases { ... }`, sequential). A 3-6-base network job
  currently takes roughly N times the single-baseline cost, sequentially,
  when the N per-base passes could run concurrently. This is very likely a
  3-6x wall-clock win on network jobs for close to zero algorithmic risk
  (no change to any estimator, purely a scheduling change) — a rare case
  where a large win doesn't require touching accuracy-sensitive code at
  all. `rayon` is the obvious tool; the per-base state (`GnssRtkIekf`,
  trajectories) already looks independent enough to parallelize directly,
  but this needs verifying against the actual mutable-state boundaries
  before assuming it's a pure drop-in.
- **No data on how this scales.** Every measurement in this project is a
  single day, at most 6 bases. Nobody has run a multi-day continuous
  mission or a network of dozens of bases (the scale a serious network-RTK
  service would eventually target) to find where memory or wall-clock
  stops scaling linearly. Not urgent today (no user need has surfaced it),
  but worth a synthetic stress test before any claim about "handles large
  jobs" ships.
- **No memory profiling has ever been done in this project.** Unknown
  whether memory footprint is trivial (a day of RINEX + Kalman state is not
  large) or has any surprises (e.g. `history: Vec<IekfSnapshot>` growing
  unbounded across the whole session for the AR-dump feature, retained for
  the full mission length rather than windowed).

**Priority for this dimension**: parallelize the per-base passes first (high
value, low risk, mechanical) before doing anything else here — there's no
reason to invest in scalability or memory work before capturing the easy win
that's already sitting in the code's own structure.

---

## 3. UX

This is where the sharpest, most consequential gap of the whole roadmap
was found this session: **the engine's best, most differentiated capability
is invisible to anyone using the actual shipped tool.**

### 3a. `gneiss-cli` cannot do network RTK at all

`gneiss-cli process` takes exactly one `--rover` and one `--base`. Every
accuracy number this project has ever measured and guarded against —
network fused fix rate 97.6%, the whole Sprint 16 analysis, all of
`check_network_benchmark.py`/`check_multignss_benchmark.py` — comes from
`eval_network_ppk`, an internal evaluation binary with a hardcoded dataset
path, not a tool an actual user can point at their own data. A user
downloading and running `gneiss-cli` today gets single-baseline RTK only,
which is the *weaker* of the two modes this project has spent most of its
recent effort validating. This is a bigger UX problem than any missing
output column: **the product's headline capability is unshippable in its
current form.**

Fixing this is plumbing, not new engineering: `eval_network_ppk.rs`'s
per-base loop, cross-base wide-lane UPD solve, and combiner-fusion logic
are all proven, tested, guarded code. The work is exposing it through
`process`'s CLI surface (`--base` needs to accept multiple paths or a
base-list file; a `--network` mode flag or auto-detection when >1 base is
given) rather than building anything new.

### 3b. The output format discards data the engine already computes

`bin/gneiss-cli/src/process.rs`'s output writer emits exactly:
`GPST-Week TOW(s) x-ecef(m) y-ecef(m) z-ecef(m) Q`. Checked
`SmoothedEpoch` (`crates/gneiss-rtk/src/post_process/combiner.rs`) — the
struct the writer iterates over already carries:

```
velocity_ecef: Option<Vector3<f64>>
attitude: Option<UnitQuaternion<f64>>
cov_position: Matrix3<f64>
std_east, std_north, std_up: f64
n_satellites: usize
```

None of this reaches the output file. Tier-1 PPK output (GrafNav's
"Flexible Export Wizard", Qinertia's export options) always includes
per-epoch uncertainty — it's the deliverable a surveyor actually needs
(a coordinate without a stated precision isn't a usable survey product).
**This is a near-zero-cost fix**: the data exists and is already computed
correctly by the smoother; it needs a few more `write!` calls, plus a call
to the already-tested `gneiss_core::coords::ecef_to_llh` for geodetic
output (surveyors work in lat/lon/height or a projected grid, essentially
never in raw ECEF). This should ship before any of the harder items below —
it's the highest value-per-hour item on this entire document.

### 3c. No QC artifact

Tier-1 tools produce a report (Leica Infinity: "detailed interactive
analysis charts and processing reports"; GrafNav 10.0 added a satellite map
view specifically for visual QC). Gneiss's only quality signal is one
console log line per run (`fix rate: 87.5%, median sep: 0.075m`). A
surveyor needs to know, per session: fix-rate over time, satellite count
and DOP history, which epochs were flagged/downgraded and why, and a
plain-language pass/fail against expected tolerances. This doesn't need a
GUI to start — a structured CSV/JSON summary alongside the trajectory file
would already close most of the gap, with a rendered chart/report as a
later, separate step.

### 3d. No batch or multi-mission processing

One rover file per invocation, no directory-of-missions mode. A production
survey workflow processes many sessions per day; requiring N manual
invocations (or a user-written wrapper script) for N sessions is exactly
the kind of friction tier-1 tools eliminate with project/batch files.

### 3e. Config ergonomics

`--config <CONFIG>` takes a raw JSON file with no visible schema, no
`--print-default-config` to get a starting point, and no validation errors
beyond whatever `serde_json` produces on a malformed file. Not audited in
detail this pass (lower priority than 3a-3d) — worth a dedicated look at
actual error messages a user would see on common mistakes (missing field,
wrong type, unknown key).

### Priority within UX, in order

1. Extend the output writer (3b) — hours of work, the data already exists.
2. Expose network processing through `gneiss-cli` (3a) — the single
   highest-leverage change on this whole document; unlocks the tool's own
   best capability for real users, largely plumbing.
3. A structured QC summary file (3c, CSV/JSON first, rendering later).
4. Batch/multi-mission processing (3d).
5. Config UX polish (3e) — real, but lower-impact than the above four.

---

## 4. Features

**What gneiss already has, at genuine tier-1-comparable depth**: multi-
constellation (GPS/GLONASS/Galileo/BeiDou/QZSS/SBAS), multi-frequency,
RTK and PPK (forward + backward + combiner smoothing — the same shape as
Qinertia's forward/backward/smoothed pipeline), a real INS/IMU tight-
coupling architecture (`swfg`'s factor-graph engine, preintegration, bias
estimation — 14 files' worth, not a stub), network multi-base fusion
*at the engine level* (see UX 3a for why this doesn't reach users yet),
precise ephemeris (SP3+CLK+satellite PCO), receiver antenna PCV, ocean
tide loading, ionosphere modeling (Klobuchar plus per-satellite estimated
states), and LAMBDA-based ambiguity resolution with partial/wide-lane
cascades. This is a substantial, non-toy feature set — the comparison
below is about the remaining distance, not starting from zero.

**Gaps, prioritized by how load-bearing they are for the PPK-software
comparison specifically:**

| Gap | Tier-1 reference | Why it matters | Effort (rough) |
|---|---|---|---|
| Network fusion not exposed via CLI | — (this is a gneiss-specific gap, not a feature tier-1 lacks) | See UX 3a; blocks every other comparison from mattering | Low (plumbing) |
| Multi-format export (KML, SBET, CSV w/ lat-lon+sigma, DXF) | GrafNav's Flexible Export Wizard names all of these explicitly | PPK output feeds downstream tools (GIS, lidar/photogrammetry georeferencing); a single custom `.pos` format is a real interoperability wall | Medium |
| Orthometric height output | Universal in surveying-grade tools | `gneiss-geodesy::geoid` already exists, complete, tested, zero callers (Sprint 13/16 finding) — this is a wiring task, not new engineering | Low |
| QC report / visual map view | Leica Infinity's charts; GrafNav 10.0's satellite map view | See UX 3c | Medium |
| Lever-arm / IMU mounting auto-estimation | Qinertia's `mechanicalEstimationPass`; Inertial Explorer's lever-arm estimation | Only matters once INS integration is a first-class delivered feature, not just internal `swfg` infrastructure — sequence after UX 3a-3c | Medium-High |
| PPP-AR (base-station-free ambiguity-resolved PPP) | Both GrafNav and Inertial Explorer support PPP as a first-class mode alongside differential | Gneiss has PPP-*adjacent* infrastructure (the legacy SWFG/`spp.rs` path, per-satellite iono states) but PPP-AR specifically (fixing PPP ambiguities without a base) hasn't been confirmed as implemented or measured this session — needs its own investigation before claiming a gap size | Unknown — investigate first |
| Network *adjustment* (least-squares multi-station network solve, distinct from RTK fusion) | Leica Infinity markets this explicitly | Different problem from network RTK fusion (that's real-time-style Kalman fusion of redundant baselines to one rover; network adjustment is a static geodetic-network least-squares solve, more relevant to control-point surveying than rover positioning) — likely genuinely out of scope for what gneiss is trying to be; noted for completeness, not recommended as a priority | High, and possibly wrong target |
| Real-time / streaming mode | All three vendors' receivers; GrafNav/Inertial Explorer also added real-time preview modes | Sprint 12f — still real, still valuable (even PPK tools increasingly offer real-time preview), but not the *blocking* gap for a post-processing-software comparison the way Sprint 16 framed it for a receiver comparison. Re-prioritized down from "largest item on the entire roadmap" to "large, real, but sequenced after the UX fixes above" in this framing. | Very High |

---

## 5. Code Quality

Current measured state (this session, `network-rtk-long-baseline` @ `ff480cd`):

- **0 compiler warnings, 0 clippy warnings** (`cargo build --workspace`,
  `cargo clippy --workspace`) — held continuously, verified dozens of times
  this session across every change.
- **`unwrap_used = "deny"` is a workspace-level clippy lint**, not just a
  documented convention — this means "no unwrap() in production" is
  *mechanically enforced*, not aspirational. Two explicit, narrow
  exceptions: `eval_swfg.rs` and `swfg/benchmark.rs`, both internal
  eval/benchmark binaries rather than shipped code, both carrying their own
  `#![allow(clippy::unwrap_used)]`.
- **21 files remain over the 500-line limit** (list: `spp.rs` 2228,
  `rtcm3/msm.rs` 1434, `ephemeris.rs` 1398, `rinex/obs.rs` 1362,
  `rtk_iekf/mod.rs` 1203, `rinex/nav.rs` 1006, `atmosphere.rs` 1005,
  `ubx.rs` 936, `rtk_iekf/update.rs` 814, `frames.rs` 699,
  `receiver_antenna.rs` 672, `eval_network_ppk.rs` 661, `ambiguity/lambda.rs`
  635, `swfg/solver.rs` 632, `rtk_iekf/formation.rs` 618,
  `swfg/engine/mod.rs` 614, `frequencies.rs` 595, `rtk_iekf/state.rs` 581,
  `receiver_pcv.rs` 549, `swfg/imu_preintegration.rs` 539,
  `swfg/pipeline/factors.rs` 536). Several of these are large mostly
  *because* of thorough co-located test coverage (CLAUDE.md's own
  convention), not unreviewable production logic — see Sprint 13's
  `rinex.rs` split entry for the reasoning on reading this number
  honestly. `spp.rs`, `formation.rs`'s nesting, and the rest are tracked
  and deliberately not rushed (see Sprint 13).
- **Mutation testing is broken and unverified.** `cargo-mutants` is
  installed and `.cargo/mutants.toml` exists, but invoking it in this
  environment hits a broken recursive cargo alias (`alias mutants has
  unresolvable recursive definition: mutants -> mutants`) shadowing the
  real binary, and an archived doc (`docs/archive/SPRINT_ROADMAP.md`)
  separately records "42% -> 0 survivors, deferred (cargo-mutants v27
  compat)" from an earlier attempt. Net result: CLAUDE.md's own hard rule
  ("Mutation testing: 0 survivors") has **never been verified against the
  current codebase**. This is a real, honest gap, not a rounding error —
  a codebase can have 100% line coverage and still have logic bugs that
  only mutation testing would catch (this exact failure mode is why
  CLAUDE.md has the rule at all). Fixing the tooling (likely: resolve the
  alias shadowing first — that may be the whole problem, not just the
  cargo-mutants version — then re-attempt a run, possibly scoped to one
  crate at a time given the codebase's size) should happen before any
  further code-quality claims lean on "we test thoroughly," since
  thorough-looking coverage and mutation-proof coverage are different
  claims and only one of them has ever been checked here.
- **Nesting depth**: three severe (6-9 level) violations found and fixed
  this session (`rinex/obs.rs`, `ionex.rs`, `rtcm3/msm.rs`); one known,
  real depth-5 case in `rtk_iekf/formation.rs` deliberately deferred
  (core DD-formation math, higher stakes than the parser fixes) — see
  Sprint 13.
- **Duplicate-code audit**: found and removed 3 confirmed dead-duplicate
  modules this session (`HatchFilter` x2, `HelmertParams` x2 plus the
  `antex` parser x2 in the orphaned `gneiss-geodesy` crate) via a
  `pub struct`/`pub fn` name-collision grep across the workspace. That
  sweep covered public-item name collisions specifically; it has not been
  extended to look for duplicated *logic* under different names, which a
  simple name grep can't find.

**Priority for this dimension**: fix the mutation-testing tooling first —
it's the one hard rule in CLAUDE.md that has literally never been checked,
which makes every other code-quality claim in this document (and every
prior sprint) slightly less certain than it sounds. Everything else here is
incremental, known, and already tracked in Sprint 13.

---

## 6. Prioritized synthesis across all five dimensions

Combining the per-dimension priorities into one sequence, weighted toward
low-cost/high-leverage items first and toward not touching accuracy-
sensitive code without a specific reason:

1. **Extend the CLI output writer** (UX 3b) — hours, zero new engineering,
   the data already exists.
2. **Expose network fusion through `gneiss-cli`** (UX 3a) — the single
   highest-leverage item overall; makes every accuracy number this project
   has measured actually reachable by a real user.
3. **Fix the mutation-testing tooling and run it** (Code Quality) —
   foundational; every other quality claim is provisional until this
   happens at least once.
4. **Parallelize per-base passes** (Performance) — mechanical, low risk,
   real wall-clock win on exactly the network jobs item 2 just made
   reachable.
5. **Continue Sprint 16's accuracy investigation** (noise-floor / receiver-
   antenna modeling depth) — open-ended effort, but the highest-value
   *accuracy* lever currently identified.
6. **Wire `gneiss-geodesy::geoid`** for orthometric height output
   (Features) — low effort, closes a basic parity gap, independent of
   everything else.
7. **A structured QC summary artifact** (UX 3c) — medium effort, real
   value, sequence after the CLI/output work above since it should report
   on the network-mode runs too.
8. **Investigate PPP-AR status** (Features) — find out whether this is a
   real gap or already-covered-and-unmeasured before sizing any work here.
9. **Real-time/streaming architecture** (Sprint 12f) — still large, still
   eventually necessary for RTK-receiver-class parity, but sequenced after
   the above given the PPK-post-processing-software framing this document
   uses. Re-evaluate priority if the product direction shifts toward
   real-time / receiver-adjacent use cases.
10. **Remaining file-size/nesting debt, export-format breadth, batch
    processing, lever-arm estimation** — ongoing background work, already
    tracked, no single item here blocks anything else.

## What this document can't close

Repeated from Sprint 16 because it's still true and still the most
important caveat on the whole plan: certified firmware, years of field
validation across real-world conditions, 24/7 support organizations, and
(for some markets) safety certifications are not things a roadmap of
engineering sprints produces. What sections 1-6 above target is the actual
closeable gap — algorithmic, architectural, and feature-completeness — not
a claim that finishing this document's list equals being Leica.
