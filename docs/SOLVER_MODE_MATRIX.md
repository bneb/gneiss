# Solver Mode Matrix — What Actually Works

Every combination below was run against real CORS RINEX data through the
actual entry point named (not assumed from reading code). Where a cell
says BROKEN or UNREACHABLE, that's a reproduced result, not a guess.

## The independent axes

The CLI conflates these in a way that isn't obvious from `--help`, which
is half the reason this matrix exists:

| Axis | Values | Selected by |
|---|---|---|
| **Engine** | SWFG (Era 2) / rtk_iekf (Era 3, DD-RTK) | **Not a flag.** Chosen automatically: `--base` + `--base-position` (and no `--imu`, which doesn't exist as a flag) routes to rtk_iekf; anything else routes to SWFG. `--mode` does NOT select the engine. |
| **Base presence** | rover-only / rover+base | `--base` given or not |
| **Temporal** | forward-only / forward+backward (RTS smoothing) | `--enable-backward-smoothing` |
| **`--mode`** | spp / ppp / rtk | Only tunes SWFG's internal parameters (elevation mask, pipeline stages, initial-position source) for whichever engine got picked above. Does not change dispatch. `rtk-ins`/`ppp-ins` exist in `EngineConfig` but `--mode` can't reach them (string match only covers spp/ppp/rtk) — moot anyway since nothing loads IMU data (no `--imu` flag exists). |
| **Constellation** | GPS-only / +GLONASS / +Galileo | `--systems` (pre-filters observation data) for both engines; rtk_iekf *additionally* gates GLONASS at the DD-formation level via `enable_glonass`, settable only through the `GNEISS_GLONASS` env var — `--systems` alone cannot turn this on for rtk_iekf. |

## The matrix

Baseline data: P224 rover, real CORS RINEX (`datasets/cors_short_baseline/`), GLONASS+Galileo present unless noted.

| # | Engine | Base | Temporal | Constellation | Result | Evidence |
|---|---|---|---|---|---|---|
| 1 | SWFG | rover-only | forward-only | GPS-only (`--systems G`) | **Works** | 5/5 epochs, n_sat=7, no error |
| 2 | SWFG | rover-only | forward-only | +GLONASS (default) | **FIXED** (was BROKEN) | Was `OrphanVariable("IfbGlonass")` on every epoch, 0/2880 processed, whenever a GLONASS satellite is tracked with no matching ephemeris (e.g. GPS-only nav file). Fixed 2647bd3: `ensure_ifb_glonass()` now creates the variable lazily at the exact point a factor references it, instead of pre-emptively from raw satellite tracking. Real GLONASS ephemerides still get real IFB factors (verified against `multignss_2025d160`, unaffected) |
| 3 | SWFG | rover-only | forward+backward | GPS-only | **Works** | Same as #1 — `execute_post_process` falls back to the same SWFG call when there's no base |
| 4 | SWFG | rover-only | forward+backward | +GLONASS | **FIXED** (was BROKEN) | Same root cause and fix as #2 (identical underlying `SwfgEngine::process_epoch` call) |
| 5 | SWFG | rover+base | forward-only | +GLONASS | **Works, but bad accuracy** | This is the CLI's actual **default** with no flags: h_p50=100mm, v_p50=219mm, h_p95=**998mm**, on data where the properly-wired path (row 8) gets h_p50=24mm |
| 6 | SWFG | rover+base | forward+backward | any | **Not a real cell** | `run_forward_pass` unconditionally upgrades to rtk_iekf whenever base+base_position are both present and no IMU — this configuration can't execute; it silently becomes row 8/9 instead |
| 7 | rtk_iekf | rover-only | (any) | (any) | **N/A by construction** | Double-difference is definitionally base-relative; there's no rover-only DD mode, this isn't a gap |
| 8 | rtk_iekf | rover+base | forward-only | GPS-only | **Works well** | This is `eval_network_ppk`'s "Forward RTK" numbers throughout the project's own benchmark history (e.g. P181 h_p50=25mm). Reachable via `PostProcessOptions{enable_bidirectional: false}` or the eval binaries — **not exposed via any `gneiss-cli` flag** |
| 9 | rtk_iekf | rover+base | forward+backward | GPS-only | **Works well — the validated path** | `gneiss-cli --enable-backward-smoothing`: h_p50=24mm, p99 down from SWFG's 1.2m to 0.095m (this session's earlier fix) |
| 10 | rtk_iekf | rover+base | forward+backward | +GLONASS | **Works, small measured cost, unreachable from gneiss-cli** | `GNEISS_GLONASS=1` on `eval_network_ppk`: P225 fix rate 73.4%→70.2%, no crash (code ICBs leak, phase mostly absorbs it — already documented, deliberately kept opt-in). `gneiss-cli`'s `--systems GER` lets GLONASS *observations* through but `enable_glonass` still defaults false and only responds to the `GNEISS_GLONASS` env var, which nothing in `gneiss-cli` sets — **there is currently no way to opt into GLONASS participation from the CLI at all**, not even for experimentation |
| 11 | rtk-ins / ppp-ins (IMU-coupled) | — | — | — | **Unreachable** | Variants exist in `EngineConfig`; `--mode` can't select them and no `--imu`/IMU-loading path exists anywhere in `gneiss-cli`. Matches the roadmap's own note that frame-safety/IMU types aren't yet wired to estimator call sites |
| — | `ekf/` (uncommitted, main tree only) | — | — | — | **Not part of this branch or any entry point** | Untracked module found earlier this session; doesn't exist in this worktree; flagged separately, not a candidate for this matrix until someone decides its fate |

## What this means, ranked by how much it matters

1. **Row 5 is the actual product risk, still open.** A brand-new user runs `gneiss-cli process -r rover.obs -b base.obs -n nav.rnx -o out.pos` — the single most obvious command — and silently gets ~1m-class accuracy from a tool capable of 24mm, with zero indication a 5-10x-better mode exists behind an undiscoverable flag.
2. ~~Row 2/4 is a real, reproducible bug~~ **Fixed (2647bd3).** `IfbGlonass` was being added to the SWFG factor graph from a different, looser check than the one deciding which factors get built, so it could end up with zero factors connecting it whenever GLONASS ephemeris was missing/partial — an unconditional crash, not a degraded mode. Root cause and fix in the row-2 entry above.
3. **Rows 8 and 11 are still open** capability gaps, not bugs: real, working functionality (forward-only rtk_iekf; IMU coupling) that exists at the library/eval-binary level but was never exposed as CLI options.
4. **Row 10 is still open, and isn't a bug** — the current behavior (silently exclude GLONASS from AR) is the *correct*, evidence-backed choice — but it's also not a choice a `gneiss-cli` user can currently override even deliberately.

Row 2/4 closed; the rest of this document is still the map for what's left.
