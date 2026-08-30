> **Superseded.** This matrix documents solver mode combinations from early engine versions. Current post-processing uses the unified multi-base network RTK pipeline. See `docs/PROJECT_STATUS.md` and `docs/TIER1_ROADMAP.md`.

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
| **Temporal** | forward-only / forward+backward (RTS smoothing) | **Default is now forward+backward** (5afbc8a). `--single-pass` opts out to forward-only; the old opt-in `--enable-backward-smoothing` flag was renamed and inverted, not kept as an alias. |
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
| 5 | SWFG | rover+base | forward-only | +GLONASS | **No longer the default** (fixed 5afbc8a) | Was the CLI's actual default with no flags: h_p50=100mm, v_p50=219mm, h_p95=**998mm**. Now only reachable via explicit `--single-pass` (with a printed warning); the bare no-flags command now gets row 9's h_p50=24mm instead |
| 6 | SWFG | rover+base | forward+backward | any | **Not a real cell** | `run_forward_pass` unconditionally upgrades to rtk_iekf whenever base+base_position are both present and no IMU — this configuration can't execute; it silently becomes row 8/9 instead |
| 7 | rtk_iekf | rover-only | (any) | (any) | **N/A by construction** | Double-difference is definitionally base-relative; there's no rover-only DD mode, this isn't a gap |
| 8 | rtk_iekf | rover+base | forward-only | GPS-only | **Works well** | This is `eval_network_ppk`'s "Forward RTK" numbers throughout the project's own benchmark history (e.g. P181 h_p50=25mm). Reachable via `PostProcessOptions{enable_bidirectional: false}` or the eval binaries — **not exposed via any `gneiss-cli` flag** |
| 9 | rtk_iekf | rover+base | forward+backward | GPS-only | **Works well — the validated path, now the default** | `gneiss-cli` with no flags (was `--enable-backward-smoothing`, now default per row 5): h_p50=24mm, p99 down from SWFG's 1.2m to 0.095m |
| 10 | rtk_iekf | rover+base | forward+backward | +GLONASS | **Fixed and reachable (Pattern-1 round), real cost confirmed worse than documented** | `enable_glonass` was previously read from a bare env var *inside* each of the forward/backward passes, nested in an unrelated baseline-length gate that differed between the two passes — so it silently never reached baselines over 25km (P222, SLAC) and could disagree between passes near the boundary. Now a top-level `PostProcessOptions` field, applied unconditionally in both passes; `gneiss-cli --glonass` exposes it. Re-measuring with the fix in place shows P225's true cost is far worse than the old (bug-corrupted) number: v_RMS 236→**796mm**, not the "+92mm p95" previously recorded. P222 now genuinely receives GLONASS (verified) but shows zero measurable change — a real, uninvestigated difference between the two baselines, not a remaining wiring bug. Still correctly kept opt-in; the case for that is now stronger, not weaker. Full detail: docs/NETWORK_RTK_NEXT_STEPS.md. |
| 11 | rtk-ins / ppp-ins (IMU-coupled) | — | — | — | **Unreachable** | Variants exist in `EngineConfig`; `--mode` can't select them and no `--imu`/IMU-loading path exists anywhere in `gneiss-cli`. Matches the roadmap's own note that frame-safety/IMU types aren't yet wired to estimator call sites |
| — | `ekf/` (uncommitted, main tree only) | — | — | — | **Not part of this branch or any entry point** | Untracked module found earlier this session; doesn't exist in this worktree; flagged separately, not a candidate for this matrix until someone decides its fate |

## What this means, ranked by how much it matters

1. ~~Row 5 is the actual product risk~~ **Fixed (5afbc8a).** The bare no-flags command now runs the validated 4-pass pipeline by default; the old fast/inaccurate path moved behind an explicit `--single-pass` opt-out with a printed warning.
2. ~~Row 2/4 is a real, reproducible bug~~ **Fixed (2647bd3).** `IfbGlonass` was being added to the SWFG factor graph from a different, looser check than the one deciding which factors get built, so it could end up with zero factors connecting it whenever GLONASS ephemeris was missing/partial — an unconditional crash, not a degraded mode. Root cause and fix in the row-2 entry above.
3. **Rows 8 and 11 are still open** capability gaps, not bugs: real, working functionality (forward-only rtk_iekf; IMU coupling) that exists at the library/eval-binary level but was never exposed as CLI options.
4. ~~Row 10 is still open~~ **Fixed and re-measured.** `--glonass` now exists and does what it says; the bug that hid the true cost at long baselines is gone, and the corrected measurement argues for keeping it opt-in even more strongly than the original (buggy) one did.

Rows 2/4, 5, and 10 closed. Rows 8 and 11 remain: real gaps, neither urgent in the way the closed ones were.
