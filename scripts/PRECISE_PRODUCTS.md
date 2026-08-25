# Multi-GNSS Precise Products — Open Source Report

Mission: free, no-auth pipeline for multi-GNSS precise orbit/clock products
(2025-06-09 / DOY 160 / GPS week 2370 used as the test vector).
Probe date: 2026-08-24.

## VERIFIED WORKING (anonymous HTTPS, downloaded + parsed end to end)

| Source | URL pattern | Products | Scope | Latency |
|---|---|---|---|---|
| **AIUB** (CODE home server) | `https://www.aiub.unibe.ch/download/CODE_MGEX/CODE/{YYYY}/COD0MGXFIN_{YYYY}{DDD}0000_01D_05M_ORB.SP3.gz` (+`_30S_CLK.CLK.gz`, `_OSB.BIA.gz`, `_ERP.gz`) | SP3-d 5 min orbits; RINEX 3 clock 30 s | GPS+GLO+GAL+BDS+QZSS | final (~9-13 d) |
| **SIO Garner** (IGS data center) | `https://garner.ucsd.edu/pub/products/{gpsweek}/{STREAM}_{YYYY}{DDD}0000_01D_05M_ORB.SP3.gz` (+`_30S_CLK.CLK.gz`) with STREAM = `GFZ0MGXRAP` (rapid), `JAX0MGXFIN`, `IAC0MGXFIN` (finals); week dir also holds GPS-only `IGS0OPS*` combined | SP3 + RINEX clk (+ OSB biases for GFZ rapid) | multi-GNSS (verified C/E/G/J/R) | rapid ~17 h (GFZ), final (JAX/IAC) |
| **BKG mirror** | `https://igs.bkg.bund.de/root_ftp/IGS/products/{gpsweek}/IGS0OPSFIN_...15M_ORB.SP3.gz` (+`_30S_CLK.CLK.gz`) | IGS combined finals/rapids | **GPS+GLO only — no Galileo** | final/rapid |

Verified content on DOY 160:
- GFZ rapid SP3: 289 epochs x 131 sats (31 GPS, 29 GAL, plus GLO/BDS/QZSS/GEO), IGb20 frame.
- JAX MADOCA final CLK: 331 k AS records covering C/E/G/J/R at 30 s.
- Both flow through `gneiss-parsers::{sp3, precise_orbit, rinex_clk}`; E02 interpolates to |r| = 29 600 km (correct MEO radius). See `crates/gneiss-parsers/examples/precise_fetch_check.rs`.

## DEAD OR AUTH-WALLED (do not rely)

- `files.igs.org/pub/product/*` → products removed; readme points to CDDIS.
- CDDIS → directory pages return HTTP 200 but every file GET returns an HTML Earthdata login page. Auth wall confirmed live.
- BKG MGEX tree: `root_ftp/MGEX/products` empty; `root_ftp/IGS/products/mgex/` frozen at GPS week 2237 (~2022).
- WHU `igs.gnsswhu.cn`, IGN `gnss-data-portal.ign.fr` / `igs.ensg.ign.fr`, GFZ `ftp.gfz.de`, Geoscience Australia `data.geo.ga.gov.au` → TCP/DNS unreachable from this network.
- ASI `www.gsc.eur.ac.it` → DNS record gone.
- ESA `navigation-office.esa.int` alive but product paths 404; distributes via CDDIS only.

## USAGE

```bash
python scripts/fetch_precise_products.py --date 2025-06-09                 # auto tier
python scripts/fetch_precise_products.py --year 2025 --doy 160 --latency rapid --neighbors 1
```

Output: `datasets/precise/{YYYY}/{DDD}/{STREAM}_{YYYYMMDD}.sp3|.clk`
(gunzipped, header-validated, cached — re-runs never re-download valid files).
Fallback chains are tried per day; exit code != 0 only if ALL sources fail.
21 unit tests in `scripts/test_fetch_precise_products.py` (network mocked).

## SELF RED TEAM

1. **Single-host concentration**: both primary multi-GNSS streams depend on one host each (aiub.unibe.ch, garner.ucsd.edu). Garner has been an IGS DC for decades (low risk); AIUB's `download.` host had intermittent DNS *during testing* — the fallback chain absorbed it (JAX MADOCA final is nearly as good for orbits; clocks differ slightly).
2. **No open combined multi-GNSS RAPID exists** (IGS0MGXRAP not operational). We take GFZ's single-AC rapid. If GFZ has a bad day, next-day PPK falls back to CODE final-quality orbits or GPS-only IGS rapid — quality cliff, flag it in logs.
3. **TLS**: garner serves a chain macOS python rejects (self-signed root); AIUB redirect target has flaky DNS. The script verifies first, falls back to unverified with a loud warning. Risk: MITM could inject bad ephemeris — mitigated by format validation + the fact that spoofed orbits produce obvious residuals; consider pinning or checksums vs. a second source if paranoid.
4. **Day-boundary**: files cover [00:00, 24:00]; interpolation near edges needs neighbor days → use `--neighbors 1`. Lagrange extrapolation outside data range explodes (verified: garbage positions when querying Sunday TOW against Monday-start data) — engine must reject out-of-range queries, not clamp.
5. **Rate limits**: none observed (garner/AIUB happily served repeated fetches). Be polite anyway: caching means N runs == 1 download per file; keep UA string so abuse is traceable. Garner may throttle abusive IPs like other .edu hosts.
6. **License**: IGS/CODE/GFZ products are openly redistributed ("open access, cite the analysis centers"). Not legal advice; attribution line recommended in engine docs: "orbits/clocks courtesy of CODE(AIUB)/GFZ via IGS".
7. **Naming churn**: IGS long-product names changed once already (igl→IGS0OPSFIN era). If downloads start failing en masse, first suspect a naming-scheme migration — re-probe the week-dir listing rather than debugging code.
8. **Clock-vs-orbit consistency**: SP3 clocks are 5-min; use the RINEX CLK (30 s) for the clock term as the engine intends; mixing SP3 orbit + RINEX clock from different ACs adds ~cm-level inconsistency (acceptable for >20 km baselines, avoid for PPP).
