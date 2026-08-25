# RTKLIB 2.4.3 b34 — Satellite Position/Clock Computation & Ionosphere Handling
### Reverse-engineering notes for building a competing GNSS engine

Source audited: `/tmp/rtklib_src/src/` (read-only). Primary files: `preceph.c` (634 lines), `ephemeris.c` (767 lines), `rtkcmn.c`, `pntpos.c`, plus `ionex.c`, `sbas.c`, `qzslex.c`, `rtkpos.c` for cross-references.

**Structural note (important):** in b34 the *broadcast* ephemeris functions the task asked about under preceph.c — `eph2pos`, `eph2clk`, `geph2pos`, `seph2pos`, `satpos`, `satposs` — live in **ephemeris.c** (they were moved there in 2010, see ephemeris.c header history line 30). `preceph.c` owns SP3 precise ephemeris I/O + interpolation, satellite antenna offset, and DCB files. There is **no** `readephdeps()` in this version; `readsp3()` is at preceph.c:253. Also: **`satrange()` and `satazon()` do not exist anywhere in this codebase** (verified by grep) — their presumed functionality is decomposed into `geodist` + `satazel` + `prange/ionocorr/tropcorr` (pntpos.c) and `satantoff` (preceph.c).

---

## 1. Broadcast ephemeris → position & clock (`ephemeris.c`)

### 1.1 Physical constants used (ephemeris.c:63-88)
| Constant | Value | Line |
|---|---|---|
| RE_GLO | 6378136.0 m | 63 |
| MU_GPS | 3.9860050E14 | 64 |
| MU_GLO | 3.9860044E14 | 65 |
| MU_GAL / MU_CMP | 3.986004418E14 | 66-67 |
| J2_GLO | 1.0826257E-3 | 68 |
| OMGE (GPS) | 7.2921151467E-5 rad/s (rtklib.h:64) | — |
| OMGE_GAL | 7.2921151467E-5 | 71 |
| OMGE_CMP | 7.292115E-5 | 72 |
| RTOL_KEPLER | 1E-14 relative tolerance | 79 |
| MAX_ITER_KEPLER | 30 | 88 |
| TSTEP (GLO RK4 step) | 60.0 s | 78 |
| ERREPH_GLO | 5.0 m (→ var = 25 m²) | 77 |
| STD_BRDCCLK | 30.0 m (fallback clock sigma) | 86 |
| SIN_5 / COS_5 | −0.0871557427476582 / 0.9961946980917456 (BDS GEO tilt) | 74-75 |

### 1.2 `eph2pos()` — GPS/Galileo/QZS/BDS-MEO Keplerian (ephemeris.c:181-250)
- Guard: `A<=0` → zeros out everything (190).
- `tk = time − toe`; μ and ωe selected per system via `satsys` switch (196-200) — note BDS uses its own slightly different ωe.
- Mean motion perturbed by `deln`: `M = M0 + (√(μ/A³)+deln)·tk` (201).
- **Kepler solved by Newton–Raphson**, not the usual fixed-point iteration: `E -= (E−e·sinE−M)/(1−e·cosE)`, iterate until `|ΔE|<1E-14` or 30 iters; overflow → trace error + zeroed output (203-209). (History note line 45: switched from fixed-point to Newton in 2013.)
- Argument of latitude/radius/inclination harmonic terms cus,cuc,crs,crc,cis,cic applied at 2u (217-220).
- Right ascension → ECEF: `Ω = OMG0 + (OMGd−ωe)·tk − ωe·toes` (236). The `−ωe·toes` term anchors the frame at ephemeris epoch.
- **BeiDou GEO special case (PRN ≤ 5)** (223-234): computes orbital plane coords `(xg,yg,zg)` with `Ω = OMG0+OMGd·tk−ωe·toes`, then rotates about the Y axis by **−5°** (SIN_5/COS_5 constants) combined with the `ωe·tk` rotation:
  ```
  rs[0]= xg·cos(ωe·tk) + yg·sin(ωe·tk)·cos5 + zg·sin(ωe·tk)·sin5
  rs[1]=−xg·sin(ωe·tk) + yg·cos(ωe·tk)·cos5 + zg·cos(ωe·tk)·sin5
  rs[2]=−yg·sin5 + zg·cos5
  ```
- Clock: polynomial about **toc** (`tk=timediff(time,toc)` recomputed, line 242): `dts = af0 + af1·tk + af2·tk²`, then
  - **Relativistic correction** (246): `dts -= 2√(μA)·e·sinE/c²` — the F·e·√A form with sinE from the Kepler solution.
  - **TGD is deliberately NOT included** (doc lines 177-179); applied later at measurement level (§4.3).
- Variance: `var_uraeph(sva)` (91-98): URA index table `{2.4,3.4,4.85,6.85,9.65,13.65,24,48,96,192,384,768,1536,3072,6144}` m, squared; index <0 or >15 → 6144².

### 1.3 `eph2clk()` (154-167)
Fixed-point iteration ×2: `t -= f0+f1·t+f2·t²`, returns the same polynomial. Used **only** by `satposs` to refine transmission time — no relativity here (correct, since relativity is a function of orbit geometry added in eph2pos).

### 1.4 GLONASS: `geph2pos()` (302-335) + `glorbit()`/`deq()` (252-280)
Not Keplerian — **RK4 numerical integration of the GLONASS ICD equations in ECEF (PZ-90)**:
- State vector {x,y,z,vx,vy,vz} initialized from broadcast pos/vel at toe; integrated from `t = time−toe` in fixed steps of **TSTEP = 60 s** (sign flips for backward integration; last partial step truncated; loop until |t|≤1e-9) (328-331).
- Acceleration model `deq()` (252-268):
  - Central term −μ/r³ with MU_GLO;
  - **J2 zonal harmonic** term `a = 1.5·J2_GLO·MU_GLO·RE_GLO²/r⁵` with the `(1−5z²/r²)` factor (z-equation gets `c−2a`) — i.e., full oblateness perturbation;
  - **Centrifugal** `+ωe²·x,y`;
  - **Coriolis** `±2ωe·v` coupling x↔y;
  - **Luni-solar acceleration taken directly from the broadcast message field `geph.acc[3]`** — RTKLIB does *not* compute Sun/Moon gravity itself for GLO (comment cites ICD ref [2] A.3.1.2 *"with bug fix for xdot[4],xdot[5]"* — i.e., they corrected a known erratum in the ICD's printed equations).
- Clock `geph2clk` (288-301): `τn` stored positive, returned negated per ICD sign convention: iterates twice `t -= (−taun + gamn·t)`, returns `−taun + gamn·t`. **No relativistic correction** (absorbed into the GLONASS broadcast values by design).
- Variance hard-coded `SQR(ERREPH_GLO=5.0)`.

### 1.5 SBAS/GEO: `seph2pos()` (367-383)
Trivial second-order Taylor: `r = pos + vel·t + acc·t²/2`, clock `af0 + af1·t`, **no relativity correction**, variance from the GPS URA table. This serves both WAAS/EGNOS GEOs and (via LEX path separately) QZS. Note QZSS geostationary satellites are handled through the ordinary `eph_t` Keplerian path, not seph.

### 1.6 Ephemeris selection (seleph 385-454)
- Validity window on |toe−teph|: `MAXDTOE` values (+1 s slack) — rtklib.h:200-205:

| System | MAXDTOE |
|---|---|
| GPS / QZSS | 7200 s |
| Galileo | 10800 s |
| BeiDou | 21600 s |
| GLONASS | 1800 s |
| SBAS | 360 s |

- Among candidates inside the window, picks the one with **toe closest to teph** (not newest!) unless an exact `iode` match is requested (used by SBAS/SSR paths). No health-bit filtering during selection — health handled downstream (`svh`, `satexclude`).

---

## 2. Dispatch and multi-source handling

### 2.1 `satpos()` (ephemeris.c:673-693)
Switch on `ephopt`: BRDC→ephpos; SBAS→satpos_sbas; SSRAPC/SSRCOM→satpos_ssr (opt 0/1 selects APC vs CoM SSR convention); PREC→peph2pos(time,sat,nav,**opt=1**) (always antenna phase center); LEX→lexeph2pos. Any failure sets `*svh=-1` ("correction not available") and returns 0.

### 2.2 `satposs()` (ephemeris.c:718-767)
Per observation:
1. Grab the first nonzero pseudorange across NFREQ (any frequency works).
2. **Transmission time**: `time = obs.time − P/c` (receiver-clock-included receive tag), then **subtract satellite clock bias computed from BROADCAST ephemeris** regardless of the chosen ephopt (`ephclk()`, line 743-747) — a two-stage bootstrapping: you need *some* clock to compute positions even when using precise products.
3. Call `satpos()` at that transmission time.
4. **Fallback chain** (755-760): if the returned `dts==0.0` (no precise clock available), re-fetch the broadcast clock via `ephclk`, zero the drift output, and inflate variance to `SQR(STD_BRDCCLK=30 m)`.

**GPS–GLO time conversion:** there is *no* runtime `glo2gps()` transform in satpos/satposs. Instead:
- At RINEX decode time GLONASS ephem tags are converted UTC→GPST: rinex.c:1128-1129 `geph->toe=utc2gpst(toc); geph->tof=utc2gpst(tof)` (tod/tow parsed as UTC per v2/v3 conventions).
- SP3 files tagged "UTC" are likewise converted (preceph.c:137 `if (!strcmp(tsys,"UTC")) time=utc2gpst(time)`).
- The residual GLONASS−GPS system time offset is **estimated as an extra receiver clock state**: SPP estimates x[4] (pntpos.c rescode:256, estpos:350 `sol->dtr[1] /* glo-gps time offset */`), plus x[5]/x[6] GAL/BDS offsets; missing systems get pseudo-measurement constraint rows with variance 0.01 to keep the normal matrix invertible (rescode:269-275). RTK outputs `dtr[1]` as "receiver glonass-gps time offset" (rtkpos.c:1727).

### 2.3 SSR (RTK-SSR/PPP-SSR): `satpos_ssr()` (560-656)
- Requires both orbit & clock corrections present and `iod[0]==iod[1]` (consistency check, 581).
- Age limits: **MAXAGESSR = 90 s** orbit/clock (history: was 70 before Nov 2013), high-rate clock **MAXAGESSR_HRCLK = 10 s** (85).
- If update interval ≥1 s, evaluation time shifted to interval midpoint (`t -= udi/2`, 598-599).
- Orbit correction extrapolated linearly along radial/along/cross unit vectors built from the *broadcast* position/velocity; clock correction quadratic poly `dclk0+dclk1·t+dclk2·t²` plus optional HR clock.
- Sanity limits: `norm(deph)>10 m` (MAXECORSSR) or `|dclk|>1e-6·c ≈ 299.79 m` (MAXCCORSSR) rejected.
- **Clock subtlety** (619-628, 646-647): for GPS/GAL/QZS/CMP the satellite clock is *recomputed from the broadcast polynomial* (not the relativistic value returned by ephpos!), then relativity is re-applied as `-2·r·v/c²`, then SSR dclk/c added — comment cites RTCM SSR eq. 3.12-7: `t_corr = t_sv − (dts_brdc + dclk/c)`.
- Antenna offset optionally applied (opt distinguishes SSR APC vs CoM products).
- Variance `var_urassr(ura)` (100-107): `std = (3^((ura>>3)&7)·(1+(ura&7)/4)−1)·1E-3 m`; ura≥63 → 5.4665 m; ura≤0 → default 0.15 m.

---

## 3. Precise ephemeris (`preceph.c`)

### 3.1 Constants (48-51)
```
NMAX       10      /* Lagrange interpolation ORDER (11 nodes!) */
MAXDTE     900.0   /* max |Δt| outside precise eph span (s) */
EXTERR_CLK 1E-3    /* clock extrapolation error growth (m/s) */
EXTERR_EPH 5E-7    /* orbit extrapolation error growth (m/s^2) */
```

### 3.2 `readsp3()` / SP3 ingestion (253-294 + helpers)
- Wildcard expansion up to MAXEXFILE=1024 files; extensions `.sp3/.SP3/.eph*/.EPH*`.
- Header: 22 lines; time-system flag read (line 91); if "UTC" → convert each epoch to GPST (137).
- Record scaling (179-195): position km→m (`×1000`), clock µs→s (`×1E-6`); velocity dm/s→m/s (`×0.1`), rate `×1E-10`; formal std-dev decoded via base-exponent: `pow(bfact,std_exp)` then `×1E-3` (pos m) / `×1E-12` (clock s); velocities `×1E-7/1E-16`.
- Sentinel `999999.999999` treated as invalid; predicted-vs-observed flags read from columns 75 ('P' clock pred) / 79 (orbit pred) with `opt` bitmask filter: 1 = observed only, 2 = predicted only, 4 = don't merge duplicates.
- PRN mapping: `'G'/space`→GPS, R/E/J/C/L codes; **QZS prn += 192, SBAS prn += 100**.
- `combpeph()` (211-238): sorts by time (ties broken by file index), merges records within 1 ns of each other, filling only slots that are still zero — allows combining multiple files covering different satellites of the same epochs.
- Epoch kept only if ≥1 valid position record existed (`v` flag).

### 3.3 `pephpos()` — THE clever function (395-481)
- Binary search over sorted peph array for bracketing epoch; reject requests >900 s outside data span (405-407).
- **Orbit: order-10 Lagrange interpolation over NMAX+1 = 11 nodes** via Neville's algorithm (`interppol`, 383-393). Window start centered: `i = index − (NMAX+1)/2`, clamped into `[0, ne−NMAX−1]` (419-420). With standard 15-min SP3 spacing the window spans ±75 min around the request.
- **Earth-rotation pre-correction (the classic trick)**, lines 435-441 ("correciton for earh rotation ver.2.4.0"): every node is rotated about Z by `ωe·Δt_j` (Δt_j = node_time − request_time) *before* interpolation:
  ```
  sinl=sin(OMGE·t[j]); cosl=cos(OMGE·t[j]);
  p[0][j]= cosl·x − sinl·y;  p[1][j]= sinl·x + cosl·y;  p[2][j]=z;
  ```
  Rationale: raw ECEF coordinates of the same satellite differ between epochs mostly by Earth rotation (~15 m/s ground-track drift in ECEF); interpolating them directly produces large errors because the trajectory in ECEF is not smooth over long windows. Rotating all nodes into the ECEF frame of the *request* epoch makes the interpoland smooth again. A competitor who interpolates raw SP3 xyz loses ~tens of meters.
- All 11 nodes must be non-zero for that sat, else "prec ephem outage" → failure (424-427).
- Orbit std: norm of the xyz formal sigmas of the nearest node, inflated **quadratically** when extrapolating beyond the window ends: `std += EXTERR_EPH·Δt²/2` (450-453).
- **Clock: plain LINEAR interpolation** between the two bracketing samples (455-479):
  `dts = (c[1]·t[0] − c[0]·t[1])/(t[0]−t[1])`; outside bracketing hold endpoint value; sigma inflated linearly `EXTERR_CLK·|Δt from nearer epoch|`.

### 3.4 About "16.7 ms quantization" / clock jitter handling
There is **no dedicated quantization/jitter handler** in b34 (nothing named or commented that way). What exists instead:
- SP3 clocks are stored at ns resolution and interpolated **linearly** between epochs (30 s clocks → 5 min orbits files; 15 min orbits-only files use the pos[][3] slot).
- The sub-sample curvature that linear interpolation misses is partially reconstructed *afterwards* by the velocity-based relativistic term in `peph2pos` (below) — since the periodic part of the satellite clock IS the relativistic term `−2r·v/c²` (±23 ns peak for GPS), adding it back using interpolated r,v restores dynamics that linear clock interpolation cannot.
- Residual jitter is absorbed stochastically: EXTERR_CLK sigma growth + PPP measurement noise models, not deterministic modeling.

### 3.5 `pephclk()` (483-529)
Same linear scheme against a separate CLK file (`nav->pclk`). Deliberate asymmetry: if `nav->nc<2` or request outside span → **returns 1 (success)** leaving dts untouched, so the clock embedded in the SP3 position file is used; but a mid-span outage (both neighbors zero) → returns 0 (hard failure).

### 3.6 `satantoff()` — satellite PCO (541-578)
- Builds satellite-fixed axes without yaw telemetry: `ez = −r̂` (radial/nadir), sun unit vector es, `ey = ez×es` normalized, `ex = ey×ez` — a **sun-pointing-y approximation of nominal yaw steering** (no noon/midnight maneuver modeling — known RTKLIB gap).
- Frequency pair: j=L1(idx0), k=L2(idx1), except **k=2 (L5/E5a) for GALILEO and SBAS** when NFREQ≥3 (564).
- Offsets combined into the **iono-free LC**: `γ=(λk/λj)²; C1=γ/(γ−1); C2=−1/(γ−1); dant = C1·PCO_L1 + C2·PCO_Lk` (568-577) so the result is consistent with IF-combined observations.
- PCVs come from ANTEX via `readsap()` (303-319). Nadir-dependent PCV is NOT applied here (only offsets); `antmodel_s()` (rtkcmn.c:3529) applies nadir PCV elsewhere (ppp.c), interpolated on a 5° grid (`interpvar`, 3487-3493, 19 bins 0-90°).

### 3.7 `peph2pos()` (597-634)
1. Position+clock at t (pephpos then pephclk override).
2. Position+clock again at **t + 1 ms** (`tt=1E-3`); velocity by finite difference `rs[3..5]=(rst−rss)/tt`.
3. `satantoff` added if opt=1.
4. **Relativistic clock correction applied post-interpolation** (623-627):
   ```
   dts[0] = dtss[0] − 2·dot(r,v)/c²      // periodic relativistic signature restored
   dts[1] = (dtst−dtss)/tt               // drift by finite difference
   ```
5. Output variance = orbit var + clock var.

### 3.8 DCB files (`readdcb`, 321-381)
Legacy ASCII DCB (P1-P2, P1-C1, P2-C2 sections), columns 26-34, ns→m via `×1E-9·CLIGHT`, stored `nav->cbias[sat][0..2]`. These feed `prange()` (§4.3).

---

## 4. Measurement-side functions (`rtkcmn.c` + `pntpos.c`)

### 4.1 `geodist()` — Sagnac/earth rotation (rtkcmn.c:3199-3209)
```c
if (norm(rs)<RE_WGS84) return -1;          // sanity: sat above earth radius
e = (rs−rr)/|rs−rr|;                        // LOS vector WITHOUT Sagnac rotation
return r + OMGE·(rs_x·rr_y − rs_y·rr_x)/CLIGHT;
```
Sign convention: the Sagnac term is **added** to the geometric range (positive for typical east-of-sat geometries). Equivalent to rotating receiver coords into transmission-time frame. The unit vector e is *not* corrected for the rotation (small effect, ~µrad).

### 4.2 `satazel()` (3218-3230)
ENU conversion; `az = atan2(E,N)` wrapped to [0,2π); `el = asin(U)`; degenerate horizontal component (<1e-12) → az=0; guard `pos[2] > −RE_WGS84`. Returns elevation.

### 4.3 TGD / code-bias application — `prange()` + `gettgd()` (pntpos.c:47-116)
Broadcast TGD never touches the satellite clock; instead pseudoranges are corrected:
- `gettgd`: returns `CLIGHT·eph.tgd[0]` (first TGD slot: GPS/QZS TGD, GAL BGD(E5a/E1), BDS TGD(B1/B3)).
- γ = (f1/f2)² (f2 = L5 for GAL/SBS).
- If no P1-P2 DCB product: `P1_P2 = (1−γ)·TGD` (line 93-95) — converts TGD into the P1-P2 bias domain.
- Single-frequency: `PC = P1 − P1_P2/(1−γ)` (109) — removes P1-P2 DCB mapped onto P1.
- IFLC option: `PC = (γP1 − P2)/(γ−1)` (103).
- SBAS mode: `PC -= P1_C1` (111) because SBAS clocks are referenced to C1.
- C/A→P code conversions via P1_C1/P2_C2 DCBs when obs is C1/C2.
- Bias variance contribution `SQR(ERR_CBIAS=0.3 m)`.

### 4.4 `ionmodel()` — Klobuchar (rtkcmn.c:3275-3315)
All angles in semi-circles (x/π), exactly per IS-GPS-200:
- Default coefficient set substituted if input all-zero (2004/1/1 values listed inline, 3278-3281).
- Earth-center angle `ψ = 0.0137/(el/π + 0.11) − 0.022`.
- IPP latitude clamped to ±0.416 semicircles (=±74.88°); longitude `λ = λ_u/π + ψ·sin(az)/cos(φ·π)`.
- Geomagnetic tweak: `φ += 0.064·cos((λ−1.617)·π)` (meridian ≈ −11.6° E).
- Local time `tt = 43200·λ + seconds_of_week mod 86400`.
- Obliquity factor `f = 1 + 16·(0.53 − el/π)³`.
- Amplitude floored at 0; period floored at **72000 s (20 h)**; phase centered at **50400 s (14:00 local)**; cosine argument `x = 2π(tt−50400)/per`.
- `delay = c·f·(|x|<1.57 ? 5E-9 + amp·(1 − x²/2 + x⁴/24) : 5E-9)` — night-time floor 5 ns.
- Guards: height < −1000 m or el ≤ 0 → 0 delay.

Related: `ionmapf` (3322-3326) single-layer mapping at HION = 350000 m (rtklib.h:69): `1/cos(asin((Re+h)/(Re+HION)·cos el))`. `ionppp` (3338-3358) rigorous spherical-trig pierce point incl. polar longitude-branch fix, slant factor `1/√(1−rp²)`.

### 4.5 `ionocorr()` — actually lives in pntpos.c:129-163
Option dispatch (IONOOPT_* enum, rtklib.h:331-339):

| opt | Source | Function / notes |
|---|---|---|
| 0 OFF | none | ion=0, var=SQR(ERR_ION=5.0 m) |
| 1 BRDC | nav→ion_gps Klobuchar | ionmodel; var=SQR(ion·ERR_BRDCI=0.5) |
| 2 SBAS | SBAS Message 18/26 grid | sbsioncorr (sbas.c:667) |
| 3 IFLC | dual-freq IF combo | handled in prange, not here |
| 4 EST | estimated per-sat | (RTK Kalman states) |
| 5 TEC | IONEX files | iontec(opt=1) (ionex.c:429) |
| 6 QZS | nav→ion_qzs Klobuchar | ionmodel with QZSS coefficients |
| 7 LEX | QZSS LEX MADOCA-style surface | lexioncorr (qzslex.c:600) |
| 8 STEC | slant TEC | (external) |

Two details worth copying:
- **First LSQ iteration forces BRDC iono + SAAS tropo** (`iter>0 ? opt->ionoopt : IONOOPT_BRDC`, rescode:237-246) so iteration 0 doesn't depend on the (bad) a-priori position.
- **Non-L1 frequency scaling** (rescode:240-243): `dion *= (λ_satL1/λ_carr[0])²` where λ_carr[0]=c/FREQ1 — correctly scales Klobuchar L1 delays for BeiDou B1 (1561.098 MHz) and GLONASS FDMA channels.

### 4.6 SBAS iono: `sbsioncorr()` (sbas.c:667-724)
- Shell re = 6378.1363 km, hion = 350 km; IPP via ionppp; 4 surrounding IGPs searched (`searchigp`).
- Bilinear weights w={ws,wn,es,en}; three-corner fallback configurations with negative-weight guard → fail.
- Delay scaled by obliquity fp; variance `Σ wᵢ·varicorr(GIVE)·9E-8·|age|` — **GIVE variance table** (sbas.c:113-116) `{0.0084,0.0333,...,187.0826}` m² with time-degradation coefficient 9e-8 s⁻¹.
- UDRE variance table at sbas.c:105-107; fast-correction degradation coefficients `degfcorr` (120-127).

### 4.7 IONEX TEC: `iontec()`/`iondelay()`/`interptec()` (ionex.c)
- Needs two maps bracketing the epoch (typically 2 h apart); delay computed on both maps then **linearly blended in time**; nearest map used if only one side valid; MIN_EL=0 rad, MIN_HGT=−1000 m; no-data variance VAR_NOTEC=(30 m)² (24-26).
- Per layer: pierce point, optional modified-SLM (`rp = rb/(rb+h)·sin(0.9782·zd)`, opt bit 2), optional **sun-fixed frame earth-rotation compensation** `lon += 2π·Δt/86400` (opt bit 1 — ionocorr passes opt=1, so ON).
- TECU→L1 metres: `40.30E16/FREQ1²`.
- Spatial: **bilinear in lat/lon** requiring all four corners >0; else nearest-quadrant corner; else mean of available corners (356-372).

### 4.8 LEX iono: `lexioncorr()` (qzslex.c:600-669)
Regional (Japan) quadratic/cubic surface about anchor `lexion.pos0`: slant delay = `F·Σ_{n≤2,m≤1} Enm·Δlat^n·Δlon^m`, F = 1/√(1−rp²) with hion=350 km; validity `|Δt| ≤ lexion.tspan`; geographic coverage check compiled out (#if 0).

### 4.9 `tropmodel()` — Saastamoinen (rtkcmn.c:3367-3387)
Standard atmosphere: `p = 1013.25·(1−2.2557E-5·h)^5.2568 hPa`; `T = 15 − 6.5E-3·h + 273.16 K`; partial water pressure `e = 6.108·H·exp((17.15T−4684)/(T−38.45))` with **relative humidity hard-default 0.7** (REL_HUMI, pntpos.c:36).
- Zenith hydrostatic: `0.0022768·p / (1 − 0.00266·cos2φ − 0.00028·h_km)`
- Zenith wet: `0.002277·(1255/T + 0.05)·e`
- **Mapping for both: plain 1/cos z** (no Niell split at this level!). Guards: h∈[−100, 10000] m, el>0.

### 4.10 `tropmapf()` — Niell Mapping Function default, GMF optional (3456-3485)
- Default build: **NMF** (`nmf()`, 3399-3442): Niell 1996 coefficient table at latitudes 15/30/45/60/75° (hydrostatic avg+amplitude, wet), linear interp between lat rows (`interpc`), annual variation `cos(2π(doy−28)/365.25)` with half-year phase flip for southern hemisphere; ellipsoidal height correction `aht={2.53E-5, 5.49E-3, 1.14E-3}` applied as `dm=(1/sinel − mf_aht(el))·h_km`; continued-fraction map `mapf(el,a,b,c)` (3394-3398). Returns dry MF, optionally wet MF.
- Compiled with `IERS_MODEL`: calls Fortran **GMF** (`gmf_`) with geoid-corrected height.
- Valid range: height ∈ [−1000, 20000] m, el > 0.

### 4.11 SBAS/MOPS troposphere: `sbstropcorr()` (sbas.c:753-782)
Latitude-dependent met tables (`getmet`, 15-75° rows, 10 params), seasonal modulation `cos(2π(doy−28 N-hemi /211 S-hemi)/365.25)`, zenith hydrostatic/wet via k1=77.604, k2=382000, Rd=287.054, gm=9.784; exponential height scaling; **mapping `m = 1.001/√(0.002001+sin²el)`**; σ = 0.12 m·m. Zenith values cached per receiver position (>1 m move triggers refresh).

### 4.12 Measurement error models — `varerr()`
**SPP code (pntpos.c:39-46):**
```c
fact = SYS_GLO?1.5 : SYS_SBS?3.0 : 1.0        // EFACT_GLO/SBS, others 1.0
var  = err[0]²·(err[1]² + err[2]²/sin(el))    // defaults err={100, 0.003, 0.003, 0, 1}
if IFLC: var ×= 9                              // ×3 sigma for IF combination
var  ×= fact²
```
System UERE factors (rtklib.h:88-93): **GPS/GAL/QZS/BDS = 1.0, GLONASS = 1.5, SBAS = 3.0**. With defaults, σ_code ≈ 100·0.003·√(1+1/sin el) ≈ 0.42 m @zenith … 1.2 m @5°.

**RTK double-difference (rtkpos.c:367-393):**
```c
c = err[3]·baseline/10km            // baseline-dependent (default err[3]=0)
d = (c·sclkstab·dt)²                // reference-rover clock stability term, sclkstab=5E-12
a,b = exterr tables if enabled, else fact·err[1], fact·err[2]
var = 2·(a² + b²/sin²el + c²) + d²  // ×2 for single-diff→double-diff
IFLC: ×3; fact includes eratio[code]≈100 and EFACT_GLO/SBS
extended per-signal error tables indexed sys∈{GPS→0,GLO→1,GAL→2}, per freq/code
```

### 4.13 Small utilities worth noting
- `satwavelen` (3162-3190): GLONASS FDMA `λ=c/(FREQn_GLO + DFRQn_GLO·k)`, DFRQ1=562.5 kHz, DFRQ2=437.5 kHz (rtklib.h:80,82); **BeiDou quirk: freq index order is B1(1561.098), B3(1268.52), B2(1207.14)** — index 1 is B3, not B2.
- `testsnr` (501-516): SNR mask, 9 knots every 10° from 5°, linear interpolation.
- `satexclude` (474-491): svh<0 → excluded; user exclude/include lists; navsys filter; **QZS health bit 0 masked (`svh&=0xFE`) so LEX health doesn't kill positioning**.
- `dops` (3242-3266), chi-square + max-GDOP solution validation (pntpos.c valsol:279-307, default maxgdop=30).

---

## 5. Unusual / clever items a competitor might miss

1. **Earth-rotation pre-rotation of SP3 nodes before polynomial interpolation** (preceph.c:435-441). Interpolating raw ECEF positions across 15-min epochs is wrong by tens of meters; rotating nodes by ωe·Δt into the target epoch's frame first fixes it. Single most impactful trick in the file.
2. **Post-interpolation relativistic clock repair** (preceph.c:624-626): linearly interpolated SP3 clocks lack the ±23 ns periodic term; recomputing `−2r·v/c²` from interpolated position/finite-differenced velocity restores it. Combined effect: cheap interpolation + correct dynamics.
3. **Graceful clock fallback ladder** (ephemeris.c:755-760 + preceph.c:491-496): separate CLK file → SP3-embedded clock → broadcast clock with σ=30 m, never hard-failing for lack of precise clock while keeping orbit quality.
4. **Extrapolation inflation instead of rejection**: EXTERR_CLK (1 mm/s per second) and EXTERR_EPH (5e-7 m/s², quadratic) let the engine keep working through short data gaps with honestly degraded sigmas feeding the estimator.
5. **Kepler solver discipline**: Newton with 1e-14 tolerance + 30-iter cap and explicit failure path (zeroed output) rather than silently returning garbage (ephemeris.c:203-209).
6. **BDS GEO −5° tilted-orbit rotation** implemented inline with precomputed SIN_5/COS_5 (ephemeris.c:223-234) — frequently forgotten in homebrew engines.
7. **GLONASS ICD equation bugfix** comment in deq() (ephemeris.c:260) — the ICD's published xdot[4]/xdot[5] equations contain errata; copy RTKLIB's corrected form (Coriolis signs).
8. **Satellite PCO combined into iono-free LC** (preceph.c:568-577), including the L1/L5 pairing switch for Galileo/SBAS — needed for consistency when processing IF combinations.
9. **First-iteration forcing of broadcast Klobuchar + Saastamoinen** in SPP (pntpos.c:237-246) decouples convergence from the initial guess.
10. **Rank-deficiency guards**: zero-variance 0.01 pseudo-rows for absent system clocks (pntpos.c:270-275); QZS LEX health bit masking (rtkcmn.c:485).
11. **Known gaps** (opportunities for a competitor): no yaw-attitude model for eclipse/noon-turns in satantoff; no stochastic clock jitter model for SP3; no ISC/or modern multi-frequency TGD usage beyond tgd[0..1]; BeiDou freq-index ordering (B3 as index 1) is a trap; tropospheric mapping in tropmodel is plain 1/cos(z) (Saastamoinen zenith × cosecant), with NMF only available via tropmapf for estimation modes.

---

## 6. Five most adoptable techniques

1. **Frame-consistent polynomial orbit interpolation** — rotate SP3 nodes by ωe·Δt before an order-10 Neville/Lagrange interpolation, clamp the 11-node window at array edges, require all nodes valid, and inflate σ quadratically outside the window (pephpos). Directly transplantable and testable.
2. **Velocity-based relativistic clock reconstruction** — after any clock interpolation, add `−2r·v/c²` using finite-difference velocities (1 ms step everywhere, both broadcast and precise paths). One formula fixes two interpolation sins at once.
3. **Layered clock fallback with honest variances** — precise clock → SP3 clock → broadcast clock, tagging fallbacks with σ=30 m so the estimator down-weights automatically; never reject for missing precise clock alone.
4. **Measurement-level TGD/DCB architecture** — keep broadcast clocks free of group delay; apply TGD-derived biases where observations enter (prange), with the `(1−γ)·TGD` P1-P2 synthesis fallback when no DCB product exists.
5. **Uniform variance plumbing end-to-end** — every correction source (URA/SSR URA/GIVE/UDRE/Klobuchar factor 0.5/Saastamoinen 0.3 m/IONEX RMS/no-TEC 30 m) returns a variance that flows into one weighted least-squares stack, with per-system UERE factors (GLO ×1.5, SBAS ×3.0) and elevation terms `σa²+σb²/sin²el`. Cheap, coherent, and immediately reusable in an EKF.
