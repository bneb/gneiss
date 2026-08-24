# OHLN End-of-Day Anomaly — Investigation Log

## Symptom
OHLN (16.5 km baseline, sea-level shore station) shows a single sustained
vertical excursion in the last ~37 minutes of the GPS day (TOW 429780–
431970, i.e. 23:23–24:00 GPS = 15:23–16:00 PDT): smoothed vertical error
ramps smoothly from +0.1 m to +3.0 m and back partway. This one window
accounts for essentially ALL of OHLN's vRMS (~433 mm; 74 epochs at ~2.5 m
over 2880 → √(74·6.25/2880) ≈ 0.40 m).

## Established facts (measured, not assumed)
1. **Forward pass diverges, not backward.** Forward-dump shows the ramp;
   combiner follows it. Both passes claim q=1 (fixed) throughout.
2. **Displacement is continuous** (+14→+90→−57 mm/epoch), NOT staircase
   AR fix-jumps, despite static-Q lock (1e-9) active. One discrete step
   of +287 mm at TOW 429930 embedded mid-ramp.
3. **Receiver clocks ruled out.** Base-minus-rover pseudorange single-
   differences (scripts/receiver_clock_probe.py, byte-verified parser):
   OHLN relative clock mean +1.6 µs, drift −0.0075 ns/s, residual RMS
   1.49 µs, max epoch jump 5.9 µs — statistically identical to clean-
   control P181 (RMS 1.08 µs). No steps, no EOD signature.
4. **No slip storm visible.** Geometry-free jumps >0.5 cyc: OHLN has 11
   all day (cleanest file; CAPO 242, P181 140). LLI flags: zero at OHLN
   vs hundreds elsewhere — suspected TPS NET-G3A converter quirk, so
   LLI-based screening is blind on this stream, but independent GF and
   SNR checks show nothing anomalous. SNR at EOD 44.6 dBHz (slightly
   better than midday).
5. **No reference-satellite switch** during the window (G25 holds);
   sat count dips only 8↔7.
6. **Static water-multipath geometry insufficient.** Midday had MORE
   low-elevation over-bay satellites (G6@5°/az268, G22@10°/az164,
   G30@20°/az259) than the blowup window (only G21 marginal) — yet no
   midday blowup.
7. **Cascade AR never fixes on OHLN**: 1,232 veto→cascade-fail plus 681
   float→cascade-fail across the day (0 cascade fixes). Veto storms
   cluster at hours 07-09, 13, 18, 20-21 — hours BEFORE the blowup.

## Site context
- OHLN: "Ohlone Park" BARD station (UC Berkeley Seismo Lab), San Pablo
  Bay shoreline, ellipsoid height −0.5 m. Antenna in 2020: ASH701945B_M
  SCIT (Ashtech choke ring) — the ONLY non-Trimble antenna in the
  benchmark set. NGS coord file describes a Mar-2024 SEPCHOKE install,
  NOT the 2020-era hardware.
- P224 rover: Sibley Volcanic Preserve, Berkeley Hills, height 407 m —
  dead-center in the typical SF-Bay marine-layer inversion band
  (300–500 m). Baseline crosses the bay northward.
- Blowup clock: GPS 23:23 ≈ 16:23 PDT = late-afternoon marine surge.

## Working hypothesis (unproven)
Marine-layer wet-delay decorrelation: base column inside moist layer,
rover column near/above inversion top. Our filter estimates ONE rover-
side ZWD scalar mapped by ROVER mapping functions — a BASE-side residual
with DIFFERENT mapping functions cannot be represented and injects
elevation-dependent bias into every DD pair. During afternoon surge this
grows coherently; with ambiguities confidently held (q=1) and innovation
gate at 500 cycles far too loose to trip on dm-scale biases, the locked-
Q solution is dragged vertically. Uniquely severe at OHLN because only
it straddles the inversion at sea level on a bay-crossing ray.

Alternative candidates still open: same-cycle dual-frequency slips
(invisible to both LLI and GF checks) — needs sat-pair DD slip detection
to rule out; tidal/wind-driven water-surface multipath variation.

## Ruled out (this session)
- Receiver clock steps/drift (probe above)
- Ephemeris staleness (nav ToC covers full day, uploads every 2 h)
- Unhealthy SVs (all 474 records health=0)
- Seismic activity (max M2.4 ~100 km away, USGS catalog 2020-05-14)
- Reference-satellite switching
- Static (time-invariant) water multipath geometry

## Tooling built (committed)
- scripts/receiver_clock_probe.py — RINEX2 nav+obs parser with
  fixed-width cell extraction (LLI/SSI glued-digit safe), week-aligned
  ephemeris selection, Kepler propagation, per-epoch inter-receiver
  clock via median pseudorange single-difference.
- scripts/test_receiver_clock_probe.py — red/green suite asserting
  hand-transcribed byte-level ground truth (types lists incl. 20-type
  wrapped header, exact first-epoch observables, nav af0 to 1e-12,
  2880-epoch grid integrity, <500 µs physics gate on clock-term spread;
  achieved worst=10.3 µs).
- eval_network_ppk WL_DUMP env: per-epoch tow,h,v,q,sep,nsat CSVs per
  base/pass for cross-base correlation studies.
- ar-decision debug line (tow-tagged) in process_epoch veto/cascade path.

## Next actions
1. Sat-pair double-difference slip detector (same-cycle L1+L2 slips are
   invisible to GF; must difference phase against DD range model with
   reference-SV differencing to cancel base-coordinate/orbit error).
2. Quantify marine-layer hypothesis: fit per-base elevation-dependent
   residual bias vs PDT hour; predict OHLN-only signature.
3. Engineering fix candidates once driver confirmed:
   - Two-station tropo (base-side ZWD from network or model inversion)
   - Per-base cascade threshold tuning (cascade 0-for-1913 today)
   - Innovation-gate scaling by baseline length
