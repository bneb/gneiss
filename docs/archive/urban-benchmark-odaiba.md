> **Superseded.** Historical urban canyon benchmark on Odaiba dataset from early PPP evaluation. See `docs/PROJECT_STATUS.md`.

Urban PPP Benchmark — Odaiba, Tokyo
====================================
Date: 2018-12-20 (GPS Week 2032, Day 354)
Receiver: u-blox F9P (dual-frequency GPS L1/L2)
Environment: Urban canyon (Odaiba, Tokyo)
Precise products: CODE MGEX SP3/CLK/BIA
Mode: GPS-only PPP (float, no AR — insufficient sats)

Results (6206 epochs, 10 Hz reference):
  All epochs:
    50% (median): 26.0 m
    68%: 68.4 m
    95%: 318.3 m
    RMS: 20,064 m (extreme outliers dominate)

  Excluding >500m outliers (2.7% of epochs):
    50%: 25.5 m
    68%: 61.3 m
    95%: 244.6 m
    RMS: 106.7 m

  Excluding >100m outliers (24.3% of epochs):
    50%: 22.9 m
    95%: 85.4 m
    RMS: 38.0 m

  Convergence: DIVERGING (first half 20.6m → last half 97.4m, 4.74x)

Analysis:
  - No coasting events — adaptive PR threshold working
  - AR failed (insufficient dual-frequency sats in urban canyon)
  - Filter degrades over time due to satellite dropouts/reacquisitions
  - 24.3% of epochs >100m error (likely deep canyon, <4 sats visible)
  - 2.7% extreme outliers >500m (likely <3 sats visible)
  - Expected improvement: multi-constellation (GPS+QZSS for Tokyo), INS coupling

Comparison with IGS open-sky stations (GPS-only, with AR):
  ALIC: 0.84m median — 31x better than urban
  CEDU: 0.73m median — 35x better than urban
