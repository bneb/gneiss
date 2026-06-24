# RALPH Accuracy Drive — Final Report

## 2026-06-22 23:30 UTC — LOOP TERMINATED

### Reason: Cannot produce valid PPP benchmarks on available datasets

### Ground Truth Audit
| Dataset | "Truth" | Actual | Independent? | PPP Viable? |
|---------|---------|--------|-------------|-------------|
| Odaiba | SPAN-CPT ✅ | reference.csv | ✅ Yes | ❌ <4 sats 90%+ |
| Shinjuku | SPAN-CPT ✅ | reference.csv | ✅ Yes | ❌ <4 sats most epochs |
| f9p | RTKLIB RTK | rover_ppk.pos | ❌ Competing solver | ✅ Open sky |

### Best Validated Results
| Dataset | Gneiss PPP | Gneiss RTK | Notes |
|---------|-----------|-----------|-------|
| Odaiba | 330m (urban canyon) | — | SPAN-CPT truth, not enough sats for PPP |
| Shinjuku | 10km (urban canyon) | — | SPAN-CPT truth, not enough sats for PPP |
| f9p | 2.44m | TBD | No independent truth available |

### Production Fixes Shipped (873 tests pass)
| Fix | Impact |
|-----|--------|
| process_noise_cb: 1e6→1.0 | Single most impactful — prevents clock variance explosion |
| Clock model: φ=1.0 random walk | Preserves temporal correlation for PPP convergence |
| Elevation mask: 15°→5° | 3× more satellites in urban environments |
| InsufficientSatellites: coast+decouple | State preserved through gaps |
| State history: coasted epochs saved | Smoother can bridge gaps |
| UDUC AR auto-enable: removed | Eliminated 18km regression with SP3/CLK |
| Automatic clock calibration | Data-driven process noise from first 200 epochs |

### Blockers
1. **No PPP-viable dataset with independent ground truth.** Odaiba/Shinjuku have truth but insufficient satellites. f9p has satellites but no truth.
2. **PPP requires dual-frequency measurements from ≥4 satellites** — not available in urban canyons.
3. **f9p needs surveyed position or trusted RTK solution** as ground truth.

### Next Step
Deploy f9p receiver at a surveyed location, or acquire open-sky dataset with known ground truth. Then PPP benchmarks become meaningful.
