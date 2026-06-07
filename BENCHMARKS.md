# Gneiss Comprehensive Benchmarks

This document empirically maps the performance of Gneiss across varying datasets compared against RTKLIB (demo5).

## RTK Kinematic Accuracy Benchmarks

### Shinjuku (u-blox automotive)
- **Gneiss**: 1.454 m (Hz50), 1.314 m (Vt50)
- **RTKLIB**: 1.670 m (Hz50), 5.394 m (Vt50)
- 🏆 **Gneiss WINS by 12.9%** (Horizontal 50th percentile)

### Odaiba (u-blox automotive)
- **Gneiss**: 1.271 m (Hz50), 2.465 m (Vt50)
- **RTKLIB**: 2.239 m (Hz50), 6.413 m (Vt50)
- 🏆 **Gneiss WINS by 43.2%** (Horizontal 50th percentile)

### PPP (u-blox f9p)
- **Gneiss**: 0.424 m (Hz50)
- **RTKLIB**: Failed to produce solutions.
- 🏆 **Gneiss WINS**

*(Note: On the GSDC smartphone dataset, RTKLIB still currently has an edge due to different cycle slip handling for smartphone observables, but we successfully achieved our goal of beating RTKLIB on the core u-blox automotive datasets).*

## SPP Global Fallback
| Metric | Gneiss SPP | RTKLIB SPP |
|---|---|---|
| Horiz 50th % | **3.297 m** | 3.520 m |
| Horiz 95th % | 8.784 m | **7.490 m** |

*Benchmarks last updated: June 2026*
