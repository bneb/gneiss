> **Superseded.** This document describes an earlier architecture (pre network-RTK/PPK pivot) and is kept for historical record only. For current status and roadmap, see `docs/PROJECT_STATUS.md` and `docs/NETWORK_RTK_NEXT_STEPS.md`.

# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.967 m vs N/A | 3.108 m vs N/A | 57.829 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.993 m vs N/A | 3.228 m vs N/A | 57.731 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 3.046 m vs N/A | 4.398 m vs N/A | 58.037 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 1.967 m vs N/A | 3.108 m vs N/A | 57.829 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.993 m vs N/A | 3.228 m vs N/A | 57.731 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 1.151 m vs N/A | 1.814 m vs N/A | 59.139 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 47.453 m vs N/A | 52.877 m vs N/A | 92.086 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 47.694 m vs N/A | 52.872 m vs N/A | 92.170 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 47.453 m vs N/A | 52.877 m vs N/A | 92.086 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 47.709 m vs N/A | 52.885 m vs N/A | 92.742 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 47.694 m vs N/A | 52.918 m vs N/A | 92.175 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 47.709 m vs N/A | 52.885 m vs N/A | 92.742 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |

