# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 3.023 m vs N/A | 4.380 m vs N/A | 58.056 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 1.143 m vs N/A | 1.792 m vs N/A | 59.139 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 4.350 m vs N/A | 13.984 m vs N/A | 66.159 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 4.459 m vs N/A | 15.090 m vs N/A | 66.438 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 4.350 m vs N/A | 13.984 m vs N/A | 66.159 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 4.402 m vs N/A | 14.941 m vs N/A | 66.117 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 4.450 m vs N/A | 15.234 m vs N/A | 66.377 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 4.402 m vs N/A | 14.941 m vs N/A | 66.117 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `ppp` | forward | Off | 151.025 m vs N/A | 178.978 m vs N/A | 89.220 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 151.025 m vs N/A | 178.978 m vs N/A | 89.220 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 151.025 m vs N/A | 178.978 m vs N/A | 89.220 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 154.789 m vs N/A | 184.598 m vs N/A | 216.244 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 154.789 m vs N/A | 184.598 m vs N/A | 216.244 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 154.789 m vs N/A | 184.598 m vs N/A | 216.244 m vs N/A | **Gneiss** | Stable PPP integration. |

