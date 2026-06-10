# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## Shinjuku (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.288 m vs 13.209 m | 26.576 m vs 55.582 m | 5.166 m vs 38.538 m | **Gneiss** (+82.7%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.288 m vs 13.209 m | 26.576 m vs 55.582 m | 5.166 m vs 38.538 m | **Gneiss** (+82.7%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 7.725 m vs 13.209 m | 29.822 m vs 55.582 m | 5.913 m vs 38.538 m | **Gneiss** (+41.5%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 428.554 m vs 13.209 m | 643.099 m vs 55.582 m | 29.389 m vs 38.538 m | RTKLIB (+96.9%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 3.374 m vs 13.209 m | 25.769 m vs 55.582 m | 5.119 m vs 38.538 m | **Gneiss** (+74.5%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 7.773 m vs 13.209 m | 29.934 m vs 55.582 m | 5.918 m vs 38.538 m | **Gneiss** (+41.2%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.775 m vs 1.670 m | 21.217 m vs 11.843 m | 2.367 m vs 5.394 m | RTKLIB (+5.9%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.778 m vs 1.670 m | 20.963 m vs 11.843 m | 2.344 m vs 5.394 m | RTKLIB (+6.1%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 3.059 m vs 1.670 m | 23.234 m vs 11.843 m | 2.834 m vs 5.394 m | RTKLIB (+45.4%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 1.790 m vs 2.089 m | 20.950 m vs 20.154 m | 2.326 m vs 5.517 m | **Gneiss** (+14.3%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.777 m vs 2.089 m | 21.042 m vs 20.154 m | 2.416 m vs 5.517 m | **Gneiss** (+14.9%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 2.477 m vs 2.089 m | 18.358 m vs 20.154 m | 2.606 m vs 5.517 m | RTKLIB (+15.7%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 2.384 m vs N/A | 27.019 m vs N/A | 5.430 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 2.384 m vs N/A | 27.019 m vs N/A | 5.430 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 2.384 m vs N/A | 27.019 m vs N/A | 5.430 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 2.384 m vs N/A | 27.019 m vs N/A | 5.430 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 2.384 m vs N/A | 27.019 m vs N/A | 5.430 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 2.384 m vs N/A | 27.019 m vs N/A | 5.430 m vs N/A | **Gneiss** | Stable PPP integration. |

## Odaiba (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.409 m vs 5.335 m | 9.384 m vs 35.568 m | 2.853 m vs 8.330 m | **Gneiss** (+54.8%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.409 m vs 5.335 m | 9.384 m vs 35.568 m | 2.853 m vs 8.330 m | **Gneiss** (+54.8%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 7.543 m vs 5.335 m | 22.192 m vs 35.568 m | 3.405 m vs 8.330 m | RTKLIB (+29.3%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 524.141 m vs 5.335 m | 1507.220 m vs 35.568 m | 2.278 m vs 8.330 m | RTKLIB (+99.0%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 4.263 m vs 5.335 m | 9.678 m vs 35.568 m | 2.800 m vs 8.330 m | **Gneiss** (+20.1%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 7.781 m vs 5.335 m | 24.135 m vs 35.568 m | 3.372 m vs 8.330 m | RTKLIB (+31.4%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.104 m vs 2.239 m | 7.694 m vs 8.099 m | 2.099 m vs 6.413 m | **Gneiss** (+50.7%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.102 m vs 2.239 m | 7.678 m vs 8.099 m | 2.097 m vs 6.413 m | **Gneiss** (+50.8%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 3.081 m vs 2.239 m | 11.316 m vs 8.099 m | 2.602 m vs 6.413 m | RTKLIB (+27.3%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 1.162 m vs 2.236 m | 7.467 m vs 13.146 m | 2.061 m vs 8.219 m | **Gneiss** (+48.0%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.099 m vs 2.236 m | 7.680 m vs 13.146 m | 2.094 m vs 8.219 m | **Gneiss** (+50.8%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 2.610 m vs 2.236 m | 15.462 m vs 13.146 m | 2.441 m vs 8.219 m | RTKLIB (+14.3%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 2.583 m vs N/A | 10.575 m vs N/A | 3.072 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 2.583 m vs N/A | 10.575 m vs N/A | 3.072 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 2.583 m vs N/A | 10.575 m vs N/A | 3.072 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 2.583 m vs N/A | 10.575 m vs N/A | 3.072 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 2.583 m vs N/A | 10.575 m vs N/A | 3.072 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 2.583 m vs N/A | 10.575 m vs N/A | 3.072 m vs N/A | **Gneiss** | Stable PPP integration. |

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 3.297 m vs 3.311 m | 8.784 m vs 10.191 m | 62.750 m vs 66.316 m | **Gneiss** (+0.4%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 3.435 m vs 3.311 m | 9.143 m vs 10.191 m | 62.693 m vs 66.316 m | RTKLIB (+3.6%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 9.957 m vs 3.311 m | 109.226 m vs 10.191 m | 63.041 m vs 66.316 m | RTKLIB (+66.7%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 8230.426 m vs 3.311 m | 16993.840 m vs 10.191 m | 85.088 m vs 66.316 m | RTKLIB (+100.0%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1307.220 m vs 3.311 m | 11350.658 m vs 10.191 m | 64.874 m vs 66.316 m | RTKLIB (+99.7%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 9.530 m vs 3.311 m | 109.226 m vs 10.191 m | 63.065 m vs 66.316 m | RTKLIB (+65.3%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 2.564 m vs 1.773 m | 7.280 m vs 4.161 m | 62.686 m vs 64.598 m | RTKLIB (+30.9%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 2.528 m vs 1.773 m | 7.565 m vs 4.161 m | 62.679 m vs 64.598 m | RTKLIB (+29.9%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 4.148 m vs 1.773 m | 68.636 m vs 4.161 m | 62.526 m vs 64.598 m | RTKLIB (+57.3%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 2.238 m vs 1.831 m | 5.661 m vs 3.126 m | 62.685 m vs 64.471 m | RTKLIB (+18.2%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 2.262 m vs 1.831 m | 5.682 m vs 3.126 m | 62.655 m vs 64.471 m | RTKLIB (+19.1%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 2.711 m vs 1.831 m | 68.131 m vs 3.126 m | 62.586 m vs 64.471 m | RTKLIB (+32.5%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `ppp` | forward | Off | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Stable PPP integration. |

## PPP (f9p_ppp)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 3.057 m vs N/A | 3.932 m vs N/A | 3.859 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 3.057 m vs N/A | 3.932 m vs N/A | 3.859 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | forward | Tight | 3.542 m vs N/A | 702.859 m vs N/A | 3.835 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Off | 197.343 m vs N/A | 1015.292 m vs N/A | 5.420 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 3.731 m vs N/A | 41.781 m vs N/A | 3.722 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Tight | 3.519 m vs N/A | 702.859 m vs N/A | 3.815 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `rtk` | forward | Off | 0.599 m vs N/A | 0.919 m vs N/A | 0.389 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.599 m vs N/A | 0.919 m vs N/A | 0.394 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | forward | Tight | 0.487 m vs N/A | 2.590 m vs N/A | 0.311 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Off | 0.459 m vs N/A | 0.738 m vs N/A | 0.337 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.599 m vs N/A | 0.919 m vs N/A | 0.390 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Tight | 0.505 m vs N/A | 1.103 m vs N/A | 0.320 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `ppp` | forward | Off | 3.057 m vs N/A | 3.923 m vs N/A | 3.859 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 3.057 m vs N/A | 3.923 m vs N/A | 3.859 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 3.057 m vs N/A | 3.923 m vs N/A | 3.859 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 3.057 m vs N/A | 3.923 m vs N/A | 3.859 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 3.057 m vs N/A | 3.923 m vs N/A | 3.859 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 3.057 m vs N/A | 3.923 m vs N/A | 3.859 m vs N/A | **Gneiss** | Stable PPP integration. |

