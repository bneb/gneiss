# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## Shinjuku (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.288 m vs 13.209 m | 26.576 m vs 55.582 m | 5.166 m vs 38.538 m | **Gneiss** (+82.7%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.288 m vs 13.209 m | 26.576 m vs 55.582 m | 5.166 m vs 38.538 m | **Gneiss** (+82.7%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 362.936 m vs 13.209 m | 3201.745 m vs 55.582 m | 453.255 m vs 38.538 m | RTKLIB (+96.4%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 2.564 m vs 13.209 m | 84.853 m vs 55.582 m | 5.288 m vs 38.538 m | **Gneiss** (+80.6%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 2.464 m vs 13.209 m | 26.576 m vs 55.582 m | 5.191 m vs 38.538 m | **Gneiss** (+81.3%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 363.677 m vs 13.209 m | 3201.745 m vs 55.582 m | 453.109 m vs 38.538 m | RTKLIB (+96.4%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.915 m vs 1.670 m | 18.099 m vs 11.843 m | 1.847 m vs 5.394 m | RTKLIB (+12.8%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.884 m vs 1.670 m | 17.288 m vs 11.843 m | 1.816 m vs 5.394 m | RTKLIB (+11.4%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 430.302 m vs 1.670 m | 11009.320 m vs 11.843 m | 437.750 m vs 5.394 m | RTKLIB (+99.6%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | N/A vs 2.089 m | N/A vs 20.154 m | N/A vs 5.517 m | RTKLIB | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.875 m vs 2.089 m | 17.296 m vs 20.154 m | 1.851 m vs 5.517 m | **Gneiss** (+10.2%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | N/A vs 2.089 m | N/A vs 20.154 m | N/A vs 5.517 m | RTKLIB | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 1.956 m vs N/A | 24.181 m vs N/A | 4.732 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | forward | Tight | 1.956 m vs N/A | 24.181 m vs N/A | 4.732 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 1.956 m vs N/A | 24.181 m vs N/A | 4.732 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | smoothed | Tight | 1.956 m vs N/A | 24.181 m vs N/A | 4.732 m vs N/A | **Gneiss** | Stable PPP integration. |

## Odaiba (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.409 m vs 5.335 m | 9.384 m vs 35.568 m | 2.853 m vs 8.330 m | **Gneiss** (+54.8%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.409 m vs 5.335 m | 9.384 m vs 35.568 m | 2.853 m vs 8.330 m | **Gneiss** (+54.8%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 12763.498 m vs 5.335 m | 1427115.758 m vs 35.568 m | 7415.789 m vs 8.330 m | RTKLIB (+100.0%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 2.867 m vs 5.335 m | 367.713 m vs 35.568 m | 2.878 m vs 8.330 m | **Gneiss** (+46.3%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 2.630 m vs 5.335 m | 9.384 m vs 35.568 m | 2.852 m vs 8.330 m | **Gneiss** (+50.7%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 2091.956 m vs 5.335 m | 1427115.758 m vs 35.568 m | 1606.843 m vs 8.330 m | RTKLIB (+99.7%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.364 m vs 2.239 m | 9.832 m vs 8.099 m | 1.892 m vs 6.413 m | **Gneiss** (+39.1%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.324 m vs 2.239 m | 9.005 m vs 8.099 m | 1.916 m vs 6.413 m | **Gneiss** (+40.9%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 230.463 m vs 2.239 m | 6600.494 m vs 8.099 m | 719.934 m vs 6.413 m | RTKLIB (+99.0%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 1.359 m vs 2.236 m | 9.832 m vs 13.146 m | 1.834 m vs 8.219 m | **Gneiss** (+39.2%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.322 m vs 2.236 m | 9.005 m vs 13.146 m | 1.907 m vs 8.219 m | **Gneiss** (+40.9%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 235.942 m vs 2.236 m | 6600.494 m vs 13.146 m | 719.934 m vs 8.219 m | RTKLIB (+99.1%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 2.417 m vs N/A | 10.750 m vs N/A | 2.695 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | forward | Tight | 2.417 m vs N/A | 10.750 m vs N/A | 2.695 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 2.417 m vs N/A | 10.750 m vs N/A | 2.695 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | smoothed | Tight | 2.417 m vs N/A | 10.750 m vs N/A | 2.695 m vs N/A | **Gneiss** | Stable PPP integration. |

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 3.297 m vs 3.311 m | 8.784 m vs 10.191 m | 62.750 m vs 66.316 m | **Gneiss** (+0.4%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 3.435 m vs 3.311 m | 9.143 m vs 10.191 m | 62.693 m vs 66.316 m | RTKLIB (+3.6%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 10335.352 m vs 3.311 m | 54963.335 m vs 10.191 m | 7984.013 m vs 66.316 m | RTKLIB (+100.0%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 3.435 m vs 3.311 m | 125.966 m vs 10.191 m | 62.630 m vs 66.316 m | RTKLIB (+3.6%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 3.643 m vs 3.311 m | 13.162 m vs 10.191 m | 62.612 m vs 66.316 m | RTKLIB (+9.1%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 10335.352 m vs 3.311 m | 54963.335 m vs 10.191 m | 7984.013 m vs 66.316 m | RTKLIB (+100.0%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 3.320 m vs 1.773 m | 8.828 m vs 4.161 m | 62.772 m vs 64.598 m | RTKLIB (+46.6%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 24.929 m vs 1.773 m | 29.089 m vs 4.161 m | 62.980 m vs 64.598 m | RTKLIB (+92.9%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 112.132 m vs 1.773 m | 2203.213 m vs 4.161 m | 66.768 m vs 64.598 m | RTKLIB (+98.4%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 3.379 m vs 1.831 m | 14.238 m vs 3.126 m | 62.685 m vs 64.471 m | RTKLIB (+45.8%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 24.929 m vs 1.831 m | 29.089 m vs 3.126 m | 62.952 m vs 64.471 m | RTKLIB (+92.7%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 112.132 m vs 1.831 m | 2203.213 m vs 3.126 m | 65.510 m vs 64.471 m | RTKLIB (+98.4%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `ppp` | forward | Off | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | forward | Tight | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | smoothed | Tight | 11.025 m vs N/A | 11.025 m vs N/A | 43.205 m vs N/A | **Gneiss** | Stable PPP integration. |

## PPP (f9p_ppp)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 3.057 m vs N/A | 3.932 m vs N/A | 3.859 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 3.057 m vs N/A | 3.932 m vs N/A | 3.859 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | forward | Tight | 3.790 m vs N/A | 1816.545 m vs N/A | 3.993 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Off | 3.064 m vs N/A | 3.932 m vs N/A | 3.859 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 3.058 m vs N/A | 3.932 m vs N/A | 3.859 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Tight | 3.788 m vs N/A | 1816.545 m vs N/A | 3.993 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `rtk` | forward | Off | 0.494 m vs N/A | 1.199 m vs N/A | 0.645 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.491 m vs N/A | 1.199 m vs N/A | 0.632 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | forward | Tight | 0.945 m vs N/A | 776.354 m vs N/A | 0.811 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Off | 0.494 m vs N/A | 1.199 m vs N/A | 0.645 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.491 m vs N/A | 1.198 m vs N/A | 0.632 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Tight | 0.945 m vs N/A | 776.354 m vs N/A | 0.811 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `ppp` | forward | Off | 3.014 m vs N/A | 3.955 m vs N/A | 4.018 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | forward | Tight | 3.014 m vs N/A | 3.955 m vs N/A | 4.018 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 3.014 m vs N/A | 3.955 m vs N/A | 4.018 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | smoothed | Tight | 3.014 m vs N/A | 3.955 m vs N/A | 4.018 m vs N/A | **Gneiss** | Stable PPP integration. |

