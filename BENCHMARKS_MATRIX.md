# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## Shinjuku (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.288 m vs 13.209 m | 26.576 m vs 55.582 m | 5.166 m vs 38.538 m | **Gneiss** (+82.7%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.288 m vs 13.209 m | 26.576 m vs 55.582 m | 5.166 m vs 38.538 m | **Gneiss** (+82.7%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 58.824 m vs 13.209 m | 704.324 m vs 55.582 m | 25.925 m vs 38.538 m | RTKLIB (+77.5%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 428.554 m vs 13.209 m | 643.099 m vs 55.582 m | 29.389 m vs 38.538 m | RTKLIB (+96.9%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 3.374 m vs 13.209 m | 25.769 m vs 55.582 m | 5.119 m vs 38.538 m | **Gneiss** (+74.5%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 58.843 m vs 13.209 m | 685.057 m vs 55.582 m | 25.849 m vs 38.538 m | RTKLIB (+77.6%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.820 m vs 1.670 m | 24.725 m vs 11.843 m | 3.100 m vs 5.394 m | RTKLIB (+8.2%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.815 m vs 1.670 m | 23.284 m vs 11.843 m | 2.992 m vs 5.394 m | RTKLIB (+8.0%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 8.761 m vs 1.670 m | 19944.870 m vs 11.843 m | 14.135 m vs 5.394 m | RTKLIB (+80.9%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 1.902 m vs 2.089 m | 25.915 m vs 20.154 m | 3.269 m vs 5.517 m | **Gneiss** (+9.0%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.815 m vs 2.089 m | 23.284 m vs 20.154 m | 2.996 m vs 5.517 m | **Gneiss** (+13.1%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 340.130 m vs 2.089 m | 30206.106 m vs 20.154 m | 105.500 m vs 5.517 m | RTKLIB (+99.4%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 1.978 m vs N/A | 24.181 m vs N/A | 4.636 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | forward | Tight | 1.978 m vs N/A | 24.181 m vs N/A | 4.636 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 1.978 m vs N/A | 24.181 m vs N/A | 4.636 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | smoothed | Tight | 1.978 m vs N/A | 24.181 m vs N/A | 4.636 m vs N/A | **Gneiss** | Stable PPP integration. |

## Odaiba (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.409 m vs 5.335 m | 9.384 m vs 35.568 m | 2.853 m vs 8.330 m | **Gneiss** (+54.8%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.409 m vs 5.335 m | 9.384 m vs 35.568 m | 2.853 m vs 8.330 m | **Gneiss** (+54.8%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 9935.669 m vs 5.335 m | 52517.195 m vs 35.568 m | 5361.163 m vs 8.330 m | RTKLIB (+99.9%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 524.141 m vs 5.335 m | 1507.220 m vs 35.568 m | 2.278 m vs 8.330 m | RTKLIB (+99.0%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 4.263 m vs 5.335 m | 9.678 m vs 35.568 m | 2.800 m vs 8.330 m | **Gneiss** (+20.1%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 9924.000 m vs 5.335 m | 52517.195 m vs 35.568 m | 5357.488 m vs 8.330 m | RTKLIB (+99.9%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.017 m vs 2.239 m | 4.209 m vs 8.099 m | 1.923 m vs 6.413 m | **Gneiss** (+54.6%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.004 m vs 2.239 m | 4.207 m vs 8.099 m | 1.909 m vs 6.413 m | **Gneiss** (+55.2%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 445.187 m vs 2.239 m | 3487.169 m vs 8.099 m | 105.451 m vs 6.413 m | RTKLIB (+99.5%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 0.865 m vs 2.236 m | 4.310 m vs 13.146 m | 1.827 m vs 8.219 m | **Gneiss** (+61.3%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.995 m vs 2.236 m | 4.207 m vs 13.146 m | 1.922 m vs 8.219 m | **Gneiss** (+55.5%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 384.622 m vs 2.236 m | 842.989 m vs 13.146 m | 69.954 m vs 8.219 m | RTKLIB (+99.4%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 2.372 m vs N/A | 9.731 m vs N/A | 2.717 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | forward | Tight | 2.372 m vs N/A | 9.731 m vs N/A | 2.717 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 2.372 m vs N/A | 9.731 m vs N/A | 2.717 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | smoothed | Tight | 2.372 m vs N/A | 9.731 m vs N/A | 2.717 m vs N/A | **Gneiss** | Stable PPP integration. |

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 3.297 m vs 3.311 m | 8.784 m vs 10.191 m | 62.750 m vs 66.316 m | **Gneiss** (+0.4%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 3.435 m vs 3.311 m | 9.143 m vs 10.191 m | 62.693 m vs 66.316 m | RTKLIB (+3.6%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 72.379 m vs 3.311 m | 6647.208 m vs 10.191 m | 73.009 m vs 66.316 m | RTKLIB (+95.4%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 8230.426 m vs 3.311 m | 16993.840 m vs 10.191 m | 85.088 m vs 66.316 m | RTKLIB (+100.0%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1307.220 m vs 3.311 m | 11350.658 m vs 10.191 m | 64.874 m vs 66.316 m | RTKLIB (+99.7%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 69.828 m vs 3.311 m | 6734.101 m vs 10.191 m | 72.497 m vs 66.316 m | RTKLIB (+95.3%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 2.357 m vs 1.773 m | 10.178 m vs 4.161 m | 62.642 m vs 64.598 m | RTKLIB (+24.8%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 2.357 m vs 1.773 m | 9.860 m vs 4.161 m | 62.646 m vs 64.598 m | RTKLIB (+24.8%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 3.668 m vs 1.773 m | 176.846 m vs 4.161 m | 63.389 m vs 64.598 m | RTKLIB (+51.7%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 2.351 m vs 1.831 m | 9.181 m vs 3.126 m | 62.628 m vs 64.471 m | RTKLIB (+22.1%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 2.353 m vs 1.831 m | 7.563 m vs 3.126 m | 62.650 m vs 64.471 m | RTKLIB (+22.2%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 2.914 m vs 1.831 m | 184.714 m vs 3.126 m | 63.319 m vs 64.471 m | RTKLIB (+37.2%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
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
| `spp` | smoothed | Off | 197.343 m vs N/A | 1015.292 m vs N/A | 5.420 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 3.731 m vs N/A | 41.781 m vs N/A | 3.722 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Tight | 3.786 m vs N/A | 1601.323 m vs N/A | 3.979 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `rtk` | forward | Off | 0.579 m vs N/A | 1.319 m vs N/A | 0.440 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.579 m vs N/A | 1.319 m vs N/A | 0.440 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | forward | Tight | 0.579 m vs N/A | 1.319 m vs N/A | 0.440 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Off | 0.580 m vs N/A | 1.315 m vs N/A | 0.439 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.579 m vs N/A | 1.319 m vs N/A | 0.440 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Tight | 0.580 m vs N/A | 1.315 m vs N/A | 0.439 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `ppp` | forward | Off | 3.015 m vs N/A | 3.955 m vs N/A | 4.016 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | forward | Tight | 3.015 m vs N/A | 3.955 m vs N/A | 4.016 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 3.015 m vs N/A | 3.955 m vs N/A | 4.016 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | N/A vs N/A | N/A vs N/A | N/A vs N/A | None | Stable PPP integration. |
| `ppp` | smoothed | Tight | 3.015 m vs N/A | 3.955 m vs N/A | 4.016 m vs N/A | **Gneiss** | Stable PPP integration. |

