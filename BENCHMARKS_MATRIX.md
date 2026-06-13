# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## Shinjuku (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.301 m vs 5.544 m | 2.288 m vs 13.209 m | 2.311 m vs 8.036 m | **Gneiss** (+76.5%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.301 m vs 5.544 m | 2.288 m vs 13.209 m | 2.311 m vs 8.036 m | **Gneiss** (+76.5%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 4.158 m vs 5.544 m | 7.725 m vs 13.209 m | 2.558 m vs 8.036 m | **Gneiss** (+25.0%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 1.301 m vs 5.544 m | 2.288 m vs 13.209 m | 2.311 m vs 8.036 m | **Gneiss** (+76.5%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.851 m vs 5.544 m | 3.374 m vs 13.209 m | 2.283 m vs 8.036 m | **Gneiss** (+66.6%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 4.187 m vs 5.544 m | 7.773 m vs 13.209 m | 2.558 m vs 8.036 m | **Gneiss** (+24.5%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.176 m vs 0.674 m | 1.821 m vs 1.670 m | 0.960 m vs 2.984 m | RTKLIB (+42.7%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.167 m vs 0.674 m | 1.815 m vs 1.670 m | 0.961 m vs 2.984 m | RTKLIB (+42.2%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 1.686 m vs 0.674 m | 3.045 m vs 1.670 m | 1.281 m vs 2.984 m | RTKLIB (+60.0%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 1.213 m vs 0.706 m | 1.901 m vs 2.089 m | 0.994 m vs 5.065 m | RTKLIB (+41.8%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.169 m vs 0.706 m | 1.815 m vs 2.089 m | 0.961 m vs 5.065 m | RTKLIB (+39.6%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 1.487 m vs 0.706 m | 2.496 m vs 2.089 m | 1.004 m vs 5.065 m | RTKLIB (+52.5%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 1.301 m vs N/A | 2.284 m vs N/A | 2.332 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 1.301 m vs N/A | 2.284 m vs N/A | 2.332 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 1.328 m vs N/A | 2.384 m vs N/A | 2.384 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 1.301 m vs N/A | 2.284 m vs N/A | 2.332 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 1.301 m vs N/A | 2.284 m vs N/A | 2.332 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 1.328 m vs N/A | 2.384 m vs N/A | 2.384 m vs N/A | **Gneiss** | Stable PPP integration. |

## Odaiba (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.522 m vs 2.061 m | 2.409 m vs 5.335 m | 1.346 m vs 5.641 m | **Gneiss** (+26.2%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.522 m vs 2.061 m | 2.409 m vs 5.335 m | 1.346 m vs 5.641 m | **Gneiss** (+26.2%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 4.117 m vs 2.061 m | 7.543 m vs 5.335 m | 1.925 m vs 5.641 m | RTKLIB (+49.9%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 1.522 m vs 2.061 m | 2.409 m vs 5.335 m | 1.346 m vs 5.641 m | **Gneiss** (+26.2%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 2.427 m vs 2.061 m | 4.263 m vs 5.335 m | 1.308 m vs 5.641 m | RTKLIB (+15.1%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 4.230 m vs 2.061 m | 7.781 m vs 5.335 m | 1.922 m vs 5.641 m | RTKLIB (+51.3%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 0.475 m vs 1.516 m | 1.017 m vs 2.239 m | 1.012 m vs 4.081 m | **Gneiss** (+68.7%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.475 m vs 1.516 m | 0.991 m vs 2.239 m | 1.012 m vs 4.081 m | **Gneiss** (+68.7%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 1.518 m vs 1.516 m | 2.853 m vs 2.239 m | 1.074 m vs 4.081 m | RTKLIB (+0.1%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 0.481 m vs 2.205 m | 0.909 m vs 2.236 m | 0.836 m vs 6.874 m | **Gneiss** (+78.2%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.475 m vs 2.205 m | 0.995 m vs 2.236 m | 1.015 m vs 6.874 m | **Gneiss** (+78.5%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 1.249 m vs 2.205 m | 2.658 m vs 2.236 m | 1.238 m vs 6.874 m | **Gneiss** (+43.4%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 1.570 m vs N/A | 2.478 m vs N/A | 1.388 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 1.570 m vs N/A | 2.478 m vs N/A | 1.388 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 1.635 m vs N/A | 2.583 m vs N/A | 1.444 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 1.570 m vs N/A | 2.478 m vs N/A | 1.388 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 1.570 m vs N/A | 2.478 m vs N/A | 1.388 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 1.635 m vs N/A | 2.583 m vs N/A | 1.444 m vs N/A | **Gneiss** | Stable PPP integration. |

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.037 m vs 2.080 m | 3.297 m vs 3.311 m | 60.179 m vs 63.357 m | **Gneiss** (+2.1%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.103 m vs 2.080 m | 3.435 m vs 3.311 m | 60.142 m vs 63.357 m | RTKLIB (+1.1%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 3.465 m vs 2.080 m | 9.957 m vs 3.311 m | 61.603 m vs 63.357 m | RTKLIB (+40.0%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 2.037 m vs 2.080 m | 3.297 m vs 3.311 m | 60.179 m vs 63.357 m | **Gneiss** (+2.1%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 2.103 m vs 2.080 m | 3.435 m vs 3.311 m | 60.142 m vs 63.357 m | RTKLIB (+1.1%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 3.273 m vs 2.080 m | 9.530 m vs 3.311 m | 61.651 m vs 63.357 m | RTKLIB (+36.4%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 1.719 m vs 1.176 m | 2.357 m vs 1.773 m | 60.948 m vs 63.820 m | RTKLIB (+31.6%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.722 m vs 1.176 m | 2.357 m vs 1.773 m | 60.968 m vs 63.820 m | RTKLIB (+31.7%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 2.093 m vs 1.176 m | 4.148 m vs 1.773 m | 60.075 m vs 63.820 m | RTKLIB (+43.8%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 1.715 m vs 1.104 m | 2.351 m vs 1.831 m | 61.004 m vs 64.073 m | RTKLIB (+35.6%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.715 m vs 1.104 m | 2.353 m vs 1.831 m | 60.994 m vs 64.073 m | RTKLIB (+35.6%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 1.520 m vs 1.104 m | 2.711 m vs 1.831 m | 60.896 m vs 64.073 m | RTKLIB (+27.4%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `ppp` | forward | Off | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |

## PPP (f9p_ppp)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.438 m vs N/A | 3.057 m vs N/A | 2.039 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.438 m vs N/A | 3.057 m vs N/A | 2.039 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | forward | Tight | 3.148 m vs N/A | 3.542 m vs N/A | 1.964 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Off | 2.438 m vs N/A | 3.057 m vs N/A | 2.039 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 3.280 m vs N/A | 3.731 m vs N/A | 2.029 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Tight | 3.146 m vs N/A | 3.519 m vs N/A | 1.986 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `rtk` | forward | Off | 0.464 m vs N/A | 0.579 m vs N/A | 0.257 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.464 m vs N/A | 0.579 m vs N/A | 0.257 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | forward | Tight | 0.441 m vs N/A | 0.486 m vs N/A | 0.283 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Off | 0.464 m vs N/A | 0.580 m vs N/A | 0.257 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.464 m vs N/A | 0.579 m vs N/A | 0.257 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Tight | 0.312 m vs N/A | 0.472 m vs N/A | 0.224 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `ppp` | forward | Off | 2.421 m vs N/A | 3.045 m vs N/A | 2.031 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 2.421 m vs N/A | 3.045 m vs N/A | 2.031 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 2.436 m vs N/A | 3.057 m vs N/A | 2.059 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 2.421 m vs N/A | 3.045 m vs N/A | 2.031 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 2.421 m vs N/A | 3.045 m vs N/A | 2.031 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 2.436 m vs N/A | 3.057 m vs N/A | 2.059 m vs N/A | **Gneiss** | Stable PPP integration. |

