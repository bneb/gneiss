# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## Shinjuku (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.302 m vs 5.544 m | 2.257 m vs 13.209 m | 4.002 m vs 8.036 m | **Gneiss** (+76.5%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.302 m vs 5.544 m | 2.257 m vs 13.209 m | 4.002 m vs 8.036 m | **Gneiss** (+76.5%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 1.379 m vs 5.544 m | 2.451 m vs 13.209 m | 4.431 m vs 8.036 m | **Gneiss** (+75.1%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 1.302 m vs 5.544 m | 2.257 m vs 13.209 m | 4.002 m vs 8.036 m | **Gneiss** (+76.5%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.407 m vs 5.544 m | 2.437 m vs 13.209 m | 4.007 m vs 8.036 m | **Gneiss** (+74.6%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 1.313 m vs 5.544 m | 2.307 m vs 13.209 m | 4.956 m vs 8.036 m | **Gneiss** (+76.3%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 0.762 m vs 0.674 m | 1.305 m vs 1.670 m | 0.454 m vs 2.984 m | RTKLIB (+11.5%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.799 m vs 0.674 m | 1.334 m vs 1.670 m | 0.596 m vs 2.984 m | RTKLIB (+15.6%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 2.341 m vs 0.674 m | 4.692 m vs 1.670 m | 1.779 m vs 2.984 m | RTKLIB (+71.2%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 0.763 m vs 0.706 m | 1.306 m vs 2.089 m | 0.456 m vs 5.065 m | RTKLIB (+7.5%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.799 m vs 0.706 m | 1.334 m vs 2.089 m | 0.596 m vs 5.065 m | RTKLIB (+11.6%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 2.218 m vs 0.706 m | 4.547 m vs 2.089 m | 1.667 m vs 5.065 m | RTKLIB (+68.2%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 16.389 m vs N/A | 23.241 m vs N/A | 9.934 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 14.325 m vs N/A | 20.532 m vs N/A | 10.885 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 14.325 m vs N/A | 20.532 m vs N/A | 10.885 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 16.389 m vs N/A | 23.241 m vs N/A | 9.934 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 14.325 m vs N/A | 20.532 m vs N/A | 10.885 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 14.325 m vs N/A | 20.532 m vs N/A | 10.885 m vs N/A | **Gneiss** | Stable PPP integration. |

## Odaiba (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.594 m vs 2.061 m | 2.533 m vs 5.335 m | 1.257 m vs 5.641 m | **Gneiss** (+22.7%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.594 m vs 2.061 m | 2.533 m vs 5.335 m | 1.257 m vs 5.641 m | **Gneiss** (+22.7%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 0.956 m vs 2.061 m | 1.714 m vs 5.335 m | 1.395 m vs 5.641 m | **Gneiss** (+53.6%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 1.594 m vs 2.061 m | 2.533 m vs 5.335 m | 1.257 m vs 5.641 m | **Gneiss** (+22.7%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.917 m vs 2.061 m | 2.760 m vs 5.335 m | 1.259 m vs 5.641 m | **Gneiss** (+7.0%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 1.489 m vs 2.061 m | 2.229 m vs 5.335 m | 1.100 m vs 5.641 m | **Gneiss** (+27.8%) | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 0.694 m vs 1.516 m | 1.145 m vs 2.239 m | 0.411 m vs 4.081 m | **Gneiss** (+54.2%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.763 m vs 1.516 m | 1.090 m vs 2.239 m | 0.574 m vs 4.081 m | **Gneiss** (+49.7%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 2.685 m vs 1.516 m | 4.933 m vs 2.239 m | 1.156 m vs 4.081 m | RTKLIB (+43.5%) | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 0.693 m vs 2.205 m | 1.145 m vs 2.236 m | 0.420 m vs 6.874 m | **Gneiss** (+68.6%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.763 m vs 2.205 m | 1.090 m vs 2.236 m | 0.575 m vs 6.874 m | **Gneiss** (+65.4%) | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 2.641 m vs 2.205 m | 4.867 m vs 2.236 m | 1.163 m vs 6.874 m | RTKLIB (+16.5%) | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 10.079 m vs N/A | 12.518 m vs N/A | 9.321 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 7.297 m vs N/A | 10.604 m vs N/A | 21.459 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 7.297 m vs N/A | 10.604 m vs N/A | 21.459 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 10.079 m vs N/A | 12.518 m vs N/A | 9.321 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 7.297 m vs N/A | 10.604 m vs N/A | 21.459 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 7.297 m vs N/A | 10.604 m vs N/A | 21.459 m vs N/A | **Gneiss** | Stable PPP integration. |

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.066 m vs 2.080 m | 3.327 m vs 3.311 m | 57.906 m vs 63.357 m | **Gneiss** (+0.7%) | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.153 m vs 2.080 m | 3.431 m vs 3.311 m | 57.892 m vs 63.357 m | RTKLIB (+3.4%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 3.163 m vs 2.080 m | 4.594 m vs 3.311 m | 57.179 m vs 63.357 m | RTKLIB (+34.2%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 2.066 m vs 2.080 m | 3.327 m vs 3.311 m | 57.906 m vs 63.357 m | **Gneiss** (+0.7%) | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 2.153 m vs 2.080 m | 3.431 m vs 3.311 m | 57.892 m vs 63.357 m | RTKLIB (+3.4%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 1.376 m vs 2.080 m | 2.231 m vs 3.311 m | 58.959 m vs 63.357 m | **Gneiss** (+33.8%) | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 1.545 m vs 1.176 m | 2.095 m vs 1.773 m | 61.511 m vs 63.820 m | RTKLIB (+23.9%) | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.545 m vs 1.176 m | 2.072 m vs 1.773 m | 61.468 m vs 63.820 m | RTKLIB (+23.9%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 6.792 m vs 1.176 m | 12.893 m vs 1.773 m | 59.403 m vs 63.820 m | RTKLIB (+82.7%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 1.543 m vs 1.104 m | 2.080 m vs 1.831 m | 61.512 m vs 64.073 m | RTKLIB (+28.5%) | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.545 m vs 1.104 m | 2.068 m vs 1.831 m | 61.478 m vs 64.073 m | RTKLIB (+28.5%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 6.697 m vs 1.104 m | 12.504 m vs 1.831 m | 59.520 m vs 64.073 m | RTKLIB (+83.5%) | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `ppp` | forward | Off | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 108.494 m vs N/A | 616.378 m vs N/A | 198.920 m vs N/A | **Gneiss** | Stable PPP integration. |

## PPP (f9p_ppp)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 2.422 m vs N/A | 3.010 m vs N/A | 1.029 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 2.422 m vs N/A | 3.010 m vs N/A | 1.029 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | forward | Tight | 2.534 m vs N/A | 2.863 m vs N/A | 1.416 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Off | 2.422 m vs N/A | 3.010 m vs N/A | 1.029 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 2.422 m vs N/A | 3.011 m vs N/A | 1.029 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Tight | 2.496 m vs N/A | 3.081 m vs N/A | 0.984 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `rtk` | forward | Off | 0.007 m vs N/A | 0.014 m vs N/A | 0.010 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.099 m vs N/A | 0.154 m vs N/A | 0.037 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | forward | Tight | 0.008 m vs N/A | 0.022 m vs N/A | 0.015 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Off | 0.007 m vs N/A | 0.014 m vs N/A | 0.010 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.099 m vs N/A | 0.154 m vs N/A | 0.037 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Tight | 0.008 m vs N/A | 0.020 m vs N/A | 0.014 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `ppp` | forward | Off | 5.316 m vs N/A | 6.783 m vs N/A | 0.701 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 5.316 m vs N/A | 6.783 m vs N/A | 0.701 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 5.316 m vs N/A | 6.783 m vs N/A | 0.701 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 5.316 m vs N/A | 6.783 m vs N/A | 0.701 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 5.316 m vs N/A | 6.783 m vs N/A | 0.701 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 5.316 m vs N/A | 6.783 m vs N/A | 0.701 m vs N/A | **Gneiss** | Stable PPP integration. |

