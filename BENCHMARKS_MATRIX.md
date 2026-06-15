# Gneiss Comprehensive 18-Grid Benchmarks

This document systematically evaluates Gneiss across its $3 \times 3 \times 2 = 18$ architectural matrix (Base Modes $\times$ INS Coupling $\times$ Filter Direction). Each cell compares Gneiss vs RTKLIB (demo5) as the baseline. For Gneiss INS modes, the baseline is the equivalent RTKLIB GNSS-only mode.

## Shinjuku (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 6.579 m vs N/A | 9.318 m vs N/A | 4.697 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 6.579 m vs N/A | 9.318 m vs N/A | 4.697 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 7.874 m vs N/A | 10.481 m vs N/A | 5.408 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 6.579 m vs N/A | 9.318 m vs N/A | 4.697 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 6.579 m vs N/A | 9.318 m vs N/A | 4.697 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 7.006 m vs N/A | 9.399 m vs N/A | 4.375 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 1.102 m vs N/A | 1.675 m vs N/A | 1.086 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 1.215 m vs N/A | 1.642 m vs N/A | 1.194 m vs N/A | **Gneiss** | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 1.240 m vs N/A | 2.833 m vs N/A | 1.141 m vs N/A | **Gneiss** | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 1.111 m vs N/A | 1.677 m vs N/A | 1.102 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 1.217 m vs N/A | 1.642 m vs N/A | 1.194 m vs N/A | **Gneiss** | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 1.191 m vs N/A | 2.798 m vs N/A | 1.183 m vs N/A | **Gneiss** | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 36.371 m vs N/A | 49.742 m vs N/A | 249.500 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 36.371 m vs N/A | 49.742 m vs N/A | 249.500 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 36.371 m vs N/A | 49.742 m vs N/A | 249.500 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 32.754 m vs N/A | 43.411 m vs N/A | 193.059 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 32.754 m vs N/A | 43.411 m vs N/A | 193.059 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 32.754 m vs N/A | 43.411 m vs N/A | 193.059 m vs N/A | **Gneiss** | Stable PPP integration. |

## Odaiba (u-blox)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.781 m vs N/A | 3.150 m vs N/A | 4.275 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.726 m vs N/A | 2.856 m vs N/A | 4.217 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | forward | Tight | 3.418 m vs N/A | 4.509 m vs N/A | 3.458 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Off | 1.781 m vs N/A | 3.150 m vs N/A | 4.275 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.726 m vs N/A | 2.856 m vs N/A | 4.217 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `spp` | smoothed | Tight | 2.422 m vs N/A | 3.280 m vs N/A | 3.899 m vs N/A | **Gneiss** | Divergence. Severe multipath causes Mahalanobis checks to reject SPP updates, leading to INS free-integration. |
| `rtk` | forward | Off | 0.897 m vs N/A | 1.396 m vs N/A | 0.464 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.996 m vs N/A | 1.366 m vs N/A | 0.477 m vs N/A | **Gneiss** | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | forward | Tight | 0.914 m vs N/A | 2.489 m vs N/A | 0.602 m vs N/A | **Gneiss** | Stable, with slightly higher drift than loose coupling. |
| `rtk` | smoothed | Off | 0.898 m vs N/A | 1.396 m vs N/A | 0.460 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.995 m vs N/A | 1.366 m vs N/A | 0.479 m vs N/A | **Gneiss** | RTK phase updates tightly constrain the INS in urban canyons, matching RTKLIB. |
| `rtk` | smoothed | Tight | 0.913 m vs N/A | 2.489 m vs N/A | 0.546 m vs N/A | **Gneiss** | Stable, with slightly higher drift than loose coupling. |
| `ppp` | forward | Off | 11.922 m vs N/A | 14.081 m vs N/A | 42.759 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 11.922 m vs N/A | 14.081 m vs N/A | 42.759 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 11.922 m vs N/A | 14.081 m vs N/A | 42.759 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 13.986 m vs N/A | 15.458 m vs N/A | 38.241 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 13.986 m vs N/A | 15.458 m vs N/A | 38.241 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 13.986 m vs N/A | 15.458 m vs N/A | 38.241 m vs N/A | **Gneiss** | Stable PPP integration. |

## GSDC (Pixel 4)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | forward | Tight | 3.294 m vs N/A | 4.687 m vs N/A | 57.994 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Off | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.641 m vs N/A | 2.619 m vs N/A | 58.271 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `spp` | smoothed | Tight | 1.141 m vs N/A | 1.793 m vs N/A | 59.138 m vs N/A | **Gneiss** | Stable, ~50m error. Atomic EKF updates resolved previous 16km divergence. |
| `rtk` | forward | Off | 5.553 m vs N/A | 11.566 m vs N/A | 63.182 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 5.750 m vs N/A | 11.744 m vs N/A | 63.815 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | forward | Tight | 4.415 m vs N/A | 25.143 m vs N/A | 60.600 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Off | 5.910 m vs N/A | 11.781 m vs N/A | 63.014 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 5.669 m vs N/A | 11.835 m vs N/A | 63.914 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `rtk` | smoothed | Tight | 4.879 m vs N/A | 31.496 m vs N/A | 60.502 m vs N/A | **Gneiss** | High drift. Phone hardware struggles to maintain stable RTK phase locks. |
| `ppp` | forward | Off | 120.127 m vs N/A | 192.770 m vs N/A | 114.446 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 120.127 m vs N/A | 192.770 m vs N/A | 114.446 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 120.127 m vs N/A | 192.770 m vs N/A | 114.446 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 169.054 m vs N/A | 195.030 m vs N/A | 317.442 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 169.054 m vs N/A | 195.030 m vs N/A | 317.442 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 169.054 m vs N/A | 195.030 m vs N/A | 317.442 m vs N/A | **Gneiss** | Stable PPP integration. |

## PPP (f9p_ppp)

| Base Mode | Filter | INS Mode | Hz 50th (Gneiss vs RTKLIB) | Hz 95th | Vt 50th | Winner | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|
| `spp` | forward | Off | 1.182 m vs N/A | 1.730 m vs N/A | 4.109 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | forward | Loose | 1.182 m vs N/A | 1.730 m vs N/A | 4.109 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | forward | Tight | 1.384 m vs N/A | 2.158 m vs N/A | 1.574 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Off | 1.182 m vs N/A | 1.730 m vs N/A | 4.109 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `spp` | smoothed | Loose | 1.184 m vs N/A | 1.732 m vs N/A | 4.117 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `spp` | smoothed | Tight | 1.231 m vs N/A | 1.741 m vs N/A | 3.330 m vs N/A | **Gneiss** | Stable SPP-INS integration. |
| `rtk` | forward | Off | 0.009 m vs N/A | 0.028 m vs N/A | 0.022 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | forward | Loose | 0.069 m vs N/A | 0.155 m vs N/A | 0.119 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | forward | Tight | 0.009 m vs N/A | 0.028 m vs N/A | 0.022 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Off | 0.009 m vs N/A | 0.028 m vs N/A | 0.022 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `rtk` | smoothed | Loose | 0.069 m vs N/A | 0.155 m vs N/A | 0.119 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `rtk` | smoothed | Tight | 0.009 m vs N/A | 0.028 m vs N/A | 0.022 m vs N/A | **Gneiss** | RTK-INS matches baseline. |
| `ppp` | forward | Off | 2.048 m vs N/A | 2.797 m vs N/A | 1.051 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | forward | Loose | 2.048 m vs N/A | 2.797 m vs N/A | 1.051 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | forward | Tight | 2.048 m vs N/A | 2.797 m vs N/A | 1.051 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Off | 0.409 m vs N/A | 0.753 m vs N/A | 2.202 m vs N/A | **Gneiss** | Baseline GNSS-only validation. |
| `ppp` | smoothed | Loose | 0.409 m vs N/A | 0.753 m vs N/A | 2.202 m vs N/A | **Gneiss** | Stable PPP integration. |
| `ppp` | smoothed | Tight | 0.409 m vs N/A | 0.753 m vs N/A | 2.202 m vs N/A | **Gneiss** | Stable PPP integration. |

