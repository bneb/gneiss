# Gneiss Comprehensive Benchmarks

This document empirically maps the performance of Gneiss across varying modes and datasets.

## GSDC (Pixel 4)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 3.327 m | 8.648 m | 60.507 m |
| `spp-ins` | 7.887 m | 75.339 m | 60.687 m |
| `spp-ins-loosely-coupled` | 3.431 m | 9.157 m | 60.499 m |
| `rtk` | 2.231 m | 9.799 m | 68.395 m |
| `rtk-ins` | 4.477 m | 18.638 m | 64.082 m |
| `rtk-ins-loosely-coupled` | 2.228 m | 8.528 m | 68.384 m |
| `ppp` | 616.378 m | 2988.980 m | 930.022 m |
| `ppp-ins` | 616.378 m | 2988.980 m | 930.022 m |

## Shinjuku (UrbanNav)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 2.257 m | 26.551 m | 7.037 m |
| `spp-ins` | 7.826 m | 29.102 m | 7.773 m |
| `spp-ins-loosely-coupled` | 3.316 m | 25.678 m | 7.084 m |
| `rtk` | 1.453 m | 20.779 m | 2.752 m |
| `rtk-ins` | 29.888 m | 106.262 m | 8.587 m |
| `rtk-ins-loosely-coupled` | 1.453 m | 22.123 m | 2.791 m |
| `ppp` | 17.582 m | 39.608 m | 30.660 m |
| `ppp-ins` | 16.623 m | 41.207 m | 51.406 m |

## PPP (f9p_ppp)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 3.010 m | 3.844 m | 1.839 m |
| `spp-ins` | 3.533 m | 2185.700 m | 1.821 m |
| `spp-ins-loosely-coupled` | 3.010 m | 3.844 m | 1.839 m |
| `rtk` | Process Failed | Process Failed | Process Failed |
| `rtk-ins` | Process Failed | Process Failed | Process Failed |
| `rtk-ins-loosely-coupled` | 0.322 m | 0.503 m | 0.210 m |
| `ppp` | 6.802 m | 8.367 m | 0.926 m |
| `ppp-ins` | 6.802 m | 8.367 m | 0.926 m |

