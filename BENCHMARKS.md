# Gneiss Comprehensive Benchmarks

This document empirically maps the performance of Gneiss across varying modes and datasets.

## GSDC (Pixel 4)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 2.066 m | 3.327 m | 57.906 m |
| `spp-ins` | 3.967 m | 7.887 m | 59.000 m |
| `spp-ins-loosely-coupled` | 2.153 m | 3.431 m | 57.892 m |
| `rtk` | 1.443 m | 1.926 m | 61.764 m |
| `rtk-ins` | 2.502 m | 4.543 m | 60.909 m |
| `rtk-ins-loosely-coupled` | 1.455 m | 1.928 m | 61.702 m |
| `ppp` | 108.494 m | 616.378 m | 198.920 m |
| `ppp-ins` | 108.494 m | 616.378 m | 198.920 m |

## Shinjuku (UrbanNav)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 1.302 m | 2.257 m | 4.002 m |
| `spp-ins` | 4.006 m | 8.073 m | 4.540 m |
| `spp-ins-loosely-coupled` | 1.843 m | 3.316 m | 4.144 m |
| `rtk` | 1.179 m | 1.634 m | 0.658 m |
| `rtk-ins` | 1.523 m | 4.570 m | 1.169 m |
| `rtk-ins-loosely-coupled` | 1.184 m | 1.635 m | 0.643 m |
| `ppp` | 16.389 m | 23.241 m | 9.934 m |
| `ppp-ins` | 14.325 m | 20.532 m | 10.885 m |

## PPP (f9p_ppp)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 2.422 m | 3.010 m | 1.029 m |
| `spp-ins` | 2.924 m | 3.378 m | 0.968 m |
| `spp-ins-loosely-coupled` | 2.422 m | 3.010 m | 1.029 m |
| `rtk` | 0.012 m | 0.212 m | 0.022 m |
| `rtk-ins` | Process Failed | Process Failed | Process Failed |
| `rtk-ins-loosely-coupled` | 0.248 m | 0.317 m | 0.126 m |
| `ppp` | 5.316 m | 6.783 m | 0.701 m |
| `ppp-ins` | 5.316 m | 6.783 m | 0.701 m |

