# Gneiss Comprehensive Benchmarks

This document empirically maps the performance of Gneiss across varying modes and datasets.

## GSDC (Pixel 4)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 2.066 m | 3.327 m | 57.906 m |
| `spp-ins` | 3.163 m | 4.594 m | 57.179 m |
| `spp-ins-loosely-coupled` | 2.153 m | 3.431 m | 57.892 m |
| `rtk` | 1.437 m | 1.935 m | 62.038 m |
| `rtk-ins` | 2.441 m | 4.459 m | 61.264 m |
| `rtk-ins-loosely-coupled` | 1.446 m | 1.948 m | 62.062 m |
| `ppp` | 108.494 m | 616.378 m | 198.920 m |
| `ppp-ins` | 108.494 m | 616.378 m | 198.920 m |

## Shinjuku (UrbanNav)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 1.302 m | 2.257 m | 4.002 m |
| `spp-ins` | 1.431 m | 2.401 m | 3.807 m |
| `spp-ins-loosely-coupled` | 1.407 m | 2.437 m | 4.007 m |
| `rtk` | 1.158 m | 1.683 m | 0.496 m |
| `rtk-ins` | 27.731 m | 47.789 m | 8.912 m |
| `rtk-ins-loosely-coupled` | 1.178 m | 1.697 m | 0.499 m |
| `ppp` | 16.389 m | 23.241 m | 9.934 m |
| `ppp-ins` | 14.325 m | 20.532 m | 10.885 m |

## PPP (f9p_ppp)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 2.422 m | 3.010 m | 1.029 m |
| `spp-ins` | 2.534 m | 2.863 m | 1.416 m |
| `spp-ins-loosely-coupled` | 2.422 m | 3.010 m | 1.029 m |
| `rtk` | 0.009 m | 0.020 m | 0.016 m |
| `rtk-ins` | 0.009 m | 0.024 m | 0.013 m |
| `rtk-ins-loosely-coupled` | 0.216 m | 0.300 m | 0.104 m |
| `ppp` | 5.316 m | 6.783 m | 0.701 m |
| `ppp-ins` | 5.316 m | 6.783 m | 0.701 m |

