# Gneiss Comprehensive Benchmarks

This document empirically maps the performance of Gneiss across varying modes and datasets.

## GSDC (Pixel 4)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 3.297 m | 8.784 m | 62.750 m |
| `spp-ins` | 72.379 m | 6647.208 m | 73.009 m |
| `spp-ins-loosely-coupled` | 3.435 m | 9.143 m | 62.693 m |
| `rtk` | 2.357 m | 10.178 m | 62.642 m |
| `rtk-ins` | 3.668 m | 176.846 m | 63.389 m |
| `rtk-ins-loosely-coupled` | 2.357 m | 9.860 m | 62.646 m |
| `ppp` | 11.025 m | 11.025 m | 43.205 m |
| `ppp-ins` | 11.025 m | 11.025 m | 43.205 m |

## Shinjuku (UrbanNav)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 428.554 m | 643.099 m | 29.389 m |
| `spp-ins` | 58.843 m | 685.057 m | 25.849 m |
| `spp-ins-loosely-coupled` | 3.374 m | 25.769 m | 5.119 m |
| `rtk` | 1.789 m | 24.058 m | 2.696 m |
| `rtk-ins` | 3.862 m | 20.961 m | 3.179 m |
| `rtk-ins-loosely-coupled` | 1.790 m | 23.881 m | 2.725 m |
| `ppp` | 1.978 m | 24.181 m | 4.636 m |
| `ppp-ins` | 1.978 m | 24.181 m | 4.636 m |

## PPP (f9p_ppp)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 3.057 m | 3.932 m | 3.859 m |
| `spp-ins` | 3.790 m | 1816.545 m | 3.993 m |
| `spp-ins-loosely-coupled` | 3.057 m | 3.932 m | 3.859 m |
| `rtk` | 0.579 m | 1.319 m | 0.440 m |
| `rtk-ins` | 0.579 m | 1.319 m | 0.440 m |
| `rtk-ins-loosely-coupled` | 0.579 m | 1.319 m | 0.440 m |
| `ppp` | 3.015 m | 3.955 m | 4.016 m |
| `ppp-ins` | 3.015 m | 3.955 m | 4.016 m |

