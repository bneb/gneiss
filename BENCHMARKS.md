# Gneiss Comprehensive Benchmarks

## GSDC (Pixel 4)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 1.967 m | 3.421 m | 57.829 m |
| `spp-ins` | 3.046 m | 4.412 m | 58.037 m |
| `rtk` | 47.453 m | 54.382 m | 92.086 m |
| `rtk-ins` | 47.453 m | 54.382 m | 92.086 m |

> [!WARNING]
> The RTK and RTK-INS metrics for GSDC exhibit a ~45m bias shift. This is because the official NGS CORS API base coordinates (`--base-coord-api`) are utilized, revealing that the GSDC provided ground truth was erroneously aligned to the un-surveyed `APPROX POSITION XYZ` in the original RINEX header.
| `ppp` | 151.025 m | 178.978 m | 89.220 m |
| `ppp-fg` | 154.789 m | 184.598 m | 216.244 m |
| `ppp-ins-fg` | 154.789 m | 184.598 m | 216.244 m |

## Shinjuku (UrbanNav)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 5.476 m | 7.814 m | 6.356 m |
| `spp-ins` | 5.643 m | 8.476 m | 7.329 m |
| `rtk` | 1.221 m | 1.906 m | 1.447 m |
| `rtk-ins` (manual lever arm) | 1.215 m | 2.002 m | 0.776 m |
| `rtk-ins` (auto-calibrated) | 2.521 m | 6.477 m | 1.118 m |
| `ppp` | 26.495 m | 33.367 m | 229.695 m |
| `ppp-fg` | 42.033 m | 62.378 m | 451.044 m |
| `ppp-ins-fg` | 42.033 m | 62.378 m | 451.044 m |

## Odaiba (UrbanNav)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 3.914 m | 5.360 m | 3.157 m |
| `spp-ins` | 5.184 m | 6.412 m | 2.644 m |
| `rtk` | 0.900 m | 1.340 m | 0.307 m |
| `rtk-ins` | 0.686 m | 1.320 m | 0.429 m |
| `ppp` | 12.853 m | 16.104 m | 4.435 m |
| `ppp-fg` | 12.362 m | 15.446 m | 4.222 m |
| `ppp-ins-fg` | 12.362 m | 15.446 m | 4.222 m |

## PPP (f9p_ppp)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | 1.066 m | 1.549 m | 3.785 m |
| `spp-ins` | 0.668 m | 1.088 m | 1.497 m |
| `rtk` | 0.010 m | 0.029 m | 0.017 m |
| `rtk-ins` | 0.010 m | 0.029 m | 0.017 m |
| `ppp` | 0.528 m | 0.654 m | 0.198 m |
| `ppp-fg` | 0.655 m | 2.154 m | 0.366 m |
| `ppp-ins-fg` | 0.655 m | 2.154 m | 0.366 m |

## UrbanLoco (Example)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | No Truth | No Truth | No Truth |
| `spp-ins` | No Truth | No Truth | No Truth |
| `rtk` | No Truth | No Truth | No Truth |
| `rtk-ins` | No Truth | No Truth | No Truth |
| `ppp` | No Truth | No Truth | No Truth |
| `ppp-fg` | No Truth | No Truth | No Truth |
| `ppp-ins-fg` | No Truth | No Truth | No Truth |

## TEX-CUP (UT Austin)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | No Truth | No Truth | No Truth |
| `spp-ins` | No Truth | No Truth | No Truth |
| `rtk` | No Truth | No Truth | No Truth |
| `rtk-ins` | No Truth | No Truth | No Truth |
| `ppp` | No Truth | No Truth | No Truth |
| `ppp-fg` | No Truth | No Truth | No Truth |
| `ppp-ins-fg` | No Truth | No Truth | No Truth |

## WHU-Smartphone (Xiaomi)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | No Truth | No Truth | No Truth |
| `spp-ins` | No Truth | No Truth | No Truth |
| `rtk` | No Truth | No Truth | No Truth |
| `rtk-ins` | No Truth | No Truth | No Truth |
| `ppp` | No Truth | No Truth | No Truth |
| `ppp-fg` | No Truth | No Truth | No Truth |
| `ppp-ins-fg` | No Truth | No Truth | No Truth |

## smartLoc (TU Chemnitz)

| Mode | Median Horizontal | 95% Horizontal | Median Vertical |
| :--- | :--- | :--- | :--- |
| `spp` | No Truth | No Truth | No Truth |
| `spp-ins` | No Truth | No Truth | No Truth |
| `rtk` | No Truth | No Truth | No Truth |
| `rtk-ins` | No Truth | No Truth | No Truth |
| `ppp` | No Truth | No Truth | No Truth |
| `ppp-fg` | No Truth | No Truth | No Truth |
| `ppp-ins-fg` | No Truth | No Truth | No Truth |

