# Gneiss vs RTKLIB (demo5) Comparison

Side-by-side comparison on identical datasets and truth references.

## PPP (f9p_ppp)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | Failed |  |  |  |
| SPP | Gneiss | 3.010 m | 3.844 m | 1.839 m | |
| RTK Kinematic | RTKLIB | Failed |  |  |  |
| RTK Kinematic | Gneiss | 0.212 m | 0.502 m | 0.181 m | |
| RTK Kinematic (combined) | RTKLIB | Failed |  |  |  |
| RTK Kinematic (combined) | Gneiss | Failed |  |  | |
| PPP Kinematic | RTKLIB | Failed |  |  |  |
| PPP Kinematic | Gneiss | 6.783 m | 8.399 m | 0.919 m | |
| PPP Kinematic (combined) | RTKLIB | Failed |  |  |  |
| PPP Kinematic (combined) | Gneiss | 6.783 m | 8.399 m | 0.919 m | |

## GSDC (Pixel 4)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | 3.311 m | 10.191 m | 66.316 m | RTKLIB ⚠️ |
| SPP | Gneiss | 3.327 m | 8.648 m | 60.507 m | |
| RTK Kinematic | RTKLIB | 1.773 m | 4.161 m | 64.598 m | RTKLIB ⚠️ |
| RTK Kinematic | Gneiss | 2.031 m | 9.340 m | 63.252 m | |
| RTK Kinematic (combined) | RTKLIB | 1.831 m | 3.126 m | 64.471 m | RTKLIB ⚠️ |
| RTK Kinematic (combined) | Gneiss | 3.277 m | 349.678 m | 63.445 m | |
| PPP Kinematic | RTKLIB | Failed |  |  |  |
| PPP Kinematic | Gneiss | 616.378 m | 2988.980 m | 930.022 m | |
| PPP Kinematic (combined) | RTKLIB | Failed |  |  |  |
| PPP Kinematic (combined) | Gneiss | 616.378 m | 2988.980 m | 930.022 m | |

## Shinjuku (u-blox)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | 14.144 m | 67.742 m | 38.519 m | **Gneiss** ✅ |
| SPP | Gneiss | 2.257 m | 26.551 m | 7.037 m | |
| RTK Kinematic | RTKLIB | 5.730 m | 20.200 m | 9.326 m | **Gneiss** ✅ |
| RTK Kinematic | Gneiss | 1.525 m | 15.639 m | 1.772 m | |
| RTK Kinematic (combined) | RTKLIB | 5.978 m | 18.029 m | 9.112 m | **Gneiss** ✅ |
| RTK Kinematic (combined) | Gneiss | 1.521 m | 13.686 m | 1.768 m | |
| PPP Kinematic | RTKLIB | Failed |  |  |  |
| PPP Kinematic | Gneiss | 23.241 m | 57.248 m | 17.566 m | |
| PPP Kinematic (combined) | RTKLIB | Failed |  |  |  |
| PPP Kinematic (combined) | Gneiss | 23.241 m | 57.248 m | 17.566 m | |

## Odaiba (u-blox)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | 7.453 m | 48.906 m | 8.709 m | **Gneiss** ✅ |
| SPP | Gneiss | 2.533 m | 9.465 m | 2.901 m | |
| RTK Kinematic | RTKLIB | 4.014 m | 17.444 m | 5.744 m | **Gneiss** ✅ |
| RTK Kinematic | Gneiss | 1.203 m | 6.244 m | 1.784 m | |
| RTK Kinematic (combined) | RTKLIB | 3.886 m | 16.086 m | 6.134 m | **Gneiss** ✅ |
| RTK Kinematic (combined) | Gneiss | 1.201 m | 6.053 m | 1.739 m | |
| PPP Kinematic | RTKLIB | Failed |  |  |  |
| PPP Kinematic | Gneiss | 12.518 m | 21.905 m | 23.796 m | |
| PPP Kinematic (combined) | RTKLIB | Failed |  |  |  |
| PPP Kinematic (combined) | Gneiss | 12.518 m | 21.905 m | 23.796 m | |

