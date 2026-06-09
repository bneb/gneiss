# Gneiss vs RTKLIB (demo5) Comparison

Side-by-side comparison on identical datasets and truth references.

## PPP (f9p_ppp)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | Failed |  |  |  |
| SPP | Gneiss | 3.057 m | 3.932 m | 3.859 m | |
| RTK Kinematic | RTKLIB | Failed |  |  |  |
| RTK Kinematic | Gneiss | 0.430 m | 0.715 m | 0.285 m | |
| RTK Kinematic (combined) | RTKLIB | Failed |  |  |  |
| RTK Kinematic (combined) | Gneiss | Failed |  |  | |

## GSDC (Pixel 4)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | 3.311 m | 10.191 m | 66.316 m | **Gneiss** ✅ |
| SPP | Gneiss | 3.297 m | 8.784 m | 62.750 m | |
| RTK Kinematic | RTKLIB | 1.773 m | 4.161 m | 64.598 m | RTKLIB ⚠️ |
| RTK Kinematic | Gneiss | 3.331 m | 9.437 m | 62.727 m | |
| RTK Kinematic (combined) | RTKLIB | 1.831 m | 3.126 m | 64.471 m |  |
| RTK Kinematic (combined) | Gneiss | Failed |  |  | |

## Shinjuku (u-blox)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | 13.209 m | 55.582 m | 38.538 m | **Gneiss** ✅ |
| SPP | Gneiss | 2.288 m | 26.576 m | 5.166 m | |
| RTK Kinematic | RTKLIB | 1.670 m | 11.843 m | 5.394 m |  |
| RTK Kinematic | Gneiss | Failed |  |  | |
| RTK Kinematic (combined) | RTKLIB | 2.089 m | 20.154 m | 5.517 m |  |
| RTK Kinematic (combined) | Gneiss | Failed |  |  | |

## Odaiba (u-blox)

| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-------|:--------|:--------|:--------|:-------|
| SPP | RTKLIB | 5.335 m | 35.568 m | 8.330 m | **Gneiss** ✅ |
| SPP | Gneiss | 2.409 m | 9.384 m | 2.853 m | |
| RTK Kinematic | RTKLIB | 2.239 m | 8.099 m | 6.413 m | **Gneiss** ✅ |
| RTK Kinematic | Gneiss | 1.252 m | 5.933 m | 2.807 m | |
| RTK Kinematic (combined) | RTKLIB | 2.236 m | 13.146 m | 8.219 m |  |
| RTK Kinematic (combined) | Gneiss | Failed |  |  | |

