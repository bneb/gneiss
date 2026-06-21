# Gneiss vs RTKLIB (demo5) Comparison

Last updated: 2026-06-21. Rows marked [stale] > 1 week old.

| Dataset | Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner | Freshness |
|:-----|:-----|:-------|:--------|:--------|:--------|:-------|
| GSDC (Pixel 4) | SPP | RTKLIB | 2.080 m | 3.311 m | 63.357 m | **Gneiss** | 🟡 [stale] |
| GSDC (Pixel 4) | SPP | Gneiss | 2.037 m | 3.297 m | 60.179 m | | 🟡 [stale] |
| GSDC (Pixel 4) | RTK Kinematic | RTKLIB | 1.176 m | 1.773 m | 63.820 m | RTKLIB | 🟡 [stale] |
| GSDC (Pixel 4) | RTK Kinematic | Gneiss | 8.365 m | 127.464 m | 87.220 m | | 🟡 [stale] |
| GSDC (Pixel 4) | RTK Kinematic (combined) | RTKLIB | 1.104 m | 1.831 m | 64.073 m | RTKLIB | 🟡 [stale] |
| GSDC (Pixel 4) | RTK Kinematic (combined) | Gneiss | 8.365 m | 127.464 m | 87.220 m | | 🟡 [stale] |
| GSDC (Pixel 4) | PPP Kinematic (EKF) | RTKLIB | 2.326 m | 4.253 m | 58.704 m | RTKLIB | 🟡 [stale] |
| GSDC (Pixel 4) | PPP Kinematic (EKF) | Gneiss | 147.223 m | 218.761 m | 30.284 m | | 🟡 [stale] |
| GSDC (Pixel 4) | PPP Kinematic (FG) | RTKLIB | 2.326 m | 4.253 m | 58.704 m | RTKLIB | 🟡 [stale] |
| GSDC (Pixel 4) | PPP Kinematic (FG) | Gneiss | 178.824 m | 195.867 m | 30.284 m | | 🟡 [stale] |
| Shinjuku (UrbanNav) | SPP | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) | 🟢 2026-06-19 | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | SPP | Gneiss | 1.849 m | 3.167 m | 2.446 m | | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | RTK Kinematic | RTKLIB | 2.205 m | 5.730 m | 4.115 m | | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | RTK Kinematic | Gneiss | 1.547 m | 2.603 m | 3.470 m | **Gneiss** | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | RTK Kinematic (combined) | RTKLIB | 2.809 m | 5.978 m | 4.700 m | RTKLIB | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | RTK Kinematic (combined) | Gneiss | 6.036 m | 14.981 m | 8.143 m | | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | PPP Kinematic (EKF) | RTKLIB | 2.003 m | 3.749 m | 3.235 m | RTKLIB | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | PPP Kinematic (EKF) | Gneiss | 30.267 m | 54.412 m | 223.067 m | [deprecated] | 🔴 2026-06-20 |
| Shinjuku (UrbanNav) | PPP Kinematic (IEKF) | RTKLIB | 2.003 m | 3.749 m | 3.235 m | RTKLIB | 🟢 2026-06-19 |
| Shinjuku (UrbanNav) | PPP Kinematic (IEKF) | Gneiss | 5.015 m | 7.402 m | 13.725 m | [stale] | 🔴 2026-06-19 |
| Shinjuku (UrbanNav) | PPP Kinematic (IEKF) fwd | Gneiss | 7.542 m | 41.576 m | 13.891 m | — | 🟢 2026-06-20 |
| Shinjuku (UrbanNav) | PPP Kinematic (IEKF) smooth | Gneiss | **16.661 m** | **42.556 m** | **16.090 m** | — | 🔴 2026-06-21 |
| Odaiba (UrbanNav) | SPP | RTKLIB | 2.799 m | 7.453 m | 4.885 m | **Gneiss** | 🟡 [stale] |
| Odaiba (UrbanNav) | SPP | Gneiss | 2.099 m | 2.958 m | 2.253 m | | 🟡 [stale] |
| Odaiba (UrbanNav) | RTK Kinematic | RTKLIB | 2.867 m | 4.014 m | 2.490 m | | 🟡 [stale] |
| Odaiba (UrbanNav) | RTK Kinematic | Gneiss | 1.257 m | 1.923 m | 2.569 m | **Gneiss** | 🟢 2026-06-19 |
| Odaiba (UrbanNav) | RTK Kinematic (combined) | RTKLIB | 3.031 m | 3.886 m | 4.168 m | RTKLIB | 🟡 [stale] |
| Odaiba (UrbanNav) | RTK Kinematic (combined) | Gneiss | 5.664 m | 17.037 m | 3.069 m | | 🟡 [stale] |
| Odaiba (UrbanNav) | PPP Kinematic (EKF) | RTKLIB | 3.965 m | 4.691 m | 6.089 m | RTKLIB | 🟡 [stale] |
| Odaiba (UrbanNav) | PPP Kinematic (EKF) | Gneiss | 16.549 m | 21.594 m | 14.939 m | [deprecated] | 🔴 2026-06-20 |
| Odaiba (UrbanNav) | PPP Kinematic (IEKF) | RTKLIB | 3.965 m | 4.691 m | 6.089 m | — | 🟢 2026-06-19 |
| Odaiba (UrbanNav) | PPP Kinematic (IEKF) | Gneiss | 3.404 m | 4.961 m | 5.902 m | [stale] | 🟢 2026-06-19 |
| Odaiba (UrbanNav) | PPP Kinematic (IEKF) fwd | Gneiss | 5.288 m | 20.668 m | 5.833 m | — | 🟢 2026-06-20 |
| Odaiba (UrbanNav) | PPP Kinematic (IEKF) smooth | Gneiss | 5.221 m | 35.823 m | 17.526 m | [stale, buggy smoother] | 🔴 2026-06-20 |
| Odaiba (UrbanNav) | PPP Kinematic (IEKF) smooth | Gneiss | **6.423 m** | **16.078 m** | **6.056 m** | — | 🟢 2026-06-21 |
| PPP (f9p_ppp) | SPP | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) | 🟡 [stale] |
| PPP (f9p_ppp) | SPP | Gneiss | 1.367 m | 2.008 m | 0.934 m | | 🟡 [stale] |
| PPP (f9p_ppp) | RTK Kinematic | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) | 🟡 [stale] |
| PPP (f9p_ppp) | RTK Kinematic | Gneiss | 0.246 m | 0.402 m | 0.168 m | | 🟡 [stale] |
| PPP (f9p_ppp) | RTK Kinematic (combined) | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) | 🟡 [stale] |
| PPP (f9p_ppp) | RTK Kinematic (combined) | Gneiss | 0.246 m | 0.402 m | 0.168 m | | 🟡 [stale] |
| PPP (f9p_ppp) | PPP Kinematic (EKF) | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) | 🟡 [stale] |
| PPP (f9p_ppp) | PPP Kinematic (EKF) | Gneiss | 13.137 m | 16.818 m | 16.127 m | | 🟡 [stale] |
| PPP (f9p_ppp) | PPP Kinematic (FG) | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) | 🟡 [stale] |
| PPP (f9p_ppp) | PPP Kinematic (FG) | Gneiss | 3.516 m | 5.081 m | 1.263 m | | 🟡 [stale] |
| PPP (f9p_ppp) | PPP Kinematic (IEKF) bcast | Gneiss | 1.930 m | 3.881 m | 7.275 m | — | 🟢 2026-06-20 |
| PPP (f9p_ppp) | PPP Kinematic (IEKF) precise | Gneiss | 1.492 m | 2.611 m | 5.877 m | — | 🟢 2026-06-20 |
| UrbanLoco (Example) | SPP | RTKLIB | No data |  |  |  | 🟡 [stale] |
| UrbanLoco (Example) | SPP | Gneiss | No data |  |  | | 🟡 [stale] |
| UrbanLoco (Example) | RTK Kinematic | RTKLIB | No data |  |  |  | 🟡 [stale] |
| UrbanLoco (Example) | RTK Kinematic | Gneiss | No data |  |  | | 🟡 [stale] |
| UrbanLoco (Example) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| UrbanLoco (Example) | RTK Kinematic (combined) | Gneiss | No data |  |  | | 🟡 [stale] |
| UrbanLoco (Example) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| UrbanLoco (Example) | PPP Kinematic (EKF) | Gneiss | No data |  |  | | 🟡 [stale] |
| UrbanLoco (Example) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| UrbanLoco (Example) | PPP Kinematic (FG) | Gneiss | No data |  |  | | 🟡 [stale] |
| TEX-CUP (UT Austin) | SPP | RTKLIB | No data |  |  |  | 🟡 [stale] |
| TEX-CUP (UT Austin) | SPP | Gneiss | No data |  |  | | 🟡 [stale] |
| TEX-CUP (UT Austin) | RTK Kinematic | RTKLIB | No data |  |  |  | 🟡 [stale] |
| TEX-CUP (UT Austin) | RTK Kinematic | Gneiss | No data |  |  | | 🟡 [stale] |
| TEX-CUP (UT Austin) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| TEX-CUP (UT Austin) | RTK Kinematic (combined) | Gneiss | No data |  |  | | 🟡 [stale] |
| TEX-CUP (UT Austin) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| TEX-CUP (UT Austin) | PPP Kinematic (EKF) | Gneiss | No data |  |  | | 🟡 [stale] |
| TEX-CUP (UT Austin) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| TEX-CUP (UT Austin) | PPP Kinematic (FG) | Gneiss | No data |  |  | | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | SPP | RTKLIB | No data |  |  |  | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | SPP | Gneiss | No data |  |  | | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | RTK Kinematic | RTKLIB | No data |  |  |  | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | RTK Kinematic | Gneiss | No data |  |  | | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | RTK Kinematic (combined) | Gneiss | No data |  |  | | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (EKF) | Gneiss | No data |  |  | | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (FG) | Gneiss | No data |  |  | | 🟡 [stale] |
| smartLoc (TU Chemnitz) | SPP | RTKLIB | No data |  |  |  | 🟡 [stale] |
| smartLoc (TU Chemnitz) | SPP | Gneiss | No data |  |  | | 🟡 [stale] |
| smartLoc (TU Chemnitz) | RTK Kinematic | RTKLIB | No data |  |  |  | 🟡 [stale] |
| smartLoc (TU Chemnitz) | RTK Kinematic | Gneiss | No data |  |  | | 🟡 [stale] |
| smartLoc (TU Chemnitz) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| smartLoc (TU Chemnitz) | RTK Kinematic (combined) | Gneiss | No data |  |  | | 🟡 [stale] |
| smartLoc (TU Chemnitz) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| smartLoc (TU Chemnitz) | PPP Kinematic (EKF) | Gneiss | No data |  |  | | 🟡 [stale] |
| smartLoc (TU Chemnitz) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  | 🟡 [stale] |
| smartLoc (TU Chemnitz) | PPP Kinematic (FG) | Gneiss | No data |  |  | | 🟡 [stale] |
