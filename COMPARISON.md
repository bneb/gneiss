# Gneiss vs RTKLIB (demo5) Comparison

| Dataset | Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:-----|:-----|:-------|:--------|:--------|:--------|:-------|
| GSDC (Pixel 4) | SPP | RTKLIB | 2.080 m | 3.311 m | 63.357 m | **Gneiss** |
| GSDC (Pixel 4) | SPP | Gneiss | 2.037 m | 3.297 m | 60.179 m | |
| GSDC (Pixel 4) | RTK Kinematic | RTKLIB | 1.176 m | 1.773 m | 63.820 m | RTKLIB |
| GSDC (Pixel 4) | RTK Kinematic | Gneiss | 8.365 m | 127.464 m | 87.220 m | |
| GSDC (Pixel 4) | RTK Kinematic (combined) | RTKLIB | 1.104 m | 1.831 m | 64.073 m | RTKLIB |
| GSDC (Pixel 4) | RTK Kinematic (combined) | Gneiss | 8.365 m | 127.464 m | 87.220 m | |
| GSDC (Pixel 4) | PPP Kinematic (EKF) | RTKLIB | 2.326 m | 4.253 m | 58.704 m | RTKLIB |
| GSDC (Pixel 4) | PPP Kinematic (EKF) | Gneiss | 147.223 m | 218.761 m | 30.284 m | |
| GSDC (Pixel 4) | PPP Kinematic (FG) | RTKLIB | 2.326 m | 4.253 m | 58.704 m | RTKLIB |
| GSDC (Pixel 4) | PPP Kinematic (FG) | Gneiss | 178.824 m | 195.867 m | 30.284 m | |
| Shinjuku (UrbanNav) | SPP | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) |
| Shinjuku (UrbanNav) | SPP | Gneiss | 1.849 m | 3.167 m | 2.446 m | |
| Shinjuku (UrbanNav) | RTK Kinematic | RTKLIB | 2.205 m | 5.730 m | 4.115 m | |
| Shinjuku (UrbanNav) | RTK Kinematic | Gneiss | 1.214 m | 1.339 m | 2.037 m | **Gneiss** |
| Shinjuku (UrbanNav) | RTK Kinematic (combined) | RTKLIB | 2.809 m | 5.978 m | 4.700 m | RTKLIB |
| Shinjuku (UrbanNav) | RTK Kinematic (combined) | Gneiss | 6.036 m | 14.981 m | 8.143 m | |
| Shinjuku (UrbanNav) | PPP Kinematic (EKF) | RTKLIB | 2.003 m | 3.749 m | 3.235 m | RTKLIB |
| Shinjuku (UrbanNav) | PPP Kinematic (EKF) | Gneiss | 30.267 m | 54.412 m | 223.067 m | |
| Shinjuku (UrbanNav) | PPP Kinematic (FG) | RTKLIB | 2.003 m | 3.749 m | 3.235 m | RTKLIB |
| Shinjuku (UrbanNav) | PPP Kinematic (FG) | Gneiss | 4.818 m | 7.208 m | 13.388 m | |
| Odaiba (UrbanNav) | SPP | RTKLIB | 2.799 m | 7.453 m | 4.885 m | **Gneiss** |
| Odaiba (UrbanNav) | SPP | Gneiss | 2.099 m | 2.958 m | 2.253 m | |
| Odaiba (UrbanNav) | RTK Kinematic | RTKLIB | 2.867 m | 4.014 m | 2.490 m | RTKLIB |
| Odaiba (UrbanNav) | RTK Kinematic | Gneiss | 5.664 m | 17.037 m | 3.069 m | |
| Odaiba (UrbanNav) | RTK Kinematic (combined) | RTKLIB | 3.031 m | 3.886 m | 4.168 m | RTKLIB |
| Odaiba (UrbanNav) | RTK Kinematic (combined) | Gneiss | 5.664 m | 17.037 m | 3.069 m | |
| Odaiba (UrbanNav) | PPP Kinematic (EKF) | RTKLIB | 3.965 m | 4.691 m | 6.089 m | RTKLIB |
| Odaiba (UrbanNav) | PPP Kinematic (EKF) | Gneiss | 16.549 m | 21.594 m | 14.939 m | |
| Odaiba (UrbanNav) | PPP Kinematic (FG) | RTKLIB | 3.965 m | 4.691 m | 6.089 m | — |
| Odaiba (UrbanNav) | PPP Kinematic (FG) | Gneiss | 3.299 m | 4.916 m | 6.111 m | **Tied** |
| PPP (f9p_ppp) | SPP | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) |
| PPP (f9p_ppp) | SPP | Gneiss | 1.367 m | 2.008 m | 0.934 m | |
| PPP (f9p_ppp) | RTK Kinematic | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) |
| PPP (f9p_ppp) | RTK Kinematic | Gneiss | 0.246 m | 0.402 m | 0.168 m | |
| PPP (f9p_ppp) | RTK Kinematic (combined) | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) |
| PPP (f9p_ppp) | RTK Kinematic (combined) | Gneiss | 0.246 m | 0.402 m | 0.168 m | |
| PPP (f9p_ppp) | PPP Kinematic (EKF) | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) |
| PPP (f9p_ppp) | PPP Kinematic (EKF) | Gneiss | 13.137 m | 16.818 m | 16.127 m | |
| PPP (f9p_ppp) | PPP Kinematic (FG) | RTKLIB | Failed |  |  | **Gneiss** (RTKLIB missing) |
| PPP (f9p_ppp) | PPP Kinematic (FG) | Gneiss | 3.516 m | 5.081 m | 1.263 m | |
| UrbanLoco (Example) | SPP | RTKLIB | No data |  |  |  |
| UrbanLoco (Example) | SPP | Gneiss | No data |  |  | |
| UrbanLoco (Example) | RTK Kinematic | RTKLIB | No data |  |  |  |
| UrbanLoco (Example) | RTK Kinematic | Gneiss | No data |  |  | |
| UrbanLoco (Example) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  |
| UrbanLoco (Example) | RTK Kinematic (combined) | Gneiss | No data |  |  | |
| UrbanLoco (Example) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  |
| UrbanLoco (Example) | PPP Kinematic (EKF) | Gneiss | No data |  |  | |
| UrbanLoco (Example) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  |
| UrbanLoco (Example) | PPP Kinematic (FG) | Gneiss | No data |  |  | |
| TEX-CUP (UT Austin) | SPP | RTKLIB | No data |  |  |  |
| TEX-CUP (UT Austin) | SPP | Gneiss | No data |  |  | |
| TEX-CUP (UT Austin) | RTK Kinematic | RTKLIB | No data |  |  |  |
| TEX-CUP (UT Austin) | RTK Kinematic | Gneiss | No data |  |  | |
| TEX-CUP (UT Austin) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  |
| TEX-CUP (UT Austin) | RTK Kinematic (combined) | Gneiss | No data |  |  | |
| TEX-CUP (UT Austin) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  |
| TEX-CUP (UT Austin) | PPP Kinematic (EKF) | Gneiss | No data |  |  | |
| TEX-CUP (UT Austin) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  |
| TEX-CUP (UT Austin) | PPP Kinematic (FG) | Gneiss | No data |  |  | |
| WHU-Smartphone (Xiaomi) | SPP | RTKLIB | No data |  |  |  |
| WHU-Smartphone (Xiaomi) | SPP | Gneiss | No data |  |  | |
| WHU-Smartphone (Xiaomi) | RTK Kinematic | RTKLIB | No data |  |  |  |
| WHU-Smartphone (Xiaomi) | RTK Kinematic | Gneiss | No data |  |  | |
| WHU-Smartphone (Xiaomi) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  |
| WHU-Smartphone (Xiaomi) | RTK Kinematic (combined) | Gneiss | No data |  |  | |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (EKF) | Gneiss | No data |  |  | |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  |
| WHU-Smartphone (Xiaomi) | PPP Kinematic (FG) | Gneiss | No data |  |  | |
| smartLoc (TU Chemnitz) | SPP | RTKLIB | No data |  |  |  |
| smartLoc (TU Chemnitz) | SPP | Gneiss | No data |  |  | |
| smartLoc (TU Chemnitz) | RTK Kinematic | RTKLIB | No data |  |  |  |
| smartLoc (TU Chemnitz) | RTK Kinematic | Gneiss | No data |  |  | |
| smartLoc (TU Chemnitz) | RTK Kinematic (combined) | RTKLIB | No data |  |  |  |
| smartLoc (TU Chemnitz) | RTK Kinematic (combined) | Gneiss | No data |  |  | |
| smartLoc (TU Chemnitz) | PPP Kinematic (EKF) | RTKLIB | No data |  |  |  |
| smartLoc (TU Chemnitz) | PPP Kinematic (EKF) | Gneiss | No data |  |  | |
| smartLoc (TU Chemnitz) | PPP Kinematic (FG) | RTKLIB | No data |  |  |  |
| smartLoc (TU Chemnitz) | PPP Kinematic (FG) | Gneiss | No data |  |  | |
