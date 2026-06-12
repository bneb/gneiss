# Benchmark Results

| Dataset | Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |
|:---|:---|:---|:---|:---|:---|:---|
| PPP | spp | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| PPP | spp | Gneiss | 3.057 m | 3.932 m | 3.859 m | |
| PPP | rtk_forward | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| PPP | rtk_forward | Gneiss | 2.986 m | 3.940 m | 4.017 m | |
| PPP | rtk_smoothed | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| PPP | rtk_smoothed | Gneiss | 2.986 m | 3.940 m | 4.017 m | |
| PPP | ppp_forward | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| PPP | ppp_forward | Gneiss | 3.490 m | 3.775 m | 3.115 m | |
| PPP | ppp_smoothed | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| PPP | ppp_smoothed | Gneiss | 3.045 m | 3.927 m | 3.869 m | |
| GSDC | spp | RTKLIB | 3.311 m | 10.191 m | 66.316 m | **Gneiss** ✅ |
| GSDC | spp | Gneiss | 3.297 m | 8.784 m | 62.750 m | |
| GSDC | rtk_forward | RTKLIB | 1.773 m | 4.161 m | 64.598 m | RTKLIB ⚠️ |
| GSDC | rtk_forward | Gneiss | 1.944 m | 5.569 m | 63.058 m | |
| GSDC | rtk_smoothed | RTKLIB | 1.831 m | 3.126 m | 64.471 m | RTKLIB ⚠️ |
| GSDC | rtk_smoothed | Gneiss | 3.523 m | 190.995 m | 62.961 m | |
| Shinjuku | spp | RTKLIB | 13.209 m | 55.582 m | 38.538 m | **Gneiss** ✅ |
| Shinjuku | spp | Gneiss | 2.288 m | 26.576 m | 5.166 m | |
| Shinjuku | rtk_forward | RTKLIB | 1.670 m | 11.843 m | 5.394 m | **Gneiss** ✅ |
| Shinjuku | rtk_forward | Gneiss | 1.635 m | 18.193 m | 1.988 m | |
| Shinjuku | rtk_smoothed | RTKLIB | 2.089 m | 20.154 m | 5.517 m | **Gneiss** ✅ |
| Shinjuku | rtk_smoothed | Gneiss | 1.608 m | 17.782 m | 1.992 m | |
| Shinjuku | rtk_ins_forward | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| Shinjuku | rtk_ins_forward | Gneiss | 2.716 m | 16.819 m | 2.378 m | |
| Shinjuku | rtk_ins_smoothed | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| Shinjuku | rtk_ins_smoothed | Gneiss | 5.640 m | 20.917 m | 2.421 m | |
| Odaiba | spp | RTKLIB | 5.335 m | 35.568 m | 8.330 m | **Gneiss** ✅ |
| Odaiba | spp | Gneiss | 2.409 m | 9.384 m | 2.853 m | |
| Odaiba | rtk_forward | RTKLIB | 2.239 m | 8.099 m | 6.413 m | **Gneiss** ✅ |
| Odaiba | rtk_forward | Gneiss | 1.160 m | 5.354 m | 1.516 m | |
| Odaiba | rtk_smoothed | RTKLIB | 2.236 m | 13.146 m | 8.219 m | **Gneiss** ✅ |
| Odaiba | rtk_smoothed | Gneiss | 1.129 m | 5.303 m | 1.476 m | |
| Odaiba | rtk_ins_forward | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| Odaiba | rtk_ins_forward | Gneiss | 2.312 m | 7.742 m | 2.284 m | |
| Odaiba | rtk_ins_smoothed | RTKLIB | N/A | N/A | N/A | **Gneiss** (RTKLIB missing) |
| Odaiba | rtk_ins_smoothed | Gneiss | 5.929 m | 22.493 m | 2.528 m | |
