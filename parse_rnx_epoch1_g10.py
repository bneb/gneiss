c = 299792458.0
f1 = 1575.42e6
f2 = 1227.60e6
lam1 = c / f1
lam2 = c / f2
lam_wl = c / (f1 - f2)

with open("datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.rnx") as f:
    lines = f.readlines()

in_obs = False
for i, line in enumerate(lines):
    if "END OF HEADER" in line:
        in_obs = True
        continue
    if in_obs and line.startswith(">"):
        num_sats = int(line[32:35])
        for j in range(num_sats):
            sat_line = lines[i+1+j]
            if sat_line.startswith("G10"):
                parts = []
                for k in range(3, len(sat_line), 16):
                    parts.append(sat_line[k:k+14].strip())
                print(parts)
                C1C = float(parts[0])
                L1C = float(parts[1])
                C2W = float(parts[8])
                L2W = float(parts[9])
                
                mw_raw = (L1C - L2W) - (f1 * C1C + f2 * C2W) / (f1 + f2) / lam_wl
                print("Raw MW G10:", mw_raw)
                
                # apply bias
                b_c1 = -0.3735 * 1e-9 * c
                b_c2 = 0.0 * 1e-9 * c
                b_l1 = -0.36346 * 1e-9 * c
                b_l2 = -0.49982 * 1e-9 * c
                
                L1C_corr = L1C - b_l1 / lam1
                L2W_corr = L2W - b_l2 / lam2
                C1C_corr = C1C - b_c1
                C2W_corr = C2W - b_c2
                
                mw_corr = (L1C_corr - L2W_corr) - (f1 * C1C_corr + f2 * C2W_corr) / (f1 + f2) / lam_wl
                print("Corr MW G10:", mw_corr)
        break
