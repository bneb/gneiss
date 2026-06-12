import numpy as np
with open("/Users/kevin/projects/gneiss/datasets/gsdc/Pixel4_GnssLog.20o") as f:
    last_cp = {}
    last_time = {}
    epoch_diffs = {}
    for line in f:
        if line.startswith(">"):
            parts = line.split()
            try:
                h = int(parts[4])
                m = int(parts[5])
                s = float(parts[6])
                t = h * 3600 + m * 60 + s
                epoch_diffs[t] = []
            except: pass
            continue
        if len(line) > 3 and line[0] == "G":
            sat = line[0:3]
            try:
                l1c = line[19:33].strip()
                d1c = line[35:49].strip()
                if l1c and d1c:
                    cp = float(l1c)
                    dop = float(d1c)
                    if sat in last_cp:
                        prev_cp = last_cp[sat]
                        prev_t = last_time[sat]
                        dt = t - prev_t
                        if 0 < dt < 10:
                            pred_cp = prev_cp - dop * dt
                            diff = cp - pred_cp
                            epoch_diffs[t].append((sat, diff, cp, dop))
                    last_cp[sat] = cp
                    last_time[sat] = t
            except: pass

for t, diffs in epoch_diffs.items():
    if not diffs: continue
    med_diff = np.median([d[1] for d in diffs])
    for d in diffs:
        rem = abs(d[1] - med_diff)
        if rem > 2.0 and d[2] != 0.0:
            print(f"SLIP: t={t} sat={d[0]} rem={rem:.2f} diff={d[1]:.2f} med={med_diff:.2f} cp={d[2]:.2f} dop={d[3]:.2f}")
