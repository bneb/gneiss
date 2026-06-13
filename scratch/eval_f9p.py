import sys
import numpy as np

def read_pos(path):
    times = []
    positions = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('%') or not line.strip(): continue
            parts = line.split()
            if len(parts) < 15: continue
            try:
                tow = float(parts[1])
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                times.append(tow)
                positions.append([x, y, z])
            except: pass
    return np.array(times), np.array(positions)

t1, p1 = read_pos(sys.argv[1])
t2, p2 = read_pos(sys.argv[2])

# Interpolate p2 to t1
import scipy.interpolate
interp = scipy.interpolate.interp1d(t2, p2, axis=0, bounds_error=False, fill_value=np.nan)
p2_interp = interp(t1)

valid = ~np.isnan(p2_interp[:,0])
t1 = t1[valid]
p1 = p1[valid]
p2_interp = p2_interp[valid]

errs = np.linalg.norm(p1 - p2_interp, axis=1)
print(f"Max 3D Error: {np.max(errs):.2f} m")
print(f"95th %ile 3D: {np.percentile(errs, 95):.2f} m")
print(f"Median 3D Error: {np.median(errs):.2f} m")
print(f"Mean 3D Error: {np.mean(errs):.2f} m")

# check the window 427000 to 427100 specifically
idx = (t1 > 427000) & (t1 < 427100)
if np.sum(idx) > 0:
    print(f"Window Max Error: {np.max(errs[idx]):.2f} m")

