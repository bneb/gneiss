import numpy as np

def read_ecef(path):
    times, positions = [], []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('%') or not line.strip(): continue
            parts = line.split()
            try:
                times.append(float(parts[1]))
                positions.append([float(parts[2]), float(parts[3]), float(parts[4])])
            except: pass
    return np.array(times), np.array(positions)

t1, p1 = read_ecef('scratch/test.pos')
t2, p2 = read_ecef('datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppp_rapid.pos')

# Align by closest time
import scipy.interpolate
interp = scipy.interpolate.interp1d(t2, p2, axis=0, bounds_error=False, fill_value=np.nan)
p2_interp = interp(t1)
valid = ~np.isnan(p2_interp[:,0])
t1 = t1[valid]
p1 = p1[valid]
p2_interp = p2_interp[valid]

errs = np.linalg.norm(p1 - p2_interp, axis=1)
print(f"Max 3D Error: {np.max(errs):.2f} m")

