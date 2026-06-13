import numpy as np

def read_llh(path):
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

t, p = read_llh('scratch/test.pos')
# p is in lat/lon/height. We'll just look for jumps.
dp = np.linalg.norm(p[1:] - p[:-1], axis=1)
print(f"Max epoch-to-epoch jump in pos: {np.max(dp)}")
print(f"Index of max jump: {np.argmax(dp)}")
print(f"TOW of max jump: {t[np.argmax(dp)]}")

