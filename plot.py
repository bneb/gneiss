import sys
import numpy as np

ref = np.array([-1276972.3782, -4717196.7865, 4087245.4526])

pos = []
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('%'): continue
        parts = line.split()
        if len(parts) >= 5:
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            pos.append([x, y, z])

pos = np.array(pos)
err = pos - ref
err_norm = np.linalg.norm(err, axis=1)

print("Final err:", err[-1])
print("Final err norm:", err_norm[-1])
print("Max err norm:", np.max(err_norm))
print("Mean err:", np.mean(err, axis=0))
