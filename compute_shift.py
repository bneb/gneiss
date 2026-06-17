import sys

truth_file = "datasets/gsdc/reference.csv"
sol_file = "benchmarks/matrix/gneiss_GSDC_Pixel_4_spp_forward.pos"

truth_epochs = []
with open(truth_file) as f:
    f.readline()
    for line in f:
        if not line.strip(): continue
        parts = line.strip().split(',')
        if len(parts) >= 8:
            tow = float(parts[0])
            x = float(parts[5])
            y = float(parts[6])
            z = float(parts[7])
            truth_epochs.append((tow, x, y, z))

sol_epochs = []
with open(sol_file) as f:
    for line in f:
        if line.startswith('%') or not line.strip(): continue
        parts = line.strip().split()
        if len(parts) >= 5:
            tow = float(parts[1])
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            sol_epochs.append((tow, x, y, z))

dxs, dys, dzs = [], [], []

for stow, sx, sy, sz in sol_epochs:
    best_diff = 1e9
    best_true = None
    for ttow, tx, ty, tz in truth_epochs:
        diff = abs(stow - ttow)
        if diff < 0.15 and diff < best_diff:
            best_diff = diff
            best_true = (tx, ty, tz)
    
    if best_true:
        dxs.append(sx - best_true[0])
        dys.append(sy - best_true[1])
        dzs.append(sz - best_true[2])

print(f"SPP Mean dx: {sum(dxs)/len(dxs)}")
print(f"SPP Mean dy: {sum(dys)/len(dys)}")
print(f"SPP Mean dz: {sum(dzs)/len(dzs)}")

