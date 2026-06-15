import sys, math

def read_pos(file_path):
    pos = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('%') or not line.strip(): continue
            parts = line.split()
            if len(parts) >= 5:
                pos.append((float(parts[2]), float(parts[3]), float(parts[4])))
    return pos

ref_pos = [35.67093043, 139.77189916, 22.0] # Tokyo? No, f9p_ppp_1224!
# Wait, let me get the ground truth from rover_csrs.pos!
gt_pos = []
with open('datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_csrs.pos', 'r') as f:
    for line in f:
        if line.startswith('%') or not line.strip(): continue
        parts = line.split()
        if len(parts) >= 5:
            gt_pos.append((float(parts[2]), float(parts[3]), float(parts[4])))

# Compute average of ground truth (it's static)
gt_lat = sum(p[0] for p in gt_pos) / len(gt_pos)
gt_lon = sum(p[1] for p in gt_pos) / len(gt_pos)
gt_alt = sum(p[2] for p in gt_pos) / len(gt_pos)

print(f"GT: {gt_lat}, {gt_lon}, {gt_alt}")

pos = read_pos('benchmarks/ppp17/rover.pos')

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return R * c

errors = []
for p in pos[-100:]:
    h_err = haversine(p[0], p[1], gt_lat, gt_lon)
    v_err = abs(p[2] - gt_alt)
    d3_err = math.sqrt(h_err**2 + v_err**2)
    errors.append(d3_err)

if errors:
    print(f"Last 100 epochs avg 3D error: {sum(errors)/len(errors):.3f} m")
    print(f"Min: {min(errors):.3f} m, Max: {max(errors):.3f} m")
else:
    print("No data")
