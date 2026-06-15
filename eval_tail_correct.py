import sys, math

def ecef_to_llh(x, y, z):
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2 - f)
    p = math.sqrt(x*x + y*y)
    r = math.sqrt(x*x + y*y + z*z)
    
    if p < 1e-6:
        return 90.0 if z > 0 else -90.0, 0.0, abs(z) - a * math.sqrt(1.0 - e2)
        
    lon = math.atan2(y, x)
    lat = math.atan2(z, p * (1.0 - e2))
    
    for _ in range(5):
        N = a / math.sqrt(1.0 - e2 * math.sin(lat)**2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1.0 - e2 * N / (N + h)))
        
    N = a / math.sqrt(1.0 - e2 * math.sin(lat)**2)
    alt = p / math.cos(lat) - N
    return math.degrees(lat), math.degrees(lon), alt

pos = []
with open('benchmarks/ppp17/rover.pos', 'r') as f:
    for line in f:
        if line.startswith('%') or not line.strip(): continue
        parts = line.split()
        if len(parts) >= 5:
            pos.append((float(parts[2]), float(parts[3]), float(parts[4])))

gt_lat, gt_lon, gt_alt = 40.0945373546018, -105.15534826974691, 1570.7479668448007

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
for p in pos:
    lat, lon, alt = ecef_to_llh(p[0], p[1], p[2])
    h_err = haversine(lat, lon, gt_lat, gt_lon)
    v_err = abs(alt - gt_alt)
    d3_err = math.sqrt(h_err**2 + v_err**2)
    errors.append((h_err, v_err, d3_err))

if errors:
    print("ALL EPOCHS:")
    print(f"Mean 3D: {sum(e[2] for e in errors)/len(errors):.3f} m")
    
    last = errors[-100:]
    print("\nLAST 100 EPOCHS:")
    print(f"Mean H: {sum(e[0] for e in last)/len(last):.3f} m")
    print(f"Mean V: {sum(e[1] for e in last)/len(last):.3f} m")
    print(f"Mean 3D: {sum(e[2] for e in last)/len(last):.3f} m")
