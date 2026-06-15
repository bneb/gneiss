import sys, math

def ecef_to_llh(x, y, z):
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2 - f)
    p = math.sqrt(x*x + y*y)
    r = math.sqrt(x*x + y*y + z*z)
    if p < 1e-6: return 90.0 if z > 0 else -90.0, 0.0, abs(z) - a * math.sqrt(1.0 - e2)
    lon = math.atan2(y, x)
    lat = math.atan2(z, p * (1.0 - e2))
    for _ in range(5):
        N = a / math.sqrt(1.0 - e2 * math.sin(lat)**2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1.0 - e2 * N / (N + h)))
    N = a / math.sqrt(1.0 - e2 * math.sin(lat)**2)
    alt = p / math.cos(lat) - N
    return math.degrees(lat), math.degrees(lon), alt

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return R * c

gt_lat, gt_lon, gt_alt = 40.097029227777774, -105.14726476944445, 1577.3571

with open('benchmarks/ppp17/rover.pos', 'r') as f:
    lines = [l for l in f if not l.startswith('%') and l.strip()]

for i in range(2563, 2575):
    l = lines[i]
    parts = l.split()
    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
    lat, lon, alt = ecef_to_llh(x, y, z)
    h_err = haversine(lat, lon, gt_lat, gt_lon)
    v_err = alt - gt_alt
    print(f"Epoch {i} ({parts[0]} {parts[1]}): H_ERR={h_err:.3f}m, V_ERR={v_err:.3f}m, X={x:.3f}, Y={y:.3f}, Z={z:.3f}")
