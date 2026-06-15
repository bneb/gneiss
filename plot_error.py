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

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return R * c

gt_lat = 40.097029227777774
gt_lon = -105.14726476944445
gt_alt = 1577.3571

with open('benchmarks/ppp17/rover.pos', 'r') as f:
    lines = [l for l in f if not l.startswith('%') and l.strip()]

errors = []
for l in lines:
    parts = l.split()
    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
    lat, lon, alt = ecef_to_llh(x, y, z)
    h_err = haversine(lat, lon, gt_lat, gt_lon)
    v_err = alt - gt_alt
    d3_err = math.sqrt(h_err**2 + v_err**2)
    errors.append((h_err, v_err, d3_err, x, y, z))

print(f"Epoch 0: {errors[0][0]:.3f}m, {errors[0][1]:.3f}m, {errors[0][2]:.3f}m")
print(f"Epoch 500: {errors[500][0]:.3f}m, {errors[500][1]:.3f}m, {errors[500][2]:.3f}m")
print(f"Epoch 1000: {errors[1000][0]:.3f}m, {errors[1000][1]:.3f}m, {errors[1000][2]:.3f}m")
print(f"Epoch 1500: {errors[1500][0]:.3f}m, {errors[1500][1]:.3f}m, {errors[1500][2]:.3f}m")
print(f"Epoch 2000: {errors[2000][0]:.3f}m, {errors[2000][1]:.3f}m, {errors[2000][2]:.3f}m")
print(f"Epoch 2500: {errors[2500][0]:.3f}m, {errors[2500][1]:.3f}m, {errors[2500][2]:.3f}m")
print(f"Epoch 3000: {errors[3000][0]:.3f}m, {errors[3000][1]:.3f}m, {errors[3000][2]:.3f}m")
print(f"Epoch -1: {errors[-1][0]:.3f}m, {errors[-1][1]:.3f}m, {errors[-1][2]:.3f}m")

# Wait, why was it 774m in the previous script?
# Let's print out the raw coordinates at the end
print(f"Last X,Y,Z: {errors[-1][3]:.3f}, {errors[-1][4]:.3f}, {errors[-1][5]:.3f}")
