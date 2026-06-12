import math
import sys

def parse_pos(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('%') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 15:
                # time = parts[0] + ' ' + parts[1]
                lat = float(parts[2])
                lon = float(parts[3])
                height = float(parts[4])
                data.append((lat, lon, height))
    return data

def distance(lat1, lon1, lat2, lon2):
    R = 6371000
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    a = math.sin(delta_phi/2) * math.sin(delta_phi/2) + \
        math.cos(phi1) * math.cos(phi2) * \
        math.sin(delta_lambda/2) * math.sin(delta_lambda/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return R * c

spp = parse_pos('benchmarks/GSDC_Pixel_4_spp.pos')
ppp = parse_pos('benchmarks/GSDC_Pixel_4_ppp.pos')

min_len = min(len(spp), len(ppp))
if min_len == 0:
    print("No data")
    sys.exit(0)

diffs = []
for i in range(min_len):
    d = distance(spp[i][0], spp[i][1], ppp[i][0], ppp[i][1])
    diffs.append(d)

print(f"Epochs: {min_len}")
print(f"Max distance SPP vs PPP: {max(diffs):.2f} m")
print(f"Avg distance SPP vs PPP: {sum(diffs)/min_len:.2f} m")
for i in range(0, min_len, max(1, min_len//20)):
    print(f"Epoch {i:4d}: {diffs[i]:.2f} m")

