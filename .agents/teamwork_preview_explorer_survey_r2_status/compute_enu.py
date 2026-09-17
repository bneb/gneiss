import math
import numpy as np
import datetime

def parse_time_to_tow(date_str, time_str):
    y, m, d = [int(x) for x in date_str.split('/')]
    time_parts = time_str.split(':')
    hr = int(time_parts[0])
    mn = int(time_parts[1])
    sec_float = float(time_parts[2])
    sec_int = int(sec_float)
    usec = int(round((sec_float - sec_int) * 1e6))
    if usec >= 1000000:
        sec_int += 1
        usec = 0
    dt = datetime.datetime(y, m, d, hr, mn, sec_int, usec)
    gps_epoch = datetime.datetime(1980, 1, 6)
    diff_sec = (dt - gps_epoch).total_seconds()
    tow = diff_sec % 604800.0
    return int(round(tow))

def parse_csrs(path):
    epochs = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                tow_round = parse_time_to_tow(parts[0], parts[1])
                lat_deg = float(parts[2])
                lon_deg = float(parts[3])
                h_m = float(parts[4])
                q = int(parts[5])
                
                lat_rad = math.radians(lat_deg)
                lon_rad = math.radians(lon_deg)
                a = 6378137.0
                f_inv = 298.257223563
                f_val = 1.0 / f_inv
                e2 = f_val * (2.0 - f_val)
                n_val = a / math.sqrt(1.0 - e2 * math.sin(lat_rad)**2)
                
                x = (n_val + h_m) * math.cos(lat_rad) * math.cos(lon_rad)
                y = (n_val + h_m) * math.cos(lat_rad) * math.sin(lon_rad)
                z = (n_val * (1.0 - e2) + h_m) * math.sin(lat_rad)
                
                epochs[tow_round] = {
                    'pos': np.array([x, y, z]),
                    'llh': (lat_rad, lon_rad, h_m),
                    'q': q
                }
    return epochs

def parse_ppk(path):
    epochs = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                tow_round = parse_time_to_tow(parts[0], parts[1])
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                q = int(parts[5])
                
                epochs[tow_round] = {
                    'pos': np.array([x, y, z]),
                    'q': q
                }
    return epochs

def ecef_to_enu(diff, lat, lon):
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    
    e = -sin_lon * diff[0] + cos_lon * diff[1]
    n = -sin_lat * cos_lon * diff[0] - sin_lat * sin_lon * diff[1] + cos_lat * diff[2]
    u = cos_lat * cos_lon * diff[0] + cos_lat * sin_lon * diff[1] + sin_lat * diff[2]
    return e, n, u

def print_block(title, e, n, u, h, d3):
    e = np.array(e)
    n = np.array(n)
    u = np.array(u)
    h = np.array(h)
    d3 = np.array(d3)
    print(f"\n=== {title} (N={len(e)}) ===")
    print(f"East (m):  mean={np.mean(e):+7.4f}, std={np.std(e):6.4f}, RMS={np.sqrt(np.mean(e**2)):6.4f}, min={np.min(e):+7.4f}, max={np.max(e):+7.4f}")
    print(f"North (m): mean={np.mean(n):+7.4f}, std={np.std(n):6.4f}, RMS={np.sqrt(np.mean(n**2)):6.4f}, min={np.min(n):+7.4f}, max={np.max(n):+7.4f}")
    print(f"Up (m):    mean={np.mean(u):+7.4f}, std={np.std(u):6.4f}, RMS={np.sqrt(np.mean(u**2)):6.4f}, min={np.min(u):+7.4f}, max={np.max(u):+7.4f}")
    print(f"Horizontal (m): p50={np.median(h):6.4f}, p68={np.percentile(h, 68):6.4f}, p95={np.percentile(h, 95):6.4f}, RMS={np.sqrt(np.mean(h**2)):6.4f}, max={np.max(h):6.4f}")
    print(f"3D Error (m):   p50={np.median(d3):6.4f}, p68={np.percentile(d3, 68):6.4f}, p95={np.percentile(d3, 95):6.4f}, RMS={np.sqrt(np.mean(d3**2)):6.4f}, max={np.max(d3):6.4f}")

csrs = parse_csrs('datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_csrs.pos')
ppk = parse_ppk('datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppk.pos')

print(f"Loaded CSRS: {len(csrs)} epochs, PPK: {len(ppk)} epochs")

# 1. CSRS vs PPK (Commercial Baseline: CSRS-PPP vs RTK Truth)
common_tows = sorted(set(csrs.keys()) & set(ppk.keys()))
print(f"Common CSRS and PPK epochs: {len(common_tows)}")

e_list, n_list, u_list, h_list, d3_list = [], [], [], [], []
for tow in common_tows:
    c_pos = csrs[tow]['pos']
    p_pos = ppk[tow]['pos']
    lat, lon, _ = csrs[tow]['llh']
    
    diff = c_pos - p_pos
    e, n, u = ecef_to_enu(diff, lat, lon)
    h = math.sqrt(e**2 + n**2)
    d3 = math.sqrt(e**2 + n**2 + u**2)
    
    e_list.append(e)
    n_list.append(n)
    u_list.append(u)
    h_list.append(h)
    d3_list.append(d3)

first_600_tows = sorted(csrs.keys())[:600]
matched_ppk_tows = [t for t in first_600_tows if t in ppk]
print_block("Commercial Baseline: CSRS-PPP vs RTK Truth (Matched N=583)", 
            [e_list[i] for i, tow in enumerate(common_tows) if tow in matched_ppk_tows],
            [n_list[i] for i, tow in enumerate(common_tows) if tow in matched_ppk_tows],
            [u_list[i] for i, tow in enumerate(common_tows) if tow in matched_ppk_tows],
            [h_list[i] for i, tow in enumerate(common_tows) if tow in matched_ppk_tows],
            [d3_list[i] for i, tow in enumerate(common_tows) if tow in matched_ppk_tows])

# 2. Gneiss Static PPP Solution vs CSRS-PPP (First 600 epochs)
gneiss_sol = np.array([-1276969.350, -4717194.952, 4087249.058])
first_600_tows = sorted(csrs.keys())[:600]

e_g_csrs, n_g_csrs, u_g_csrs, h_g_csrs, d3_g_csrs = [], [], [], [], []
for tow in first_600_tows:
    c_pos = csrs[tow]['pos']
    lat, lon, _ = csrs[tow]['llh']
    diff = gneiss_sol - c_pos
    e, n, u = ecef_to_enu(diff, lat, lon)
    h = math.sqrt(e**2 + n**2)
    d3 = math.sqrt(e**2 + n**2 + u**2)
    e_g_csrs.append(e)
    n_g_csrs.append(n)
    u_g_csrs.append(u)
    h_g_csrs.append(h)
    d3_g_csrs.append(d3)

print_block("Discrepancy: Current Gneiss PPP vs CSRS-PPP (N=600)", e_g_csrs, n_g_csrs, u_g_csrs, h_g_csrs, d3_g_csrs)

# 3. Gneiss Static PPP Solution vs RTK Truth (Matched N=583)
e_g_ppk, n_g_ppk, u_g_ppk, h_g_ppk, d3_g_ppk = [], [], [], [], []
matched_ppk_tows = [t for t in first_600_tows if t in ppk]
for tow in matched_ppk_tows:
    p_pos = ppk[tow]['pos']
    lat, lon, _ = csrs[tow]['llh']
    diff = gneiss_sol - p_pos
    e, n, u = ecef_to_enu(diff, lat, lon)
    h = math.sqrt(e**2 + n**2)
    d3 = math.sqrt(e**2 + n**2 + u**2)
    e_g_ppk.append(e)
    n_g_ppk.append(n)
    u_g_ppk.append(u)
    h_g_ppk.append(h)
    d3_g_ppk.append(d3)

print_block("Current Gneiss PPP vs RTK Truth (Matched N=583)", e_g_ppk, n_g_ppk, u_g_ppk, h_g_ppk, d3_g_ppk)
