#!/usr/bin/env python3
"""Evaluate all benchmark results and compute error percentiles."""
import os, math, json, re

BASE_DIR = "/Users/kevin/projects/gneiss"
BENCH_DIR = os.path.join(BASE_DIR, "benchmarks")
OUT_DIR = os.path.join(BENCH_DIR, "results")
DATASETS_DIR = os.path.join(BASE_DIR, "datasets")

# IGS ground truth (from RINEX APPROX POSITION XYZ) in ECEF
IGS_TRUTH = {
    "alic": (-4052051.7670, 4212836.2150, -2545106.0270),
    "cedu": (-3753472.1430, 3912741.0403, -3347961.0368),
    "yarr": (-2389024.5060, 5043315.3590, -3078534.1630),
    "hob2": (-3950071.2740, 2522415.2180, -4311638.5120),
    "park": (-4554255.0460, 2816652.4450, -3454059.9010),
    "pert": (-2368686.7490, 4881316.4850, -3341796.3890),
    "nklg": (6287385.7333, 1071574.6758, 39133.0395),
}

def ecef_to_llh_ref(x, y, z):
    """Convert ECEF to geodetic LLH (reference point)."""
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = 2*f - f*f
    r2 = x*x + y*y
    lon = math.atan2(y, x)
    lat = math.atan2(z, math.sqrt(r2) * (1 - e2))
    for _ in range(5):
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1 - e2 * sin_lat*sin_lat)
        lat = math.atan2(z + e2 * N * sin_lat, math.sqrt(r2))
    return lat, lon

def ecef_diff_to_enu(dx, dy, dz, ref_x, ref_y, ref_z):
    """Convert ECEF differences to ENU."""
    lat, lon = ecef_to_llh_ref(ref_x, ref_y, ref_z)
    sin_lat = math.sin(lat); cos_lat = math.cos(lat)
    sin_lon = math.sin(lon); cos_lon = math.cos(lon)
    e = -sin_lon*dx + cos_lon*dy
    n = -sin_lat*cos_lon*dx - sin_lat*sin_lon*dy + cos_lat*dz
    u = cos_lat*cos_lon*dx + cos_lat*sin_lon*dy + sin_lat*dz
    return e, n, u

def read_pos_file(path):
    """Read gneiss POS file returning list of (x, y, z, q)."""
    positions = []
    if not os.path.exists(path):
        return positions
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('%') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    x, y, z, q = float(parts[2]), float(parts[3]), float(parts[4]), int(parts[5])
                    positions.append((x, y, z, q))
                except (ValueError, IndexError):
                    continue
    return positions

def compute_percentile(sorted_vals, pct):
    if not sorted_vals:
        return None
    idx = max(0, min(int(len(sorted_vals) * pct / 100.0), len(sorted_vals) - 1))
    return sorted_vals[idx]

def evaluate_igs_station(station, pos_file):
    """Evaluate an IGS station against known truth."""
    if not os.path.exists(pos_file):
        return None

    positions = read_pos_file(pos_file)
    if not positions:
        return None

    station_lower = station.lower().replace("_ionex", "")
    tx, ty, tz = IGS_TRUTH[station_lower]

    hz_errors = []
    vt_errors = []
    d3_errors = []
    fixes = 0
    total = len(positions)
    coast_epochs = 0
    coast_segments = 0
    prev_fixed = True

    for x, y, z, q in positions:
        dx, dy, dz = x - tx, y - ty, z - tz
        e, n, u = ecef_diff_to_enu(dx, dy, dz, tx, ty, tz)
        hz = math.sqrt(e*e + n*n)
        vt = abs(u)
        d3 = math.sqrt(hz*hz + vt*vt)

        hz_errors.append(hz)
        vt_errors.append(vt)
        d3_errors.append(d3)

        if q == 1:
            fixes += 1

        # Coasting detection (Q >= 5 typically means coasting, Q=2 means float)
        if q >= 5:
            coast_epochs += 1
            if prev_fixed:
                coast_segments += 1
                prev_fixed = False
        else:
            prev_fixed = True

    hz_errors.sort()
    vt_errors.sort()
    d3_errors.sort()

    # Also extract coasting from log
    log_file = pos_file.replace('.pos', '.log')
    coast_log_segments = 0
    ar_fix_attempts = 0
    if os.path.exists(log_file):
        with open(log_file) as f:
            content = f.read()
            coast_log_segments = len(re.findall(r"coasting", content))
            ar_fix_attempts = len(re.findall(r"PPP Cascade AR Fixed", content))

    # Compute percentiles
    median_hz = compute_percentile(hz_errors, 50)
    p95_hz = compute_percentile(hz_errors, 95)
    median_vt = compute_percentile(vt_errors, 50)
    median_3d = compute_percentile(d3_errors, 50)
    p68_3d = compute_percentile(d3_errors, 68)
    p95_3d = compute_percentile(d3_errors, 95)

    # Also compute 68th and 95th for HZ
    p68_hz = compute_percentile(hz_errors, 68)
    p95_vt = compute_percentile(vt_errors, 95)

    fix_rate = (fixes / total * 100) if total > 0 else 0

    return {
        "station": station,
        "total_epochs": total,
        "median_hz": median_hz,
        "p68_hz": p68_hz,
        "p95_hz": p95_hz,
        "median_vt": median_vt,
        "p95_vt": p95_vt,
        "median_3d": median_3d,
        "p68_3d": p68_3d,
        "p95_3d": p95_3d,
        "fix_rate": fix_rate,
        "fixes": fixes,
        "coast_epochs": coast_epochs,
        "coast_segments": coast_segments if coast_segments > 0 else coast_log_segments,
        "ar_fix_attempts": ar_fix_attempts,
    }

def evaluate_urbannav(pos_file, ref_file):
    """Evaluate UrbanNav using gneiss-cli eval."""
    CLI = os.path.join(BASE_DIR, "target/release/gneiss-cli")
    cmd = f"{CLI} eval --solution {pos_file} --truth {ref_file}"
    try:
        proc = __import__('subprocess').run(cmd, capture_output=True, text=True, shell=True, timeout=300)
        output = proc.stdout + "\n" + proc.stderr
    except Exception as e:
        return {"error": str(e)}

    # Parse the eval table
    metrics = {}
    for line in output.splitlines():
        # Parse the table format: | Horiz 25% | 0.840 m | 1.234 m | ... |
        if line.startswith("| Horiz") and "Metric" not in line:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 7:
                metrics["hz_50"] = parts[2]
                metrics["hz_95"] = parts[5]
        elif line.startswith("| Vert"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 7:
                metrics["vt_50"] = parts[2]
        elif line.startswith("| 3D"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 7:
                metrics["3d_25"] = parts[1]
                metrics["3d_50"] = parts[2]
                metrics["3d_75"] = parts[3]
                metrics["3d_95"] = parts[5]

    m = re.search(r"Evaluated (\d+) epochs", output)
    if m:
        metrics["epochs"] = int(m.group(1))

    return metrics if metrics else {"error": "No metrics parsed"}

def parse_metric(m):
    """Parse a metric string like '0.840 m' to float."""
    if m is None:
        return None
    m = str(m).replace('m', '').replace('°', '').strip()
    try:
        return float(m)
    except (ValueError, TypeError):
        return None

def main():
    results = {}

    print("=" * 70)
    print("EVALUATING ALL BENCHMARKS")
    print("=" * 70)

    # =================================================================
    # IGS STATIONS
    # =================================================================
    print("\n--- IGS Stations (GPS PPP + AR + SP3/CLK) ---")

    for station in ["alic", "cedu", "yarr", "hob2", "park", "pert"]:
        pos_file = os.path.join(OUT_DIR, f"igs_{station}.pos")
        r = evaluate_igs_station(station, pos_file)
        if r:
            results[f"igs_{station}"] = r
            print(f"  {station.upper()}:")
            print(f"    3D: median={r['median_3d']:.3f}m, 68%={r['p68_3d']:.3f}m, 95%={r['p95_3d']:.3f}m")
            print(f"    HZ: median={r['median_hz']:.3f}m, 95%={r['p95_hz']:.3f}m")
            print(f"    Fix: {r['fix_rate']:.1f}% ({r['fixes']}/{r['total_epochs']})")
            print(f"    Coast: {r['coast_epochs']} epochs, {r['coast_segments']} segments")
        else:
            print(f"  {station.upper()}: [FAIL]")

    # =================================================================
    # NKLG WITH IONEX
    # =================================================================
    print("\n--- IGS NKLG (GPS PPP + AR + IONEX + SP3/CLK) ---")

    pos_file = os.path.join(OUT_DIR, "igs_nklg_ionex.pos")
    r = evaluate_igs_station("nklg_ionex", pos_file)
    if r:
        results["igs_nklg_ionex"] = r
        print(f"  NKLG (IONEX):")
        print(f"    3D: median={r['median_3d']:.3f}m, 68%={r['p68_3d']:.3f}m, 95%={r['p95_3d']:.3f}m")
        print(f"    HZ: median={r['median_hz']:.3f}m, 95%={r['p95_hz']:.3f}m")
        print(f"    Fix: {r['fix_rate']:.1f}% ({r['fixes']}/{r['total_epochs']})")
        print(f"    Coast: {r['coast_epochs']} epochs, {r['coast_segments']} segments")
    else:
        print("  NKLG (IONEX): [FAIL]")

    # =================================================================
    # ODAIBA
    # =================================================================
    print("\n--- Odaiba (GPS+Galileo PPP + SP3/CLK/BIA) ---")

    odaiba_pos = os.path.join(OUT_DIR, "urbannav_odaiba_ppp.pos")
    odaiba_ref = os.path.join(DATASETS_DIR, "urbannav", "tokyo", "Tokyo_Data", "Odaiba", "reference.csv")

    if os.path.exists(odaiba_pos) and os.path.exists(odaiba_ref):
        metrics = evaluate_urbannav(odaiba_pos, odaiba_ref)
        if metrics and "error" not in metrics:
            positions = read_pos_file(odaiba_pos)
            fixes = sum(1 for _, _, _, q in positions if q == 1)
            total = len(positions)

            # Coasting from POS
            coast_epochs = sum(1 for _, _, _, q in positions if q >= 5)
            prev_fixed = True
            coast_segments = 0
            for _, _, _, q in positions:
                if q >= 5:
                    if prev_fixed:
                        coast_segments += 1
                        prev_fixed = False
                else:
                    prev_fixed = True

            hz_50 = parse_metric(metrics.get("hz_50"))
            hz_95 = parse_metric(metrics.get("hz_95"))
            vt_50 = parse_metric(metrics.get("vt_50"))
            d3_50 = parse_metric(metrics.get("3d_50"))
            d3_68 = parse_metric(metrics.get("3d_75"))  # Approximate 68% with 75th
            d3_95 = parse_metric(metrics.get("3d_95"))

            results["urbannav_odaiba"] = {
                "median_hz": hz_50,
                "p95_hz": hz_95,
                "median_vt": vt_50,
                "median_3d": d3_50,
                "p68_3d": d3_68,
                "p95_3d": d3_95,
                "fix_rate": (fixes / total * 100) if total > 0 else 0,
                "fixes": fixes,
                "total_epochs": total,
                "coast_epochs": coast_epochs,
                "coast_segments": coast_segments,
            }

            print(f"    3D: median={d3_50:.3f}m, 68%={d3_68:.3f}m, 95%={d3_95:.3f}m")
            print(f"    HZ: median={hz_50:.3f}m, 95%={hz_95:.3f}m")
            print(f"    Fix: {results['urbannav_odaiba']['fix_rate']:.1f}% ({fixes}/{total})")
            print(f"    Coast: {coast_epochs} epochs, {coast_segments} segments")
        else:
            print(f"  [FAIL] Could not evaluate Odaiba: {metrics}")
    else:
        print("  [SKIP] Missing Odaiba files")

    # =================================================================
    # SHINJUKU
    # =================================================================
    print("\n--- Shinjuku (GPS+Galileo PPP + SP3/CLK/BIA) ---")

    shinjuku_pos = os.path.join(OUT_DIR, "urbannav_shinjuku_ppp.pos")
    shinjuku_ref = os.path.join(DATASETS_DIR, "urbannav", "tokyo", "Tokyo_Data", "Shinjuku", "reference.csv")

    if os.path.exists(shinjuku_pos) and os.path.exists(shinjuku_ref):
        metrics = evaluate_urbannav(shinjuku_pos, shinjuku_ref)
        if metrics and "error" not in metrics:
            positions = read_pos_file(shinjuku_pos)
            fixes = sum(1 for _, _, _, q in positions if q == 1)
            total = len(positions)

            coast_epochs = sum(1 for _, _, _, q in positions if q >= 5)
            prev_fixed = True
            coast_segments = 0
            for _, _, _, q in positions:
                if q >= 5:
                    if prev_fixed:
                        coast_segments += 1
                        prev_fixed = False
                else:
                    prev_fixed = True

            hz_50 = parse_metric(metrics.get("hz_50"))
            hz_95 = parse_metric(metrics.get("hz_95"))
            vt_50 = parse_metric(metrics.get("vt_50"))
            d3_50 = parse_metric(metrics.get("3d_50"))
            d3_68 = parse_metric(metrics.get("3d_75"))  # Approximate 68% with 75th
            d3_95 = parse_metric(metrics.get("3d_95"))

            results["urbannav_shinjuku"] = {
                "median_hz": hz_50,
                "p95_hz": hz_95,
                "median_vt": vt_50,
                "median_3d": d3_50,
                "p68_3d": d3_68,
                "p95_3d": d3_95,
                "fix_rate": (fixes / total * 100) if total > 0 else 0,
                "fixes": fixes,
                "total_epochs": total,
                "coast_epochs": coast_epochs,
                "coast_segments": coast_segments,
            }

            print(f"    3D: median={d3_50:.3f}m, 68%={d3_68:.3f}m, 95%={d3_95:.3f}m")
            print(f"    HZ: median={hz_50:.3f}m, 95%={hz_95:.3f}m")
            print(f"    Fix: {results['urbannav_shinjuku']['fix_rate']:.1f}% ({fixes}/{total})")
            print(f"    Coast: {coast_epochs} epochs, {coast_segments} segments")
        else:
            print(f"  [FAIL] Could not evaluate Shinjuku: {metrics}")
    else:
        print("  [SKIP] Missing Shinjuku files")

    # =================================================================
    # SUMMARY TABLE
    # =================================================================
    print("\n" + "=" * 110)
    print("COMPREHENSIVE BENCHMARK MATRIX - SUMMARY")
    print("=" * 110)
    print(f"\n{'Dataset':<18} {'Mode':<32} {'3D 50%':<12} {'3D 68%':<12} {'3D 95%':<12} {'Coasts':<10} {'Fix%':<8}")
    print("-" * 104)

    for key in sorted(results.keys()):
        r = results[key]
        if key == "igs_nklg_ionex":
            ds = "IGS NKLG"
            mode = "GPS PPP+AR+IONEX+SP3/CLK"
        elif key.startswith("igs_"):
            st = key[4:].upper()
            ds = f"IGS {st}"
            mode = "GPS PPP+AR+SP3/CLK"
        elif key == "urbannav_odaiba":
            ds = "Odaiba"
            mode = "GPS+Galileo PPP+SP3/CLK/BIA"
        elif key == "urbannav_shinjuku":
            ds = "Shinjuku"
            mode = "GPS+Galileo PPP+SP3/CLK/BIA"
        else:
            ds = key
            mode = ""

        d3_50s = f"{r['median_3d']:.3f}m" if r.get('median_3d') is not None else "N/A"
        d3_68s = f"{r['p68_3d']:.3f}m" if r.get('p68_3d') is not None else "N/A"
        d3_95s = f"{r['p95_3d']:.3f}m" if r.get('p95_3d') is not None else "N/A"
        coasts = str(r.get('coast_epochs', 'N/A'))
        fix = f"{r['fix_rate']:.1f}%" if r.get('fix_rate', 0) > 0 else "0%"

        print(f"{ds:<18} {mode:<32} {d3_50s:<12} {d3_68s:<12} {d3_95s:<12} {coasts:<10} {fix:<8}")

    # Save to JSON
    serializable = {}
    for key, r in results.items():
        serializable[key] = {
            k: (round(v, 3) if isinstance(v, float) else v)
            for k, v in r.items()
        }
    json_path = os.path.join(OUT_DIR, "benchmark_results.json")
    with open(json_path, 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f"\nResults saved to {json_path}")

    return results

if __name__ == "__main__":
    results = main()
