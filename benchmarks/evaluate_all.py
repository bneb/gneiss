#!/usr/bin/env python3
"""Evaluate all benchmark results and produce summary."""
import os
import math
import json
import subprocess
import re

CLI_BIN = "/Users/kevin/projects/gneiss/target/release/gneiss-cli"
OUT_DIR = "/Users/kevin/projects/gneiss/benchmarks/results"
DATASETS_DIR = "/Users/kevin/projects/gneiss/datasets"

# IGS ground truth (from RINEX APPROX POSITION XYZ)
IGS_TRUTH = {
    "alic": (-4052051.7670, 4212836.2150, -2545106.0270),
    "cedu": (-3753472.1430, 3912741.0403, -3347961.0368),
    "yarr": (-2389024.5060, 5043315.3590, -3078534.1630),
    "hob2": (-3950071.2740, 2522415.2180, -4311638.5120),
    "park": (-4554255.0460, 2816652.4450, -3454059.9010),
    "pert": (-2368686.7490, 4881316.4850, -3341796.3890),
    "nklg": (6287385.7333, 1071574.6758, 39133.0395),
}

def ecef_to_enu(dx, dy, dz, ref_x, ref_y, ref_z):
    """Convert ECEF difference to ENU using WGS84 ellipsoid."""
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = 2*f - f*f

    r2 = ref_x*ref_x + ref_y*ref_y
    lon = math.atan2(ref_y, ref_x)
    lat = math.atan2(ref_z, r2 * (1 - e2))
    for _ in range(5):
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1 - e2 * sin_lat*sin_lat)
        lat = math.atan2(ref_z + e2 * N * sin_lat, math.sqrt(r2))

    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)

    e = -sin_lon*dx + cos_lon*dy
    n = -sin_lat*cos_lon*dx - sin_lat*sin_lon*dy + cos_lat*dz
    u = cos_lat*cos_lon*dx + cos_lat*sin_lon*dy + sin_lat*dz
    return e, n, u

def read_pos_file(path):
    """Read gneiss POS file returning list of (x, y, z, q) tuples."""
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

def evaluate_igs(station):
    """Evaluate IGS station PPP-AR result."""
    pos_file = os.path.join(OUT_DIR, f"igs_{station}.pos")
    if not os.path.exists(pos_file):
        return None

    positions = read_pos_file(pos_file)
    if not positions:
        return None

    tx, ty, tz = IGS_TRUTH[station]

    hz_errors = []
    vt_errors = []
    d3_errors = []
    fixes = 0
    total_epochs = len(positions)

    for x, y, z, q in positions:
        dx, dy, dz = x - tx, y - ty, z - tz
        e, n, u = ecef_to_enu(dx, dy, dz, tx, ty, tz)
        hz = math.sqrt(e*e + n*n)
        vt = abs(u)
        d3 = math.sqrt(hz*hz + vt*vt)

        hz_errors.append(hz)
        vt_errors.append(vt)
        d3_errors.append(d3)
        if q == 1:
            fixes += 1

    hz_errors.sort()
    vt_errors.sort()
    d3_errors.sort()

    # Count coasting events from pos file Q flag
    coast_count = sum(1 for _, _, _, q in positions if q >= 5)

    # Count AR fixes from log
    log_file = os.path.join(OUT_DIR, f"igs_{station}.log")
    ar_fix_logs = 0
    if os.path.exists(log_file):
        with open(log_file) as f:
            content = f.read()
            ar_fix_logs = len(re.findall(r"PPP Cascade AR Fixed", content))
            coast_logs = len(re.findall(r"coasting", content))

    return {
        "median_3d": compute_percentile(d3_errors, 50),
        "p68_3d": compute_percentile(d3_errors, 68),
        "p95_3d": compute_percentile(d3_errors, 95),
        "median_hz": compute_percentile(hz_errors, 50),
        "p95_hz": compute_percentile(hz_errors, 95),
        "median_vt": compute_percentile(vt_errors, 50),
        "fix_rate": fixes / total_epochs * 100 if total_epochs > 0 else 0,
        "fixes": fixes,
        "total_epochs": total_epochs,
        "coast_epochs": coast_count,
        "ar_fix_attempts": ar_fix_logs,
        "coast_segments": coast_logs,
    }

def parse_eval_output(text):
    """Parse gneiss-cli eval output."""
    metrics = {}
    for line in text.splitlines():
        if line.startswith("| Horiz"):
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
                metrics["3d_50"] = parts[2]
                metrics["3d_75"] = parts[3]
                metrics["3d_95"] = parts[5]

    m = re.search(r"Evaluated (\d+) epochs", text)
    if m:
        metrics["epochs"] = int(m.group(1))
    return metrics

def evaluate_urbannav(pos_file, csv_file):
    """Evaluate using gneiss-cli eval."""
    cmd = [CLI_BIN, "eval", "--solution", pos_file, "--truth", csv_file]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        output = proc.stdout + "\n" + proc.stderr
        return parse_eval_output(output)
    except Exception as e:
        return {"error": str(e)}

def parse_metric(m):
    """Parse a metric string like '0.840 m' to float."""
    m = m.replace('m', '').replace('°', '').strip()
    try:
        return float(m)
    except (ValueError, TypeError):
        return None

def main():
    results = {}

    # Evaluate IGS stations
    print("=" * 70)
    print("EVALUATING IGS STATIONS")
    print("=" * 70)

    for station in ["alic", "cedu", "yarr", "hob2", "park", "pert"]:
        r = evaluate_igs(station)
        if r:
            results[f"igs_{station}"] = r
            print(f"\n{station.upper()}:")
            print(f"  3D: median={r['median_3d']:.3f}m, p68={r['p68_3d']:.3f}m, p95={r['p95_3d']:.3f}m")
            print(f"  HZ: median={r['median_hz']:.3f}m, p95={r['p95_hz']:.3f}m")
            print(f"  Fix rate: {r['fix_rate']:.1f}% ({r['fixes']}/{r['total_epochs']})")
            print(f"  Coast epochs: {r['coast_epochs']}")

    # Evaluate NKLG with IONEX
    print("\n" + "=" * 70)
    print("EVALUATING NKLG WITH IONEX")
    print("=" * 70)
    r = evaluate_igs("nklg_ionex")
    if r:
        results["igs_nklg_ionex"] = r
        print(f"\nNKLG (IONEX):")
        print(f"  3D: median={r['median_3d']:.3f}m, p68={r['p68_3d']:.3f}m, p95={r['p95_3d']:.3f}m")
        print(f"  HZ: median={r['median_hz']:.3f}m, p95={r['p95_hz']:.3f}m")
        print(f"  Fix rate: {r['fix_rate']:.1f}%")
        print(f"  Coast epochs: {r['coast_epochs']}")

    # Evaluate Odaiba
    print("\n" + "=" * 70)
    print("EVALUATING ODAIBA")
    print("=" * 70)
    odaiba_pos = os.path.join(OUT_DIR, "urbannav_odaiba_ppp.pos")
    odaiba_ref = os.path.join(DATASETS_DIR, "urbannav", "tokyo", "Tokyo_Data", "Odaiba", "reference.csv")
    metrics = evaluate_urbannav(odaiba_pos, odaiba_ref)
    if metrics and "error" not in metrics:
        # Evaluate POS file
        positions = read_pos_file(odaiba_pos)
        fixes = sum(1 for _, _, _, q in positions if q == 1)
        coast_epochs = sum(1 for _, _, _, q in positions if q >= 5)

        # Also get coasting from log
        log_file = os.path.join(OUT_DIR, "urbannav_odaiba_ppp.log")
        coast_segments = 0
        if os.path.exists(log_file):
            with open(log_file) as f:
                coast_segments = len(re.findall(r"coasting", f.read()))

        d3_50 = parse_metric(metrics.get("3d_50"))
        d3_75 = parse_metric(metrics.get("3d_75"))
        d3_95 = parse_metric(metrics.get("3d_95"))
        hz_50 = parse_metric(metrics.get("hz_50"))
        hz_95 = parse_metric(metrics.get("hz_95"))
        vt_50 = parse_metric(metrics.get("vt_50"))

        results["urbannav_odaiba"] = {
            "median_3d": d3_50,
            "p68_3d": d3_75,  # Approximate 68% with 75th percentile
            "p95_3d": d3_95,
            "median_hz": hz_50,
            "p95_hz": hz_95,
            "median_vt": vt_50,
            "fix_rate": fixes / len(positions) * 100 if positions else 0,
            "fixes": fixes,
            "total_epochs": metrics.get("epochs", 0),
            "coast_epochs": coast_epochs,
            "coast_segments": coast_segments,
        }
        print(f"\nOdaiba:")
        print(f"  3D: median={d3_50:.3f}m, p68={d3_75:.3f}m, p95={d3_95:.3f}m")
        print(f"  HZ: median={hz_50:.3f}m, p95={hz_95:.3f}m")
        print(f"  Fix rate: {fixes}/{len(positions)}")
        print(f"  Coasts: {coast_epochs} epochs, {coast_segments} segments")

    # Evaluate Shinjuku
    print("\n" + "=" * 70)
    print("EVALUATING SHINJUKU")
    print("=" * 70)
    shinjuku_pos = os.path.join(OUT_DIR, "urbannav_shinjuku_ppp.pos")
    shinjuku_ref = os.path.join(DATASETS_DIR, "urbannav", "tokyo", "Tokyo_Data", "Shinjuku", "reference.csv")
    metrics = evaluate_urbannav(shinjuku_pos, shinjuku_ref)
    if metrics and "error" not in metrics:
        positions = read_pos_file(shinjuku_pos)
        fixes = sum(1 for _, _, _, q in positions if q == 1)
        coast_epochs = sum(1 for _, _, _, q in positions if q >= 5)

        log_file = os.path.join(OUT_DIR, "urbannav_shinjuku_ppp.log")
        coast_segments = 0
        if os.path.exists(log_file):
            with open(log_file) as f:
                coast_segments = len(re.findall(r"coasting", f.read()))

        d3_50 = parse_metric(metrics.get("3d_50"))
        d3_75 = parse_metric(metrics.get("3d_75"))
        d3_95 = parse_metric(metrics.get("3d_95"))
        hz_50 = parse_metric(metrics.get("hz_50"))
        hz_95 = parse_metric(metrics.get("hz_95"))
        vt_50 = parse_metric(metrics.get("vt_50"))

        results["urbannav_shinjuku"] = {
            "median_3d": d3_50,
            "p68_3d": d3_75,
            "p95_3d": d3_95,
            "median_hz": hz_50,
            "p95_hz": hz_95,
            "median_vt": vt_50,
            "fix_rate": fixes / len(positions) * 100 if positions else 0,
            "fixes": fixes,
            "total_epochs": metrics.get("epochs", 0),
            "coast_epochs": coast_epochs,
            "coast_segments": coast_segments,
        }
        print(f"\nShinjuku:")
        print(f"  3D: median={d3_50:.3f}m, p68={d3_75:.3f}m, p95={d3_95:.3f}m")
        print(f"  HZ: median={hz_50:.3f}m, p95={hz_95:.3f}m")
        print(f"  Fix rate: {fixes}/{len(positions)}")
        print(f"  Coasts: {coast_epochs} epochs, {coast_segments} segments")

    # ======================================================================
    # SUMMARY TABLE
    # ======================================================================
    print("\n" + "=" * 100)
    print("COMPREHENSIVE BENCHMARK MATRIX - SUMMARY")
    print("=" * 100)
    print(f"\n{'Dataset':<18} {'Mode':<32} {'3D 50%':<10} {'3D 68%':<10} {'3D 95%':<10} {'Coasts':<10} {'Fix%':<8}")
    print("-" * 100)

    for key in sorted(results.keys()):
        r = results[key]
        if key == "igs_nklg_ionex":
            ds = "IGS NKLG"
            mode = "GPS PPP+AR+IONEX"
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

        d3_50s = f"{r['median_3d']:.3f}m" if r.get('median_3d') else "N/A"
        d3_68s = f"{r['p68_3d']:.3f}m" if r.get('p68_3d') else "N/A"
        d3_95s = f"{r['p95_3d']:.3f}m" if r.get('p95_3d') else "N/A"
        coasts = str(r.get('coast_epochs', 'N/A'))
        fix = f"{r['fix_rate']:.1f}%" if r.get('fix_rate', 0) > 0 else "0%"

        print(f"{ds:<18} {mode:<32} {d3_50s:<10} {d3_68s:<10} {d3_95s:<10} {coasts:<10} {fix:<8}")

    # Save JSON
    serializable = {}
    for key, r in results.items():
        serializable[key] = {
            k: (round(v, 3) if isinstance(v, float) else v)
            for k, v in r.items()
        }
    with open(os.path.join(OUT_DIR, "benchmark_results.json"), 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f"\nResults saved to {os.path.join(OUT_DIR, 'benchmark_results.json')}")

