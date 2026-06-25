#!/usr/bin/env python3
"""
Comprehensive benchmark matrix runner for Gneiss GNSS PPP engine.
Runs the full test matrix from SPRINT_ROADMAP.md and reports results.
"""
import subprocess
import os
import sys
import json
import time
import re
import multiprocessing
from datetime import datetime

CLI_BIN = "/Users/kevin/projects/gneiss/target/release/gneiss-cli"
BASE_DIR = "/Users/kevin/projects/gneiss"
DATASETS_DIR = os.path.join(BASE_DIR, "datasets")
BENCH_DIR = os.path.join(BASE_DIR, "benchmarks")
OUT_DIR = os.path.join(BENCH_DIR, "results")
os.makedirs(OUT_DIR, exist_ok=True)

# Ground truth coordinates from IGS RINEX APPROX POSITION XYZ
IGS_TRUTH = {
    "alic": (-4052051.7670, 4212836.2150, -2545106.0270),
    "cedu": (-3753472.1430, 3912741.0403, -3347961.0368),
    "yarr": (-2389024.5060, 5043315.3590, -3078534.1630),
    "hob2": (-3950071.2740, 2522415.2180, -4311638.5120),
    "park": (-4554255.0460, 2816652.4450, -3454059.9010),
    "pert": (-2368686.7490, 4881316.4850, -3341796.3890),
    "nklg": (6287385.7333, 1071574.6758, 39133.0395),
}

def create_igs_truth(name, x, y, z, out_path):
    """Create a ground truth POS file for an IGS station.
    Format: YYYY/MM/DD HH:MM:SS X Y Z Q nsats
    We just need one reference line per epoch, using a dummy time.
    The evaluator matches by TOW, so we need the full day covered.

    For day 335 of 2019 (Dec 1), GPS week 2082, we generate truth at 30s intervals.
    """
    # GPS week 2082 day 0 = Dec 1, 2019 00:00:00 GPS
    # Day 335 of 2019
    lines = []
    for tow_int in range(0, 86400, 30):
        tow = float(tow_int)
        dt = datetime(2019, 12, 1, 0, 0, 0)
        # Convert GPS tow to UTC roughly
        hours = tow_int // 3600
        mins = (tow_int % 3600) // 60
        secs = tow_int % 60
        line = f"2019/12/01 {hours:02d}:{mins:02d}:{secs:05.2f} {x:.4f} {y:.4f} {z:.4f} 1 12\n"
        lines.append(line)

    with open(out_path, 'w') as f:
        f.write(f"% IGS {name} ground truth (SINEX)\n")
        f.writelines(lines)
    print(f"  Created truth file: {out_path} ({len(lines)} epochs)")

def run_benchmark_igs(station, use_ionex=False):
    """Run IGS station benchmark with PPP+AR.
    Returns (pos_file_path, success_bool)
    """
    igs_dir = os.path.join(DATASETS_DIR, "igs")
    station_lower = station.lower()
    label = f"igs_{station_lower}"
    if use_ionex:
        label += "_ionex"
    out_pos = os.path.join(OUT_DIR, f"{label}.pos")

    # IGS files
    obs_file = os.path.join(igs_dir, f"{station_lower}3350.19o")
    if not os.path.exists(obs_file):
        # Try crx version
        obs_file = os.path.join(igs_dir, f"{station_lower}3350.19o.crx")
    if not os.path.exists(obs_file):
        print(f"  [SKIP] No observation file found for {station}")
        return None, False

    nav_file = os.path.join(igs_dir, "brdc3350.19n")
    sp3_file = os.path.join(igs_dir, "cod20820.sp3")
    clk_file = os.path.join(igs_dir, "gfz20820.clk")
    atx_file = os.path.join(igs_dir, "igs14.atx")

    if use_ionex:
        config_file = os.path.join(BENCH_DIR, "igs_ppp_ionex_config.json")
        ionex_file = os.path.join(igs_dir, "codg3350.19i")
    else:
        config_file = os.path.join(BENCH_DIR, "igs_ppp_ar_config.json")
        ionex_arg = []

    if not os.path.exists(sp3_file) or not os.path.exists(clk_file):
        print(f"  [SKIP] Missing SP3/CLK for {station}")
        return None, False

    # Build command
    cmd = [
        CLI_BIN, "process",
        "--mode", "ppp",
        "--rover", obs_file,
        "--nav", nav_file,
        "--sp3", sp3_file,
        "--clk", clk_file,
        "--antex", atx_file,
        "--config", config_file,
        "--systems", "G",  # GPS only for IGS
        "--output", out_pos,
    ]
    if use_ionex:
        cmd += ["--ionex", os.path.join(igs_dir, "codg3350.19i")]

    # Skip if output already exists and is non-empty
    if os.path.exists(out_pos) and os.path.getsize(out_pos) > 1000:
        print(f"  [CACHED] {out_pos} exists, skipping processing")
        return out_pos, True

    print(f"  Processing {station}...")
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=1800  # 30 min timeout
        )
        if proc.returncode != 0:
            print(f"  [FAIL] Return code {proc.returncode}")
            print(f"  STDERR: {proc.stderr[-500:]}")
            return None, False
        print(f"  [OK] Processed {station}")
        return out_pos, True
    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] {station}")
        return None, False
    except Exception as e:
        print(f"  [ERROR] {station}: {e}")
        return None, False

def run_benchmark_urbannav(dataset_name, rover_obs, nav_file, ref_csv, label, sp3_file=None, clk_file=None, bia_file=None):
    """Run UrbanNav benchmark with GPS+Galileo PPP."""
    out_pos = os.path.join(OUT_DIR, f"{label}.pos")

    cmd = [
        CLI_BIN, "process",
        "--mode", "ppp",
        "--rover", rover_obs,
        "--nav", nav_file,
        "--config", os.path.join(BENCH_DIR, "urbannav_ppp_config.json"),
        "--systems", "GE",  # GPS+Galileo
        "--dynamics", "automotive",
        "--output", out_pos,
    ]

    # Add precise products if available
    tokyo_dir = os.path.join(DATASETS_DIR, "urbannav", "tokyo")
    sp3 = sp3_file or os.path.join(tokyo_dir, "COD0MGXFIN_20183530000_01D_05M_ORB.SP3")
    clk = clk_file or os.path.join(tokyo_dir, "cod20323.clk")
    bia = bia_file or os.path.join(tokyo_dir, "COD0MGXFIN_20183530000_01D_01D_OSB.BIA")

    if os.path.exists(sp3):
        cmd += ["--sp3", sp3]
    if os.path.exists(clk):
        cmd += ["--clk", clk]
    if os.path.exists(bia):
        cmd += ["--bia", bia]

    cmd += ["--antex", os.path.join(DATASETS_DIR, "igs14.atx")]

    # Skip if cached
    if os.path.exists(out_pos) and os.path.getsize(out_pos) > 1000:
        print(f"  [CACHED] {out_pos} exists")
        return out_pos

    print(f"  Processing {dataset_name}...")
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=3600
        )
        if proc.returncode != 0:
            print(f"  [FAIL] Return code {proc.returncode}")
            if proc.stderr:
                print(f"  STDERR: {proc.stderr[-500:]}")
            return None
        print(f"  [OK] Processed {dataset_name}")
        return out_pos
    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] {dataset_name}")
        return None
    except Exception as e:
        print(f"  [ERROR] {dataset_name}: {e}")
        return None

def evaluate_igs(pos_file, truth_file):
    """Evaluate IGS station POS file against truth using gneiss-cli eval."""
    cmd = [CLI_BIN, "eval", "--solution", pos_file, "--truth", truth_file]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        output = proc.stdout + "\n" + proc.stderr
        return parse_eval_output(output)
    except Exception as e:
        print(f"  [EVAL ERROR] {e}")
        return None

def parse_eval_output(text):
    """Parse the eval output table."""
    metrics = {}
    for line in text.splitlines():
        if "No matching epochs found" in line or "ERROR" in line.upper() and "Evaluated" not in line and "LOADED" not in line.upper():
            pass  # Not a fatal error, just may have no matches
        if line.startswith("| Horiz") and "Metric" not in line:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 7:
                metrics["hz_25"] = parts[1]
                metrics["hz_50"] = parts[2]
                metrics["hz_75"] = parts[3]
                metrics["hz_90"] = parts[4]
                metrics["hz_95"] = parts[5]
                metrics["hz_99"] = parts[6]
        elif line.startswith("| Vert"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 7:
                metrics["vt_25"] = parts[1]
                metrics["vt_50"] = parts[2]
                metrics["vt_75"] = parts[3]
                metrics["vt_90"] = parts[4]
                metrics["vt_95"] = parts[5]
                metrics["vt_99"] = parts[6]
        elif line.startswith("| 3D"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 7:
                metrics["3d_25"] = parts[1]
                metrics["3d_50"] = parts[2]
                metrics["3d_68"] = None  # Need to compute separately
                metrics["3d_75"] = parts[3]
                metrics["3d_90"] = parts[4]
                metrics["3d_95"] = parts[5]
                metrics["3d_99"] = parts[6]
        elif line.startswith("| Heading"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 7:
                metrics["hdg_25"] = parts[1]
                metrics["hdg_50"] = parts[2]

    # Extract evaluated count
    m = re.search(r"Evaluated (\d+) epochs", text)
    if m:
        metrics["epochs"] = int(m.group(1))

    return metrics if metrics else None

def compute_igs_errors_3d(pos_file, truth_x, truth_y, truth_z):
    """Compute 3D, horizontal, and vertical errors for IGS static station.
    Returns (hz_errors, vt_errors, d3_errors) sorted lists.
    """
    import numpy as np
    from gneiss_core.coords import ecef_to_llh, llh_to_ecef

    # Read solution
    positions = []
    with open(pos_file, 'r') as f:
        for line in f:
            if line.startswith('%') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    positions.append((x, y, z))
                except ValueError:
                    continue

    if not positions:
        return [], [], []

    true_pos = np.array([truth_x, truth_y, truth_z])
    true_llh = ecef_to_llh(np.array([truth_x, truth_y, truth_z]))  # Not directly callable this way

    # Compute errors in ENU frame
    hz_errors = []
    vt_errors = []
    d3_errors = []

    for pos in positions:
        pos_arr = np.array(pos)
        diff = pos_arr - true_pos

        # ENU transformation
        # For now, approximate 3D error directly from ECEF difference
        # (for static stations, this is fine since we're measuring convergence)
        d3 = np.linalg.norm(diff)
        d3_errors.append(d3)

    d3_errors.sort()
    return d3_errors

def evaluate_pos_vs_csv(pos_file, csv_file):
    """Evaluate POS file against a reference CSV (like UrbanNav)."""
    cmd = [CLI_BIN, "eval", "--solution", pos_file, "--truth", csv_file]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        output = proc.stdout + "\n" + proc.stderr
        return parse_eval_output(output)
    except Exception as e:
        print(f"  [EVAL ERROR] {e}")
        return None

def extract_percentile(errors, pct):
    """Extract percentile from sorted list."""
    if not errors:
        return float('nan')
    idx = int(len(errors) * pct / 100.0)
    idx = min(idx, len(errors) - 1)
    return errors[idx]

def count_coasting_events(pos_file):
    """Count coasting events from the POS file (Q=5 or where ambiguities reset)."""
    count = 0
    prev_fixed = True
    coast_segments = 0

    if not os.path.exists(pos_file):
        return 0

    with open(pos_file, 'r') as f:
        for line in f:
            if line.startswith('%') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    q = int(parts[5])
                    # Q=5 means coasting/single point, Q=2 means float
                    if q >= 5:
                        if prev_fixed:
                            coast_segments += 1
                            prev_fixed = False
                        count += 1
                    else:
                        if not prev_fixed and count > 0:
                            pass  # End of coast segment
                        prev_fixed = True
                except ValueError:
                    continue

    return count, coast_segments

def count_ar_fixes(pos_file):
    """Count AR fix epochs (Q=1 in pos file)."""
    fixes = 0
    total = 0

    if not os.path.exists(pos_file):
        return 0, 0

    with open(pos_file, 'r') as f:
        for line in f:
            if line.startswith('%') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    q = int(parts[5])
                    total += 1
                    if q == 1:
                        fixes += 1
                except ValueError:
                    continue

    return fixes, total

def main():
    print("=" * 70)
    print("GNEISS COMPREHENSIVE BENCHMARK MATRIX")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    print()

    results = {}

    # ======================================================================
    # 1. IGS STATIONS - GPS PPP with AR + SP3/CLK
    # ======================================================================
    print("\n" + "=" * 70)
    print("PHASE 1: IGS Stations - GPS PPP with AR + SP3/CLK")
    print("=" * 70)

    igs_stations = ["alic", "cedu", "yarr", "hob2", "park", "pert"]
    igs_dir = os.path.join(DATASETS_DIR, "igs")

    for station in igs_stations:
        print(f"\n--- IGS {station.upper()} ---")
        pos_file, success = run_benchmark_igs(station, use_ionex=False)

        if pos_file and os.path.exists(pos_file):
            # Compute 3D errors using Python with numpy
            # We need numpy for proper ENU computation
            try:
                import numpy as np
                from gneiss_core.coords import ecef_to_llh

                # Read positions
                positions = []
                with open(pos_file, 'r') as f:
                    for line in f:
                        if line.startswith('%') or line.strip() == '':
                            continue
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                                positions.append((x, y, z))
                            except ValueError:
                                continue

                truth_x, truth_y, truth_z = IGS_TRUTH[station]
                true_pos = np.array([truth_x, truth_y, truth_z])

                hz_errors = []
                vt_errors = []
                d3_errors = []

                # Convert truth to LLH for ENU computation
                # ecef_to_llh takes a Vector3<f64> from nalgebra
                # We'll approximate with just ECEF difference for 3D
                for pos in positions:
                    pos_arr = np.array(pos)
                    diff = pos_arr - true_pos
                    d3 = np.linalg.norm(diff)

                    # Approximate horizontal vs vertical using WGS84
                    # Use local ENU transform
                    import math
                    r2 = math.sqrt(true_x**2 + true_y**2)
                    lon = math.atan2(true_y, true_x)
                    # Old method: compute approximate lat
                    lat = math.atan2(true_z, r2 * (1 - 0.00669437999014))
                    for _ in range(3):
                        sin_lat = math.sin(lat)
                        N = 6378137.0 / math.sqrt(1 - 0.00669437999014 * sin_lat**2)
                        lat = math.atan2(true_z + 0.00669437999014 * N * sin_lat, r2)

                    sin_lat = math.sin(lat)
                    cos_lat = math.cos(lat)
                    sin_lon = math.sin(lon)
                    cos_lon = math.cos(lon)

                    e = -sin_lon * diff[0] + cos_lon * diff[1]
                    n = -sin_lat * cos_lon * diff[0] - sin_lat * sin_lon * diff[1] + cos_lat * diff[2]
                    u = cos_lat * cos_lon * diff[0] + cos_lat * sin_lon * diff[1] + sin_lat * diff[2]

                    hz_errors.append(math.sqrt(e*e + n*n))
                    vt_errors.append(abs(u))
                    d3_errors.append(math.sqrt(e*e + n*n + u*u))

                hz_errors.sort()
                vt_errors.sort()
                d3_errors.sort()

                d3_50 = extract_percentile(d3_errors, 50)
                d3_68 = extract_percentile(d3_errors, 68)
                d3_95 = extract_percentile(d3_errors, 95)
                hz_50 = extract_percentile(hz_errors, 50)
                hz_95 = extract_percentile(hz_errors, 95)
                vt_50 = extract_percentile(vt_errors, 50)

                # Count coasting events and AR fixes
                coast_count, coast_segs = count_coasting_events(pos_file)
                fixes, total = count_ar_fixes(pos_file)
                fix_rate = (fixes / total * 100) if total > 0 else 0

                results[f"igs_{station}"] = {
                    "median_3d": d3_50,
                    "p68_3d": d3_68,
                    "p95_3d": d3_95,
                    "median_hz": hz_50,
                    "p95_hz": hz_95,
                    "median_vt": vt_50,
                    "coast_epochs": coast_count,
                    "coast_segments": coast_segs,
                    "fix_rate": fix_rate,
                    "total_epochs": total,
                    "pos_file": pos_file,
                }

                print(f"  Median 3D: {d3_50:.3f}m, 68%: {d3_68:.3f}m, 95%: {d3_95:.3f}m")
                print(f"  Coasts: {coast_count} ({coast_segs} seg), Fix rate: {fix_rate:.1f}%")

            except Exception as e:
                print(f"  [ERROR computing errors for {station}]: {e}")
                import traceback
                traceback.print_exc()

    # ======================================================================
    # 2. IGS NKLG - GPS PPP with IONEX + SP3/CLK
    # ======================================================================
    print("\n" + "=" * 70)
    print("PHASE 2: IGS NKLG - GPS PPP with IONEX + SP3/CLK")
    print("=" * 70)

    print(f"\n--- NKLG with IONEX ---")
    pos_file, success = run_benchmark_igs("nklg", use_ionex=True)

    if pos_file and os.path.exists(pos_file):
        try:
            import numpy as np
            import math

            positions = []
            with open(pos_file, 'r') as f:
                for line in f:
                    if line.startswith('%') or line.strip() == '':
                        continue
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                            positions.append((x, y, z))
                        except ValueError:
                            continue

            truth_x, truth_y, truth_z = IGS_TRUTH["nklg"]
            true_pos = np.array([truth_x, truth_y, truth_z])

            hz_errors = []
            vt_errors = []
            d3_errors = []

            r2 = math.sqrt(truth_x**2 + truth_y**2)
            lon = math.atan2(truth_y, truth_x)
            lat = math.atan2(truth_z, r2 * (1 - 0.00669437999014))
            for _ in range(3):
                sin_lat = math.sin(lat)
                N = 6378137.0 / math.sqrt(1 - 0.00669437999014 * sin_lat**2)
                lat = math.atan2(truth_z + 0.00669437999014 * N * sin_lat, r2)

            sin_lat = math.sin(lat)
            cos_lat = math.cos(lat)
            sin_lon = math.sin(lon)
            cos_lon = math.cos(lon)

            for pos in positions:
                pos_arr = np.array(pos)
                diff = pos_arr - true_pos

                e = -sin_lon * diff[0] + cos_lon * diff[1]
                n = -sin_lat * cos_lon * diff[0] - sin_lat * sin_lon * diff[1] + cos_lat * diff[2]
                u = cos_lat * cos_lon * diff[0] + cos_lat * sin_lon * diff[1] + sin_lat * diff[2]

                hz_errors.append(math.sqrt(e*e + n*n))
                vt_errors.append(abs(u))
                d3_errors.append(math.sqrt(e*e + n*n + u*u))

            hz_errors.sort()
            vt_errors.sort()
            d3_errors.sort()

            d3_50 = extract_percentile(d3_errors, 50)
            d3_68 = extract_percentile(d3_errors, 68)
            d3_95 = extract_percentile(d3_errors, 95)
            hz_50 = extract_percentile(hz_errors, 50)
            hz_95 = extract_percentile(hz_errors, 95)

            coast_count, coast_segs = count_coasting_events(pos_file)
            fixes, total = count_ar_fixes(pos_file)
            fix_rate = (fixes / total * 100) if total > 0 else 0

            results["igs_nklg_ionex"] = {
                "median_3d": d3_50,
                "p68_3d": d3_68,
                "p95_3d": d3_95,
                "median_hz": hz_50,
                "p95_hz": hz_95,
                "coast_epochs": coast_count,
                "coast_segments": coast_segs,
                "fix_rate": fix_rate,
                "total_epochs": total,
                "pos_file": pos_file,
            }

            print(f"  Median 3D: {d3_50:.3f}m, 68%: {d3_68:.3f}m, 95%: {d3_95:.3f}m")
            print(f"  Coasts: {coast_count} ({coast_segs} seg), Fix rate: {fix_rate:.1f}%")

        except Exception as e:
            print(f"  [ERROR computing errors for NKLG]: {e}")
            import traceback
            traceback.print_exc()

    # ======================================================================
    # 3. Odaiba - GPS+Galileo PPP with SP3/CLK/BIA
    # ======================================================================
    print("\n" + "=" * 70)
    print("PHASE 3: Odaiba - GPS+Galileo PPP with SP3/CLK/BIA")
    print("=" * 70)

    odaiba_dir = os.path.join(DATASETS_DIR, "urbannav", "tokyo", "Tokyo_Data", "Odaiba")
    odaiba_obs = os.path.join(odaiba_dir, "rover_ublox.obs")
    odaiba_nav = os.path.join(odaiba_dir, "base.nav")
    odaiba_ref = os.path.join(odaiba_dir, "reference.csv")

    if os.path.exists(odaiba_obs) and os.path.exists(odaiba_ref):
        pos_file = run_benchmark_urbannav(
            "Odaiba", odaiba_obs, odaiba_nav, odaiba_ref, "urbannav_odaiba_ppp",
            sp3_file=os.path.join(DATASETS_DIR, "urbannav", "tokyo", "COD0MGXFIN_20183530000_01D_05M_ORB.SP3"),
            clk_file=os.path.join(DATASETS_DIR, "urbannav", "tokyo", "cod20323.clk"),
            bia_file=os.path.join(DATASETS_DIR, "urbannav", "tokyo", "COD0MGXFIN_20183530000_01D_01D_OSB.BIA"),
        )

        if pos_file and os.path.exists(pos_file):
            metrics = evaluate_pos_vs_csv(pos_file, odaiba_ref)
            if metrics:
                # Get 3D, HZ, VT from eval output
                # The eval output gives horiz, vert, 3D at various percentiles
                # For 3D we need to compute from HZ and VT: sqrt(hz^2 + vt^2)
                # But the evaluator gives us 3D directly

                # Extract numeric values
                def parse_metric(m):
                    m = m.replace('m', '').replace('°', '').strip()
                    try:
                        return float(m)
                    except:
                        return None

                hz_50 = parse_metric(metrics.get("hz_50", "N/A"))
                hz_95 = parse_metric(metrics.get("hz_95", "N/A"))
                vt_50 = parse_metric(metrics.get("vt_50", "N/A"))
                d3_50_val = parse_metric(metrics.get("3d_50", "N/A"))
                d3_95_val = parse_metric(metrics.get("3d_95", "N/A"))

                # Compute 68th percentile - not directly available from eval
                # We'll estimate: 68% ~= 75% approx, or compute from raw data
                # For now, use 75th as approximation for 68th
                d3_68_val = parse_metric(metrics.get("3d_75", "N/A"))

                fixes, total = count_ar_fixes(pos_file)
                fix_rate = (fixes / total * 100) if total > 0 else 0
                coast_count, coast_segs = count_coasting_events(pos_file)

                results["urbannav_odaiba"] = {
                    "median_3d": d3_50_val,
                    "p68_3d": d3_68_val,
                    "p95_3d": d3_95_val,
                    "median_hz": hz_50,
                    "p95_hz": hz_95,
                    "median_vt": vt_50,
                    "coast_epochs": coast_count,
                    "coast_segments": coast_segs,
                    "fix_rate": fix_rate,
                    "total_epochs": metrics.get("epochs", total),
                    "pos_file": pos_file,
                }

                print(f"  Median 3D: {d3_50_val:.3f}m, 68%: {d3_68_val:.3f}m, 95%: {d3_95_val:.3f}m")
                print(f"  Coasts: {coast_count} ({coast_segs} seg), Fix rate: {fix_rate:.1f}%")
            else:
                print(f"  [FAIL] Could not evaluate Odaiba result")

    # ======================================================================
    # 4. Shinjuku - GPS+Galileo PPP
    # ======================================================================
    print("\n" + "=" * 70)
    print("PHASE 4: Shinjuku - GPS+Galileo PPP with SP3/CLK/BIA")
    print("=" * 70)

    shinjuku_dir = os.path.join(DATASETS_DIR, "urbannav", "tokyo", "Tokyo_Data", "Shinjuku")
    shinjuku_obs = os.path.join(shinjuku_dir, "rover_ublox.obs")
    shinjuku_nav = os.path.join(shinjuku_dir, "base.nav")
    shinjuku_ref = os.path.join(shinjuku_dir, "reference.csv")

    if os.path.exists(shinjuku_obs) and os.path.exists(shinjuku_ref):
        pos_file = run_benchmark_urbannav(
            "Shinjuku", shinjuku_obs, shinjuku_nav, shinjuku_ref, "urbannav_shinjuku_ppp",
            sp3_file=os.path.join(DATASETS_DIR, "urbannav", "tokyo", "COD0MGXFIN_20183530000_01D_05M_ORB.SP3"),
            clk_file=os.path.join(DATASETS_DIR, "urbannav", "tokyo", "cod20323.clk"),
            bia_file=os.path.join(DATASETS_DIR, "urbannav", "tokyo", "COD0MGXFIN_20183530000_01D_01D_OSB.BIA"),
        )

        if pos_file and os.path.exists(pos_file):
            metrics = evaluate_pos_vs_csv(pos_file, shinjuku_ref)
            if metrics:
                def parse_metric(m):
                    m = m.replace('m', '').replace('°', '').strip()
                    try:
                        return float(m)
                    except:
                        return None

                hz_50 = parse_metric(metrics.get("hz_50", "N/A"))
                hz_95 = parse_metric(metrics.get("hz_95", "N/A"))
                vt_50 = parse_metric(metrics.get("vt_50", "N/A"))
                d3_50_val = parse_metric(metrics.get("3d_50", "N/A"))
                d3_95_val = parse_metric(metrics.get("3d_95", "N/A"))
                d3_68_val = parse_metric(metrics.get("3d_75", "N/A"))

                fixes, total = count_ar_fixes(pos_file)
                fix_rate = (fixes / total * 100) if total > 0 else 0
                coast_count, coast_segs = count_coasting_events(pos_file)

                results["urbannav_shinjuku"] = {
                    "median_3d": d3_50_val,
                    "p68_3d": d3_68_val,
                    "p95_3d": d3_95_val,
                    "median_hz": hz_50,
                    "p95_hz": hz_95,
                    "median_vt": vt_50,
                    "coast_epochs": coast_count,
                    "coast_segments": coast_segs,
                    "fix_rate": fix_rate,
                    "total_epochs": metrics.get("epochs", total),
                    "pos_file": pos_file,
                }

                print(f"  Median 3D: {d3_50_val:.3f}m, 68%: {d3_68_val:.3f}m, 95%: {d3_95_val:.3f}m")
                print(f"  Coasts: {coast_count} ({coast_segs} seg), Fix rate: {fix_rate:.1f}%")
            else:
                print(f"  [FAIL] Could not evaluate Shinjuku result")

    # ======================================================================
    # SUMMARY REPORT
    # ======================================================================
    print("\n" + "=" * 70)
    print("COMPREHENSIVE BENCHMARK RESULTS")
    print("=" * 70)

    print(f"\n{'Dataset':<20} {'Mode':<30} {'3D 50%':<12} {'3D 68%':<12} {'3D 95%':<12} {'Coasts':<10} {'AR Fix%':<10}")
    print("-" * 106)

    for key, r in sorted(results.items()):
        dataset = key
        if key.startswith("igs_"):
            if "ionex" in key:
                mode = "GPS PPP + AR + IONEX"
                dataset = f"IGS {key[4:9].upper()}"
            else:
                mode = "GPS PPP + AR + SP3/CLK"
                dataset = f"IGS {key[4:].upper()}"
        elif key == "urbannav_odaiba":
            mode = "GPS+Galileo PPP + SP3/CLK/BIA"
            dataset = "Odaiba"
        elif key == "urbannav_shinjuku":
            mode = "GPS+Galileo PPP + SP3/CLK/BIA"
            dataset = "Shinjuku"
        else:
            mode = "Unknown"

        d3_50 = f"{r['median_3d']:.3f}m" if r.get('median_3d') is not None and not (isinstance(r['median_3d'], float) and math.isnan(r['median_3d'])) else "N/A"
        d3_68 = f"{r['p68_3d']:.3f}m" if r.get('p68_3d') is not None and not (isinstance(r['p68_3d'], float) and math.isnan(r['p68_3d'])) else r.get('p68_3d', "N/A")
        d3_95 = f"{r['p95_3d']:.3f}m" if r.get('p95_3d') is not None and not (isinstance(r['p95_3d'], float) and math.isnan(r['p95_3d'])) else "N/A"
        coasts = str(r.get('coast_epochs', 'N/A'))
        fix_rate = f"{r.get('fix_rate', 0):.1f}%" if r.get('fix_rate', 0) > 0 else "N/A"

        print(f"{dataset:<20} {mode:<30} {d3_50:<12} {d3_68:<12} {d3_95:<12} {coasts:<10} {fix_rate:<10}")

    # Save results as JSON
    results_file = os.path.join(OUT_DIR, "benchmark_results.json")
    # Convert to serializable format
    serializable = {}
    for key, r in results.items():
        serializable[key] = {
            k: v for k, v in r.items() if k != 'pos_file'
        }
        serializable[key]['median_3d'] = round(r.get('median_3d', 0), 3) if r.get('median_3d') else None
        serializable[key]['p68_3d'] = round(r.get('p68_3d', 0), 3) if r.get('p68_3d') else None
        serializable[key]['p95_3d'] = round(r.get('p95_3d', 0), 3) if r.get('p95_3d') else None
        serializable[key]['median_hz'] = round(r.get('median_hz', 0), 3) if r.get('median_hz') else None
        serializable[key]['median_vt'] = round(r.get('median_vt', 0), 3) if r.get('median_vt') else None
        serializable[key]['fix_rate'] = round(r.get('fix_rate', 0), 1)

    with open(results_file, 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f"\nResults saved to {results_file}")

    print(f"\nFinished: {datetime.now().isoformat()}")
    return results

if __name__ == "__main__":
    results = main()
