#!/usr/bin/env python3
"""Run RTKLIB demo5 on the same datasets as Gneiss and compare results."""
import os
import subprocess
import re
import sys

RTKLIB_BIN = "/tmp/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp"
GNEISS_CLI = "target/release/gneiss-cli"
OUT_DIR = "benchmarks/rtklib_comparison"

# RTKLIB modes: -p flag
# 0 = single (SPP)
# 2 = kinematic (RTK)
# 3 = static (RTK static)
# 6 = ppp-kinematic
# 7 = ppp-static

DATASETS = {
    "PPP (f9p_ppp)": {
        "rover": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs",
        "base": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/tmg23590.20o",
        "nav": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav",
        "truth": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppk.pos",
        "conf": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ppk.conf",
        "gneiss_config": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json",
        "is_static": False,
    },
    "GSDC (Pixel 4)": {
        "rover": "datasets/gsdc/Pixel4_GnssLog.20o",
        "base": "datasets/gsdc/p2221350.20o",
        "nav": "datasets/gsdc/rover.nav",
        "truth": "datasets/gsdc/reference.csv",
        "conf": None,
        "gneiss_config": "datasets/gsdc/gsdc_config.json",
        "is_static": True,  # GSDC is a static dataset
    },
    "Shinjuku (u-blox)": {
        "rover": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
        "base": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
        "nav": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav",
        "truth": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv",
        "conf": None,
        "gneiss_config": "datasets/urbannav/tokyo/tokyo_config.json",
        "is_static": False,
    },
    "Odaiba (u-blox)": {
        "rover": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/rover_ublox.obs",
        "base": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base_trimble.obs",
        "nav": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base.nav",
        "truth": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/reference.csv",
        "conf": None,
        "gneiss_config": "datasets/urbannav/tokyo/tokyo_config.json",
        "is_static": False,
    },
}

# Comparison matrix: each entry = (label, rtklib_args, gneiss_mode, gneiss_flags)
COMPARISONS = [
    ("SPP", ["-p", "0"], "spp", []),
    ("RTK Kinematic", ["-p", "2", "-h"], "rtk", []),
    ("RTK Kinematic (combined)", ["-p", "2", "-c", "-h"], "rtk", ["--enable-backward-smoothing"]),
    ("PPP Kinematic", ["-p", "6"], "ppp", []),
    ("PPP Kinematic (combined)", ["-p", "6", "-c"], "ppp", ["--enable-backward-smoothing"]),
]


def parse_eval_metrics(text):
    """Parse gneiss-cli eval output for metric lines."""
    metrics = {}
    for line in text.splitlines():
        if line.startswith("| Horiz"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) > 4:
                metrics["hz_50"] = parts[2]
                metrics["hz_95"] = parts[3]
                metrics["hz_99"] = parts[4]
        elif line.startswith("| Vert"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) > 4:
                metrics["vt_50"] = parts[2]
                metrics["vt_95"] = parts[3]
                metrics["vt_99"] = parts[4]
        elif line.startswith("| 3D"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) > 4:
                metrics["3d_50"] = parts[2]
                metrics["3d_95"] = parts[3]
    return metrics


def extract_float(s):
    """Extract float from string like '0.501 m'"""
    m = re.search(r"([\d.]+)", s)
    if m:
        return float(m.group(1))
    return float("inf")


def run_rtklib(dataset_name, ds_config, rtklib_args, label):
    """Run RTKLIB rnx2rtkp with given arguments."""
    safe_name = dataset_name.replace(" ", "_").replace("(", "").replace(")", "")
    safe_label = label.replace(" ", "_").replace("(", "").replace(")", "")
    out_file = os.path.join(OUT_DIR, f"rtklib_{safe_name}_{safe_label}.pos")

    cmd = [RTKLIB_BIN] + rtklib_args + ["-e", "-t"]

    # Use config file if available
    if ds_config.get("conf"):
        cmd += ["-k", ds_config["conf"]]

    cmd += [ds_config["rover"]]

    # For RTK modes, add base station
    if any(p in rtklib_args for p in ["1", "2", "3", "4"]) and ds_config.get("base"):
        cmd += [ds_config["base"]]

    # Add nav file
    cmd += [ds_config["nav"]]

    print(f"    RTKLIB cmd: {' '.join(cmd)}")

    try:
        with open(out_file, "w") as f:
            proc = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, timeout=600)
        if proc.returncode != 0:
            print(f"    RTKLIB failed (exit {proc.returncode})")
            return None, out_file
    except subprocess.TimeoutExpired:
        print(f"    RTKLIB timeout")
        return None, out_file

    # Count output lines
    with open(out_file) as f:
        lines = [l for l in f if not l.startswith("%")]
    print(f"    RTKLIB produced {len(lines)} solution epochs")

    if len(lines) < 10:
        print(f"    RTKLIB produced too few solutions")
        return None, out_file

    return out_file, out_file


def run_gneiss(dataset_name, ds_config, gneiss_mode, gneiss_flags):
    """Run Gneiss with given mode."""
    safe_name = dataset_name.replace(" ", "_").replace("(", "").replace(")", "")
    out_file = os.path.join(OUT_DIR, f"gneiss_{safe_name}_{gneiss_mode}.pos")

    cmd = [
        GNEISS_CLI, "process",
        "--mode", gneiss_mode,
        "--rover", ds_config["rover"],
        "--output", out_file,
    ]
    if ds_config.get("base"):
        cmd += ["--base", ds_config["base"]]
    if ds_config.get("nav"):
        cmd += ["--nav", ds_config["nav"]]
    if ds_config.get("gneiss_config"):
        cmd += ["--config", ds_config["gneiss_config"]]
    cmd += gneiss_flags

    print(f"    Gneiss cmd: {' '.join(cmd)}")

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            print(f"    Gneiss failed (exit {proc.returncode})")
            return None
    except subprocess.TimeoutExpired:
        print(f"    Gneiss timeout")
        return None

    return out_file


def evaluate(sol_file, truth_file):
    """Run gneiss-cli eval and return metrics."""
    if not sol_file or not os.path.exists(sol_file):
        return None

    cmd = [GNEISS_CLI, "eval", "--solution", sol_file, "--truth", truth_file]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    output = proc.stdout + "\n" + proc.stderr
    if "No matching epochs found" in output:
        return {"error": "Mismatch"}
    metrics = parse_eval_metrics(output)
    if not metrics:
        return {"error": "Parse Error"}
    return metrics


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Build Gneiss first
    print("Building Gneiss...")
    subprocess.run(["cargo", "build", "--release", "--bin", "gneiss-cli"], check=True)

    all_results = {}

    for ds_name, ds_config in DATASETS.items():
        print(f"\n{'='*60}")
        print(f"  Dataset: {ds_name}")
        print(f"{'='*60}")

        all_results[ds_name] = {}

        for label, rtklib_args, gneiss_mode, gneiss_flags in COMPARISONS:
            print(f"\n  --- {label} ---")

            # Check if this is an RTK mode and whether base station is available
            needs_base = any(p in rtklib_args for p in ["1", "2", "3", "4"])
            if needs_base and not ds_config.get("base"):
                print(f"    Skipping (no base station for this dataset)")
                continue

            # Run RTKLIB
            print(f"  Running RTKLIB ({label})...")
            rtklib_sol, _ = run_rtklib(ds_name, ds_config, rtklib_args, label)
            rtklib_metrics = evaluate(rtklib_sol, ds_config["truth"])

            # Run Gneiss
            print(f"  Running Gneiss ({gneiss_mode})...")
            gneiss_sol = run_gneiss(ds_name, ds_config, gneiss_mode, gneiss_flags)
            gneiss_metrics = evaluate(gneiss_sol, ds_config["truth"])

            all_results[ds_name][label] = {
                "rtklib": rtklib_metrics,
                "gneiss": gneiss_metrics,
            }

            # Print comparison
            if rtklib_metrics and "error" not in rtklib_metrics:
                print(f"    RTKLIB  -> Hz50: {rtklib_metrics.get('hz_50', 'N/A')}, Hz95: {rtklib_metrics.get('hz_95', 'N/A')}, Vt50: {rtklib_metrics.get('vt_50', 'N/A')}")
            else:
                print(f"    RTKLIB  -> {rtklib_metrics}")

            if gneiss_metrics and "error" not in gneiss_metrics:
                print(f"    Gneiss  -> Hz50: {gneiss_metrics.get('hz_50', 'N/A')}, Hz95: {gneiss_metrics.get('hz_95', 'N/A')}, Vt50: {gneiss_metrics.get('vt_50', 'N/A')}")
            else:
                print(f"    Gneiss  -> {gneiss_metrics}")

            # Compare
            if (rtklib_metrics and "error" not in rtklib_metrics and
                gneiss_metrics and "error" not in gneiss_metrics):
                r_hz = extract_float(rtklib_metrics.get("hz_50", "inf"))
                g_hz = extract_float(gneiss_metrics.get("hz_50", "inf"))
                if g_hz < r_hz:
                    pct = (1 - g_hz / r_hz) * 100 if r_hz > 0 else 0
                    print(f"    >>> Gneiss WINS by {pct:.1f}% (Hz50)")
                elif r_hz < g_hz:
                    pct = (1 - r_hz / g_hz) * 100 if g_hz > 0 else 0
                    print(f"    >>> RTKLIB WINS by {pct:.1f}% (Hz50) -- NEEDS IMPROVEMENT")
                else:
                    print(f"    >>> TIE (Hz50)")

    # Generate comparison markdown
    generate_comparison_md(all_results)

    return all_results


def generate_comparison_md(all_results):
    md = "# Gneiss vs RTKLIB (demo5) Comparison\n\n"
    md += "Side-by-side comparison on identical datasets and truth references.\n\n"

    for ds_name, modes in all_results.items():
        md += f"## {ds_name}\n\n"
        md += "| Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |\n"
        md += "|:-----|:-------|:--------|:--------|:--------|:-------|\n"

        for label, data in modes.items():
            r = data["rtklib"]
            g = data["gneiss"]

            r_hz50 = r.get("hz_50", "N/A") if r and "error" not in r else (r.get("error", "Failed") if r else "Failed")
            r_hz95 = r.get("hz_95", "N/A") if r and "error" not in r else ""
            r_vt50 = r.get("vt_50", "N/A") if r and "error" not in r else ""

            g_hz50 = g.get("hz_50", "N/A") if g and "error" not in g else (g.get("error", "Failed") if g else "Failed")
            g_hz95 = g.get("hz_95", "N/A") if g and "error" not in g else ""
            g_vt50 = g.get("vt_50", "N/A") if g and "error" not in g else ""

            # Determine winner
            winner = ""
            if (r and "error" not in r and g and "error" not in g):
                r_val = extract_float(r.get("hz_50", "inf"))
                g_val = extract_float(g.get("hz_50", "inf"))
                if g_val < r_val:
                    winner = "**Gneiss** ✅"
                elif r_val < g_val:
                    winner = "RTKLIB ⚠️"
                else:
                    winner = "Tie"

            md += f"| {label} | RTKLIB | {r_hz50} | {r_hz95} | {r_vt50} | {winner} |\n"
            md += f"| {label} | Gneiss | {g_hz50} | {g_hz95} | {g_vt50} | |\n"

        md += "\n"

    with open("COMPARISON.md", "w") as f:
        f.write(md)
    print(f"\nWrote COMPARISON.md")


if __name__ == "__main__":
    results = main()
