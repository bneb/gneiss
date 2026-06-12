#!/usr/bin/env python3
"""Unified Benchmark Orchestrator for Gneiss and RTKLIB comparisons."""
import os
import subprocess
import re
import argparse
import sys

RTKLIB_BIN = "/tmp/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp"
GNEISS_CLI = "target/release/gneiss-cli"
OUT_DIR_GNEISS = "benchmarks"
OUT_DIR_RTKLIB = "benchmarks/rtklib_comparison"
OUT_DIR_MATRIX = "benchmarks/matrix"

# Master Datasets dictionary
DATASETS = {
    "GSDC (Pixel 4)": {
        "rover": "datasets/gsdc/Pixel4_GnssLog.20o",
        "base": "datasets/gsdc/p2221350.20o",
        "nav": "datasets/gsdc/rover.nav",
        "truth": "datasets/gsdc/reference.csv",
        "gneiss_config": "datasets/gsdc/gsdc_config.json",
        "is_static": True,
    },
    "Shinjuku (UrbanNav)": {
        "rover": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
        "base": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
        "nav": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav",
        "truth": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv",
        "gneiss_config": "datasets/urbannav/tokyo/tokyo_config.json",
        "is_static": False,
    },
    "Odaiba (UrbanNav)": {
        "rover": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/rover_ublox.obs",
        "base": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base_trimble.obs",
        "nav": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base.nav",
        "truth": "datasets/urbannav/tokyo/Tokyo_Data/Odaiba/reference.csv",
        "gneiss_config": "datasets/urbannav/tokyo/tokyo_config.json",
        "is_static": False,
    },
    "PPP (f9p_ppp)": {
        "rover": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs",
        "base": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/tmg23590.20o",
        "nav": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav",
        "truth": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppk.pos",
        "sp3": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_05M_ORB.SP3",
        "clk": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_30S_CLK.CLK",
        "gneiss_config": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json",
        "conf": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ppk.conf",
        "is_static": False,
    },
    "UrbanLoco (Example)": {
        "rover": "datasets/urbanloco/CA/rover.obs",
        "base": "datasets/urbanloco/CA/base.obs",
        "nav": "datasets/urbanloco/CA/base.nav",
        "truth": "datasets/urbanloco/CA/ground_truth.csv",
        "gneiss_config": None,
        "is_static": False,
    },
    "TEX-CUP (UT Austin)": {
        "rover": "datasets/tex_cup/rover.obs",
        "base": "datasets/tex_cup/base.obs",
        "nav": "datasets/tex_cup/nav.nav",
        "truth": "datasets/tex_cup/truth.csv",
        "gneiss_config": None,
        "is_static": False,
    },
    "WHU-Smartphone (Xiaomi)": {
        "rover": "datasets/whu_smartphone/xiaomi/rover.obs",
        "base": "datasets/whu_smartphone/base/base.obs",
        "nav": "datasets/whu_smartphone/base/brdc.nav",
        "truth": "datasets/whu_smartphone/xiaomi/truth.csv",
        "gneiss_config": None,
        "is_static": False,
    },
    "smartLoc (TU Chemnitz)": {
        "rover": "datasets/smartloc/frankfurt/rover.obs",
        "base": "datasets/smartloc/frankfurt/base.obs",
        "nav": "datasets/smartloc/frankfurt/nav.nav",
        "truth": "datasets/smartloc/frankfurt/truth.csv",
        "gneiss_config": None,
        "is_static": False,
    }
}

# Shared Utilities
def extract_float(s):
    m = re.search(r"([\d.]+)", str(s))
    if m: return float(m.group(1))
    return float("inf")

def parse_eval_metrics(text):
    metrics = {}
    for line in text.splitlines():
        if line.startswith("| Horiz"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) > 4:
                metrics["hz_50"] = parts[2]
                metrics["hz_95"] = parts[3]
                if len(parts) > 5: metrics["hz_99"] = parts[4]
        elif line.startswith("| Vert"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) > 4:
                metrics["vt_50"] = parts[2]
                metrics["vt_95"] = parts[3]
                if len(parts) > 5: metrics["vt_99"] = parts[4]
    return metrics

def evaluate(sol_file, truth_file, dry_run=False):
    if dry_run:
        return {"hz_50": "0.0m", "hz_95": "0.0m", "vt_50": "0.0m"}
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

# Suite 1: Gneiss Only
def run_gneiss_suite(dry_run=False, eval_only=False):
    modes = ["spp", "ppp", "ppp-fg", "ppp-ins-fg"]
    results = {ds: {} for ds in DATASETS}
    os.makedirs(OUT_DIR_GNEISS, exist_ok=True)
    
    for ds_name, config in DATASETS.items():
        print(f"\n=== Gneiss Suite: {ds_name} ===")
        for mode in modes:
            out_file = f"{OUT_DIR_GNEISS}/{ds_name.replace(' ', '_').replace('(', '').replace(')', '')}_{mode}.pos"
            
            if not eval_only:
                cmd = [GNEISS_CLI, "process", "--mode", mode, "--rover", config["rover"], "--output", out_file]
                if config.get("base"): cmd.extend(["--base", config["base"]])
                if config.get("nav"): cmd.extend(["--nav", config["nav"]])
                if config.get("sp3"): cmd.extend(["--sp3", config["sp3"]])
                if config.get("clk"): cmd.extend(["--clk", config["clk"]])
                if config.get("gneiss_config"): cmd.extend(["--config", config["gneiss_config"]])
                if mode in ["ppp-fg", "ppp-ins-fg"]: cmd.append("--enable-backward-smoothing")
                
                if dry_run:
                    print(f"[DRY] {' '.join(cmd)}")
                else:
                    try:
                        subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                    except subprocess.TimeoutExpired:
                        results[ds_name][mode] = "Timeout"
                        continue

            if not os.path.exists(config["truth"]):
                results[ds_name][mode] = "No Truth"
                continue
                
            metrics = evaluate(out_file, config["truth"], dry_run)
            results[ds_name][mode] = metrics
            
    # Generate Markdown
    md = "# Gneiss Comprehensive Benchmarks\n\n"
    for ds_name, modes_data in results.items():
        md += f"## {ds_name}\n\n| Mode | Median Horizontal | 95% Horizontal | Median Vertical |\n| :--- | :--- | :--- | :--- |\n"
        for mode in modes:
            res = modes_data.get(mode, "N/A")
            if isinstance(res, dict):
                md += f"| `{mode}` | {res.get('hz_50','N/A')} | {res.get('hz_95','N/A')} | {res.get('vt_50','N/A')} |\n"
            else:
                md += f"| `{mode}` | {res} | {res} | {res} |\n"
        md += "\n"
    with open("BENCHMARKS.md", "w") as f: f.write(md)
    print("Wrote BENCHMARKS.md")

# Suite 2: RTKLIB Comparison
def run_rtklib_suite(dry_run=False, eval_only=False):
    comparisons = [
        ("SPP", ["-p", "0"], "spp", []),
        ("RTK Kinematic", ["-p", "2", "-h"], "rtk", []),
        ("RTK Kinematic (combined)", ["-p", "2", "-c", "-h"], "rtk", ["--enable-backward-smoothing"]),
        ("PPP Kinematic (EKF)", ["-p", "7", "-h"], "ppp", []),
        ("PPP Kinematic (FG)", ["-p", "7", "-h"], "ppp-fg", []),
    ]
    os.makedirs(OUT_DIR_RTKLIB, exist_ok=True)
    all_results = {}
    
    for ds_name, config in DATASETS.items():
        print(f"\n=== RTKLIB Compare: {ds_name} ===")
        all_results[ds_name] = {}
        for label, rtklib_args, gneiss_mode, gneiss_flags in comparisons:
            needs_base = any(p in rtklib_args for p in ["1", "2", "3", "4"])
            if needs_base and not config.get("base"): continue
            
            safe_name = ds_name.replace(" ", "_").replace("(", "").replace(")", "")
            safe_label = label.replace(" ", "_").replace("(", "").replace(")", "")
            r_out = os.path.join(OUT_DIR_RTKLIB, f"rtklib_{safe_name}_{safe_label}.pos")
            g_out = os.path.join(OUT_DIR_RTKLIB, f"gneiss_{safe_name}_{gneiss_mode}_{safe_label}.pos")
            
            if not eval_only:
                # RTKLIB
                cmd_r = [RTKLIB_BIN] + rtklib_args + ["-e", "-t"]
                if config.get("conf"): cmd_r += ["-k", config["conf"]]
                cmd_r += [config["rover"]]
                if needs_base: cmd_r += [config["base"]]
                cmd_r += [config["nav"]]
                if dry_run: print(f"[DRY] {' '.join(cmd_r)}")
                else:
                    with open(r_out, "w") as f: subprocess.run(cmd_r, stdout=f, stderr=subprocess.PIPE, text=True, timeout=600)
                
                # Gneiss
                cmd_g = [GNEISS_CLI, "process", "--mode", gneiss_mode, "--rover", config["rover"], "--output", g_out]
                if config.get("base"): cmd_g += ["--base", config["base"]]
                if config.get("nav"): cmd_g += ["--nav", config["nav"]]
                if config.get("gneiss_config"): cmd_g += ["--config", config["gneiss_config"]]
                if config.get("sp3"): cmd_g += ["--sp3", config["sp3"]]
                if config.get("clk"): cmd_g += ["--clk", config["clk"]]
                cmd_g += gneiss_flags
                if dry_run: print(f"[DRY] {' '.join(cmd_g)}")
                else: subprocess.run(cmd_g, capture_output=True, text=True, timeout=600)

            all_results[ds_name][label] = {
                "rtklib": evaluate(r_out, config["truth"], dry_run),
                "gneiss": evaluate(g_out, config["truth"], dry_run)
            }

    md = "# Gneiss vs RTKLIB (demo5) Comparison\n\n| Dataset | Mode | Engine | Hz 50th | Hz 95th | Vt 50th | Winner |\n|:-----|:-----|:-------|:--------|:--------|:--------|:-------|\n"
    for ds_name, modes in all_results.items():
        for label, data in modes.items():
            r = data["rtklib"]
            g = data["gneiss"]
            r_hz50 = r.get("hz_50", "N/A") if r and "error" not in r else "Failed"
            g_hz50 = g.get("hz_50", "N/A") if g and "error" not in g else "Failed"
            
            winner = ""
            if r and "error" not in r and g and "error" not in g:
                r_val, g_val = extract_float(r_hz50), extract_float(g_hz50)
                if g_val < r_val: winner = "**Gneiss**"
                elif r_val < g_val: winner = "RTKLIB"
                else: winner = "Tie"
            elif g and "error" not in g: winner = "**Gneiss** (RTKLIB missing)"
            
            md += f"| {ds_name} | {label} | RTKLIB | {r_hz50} | {r.get('hz_95','') if r else ''} | {r.get('vt_50','') if r else ''} | {winner} |\n"
            md += f"| {ds_name} | {label} | Gneiss | {g_hz50} | {g.get('hz_95','') if g else ''} | {g.get('vt_50','') if g else ''} | |\n"

    with open("COMPARISON.md", "w") as f: f.write(md)
    print("Wrote COMPARISON.md")

# Suite 3: Matrix
def run_matrix_suite(dry_run=False, eval_only=False):
    # For brevity, this merges the 18-grid logic
    print("Matrix suite triggered. Logic consolidated into COMPARISON.md for now.")
    # You can expand this with the full 18-grid logic if needed.
    run_rtklib_suite(dry_run, eval_only)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gneiss Unified Benchmark Orchestrator")
    parser.add_argument("--suite", choices=["gneiss", "rtklib", "matrix"], default="gneiss", help="Benchmark suite to run")
    parser.add_argument("--dataset", type=str, help="Filter by dataset name (substring match)")
    parser.add_argument("--dry-run", action="store_true", help="Print commands instead of executing")
    parser.add_argument("--eval-only", action="store_true", help="Skip processing, just evaluate existing .pos files")
    args = parser.parse_args()

    if args.dataset:
        DATASETS = {k: v for k, v in DATASETS.items() if args.dataset.lower() in k.lower()}

    if not args.dry_run and not args.eval_only:
        print("Building Gneiss...")
        subprocess.run(["cargo", "build", "--release", "--bin", "gneiss-cli"], check=True)

    if args.suite == "gneiss":
        run_gneiss_suite(args.dry_run, args.eval_only)
    elif args.suite == "rtklib":
        run_rtklib_suite(args.dry_run, args.eval_only)
    elif args.suite == "matrix":
        # Matrix delegates to run_full_matrix or an expanded version here
        run_matrix_suite(args.dry_run, args.eval_only)

    print("Done!")
