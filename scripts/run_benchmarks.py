#!/usr/bin/env python3
import os
import subprocess
import json

MODES = [
    "spp",
    "spp-ins",
    "spp-ins-loosely-coupled",
    "rtk",
    "rtk-ins",
    "rtk-ins-loosely-coupled",
    "ppp",
    "ppp-ins",
]

DATASETS = {
    "GSDC (Pixel 4)": {
        "rover": "datasets/gsdc/Pixel4_GnssLog.20o",
        "base": "datasets/gsdc/p2221350.20o",
        "nav": "datasets/gsdc/rover.nav",
        "truth": "datasets/gsdc/reference.csv",
        "extra_args": ["--config", "datasets/gsdc/gsdc_config.json"],
    },
    "Shinjuku (UrbanNav)": {
        "rover": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
        "base": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
        "nav": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav",
        "truth": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv",
        "extra_args": ["--enable-backward-smoothing"],
    },
    "PPP (f9p_ppp)": {
        "rover": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs",
        "base": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/tmg23590.20o",
        "nav": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav",
        "truth": "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppk.pos",
        "extra_args": ["--config", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json"],
    }
}

CLI_BIN = "target/release/gneiss-cli"

def build_cli():
    print("Building gneiss-cli...")
    subprocess.run(["cargo", "build", "--release", "--bin", "gneiss-cli"], check=True)

def parse_eval_output(stdout, stderr):
    # gneiss-cli eval logs its output to stderr via tracing
    # Looking for lines like:
    # Median Horizontal Error: 0.028 m
    # 95th Percentile Horizontal Error: 0.640 m
    # Median Vertical Error: 0.012 m
    
    # In case it's in stdout or stderr
    text = stdout + "\n" + stderr
    metrics = {
        "hz_median": "N/A",
        "hz_95": "N/A",
        "vt_median": "N/A"
    }
    
    if "No matching epochs found" in text:
        return "Mismatch"
    
    for line in text.splitlines():
        if line.startswith("| Horiz   |"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) > 4:
                metrics["hz_median"] = parts[2]
                metrics["hz_95"] = parts[3]
        elif line.startswith("| Vert    |"):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) > 4:
                metrics["vt_median"] = parts[2]
                
    if metrics["hz_median"] == "N/A":
        # Check if there's any generic error
        if "ERROR" in text or "WARN" in text:
            return "Failed/Error"
            
    return metrics

def run_benchmark():
    results = {ds: {} for ds in DATASETS}
    
    for ds_name, config in DATASETS.items():
        print(f"\n=============================")
        print(f"Running Dataset: {ds_name}")
        print(f"=============================")
        
        for mode in MODES:
            print(f"  Mode: {mode}")
            
            # 1. Process
            out_file = f"benchmarks/{ds_name.replace(' ', '_').replace('(', '').replace(')', '')}_{mode}.pos"
            
            cmd = [CLI_BIN, "process", "--mode", mode, "--rover", config["rover"], "--output", out_file]
            if config["base"]:
                cmd.extend(["--base", config["base"]])
            if config["nav"]:
                cmd.extend(["--nav", config["nav"]])
            cmd.extend(config["extra_args"])
            
            try:
                # Run process
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                if proc.returncode != 0:
                    print(f"    -> Process Failed")
                    results[ds_name][mode] = "Process Failed"
                    continue
            except subprocess.TimeoutExpired:
                print(f"    -> Timeout")
                results[ds_name][mode] = "Timeout"
                continue
                
            # 2. Evaluate
            if not os.path.exists(config["truth"]):
                print(f"    -> Missing truth file: {config['truth']}")
                results[ds_name][mode] = "No Truth"
                continue
                
            eval_cmd = [CLI_BIN, "eval", "--solution", out_file, "--truth", config["truth"]]
            eval_proc = subprocess.run(eval_cmd, capture_output=True, text=True)
            
            metrics = parse_eval_output(eval_proc.stdout, eval_proc.stderr)
            if isinstance(metrics, dict):
                print(f"    -> {metrics['hz_median']} Hz, {metrics['vt_median']} Vt")
            else:
                print(f"    -> {metrics}")
            
            results[ds_name][mode] = metrics
            
    return results

def generate_markdown(results):
    md = "# Gneiss Comprehensive Benchmarks\n\n"
    md += "This document empirically maps the performance of Gneiss across varying modes and datasets.\n\n"
    
    for ds_name, modes_data in results.items():
        md += f"## {ds_name}\n\n"
        md += "| Mode | Median Horizontal | 95% Horizontal | Median Vertical |\n"
        md += "| :--- | :--- | :--- | :--- |\n"
        
        for mode in MODES:
            res = modes_data.get(mode, "N/A")
            if isinstance(res, dict):
                md += f"| `{mode}` | {res['hz_median']} | {res['hz_95']} | {res['vt_median']} |\n"
            else:
                md += f"| `{mode}` | {res} | {res} | {res} |\n"
        
        md += "\n"
        
    with open("BENCHMARKS.md", "w") as f:
        f.write(md)
    print("Wrote BENCHMARKS.md")

if __name__ == "__main__":
    os.makedirs("benchmarks", exist_ok=True)
    build_cli()
    results = run_benchmark()
    generate_markdown(results)
