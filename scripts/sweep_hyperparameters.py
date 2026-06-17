import json
import subprocess
import itertools
import os
import re

CONFIG_FILE = "tokyo_tc_config.json"
OUTPUT_FILE = "sweep_results.csv"

# The EKF tuning parameters to sweep (Reduced search space for speed)
PARAM_GRID = {
    "sigma_ab": [1e-4, 1e-5, 1e-6],
    "sigma_gb": [1e-5, 1e-6, 1e-7]
}

def load_config():
    with open(CONFIG_FILE, 'r') as f:
        return json.load(f)

def save_config(config):
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f, indent=4)

def run_benchmark():
    # Run process
    cmd_process = [
        "target/release/gneiss-cli", "process",
        "--mode", "rtk-ins",
        "--rover", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
        "--base", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
        "--nav", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav",
        "--output", "out_shinjuku_ekf.pos",
        "--config", CONFIG_FILE
    ]
    subprocess.run(cmd_process, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Run eval
    cmd_eval = [
        "target/release/gneiss-cli", "eval",
        "--solution", "out_shinjuku_ekf.pos",
        "--truth", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv"
    ]
    res = subprocess.run(cmd_eval, capture_output=True, text=True)
    
    # Parse 95th percentile horizontal error
    for line in res.stdout.split('\n'):
        if line.startswith("| Horiz"):
            parts = line.split('|')
            if len(parts) >= 7:
                p95 = float(parts[6].strip().replace('m', '').strip())
                return p95
    return float('inf')

def main():
    keys, values = zip(*PARAM_GRID.items())
    combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    print(f"Total combinations to test: {len(combinations)}")
    base_config = load_config()

    with open(OUTPUT_FILE, 'w') as f:
        header = list(keys) + ["p95_horiz_error"]
        f.write(",".join(header) + "\n")

    best_p95 = float('inf')
    best_combo = None

    for i, combo in enumerate(combinations):
        print(f"[{i+1}/{len(combinations)}] Testing {combo}...")
        
        # Inject combo into config
        for k, v in combo.items():
            base_config["tuning"][k] = v
        
        save_config(base_config)
        
        p95 = run_benchmark()
        print(f"  -> 95th Percentile Horiz Error: {p95}m")

        with open(OUTPUT_FILE, 'a') as f:
            row = [str(combo[k]) for k in keys] + [str(p95)]
            f.write(",".join(row) + "\n")

        if p95 < best_p95:
            best_p95 = p95
            best_combo = combo

    print("\n===============================")
    print(f"BEST CONFIG (95th %: {best_p95}m):")
    for k, v in best_combo.items():
        print(f"  {k}: {v}")
    print("===============================")

if __name__ == "__main__":
    main()
