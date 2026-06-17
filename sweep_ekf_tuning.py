import json
import subprocess
import itertools
import os
import concurrent.futures

CONFIG_FILE = "tokyo_tc_config.json"
OUTPUT_FILE = "sweep_lever_arm_results.csv"

PARAM_GRID = {
    "sigma_v": [0.05, 0.1, 0.5],
    "sigma_phi": [0.01, 0.05, 0.1],
    "sigma_ab": [1e-6, 1e-5],
    "nhc_sigma": [0.5, 1.0, 2.0, 5.0]
}

def load_config():
    with open(CONFIG_FILE, 'r') as f:
        return json.load(f)

def save_config(config, filename):
    with open(filename, 'w') as f:
        json.dump(config, f, indent=4)

def run_benchmark(combo, i):
    # Create unique config and output files for parallel runs
    conf_filename = f"config_worker_{i}.json"
    out_filename = f"out_shinjuku_ekf_{i}.pos"

    base_config = load_config()
    base_config["enable_nhc"] = True
    
    # We found a reasonable GNSS lever arm
    base_config["imu_to_antenna_lever_arm"] = [1.0, 0.0, -0.5]
    
    # Apply tuning parameters
    if "tuning" not in base_config:
        base_config["tuning"] = {}
        
    base_config["tuning"]["sigma_v"] = combo["sigma_v"]
    base_config["tuning"]["sigma_phi"] = combo["sigma_phi"]
    base_config["tuning"]["sigma_ab"] = combo["sigma_ab"]
    base_config["tuning"]["nhc_sigma_lateral"] = combo["nhc_sigma"]
    base_config["tuning"]["nhc_sigma_vertical"] = combo["nhc_sigma"]
    
    save_config(base_config, conf_filename)

    cmd_process = [
        "target/release/gneiss-cli", "process",
        "--mode", "rtk-ins",
        "--rover", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
        "--base", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
        "--nav", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav",
        "--output", out_filename,
        "--config", conf_filename
    ]
    subprocess.run(cmd_process, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd_eval = [
        "target/release/gneiss-cli", "eval",
        "--solution", out_filename,
        "--truth", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv"
    ]
    res = subprocess.run(cmd_eval, capture_output=True, text=True)
    
    p95 = float('inf')
    for line in res.stdout.split('\n'):
        if line.startswith("| Horiz"):
            parts = line.split('|')
            if len(parts) >= 7:
                p95 = float(parts[6].strip().replace('m', '').strip())
                break
                
    # Cleanup temporary files
    if os.path.exists(conf_filename): os.remove(conf_filename)
    if os.path.exists(out_filename): os.remove(out_filename)
    
    return combo, p95

def main():
    keys, values = zip(*PARAM_GRID.items())
    combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    print(f"Total combinations to test: {len(combinations)}")
    
    with open(OUTPUT_FILE, 'w') as f:
        header = list(keys) + ["p95_horiz_error"]
        f.write(",".join(header) + "\n")

    best_p95 = float('inf')
    best_combo = None

    # Run in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(run_benchmark, combo, i): combo for i, combo in enumerate(combinations)}
        for future in concurrent.futures.as_completed(futures):
            combo = futures[future]
            try:
                result_combo, p95 = future.result()
                print(f"[{combo}] -> 95th Percentile Horiz Error: {p95}m")
                
                with open(OUTPUT_FILE, 'a') as f:
                    row = [str(result_combo[k]) for k in keys] + [str(p95)]
                    f.write(",".join(row) + "\n")

                if p95 < best_p95:
                    best_p95 = p95
                    best_combo = result_combo
            except Exception as exc:
                print(f"{combo} generated an exception: {exc}")

    print("\n===============================")
    print(f"BEST CONFIG (95th %: {best_p95}m):")
    for k, v in best_combo.items():
        print(f"  {k}: {v}")
    print("===============================")

if __name__ == "__main__":
    main()
