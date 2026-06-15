#!/usr/bin/env python3
import json
import random
import subprocess
import os
import re

def random_params():
    return {
        "snr_a": random.uniform(0.001, 3.0),
        "snr_b": random.uniform(10.0, 500.0),
        "phase_outlier_ratio_thresh": random.uniform(2.0, 15.0),
        "doppler_outlier_ratio_mult": random.uniform(1.0, 10.0),
        "loosely_coupled_mahalanobis_sq": random.uniform(10.0, 1000.0),
    }

def run_cmd(cmd):
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    return result.stdout, result.stderr

def extract_h50_v50(eval_output):
    h50 = float('inf')
    v50 = float('inf')
    for line in eval_output.split('\n'):
        if line.startswith('| Horiz'):
            parts = line.split('|')
            if len(parts) > 3:
                h50 = float(parts[2].strip().replace('m', '').strip())
        if line.startswith('| Vert'):
            parts = line.split('|')
            if len(parts) > 3:
                v50 = float(parts[2].strip().replace('m', '').strip())
    return h50, v50

def sweep_dataset(name, base_config_path, process_cmd_template, eval_cmd_template, iterations=30):
    print(f"\n======================================")
    print(f"Sweeping {name}")
    print(f"======================================")
    
    with open(base_config_path, 'r') as f:
        base_config = json.load(f)
        
    best_score = float('inf')
    best_params = None
    
    for i in range(iterations):
        params = random_params()
        config = base_config.copy()
        
        # Inject tuning parameters
        if "tuning" not in config:
            config["tuning"] = {}
            
        for k, v in params.items():
            config["tuning"][k] = v
            
        temp_config_path = f"sweep_temp_{name}.json"
        with open(temp_config_path, 'w') as f:
            json.dump(config, f, indent=4)
            
        process_cmd = process_cmd_template.replace("{config}", temp_config_path)
        print(f"[{i+1}/{iterations}] Running process...")
        run_cmd(process_cmd)
        
        eval_cmd = eval_cmd_template
        stdout, _ = run_cmd(eval_cmd)
        
        h50, v50 = extract_h50_v50(stdout)
        score = h50 + v50 * 0.5  # Weight H50 more, but factor in V50
        
        print(f"[{i+1}/{iterations}] Score: {score:.4f} (H50={h50:.4f}, V50={v50:.4f}) | Params: a={params['snr_a']:.3f}, b={params['snr_b']:.1f}, phase_thresh={params['phase_outlier_ratio_thresh']:.1f}")
        
        if score < best_score:
            best_score = score
            best_params = params
            print(f"  *** NEW BEST! ***")
            
        # Clean up temp files
        if os.path.exists(temp_config_path):
            os.remove(temp_config_path)
            
    print(f"\n[DONE] Best {name} Score: {best_score:.4f}")
    print(f"Best Params: {json.dumps(best_params, indent=4)}")
    return best_params

if __name__ == '__main__':
    # Build release first
    print("Building release...")
    run_cmd("cargo build --release")
    
    # 1. F9P PPP
    best_f9p = sweep_dataset(
        name="F9P_PPP",
        base_config_path="datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json",
        process_cmd_template="target/release/gneiss-cli process --mode rtk --rover datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs --base datasets/rtkexplorer/sample_1/f9p_ppp_1224/tmg23590.20o --nav datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav --output sweep_f9p.pos --config {config}",
        eval_cmd_template="target/release/gneiss-cli eval --solution sweep_f9p.pos --truth datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppk.pos",
        iterations=30
    )
    
    # 2. Odaiba
    best_odaiba = sweep_dataset(
        name="Odaiba_UrbanNav",
        base_config_path="datasets/urbannav/tokyo/tokyo_config.json",
        process_cmd_template="target/release/gneiss-cli process --mode rtk-ins --rover datasets/urbannav/tokyo/Tokyo_Data/Odaiba/rover_ublox.obs --base datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base_trimble.obs --nav datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base.nav --output sweep_odaiba.pos --config {config}",
        eval_cmd_template="target/release/gneiss-cli eval --solution sweep_odaiba.pos --truth datasets/urbannav/tokyo/Tokyo_Data/Odaiba/reference.csv",
        iterations=30
    )
    
    # 3. GSDC Pixel 4
    best_gsdc = sweep_dataset(
        name="GSDC_Pixel4",
        base_config_path="datasets/gsdc/gsdc_config.json",
        process_cmd_template="target/release/gneiss-cli process --mode rtk --rover datasets/gsdc/Pixel4_GnssLog.20o --base datasets/gsdc/p2221350.20o --nav datasets/gsdc/Pixel4_GnssLog.nav --output sweep_gsdc.pos --config {config}",
        eval_cmd_template="target/release/gneiss-cli eval --solution sweep_gsdc.pos --truth datasets/gsdc/reference.csv",
        iterations=30
    )
    
    print("\n\n====== FINAL RESULTS ======")
    print("F9P_PPP:", json.dumps(best_f9p, indent=4))
    print("Odaiba_UrbanNav:", json.dumps(best_odaiba, indent=4))
    print("GSDC_Pixel4:", json.dumps(best_gsdc, indent=4))
    
    with open('sweep_results_best.json', 'w') as f:
        json.dump({
            "f9p": best_f9p,
            "odaiba": best_odaiba,
            "gsdc": best_gsdc
        }, f, indent=4)
