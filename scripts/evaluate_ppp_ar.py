#!/usr/bin/env python3
import sys
import numpy as np
import argparse

def evaluate_ppp_ar(pos_file, true_x, true_y, true_z):
    true_pos = np.array([true_x, true_y, true_z])
    
    fixes = 0
    total = 0
    errors = []

    with open(pos_file, 'r') as f:
        for line in f:
            if line.startswith('%'): continue
            parts = line.split()
            # Expecting RTKLIB pos format or similar: 
            # Date Time X Y Z Q nsats ...
            if len(parts) < 6: continue
            
            try:
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                q = int(parts[5])
            except ValueError:
                continue

            total += 1
            pos = np.array([x, y, z])
            err = np.linalg.norm(pos - true_pos)
            errors.append(err)
            
            # Q=1 means FIXED, Q=2 means FLOAT, Q=5 means SPP
            if q == 1:
                fixes += 1

    if total == 0:
        print("No valid epochs found in file.")
        return

    fix_rate = (fixes / total) * 100
    errors_np = np.array(errors)
    
    # Calculate Final 3D Error (average of last 10 epochs)
    final_err = np.mean(errors_np[-10:]) if len(errors_np) >= 10 else np.mean(errors_np)
    
    # Calculate RMS 3D Error over the entire file
    rms_err = np.sqrt(np.mean(errors_np**2))

    print(f"Total Epochs: {total}")
    print(f"Fixed Epochs: {fixes}")
    print(f"Fix Rate: {fix_rate:.1f}%")
    print(f"Final 3D Error: {final_err:.3f} m")
    print(f"RMS 3D Error: {rms_err:.3f} m")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate PPP-AR performance from a .pos file.")
    parser.add_argument("pos_file", help="Path to the output .pos file.")
    parser.add_argument("--true-x", type=float, required=True, help="True ECEF X coordinate.")
    parser.add_argument("--true-y", type=float, required=True, help="True ECEF Y coordinate.")
    parser.add_argument("--true-z", type=float, required=True, help="True ECEF Z coordinate.")
    
    args = parser.parse_args()
    evaluate_ppp_ar(args.pos_file, args.true_x, args.true_y, args.true_z)
