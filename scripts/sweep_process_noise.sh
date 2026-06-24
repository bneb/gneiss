#!/bin/bash
# Phase 4: Automated process noise calibration sweep
# Sweeps process_noise_cb, process_noise_cd, process_noise_amb_float
# and evaluates PPP accuracy against IGS ground truth.
# Usage: ./scripts/sweep_process_noise.sh

set -e
BINARY=./target/release/gneiss-cli
ROVER=datasets/igs/suth3350.19o
NAV=datasets/igs/brdc3350.19n
SP3=datasets/igs/cod20820.sp3
CLK=datasets/igs/gfz20820.clk
TRUTH_X=5041274.8417
TRUTH_Y=1916054.1204
TRUTH_Z=-3397075.9716
EPOCHS=500
RESULTS=/tmp/sweep_results.csv

echo "cb,cd,amb_float,hz50,min" > $RESULTS

for cb in 0.1 1.0 10.0 100.0; do
for cd in 1.0 10.0 100.0 1000.0; do
for amb in 1e-9 1e-8 1e-7 1e-6 1e-5; do
    CONFIG=$(mktemp)
    echo "{\"mode\":\"Ppp\",\"process_noise_cb\":$cb,\"process_noise_cd\":$cd,\"process_noise_amb_float\":$amb,\"min_snr_dbhz\":0,\"enable_backward_smoothing\":false}" > $CONFIG
    OUT=/tmp/sweep_${cb}_${cd}_${amb}.pos

    $BINARY process --mode ppp --rover $ROVER --nav $NAV --output $OUT \
        --sp3 $SP3 --clk $CLK --config $CONFIG --max-epochs $EPOCHS 2>/dev/null

    # Evaluate
    python3 -c "
import math
truth = ($TRUTH_X, $TRUTH_Y, $TRUTH_Z)
errors = []
with open('$OUT') as f:
    for line in f:
        if line.startswith('%'): continue
        p = line.split()
        if len(p) >= 5:
            dx = float(p[2]) - truth[0]
            dy = float(p[3]) - truth[1]
            dz = float(p[4]) - truth[2]
            errors.append(math.sqrt(dx*dx + dy*dy))
errors.sort()
n = len(errors)
if n > 0:
    print(f'{errors[n//2]:.4f},{min(errors):.4f}')
" >> $RESULTS

    rm -f $CONFIG $OUT
done
done
done

echo "Sweep complete. Results in $RESULTS"
sort -t, -k4 -n $RESULTS | head -10