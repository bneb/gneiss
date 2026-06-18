#!/usr/bin/env bash
# Reproduce the Shinjuku RTK benchmark: gneiss vs RTKLIB demo5
# Requires: gneiss-cli built in release mode, dataset downloaded
# Dataset: UrbanNav Tokyo Shinjuku (https://github.com/IPNL-POLYU/UrbanNavDataset)
# RTKLIB demo5: https://github.com/rtklibexplorer/RTKLIB

set -euo pipefail

ROVER="datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs"
BASE="datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs"
NAV="datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav"
TRUTH="datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv"
OUTDIR="/tmp/shinjuku_rtk_benchmark"

mkdir -p "$OUTDIR"

echo "=== Gneiss Shinjuku RTK ==="
echo "Rover: $ROVER"
echo "Base:  $BASE"

# Build if needed
cargo build --release -p gneiss-cli 2>/dev/null

# Run gneiss
echo ""
echo "Running gneiss RTK..."
./target/release/gneiss-cli process \
    --rover "$ROVER" \
    --base "$BASE" \
    --nav "$NAV" \
    --mode rtk \
    --output "$OUTDIR/gneiss_rtk.pos"

# Evaluate
echo ""
echo "=== Gneiss RTK Results ==="
./target/release/gneiss-cli eval \
    --solution "$OUTDIR/gneiss_rtk.pos" \
    --truth "$TRUTH"

# RTKLIB comparison (if available)
RTKLIB_BIN="${RTKLIB_BIN:-rtkrcv}"
if command -v "$RTKLIB_BIN" &>/dev/null; then
    echo ""
    echo "=== RTKLIB Comparison ==="
    echo "RTKLIB binary: $RTKLIB_BIN"
    echo "To complete the comparison, run RTKLIB on the same rover/base/nav files"
    echo "and evaluate with: gneiss-cli eval --solution rtklib_output.pos --truth $TRUTH"
else
    echo ""
    echo "RTKLIB not found in PATH."
    echo "Install RTKLIB demo5 from https://github.com/rtklibexplorer/RTKLIB"
    echo "Expected RTKLIB result: 2.21m horizontal median (from COMPARISON.md)"
fi

echo ""
echo "Expected gneiss result: ~1.34m horizontal median (from COMPARISON.md)"
echo "Output: $OUTDIR/gneiss_rtk.pos"
