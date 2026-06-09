#!/bin/bash
set -e

datasets=("Odaiba" "Shinjuku" "GSDC" "PPP")

for d in "${datasets[@]}"; do
  echo "========================================"
  echo "Evaluating $d"
  echo "========================================"
  
  if [ "$d" == "GSDC" ]; then
    truth="datasets/gsdc/reference.csv"
    
    echo "--- SPP ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_spp.pos --truth $truth | grep "Horiz"
    echo "RTKLIB:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/rtklib_GSDC_Pixel_4_SPP.pos --truth $truth | grep "Horiz"
    
    echo "--- RTK Forward ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_forward.pos --truth $truth | grep "Horiz"
    echo "RTKLIB:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/rtklib_GSDC_Pixel_4_RTK_Kinematic.pos --truth $truth | grep "Horiz"

    echo "--- RTK Smoothed ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_smoothed.pos --truth $truth | grep "Horiz"
    echo "RTKLIB:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/rtklib_GSDC_Pixel_4_RTK_Kinematic_combined.pos --truth $truth | grep "Horiz"
    
  elif [ "$d" == "PPP" ]; then
    truth="datasets/rtkexplorer/sample_1/ground_truth.pos"

    echo "--- SPP ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_spp.pos --truth $truth | grep "Horiz"
    
    echo "--- RTK Forward ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_forward.pos --truth $truth | grep "Horiz"

    echo "--- RTK Smoothed ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_smoothed.pos --truth $truth | grep "Horiz"

  else
    truth="datasets/urbannav/tokyo/Tokyo_Data/$d/reference.csv"

    echo "--- SPP ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_spp.pos --truth $truth | grep "Horiz"
    echo "RTKLIB:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/rtklib_${d}_u-blox_SPP.pos --truth $truth | grep "Horiz"
    
    echo "--- RTK Forward ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_forward.pos --truth $truth | grep "Horiz"
    echo "RTKLIB:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/rtklib_${d}_u-blox_RTK_Kinematic.pos --truth $truth | grep "Horiz"

    echo "--- RTK Smoothed ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_smoothed.pos --truth $truth | grep "Horiz"
    echo "RTKLIB:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/rtklib_${d}_u-blox_RTK_Kinematic_combined.pos --truth $truth | grep "Horiz"

    echo "--- RTK+INS Forward ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_ins_forward.pos --truth $truth | grep "Horiz"
    
    echo "--- RTK+INS Smoothed ---"
    echo "Gneiss:"
    cargo run --release --quiet --bin gneiss-cli -- eval --solution benchmarks/rtklib_comparison/gneiss_${d}_rtk_ins_smoothed.pos --truth $truth | grep "Horiz"
  fi
  echo ""
done
