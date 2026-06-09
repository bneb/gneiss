#!/bin/bash
set -e

echo "Building Gneiss..."
cargo build --release

datasets=("Odaiba" "Shinjuku" "GSDC" "PPP")

for d in "${datasets[@]}"; do
  echo "Evaluating $d..."
  
  if [ "$d" == "GSDC" ]; then
    # SPP
    target/release/gneiss-cli process --mode spp --rover datasets/gsdc/Pixel4_GnssLog.20o --base datasets/gsdc/p2221350.20o --nav datasets/gsdc/Pixel4_GnssLog.nav --output benchmarks/rtklib_comparison/gneiss_${d}_spp.pos --config datasets/gsdc/gsdc_config.json
    
    # RTK Forward
    target/release/gneiss-cli process --mode rtk --rover datasets/gsdc/Pixel4_GnssLog.20o --base datasets/gsdc/p2221350.20o --nav datasets/gsdc/Pixel4_GnssLog.nav --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_forward.pos --config datasets/gsdc/gsdc_config.json

    # RTK Smoothed
    target/release/gneiss-cli process --mode rtk --enable-backward-smoothing --rover datasets/gsdc/Pixel4_GnssLog.20o --base datasets/gsdc/p2221350.20o --nav datasets/gsdc/Pixel4_GnssLog.nav --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_smoothed.pos --config datasets/gsdc/gsdc_config.json
    
  elif [ "$d" == "PPP" ]; then
    # SPP
    target/release/gneiss-cli process --mode spp --rover datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs --base datasets/rtkexplorer/sample_1/base.rtcm3 --nav datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav --output benchmarks/rtklib_comparison/gneiss_${d}_spp.pos --config datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json
    
    # RTK Forward
    target/release/gneiss-cli process --mode rtk --rover datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs --base datasets/rtkexplorer/sample_1/base.rtcm3 --nav datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_forward.pos --config datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json

    # RTK Smoothed
    target/release/gneiss-cli process --mode rtk --enable-backward-smoothing --rover datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs --base datasets/rtkexplorer/sample_1/base.rtcm3 --nav datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_smoothed.pos --config datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json
    
  else
    # SPP
    target/release/gneiss-cli process --mode spp --rover datasets/urbannav/tokyo/Tokyo_Data/$d/rover_ublox.obs --output benchmarks/rtklib_comparison/gneiss_${d}_spp.pos --base datasets/urbannav/tokyo/Tokyo_Data/$d/base_trimble.obs --nav datasets/urbannav/tokyo/Tokyo_Data/$d/base.nav --config datasets/urbannav/tokyo/tokyo_config.json

    # RTK Forward
    target/release/gneiss-cli process --mode rtk --rover datasets/urbannav/tokyo/Tokyo_Data/$d/rover_ublox.obs --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_forward.pos --base datasets/urbannav/tokyo/Tokyo_Data/$d/base_trimble.obs --nav datasets/urbannav/tokyo/Tokyo_Data/$d/base.nav --config datasets/urbannav/tokyo/tokyo_config.json

    # RTK Smoothed
    target/release/gneiss-cli process --mode rtk --enable-backward-smoothing --rover datasets/urbannav/tokyo/Tokyo_Data/$d/rover_ublox.obs --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_smoothed.pos --base datasets/urbannav/tokyo/Tokyo_Data/$d/base_trimble.obs --nav datasets/urbannav/tokyo/Tokyo_Data/$d/base.nav --config datasets/urbannav/tokyo/tokyo_config.json

    # RTK+INS Forward
    target/release/gneiss-cli process --mode rtk-ins --rover datasets/urbannav/tokyo/Tokyo_Data/$d/rover_ublox.obs --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_ins_forward.pos --base datasets/urbannav/tokyo/Tokyo_Data/$d/base_trimble.obs --nav datasets/urbannav/tokyo/Tokyo_Data/$d/base.nav --config datasets/urbannav/tokyo/tokyo_config.json

    # RTK+INS Smoothed
    target/release/gneiss-cli process --mode rtk-ins --enable-backward-smoothing --rover datasets/urbannav/tokyo/Tokyo_Data/$d/rover_ublox.obs --output benchmarks/rtklib_comparison/gneiss_${d}_rtk_ins_smoothed.pos --base datasets/urbannav/tokyo/Tokyo_Data/$d/base_trimble.obs --nav datasets/urbannav/tokyo/Tokyo_Data/$d/base.nav --config datasets/urbannav/tokyo/tokyo_config.json
  fi
done
