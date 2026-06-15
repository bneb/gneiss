#!/bin/bash
cargo build --release

# Ensure IMU fusion is enabled for Shinjuku benchmark
cat datasets/urbannav/tokyo/tokyo_config.json | sed 's/"enable_imu_fusion": false/"enable_imu_fusion": true/' > tokyo_tc_config.json

target/release/gneiss-cli process \
  --mode rtk-ins \
  --rover datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs \
  --base datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs \
  --nav datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav \
  --output out_shinjuku.pos \
  --config tokyo_tc_config.json

target/release/gneiss-cli eval \
  --solution out_shinjuku.pos \
  --truth datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv
