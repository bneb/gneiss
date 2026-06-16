#!/bin/bash
cargo build --release

# Ensure IMU fusion is enabled for Shinjuku benchmark
cat datasets/urbannav/tokyo/tokyo_config.json | sed 's/"enable_imu_fusion": false/"enable_imu_fusion": true/' | sed 's/"mode": "Rtk"/"mode": "Rtk", "export_gnn_dataset_path": "shinjuku_gnn_dataset.csv"/' > tokyo_tc_config.json

echo "Running EKF Baseline..."
target/release/gneiss-cli process \
  --mode rtk-ins \
  --rover datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs \
  --base datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs \
  --nav datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav \
  --output out_shinjuku_ekf.pos \
  --config tokyo_tc_config.json

echo "Running FGO Backend..."
target/release/gneiss-cli process \
  --mode rtk-ins-fg \
  --rover datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs \
  --base datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs \
  --nav datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav \
  --output out_shinjuku_fgo.pos \
  --config tokyo_tc_config.json

echo "=== EKF Evaluation ==="
target/release/gneiss-cli eval \
  --solution out_shinjuku_ekf.pos \
  --truth datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv

echo "=== FGO Evaluation ==="
target/release/gneiss-cli eval \
  --solution out_shinjuku_fgo.pos \
  --truth datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv
