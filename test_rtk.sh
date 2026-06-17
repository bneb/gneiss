#!/bin/bash
set -e
cp datasets/urbannav/tokyo/tokyo_config.json tokyo_lc_config.json

target/release/gneiss-cli process \
    --mode rtk \
    --rover datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs \
    --base datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs \
    --nav datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav \
    --config tokyo_lc_config.json \
    --output out_shinjuku_rtk.pos

target/release/gneiss-cli eval \
    --solution out_shinjuku_rtk.pos \
    --truth datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv
