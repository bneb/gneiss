#!/usr/bin/env bash
for dataset in f9p_ppp GSDC_Pixel_4 Shinjuku Odaiba; do
    echo "======================================"
    echo "DATASET: $dataset"
    echo "======================================"
    if [[ "$dataset" == "f9p_ppp" ]]; then
        truth="datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppk.pos"
    elif [[ "$dataset" == "GSDC_Pixel_4" ]]; then
        truth="datasets/gsdc/reference.csv"
    elif [[ "$dataset" == "Shinjuku" ]]; then
        truth="datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv"
    elif [[ "$dataset" == "Odaiba" ]]; then
        truth="datasets/urbannav/tokyo/Tokyo_Data/Odaiba/reference.csv"
    fi
    for mode in spp rtk RTK_Kinematic RTK_Kinematic_combined SPP; do
        for f in benchmarks/rtklib_comparison/*${dataset}*${mode}.pos benchmarks/rtklib_comparison/*${mode}*${dataset}*.pos; do
            if [[ -f "$f" ]]; then
                echo "--> $f"
                target/release/gneiss-cli eval --solution "$f" --truth "$truth" | grep -A 3 "| Horiz"
            fi
        done
    done
done
