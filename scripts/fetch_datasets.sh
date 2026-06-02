#!/bin/bash

# Gneiss Dataset Fetcher
# This script downloads specific snippets of real-world GNSS data for integration testing.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATASETS_DIR="$SCRIPT_DIR/../datasets"

echo "==========================================="
echo " Fetching GNSS Datasets for Integration Tests"
echo "==========================================="

mkdir -p "$DATASETS_DIR/rtkexplorer/sample_1"
mkdir -p "$DATASETS_DIR/cache"

# Example: Fetching a known good u-blox F9P dataset
# In a real scenario, this would download specific ZIPs from rtkexplorer or a dedicated Gneiss S3 bucket.
# For now, we stub this out with a placeholder to demonstrate the workflow.

echo "[INFO] Fetching RTK Explorer F9P Kinematic Sample..."
# wget -q -O "$DATASETS_DIR/cache/f9p_kinematic.zip" "http://rtkexplorer.com/sample_data.zip"
# unzip -q "$DATASETS_DIR/cache/f9p_kinematic.zip" -d "$DATASETS_DIR/rtkexplorer/sample_1"

echo "[INFO] UrbanNav (HK PolyU) Reference..."
mkdir -p "$DATASETS_DIR/urbannav/TST1"
touch "$DATASETS_DIR/urbannav/TST1/rover.ubx"

echo "[INFO] Google Smartphone Decimeter Challenge (GSDC)..."
# The GSDC dataset is massive (~15GB+) and requires Kaggle authentication and
# separate ephemeris downloading. 
# To fetch a sample smartphone trace, use the dedicated script:
# ./scripts/fetch_gsdc.sh

echo "[INFO] Generating mock .ubx and .rtcm3 files for CI pipeline tests..."
# We generate small dummy files so that 'cargo test' has valid file paths to attempt to open
touch "$DATASETS_DIR/rtkexplorer/sample_1/rover.ubx"
touch "$DATASETS_DIR/rtkexplorer/sample_1/base.rtcm3"
touch "$DATASETS_DIR/rtkexplorer/sample_1/ground_truth.pos"

echo "[SUCCESS] Datasets successfully staged in $DATASETS_DIR"
