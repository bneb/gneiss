#!/bin/bash

# Gneiss: Google Smartphone Decimeter Challenge (GSDC) Fetcher
# 
# The GSDC dataset is hosted on Kaggle and is massive (~15GB+).
# It provides highly degraded smartphone GNSS logs (Broadcom/Qualcomm)
# which are used for testing our Adaptive EKF and MAD RAIM.
#
# CRITICAL REQUIREMENTS:
# 1. Kaggle Authentication: You must have the `kaggle` CLI installed and authenticated.
#    - pip install kaggle
#    - Place your kaggle.json token in ~/.kaggle/
# 2. Rules Acceptance: You MUST manually accept the competition rules on the Kaggle website
#    before the API will authorize any downloads.
#    https://www.kaggle.com/c/google-smartphone-decimeter-challenge/rules
# 3. Ephemeris/CORS Missing: The GSDC dataset ONLY provides rover observation files.
#    It does NOT provide satellite ephemeris (.nav) or base station (.obs) files.
#    To process these files, you must independently fetch broadcast ephemeris (BRDC) 
#    from NASA CDDIS/IGS for the specific date of the drive.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATASETS_DIR="$SCRIPT_DIR/../datasets/gsdc"

echo "=========================================================="
echo " Fetching Google Smartphone Decimeter Challenge (GSDC)"
echo "=========================================================="

if ! command -v kaggle &> /dev/null; then
    echo "[ERROR] kaggle CLI is not installed or not in PATH."
    echo "Please run: pip install kaggle"
    exit 1
fi

mkdir -p "$DATASETS_DIR"

# Example: Fetching a single 20-minute drive from 2020-05-14 to prevent filling the hard drive
DRIVE="train/2020-05-14-US-MTV-1/Pixel4"

echo "[INFO] Downloading Rover Observation for Drive: $DRIVE..."
if kaggle competitions download -c google-smartphone-decimeter-challenge -f "$DRIVE/supplemental/Pixel4_GnssLog.20o" -p "$DATASETS_DIR"; then
    unzip -o "$DATASETS_DIR/Pixel4_GnssLog.20o.zip" -d "$DATASETS_DIR"
    rm "$DATASETS_DIR/Pixel4_GnssLog.20o.zip"
else
    echo "[ERROR] Failed to download. Did you accept the competition rules on Kaggle?"
    exit 1
fi

echo "[INFO] Downloading Ground Truth for Drive: $DRIVE..."
kaggle competitions download -c google-smartphone-decimeter-challenge -f "$DRIVE/ground_truth.csv" -p "$DATASETS_DIR"

echo ""
echo "[SUCCESS] GSDC sample drive staged in $DATASETS_DIR"
echo ""
echo "----------------------------------------------------------"
echo " ATTENTION: EPHEMERIS REQUIRED"
echo "----------------------------------------------------------"
echo "The GSDC dataset does not include satellite navigation files."
echo "Before processing 'Pixel4_GnssLog.20o' through Gneiss, you must:"
echo "1. Fetch the daily broadcast ephemeris for 2020-05-14."
echo "   (e.g., from NASA CDDIS FTP: https://cddis.nasa.gov/archive/gnss/data/daily/2020/brdc/)"
echo "2. Place the ephemeris file in the $DATASETS_DIR folder as 'rover.nav'."
echo "----------------------------------------------------------"
