#!/bin/bash
set -e

# Download and compile RNXCMP (crx2rnx) from GSI

INSTALL_DIR="/Users/kevin/.cargo/bin"
RNXCMP_VERSION="4.2.0"
RNXCMP_URL="https://terras.gsi.go.jp/ja/crx2rnx/RNXCMP_${RNXCMP_VERSION}_src.tar.gz"

echo "[INFO] Downloading RNXCMP v${RNXCMP_VERSION}..."
curl -sL -O "$RNXCMP_URL"

echo "[INFO] Extracting..."
tar -xzf "RNXCMP_${RNXCMP_VERSION}_src.tar.gz"
cd "RNXCMP_${RNXCMP_VERSION}_src/source"

echo "[INFO] Compiling..."
gcc -O2 -o crx2rnx crx2rnx.c
gcc -O2 -o rnx2crx rnx2crx.c

echo "[INFO] Installing to $INSTALL_DIR (may require sudo)..."
cp crx2rnx "$INSTALL_DIR/"
cp rnx2crx "$INSTALL_DIR/"

echo "[INFO] Cleanup..."
cd ../..
rm -rf "RNXCMP_${RNXCMP_VERSION}_src" "RNXCMP_${RNXCMP_VERSION}_src.tar.gz"

echo "[INFO] crx2rnx successfully installed."
crx2rnx -h | head -n 5
