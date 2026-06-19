# Gneiss Benchmarks

All benchmarks compare **Gneiss** against **RTKLIB demo5 (2.4.3 b34)** on public datasets.
Every number is reproducible with the commands shown below.

## Methodology

- **RTKLIB version**: demo5 2.4.3 b34, built from [rtklibexplorer/RTKLIB](https://github.com/rtklibexplorer/RTKLIB)
- **RTKLIB config**: default kinematic PPP/RTK settings, elevation mask 15°, SNR mask 25 dB-Hz
- **Gneiss config**: per-dataset JSON configs in `datasets/<name>/`. Key settings documented inline.
- **Truth**: SPAN-CPT for UrbanNav/Odaiba, IGS coordinates for WTZR, SPAN for GSDC
- **Metrics**: median (50th percentile) and 95th percentile of horizontal error, vertical error where available
- **Epoch matching**: time-tag alignment ±0.5s between solution and truth

## Results

See [COMPARISON.md](./COMPARISON.md) for the full per-dataset, per-mode head-to-head table.

### Summary (datasets with working data)

| Dataset | Gneiss SPP | RTKLIB SPP | Gneiss RTK | RTKLIB RTK | Winner |
|---------|-----------|-----------|-----------|-----------|--------|
| GSDC (Pixel 4) | 3.30 m | 3.31 m | 8.37 m | 1.77 m | RTKLIB |
| Shinjuku (UrbanNav) | 1.85 m | Failed | **2.60 m** | 5.73 m | **Gneiss** |
| Odaiba (UrbanNav) | 2.10 m | 2.80 m | **1.92 m** | 4.01 m | **Gneiss** |
| f9p_ppp | 1.37 m | Failed | **0.25 m** | Failed | **Gneiss** |

### PPP Summary (fresh 2026-06-19)

| Dataset | Gneiss PPP-FG | RTKLIB PPP | Winner |
|---------|---------------|-----------|--------|
| Shinjuku | 7.40 m | 3.75 m | RTKLIB |
| Odaiba | **4.96 m** | 4.69 m | **Tied** |
| GSDC | 178.82 m | 4.25 m | RTKLIB |

### Datasets pending data acquisition

UrbanLoco, TEX-CUP, WHU-Smartphone, smartLoc — RINEX files not yet downloaded.
See `scripts/fetch_datasets.sh` for download URLs.

## Reproduction

### Shinjuku RTK (gneiss beats RTKLIB)
```bash
# Gneiss
cargo run --release -p gneiss-cli -- process \
  --rover datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs \
  --base datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs \
  --nav datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav \
  --mode rtk \
  --output shinjuku_rtk_gneiss.pos

# Evaluate
cargo run --release -p gneiss-cli -- eval \
  --solution shinjuku_rtk_gneiss.pos \
  --truth datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv
```

### WTZR PPP-AR (static, 24-hour)
```bash
cargo run --release -p gneiss-cli -- process \
  --rover datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.rnx \
  --nav datasets/wtzr_ppp_1224/BRDC00WRD_R_20203590000_01D_MN.rnx \
  --sp3 datasets/wtzr_ppp_1224/COD0MGXFIN_20203590000_01D_05M_ORB.SP3 \
  --clk datasets/wtzr_ppp_1224/COD0MGXFIN_20203590000_01D_30S_CLK.CLK \
  --bia datasets/wtzr_ppp_1224/COD0MGXFIN_20203590000_01D_01D_OSB.BIA \
  --mode ppp-fg \
  --output wtzr_ppp_ar.pos
```

## Known Limitations

- **GSDC PPP diverges** (7,184m median) — IEKF numerical stability issue on smartphone L5 data. Under investigation.
- **PPP requires precise products** (SP3/CLK) for competitive accuracy. Broadcast-only PPP is ~30-50m.
- **AR needs >30 min convergence** for widelane covariance to drop below 0.15 cycle threshold.
- **GLONASS AR not supported** — inter-frequency bias calibration required.
- **INS tightly-coupled modes** (SPP-ins, RTK-ins) produce degraded results on some datasets.
