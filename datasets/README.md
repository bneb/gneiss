# Gneiss Datasets

This directory contains real-world GNSS/INS datasets used for integration testing, algorithm benchmarking, and Post-Processed Kinematic (PPK) validation.

Because raw `.ubx` and `.rtcm3` observation files can be extremely large, they are **not committed** to the git repository. 
Instead, we provide automated shell scripts in the `../scripts/` directory to fetch the exact test epochs required from public repositories.

## Directory Structure

*   `rtkexplorer/`: Sample datasets from the RTK Explorer community, primarily focusing on u-blox ZED-F9P dual-frequency kinematic and static scenarios. Excellent for validating the core LAMBDA integer resolution.
*   `urbannav/`: Complex, multi-sensor datasets recorded in harsh urban environments (high multipath, severe signal blockage) with SPAN-CPT ground truth. Excellent for stress-testing the Extended Kalman Filter (EKF).
*   `cache/`: Temporary directory for unzipping and processing large intermediate files.

## Recommended Data Sources

To rigorously test the engine, we recommend leveraging these public datasets:

### 1. High-End Commercial & Autonomous Benchmarks
*   **UrbanNav (HK PolyU):** Multi-sensor (LiDAR, IMU, GNSS) data from dense urban canyons. Includes NovAtel SPAN ground truth. [GitHub Repository](https://github.com/IPNL-POLYU/UrbanNavDataset)
*   **4Seasons:** Autonomous driving data across highway, tunnel, and parking environments using Septentrio mosaic-X5. [Website](https://www.4seasons-dataset.com/)
*   **LUCOOP (Leibniz Univ):** High-rate (10Hz) raw GNSS synchronized with multiple IMUs and LiDAR.

### 2. Smartphone & Low-Cost Hardware
*   **Google Smartphone Decimeter Challenge (GSDC):** Massive collection of smartphone GNSS data (Broadcom/Qualcomm) with NovAtel SPAN ground truth. Ideal for testing multipath mitigation and cycle slip detection. [Kaggle](https://www.kaggle.com/c/google-smartphone-decimeter-challenge)
*   **RTK Explorer (rtklibexplorer):** Curated datasets for low-cost receivers (u-blox F9P/M8T). [Blog/Resources](https://rtklibexplorer.wordpress.com/)

### 3. Geodetic & Atmospheric Reference
*   **NASA CDDIS:** Archive for space geodesy, daily/hourly RINEX data from global permanent stations. [Website](https://cddis.nasa.gov/)
*   **International GNSS Service (IGS):** Global network providing raw data and precise products (ephemerides, clocks). [Website](https://www.igs.org/)

## Fetching Data
To populate the standard community and academic datasets, run the main ingestion script:
```bash
./scripts/fetch_datasets.sh
```

**Note on Google Smartphone Decimeter Challenge (GSDC):**
The GSDC dataset requires a Kaggle account to download and does not provide satellite ephemeris natively. We have provided a dedicated script to fetch a sample smartphone trace from this competition. See the script for authentication requirements:
```bash
./scripts/fetch_gsdc.sh
```

## Adding New Datasets
When adding a new dataset for a specific integration test:
1. Create a dedicated folder (e.g., `datasets/my_custom_scenario/`).
2. Provide a `.json` or `.toml` metadata file detailing the known Base Station coordinates, antenna heights, and the expected ground-truth trajectory.
3. Update `scripts/fetch_datasets.sh` if the data is publicly hosted.
