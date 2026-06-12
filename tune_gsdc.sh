#!/bin/bash
set -e

cargo build --release

truth="datasets/gsdc/reference.csv"

run_gsdc() {
  local slip=$1
  local cn0=$2
  local var=$3
  
  cat <<EOF > datasets/gsdc/test_config.json
{
    "enable_imu_fusion": false,
    "enable_backward_smoothing": true,
    "enable_nhc": false,
    "use_carrier_phase": true,
    "use_doppler_smoothing": true,
    "hatch_filter_length": 100,
    "min_elevation_deg": 15.0,
    "min_cn0_dbhz": $cn0,
    "disable_glonass": false,
    "disable_l2": false,
    "chi_square_pr_threshold": 30.0,
    "chi_square_cp_threshold": 1000000.0,
    "raim_pseudorange_outlier_m": 500.0,
    "nominal_snr_dbhz": 35.0,
    "imu_to_antenna_lever_arm": [0.0, 0.0, 0.0],
    "max_base_age_s": 31.0,
    "gnss_process_noise_var": $var,
    "initial_ambiguity_variance": 10000.0,
    "lambda_min_ratio": 3.0,
    "lambda_min_subset": 4,
    "ar_min_epoch_count": 5,
    "ar_min_lock": 3,
    "spp_consistency_threshold_m": 500.0,
    "max_ephemeris_age_s": 1800.0,
    "doppler_slip_threshold_cycles": $slip,
    "dynamics_model": "automotive",
    "max_reject_count": 10,
    "base_position": [-2689639.5060, -4290438.6360, 3865050.9560]
}
EOF

  target/release/gneiss-cli process --mode rtk --enable-backward-smoothing --rover datasets/gsdc/Pixel4_GnssLog.20o --base datasets/gsdc/p2221350.20o --nav datasets/gsdc/Pixel4_GnssLog.nav --output gsdc_test.pos --config datasets/gsdc/test_config.json > /dev/null 2>&1
  
  echo -n "Slip: $slip, CN0: $cn0, Var: $var | "
  target/release/gneiss-cli eval --solution gsdc_test.pos --truth $truth | grep "Horiz"
}

for slip in 15.0 50.0 100.0; do
  for cn0 in 15.0 20.0 25.0; do
    for var in 10.0 100.0 500.0; do
      run_gsdc $slip $cn0 $var
    done
  done
done
