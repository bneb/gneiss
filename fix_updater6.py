import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

# test_compute_scalar_thresholds should set pr_abs_thresh = 50.0, dop_abs_thresh = 60.0, phase_outlier_ratio_thresh=2.0
content = content.replace(
    "let tuning = EkfTuningConfig::default();",
    "let mut tuning = EkfTuningConfig::default(); tuning.pr_abs_thresh = 50.0; tuning.dop_abs_thresh = 60.0; tuning.phase_outlier_ratio_thresh = 2.0; tuning.doppler_outlier_ratio_mult = 3.0;"
)
# We also need to update the assertions!
# PR (meas=0): thresh = max_inn, abs_thresh = 50.0
content = content.replace("assert_eq!(abs_thresh_pr, 40.0);", "assert_eq!(abs_thresh_pr, 50.0);")
content = content.replace("assert_eq!(abs_thresh_cp, 1000000.0);", "assert_eq!(abs_thresh_cp, 1000000.0);")
content = content.replace("assert_eq!(abs_thresh_doppler, tuning.dop_abs_thresh);", "assert_eq!(abs_thresh_doppler, 60.0);")

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
