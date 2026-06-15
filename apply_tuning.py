import json

best_f9p = {
    "snr_a": 2.1030935228290484,
    "snr_b": 396.44303052639924,
    "phase_outlier_ratio_thresh": 7.191519122560652,
    "doppler_outlier_ratio_mult": 4.8235811218000535,
    "loosely_coupled_mahalanobis_sq": 175.83997002496068
}
best_odaiba = {
    "snr_a": 0.15719592697508833,
    "snr_b": 144.48196248109252,
    "phase_outlier_ratio_thresh": 14.584475964582582,
    "doppler_outlier_ratio_mult": 5.809658542941578,
    "loosely_coupled_mahalanobis_sq": 880.1188109642901
}
best_gsdc = {
    "snr_a": 2.913449293357578,
    "snr_b": 408.75325641047357,
    "phase_outlier_ratio_thresh": 5.1815692422278,
    "doppler_outlier_ratio_mult": 9.656904482353903,
    "loosely_coupled_mahalanobis_sq": 736.4084756538158
}

configs = [
    ("datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json", best_f9p),
    ("datasets/urbannav/tokyo/tokyo_config.json", best_odaiba),
    ("datasets/gsdc/gsdc_config.json", best_gsdc)
]

for path, tuning in configs:
    with open(path, 'r') as f:
        config = json.load(f)
    config["tuning"] = tuning
    with open(path, 'w') as f:
        json.dump(config, f, indent=4)
        
print("Updated configurations.")
