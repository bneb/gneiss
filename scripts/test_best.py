import json
import subprocess

with open("datasets/urbannav/tokyo/tokyo_config.json", "r") as f:
    config = json.load(f)

config["enable_imu_fusion"] = True
config["mode"] = "Rtk"
config["export_gnn_dataset_path"] = "shinjuku_gnn_dataset.csv"
config["tuning"]["huber_threshold_tightly"] = 4.0

with open("tokyo_tc_config.json", "w") as f:
    json.dump(config, f, indent=4)

cmd_process = [
    "target/release/gneiss-cli", "process",
    "--mode", "rtk-ins",
    "--rover", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
    "--base", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
    "--nav", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav",
    "--output", "out_shinjuku_ekf.pos",
    "--config", "tokyo_tc_config.json"
]
subprocess.run(cmd_process, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

cmd_eval = [
    "target/release/gneiss-cli", "eval",
    "--solution", "out_shinjuku_ekf.pos",
    "--truth", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv"
]
subprocess.run(cmd_eval)
