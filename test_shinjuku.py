import subprocess
import os

ds_config = {
    "rover": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs",
    "base": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base_trimble.obs",
    "nav": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav",
    "truth": "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv",
    "conf": "datasets/urbannav/tokyo/tokyo_config.json"
}

cmd = [
    "target/release/gneiss-cli", "process",
    "--mode", "rtk",
    "--rover", ds_config["rover"],
    "--base", ds_config["base"],
    "--nav", ds_config["nav"],
    "--config", ds_config["conf"],
    "--enable-backward-smoothing",
    "--output", "shinjuku_test.pos"
]

print("Running command...")
proc = subprocess.run(cmd, capture_output=True, text=True)
if proc.returncode != 0:
    print("FAILED!", proc.returncode)
    print(proc.stderr)
else:
    print("Success! Evaluating...")
    eval_cmd = ["target/release/gneiss-cli", "eval", "--solution", "shinjuku_test.pos", "--truth", ds_config["truth"]]
    eval_proc = subprocess.run(eval_cmd, capture_output=True, text=True)
    print(eval_proc.stdout)
