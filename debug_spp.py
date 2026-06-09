import subprocess
import sys

cmd = [
    "cargo", "run", "--release", "--bin", "gneiss-cli", "--",
    "process", "--mode", "spp", 
    "--rover", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/rover_ublox.obs", 
    "--output", "/tmp/shinjuku_spp_test.pos", 
    "--nav", "datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/base.nav", 
    "--config", "datasets/urbannav/tokyo/tokyo_config.json"
]

proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env={"RUST_LOG": "debug", "PATH": "/opt/homebrew/bin:/usr/bin:/bin"})
lines = []
for i, line in enumerate(proc.stdout):
    lines.append(line)
    if i >= 1000:
        proc.kill()
        break

with open("/tmp/spp_debug_py.log", "w") as f:
    f.writelines(lines)
