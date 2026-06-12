import subprocess
cmd = ["cargo", "run", "--release", "--bin", "gneiss-cli", "--", "process", "--mode", "ppp", "--rover", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs", "--nav", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav", "--output", "out_ppp.pos", "--sp3", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_05M_ORB.SP3", "--clk", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_30S_CLK.CLK"]
proc = subprocess.run(cmd, env={"RUST_LOG": "trace", **__import__('os').environ}, capture_output=True, text=True)
with open("ppp_trace.log", "w") as f:
    f.write(proc.stderr)
