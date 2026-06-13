import subprocess
import os

def run():
    print("Running PPP benchmark...")
    env = os.environ.copy()
    subprocess.run(["target/release/gneiss-cli", "process", "--mode", "ppp", 
                    "--rover", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs", 
                    "--nav", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav", 
                    "--config", "datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json", 
                    "--output", "scratch/ppp_only_test.pos"], env=env, check=True)
    print("Done")

if __name__ == "__main__":
    run()
