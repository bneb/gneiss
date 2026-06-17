import sys

time_target = "2020 12 24 22 13"
in_time = False

with open("datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs") as f:
    for line in f:
        if line.startswith(">"):
            if time_target in line:
                in_time = True
                print("\n" + line.strip())
            elif in_time:
                # Stop after a few epochs
                print("\n" + line.strip())
                if "2020 12 24 22 13 10" in line:
                    break
        elif in_time:
            if line.startswith("G04"):
                print(line.strip())
