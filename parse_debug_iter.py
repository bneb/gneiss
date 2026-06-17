with open("/tmp/debug.log") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "vel=[0.42, -0.29, 0.11]" in line:
        for j in range(max(0, i-50), min(len(lines), i+1)):
            if "GN Iter" in lines[j]:
                print(lines[j].strip())
        break
