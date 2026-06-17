with open("/tmp/debug.log") as f:
    for line in f:
        if "PPP Epoch" in line:
            print(line.strip())
            break
