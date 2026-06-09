import json
with open('out.log', 'r') as f:
    for line in f:
        if "Smoother at k=5000" in line:
            break
