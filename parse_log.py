import re

with open('out.log', 'r') as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "Smoother at k=5000." in line:
        print(line.strip())
        print(lines[i+1].strip())
        print(lines[i+2].strip())
        print(lines[i+3].strip())
        break
