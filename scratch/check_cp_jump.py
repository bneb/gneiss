import csv
with open('data/gsdc/Pixel4/2020-05-14-US-MTV-1/GNSS_Log.txt') as f:
    for line in f:
        if line.startswith('Raw'):
            parts = line.split(',')
            if len(parts) > 20:
                print(parts)
                break
