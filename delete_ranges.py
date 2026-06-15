with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    lines = f.readlines()

ranges_to_delete = [
    (35, 62),
    (65, 106),
    (352, 354),
    (371, 379),
    (381, 398),
    (400, 406),
    (408, 418),
    (707, 726), # includes the const DD_CROSS_CORRELATION_SCALE
]

lines_to_keep = []
for i, line in enumerate(lines):
    line_num = i + 1
    delete = False
    for start, end in ranges_to_delete:
        if start <= line_num <= end:
            delete = True
            break
    if not delete:
        lines_to_keep.append(line)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.writelines(lines_to_keep)

print("Deleted ranges.")
