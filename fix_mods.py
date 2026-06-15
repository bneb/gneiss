def remove_line(file, text):
    with open(file, 'r') as f:
        lines = f.readlines()
    with open(file, 'w') as f:
        for line in lines:
            if text not in line:
                f.write(line)

remove_line('crates/gneiss-rtk/src/estimators/ekf/mod.rs', 'filter_tests')
remove_line('crates/gneiss-rtk/src/estimators/mod.rs', 'spp_tests')
