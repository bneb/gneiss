import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

content = content.replace("""    let lower_bound = match meas_type {
        0 => max_innovation * max_innovation * 25.0, // Phase
        1 | 2 => max_innovation * max_innovation * 225.0, // Doppler
        3 => max_innovation * max_innovation * 625.0, // Pseudorange
        _ => max_innovation * max_innovation * 25.0,
    };""", """    let lower_bound = match meas_type {
        1 | 2 => max_innovation * max_innovation * 225.0, // Doppler
        3 => max_innovation * max_innovation * 625.0, // Pseudorange
        _ => max_innovation * max_innovation * 25.0, // Phase and default
    };""")

content = content.replace("""    let upper_bound = match meas_type {
        0 => max_innovation * max_innovation * 625.0, // Phase
        1 | 2 => max_innovation * max_innovation * 2500.0, // Doppler
        3 => max_innovation * max_innovation * 10000.0, // Pseudorange
        _ => max_innovation * max_innovation * 625.0,
    };""", """    let upper_bound = match meas_type {
        1 | 2 => max_innovation * max_innovation * 2500.0, // Doppler
        3 => max_innovation * max_innovation * 10000.0, // Pseudorange
        _ => max_innovation * max_innovation * 625.0, // Phase and default
    };""")

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
