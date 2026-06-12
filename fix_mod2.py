import re
with open('crates/gneiss-rtk/src/factor_graph/mod.rs', 'r') as f:
    text = f.read()

target = """            let err = 0.5 * (res.transpose() * info.clone() * res.clone())[0];
            current_error += err;
            if err > 1_000_000.0 {
                println!("HUGE ERROR in Factor {}: {}", i, err);
            }"""

replacement = """            let err = 0.5 * (res.transpose() * info.clone() * res.clone())[0];
            current_error += err;
            if _iter == 0 && i < 5 {
                println!("Factor {} initial err: {}, res: {:?}", i, err, res.as_slice());
            }"""

if target in text:
    text = text.replace(target, replacement)
    with open('crates/gneiss-rtk/src/factor_graph/mod.rs', 'w') as f:
        f.write(text)
    print("Replaced!")
else:
    print("Target not found!")
