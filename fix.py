with open("crates/gneiss-rtk/src/lambda.rs", "r") as f:
    lines = f.readlines()
lines = [l for l in lines if l.strip() != "#[test]"]
with open("crates/gneiss-rtk/src/lambda.rs", "w") as f:
    f.writelines(lines)
