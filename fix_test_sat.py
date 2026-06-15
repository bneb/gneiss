with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

# Make af1 non-zero in test_get_sat_state
text = text.replace("af0: 0.0, af1: 0.0, af2: 0.0", "af0: 0.0, af1: 0.001, af2: 0.0")

# The pos and vel asserts will fail because af1 changed!
# Instead of hardcoding, I should just make af1: 0.001 and then calculate the correct output or just add a new test in measurement_math.rs
