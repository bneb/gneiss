with open("crates/gneiss-rtk/src/engine/processor.rs", "r") as f:
    text = f.read()

text = text.replace("MeasurementEnvironment, MeasurementModel}", "MeasurementEnvironment}")

with open("crates/gneiss-rtk/src/engine/processor.rs", "w") as f:
    f.write(text)


with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    text = f.read()

text = text.replace("crate::engine::measurement::build_dense_covariance_matrix", "crate::engine::measurement_math::build_dense_covariance_matrix")
text = text.replace("crate::engine::measurement::MeasurementModel", "crate::engine::measurement_math::MeasurementModel")

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(text)


with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

text = text.replace("env.tuning.snr_a", "snr_a").replace("env.tuning.snr_b", "snr_b")

if "pub use crate::engine::measurement_math::*;" not in text:
    text = "pub use crate::engine::measurement_math::*;\n" + text

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(text)
