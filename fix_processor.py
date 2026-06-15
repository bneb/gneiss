with open("crates/gneiss-rtk/src/engine/processor.rs", "r") as f:
    text = f.read()

# Fix destructuring
text = text.replace(
    "if let Some(MeasurementModel { z_vec: z_safe, h_mat: h_safe, r_mat: r_safe, meas_types: type_safe }) = crate::engine::measurement::build_measurement_model(",
    "if let Some((z_safe, h_safe, r_safe, type_safe)) = crate::engine::measurement::build_measurement_model("
)

# Fix updater call missing is_tightly_coupled
text = text.replace(
    "crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, self.config.spp_consistency_threshold_m, None,  &self.config.tuning)",
    "crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, self.config.spp_consistency_threshold_m, None, self.config.mode.is_tightly_coupled(), &self.config.tuning)"
)

text = text.replace(
    "Some(type_stripped)\n                        }.as_deref(),  &self.config.tuning)",
    "Some(type_stripped)\n                        }.as_deref(), self.config.mode.is_tightly_coupled(), &self.config.tuning)"
)

with open("crates/gneiss-rtk/src/engine/processor.rs", "w") as f:
    f.write(text)

