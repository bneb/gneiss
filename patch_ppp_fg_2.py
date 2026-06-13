with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    content = f.read()

content = content.replace("gneiss_core::ephemeris::Constellation", "gneiss_core::sat::Constellation")
content = content.replace("use nalgebra::{DMatrix, DVector, Vector3};", "use nalgebra::{DMatrix, DVector, Vector3, UnitQuaternion};")
content = content.replace("apply_state_vector(state, &x_i);", "apply_state_vector(state, &x_i, state.covariance.clone());")

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(content)

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    content = f.read()

content = content.replace("gneiss_core::ephemeris::Constellation", "gneiss_core::sat::Constellation")

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(content)

