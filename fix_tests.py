import re

with open("crates/gneiss-rtk/src/engine/tests_updater.rs", "r") as f:
    updater = f.read()
updater = updater.replace("apply_state_correction, apply_state_vector, extract_state_vector, joseph_update", "apply_state_correction, joseph_update")
with open("crates/gneiss-rtk/src/engine/tests_updater.rs", "w") as f:
    f.write(updater)

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    ppp_fg = f.read()
ppp_fg = ppp_fg.replace("let h = build_h_row(&los, 4.0, Some(18), 19);", "let h = build_h_row(&los, 4.0, Some(18), 19, gneiss_core::sat::Constellation::Gps);")
ppp_fg = ppp_fg.replace("let h2 = build_h_row(&los, 4.0, None, 16);", "let h2 = build_h_row(&los, 4.0, None, 16, gneiss_core::sat::Constellation::Gps);")
ppp_fg = ppp_fg.replace("apply_state_vector(&mut state2, &x);", "apply_state_vector(&mut state2, &x, state2.covariance.clone());")
with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(ppp_fg)

with open("crates/gneiss-rtk/src/engine/predictor.rs", "r") as f:
    predictor = f.read()
predictor = predictor.replace("use crate::engine::{EngineConfig, DynamicsModel, TuningParams};", "use crate::engine::{EngineConfig, DynamicsModel, EngineMode};")
predictor = predictor.replace("tuning: TuningParams::default(),", "")
predictor = predictor.replace("elevation_mask: 0.2,", "")
predictor = predictor.replace("snr_mask: 30.0,", "")
predictor = predictor.replace("let mut state = RtkState::new();", "let mut state = RtkState::new(gneiss_core::time::GpsTime::new(2000, 0.0), gneiss_core::coords::Coordinate::new(nalgebra::Vector3::zeros(), gneiss_core::coords::Datum::WGS84, gneiss_core::coords::Frame::ECEF, gneiss_core::time::GpsTime::new(2000, 0.0)), 10.0);")
with open("crates/gneiss-rtk/src/engine/predictor.rs", "w") as f:
    f.write(predictor)

