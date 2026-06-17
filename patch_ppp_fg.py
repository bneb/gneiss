import re

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'r') as f:
    content = f.read()

# Replace find_worst_outlier
old_find_worst_outlier = """    fn find_worst_outlier(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>) -> Option<gneiss_core::sat::SatelliteId> {
        let final_meas = self.build_measurements(state, sats, x_i, self.max_iterations);
        let mut worst_sat = None;
        let mut max_norm = 15.0;
        for m in &final_meas {
            if m.is_phase {
                let norm = m.res.abs() / m.raw_var.sqrt();
                if norm > max_norm { max_norm = norm; worst_sat = m.sat; }
            }
        }
        worst_sat
    }"""

new_find_worst_outlier = """    fn find_worst_outlier(&self, state: &RtkState, sats: &[ProcessedSat], x_i: &DVector<f64>) -> Option<gneiss_core::sat::SatelliteId> {
        let final_meas = self.build_measurements(state, sats, x_i, self.max_iterations);
        find_worst_outlier_sat(&final_meas)
    }"""

content = content.replace(old_find_worst_outlier, new_find_worst_outlier)

find_worst_outlier_sat_func = """
fn find_worst_outlier_sat(meas: &[Measurement]) -> Option<gneiss_core::sat::SatelliteId> {
    let mut worst_sat = None;
    let mut max_norm = 15.0;
    for m in meas {
        if m.is_phase {
            let norm = m.res.abs() / m.raw_var.sqrt();
            if norm > max_norm { max_norm = norm; worst_sat = m.sat; }
        }
    }
    worst_sat
}
"""

content = content.replace("#[cfg(test)]\nmod nan_tests", find_worst_outlier_sat_func + "\n#[cfg(test)]\nmod nan_tests")

tests = """
#[cfg(test)]
mod mutant_killer_tests {
    use super::*;
    use nalgebra::{DMatrix, DVector, Vector3};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use gneiss_core::sat::{SatelliteId, Constellation};
    use crate::engine::updater::MeasurementType;

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0)),
            1.0
        )
    }

    #[test]
    fn test_find_worst_outlier() {
        let sat_id1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_id2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let sat_id3 = SatelliteId { constellation: Constellation::Gps, prn: 3 };
        
        let meas1 = Measurement { res: 14.0, raw_var: 1.0, is_phase: true, sat: Some(sat_id1), row_idx: 0, typ: MeasurementType::CarrierPhase };
        let meas2 = Measurement { res: 16.0, raw_var: 1.0, is_phase: true, sat: Some(sat_id2), row_idx: 1, typ: MeasurementType::CarrierPhase };
        let meas3 = Measurement { res: 17.0, raw_var: 1.0, is_phase: true, sat: Some(sat_id3), row_idx: 2, typ: MeasurementType::CarrierPhase };
        let meas4 = Measurement { res: 100.0, raw_var: 1.0, is_phase: false, sat: Some(sat_id1), row_idx: 3, typ: MeasurementType::Pseudorange };

        assert_eq!(find_worst_outlier_sat(&[meas1.clone()]), None);
        assert_eq!(find_worst_outlier_sat(&[meas2.clone()]), Some(sat_id2));
        assert_eq!(find_worst_outlier_sat(&[meas2.clone(), meas3.clone()]), Some(sat_id3));
        assert_eq!(find_worst_outlier_sat(&[meas4.clone()]), None);
    }
}
"""

content = content + tests

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'w') as f:
    f.write(content)

