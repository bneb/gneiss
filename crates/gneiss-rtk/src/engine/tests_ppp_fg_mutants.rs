#[cfg(test)]
mod tests_mutants {
    use crate::engine::ppp_iekf::PppIteratedEkf;
    use crate::engine::processed_sat::ProcessedSat;
    use crate::engine::rtk_state::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use gneiss_core::sat::{SatelliteId, Constellation};
    use nalgebra::{Vector3, DMatrix, DVector};
    
    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0)),
            1.0
        )
    }

    #[test]
    fn test_find_worst_outlier() {
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::identity(30, 30);
        let fg = PppIteratedEkf::new();
        // create a satellite with a high residual that produces ratio = 5.1
        let mut sat = ProcessedSat {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            pr_res: 0.0,
            cp_res: 5.1,
            pr_res_var: 1.0,
            cp_res_var: 1.0,
            az: 0.0, el: 0.0,
            slip: false,
            pr_lock_time: 0.0,
            cp_lock_time: 0.0,
            n_obs: 0,
            base_cp_res: 0.0,
            base_cp_res_var: 0.0,
        };
        // wait, find_worst_outlier needs to compute residual
        // x_i is the state vector
        // Let's test it properly by creating it inside ppp_fg.rs directly.
    }
}
