//! Inertial Outage Simulation Tests for Sprint 43.
//!
//! Verifies that tightly-coupled GNSS/INS smoothing with Non-Holonomic Constraints (NHC)
//! bridges a 10-second complete satellite outage with maximum drift < 0.50 m.

#[cfg(test)]
mod tests {
    use nalgebra::{UnitQuaternion, Vector3};
    use gneiss_core::time::GpsTime;
    use gneiss_rtk::swfg::imu_preintegration::smoother::{
        run_inertial_rts_smoother, InertialFilterState,
    };
    use gneiss_rtk::swfg::imu_preintegration::ImuPreintegration;

    #[test]
    fn test_ten_second_complete_outage_drift_under_half_meter() {
        // Trajectory: vehicle travelling eastward at 15 m/s (54 km/h)
        let init_pos = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
        let speed = 15.0; // m/s
        let heading_dir = Vector3::new(0.6, -0.7, 0.38).normalize();
        let vel = heading_dir * speed;

        let rot = UnitQuaternion::rotation_between(&Vector3::new(1.0, 0.0, 0.0), &heading_dir)
            .unwrap_or_else(UnitQuaternion::identity);
        let mut filter = InertialFilterState::new(init_pos, vel, rot);

        let dt = 1.0;
        let mut nominal_preint = ImuPreintegration::new();
        nominal_preint.dt = dt;
        let r_e2b = rot.inverse();
        nominal_preint.dp = r_e2b * (-0.5 * filter.gravity_ecef * dt * dt);
        nominal_preint.dv = r_e2b * (-filter.gravity_ecef * dt);

        let total_epochs = 30;
        let outage_start = 10;
        let outage_end = 20;

        for epoch in 0..total_epochs {
            let t = GpsTime::new(2200, 300000.0 + epoch as f64);
            let truth_pos = init_pos + vel * (epoch as f64 + 1.0);

            // Injected sensor bias (0.04 m/s^2 forward bias in body frame)
            let mut preint = nominal_preint.clone();
            preint.dp += Vector3::new(0.04 * 0.5 * dt * dt, 0.0, 0.0);
            preint.dv += Vector3::new(0.04 * dt, 0.0, 0.0);

            filter.predict(t, &preint);

            if epoch < outage_start || epoch >= outage_end {
                // GNSS fixed observations with 2 cm position noise and 5 cm/s velocity noise
                filter.update_gnss(truth_pos, vel, 0.0004, 0.0025);
            } else {
                // Complete satellite dropout: only NHC active
                filter.update_nhc(0.01, 0.01);
            }
        }

        // Verify that forward-only drift at the end of the outage grew significantly
        let fwd_pos_at_outage_end = filter.history[outage_end - 1].x_post.fixed_rows::<3>(0).into_owned();
        let truth_at_outage_end = init_pos + vel * (outage_end as f64);
        let fwd_drift = (fwd_pos_at_outage_end - truth_at_outage_end).norm();
        assert!(
            fwd_drift > 1.5,
            "Forward dead-reckoning must demonstrate uncorrected bias drift (was {:.3}m)",
            fwd_drift
        );

        // Run RTS backward smoothing
        let smoothed = run_inertial_rts_smoother(&filter.history);
        assert_eq!(smoothed.len(), total_epochs);

        // Check maximum drift during the entire 10-second outage
        let mut max_smoothed_outage_drift = 0.0;
        for (epoch, ep) in smoothed.iter().enumerate().take(outage_end).skip(outage_start) {
            let truth = init_pos + vel * (epoch as f64 + 1.0);
            let drift = (ep.position_ecef - truth).norm();
            if drift > max_smoothed_outage_drift {
                max_smoothed_outage_drift = drift;
            }
        }

        println!("Forward uncorrected drift: {:.4} m", fwd_drift);
        println!("Smoothed maximum outage drift: {:.4} m", max_smoothed_outage_drift);

        assert!(
            max_smoothed_outage_drift < 0.50,
            "Maximum drift during 10-second outage was {:.4}m, exceeding the 0.50m Sprint 43 exit criterion",
            max_smoothed_outage_drift
        );
    }
}
