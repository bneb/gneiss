//! Tests for the bidirectional RTS inertial smoother.
//!
//! Kalman updates are checked against the textbook form.  Straight after
//! `InertialFilterState::new` the state is P = 0.01 I, so with R = 0.01 I the
//! gain is K = P (P + R)^-1 = 0.5 I exactly; after a `predict` the blocks are
//! no longer equal and the test derives K block by block.
#![allow(clippy::unwrap_used)]

use super::*;

    #[test]
    fn test_smoother_empty_snapshots() {
        let res = run_inertial_rts_smoother(&[]);
        assert!(res.is_empty());
    }

    #[test]
    fn test_outage_bridging_eliminates_forward_drift() {
        let init_pos = Vector3::new(4_000_000.0, 0.0, 4_000_000.0);
        let init_vel = Vector3::new(10.0, 0.0, 0.0);
        let mut filter = InertialFilterState::new(init_pos, init_vel, UnitQuaternion::identity());

        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;
        preint.dp = -0.5 * filter.gravity_ecef * preint.dt * preint.dt;
        preint.dv = -filter.gravity_ecef * preint.dt;

        // 1. Initial 3 epochs with GNSS fixes
        for i in 0..3 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            filter.predict(t, &preint);
            let gnss_p = init_pos + init_vel * (i as f64 + 1.0);
            filter.update_gnss(gnss_p, init_vel, 0.0004, 0.0025);
        }

        // 2. 10-second complete GNSS outage (epochs 3..13)
        // With constant velocity + small synthetic bias drift
        for i in 3..13 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            let mut biased_preint = preint.clone();
            // Inject 0.05 m/s^2 bias drift
            biased_preint.dp.x += 0.05 * 0.5;
            biased_preint.dv.x += 0.05;
            filter.predict(t, &biased_preint);
            filter.update_nhc(0.01, 0.01);
        }

        let forward_drift_at_end = (filter.x.fixed_rows::<3>(0).into_owned() - (init_pos + init_vel * 13.0)).norm();
        assert!(forward_drift_at_end > 1.0, "Forward pass must exhibit uncorrected bias drift without smoothing");

        // 3. Post-outage GNSS re-lock (epochs 13..15)
        for i in 13..16 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            filter.predict(t, &preint);
            let gnss_p = init_pos + init_vel * (i as f64 + 1.0);
            filter.update_gnss(gnss_p, init_vel, 0.0004, 0.0025);
        }

        // 4. Run RTS backward smoothing
        let smoothed = run_inertial_rts_smoother(&filter.history);
        assert_eq!(smoothed.len(), 16);

        // Verify maximum drift across all 10 outage epochs stays < 0.50 m
        let mut max_outage_drift = 0.0;
        for (i, ep) in smoothed.iter().enumerate().take(13).skip(3) {
            let truth_p = init_pos + init_vel * (i as f64 + 1.0);
            let err = (ep.position_ecef - truth_p).norm();
            if err > max_outage_drift {
                max_outage_drift = err;
            }
        }
        assert!(
            max_outage_drift < 0.50,
            "Maximum drift during 10-second outage was {:.4}m, which exceeds the 0.50m threshold",
            max_outage_drift
        );
    }

    #[test]
    fn test_stationary_zupt_locks_velocity_and_position() {
        let init_pos = Vector3::new(4_000_000.0, 0.0, 4_000_000.0);
        let init_vel = Vector3::zeros();
        let mut filter = InertialFilterState::new(init_pos, init_vel, UnitQuaternion::identity());

        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;
        preint.dp = -0.5 * filter.gravity_ecef * preint.dt * preint.dt;
        preint.dv = -filter.gravity_ecef * preint.dt;

        // Simulate 10 seconds of complete outage while vehicle is stationary (ZUPT active)
        for i in 0..10 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            let mut biased = preint.clone();
            // Inject small sensor noise / bias
            biased.dp.x += 0.02 * 0.5;
            biased.dv.x += 0.02;
            filter.predict(t, &biased);
            filter.update_zupt(0.0001); // 1 cm/s 1-sigma
        }

        let smoothed = run_inertial_rts_smoother(&filter.history);
        for ep in smoothed {
            let drift = (ep.position_ecef - init_pos).norm();
            assert!(
                drift < 0.20,
                "Position drift {:.4}m exceeded 0.20m during stationary ZUPT outage",
                drift
            );
            assert!(
                ep.velocity_ecef.norm() < 0.02,
                "Velocity was not locked during stationary ZUPT outage"
            );
        }

    }

    // ---- seeding and prediction -------------------------------------------

    /// Receiver on the line x = z, so `normalize` is exactly (1, 0, 1)/sqrt(2)
    /// and the gravity vector is exactly 9.80665 / sqrt(2) along each axis.
    fn seed() -> InertialFilterState {
        InertialFilterState::new(
            Vector3::new(4_000_000.0, 0.0, 4_000_000.0),
            Vector3::new(3.0, 0.0, 0.0),
            UnitQuaternion::identity(),
        )
    }

    #[test]
    fn initial_state_matches_the_documented_seed() {
        let f = seed();
        assert_eq!(f.x.fixed_rows::<3>(0).into_owned(), Vector3::new(4_000_000.0, 0.0, 4_000_000.0));
        assert_eq!(f.x.fixed_rows::<3>(3).into_owned(), Vector3::new(3.0, 0.0, 0.0));
        // 10 cm position and 10 cm/s velocity uncertainty on every axis.
        for i in 0..6 {
            assert!((f.p[(i, i)] - 0.01).abs() < 1e-15, "P[{i},{i}] = {}", f.p[(i, i)]);
        }
        let g = 9.80665 / std::f64::consts::SQRT_2;
        let want = Vector3::new(-g, 0.0, -g);
        assert!((f.gravity_ecef - want).norm() < 1e-12, "gravity = {:?}, want {want:?}", f.gravity_ecef);
        assert!(f.history.is_empty());
    }

    #[test]
    fn predict_applies_the_constant_gravity_kinematics() {
        // With an identity attitude and an empty preintegration (dp = dv = 0):
        //   p' = p + v dt + 0.5 g dt^2,   v' = v + g dt
        // dt = 2 s, v = (3, 0, 0), g = (-g0, 0, -g0) with g0 = 9.80665/sqrt2:
        //   v' = (3 - 2 g0, 0, -2 g0),  p' = p + (6, 0, 0) + 2 g
        let mut f = seed();
        let mut preint = ImuPreintegration::new();
        preint.dt = 2.0;
        let f_mat = f.predict(GpsTime::new(2000, 100.0), &preint);
        let g0 = 9.80665 / std::f64::consts::SQRT_2;
        let want_v = Vector3::new(3.0 - 2.0 * g0, 0.0, -2.0 * g0);
        // 0.5 * dt^2 = 2, so the gravity term contributes -2 * g0 per axis.
        let want_p = Vector3::new(4_000_000.0 + 6.0 - 2.0 * g0, 0.0, 4_000_000.0 - 2.0 * g0);
        let got_v = f.x.fixed_rows::<3>(3).into_owned();
        let got_p = f.x.fixed_rows::<3>(0).into_owned();
        assert!((got_v - want_v).norm() < 1e-9, "v' = {got_v:?}, want {want_v:?}");
        assert!((got_p - want_p).norm() < 1e-9, "p' = {got_p:?}, want {want_p:?}");

        // The transition matrix couples position to velocity with dt.
        for i in 0..3 {
            assert!((f_mat[(i, i)] - 1.0).abs() < 1e-15 && (f_mat[(i, i + 3)] - 2.0).abs() < 1e-15);
            assert!((f_mat[(i + 3, i + 3)] - 1.0).abs() < 1e-15, "velocity has no position coupling");
        }
        assert_eq!(f.history.len(), 1, "every predict must be recorded for the backward pass");
        assert!(!f.history[0].is_gnss_available, "predict alone is not a GNSS fix");
    }

    // ---- measurement updates ----------------------------------------------

    #[test]
    fn gnss_update_uses_the_exact_block_kalman_gain() {
        // After one predict with dt = 1 from P = 0.01 I, F = [[I, I], [0, I]]
        // and Q = diag(0.01 I, 0.05 I) give
        //   P = 1/100 [[3, 1], [1, 6]].
        // With R = 0.01 I the innovation covariance keeps the off-diagonal
        // block, so S = 1/100 [[4, 1], [1, 7]] and
        //   det = (4*7 - 1)/100^2 = 27/10 000,
        //   K = P S^-1 = 1/27 [[20, 1], [1, 23]],
        //   P' = (I - K) P = 1/2700 [[20, 1], [1, 23]].
        let mut f = seed();
        f.predict(GpsTime::new(2000, 100.0), &ImuPreintegration { dt: 1.0, ..ImuPreintegration::new() });
        let before = f.x;
        let z_pos = Vector3::new(4_000_010.0, 5.0, 4_000_000.0);
        f.update_gnss(z_pos, Vector3::zeros(), 0.01, 0.01);

        let inn_p = z_pos - before.fixed_rows::<3>(0).into_owned();
        let inn_v = Vector3::zeros() - before.fixed_rows::<3>(3).into_owned();
        let want_p = before.fixed_rows::<3>(0).into_owned() + inn_p * (20.0 / 27.0) + inn_v * (1.0 / 27.0);
        let got_p = f.x.fixed_rows::<3>(0).into_owned();
        assert!((got_p - want_p).norm() < 1e-9, "p = {got_p:?}, want {want_p:?}");

        for (r, c, num) in [(0, 0, 20), (0, 3, 1), (3, 0, 1), (3, 3, 23)] {
            let want = f64::from(num) / 2700.0;
            assert!((f.p[(r, c)] - want).abs() < 1e-12, "P'[{r},{c}] = {}, want {want}", f.p[(r, c)]);
        }
        assert!(f.history.last().unwrap().is_gnss_available, "the update must mark the epoch as fixed");
    }

    #[test]
    fn zupt_halves_the_velocity_when_state_and_measurement_agree() {
        // p_v = 0.01, R = 0.01 -> K = 0.5, so v' = v - 0.5 v = v / 2.
        let mut f = seed();
        f.x.fixed_rows_mut::<3>(3).copy_from(&Vector3::new(2.0, -4.0, 6.0));
        f.update_zupt(0.01);
        let got = f.x.fixed_rows::<3>(3).into_owned();
        assert!((got - Vector3::new(1.0, -2.0, 3.0)).norm() < 1e-9, "v' = {got:?}, want (1, -2, 3)");
        // p_vv' = p_v - K p_v = 0.005
        for i in 3..6 {
            assert!((f.p[(i, i)] - 0.005).abs() < 1e-12, "P'[{i},{i}] = {}", f.p[(i, i)]);
        }
    }

    #[test]
    fn nhc_leaves_forward_motion_alone_and_damps_lateral_motion() {
        // With an identity attitude the body frame is the ECEF frame, so the
        // constraint is simply v_y = v_z = 0 and v_x is unconstrained.
        let mut f = seed();
        f.x.fixed_rows_mut::<3>(3).copy_from(&Vector3::new(5.0, 1.0, 1.0));
        f.update_nhc(0.01, 0.01);
        let v = f.x.fixed_rows::<3>(3).into_owned();
        assert!((v.x - 5.0).abs() < 1e-9, "forward speed must be untouched, got {v:?}");
        // p_vy = p_vz = 0.01, R = 0.01 -> K = 0.5, so both halve.
        assert!((v.y - 0.5).abs() < 1e-9 && (v.z - 0.5).abs() < 1e-9, "lateral/vertical must halve, got {v:?}");
    }

    #[test]
    fn covariance_stays_symmetric_and_positive_definite_across_a_full_cycle() {
        // P is built from F P F' + Q and from (I - K H) P, both of which are
        // symmetric by construction; an asymmetric or indefinite P would mean
        // a transposed term somewhere.
        let mut f = seed();
        let mut preint = ImuPreintegration::new();
        preint.dt = 0.5;
        for i in 0..25 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            f.predict(t, &preint);
            match i % 3 {
                0 => f.update_gnss(Vector3::new(4_000_000.0 + i as f64, 0.0, 4_000_000.0), Vector3::zeros(), 4e-4, 2.5e-3),
                1 => f.update_nhc(0.01, 0.01),
                _ => f.update_zupt(1e-4),
            }
            assert!((f.p - f.p.transpose()).norm() < 1e-15, "P lost symmetry at epoch {i}");
            let e = nalgebra::linalg::SymmetricEigen::new(f.p).eigenvalues;
            assert!(e.min() > 0.0, "P lost positive definiteness at epoch {i}: min eig {}", e.min());
        }
        assert_eq!(f.history.len(), 25);
    }

    // ---- backward pass -----------------------------------------------------

    #[test]
    fn a_single_snapshot_smoother_is_a_passthrough() {
        let mut f = seed();
        f.predict(GpsTime::new(2000, 100.0), &ImuPreintegration { dt: 1.0, ..ImuPreintegration::new() });
        f.update_gnss(Vector3::new(4_000_001.0, 0.0, 4_000_000.0), Vector3::zeros(), 0.01, 0.01);
        let out = run_inertial_rts_smoother(&f.history);
        assert_eq!(out.len(), 1);
        let last = f.history[0].clone();
        assert!((out[0].position_ecef - last.x_post.fixed_rows::<3>(0).into_owned()).norm() < 1e-12);
        assert!((out[0].velocity_ecef - last.x_post.fixed_rows::<3>(3).into_owned()).norm() < 1e-12);
        assert!((out[0].cov_position - last.p_post.fixed_view::<3, 3>(0, 0).into_owned()).norm() < 1e-12);
    }

    #[test]
    fn rts_smoothing_never_increases_the_covariance() {
        // The RTS gain is a projection, so P_s = P + C dP C' with
        // dP = P_s(next) - P_pred(next) must stay no worse than P_post here.
        let mut f = seed();
        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;
        preint.dp = Vector3::zeros();
        preint.dv = Vector3::zeros();
        for i in 0..6 {
            let t = GpsTime::new(2000, 100.0 + i as f64);
            f.predict(t, &preint);
            f.update_gnss(Vector3::new(4_000_000.0 + 0.5 * (i as f64 + 1.0), 0.0, 4_000_000.0), Vector3::zeros(), 4e-4, 2.5e-3);
        }
        let out = run_inertial_rts_smoother(&f.history);
        assert_eq!(out.len(), 6);
        for (i, ep) in out.iter().enumerate() {
            let forward = f.history[i].p_post;
            let diag = |m: Cov6, off: usize| (0..3).map(|k| m[(off + k, off + k)]).sum::<f64>();
            let trace_f: f64 = diag(forward, 0) + diag(forward, 3);
            let trace_s: f64 = ep.cov_position.diagonal().sum() + ep.cov_velocity.diagonal().sum();
            assert!(trace_s <= trace_f + 1e-9, "epoch {i}: smoothed trace {trace_s} exceeds forward {trace_f}");
        }
    }
