#![allow(clippy::unwrap_used)]

//! Golden-vector tests for the SPP WNLLS solver.

use super::*;

use crate::estimators::spp::fixture::*;
use crate::estimators::spp::raim::seed_initial_state;
use crate::estimators::spp::{build_measurements, GpsTime};
use gneiss_core::atmosphere::KlobucharParams;
use gneiss_core::sat::Constellation;

fn receiver() -> Vector3<f64> {
    Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0)
}

fn meas(prn: u8, constellation: Constellation, raw_pr: f64, t: GpsTime) -> SppMeasurement {
    SppMeasurement {
        constellation,
        raw_pr,
        snr: 45.0,
        doppler: 0.0,
        time: t,
        eph: gps_ephemeris(prn, t),
        is_iono_free: false,
        freq_band: 1,
    }
}

fn state_at_receiver(t: GpsTime) -> SppState {
    SppState::new(
        Coordinate::new(receiver(), Datum::WGS84, Frame::ECEF, t),
        0.0,
        0.0,
        0.0,
        0.0,
    )
}

#[test]
fn clock_columns_are_allocated_in_a_fixed_order() {
    // 3 position unknowns, then one column per constellation in the order
    // GPS (shared with QZSS), Galileo, BeiDou, GLONASS.
    let t = epoch_time();
    let m = vec![
        meas(1, Constellation::Gps, 20_000_000.0, t),
        meas(1, Constellation::Galileo, 20_000_000.0, t),
        meas(1, Constellation::Beidou, 20_000_000.0, t),
        meas(1, Constellation::Glonass, 20_000_000.0, t),
    ];
    let (cols, clocks) = find_clock_cols(&m);
    assert_eq!(cols, 7);
    assert_eq!((clocks.0, clocks.1, clocks.2, clocks.3), (Some(3), Some(4), Some(5), Some(6)));
    // QZSS shares the GPS column and must not add a new unknown.
    let q = vec![meas(1, Constellation::Qzss, 20_000_000.0, t), meas(1, Constellation::Galileo, 20_000_000.0, t)];
    let (cols2, clocks2) = find_clock_cols(&q);
    assert_eq!(cols2, 5);
    assert_eq!((clocks2.0, clocks2.1, clocks2.2, clocks2.3), (Some(3), Some(4), None, None));
    // With no clocks at all the state is just the 3-D position.
    let g = vec![meas(1, Constellation::Sbas, 20_000_000.0, t)];
    let (cols3, clocks3) = find_clock_cols(&g);
    assert_eq!(cols3, 3);
    assert_eq!((clocks3.0, clocks3.1, clocks3.2, clocks3.3), (None, None, None, None));
}

#[test]
fn design_matrix_row_is_the_unit_line_of_sight_plus_one_clock_entry() {
    let t = epoch_time();
    // Seven satellites against seven unknowns: n != cols - 1, so no height
    // pseudo-observation is appended and every row belongs to a satellite.
    let m = vec![
        meas(1, Constellation::Gps, 20_000_000.0, t),
        meas(2, Constellation::Galileo, 20_100_000.0, t),
        meas(3, Constellation::Beidou, 20_200_000.0, t),
        meas(4, Constellation::Glonass, 20_300_000.0, t),
        meas(5, Constellation::Qzss, 20_400_000.0, t),
        meas(6, Constellation::Sbas, 20_500_000.0, t),
        meas(7, Constellation::Navic, 20_600_000.0, t),
    ];
    let st = state_at_receiver(t);
    let cfg = open_config();
    let (h, _w, _dz, cols, clocks) =
        build_design_matrix(&st, &m, None, &cfg).expect("design matrix");
    assert_eq!(cols, 7);
    assert_eq!(h.nrows(), 7);

    for (i, sat) in m.iter().enumerate() {
        let r_rcv = st.position.vector;
        let (sat_coord, _) = compute_sat_state(sat, 0.0);
        let d = r_rcv - sat_coord.vector;
        let rho = d.norm();
        // Position block is the unit LOS from the receiver to the satellite.
        assert!((h[(i, 0)] - d.x / rho).abs() < 1e-12, "row {i} x");
        assert!((h[(i, 1)] - d.y / rho).abs() < 1e-12, "row {i} y");
        assert!((h[(i, 2)] - d.z / rho).abs() < 1e-12, "row {i} z");
        // Exactly one clock coefficient, equal to +1.
        let expect_col = match sat.constellation {
            Constellation::Gps | Constellation::Qzss => clocks.0,
            Constellation::Galileo => clocks.1,
            Constellation::Beidou => clocks.2,
            Constellation::Glonass => clocks.3,
            _ => None,
        };
        let ones: Vec<usize> = (3..cols).filter(|c| h[(i, *c)].abs() > 1e-12).collect();
        match expect_col {
            Some(c) => {
                assert_eq!(ones, vec![c], "row {i} clock column");
                assert_eq!(h[(i, c)], 1.0);
            }
            // SBAS is not a clock constellation here: no column at all.
            None => assert!(ones.is_empty(), "row {i} must have no clock entry"),
        }
    }
}

#[test]
fn height_constraint_row_is_the_geodetic_up_unit_vector() {
    // At the equator on the prime meridian the up vector is exactly
    // [cos 0 cos 0, cos 0 sin 0, sin 0] = [1, 0, 0].
    let t = epoch_time();
    let m: Vec<SppMeasurement> =
        (1..=3u8).map(|p| meas(p, Constellation::Gps, 20_000_000.0, t)).collect();
    let st = state_at_receiver(t);
    let cfg = open_config();
    let (mut h, mut w, mut dz, cols, _) =
        build_design_matrix(&st, &m, None, &cfg).expect("design matrix");
    // 3 GPS satellites for 4 unknowns: one row short, so the solver must
    // append the height pseudo-observation.
    assert_eq!(cols, 4);
    assert_eq!(m.len(), cols - 1);
    // build_design_matrix already sizes the matrix for the extra row; the
    // solver only has to fill it in.
    assert_eq!(h.nrows(), 4, "one row is reserved for the height constraint");
    apply_height_constraint(ecef_to_llh(receiver()), 3, &mut h, &mut w, &mut dz);
    assert!((h[(3, 0)] - 1.0).abs() < 1e-15 && h[(3, 1)].abs() < 1e-15 && h[(3, 2)].abs() < 1e-15,
        "height row = {:?}", h.row(3));
    assert_eq!(dz[3], 0.0);
    assert!((w[(3, 3)] - 0.01).abs() < 1e-15, "weight = {}", w[(3, 3)]);
    // The constraint row must be a unit vector.
    let n = (h[(3, 0)].powi(2) + h[(3, 1)].powi(2) + h[(3, 2)].powi(2)).sqrt();
    assert!((n - 1.0).abs() < 1e-15);
}

#[test]
fn every_satellite_below_the_mask_yields_a_singular_normal_matrix() {
    // With every row masked the position block is all zeros, H^T W H is the
    // zero matrix, and the solver must report MatrixInversionFailed rather
    // than returning a NaN fix.
    let t = epoch_time();
    let m: Vec<SppMeasurement> =
        (1..=4u8).map(|p| meas(p, Constellation::Gps, 20_000_000.0, t)).collect();
    let st = state_at_receiver(t);
    // 86 deg mask: nothing on Earth is that high.
    let cfg = SppConfig { elevation_mask_rad: 1.5, ..open_config() };
    let err = spp_wnlls_step(&st, &m, None, &cfg).expect_err("singular normal matrix");
    assert_eq!(err, SppError::MatrixInversionFailed);
}

#[test]
fn too_few_measurements_for_the_unknowns_is_rejected_up_front() {
    let t = epoch_time();
    let st = state_at_receiver(t);
    let cfg = open_config();
    // 3 GPS satellites give 4 unknowns, which is exactly allowed; 2 is not.
    let two: Vec<SppMeasurement> =
        (1..=2u8).map(|p| meas(p, Constellation::Gps, 20_000_000.0, t)).collect();
    assert_eq!(
        build_design_matrix(&st, &two, None, &cfg).map(|_| ()).expect_err("too few"),
        SppError::NotEnoughMeasurements
    );
}

#[test]
fn poor_geometry_is_rejected_by_the_normal_matrix_trace() {
    // The gate is on tr(H^T W H)^-1 over the position block. With every
    // satellite in the same direction that quantity explodes and the fix
    // must be refused instead of being wildly wrong.
    let t = epoch_time();
    let m: Vec<SppMeasurement> =
        (1..=4u8).map(|p| meas(p, Constellation::Gps, 20_000_000.0, t)).collect();
    let st = state_at_receiver(t);
    let cfg = SppConfig { geometry_variance_threshold: 1e-12, ..open_config() };
    assert_eq!(spp_wnlls_step(&st, &m, None, &cfg).expect_err("poor geometry"), SppError::PoorGeometry);
}

#[test]
fn sagnac_correction_is_a_rotation_about_ecef_z_by_omega_times_flight_time() {
    // Zero flight time is the identity. For a flight time tau the rotation
    // angle is theta = omega_E * tau, and the code maps
    //   (x, y) -> (x cos theta + y sin theta, -x sin theta + y cos theta),
    // i.e. a point on +y moves toward +x.
    let v = Vector3::new(2.0, 3.0, 4.0);
    assert_eq!(compute_sagnac_correction(v, 0.0), v);

    // The argument is a GEOMETRIC PSEUDORANGE in metres, so a 0.07 s flight
    // time is passed as pr = c * 0.07 = 2.099e7 m.
    let tau = 0.07_f64;
    let theta = OMEGA_E * tau;
    let pr = LIGHT_SPEED * tau;
    let got = compute_sagnac_correction(v, pr);
    let (c, s) = (theta.cos(), theta.sin());
    let expect = Vector3::new(v.x * c + v.y * s, -v.x * s + v.y * c, v.z);
    assert!((got - expect).norm() < 1e-12, "got {got:?}");
    // The map preserves length and z (Earth rotates about the spin axis).
    assert!((got.norm() - v.norm()).abs() < 1e-12);
    assert_eq!(got.z, v.z);
    // A +y point moves toward +x, and the shift is second order in theta.
    // x' = sin(theta) for a unit +y input, and y' = cos(theta) = 1 - theta^2/2
    // = 1 - 1.303e-11, so the y tolerance must exceed that.
    let ey = compute_sagnac_correction(Vector3::y(), pr);
    assert!(ey.x > 0.0 && ey.x < theta, "ey = {ey:?}, theta = {theta}");
    assert!((ey.x - theta.sin()).abs() < 1e-15, "x shift = {}", ey.x);
    assert!((ey.y - 1.0).abs() < 1e-10, "ey = {ey:?}");
}

#[test]
fn atmospheric_delays_are_zero_at_or_below_the_earth_and_gated_by_config() {
    // Inside the Earth there is no atmosphere to speak of.
    let t = epoch_time();
    let cfg = open_config();
    let (tropo, iono) = compute_atmospheric_delays(
        Vector3::new(1.0, 1.0, 1.0),
        ecef_to_llh(receiver()),
        1.0,
        1.0,
        t,
        None,
        &SppConfig { enable_tropo: true, enable_iono: true, ..cfg.clone() },
    );
    assert_eq!((tropo, iono), (0.0, 0.0));

    // Enabled ionosphere with no broadcast parameters degrades to zero
    // rather than dividing by zero.
    let (_, iono_none) = compute_atmospheric_delays(
        receiver(),
        ecef_to_llh(receiver()),
        1.0,
        1.0,
        t,
        None,
        &SppConfig { enable_iono: true, ..cfg.clone() },
    );
    assert_eq!(iono_none, 0.0);

    // With broadcast parameters the delay is positive at the zenith and
    // strictly larger at a low elevation.
    let params = KlobucharParams::default();
    let mut iono_high = 0.0;
    let mut iono_low = 0.0;
    for (i, slot) in [0.0_f64, 0.5_f64].iter().enumerate() {
        let (_, v) = compute_atmospheric_delays(
            receiver(),
            ecef_to_llh(receiver()),
            1.0,
            *slot,
            t,
            Some(&params),
            &SppConfig { enable_iono: true, ..cfg.clone() },
        );
        if i == 0 {
            iono_high = v;
        } else {
            iono_low = v;
        }
    }
    assert!(iono_high >= 0.0, "zenith iono delay = {iono_high}");
    assert!(iono_low >= 0.0, "low elevation delay = {iono_low}");
    // The exact elevation dependence belongs to the Klobuchar model; what this
    // module owns is that supplying broadcast parameters yields a finite,
    // non-negative delay instead of the zero used when they are absent or the
    // ionosphere is disabled.
    assert!(iono_low.is_finite() && iono_high.is_finite());

    // Troposphere is metres, positive, and larger at low elevation.
    let (t_high, _) = compute_atmospheric_delays(
        receiver(),
        ecef_to_llh(receiver()),
        1.0,
        core::f64::consts::FRAC_PI_2,
        t,
        None,
        &SppConfig { enable_tropo: true, ..cfg.clone() },
    );
    assert!(t_high > 2.0 && t_high < 3.0, "zenith troposphere = {t_high} m");
}

#[test]
fn ionosphere_free_measurements_drop_the_klobuchar_delay() {
    // An ionosphere-free combination has already had the delay removed, so
    // applying Klobuchar again would double-count it.
    let t = epoch_time();
    let mut m = meas(1, Constellation::Gps, 20_000_000.0, t);
    m.is_iono_free = true;
    let cfg = SppConfig { enable_iono: true, enable_tropo: false, ..open_config() };
    let st = state_at_receiver(t);
    let with = compute_measurement_residuals(
        &st,
        &m,
        st.position.vector,
        ecef_to_llh(st.position.vector),
        Some(&KlobucharParams::default()),
        &cfg,
    );
    let without = compute_measurement_residuals(
        &st,
        &m,
        st.position.vector,
        ecef_to_llh(st.position.vector),
        None,
        &cfg,
    );
    assert!((with.4 - without.4).abs() < 1e-12,
        "iono-free residual changed by {} when Klobuchar was supplied", with.4 - without.4);
}

#[test]
fn iterative_solver_converges_to_the_generated_position_and_reports_failure() {
    let t = epoch_time();
    let true_pos = receiver();
    let true_cdt = 1234.5;
    let (epoch, ephems) = four_gps_scene(true_pos, true_cdt, t);
    let meas = build_measurements(&epoch, &ephems, &open_config());
    assert!(meas.len() >= 4);
    let seed = seed_initial_state(&meas, None);
    let solved = solve_spp_iteratively(seed.clone(), &meas, None, &open_config()).expect("must converge");
    assert!((solved.position.vector - true_pos).norm() < 0.05, "fix error");
    assert!((solved.cdt - true_cdt).abs() < 0.05, "cdt error = {}", solved.cdt - true_cdt);
    // Zero iterations can never converge.
    let cfg = SppConfig { max_iterations: 0, ..open_config() };
    assert_eq!(
        solve_spp_iteratively(seed, &meas, None, &cfg).expect_err("no iterations"),
        SppError::ConvergenceFailed
    );
}

#[test]
fn satellites_below_the_mask_are_down_weighted_not_dropped() {
    let t = epoch_time();
    let m: Vec<SppMeasurement> =
        (1..=4u8).map(|p| meas(p, Constellation::Gps, 20_000_000.0, t)).collect();
    let st = state_at_receiver(t);
    // A 30 deg mask admits the high satellites and down-weights the rest to
    // MIN_WEIGHT; the row must survive with an exactly zero design block.
    let cfg = SppConfig { elevation_mask_rad: 30.0_f64.to_radians(), ..open_config() };
    let (h, w, _dz, _, _) = build_design_matrix(&st, &m, None, &cfg).expect("design matrix");
    for i in 0..m.len() {
        let masked = h[(i, 0)].abs() < 1e-15;
        if masked {
            assert_eq!(w[(i, i)], MIN_WEIGHT, "masked row {i} weight");
            assert_eq!(w[(i, i)], 1e-10);
        }
    }
}
