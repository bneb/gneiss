use crate::engine::config::EkfTuningConfig;
use crate::engine::measurement::{compute_innovations, MeasurementEnvironment};
use crate::filter::{DdObservation, RtkState};
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

#[test]
fn test_measurement_model_against_rtklib_golden_data() {
    let time = GpsTime::new(2137, 422922.0);
    let mut state = RtkState::new(
        time,
        Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ),
        10.0,
    );
    state.velocity = Vector3::new(10.0, -5.0, 2.0);

    let ref_sat = SatelliteId {
        constellation: Constellation::Gps,
        prn: 1,
    };
    let rov_sat1 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 2,
    };
    let rov_sat2 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 3,
    };

    state.add_ambiguity(ref_sat, 1, 5.0, 100.0);
    state.add_ambiguity(rov_sat1, 1, 10.0, 100.0);
    state.add_ambiguity(rov_sat2, 1, 15.0, 100.0);

    state.windup.insert(ref_sat, 0.0);
    state.windup.insert(rov_sat1, 0.0);
    state.windup.insert(rov_sat2, 0.0);

    let ref_rover = DdObservation {
        sat: ref_sat,
        pr_l1: 20000000.0,
        pr_l2: Some(20000001.0),
        cp_l1: Some(100000000.0),
        cp_l2: Some(80000000.0),
        doppler: 100.0,
        snr: 45.0,
        locktime: Some(100),
    };
    let ref_base = DdObservation {
        sat: ref_sat,
        pr_l1: 20005000.0,
        pr_l2: Some(20005001.0),
        cp_l1: Some(100020000.0),
        cp_l2: Some(80016000.0),
        doppler: 10.0,
        snr: 45.0,
        locktime: Some(100),
    };

    let rov1_rover = DdObservation {
        sat: rov_sat1,
        pr_l1: 21000000.0,
        pr_l2: Some(21000001.0),
        cp_l1: Some(105000000.0),
        cp_l2: Some(84000000.0),
        doppler: -50.0,
        snr: 45.0,
        locktime: Some(100),
    };
    let rov1_base = DdObservation {
        sat: rov_sat1,
        pr_l1: 21005000.0,
        pr_l2: Some(21005001.0),
        cp_l1: Some(105020000.0),
        cp_l2: Some(84016000.0),
        doppler: 10.0,
        snr: 45.0,
        locktime: Some(100),
    };

    let rov2_rover = DdObservation {
        sat: rov_sat2,
        pr_l1: 22000000.0,
        pr_l2: Some(22000001.0),
        cp_l1: Some(110000000.0),
        cp_l2: Some(88000000.0),
        doppler: -20.0,
        snr: 45.0,
        locktime: Some(100),
    };
    let rov2_base = DdObservation {
        sat: rov_sat2,
        pr_l1: 22005000.0,
        pr_l2: Some(22005001.0),
        cp_l1: Some(110020000.0),
        cp_l2: Some(88016000.0),
        doppler: 10.0,
        snr: 45.0,
        locktime: Some(100),
    };

    let matched_obs = vec![(rov1_rover, rov1_base), (rov2_rover, rov2_base)];

    let eph_ref = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
        sat: ref_sat,
        toe: time,
        toc: time,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        crs: 0.0,
        crc: 0.0,
        cuc: 0.0,
        cus: 0.0,
        cic: 0.0,
        cis: 0.0,
        m0: 0.0,
        e: 0.01,
        sqrt_a: 5153.6,
        delta_n: 0.0,
        omega0: 0.0,
        omega_dot: 0.0,
        i0: 1.0,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 0,
        iodc: 0,
    });

    let eph_rov1 = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
        sat: rov_sat1,
        toe: time,
        toc: time,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        crs: 0.0,
        crc: 0.0,
        cuc: 0.0,
        cus: 0.0,
        cic: 0.0,
        cis: 0.0,
        m0: 1.0,
        e: 0.01,
        sqrt_a: 5153.6,
        delta_n: 0.0,
        omega0: 0.5,
        omega_dot: 0.0,
        i0: 1.0,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 0,
        iodc: 0,
    });

    let eph_rov2 = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
        sat: rov_sat2,
        toe: time,
        toc: time,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        crs: 0.0,
        crc: 0.0,
        cuc: 0.0,
        cus: 0.0,
        cic: 0.0,
        cis: 0.0,
        m0: 2.0,
        e: 0.01,
        sqrt_a: 5153.6,
        delta_n: 0.0,
        omega0: 1.0,
        omega_dot: 0.0,
        i0: 1.0,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 0,
        iodc: 0,
    });

    let ephemerides = vec![eph_ref, eph_rov1, eph_rov2];
    let base_coord = Coordinate::new(
        Vector3::new(1005.0, 2005.0, 3005.0),
        Datum::WGS84,
        Frame::ECEF,
        time,
    );

    let config = crate::engine::EngineConfig::default();
    let env = MeasurementEnvironment {
        ephemerides: &ephemerides,
        base_coord: &base_coord,
        base_time: base_coord.epoch,
        lever_arm: Vector3::zeros(),
        omega_b: Vector3::zeros(),
        tuning: &config.tuning,
        gnn_variances: std::collections::HashMap::new(),
        klobuchar_params: None,
    };
    let updates = compute_innovations(&mut state, &matched_obs, &ref_rover, &ref_base, &env)
        .unwrap();
    let z = updates.z;
    let r = updates.r;

    println!("Z: {:?}", z);

    // Lock in the golden Z vector (updated for iterative ecef_to_llh refinement)
    assert!((z[0] - 2.95577).abs() < 1e-3, "z[0]={}", z[0]);
    assert!((z[1] - 2.95627).abs() < 1e-3, "z[1]={}", z[1]);
    assert!((z[2] - -2.04577).abs() < 1e-3, "z[2]={}", z[2]);
    assert!((z[3] - 19.36016).abs() < 1e-3, "z[3]={}", z[3]);
    assert!((z[4] - 9.97132).abs() < 1e-3, "z[4]={}", z[4]);
    assert!((z[5] - 9.97206).abs() < 1e-3, "z[5]={}", z[5]);
    assert!((z[6] - -0.03096).abs() < 1e-3, "z[6]={}", z[6]);
    assert!((z[7] - 16.01613).abs() < 1e-3, "z[7]={}", z[7]);

    // Lock in the golden R diagonal
    assert!(r[0] >= 16.0);
    assert!(r[1] >= 16.0);
    assert!(r[2] >= 0.0001);
    assert!(r[3] >= 0.1);
    assert!(r[4] >= 16.0);
    assert!(r[5] >= 16.0);
    assert!(r[6] >= 0.0001);
    assert!(r[7] >= 0.1);
}

#[test]
fn test_compute_innovations() {
    let time = GpsTime::new(2137, 422922.0);
    let mut state = RtkState::new(
        time,
        Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ),
        10.0,
    );

    let sat1 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 1,
    };
    let sat2 = SatelliteId {
        constellation: Constellation::Gps,
        prn: 2,
    };

    let rov_sat = DdObservation {
        sat: sat1,
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: Some(100000000.0),
        cp_l2: Some(80000000.0),
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };
    let base_sat = DdObservation {
        sat: sat1,
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: Some(100000000.0),
        cp_l2: Some(80000000.0),
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };
    let rov_ref = DdObservation {
        sat: sat2,
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: Some(100000000.0),
        cp_l2: Some(80000000.0),
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };
    let ref_base = DdObservation {
        sat: sat2,
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: Some(100000000.0),
        cp_l2: Some(80000000.0),
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };

    let eph_ref = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
        sat: sat2,
        toe: time,
        toc: time,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        crs: 0.0,
        crc: 0.0,
        cuc: 0.0,
        cus: 0.0,
        cic: 0.0,
        cis: 0.0,
        m0: 0.0,
        e: 0.01,
        sqrt_a: 5153.6,
        delta_n: 0.0,
        omega0: 0.0,
        omega_dot: 0.0,
        i0: 1.0,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 0,
        iodc: 0,
    });

    let eph_rov = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
        sat: sat1,
        toe: time,
        toc: time,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        crs: 0.0,
        crc: 0.0,
        cuc: 0.0,
        cus: 0.0,
        cic: 0.0,
        cis: 0.0,
        m0: 1.0,
        e: 0.01,
        sqrt_a: 5153.6,
        delta_n: 0.0,
        omega0: 0.5,
        omega_dot: 0.0,
        i0: 1.0,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 0,
        iodc: 0,
    });

    state.add_ambiguity(sat1, 1, 0.0, 1.0);
    state.add_ambiguity(sat2, 1, 0.0, 1.0);
    state.add_ambiguity(sat1, 2, 0.0, 1.0);
    state.add_ambiguity(sat2, 2, 0.0, 1.0);

    let tuning = EkfTuningConfig::default();
    let base_coord = Coordinate::new(
        Vector3::zeros(),
        Datum::WGS84,
        Frame::ECEF,
        GpsTime::new(0, 0.0),
    );
    let env = MeasurementEnvironment {
        ephemerides: &[eph_ref, eph_rov],
        base_coord: &base_coord,
        base_time: GpsTime::new(0, 0.0),
        lever_arm: Vector3::zeros(),
        omega_b: Vector3::zeros(),
        tuning: &tuning,
        gnn_variances: std::collections::HashMap::new(),
        klobuchar_params: None,
    };

    let group = vec![(rov_sat, base_sat)];

    let res = compute_innovations(&mut state, &group, &rov_ref, &ref_base, &env);
    assert!(res.is_some());
    let updates = res.unwrap();
    let (z_vals, h_rows, r_vals, meas_type) = (updates.z, updates.h, updates.r, updates.mt);
    assert!(z_vals.len() > 0);
    assert_eq!(h_rows.len(), z_vals.len());
    assert_eq!(r_vals.len(), z_vals.len());
    assert_eq!(meas_type.len(), z_vals.len());

    println!("h_rows[0] = {:?}", h_rows[0]);
    assert_eq!(h_rows[0][0], 0.993379331940329);
}
