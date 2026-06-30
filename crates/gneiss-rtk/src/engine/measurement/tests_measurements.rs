use crate::engine::config::EkfTuningConfig;
use crate::engine::measurement::{
    compute_dd_doppler, compute_dd_pseudorange, get_sat_state, compute_atmospheric_delays,
    DdComponents, DdContext, DdMeasurementContext, MeasurementEnvironment, SatState,
};
use crate::engine::measurement::types::UpdateGeometry;
use crate::filter::{DdObservation, RtkState};
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

#[test]
fn test_get_sat_state() {
    let time = GpsTime::new(2137, 422922.0);
    let rx_pos = Vector3::new(1000.0, 2000.0, 3000.0);
    let sat = SatelliteId {
        constellation: Constellation::Gps,
        prn: 1,
    };
    let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
        sat,
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
        omega0: 0.0,
        omega_dot: 0.0,
        i0: 1.0,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 0,
        iodc: 0,
    });

    let pr = 20000000.0;
    let (pos, vel) = get_sat_state(&eph, pr, 0.0, time, rx_pos);

    // Assert non-zero output
    println!("pos: {:?}", pos);
    println!("vel: {:?}", vel);
    assert!((pos.x - 5041617.577444584).abs() < 1e-6);
    assert!((pos.y - 17749192.669882875).abs() < 1e-6);
    assert!((pos.z - 18906563.91687504).abs() < 1e-6);
    assert!((vel.x - (-2062.128434083682)).abs() < 1e-6);
    assert!((vel.y - (-1216.6504200626093)).abs() < 1e-6);
    assert!((vel.z - 1738.0997934704972).abs() < 1e-6);

    let (pos0, _vel0) = get_sat_state(&eph, 0.0, 0.0, time, rx_pos);
    assert!((pos.x - pos0.x).abs() > 0.0);
}

#[test]
fn test_compute_atmospheric_delays() {
    let state_time = GpsTime::new(2137, 422922.0);
    let pos_apc = Vector3::new(1000.0, 2000.0, 3000.0);
    let base_coord_vec = Vector3::new(1005.0, 2005.0, 3005.0);
    let sat_vec_rov = Vector3::new(15000000.0, 20000000.0, 30000000.0);
    let ref_sat_vec_rov = Vector3::new(-15000000.0, 20000000.0, -30000000.0);
    let sat_vec_bas = Vector3::new(15000005.0, 20000005.0, 30000005.0);
    let ref_sat_vec_bas = Vector3::new(-15000005.0, 20000005.0, -30000005.0);
    let sat_f1 = 1575.42e6;
    let sat_f2 = 1227.60e6;
    let ref_f1 = 1575.42e6;
    let ref_f2 = 1227.60e6;

    let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(
        state_time,
        pos_apc,
        base_coord_vec,
        sat_vec_rov,
        ref_sat_vec_rov,
        sat_vec_bas,
        ref_sat_vec_bas,
        sat_f1,
        sat_f2,
        ref_f1,
        ref_f2,
        None,
    );
    assert!((tropo_dd - 0.0).abs() < 1e-6);
    assert!((iono_dd_l1 - (-0.00020734960295598626)).abs() < 1e-6);
    assert!((iono_dd_l2 - (-0.0003414932766467871)).abs() < 1e-6);
}

#[test]
fn test_compute_dd_pseudorange() {
    let mut rov_sat = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };
    let mut base_sat = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };
    let mut rov_ref = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };
    let mut ref_base = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 0.0,
        snr: 45.0,
        locktime: None,
    };

    let sat_state = SatState {
        rov_pos: Vector3::new(20000000.0, 0.0, 0.0),
        rov_vel: Vector3::zeros(),
        bas_pos: Vector3::new(0.0, 20000000.0, 0.0),
        bas_vel: Vector3::zeros(),
        f1: 1575.42e6,
        f2: 1227.60e6,
    };
    let ref_state = SatState {
        rov_pos: Vector3::new(20000000.0, 0.0, 0.0),
        rov_vel: Vector3::zeros(),
        bas_pos: Vector3::new(0.0, 20000000.0, 0.0),
        bas_vel: Vector3::zeros(),
        f1: 1575.42e6,
        f2: 1227.60e6,
    };

    let ctx = DdContext {
        rov_sat: &mut rov_sat,
        base_sat: &mut base_sat,
        rov_ref: &mut rov_ref,
        ref_base: &mut ref_base,
        sat_state: &sat_state,
        ref_state: &ref_state,
    };

    let tuning = EkfTuningConfig::default();
    let base_coord = Coordinate::new(
        Vector3::zeros(),
        Datum::WGS84,
        Frame::ECEF,
        GpsTime::new(0, 0.0),
    );
    let env = MeasurementEnvironment {
        ephemerides: &[],
        base_coord: &base_coord,
        base_time: GpsTime::new(0, 0.0),
        lever_arm: Vector3::zeros(),
        omega_b: Vector3::zeros(),
        tuning: &tuning,
        gnn_variances: std::collections::HashMap::new(),
        klobuchar_params: None,
    };

    let ugeom = UpdateGeometry {
        comp_dd: 0.0,
        h_r: Vector3::new(1.0, 0.0, 0.0),
        h_att: Vector3::zeros(),
        h_zwd: 0.0,
        state_size: 22,
    };
    let comps = DdComponents {
        comp_pr_dd: 0.0,
        iono_dd_l1: 0.0,
        iono_dd_l2: 0.0,
        var_factor: 1.0,
        ref_var_factor: 1.0,
        h_zwd: 0.0,
    };
    let mctx = DdMeasurementContext {
        ctx: &ctx,
        geom: &ugeom,
        comps: &comps,
        env: &env,
    };
    let updates = compute_dd_pseudorange(&mctx);
    assert_eq!(updates.len(), 2);
    assert_eq!(updates[0].z, 0.0);
    assert_eq!(updates[1].z, 0.0);
}

#[test]
fn test_compute_dd_doppler() {
    let time = GpsTime::new(2137, 422922.0);
    let state = RtkState::new(
        time,
        Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ),
        10.0,
    );

    let mut rov_sat = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };
    let mut base_sat = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };
    let mut rov_ref = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };
    let mut ref_base = DdObservation {
        sat: SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        },
        pr_l1: 20000000.0,
        pr_l2: Some(20000000.0),
        cp_l1: None,
        cp_l2: None,
        doppler: 100.0,
        snr: 45.0,
        locktime: None,
    };

    let sat_state = SatState {
        rov_pos: Vector3::new(20000000.0, 0.0, 0.0),
        rov_vel: Vector3::zeros(),
        bas_pos: Vector3::new(0.0, 20000000.0, 0.0),
        bas_vel: Vector3::zeros(),
        f1: 1575.42e6,
        f2: 1227.60e6,
    };
    let ref_state = SatState {
        rov_pos: Vector3::new(20000000.0, 0.0, 0.0),
        rov_vel: Vector3::zeros(),
        bas_pos: Vector3::new(0.0, 20000000.0, 0.0),
        bas_vel: Vector3::zeros(),
        f1: 1575.42e6,
        f2: 1227.60e6,
    };

    let ctx = DdContext {
        rov_sat: &mut rov_sat,
        base_sat: &mut base_sat,
        rov_ref: &mut rov_ref,
        ref_base: &mut ref_base,
        sat_state: &sat_state,
        ref_state: &ref_state,
    };

    let tuning = EkfTuningConfig::default();
    let base_coord = Coordinate::new(
        Vector3::zeros(),
        Datum::WGS84,
        Frame::ECEF,
        GpsTime::new(0, 0.0),
    );
    let env = MeasurementEnvironment {
        ephemerides: &[],
        base_coord: &base_coord,
        base_time: GpsTime::new(0, 0.0),
        lever_arm: Vector3::zeros(),
        omega_b: Vector3::zeros(),
        tuning: &tuning,
        gnn_variances: std::collections::HashMap::new(),
        klobuchar_params: None,
    };

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    let update = compute_dd_doppler(
        &ctx,
        Vector3::zeros(),
        Vector3::zeros(),
        Vector3::new(1.0, 0.0, 0.0),
        &r_b_e_rot,
        &state.velocity,
        &env.omega_b,
        &env.lever_arm,
        env.tuning.dop_base_var,
        22,
        1.0,
        1.0,
    );
    let u = update.unwrap();
    assert!(!u.z.is_nan());
    assert!(u.z.is_finite(), "doppler innovation should be finite");
    assert!(
        u.z.abs() < 1000.0,
        "doppler innovation abs should be < 1000 Hz"
    );
    assert!(u.r > 0.0, "doppler variance should be positive");
}
