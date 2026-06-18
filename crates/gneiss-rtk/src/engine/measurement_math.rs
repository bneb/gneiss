use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

pub struct WindupUpdates {
    pub w_sat: f64,
    pub w_ref: f64,
    pub w_bas_sat: f64,
    pub w_bas_ref: f64,
}

pub fn get_sat_state(
    eph: &Ephemeris,
    pr: f64,
    t_rx: GpsTime,
    rx_pos: Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    let tau_pr = pr / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    assert!(tau_pr.abs() < 1000.0, "tau_pr must be < 1000.0s");
    let t_tx_nom = GpsTime::new(t_rx.week, t_rx.tow - tau_pr);
    let (_, _, dt_s, _) = eph.position(t_tx_nom);
    let t_tx_true = GpsTime::new(t_rx.week, t_rx.tow - tau_pr - dt_s);
    let (raw_vec, raw_vel, _, _) = eph.position(t_tx_true);

    let mut sat_pos = raw_vec;
    let mut sat_vel = raw_vel;
    for _ in 0..2 {
        let geometric_range = (sat_pos - rx_pos).norm();
        let true_tau = geometric_range / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let theta = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * true_tau;
        let cos_t = f64::cos(theta);
        let sin_t = f64::sin(theta);
        sat_pos = nalgebra::Vector3::new(
            raw_vec.x * cos_t + raw_vec.y * sin_t,
            -raw_vec.x * sin_t + raw_vec.y * cos_t,
            raw_vec.z,
        );
        sat_vel = nalgebra::Vector3::new(
            raw_vel.x * cos_t + raw_vel.y * sin_t,
            -raw_vel.x * sin_t + raw_vel.y * cos_t,
            raw_vel.z,
        );
    }
    (sat_pos, sat_vel)
}

#[allow(clippy::too_many_arguments)]
pub fn compute_atmospheric_delays(
    state_time: GpsTime,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    sat_vec_rov: Vector3<f64>,
    ref_sat_vec_rov: Vector3<f64>,
    sat_vec_bas: Vector3<f64>,
    ref_sat_vec_bas: Vector3<f64>,
    sat_f1: f64,
    sat_f2: f64,
    ref_f1: f64,
    ref_f2: f64,
) -> (f64, f64, f64) {
    let tropo_params = gneiss_core::atmosphere::TropoParams::default();
    let iono_params = gneiss_core::atmosphere::KlobucharParams::default();

    let base_llh = gneiss_core::coords::ecef_to_llh(base_coord_vec);
    let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);

    let (az_rov_sat, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_vec_rov);
    let (az_rov_ref, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, pos_apc, ref_sat_vec_rov);
    let (az_bas_sat, el_bas_sat) =
        gneiss_core::coords::az_el(base_llh, base_coord_vec, sat_vec_bas);
    let (az_bas_ref, el_bas_ref) =
        gneiss_core::coords::az_el(base_llh, base_coord_vec, ref_sat_vec_bas);

    let tropo_rov_sat = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        rov_llh,
        el_rov_sat,
    );
    let tropo_rov_ref = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        rov_llh,
        el_rov_ref,
    );
    let tropo_bas_sat = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        base_llh,
        el_bas_sat,
    );
    let tropo_bas_ref = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        base_llh,
        el_bas_ref,
    );
    let tropo_dd = (tropo_rov_sat - tropo_rov_ref) - (tropo_bas_sat - tropo_bas_ref);

    let iono_rov_sat = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
        &iono_params,
        rov_llh,
        az_rov_sat,
        el_rov_sat,
        state_time,
    );
    let iono_rov_ref = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
        &iono_params,
        rov_llh,
        az_rov_ref,
        el_rov_ref,
        state_time,
    );
    let iono_bas_sat = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
        &iono_params,
        base_llh,
        az_bas_sat,
        el_bas_sat,
        state_time,
    );
    let iono_bas_ref = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
        &iono_params,
        base_llh,
        az_bas_ref,
        el_bas_ref,
        state_time,
    );
    let iono_dd_l1 = (iono_rov_sat - iono_rov_ref) - (iono_bas_sat - iono_bas_ref);

    let f_ratio_sat_l2 = (sat_f1 / sat_f2).powi(2);
    let f_ratio_ref_l2 = (ref_f1 / ref_f2).powi(2);
    let iono_dd_l2 = (iono_rov_sat * f_ratio_sat_l2 - iono_rov_ref * f_ratio_ref_l2)
        - (iono_bas_sat * f_ratio_sat_l2 - iono_bas_ref * f_ratio_ref_l2);

    (tropo_dd, iono_dd_l1, iono_dd_l2)
}

#[allow(clippy::too_many_arguments)]
pub fn compute_phase_windup(
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    sun_pos: Vector3<f64>,
    rov_pos_sat: Vector3<f64>,
    rov_pos_ref: Vector3<f64>,
    bas_pos_sat: Vector3<f64>,
    bas_pos_ref: Vector3<f64>,
    prev_w_sat: f64,
    prev_w_ref: f64,
    prev_w_bas_sat: f64,
    prev_w_bas_ref: f64,
) -> WindupUpdates {
    let w_sat = gneiss_core::windup::phase_windup(rov_pos_sat, sun_pos, pos_apc, prev_w_sat);
    let w_ref = gneiss_core::windup::phase_windup(rov_pos_ref, sun_pos, pos_apc, prev_w_ref);
    let w_bas_sat =
        gneiss_core::windup::phase_windup(bas_pos_sat, sun_pos, base_coord_vec, prev_w_bas_sat);
    let w_bas_ref =
        gneiss_core::windup::phase_windup(bas_pos_ref, sun_pos, base_coord_vec, prev_w_bas_ref);

    WindupUpdates {
        w_sat,
        w_ref,
        w_bas_sat,
        w_bas_ref,
    }
}

pub fn range_attitude_jacobian(lever_ecef: &Vector3<f64>, h_r: &Vector3<f64>) -> Vector3<f64> {
    lever_ecef.cross(h_r)
}

pub fn doppler_attitude_jacobian(
    r_b_e: &nalgebra::Matrix3<f64>,
    omega_b: &Vector3<f64>,
    lever_arm: &Vector3<f64>,
    h_r: &Vector3<f64>,
) -> Vector3<f64> {
    let a = r_b_e * omega_b.cross(lever_arm);
    a.cross(h_r)
}

pub struct VarianceFactors {
    pub snr_rov_sat: f64,
    pub snr_rov_ref: f64,
    pub el_rov_sat: f64,
    pub el_rov_ref: f64,
    pub el_bas_sat: f64,
    pub el_bas_ref: f64,
    pub snr_a: f64,
    pub snr_b: f64,
    pub gnn_var_sat: Option<f64>,
    pub gnn_var_ref: Option<f64>,
}

pub fn compute_variance_factors(v: &VarianceFactors) -> (f64, f64) {
    let ref_var = if let Some(var) = v.gnn_var_ref {
        var
    } else {
        gneiss_core::variance::observation_variance(v.snr_rov_ref, v.el_rov_ref, v.snr_a, v.snr_b)
            + gneiss_core::variance::elevation_variance_scale(v.el_bas_ref)
    };

    let sat_var = if let Some(var) = v.gnn_var_sat {
        var
    } else {
        gneiss_core::variance::observation_variance(v.snr_rov_sat, v.el_rov_sat, v.snr_a, v.snr_b)
            + gneiss_core::variance::elevation_variance_scale(v.el_bas_sat)
    };

    (sat_var + ref_var, ref_var)
}

pub fn compute_zwd_mapping(el_rov_sat: f64, el_rov_ref: f64, zwd: f64) -> (f64, f64) {
    let m_w_rov_sat = 1.0 / el_rov_sat.max(0.001).sin();
    let m_w_rov_ref = 1.0 / el_rov_ref.max(0.001).sin();
    let h_zwd = m_w_rov_sat - m_w_rov_ref;
    let zwd_dd = h_zwd * zwd;
    (h_zwd, zwd_dd)
}

pub fn compute_geometric_dd(
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    sat_vec_rov: Vector3<f64>,
    ref_sat_vec_rov: Vector3<f64>,
    sat_vec_bas: Vector3<f64>,
    ref_sat_vec_bas: Vector3<f64>,
) -> f64 {
    ((pos_apc - sat_vec_rov).norm() - (pos_apc - ref_sat_vec_rov).norm())
        - ((base_coord_vec - sat_vec_bas).norm() - (base_coord_vec - ref_sat_vec_bas).norm())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::ephemeris::Ephemeris;
    use gneiss_core::sat::SatelliteId;
    use gneiss_core::time::GpsTime;
    use nalgebra::{Matrix3, Vector3};

    #[test]
    fn test_compute_zwd_mapping() {
        let el_sat = 45.0_f64.to_radians();
        let el_ref = 60.0_f64.to_radians();
        let zwd = 0.1;
        let (h_zwd, zwd_dd) = compute_zwd_mapping(el_sat, el_ref, zwd);

        let m_w_sat = 1.0 / el_sat.sin();
        let m_w_ref = 1.0 / el_ref.sin();
        assert!((h_zwd - (m_w_sat - m_w_ref)).abs() < 1e-9);
        assert!((zwd_dd - (h_zwd * zwd)).abs() < 1e-9);

        // Edge case: elevation < 0.001
        let (h_zwd_edge, _) = compute_zwd_mapping(0.0001, 0.0001, zwd);
        let max_val = 1.0 / 0.001_f64.sin();
        assert!((h_zwd_edge - (max_val - max_val)).abs() < 1e-9);
    }

    #[test]
    fn test_compute_geometric_dd() {
        let pos_apc = Vector3::new(10.0, 20.0, 30.0);
        let base_coord_vec = Vector3::new(10.0, 20.0, 30.0);
        let sat_vec_rov = Vector3::new(100.0, 200.0, 300.0);
        let ref_sat_vec_rov = Vector3::new(-100.0, -200.0, -300.0);

        let dd = compute_geometric_dd(
            pos_apc,
            base_coord_vec,
            sat_vec_rov,
            ref_sat_vec_rov,
            sat_vec_rov,
            ref_sat_vec_rov,
        );
        // Because rover and base are at the same position, DD should exactly cancel out to 0.
        assert!(dd.abs() < 1e-9);

        let pos_apc2 = pos_apc + Vector3::new(1.0, 0.0, 0.0);
        let dd2 = compute_geometric_dd(
            pos_apc2,
            base_coord_vec,
            sat_vec_rov,
            ref_sat_vec_rov,
            sat_vec_rov,
            ref_sat_vec_rov,
        );
        let d_rov_sat = (pos_apc2 - sat_vec_rov).norm();
        let d_rov_ref = (pos_apc2 - ref_sat_vec_rov).norm();
        let d_bas_sat = (base_coord_vec - sat_vec_rov).norm();
        let d_bas_ref = (base_coord_vec - ref_sat_vec_rov).norm();
        assert!((dd2 - ((d_rov_sat - d_rov_ref) - (d_bas_sat - d_bas_ref))).abs() < 1e-9);
    }

    #[test]
    fn test_compute_variance_factors() {
        let (var, ref_var) = compute_variance_factors(&VarianceFactors {
            snr_rov_sat: 40.0,
            snr_rov_ref: 45.0,
            el_rov_sat: 45.0_f64.to_radians(),
            el_rov_ref: 60.0_f64.to_radians(),
            el_bas_sat: 45.0_f64.to_radians(),
            el_bas_ref: 60.0_f64.to_radians(),
            snr_a: 100.0,
            snr_b: 100.0,
            gnn_var_sat: None,
            gnn_var_ref: None,
        });

        let expected_ref_var =
            gneiss_core::variance::observation_variance(45.0, 60.0_f64.to_radians(), 100.0, 100.0)
                + gneiss_core::variance::elevation_variance_scale(60.0_f64.to_radians());
        let expected_var =
            gneiss_core::variance::observation_variance(40.0, 45.0_f64.to_radians(), 100.0, 100.0)
                + gneiss_core::variance::elevation_variance_scale(45.0_f64.to_radians())
                + expected_ref_var;

        assert!((ref_var - expected_ref_var).abs() < 1e-9);
        assert!((var - expected_var).abs() < 1e-9);
    }

    #[test]
    fn test_range_attitude_jacobian() {
        let lever = Vector3::new(1.0, 2.0, 3.0);
        let h_r = Vector3::new(4.0, 5.0, 6.0);
        let expected = lever.cross(&h_r);
        let result = range_attitude_jacobian(&lever, &h_r);
        assert!((result - expected).norm() < 1e-9);
    }

    #[test]
    fn test_doppler_attitude_jacobian() {
        let r_b_e = Matrix3::identity();
        let omega = Vector3::new(0.1, 0.2, 0.3);
        let lever = Vector3::new(1.0, 2.0, 3.0);
        let h_r = Vector3::new(4.0, 5.0, 6.0);

        let a = r_b_e * omega.cross(&lever);
        let expected = a.cross(&h_r);
        let result = doppler_attitude_jacobian(&r_b_e, &omega, &lever, &h_r);
        assert!((result - expected).norm() < 1e-9);
    }

    #[test]
    fn test_compute_phase_windup() {
        let updates = compute_phase_windup(
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::new(1e11, 0.0, 0.0),
            Vector3::new(100.0, 0.0, 0.0),
            Vector3::new(0.0, 100.0, 0.0),
            Vector3::new(100.0, 0.0, 0.0),
            Vector3::new(0.0, 100.0, 0.0),
            0.1,
            0.2,
            0.3,
            0.4,
        );
        // Phase windup function exists, we just need to ensure the wrapper passes values identically.
        // It's tested elsewhere, but we ensure our mapping doesn't mangle it.
        let exp_w_sat = gneiss_core::windup::phase_windup(
            Vector3::new(100.0, 0.0, 0.0),
            Vector3::new(1e11, 0.0, 0.0),
            Vector3::new(10.0, 20.0, 30.0),
            0.1,
        );
        assert!(
            (updates.w_sat - exp_w_sat).abs() < 1e-9
                || (updates.w_sat.is_nan() && exp_w_sat.is_nan())
        );
    }

    #[test]
    fn test_compute_atmospheric_delays() {
        let state_time = GpsTime::new(2000, 100000.0);
        let pos_apc = Vector3::new(10.0, 20.0, 30.0);
        let base_coord_vec = Vector3::new(10.0, 20.0, 30.0);
        let sat_vec_rov = Vector3::new(100.0, 200.0, 300.0);
        let ref_sat_vec_rov = Vector3::new(-100.0, -200.0, -300.0);
        let sat_vec_bas = Vector3::new(100.0, 200.0, 300.0);
        let ref_sat_vec_bas = Vector3::new(-100.0, -200.0, -300.0);

        let sat_f1 = 1575.42e6;
        let sat_f2 = 1227.60e6;
        let ref_f1 = 1575.42e6;
        let ref_f2 = 1227.60e6;

        let (tropo, iono1, iono2) = compute_atmospheric_delays(
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
        );
        // Since rover and base are at same position and vectors are identical, double difference should be 0.
        assert!(tropo.abs() < 1e-9);
        assert!(iono1.abs() < 1e-9);
        assert!(iono2.abs() < 1e-9);
    }
}

/// DD measurements sharing a reference satellite have correlated noise equal to
/// the reference satellite's measurement variance.
pub const DD_CROSS_CORRELATION_SCALE: f64 = 1.0;

pub fn build_dense_covariance_matrix(
    r_diagonals: &[f64],
    meas_types: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> nalgebra::DMatrix<f64> {
    let mut r_mat =
        nalgebra::DMatrix::from_diagonal(&nalgebra::DVector::from_row_slice(r_diagonals));
    for i in 0..meas_types.len() {
        for j in (i + 1)..meas_types.len() {
            if meas_types[i].1 == meas_types[j].1
                && meas_types[i].0.constellation == meas_types[j].0.constellation
            {
                let cov = meas_types[i].2.min(meas_types[j].2) * DD_CROSS_CORRELATION_SCALE;
                r_mat[(i, j)] = cov;
                r_mat[(j, i)] = cov;
            }
        }
    }
    r_mat
}

#[test]
fn test_get_sat_state_catches_tau_pr_mutation() {
    let time = GpsTime::new(2137, 422922.0);
    let rx_pos = Vector3::new(1000.0, 2000.0, 3000.0);
    let sat = gneiss_core::sat::SatelliteId {
        constellation: gneiss_core::sat::Constellation::Gps,
        prn: 1,
    };

    // Use a non-zero af1 so that t_tx_nom affects dt_s
    let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
        sat,
        toe: time,
        toc: time,
        af0: 0.0,
        af1: 0.01,
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
    let (pos1, _) = get_sat_state(&eph, pr, time, rx_pos);

    // Exact assertion to catch symmetric swap of + and -
    assert!(
        (pos1.x - 5041616.189715548).abs() < 1e-4,
        "pos1.x was {}",
        pos1.x
    );
    assert!(
        (pos1.y - 17749191.855105825).abs() < 1e-4,
        "pos1.y was {}",
        pos1.y
    );
}

#[test]
fn test_get_sat_state_catches_rotation_mutation() {
    let time = GpsTime::new(2137, 422922.0);
    let rx_pos = Vector3::new(1000.0, 2000.0, 3000.0);
    let sat = gneiss_core::sat::SatelliteId {
        constellation: gneiss_core::sat::Constellation::Gps,
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

    // Use a massive pr so true_tau is huge (e.g. 100 seconds)
    // theta = 7.29e-5 * 100 = 7.29e-3 rad. cos(theta) = 0.999973.
    // If cos(theta) is divided instead of multiplied, vel x changes by raw_vel * (1/cos - cos)
    // For vel = 2000, change is 2000 * 5.4e-5 = 0.1 m/s.
    // We will assert with < 1e-9 tolerance to ensure the exact math is used.
    let pr = 3e10; // 100 seconds
    let (_, vel) = get_sat_state(&eph, pr, time, rx_pos);

    // Actually, let's just make the assertion extremely tight on the normal PR
    let pr_normal = 20000000.0;
    let (_, vel_normal) = get_sat_state(&eph, pr_normal, time, rx_pos);
    assert!((vel_normal.x - (-2062.128434083682)).abs() < 1e-12);
    assert!((vel_normal.y - (-1216.6504200626093)).abs() < 1e-12);
}

#[test]
fn test_compute_atmospheric_delays_catches_tropo_mutation() {
    let state_time = GpsTime::new(2137, 422922.0);

    let pos_apc = Vector3::new(6378137.0, 0.0, 0.0);
    let base_coord_vec = Vector3::new(6378137.0, 1000.0, 0.0);

    let sat_vec_rov = Vector3::new(26000000.0, 0.0, 0.0);
    let ref_sat_vec_rov = Vector3::new(26000000.0, 10000000.0, 0.0);

    let sat_vec_bas = Vector3::new(26000000.0, 0.0, 0.0);
    let ref_sat_vec_bas = Vector3::new(26000000.0, 10000000.0, 0.0);

    let (tropo_dd, _, _) = compute_atmospheric_delays(
        state_time,
        pos_apc,
        base_coord_vec,
        sat_vec_rov,
        ref_sat_vec_rov,
        sat_vec_bas,
        ref_sat_vec_bas,
        1.0,
        1.0,
        1.0,
        1.0,
    );

    assert!((tropo_dd - 0.0).abs() > 1e-6);
    // Specifically, let's just assert that it exactly matches the logic:
    let tropo_params = gneiss_core::atmosphere::TropoParams::default();
    let base_llh = gneiss_core::coords::ecef_to_llh(base_coord_vec);
    let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
    let (_, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_vec_rov);
    let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, pos_apc, ref_sat_vec_rov);
    let (_, el_bas_sat) = gneiss_core::coords::az_el(base_llh, base_coord_vec, sat_vec_bas);
    let (_, el_bas_ref) = gneiss_core::coords::az_el(base_llh, base_coord_vec, ref_sat_vec_bas);

    let t1 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        rov_llh,
        el_rov_sat,
    );
    let t2 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        rov_llh,
        el_rov_ref,
    );
    let t3 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        base_llh,
        el_bas_sat,
    );
    let t4 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(
        &tropo_params,
        base_llh,
        el_bas_ref,
    );

    let expected_dd = (t1 - t2) - (t3 - t4);
    assert!((tropo_dd - expected_dd).abs() < 1e-9);
}
