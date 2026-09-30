//! Shared test utilities, mathematical reference oracles, and opaque-box contract fixtures.

use core::marker::PhantomData;
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

pub use super::types::*;

pub const WGS84_A: f64 = 6378137.0;
pub const WGS84_F: f64 = 1.0 / 298.257223563;
pub const WGS84_B: f64 = WGS84_A * (1.0 - WGS84_F);
pub const SPEED_OF_LIGHT: f64 = 299792458.0;
pub const OMEGA_EARTH: f64 = 7.2921151467e-5;
pub const GPS_LEAP_SECONDS_2017: i32 = 18;
pub const WEEK_SECONDS: f64 = 604800.0;

pub fn omega_ie_ecef() -> Vector3<f64> {
    Vector3::new(0.0, 0.0, OMEGA_EARTH)
}

pub fn skew_symmetric(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
}

pub fn r_b_e_from_rpy(roll: f64, pitch: f64, yaw: f64) -> UnitQuaternion<f64> {
    UnitQuaternion::from_euler_angles(roll, pitch, yaw)
}

pub fn ecef_to_llh_analytical(p: &Vector3<f64>) -> Vector3<f64> {
    let e2 = (WGS84_A * WGS84_A - WGS84_B * WGS84_B) / (WGS84_A * WGS84_A);
    let p_dist = p.x.hypot(p.y);
    let (phi, h) = if p_dist < 1e-9 {
        let phi = if p.z >= 0.0 { core::f64::consts::FRAC_PI_2 } else { -core::f64::consts::FRAC_PI_2 };
        (phi, p.z.abs() - WGS84_B)
    } else {
        let mut phi = p.z.atan2(p_dist * (1.0 - e2));
        for _ in 0..5 {
            let sin_phi = phi.sin();
            let n = WGS84_A / (1.0 - e2 * sin_phi * sin_phi).sqrt();
            phi = (p.z + e2 * n * sin_phi).atan2(p_dist);
        }
        let sin_phi = phi.sin();
        let n = WGS84_A / (1.0 - e2 * sin_phi * sin_phi).sqrt();
        (phi, p_dist / phi.cos() - n)
    };
    let lam = p.y.atan2(p.x);
    Vector3::new(phi, lam, h)
}

pub fn ecef_to_ned_matrix_analytical(llh: &Vector3<f64>) -> Matrix3<f64> {
    let (s_lat, c_lat) = (llh.x.sin(), llh.x.cos());
    let (s_lon, c_lon) = (llh.y.sin(), llh.y.cos());
    Matrix3::new(
        -s_lat * c_lon, -s_lat * s_lon,  c_lat,
        -s_lon,          c_lon,          0.0,
        -c_lat * c_lon, -c_lat * s_lon, -s_lat,
    )
}

pub fn los_and_az_el(
    station_ecef: &Vector3<f64>,
    sat_ecef: &Vector3<f64>,
) -> (Vector3<f64>, f64, f64, f64) {
    let diff = sat_ecef - station_ecef;
    let range = diff.norm();
    let unit_los = diff / range;
    let llh = ecef_to_llh_analytical(station_ecef);
    let r_ned = ecef_to_ned_matrix_analytical(&llh);
    let ned = r_ned * unit_los;
    let el = (-ned.z).atan2(ned.x.hypot(ned.y));
    let mut az = ned.y.atan2(ned.x);
    if az < 0.0 {
        az += 2.0 * core::f64::consts::PI;
    }
    (unit_los, range, el, az)
}

pub struct LocalTangentPlane<R: ReferenceFrame> {
    origin: EcefPos<R>,
    r_ned: Matrix3<f64>,
}

impl<R: ReferenceFrame> LocalTangentPlane<R> {
    pub fn from_origin(origin: EcefPos<R>) -> Self {
        let llh = ecef_to_llh_analytical(origin.coords());
        let r_ned = ecef_to_ned_matrix_analytical(&llh);
        Self { origin, r_ned }
    }
    pub fn to_ned(&self, target: EcefPos<R>) -> SpatialVector<Ned> {
        let d = self.r_ned * (target.coords() - self.origin.coords());
        SpatialVector::from_vector(d)
    }
    pub fn to_enu(&self, target: EcefPos<R>) -> SpatialVector<Enu> {
        let d = self.r_ned * (target.coords() - self.origin.coords());
        SpatialVector::new(d.y, d.x, -d.z)
    }
    pub fn from_enu(&self, enu: SpatialVector<Enu>) -> EcefPos<R> {
        let enu_v = enu.vector();
        let ned_v = Vector3::new(enu_v.y, enu_v.x, -enu_v.z);
        let ecef_delta = self.r_ned.transpose() * ned_v;
        Point3::from_coords(self.origin.coords() + ecef_delta)
    }
    pub fn project_cov(&self, cov: &SpatialCovariance<Ecef<R>>) -> EnuCovariance {
        let cov_ned = self.r_ned * cov.matrix() * self.r_ned.transpose();
        let mut m = Matrix3::zeros();
        m[(0, 0)] = cov_ned[(1, 1)];
        m[(1, 1)] = cov_ned[(0, 0)];
        m[(2, 2)] = cov_ned[(2, 2)];
        m[(0, 1)] = cov_ned[(1, 0)];
        m[(1, 0)] = cov_ned[(0, 1)];
        m[(0, 2)] = -cov_ned[(1, 2)];
        m[(2, 0)] = -cov_ned[(2, 1)];
        m[(1, 2)] = -cov_ned[(0, 2)];
        m[(2, 1)] = -cov_ned[(2, 0)];
        SpatialCovariance::from_matrix(m)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleDiffGeometry<R: ReferenceFrame> {
    pub u_rov_sat: Vector3<f64>,
    pub u_rov_ref: Vector3<f64>,
    pub range_rov_sat_m: f64,
    pub range_rov_ref_m: f64,
    pub el_rov_sat_rad: f64,
    pub az_rov_sat_rad: f64,
    pub el_rov_ref_rad: f64,
    pub az_rov_ref_rad: f64,
    pub u_bas_sat: Vector3<f64>,
    pub u_bas_ref: Vector3<f64>,
    pub range_bas_sat_m: f64,
    pub range_bas_ref_m: f64,
    pub el_bas_sat_rad: f64,
    pub az_bas_sat_rad: f64,
    pub el_bas_ref_rad: f64,
    pub az_bas_ref_rad: f64,
    pub base_dd_range_m: f64,
    pub geometric_dd_m: f64,
    pub delta_u_rov: Vector3<f64>,
    _frame: PhantomData<R>,
}

impl<R: ReferenceFrame> DoubleDiffGeometry<R> {
    pub fn compute(
        rov_pos: EcefPos<R>,
        bas_pos: EcefPos<R>,
        sat_pos: EcefPos<R>,
        ref_sat_pos: EcefPos<R>,
    ) -> Self {
        let (u_rov_sat, r_rov_sat, el_rov_sat, az_rov_sat) = los_and_az_el(rov_pos.coords(), sat_pos.coords());
        let (u_rov_ref, r_rov_ref, el_rov_ref, az_rov_ref) = los_and_az_el(rov_pos.coords(), ref_sat_pos.coords());
        let (u_bas_sat, r_bas_sat, el_bas_sat, az_bas_sat) = los_and_az_el(bas_pos.coords(), sat_pos.coords());
        let (u_bas_ref, r_bas_ref, el_bas_ref, az_bas_ref) = los_and_az_el(bas_pos.coords(), ref_sat_pos.coords());
        let base_dd = r_bas_sat - r_bas_ref;
        let rov_dd = r_rov_sat - r_rov_ref;
        let geometric_dd = rov_dd - base_dd;
        let delta_u_rov = u_rov_sat - u_rov_ref;
        Self {
            u_rov_sat, u_rov_ref, range_rov_sat_m: r_rov_sat, range_rov_ref_m: r_rov_ref,
            el_rov_sat_rad: el_rov_sat, az_rov_sat_rad: az_rov_sat,
            el_rov_ref_rad: el_rov_ref, az_rov_ref_rad: az_rov_ref,
            u_bas_sat, u_bas_ref, range_bas_sat_m: r_bas_sat, range_bas_ref_m: r_bas_ref,
            el_bas_sat_rad: el_bas_sat, az_bas_sat_rad: az_bas_sat,
            el_bas_ref_rad: el_bas_ref, az_bas_ref_rad: az_bas_ref,
            base_dd_range_m: base_dd, geometric_dd_m: geometric_dd,
            delta_u_rov, _frame: PhantomData,
        }
    }
}
