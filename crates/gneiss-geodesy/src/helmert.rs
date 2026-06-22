use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// 14-parameter Helmert Transformation for coordinates between reference frames
/// taking into account epoch propagation (tectonic motion).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HelmertParams {
    /// Translation in X (meters)
    pub tx: f64,
    /// Translation in Y (meters)
    pub ty: f64,
    /// Translation in Z (meters)
    pub tz: f64,
    /// Rotation around X axis (milli-arcseconds)
    pub rx: f64,
    /// Rotation around Y axis (milli-arcseconds)
    pub ry: f64,
    /// Rotation around Z axis (milli-arcseconds)
    pub rz: f64,
    /// Scale factor (parts per billion)
    pub s: f64,

    /// Rate of change of tx (meters/year)
    pub dtx: f64,
    /// Rate of change of ty (meters/year)
    pub dty: f64,
    /// Rate of change of tz (meters/year)
    pub dtz: f64,
    /// Rate of change of rx (milli-arcseconds/year)
    pub drx: f64,
    /// Rate of change of ry (milli-arcseconds/year)
    pub dry: f64,
    /// Rate of change of rz (milli-arcseconds/year)
    pub drz: f64,
    /// Rate of change of scale (ppb/year)
    pub ds: f64,

    /// Reference epoch for the parameters (e.g. 2010.0)
    pub ref_epoch: f64,
}

pub trait GeodeticTransform {
    fn apply(&self, coord: gneiss_core::coords::Coordinate) -> gneiss_core::coords::Coordinate;
}

impl GeodeticTransform for HelmertParams {
    fn apply(&self, coord: gneiss_core::coords::Coordinate) -> gneiss_core::coords::Coordinate {
        if coord.frame != gneiss_core::coords::Frame::ECEF {
            // Helmert is defined for ECEF vectors
            return coord;
        }

        let obs_epoch = coord.epoch.to_fractional_year();
        let new_vector = self.transform(coord.vector, obs_epoch);

        gneiss_core::coords::Coordinate::new(new_vector, coord.datum, coord.frame, coord.epoch)
    }
}

impl HelmertParams {
    /// Computes the translation vector at the given time difference.
    fn compute_translation(&self, dt: f64) -> Vector3<f64> {
        Vector3::new(
            self.tx + self.dtx * dt,
            self.ty + self.dty * dt,
            self.tz + self.dtz * dt,
        )
    }

    /// Computes the rotation vector (in radians) at the given time difference.
    fn compute_rotation(&self, dt: f64) -> Vector3<f64> {
        let mas2rad = gneiss_core::constants::MILLIARCSEC_TO_RAD;
        Vector3::new(
            (self.rx + self.drx * dt) * mas2rad,
            (self.ry + self.dry * dt) * mas2rad,
            (self.rz + self.drz * dt) * mas2rad,
        )
    }

    /// Transforms an ECEF vector from the source frame to the target frame at a specific observation epoch.
    pub fn transform(&self, ecef: Vector3<f64>, obs_epoch: f64) -> Vector3<f64> {
        let dt = obs_epoch - self.ref_epoch;
        let t = self.compute_translation(dt);
        let r = self.compute_rotation(dt);
        let scale = 1.0 + (self.s + self.ds * dt) * 1e-9;

        // Apply Helmert Transformation (position convention)
        Vector3::new(
            t.x + scale * (ecef.x - r.z * ecef.y + r.y * ecef.z),
            t.y + scale * (r.z * ecef.x + ecef.y - r.x * ecef.z),
            t.z + scale * (-r.y * ecef.x + r.x * ecef.y + ecef.z),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;
    use nalgebra::Vector3;

    #[test]
    fn test_helmert_itrf2014_to_itrf2020() {
        // Transformation parameters from ITRF2014 to ITRF2020 at epoch 2015.0
        // Provided by IERS
        // tx, ty, tz in mm -> convert to m
        let params = HelmertParams {
            tx: -0.0014,
            ty: -0.0012,
            tz: 0.0012,
            rx: 0.0,
            ry: 0.0,
            rz: 0.0,
            s: 0.0,

            dtx: 0.0,
            dty: -0.0001,
            dtz: 0.0002,
            drx: 0.0,
            dry: 0.0,
            drz: 0.0,
            ds: 0.0,

            ref_epoch: 2015.0,
        };

        // A coordinate in ITRF2014 (approx location on Earth's surface)
        let ecef_2014 = Vector3::new(4027893.0, 307041.0, 4919475.0);

        // We want to transform to ITRF2020 at observation epoch 2025.0
        let transformed = params.transform(ecef_2014, 2025.0);

        // Expected translation at 2025.0:
        // dt = 2025.0 - 2015.0 = 10.0 years
        // T_X = tx + dt * dtx = -0.0014 + 10 * 0.0 = -0.0014 m
        // T_Y = ty + dt * dty = -0.0012 + 10 * -0.0001 = -0.0022 m
        // T_Z = tz + dt * dtz = 0.0012 + 10 * 0.0002 = 0.0032 m
        // For zero rotations and zero scale, X_new = X + T_X ...

        let expected_x = 4027893.0 - 0.0014;
        let expected_y = 307041.0 - 0.0022;
        let expected_z = 4919475.0 + 0.0032;

        assert!((transformed.x - expected_x).abs() < 1e-9);
        assert!((transformed.y - expected_y).abs() < 1e-9);
        assert!((transformed.z - expected_z).abs() < 1e-9);
    }

    #[test]
    fn test_helmert_with_rotations() {
        let params = HelmertParams {
            tx: 1.0,
            ty: 2.0,
            tz: 3.0,
            rx: 1000.0, // 1000 mas = 1 arcsec
            ry: 2000.0, // 2 arcsec
            rz: 3000.0, // 3 arcsec
            s: 10.0,    // 10 ppb
            dtx: 0.0,
            dty: 0.0,
            dtz: 0.0,
            drx: 0.0,
            dry: 0.0,
            drz: 0.0,
            ds: 0.0,
            ref_epoch: 2000.0,
        };

        let ecef = Vector3::new(6000000.0, 1000000.0, 2000000.0);
        let transformed = params.transform(ecef, 2000.0);

        let mas2rad = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);
        let rx = 1000.0 * mas2rad;
        let ry = 2000.0 * mas2rad;
        let rz = 3000.0 * mas2rad;
        let scale = 1.0 + 10.0 * 1e-9;

        let exp_x = 1.0 + scale * (ecef.x - rz * ecef.y + ry * ecef.z);
        let exp_y = 2.0 + scale * (rz * ecef.x + ecef.y - rx * ecef.z);
        let exp_z = 3.0 + scale * (-ry * ecef.x + rx * ecef.y + ecef.z);

        assert!((transformed.x - exp_x).abs() < 1e-9);
        assert!((transformed.y - exp_y).abs() < 1e-9);
        assert!((transformed.z - exp_z).abs() < 1e-9);
    }

    #[test]
    fn test_helmert_transform_negative_dt() {
        // Observation epoch before ref_epoch: dt negative
        let params = HelmertParams {
            tx: 1.0,
            ty: 2.0,
            tz: 3.0,
            rx: 1000.0,
            ry: 2000.0,
            rz: 3000.0,
            s: 10.0,
            dtx: 0.1,
            dty: 0.2,
            dtz: 0.3,
            drx: 10.0,
            dry: 20.0,
            drz: 30.0,
            ds: 1.0,
            ref_epoch: 2010.0,
        };

        let ecef = Vector3::new(4000000.0, 500000.0, 4800000.0);
        // dt = 2000.0 - 2010.0 = -10.0 years
        let transformed = params.transform(ecef, 2000.0);

        let dt = -10.0;
        let tx = 1.0 + 0.1 * dt;
        let ty = 2.0 + 0.2 * dt;
        let tz = 3.0 + 0.3 * dt;
        let mas2rad = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);
        let rx = (1000.0 + 10.0 * dt) * mas2rad;
        let ry = (2000.0 + 20.0 * dt) * mas2rad;
        let rz = (3000.0 + 30.0 * dt) * mas2rad;
        let scale = 1.0 + (10.0 + 1.0 * dt) * 1e-9;

        let exp_x = tx + scale * (ecef.x - rz * ecef.y + ry * ecef.z);
        let exp_y = ty + scale * (rz * ecef.x + ecef.y - rx * ecef.z);
        let exp_z = tz + scale * (-ry * ecef.x + rx * ecef.y + ecef.z);

        assert!((transformed.x - exp_x).abs() < 1e-9);
        assert!((transformed.y - exp_y).abs() < 1e-9);
        assert!((transformed.z - exp_z).abs() < 1e-9);
    }

    #[test]
    fn test_helmert_transform_at_ref_epoch() {
        // dt = 0 -> no rate terms contribute
        let params = HelmertParams {
            tx: 1.5,
            ty: 2.5,
            tz: 3.5,
            rx: 500.0,
            ry: 600.0,
            rz: 700.0,
            s: 5.0,
            dtx: 100.0,
            dty: 200.0,
            dtz: 300.0,
            drx: 999.0,
            dry: 999.0,
            drz: 999.0,
            ds: 999.0,
            ref_epoch: 2020.0,
        };

        let ecef = Vector3::new(5000000.0, 1000000.0, 3000000.0);
        let transformed = params.transform(ecef, 2020.0); // dt = 0

        let mas2rad = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);
        let rx = 500.0 * mas2rad;
        let ry = 600.0 * mas2rad;
        let rz = 700.0 * mas2rad;
        let scale = 1.0 + 5.0 * 1e-9;

        let exp_x = 1.5 + scale * (ecef.x - rz * ecef.y + ry * ecef.z);
        let exp_y = 2.5 + scale * (rz * ecef.x + ecef.y - rx * ecef.z);
        let exp_z = 3.5 + scale * (-ry * ecef.x + rx * ecef.y + ecef.z);

        assert!((transformed.x - exp_x).abs() < 1e-9);
        assert!((transformed.y - exp_y).abs() < 1e-9);
        assert!((transformed.z - exp_z).abs() < 1e-9);
    }

    #[test]
    fn test_helmert_serde_roundtrip() {
        use serde_json;

        let params = HelmertParams {
            tx: 0.5,
            ty: -0.3,
            tz: 0.1,
            rx: 100.0,
            ry: 200.0,
            rz: 300.0,
            s: 2.0,
            dtx: 0.01,
            dty: -0.02,
            dtz: 0.03,
            drx: 1.0,
            dry: 2.0,
            drz: 3.0,
            ds: 0.1,
            ref_epoch: 2015.0,
        };

        let json = serde_json::to_string(&params).unwrap();
        let deserialized: HelmertParams = serde_json::from_str(&json).unwrap();

        // Check all fields
        assert!((deserialized.tx - params.tx).abs() < 1e-12);
        assert!((deserialized.ty - params.ty).abs() < 1e-12);
        assert!((deserialized.tz - params.tz).abs() < 1e-12);
        assert!((deserialized.rx - params.rx).abs() < 1e-12);
        assert!((deserialized.ry - params.ry).abs() < 1e-12);
        assert!((deserialized.rz - params.rz).abs() < 1e-12);
        assert!((deserialized.s - params.s).abs() < 1e-12);
        assert!((deserialized.dtx - params.dtx).abs() < 1e-12);
        assert!((deserialized.dty - params.dty).abs() < 1e-12);
        assert!((deserialized.dtz - params.dtz).abs() < 1e-12);
        assert!((deserialized.drx - params.drx).abs() < 1e-12);
        assert!((deserialized.dry - params.dry).abs() < 1e-12);
        assert!((deserialized.drz - params.drz).abs() < 1e-12);
        assert!((deserialized.ds - params.ds).abs() < 1e-12);
        assert!((deserialized.ref_epoch - params.ref_epoch).abs() < 1e-12);
    }

    #[test]
    fn test_helmert_apply_ecef() {
        use gneiss_core::coords::{Coordinate, Datum, Frame};
        use gneiss_core::time::GpsTime;

        let params = HelmertParams {
            tx: 1.0,
            ty: 2.0,
            tz: 3.0,
            rx: 0.0,
            ry: 0.0,
            rz: 0.0,
            s: 0.0,
            dtx: 0.0,
            dty: 0.0,
            dtz: 0.0,
            drx: 0.0,
            dry: 0.0,
            drz: 0.0,
            ds: 0.0,
            ref_epoch: 2010.0,
        };

        // Use a date close to ref_epoch so dt ~ 0
        let epoch = GpsTime::from_calendar(2010, 1, 1, 0, 0, 0.0);
        let coord = Coordinate::new(
            Vector3::new(1000.0, 2000.0, 3000.0),
            Datum::ITRF2014,
            Frame::ECEF,
            epoch,
        );

        let result = params.apply(coord);

        // dt ~ 0, so just translation applies
        assert!((result.vector.x - 1001.0).abs() < 0.01);
        assert!((result.vector.y - 2002.0).abs() < 0.01);
        assert!((result.vector.z - 3003.0).abs() < 0.01);
        assert_eq!(result.datum, Datum::ITRF2014);
        assert_eq!(result.frame, Frame::ECEF);
    }

    #[test]
    fn test_helmert_apply_non_ecef_returns_unchanged() {
        use gneiss_core::coords::{Coordinate, Datum, Frame};
        use gneiss_core::time::GpsTime;

        let params = HelmertParams {
            tx: 100.0,
            ty: 200.0,
            tz: 300.0,
            rx: 0.0,
            ry: 0.0,
            rz: 0.0,
            s: 0.0,
            dtx: 0.0,
            dty: 0.0,
            dtz: 0.0,
            drx: 0.0,
            dry: 0.0,
            drz: 0.0,
            ds: 0.0,
            ref_epoch: 2000.0,
        };

        let epoch = GpsTime::from_calendar(2020, 6, 1, 12, 0, 0.0);
        let coord = Coordinate::new(
            Vector3::new(45.0, -73.0, 100.0), // LLH: lat, lon, height
            Datum::WGS84,
            Frame::LLH,
            epoch,
        );

        let result = params.apply(coord);

        // Non-ECEF frames should be returned unchanged
        assert!((result.vector.x - 45.0).abs() < 1e-12);
        assert!((result.vector.y + 73.0).abs() < 1e-12);
        assert!((result.vector.z - 100.0).abs() < 1e-12);
    }

    #[test]
    fn test_helmert_apply_enu_returns_unchanged() {
        use gneiss_core::coords::{Coordinate, Datum, Frame};
        use gneiss_core::time::GpsTime;

        let params = HelmertParams {
            tx: 100.0,
            ty: 200.0,
            tz: 300.0,
            rx: 0.0,
            ry: 0.0,
            rz: 0.0,
            s: 0.0,
            dtx: 0.0,
            dty: 0.0,
            dtz: 0.0,
            drx: 0.0,
            dry: 0.0,
            drz: 0.0,
            ds: 0.0,
            ref_epoch: 2000.0,
        };

        let epoch = GpsTime::from_calendar(2020, 6, 1, 12, 0, 0.0);
        let coord = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ENU,
            epoch,
        );

        let result = params.apply(coord);
        assert!((result.vector.x - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_helmert_clone_and_debug() {
        let params = HelmertParams {
            tx: 1.0,
            ty: 2.0,
            tz: 3.0,
            rx: 0.0,
            ry: 0.0,
            rz: 0.0,
            s: 0.0,
            dtx: 0.0,
            dty: 0.0,
            dtz: 0.0,
            drx: 0.0,
            dry: 0.0,
            drz: 0.0,
            ds: 0.0,
            ref_epoch: 2000.0,
        };

        let cloned = params.clone();
        assert!((cloned.tx - 1.0).abs() < 1e-12);

        let debug_str = format!("{:?}", params);
        assert!(debug_str.contains("tx: 1.0"));
    }
}
