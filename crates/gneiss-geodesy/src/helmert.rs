use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

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
    /// Transforms an ECEF vector from the source frame to the target frame at a specific observation epoch.
    pub fn transform(&self, ecef: Vector3<f64>, obs_epoch: f64) -> Vector3<f64> {
        let dt = obs_epoch - self.ref_epoch;

        let tx = self.tx + self.dtx * dt;
        let ty = self.ty + self.dty * dt;
        let tz = self.tz + self.dtz * dt;

        let rx_mas = self.rx + self.drx * dt;
        let ry_mas = self.ry + self.dry * dt;
        let rz_mas = self.rz + self.drz * dt;

        let s_ppb = self.s + self.ds * dt;

        // Conversions
        let mas2rad = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);
        let rx = rx_mas * mas2rad;
        let ry = ry_mas * mas2rad;
        let rz = rz_mas * mas2rad;
        let scale = 1.0 + s_ppb * 1e-9;

        // Apply Helmert Transformation (position convention)
        // X' = T + (1 + s) * R * X
        // R is the rotation matrix for small angles:
        // |   1   -rz   ry |
        // |  rz     1  -rx |
        // | -ry    rx    1 |

        let x_new = tx + scale * (ecef.x - rz * ecef.y + ry * ecef.z);
        let y_new = ty + scale * (rz * ecef.x + ecef.y - rx * ecef.z);
        let z_new = tz + scale * (-ry * ecef.x + rx * ecef.y + ecef.z);

        Vector3::new(x_new, y_new, z_new)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_helmert_itrf2014_to_itrf2020() {
        // Transformation parameters from ITRF2014 to ITRF2020 at epoch 2015.0
        // Provided by IERS
        // tx, ty, tz in mm -> convert to m
        let params = HelmertParams {
            tx: -0.0014, ty: -0.0012, tz:  0.0012,
            rx:  0.0,    ry:  0.0,    rz:  0.0,
            s:   0.0,
            
            dtx:  0.0,    dty: -0.0001, dtz:  0.0002,
            drx:  0.0,    dry:  0.0,    drz:  0.0,
            ds:   0.0,
            
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

        assert!((transformed.x - expected_x).abs() < 1e-4);
        assert!((transformed.y - expected_y).abs() < 1e-4);
        assert!((transformed.z - expected_z).abs() < 1e-4);
    }

    #[test]
    fn test_helmert_with_rotations() {
        let params = HelmertParams {
            tx: 1.0, ty: 2.0, tz: 3.0,
            rx: 1000.0, // 1000 mas = 1 arcsec
            ry: 2000.0, // 2 arcsec
            rz: 3000.0, // 3 arcsec
            s: 10.0,    // 10 ppb
            dtx: 0.0, dty: 0.0, dtz: 0.0,
            drx: 0.0, dry: 0.0, drz: 0.0, ds: 0.0,
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
}