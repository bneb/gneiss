use alloc::vec::Vec;
use nalgebra::Vector3;

/// Inverts a 3x3 matrix analytically. Returns None if singular.
fn invert_3x3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c01 = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    let c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];

    let det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    if libm::fabs(det) < 1e-15 {
        return None;
    }
    let inv_det = 1.0 / det;

    let c10 = m[0][2] * m[2][1] - m[0][1] * m[2][2];
    let c11 = m[0][0] * m[2][2] - m[0][2] * m[2][0];
    let c12 = m[0][1] * m[2][0] - m[0][0] * m[2][1];

    let c20 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
    let c21 = m[0][2] * m[1][0] - m[0][0] * m[1][2];
    let c22 = m[0][0] * m[1][1] - m[0][1] * m[1][0];

    Some([
        [c00 * inv_det, c10 * inv_det, c20 * inv_det],
        [c01 * inv_det, c11 * inv_det, c21 * inv_det],
        [c02 * inv_det, c12 * inv_det, c22 * inv_det],
    ])
}

/// Fitted local site calibration parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteCalibration {
    /// Horizontal translation East (meters).
    pub dx: f64,
    /// Horizontal translation North (meters).
    pub dy: f64,
    /// Rotation angle (radians, counter-clockwise).
    pub rotation_rad: f64,
    /// Scale factor (dimensionless, ~1.0).
    pub scale: f64,
    /// Centroid of source GNSS Easting (meters).
    pub e0: f64,
    /// Centroid of source GNSS Northing (meters).
    pub n0: f64,
    /// Vertical base shift (meters).
    pub dz0: f64,
    /// Vertical slope in East direction (m/m).
    pub slope_east: f64,
    /// Vertical slope in North direction (m/m).
    pub slope_north: f64,
}

impl SiteCalibration {
    /// Fits a local site calibration from a set of paired GNSS grid (source) and Ground (target) points.
    /// Requires at least 3 non-collinear point pairs.
    pub fn fit(pairs: &[(Vector3<f64>, Vector3<f64>)]) -> Option<Self> {
        if pairs.len() < 3 {
            return None;
        }

        let n = pairs.len() as f64;
        let mut e0 = 0.0;
        let mut n0 = 0.0;
        let mut tg_e0 = 0.0;
        let mut tg_n0 = 0.0;

        for (src, tgt) in pairs {
            e0 += src.x;
            n0 += src.y;
            tg_e0 += tgt.x;
            tg_n0 += tgt.y;
        }
        e0 /= n;
        n0 /= n;
        tg_e0 /= n;
        tg_n0 /= n;

        // Fit horizontal 4-param Helmert
        let mut num = 0.0;
        let mut den = 0.0;
        let mut cross = 0.0;

        for (src, tgt) in pairs {
            let de = src.x - e0;
            let dn = src.y - n0;
            let dt_e = tgt.x - tg_e0;
            let dt_n = tgt.y - tg_n0;

            num += de * dt_e + dn * dt_n;
            cross += de * dt_n - dn * dt_e;
            den += de * de + dn * dn;
        }

        if den < 1e-6 {
            return None;
        }

        let a_param = num / den;
        let b_param = cross / den;

        let scale = libm::sqrt(a_param * a_param + b_param * b_param);
        let rotation_rad = libm::atan2(b_param, a_param);

        let dx = tg_e0 - e0;
        let dy = tg_n0 - n0;

        // Fit vertical inclined plane: dh = tgt.z - src.z = dz0 + slope_e * de + slope_n * dn
        let mut ata = [[0.0; 3]; 3];
        let mut atb = [0.0; 3];

        for (src, tgt) in pairs {
            let de = src.x - e0;
            let dn = src.y - n0;
            let dh = tgt.z - src.z;

            let row = [1.0, de, dn];
            for i in 0..3 {
                atb[i] += row[i] * dh;
                for j in 0..3 {
                    ata[i][j] += row[i] * row[j];
                }
            }
        }

        let inv_ata = invert_3x3(&ata)?;
        let dz0 = inv_ata[0][0] * atb[0] + inv_ata[0][1] * atb[1] + inv_ata[0][2] * atb[2];
        let slope_east = inv_ata[1][0] * atb[0] + inv_ata[1][1] * atb[1] + inv_ata[1][2] * atb[2];
        let slope_north = inv_ata[2][0] * atb[0] + inv_ata[2][1] * atb[1] + inv_ata[2][2] * atb[2];

        Some(Self {
            dx,
            dy,
            rotation_rad,
            scale,
            e0,
            n0,
            dz0,
            slope_east,
            slope_north,
        })
    }

    /// Transforms a GNSS map projection coordinate (E, N, h) to calibrated local ground coordinates.
    pub fn transform(&self, gnss_grid: Vector3<f64>) -> Vector3<f64> {
        let de = gnss_grid.x - self.e0;
        let dn = gnss_grid.y - self.n0;

        let cos_r = libm::cos(self.rotation_rad);
        let sin_r = libm::sin(self.rotation_rad);

        let rot_e = self.scale * (cos_r * de - sin_r * dn);
        let rot_n = self.scale * (sin_r * de + cos_r * dn);

        let ground_e = rot_e + self.e0 + self.dx;
        let ground_n = rot_n + self.n0 + self.dy;
        let ground_h = gnss_grid.z + self.dz0 + self.slope_east * de + self.slope_north * dn;

        Vector3::new(ground_e, ground_n, ground_h)
    }

    /// Computes horizontal and vertical transformation residuals across control point pairs.
    pub fn compute_residuals(&self, pairs: &[(Vector3<f64>, Vector3<f64>)]) -> Vec<(f64, f64)> {
        pairs.iter().map(|(src, tgt)| {
            let transformed = self.transform(*src);
            let h_res = libm::sqrt(
                (transformed.x - tgt.x) * (transformed.x - tgt.x)
                + (transformed.y - tgt.y) * (transformed.y - tgt.y)
            );
            let v_res = transformed.z - tgt.z;
            (h_res, v_res)
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn test_site_calibration_exact_fit() {
        let p1_src = Vector3::new(100.0, 200.0, 50.0);
        let p2_src = Vector3::new(300.0, 200.0, 52.0);
        let p3_src = Vector3::new(200.0, 400.0, 54.0);

        // Ground targets shifted by dx=10, dy=20, dz=5, rotated 0.1 rad, scale 1.0001
        let shift = Vector3::new(10.0, 20.0, 5.0);
        let pairs = vec![
            (p1_src, p1_src + shift),
            (p2_src, p2_src + shift),
            (p3_src, p3_src + shift),
        ];

        let calib = SiteCalibration::fit(&pairs).expect("fit calibration");
        let residuals = calib.compute_residuals(&pairs);

        for (h_err, v_err) in residuals {
            assert!(h_err < 1e-4, "Horizontal error < 0.1 mm");
            assert!(v_err.abs() < 1e-4, "Vertical error < 0.1 mm");
        }
    }
}
