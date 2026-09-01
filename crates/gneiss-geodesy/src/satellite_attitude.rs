//! Satellite attitude and Antenna Phase Center (PCO) ECEF projection.
//!
//! Implements the nominal GNSS yaw-steering attitude model (Montenbruck & Gill 2000,
//! Bar-Sever 1996, IGS white paper on satellite body frames).

use nalgebra::{Matrix3, Vector3};

#[inline]
fn vec_norm(v: &Vector3<f64>) -> f64 {
    libm::sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

/// Computes the satellite body-to-ECEF rotation matrix \mathbf{R}_{\text{sat}} = [\hat{\mathbf{e}}_x, \hat{\mathbf{e}}_y, \hat{\mathbf{e}}_z].
///
/// Coordinate frame definitions:
/// - \hat{\mathbf{e}}_z: Points from satellite center of mass to Earth center (nadir).
/// - \hat{\mathbf{e}}_y: Points along solar panel rotation axis: \hat{\mathbf{e}}_z \times \hat{\mathbf{e}}_{\text{sun}} / \|\dots\|.
/// - \hat{\mathbf{e}}_x: Completes right-handed orthogonal system: \hat{\mathbf{e}}_y \times \hat{\mathbf{e}}_z.
pub fn nominal_satellite_attitude_matrix(
    sat_pos_ecef: &Vector3<f64>,
    sun_pos_ecef: &Vector3<f64>,
) -> Option<Matrix3<f64>> {
    let r_sat = vec_norm(sat_pos_ecef);
    let r_sun = vec_norm(sun_pos_ecef);
    if r_sat < 1.0 || r_sun < 1.0 {
        return None;
    }

    // e_z points nadir (towards Earth center)
    let e_z = -sat_pos_ecef / r_sat;

    // Unit vector towards the Sun
    let e_sun = sun_pos_ecef / r_sun;

    // Solar panel axis (e_y)
    let e_y_unnorm = e_z.cross(&e_sun);
    let e_y_len = vec_norm(&e_y_unnorm);
    if e_y_len < 1e-6 {
        return None; // Collinear / eclipse boundary
    }
    let e_y = e_y_unnorm / e_y_len;

    // e_x = e_y x e_z
    let e_x = e_y.cross(&e_z);

    Some(Matrix3::from_columns(&[e_x, e_y, e_z]))
}

/// Projects a satellite body-frame Phase Center Offset (PCO) into ECEF frame coordinates (meters).
pub fn project_satellite_pco_to_ecef(
    sat_pos_ecef: &Vector3<f64>,
    sun_pos_ecef: &Vector3<f64>,
    body_pco: &Vector3<f64>,
) -> Vector3<f64> {
    if let Some(r_sat) = nominal_satellite_attitude_matrix(sat_pos_ecef, sun_pos_ecef) {
        r_sat * body_pco
    } else {
        Vector3::zeros()
    }
}

/// Determines if a satellite is in the Earth's shadow cone (umbra/penumbra).
///
/// Uses the standard cylindrical shadow model (Montenbruck & Gill 2000, Section 3.4)
/// with an Earth equatorial radius of 6378.137 km plus 40 km atmospheric absorption buffer.
pub fn is_satellite_eclipsed(sat_pos_ecef: &Vector3<f64>, sun_pos_ecef: &Vector3<f64>) -> bool {
    let r_sun = vec_norm(sun_pos_ecef);
    if r_sun < 1.0 {
        return false;
    }
    let e_sun = sun_pos_ecef / r_sun;
    let proj = sat_pos_ecef.dot(&e_sun);

    // If proj >= 0, satellite is on the day-side of Earth (illuminated)
    if proj >= 0.0 {
        return false;
    }

    // Perpendicular distance to Sun-Earth axis
    let perp_vec = sat_pos_ecef - proj * e_sun;
    let d_perp = vec_norm(&perp_vec);

    let r_shadow = 6_378_137.0 + 40_000.0;
    d_perp < r_shadow
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nominal_satellite_attitude_orthonormal() {
        let sat = Vector3::new(0.0, 26_560_000.0, 0.0);
        let sun = Vector3::new(1.495e11, 0.0, 0.0);

        let rot = nominal_satellite_attitude_matrix(&sat, &sun).unwrap();

        // Must be orthonormal (R^T R = I) and det = +1
        let diff = rot.transpose() * rot - Matrix3::identity();
        let frob_norm = libm::sqrt(diff.iter().map(|&x| x * x).sum::<f64>());
        assert!(frob_norm < 1e-10);

        let det = rot[(0, 0)] * (rot[(1, 1)] * rot[(2, 2)] - rot[(1, 2)] * rot[(2, 1)])
            - rot[(0, 1)] * (rot[(1, 0)] * rot[(2, 2)] - rot[(1, 2)] * rot[(2, 0)])
            + rot[(0, 2)] * (rot[(1, 0)] * rot[(2, 1)] - rot[(1, 1)] * rot[(2, 0)]);
        assert!((det - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_project_satellite_pco_textbook_vector() {
        // Sat at [0, 26560km, 0], Sun at [1.5e11, 0, 0]
        let sat = Vector3::new(0.0, 26_560_000.0, 0.0);
        let sun = Vector3::new(1.495e11, 0.0, 0.0);

        // Body PCO: dx = 0.1m, dy = 0.2m, dz = 1.5m
        let body_pco = Vector3::new(0.1, 0.2, 1.5);
        let ecef_pco = project_satellite_pco_to_ecef(&sat, &sun, &body_pco);

        // Analytical solution:
        // e_z = (0, -1, 0)
        // e_sun = (1, 0, 0)
        // e_y = e_z x e_sun = (0, 0, 1)
        // e_x = e_y x e_z = (1, 0, 0)
        // PCO_ecef = 0.1*(1,0,0) + 0.2*(0,0,1) + 1.5*(0,-1,0) = [0.1, -1.5, 0.2]
        assert!((ecef_pco.x - 0.1).abs() < 1e-6);
        assert!((ecef_pco.y - (-1.5)).abs() < 1e-6);
        assert!((ecef_pco.z - 0.2).abs() < 1e-6);
    }

    #[test]
    fn test_is_satellite_eclipsed_geometry() {
        let sun = Vector3::new(1.495e11, 0.0, 0.0);

        // Sat 1: Day side (+x) -> NOT eclipsed
        let sat_day = Vector3::new(26_560_000.0, 0.0, 0.0);
        assert!(!is_satellite_eclipsed(&sat_day, &sun));

        // Sat 2: Night side directly behind Earth (-x, y=0, z=0) -> ECLIPSED
        let sat_shadow = Vector3::new(-26_560_000.0, 0.0, 0.0);
        assert!(is_satellite_eclipsed(&sat_shadow, &sun));

        // Sat 3: Night side but high y offset (outside shadow cylinder) -> NOT eclipsed
        let sat_outside = Vector3::new(-26_560_000.0, 10_000_000.0, 0.0);
        assert!(!is_satellite_eclipsed(&sat_outside, &sun));
    }
}
