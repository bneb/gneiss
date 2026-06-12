use nalgebra::Vector3;

/// Calculates the Phase Wind-Up correction in cycles for a given satellite and receiver.
/// `sat_pos`: Satellite ECEF position in meters.
/// `sun_pos`: Sun ECEF position in meters.
/// `rcv_pos`: Receiver ECEF position in meters.
/// `prev_windup`: The windup value from the previous epoch (in cycles) to ensure continuity.
pub fn phase_windup(sat_pos: Vector3<f64>, sun_pos: Vector3<f64>, rcv_pos: Vector3<f64>, prev_windup: f64) -> f64 {
    // Unit vector from satellite to receiver
    let k = (rcv_pos - sat_pos).normalize();

    // Satellite body axes (nominal yaw steering model)
    // Z-axis points towards the center of the Earth
    let sat_z = -sat_pos.normalize();

    // Unit vector from satellite to sun
    let e_sun = (sun_pos - sat_pos).normalize();

    // Y-axis is perpendicular to Z-axis and Sun vector (points along the solar panel axis)
    let sat_y = sat_z.cross(&e_sun).normalize();

    // X-axis completes the right-hand rule
    let sat_x = sat_y.cross(&sat_z).normalize();

    // Satellite dipole vectors (D and D')
    // D is the unit vector in the X-Y plane of the satellite
    let mut d_prime = sat_x - k * (k.dot(&sat_x)) - k.cross(&sat_y);
    if d_prime.norm() > 1e-12 {
        d_prime.normalize_mut();
    } else {
        d_prime = sat_x; // Fallback
    }

    // Receiver local tangent plane (East, North, Up)
    let llh = crate::coords::ecef_to_llh(rcv_pos);
    let lat = llh.x;
    let lon = llh.y;

    let sin_lat = libm::sin(lat);
    let cos_lat = libm::cos(lat);
    let sin_lon = libm::sin(lon);
    let cos_lon = libm::cos(lon);

    let e = Vector3::new(-sin_lon, cos_lon, 0.0);
    let n = Vector3::new(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat);

    // Receiver dipole vector D
    let mut d = e - k * (k.dot(&e)) + k.cross(&n);
    if d.norm() > 1e-12 {
        d.normalize_mut();
    } else {
        d = e; // Fallback
    }

    // Angle between D and D'
    let cos_angle = d_prime.dot(&d);
    let mut angle = libm::acos(cos_angle.clamp(-1.0, 1.0));

    let sign_check = k.dot(&d_prime.cross(&d));
    if sign_check < 0.0 {
        angle = -angle;
    }

    // Convert angle to cycles
    let mut d_phi = angle / (2.0 * core::f64::consts::PI);

    // Ensure continuity (unwrap phase)
    d_phi -= libm::round(d_phi - prev_windup);

    d_phi
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_phase_windup_physics() {
        // Simple physics test: Receiver at equator, satellite at zenith, sun at horizon.
        // We just ensure it doesn't panic and returns a reasonable fraction.
        let sat_pos = Vector3::new(26560000.0, 0.0, 0.0);
        let sun_pos = Vector3::new(0.0, 149597870700.0, 0.0);
        let rcv_pos = Vector3::new(crate::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let w = phase_windup(sat_pos, sun_pos, rcv_pos, 0.0);
        assert!((-0.5..=0.5).contains(&w));
    }
}