use crate::time::GpsTime;
use nalgebra::Vector3;

/// Calculates Julian centuries and days since J2000.0
fn julian_centuries(t: GpsTime) -> (f64, f64) {
    let days_since_gps = (t.week as f64) * 7.0 + (t.tow / crate::constants::SECONDS_PER_DAY);
    let d = days_since_gps - crate::constants::DAYS_GPS_TO_J2000;
    (d / crate::constants::DAYS_PER_JULIAN_CENTURY, d)
}

/// Rotates an ECI vector to ECEF given days since J2000.0
fn eci_to_ecef(eci: Vector3<f64>, d: f64) -> Vector3<f64> {
    let gmst = ((18.697374558 + 24.06570982441908 * d) * (core::f64::consts::PI / 12.0))
        % (2.0 * core::f64::consts::PI);
    let cos_gmst = libm::cos(gmst);
    let sin_gmst = libm::sin(gmst);
    Vector3::new(
        eci.x * cos_gmst + eci.y * sin_gmst,
        -eci.x * sin_gmst + eci.y * cos_gmst,
        eci.z,
    )
}

/// Calculates the approximate position of the Sun in ECI coordinates.
fn sun_position_eci(t_jc: f64) -> Vector3<f64> {
    let m = (357.52911 + 35999.05029 * t_jc).to_radians();
    let l_0 = (280.46646 + 36000.76983 * t_jc).to_radians();
    let lambda = l_0
        + (1.914602 - 0.004817 * t_jc) * libm::sin(m)
        + (0.019993 - 0.000101 * t_jc) * libm::sin(2.0 * m);
    let epsilon = (23.439291 - 0.0130042 * t_jc).to_radians();
    let r_au = 1.000140612 - 0.016708617 * libm::cos(m) - 0.000139589 * libm::cos(2.0 * m);
    let r_meters = r_au * crate::constants::ASTRONOMICAL_UNIT_M;
    Vector3::new(
        r_meters * libm::cos(lambda),
        r_meters * libm::sin(lambda) * libm::cos(epsilon),
        r_meters * libm::sin(lambda) * libm::sin(epsilon),
    )
}

/// Calculates the approximate position of the Sun in ECEF coordinates.
/// Accuracy is around 0.1 deg, which is sufficient for GNSS attitude and phase wind-up modeling.
pub fn sun_position_ecef(t: GpsTime) -> Vector3<f64> {
    let (t_jc, d) = julian_centuries(t);
    eci_to_ecef(sun_position_eci(t_jc), d)
}

/// Calculates the approximate position of the Moon in ECI coordinates.
fn moon_position_eci(t_jc: f64) -> Vector3<f64> {
    let fc = [
        [
            134.96340251,
            1717915923.2178,
            31.8792,
            0.051635,
            -0.00024470,
        ],
        [357.52910918, 129596581.0481, -0.5532, 0.000136, -0.00001149],
        [
            93.27209062,
            1739527262.8478,
            -12.7512,
            -0.001037,
            0.00000417,
        ],
        [
            297.85019547,
            1602961601.2090,
            -6.3706,
            0.006593,
            -0.00003169,
        ],
        [125.04455501, -6962890.2665, 7.4722, 0.007702, -0.00005939],
    ];
    let mut f = [0.0; 5];
    let tt = [
        t_jc,
        t_jc * t_jc,
        t_jc * t_jc * t_jc,
        t_jc * t_jc * t_jc * t_jc,
    ];
    for i in 0..5 {
        f[i] = fc[i][0] * 3600.0;
        for j in 0..4 {
            f[i] += fc[i][j + 1] * tt[j];
        }
        f[i] = (f[i] * (core::f64::consts::PI / (180.0 * 3600.0))) % (2.0 * core::f64::consts::PI);
    }
    let lm = 218.32 + 481267.883 * t_jc + 6.29 * libm::sin(f[0])
        - 1.27 * libm::sin(f[0] - 2.0 * f[3])
        + 0.66 * libm::sin(2.0 * f[3])
        + 0.21 * libm::sin(2.0 * f[0])
        - 0.19 * libm::sin(f[1])
        - 0.11 * libm::sin(2.0 * f[2]);
    let pm = 5.13 * libm::sin(f[2]) + 0.28 * libm::sin(f[0] + f[2])
        - 0.28 * libm::sin(f[2] - f[0])
        - 0.17 * libm::sin(f[2] - 2.0 * f[3]);
    let rm = crate::constants::WGS84_SEMI_MAJOR_AXIS_M
        / libm::sin(
            (0.9508
                + 0.0518 * libm::cos(f[0])
                + 0.0095 * libm::cos(f[0] - 2.0 * f[3])
                + 0.0078 * libm::cos(2.0 * f[3])
                + 0.0028 * libm::cos(2.0 * f[0]))
                * core::f64::consts::PI
                / 180.0,
        );

    let sinl = libm::sin(lm * core::f64::consts::PI / 180.0);
    let cosl = libm::cos(lm * core::f64::consts::PI / 180.0);
    let sinp = libm::sin(pm * core::f64::consts::PI / 180.0);
    let cosp = libm::cos(pm * core::f64::consts::PI / 180.0);
    let eps = (23.439291 - 0.0130042 * t_jc).to_radians();

    Vector3::new(
        rm * cosp * cosl,
        rm * (libm::cos(eps) * cosp * sinl - libm::sin(eps) * sinp),
        rm * (libm::sin(eps) * cosp * sinl + libm::cos(eps) * sinp),
    )
}

/// Calculates the approximate position of the Moon in ECEF coordinates.
pub fn moon_position_ecef(t: GpsTime) -> Vector3<f64> {
    let (t_jc, d) = julian_centuries(t);
    eci_to_ecef(moon_position_eci(t_jc), d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sun_position_ecef_snapshot() {
        let t1 = GpsTime::new(0, 0.0);
        let pos1 = sun_position_ecef(t1);

        let t2 = GpsTime::new(2000, 86400.0);
        let pos2 = sun_position_ecef(t2);

        assert!((pos1.x - -130361764502.90057).abs() < 1e-4);
        assert!((pos1.y - -52655996731.87905).abs() < 1e-4);
        assert!((pos1.z - -43390603069.40688).abs() < 1e-4);

        assert!((pos2.x - -10056354437.039314).abs() < 1e-4);
        assert!((pos2.y - -145030993652.9578).abs() < 1e-4);
        assert!((pos2.z - 41587067139.174736).abs() < 1e-4);
    }

    #[test]
    fn test_moon_position_ecef_snapshot() {
        let t1 = GpsTime::new(0, 0.0);
        let pos1 = moon_position_ecef(t1);

        let t2 = GpsTime::new(2000, 86400.0);
        let pos2 = moon_position_ecef(t2);

        assert!((pos1.x - -235680777.45259964).abs() < 1e-4);
        assert!((pos1.y - 286919902.5195386).abs() < 1e-4);
        assert!((pos1.z - -104315553.52824344).abs() < 1e-4);

        assert!((pos2.x - -358772193.8753879).abs() < 1e-4);
        assert!((pos2.y - 38572485.482881606).abs() < 1e-4);
        assert!((pos2.z - 78820686.42981544).abs() < 1e-4);
    }
}
