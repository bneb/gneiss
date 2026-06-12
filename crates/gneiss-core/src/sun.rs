use nalgebra::Vector3;
use crate::time::GpsTime;

/// Calculates the approximate position of the Sun in ECEF coordinates.
/// Accuracy is around 0.1 deg, which is sufficient for GNSS attitude and phase wind-up modeling.
pub fn sun_position_ecef(t: GpsTime) -> Vector3<f64> {
    // Julian centuries since J2000.0
    // J2000.0 is 2000-01-01 12:00:00 UTC = GPS Time 1000, 0.0s? No.
    // GPS Time started Jan 6 1980.
    // J2000.0 is exactly 7292.5 days after GPS epoch.
    // Weeks: 7292.5 / 7 = 1041.7857...
    let days_since_gps_epoch = (t.week as f64) * 7.0 + (t.tow / 86400.0);
    let d = days_since_gps_epoch - 7292.5; // Days since J2000.0
    
    let t_jc = d / 36525.0; // Julian centuries
    
    // Mean anomaly of the Sun
    let m = (357.52911 + 35999.05029 * t_jc).to_radians();
    
    // Mean longitude of the Sun
    let l_0 = (280.46646 + 36000.76983 * t_jc).to_radians();
    
    // Ecliptic longitude
    let lambda = l_0 + (1.914602 - 0.004817 * t_jc) * libm::sin(m) + (0.019993 - 0.000101 * t_jc) * libm::sin(2.0 * m);
    
    // Obliquity of the ecliptic
    let epsilon = (23.439291 - 0.0130042 * t_jc).to_radians();
    
    // Distance to the Sun in Astronomical Units (AU)
    let r_au = 1.000140612 - 0.016708617 * libm::cos(m) - 0.000139589 * libm::cos(2.0 * m);
    
    // Convert AU to meters
    let r_meters = r_au * 149597870700.0;
    
    // Sun position in ECI (Earth-Centered Inertial)
    let x_eci = r_meters * libm::cos(lambda);
    let y_eci = r_meters * libm::sin(lambda) * libm::cos(epsilon);
    let z_eci = r_meters * libm::sin(lambda) * libm::sin(epsilon);
    
    // Greenwich Mean Sidereal Time (GMST) to rotate ECI to ECEF
    // Simplified GMST approximation based on GPS time (this ignores precise UT1-UTC, but is enough for attitude)
    let gmst = (4.8949612128230587e-5 * d * 86400.0 + 1.7533685592333) % (2.0 * core::f64::consts::PI);
    
    let cos_gmst = libm::cos(gmst);
    let sin_gmst = libm::sin(gmst);
    
    let x_ecef = x_eci * cos_gmst + y_eci * sin_gmst;
    let y_ecef = -x_eci * sin_gmst + y_eci * cos_gmst;
    let z_ecef = z_eci;
    
    Vector3::new(x_ecef, y_ecef, z_ecef)
}

/// Calculates the approximate position of the Moon in ECEF coordinates.
pub fn moon_position_ecef(t: GpsTime) -> Vector3<f64> {
    let days_since_gps_epoch = (t.week as f64) * 7.0 + (t.tow / 86400.0);
    let d = days_since_gps_epoch - 7292.5; // Days since J2000.0
    let t_jc = d / 36525.0; // Julian centuries

    // Astronomical arguments
    let fc = [
        [ 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470 ],
        [ 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149 ],
        [  93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417 ],
        [ 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169 ],
        [ 125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939 ]
    ];
    let mut f = [0.0; 5];
    let tt = [t_jc, t_jc*t_jc, t_jc*t_jc*t_jc, t_jc*t_jc*t_jc*t_jc];
    for i in 0..5 {
        f[i] = fc[i][0] * 3600.0;
        for j in 0..4 {
            f[i] += fc[i][j+1] * tt[j];
        }
        f[i] = (f[i] * (core::f64::consts::PI / (180.0 * 3600.0))) % (2.0 * core::f64::consts::PI);
    }

    let lm = 218.32 + 481267.883 * t_jc + 6.29 * libm::sin(f[0]) - 1.27 * libm::sin(f[0] - 2.0 * f[3])
           + 0.66 * libm::sin(2.0 * f[3]) + 0.21 * libm::sin(2.0 * f[0]) - 0.19 * libm::sin(f[1]) - 0.11 * libm::sin(2.0 * f[2]);
    let pm = 5.13 * libm::sin(f[2]) + 0.28 * libm::sin(f[0] + f[2]) - 0.28 * libm::sin(f[2] - f[0])
           - 0.17 * libm::sin(f[2] - 2.0 * f[3]);
    let rm = crate::constants::WGS84_SEMI_MAJOR_AXIS_M / libm::sin((0.9508 + 0.0518 * libm::cos(f[0]) + 0.0095 * libm::cos(f[0] - 2.0 * f[3])
           + 0.0078 * libm::cos(2.0 * f[3]) + 0.0028 * libm::cos(2.0 * f[0])) * core::f64::consts::PI / 180.0);

    let sinl = libm::sin(lm * core::f64::consts::PI / 180.0);
    let cosl = libm::cos(lm * core::f64::consts::PI / 180.0);
    let sinp = libm::sin(pm * core::f64::consts::PI / 180.0);
    let cosp = libm::cos(pm * core::f64::consts::PI / 180.0);
    
    // Obliquity of the ecliptic
    let eps = (23.439291 - 0.0130042 * t_jc).to_radians();
    let sine = libm::sin(eps);
    let cose = libm::cos(eps);

    let x_eci = rm * cosp * cosl;
    let y_eci = rm * (cose * cosp * sinl - sine * sinp);
    let z_eci = rm * (sine * cosp * sinl + cose * sinp);

    let gmst = (4.8949612128230587e-5 * d * 86400.0 + 1.7533685592333) % (2.0 * core::f64::consts::PI);
    let cos_gmst = libm::cos(gmst);
    let sin_gmst = libm::sin(gmst);
    
    let x_ecef = x_eci * cos_gmst + y_eci * sin_gmst;
    let y_ecef = -x_eci * sin_gmst + y_eci * cos_gmst;
    let z_ecef = z_eci;
    
    Vector3::new(x_ecef, y_ecef, z_ecef)
}