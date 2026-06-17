fn ecef_to_llh(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let a = 6378137.0;
    let b = 6356752.314245179;
    let e_sq = 1.0 - (b * b) / (a * a);
    let p = (x * x + y * y).sqrt();
    let lon = y.atan2(x);
    let mut lat = (z / (p * (1.0 - e_sq))).atan();
    let mut alt = 0.0;
    for _ in 0..5 {
        let n = a / (1.0 - e_sq * lat.sin().powi(2)).sqrt();
        alt = p / lat.cos() - n;
        lat = (z / (p * (1.0 - e_sq * n / (n + alt)))).atan();
    }
    (lat, lon, alt)
}

fn main() {
    let pi = std::f64::consts::PI;
    let (lat, lon, _) = ecef_to_llh(-2689639.506, -4290438.636, 3865050.956);

    let az = 45.0 * pi / 180.0;
    let el = 15.0 * pi / 180.0;
    let tow = 113055.0;

    let alpha = [1.1176e-08, 2.2352e-08, -5.9605e-08, -1.1921e-07];
    let beta = [90112.0, 114688.0, -196608.0, -655360.0];

    // OLD CODE
    let f_old = 1.0 + 16.0 * f64::powf(0.53 - el / pi, 3.0);
    let phi_m_old = lat / pi + 0.064 * f64::cos(lon - 1.617);
    
    let mut t_old = 43200.0 * phi_m_old + tow;
    t_old %= 86400.0;
    if t_old < 0.0 { t_old += 86400.0; }

    let mut a_old = alpha[0] + alpha[1] * phi_m_old + alpha[2] * phi_m_old * phi_m_old + alpha[3] * phi_m_old * phi_m_old * phi_m_old;
    if a_old < 0.0 { a_old = 0.0; }

    let mut p_old = beta[0] + beta[1] * phi_m_old + beta[2] * phi_m_old * phi_m_old + beta[3] * phi_m_old * phi_m_old * phi_m_old;
    if p_old < 72000.0 { p_old = 72000.0; }

    let x_old = 2.0 * pi * (t_old - 50400.0) / p_old;
    let delay_old = if f64::abs(x_old) < 1.57 { 5e-9 + a_old * (1.0 - x_old * x_old / 2.0 + x_old * x_old * x_old * x_old / 24.0) } else { 5e-9 };

    // NEW CODE
    let el_semi = el / pi;
    let lat_semi = lat / pi;
    let lon_semi = lon / pi;

    let psi = 0.0137 / (el_semi + 0.11) - 0.022;
    let mut phi_i = lat_semi + psi * f64::cos(az);
    if phi_i > 0.416 { phi_i = 0.416; } else if phi_i < -0.416 { phi_i = -0.416; }

    let lambda_i = lon_semi + (psi * f64::sin(az)) / f64::cos(phi_i * pi);
    let phi_m_new = phi_i + 0.064 * f64::cos((lambda_i - 1.617) * pi);

    let mut t_new = 43200.0 * lambda_i + tow;
    t_new %= 86400.0;
    if t_new < 0.0 { t_new += 86400.0; }

    let mut a_new = alpha[0] + alpha[1] * phi_m_new + alpha[2] * phi_m_new * phi_m_new + alpha[3] * phi_m_new * phi_m_new * phi_m_new;
    if a_new < 0.0 { a_new = 0.0; }

    let mut p_new = beta[0] + beta[1] * phi_m_new + beta[2] * phi_m_new * phi_m_new + beta[3] * phi_m_new * phi_m_new * phi_m_new;
    if p_new < 72000.0 { p_new = 72000.0; }

    let x_new = 2.0 * pi * (t_new - 50400.0) / p_new;
    let f_new = 1.0 + 16.0 * f64::powf(0.53 - el_semi, 3.0);
    let delay_new = if f64::abs(x_new) < 1.57 { 5e-9 + a_new * (1.0 - x_new * x_new / 2.0 + x_new * x_new * x_new * x_new / 24.0) } else { 5e-9 };

    println!("OLD Delay: {} m", delay_old * f_old * 299792458.0);
    println!("NEW Delay: {} m", delay_new * f_new * 299792458.0);
}
