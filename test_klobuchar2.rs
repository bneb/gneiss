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
    let (lat, lon, alt) = ecef_to_llh(-2689639.506, -4290438.636, 3865050.956);
    println!("lat: {} rad ({} deg), lon: {} rad ({} deg), alt: {} m", lat, lat * 180.0 / pi, lon, lon * 180.0 / pi, alt);

    let az = 45.0 * pi / 180.0;
    let el = 15.0 * pi / 180.0;

    let el_semi = el / pi;
    let lat_semi = lat / pi;
    let lon_semi = lon / pi;

    let psi = 0.0137 / (el_semi + 0.11) - 0.022;
    let mut phi_i = lat_semi + psi * f64::cos(az);
    if phi_i > 0.416 { phi_i = 0.416; } else if phi_i < -0.416 { phi_i = -0.416; }

    let lambda_i = lon_semi + (psi * f64::sin(az)) / f64::cos(phi_i * pi);
    let phi_m = phi_i + 0.064 * f64::cos((lambda_i - 1.617) * pi);

    let tow = 113055.0; // from log
    let mut t = 43200.0 * lambda_i + tow;
    t %= 86400.0;
    if t < 0.0 { t += 86400.0; }

    let alpha = [1.1176e-08, 2.2352e-08, -5.9605e-08, -1.1921e-07];
    let beta = [90112.0, 114688.0, -196608.0, -655360.0];

    let mut a = alpha[0] + alpha[1] * phi_m + alpha[2] * phi_m * phi_m + alpha[3] * phi_m * phi_m * phi_m;
    if a < 0.0 { a = 0.0; }

    let mut p = beta[0] + beta[1] * phi_m + beta[2] * phi_m * phi_m + beta[3] * phi_m * phi_m * phi_m;
    if p < 72000.0 { p = 72000.0; }

    let x = 2.0 * pi * (t - 50400.0) / p;
    let f = 1.0 + 16.0 * f64::powf(0.53 - el_semi, 3.0);

    let delay = if f64::abs(x) < 1.57 {
        5e-9 + a * (1.0 - x * x / 2.0 + x * x * x * x / 24.0)
    } else {
        5e-9
    };

    println!("delay: {} m", delay * f * 299792458.0);
}
