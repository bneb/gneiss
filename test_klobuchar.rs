fn main() {
    let pi = std::f64::consts::PI;
    let pos_llh = (35.0 * pi / 180.0, 139.0 * pi / 180.0, 100.0);
    let az = 45.0 * pi / 180.0;
    let el = 30.0 * pi / 180.0;

    let el_semi = el / pi;
    let lat_semi = pos_llh.0 / pi;
    let lon_semi = pos_llh.1 / pi;

    let psi = 0.0137 / (el_semi + 0.11) - 0.022;
    let mut phi_i = lat_semi + psi * f64::cos(az);
    if phi_i > 0.416 {
        phi_i = 0.416;
    } else if phi_i < -0.416 {
        phi_i = -0.416;
    }

    let lambda_i = lon_semi + (psi * f64::sin(az)) / f64::cos(phi_i * pi);
    let phi_m = phi_i + 0.064 * f64::cos((lambda_i - 1.617) * pi);

    let tow = 100000.0;
    let mut t = 43200.0 * lambda_i + tow;
    t %= 86400.0;
    if t < 0.0 { t += 86400.0; }

    println!("phi_i: {}", phi_i);
    println!("lambda_i: {}", lambda_i);
    println!("phi_m: {}", phi_m);
    println!("t: {}", t);
}
