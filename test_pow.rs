fn main() {
    let x = -0.1_f64;
    println!("powf: {}", x.powf(3.0));
    println!("libm::pow: {}", libm::pow(x, 3.0));
}
