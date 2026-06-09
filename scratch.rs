fn main() {
    let mut m = nalgebra::DMatrix::<f64>::zeros(3, 3);
    m[(0,0)] = f64::NAN;
    println!("Trying to inverse...");
    let _ = m.pseudo_inverse(1e-12);
    println!("Done");
}
