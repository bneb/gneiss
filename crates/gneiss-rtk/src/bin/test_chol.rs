fn main() {
    use nalgebra::DMatrix;
    let mut m = DMatrix::<f64>::zeros(2, 2);
    m[(0,0)] = 2.0; m[(0,1)] = 1.0;
    m[(1,0)] = 1.0; m[(1,1)] = 2.0;
    let chol = m.clone().cholesky().unwrap();
    let inv = chol.inverse();
    println!("Inverse:\n{}", inv);
}
