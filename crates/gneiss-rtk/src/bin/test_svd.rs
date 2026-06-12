use nalgebra::DMatrix;
fn main() {
    let mut cov = DMatrix::<f64>::zeros(2, 2);
    cov[(0, 0)] = 1e-12; // Position
    cov[(1, 1)] = 1e6; // Clock bias
    let svd = cov.svd(true, true);
    println!("Singular values: {:?}", svd.singular_values);
    let inv = svd.pseudo_inverse(1e-6).unwrap();
    println!("PriorInfo pos: {}", inv[(0, 0)]);
}
