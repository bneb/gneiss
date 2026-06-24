// Quick experiment to find a matrix that triggers the fallback path
use nalgebra::dmatrix;

fn main() {
    // Try various pathological matrices
    let matrices = vec![
        ("NaN 2x2", dmatrix![f64::NAN, 0.0; 0.0, 1.0]),
        ("Inf 2x2", dmatrix![f64::INFINITY, 0.0; 0.0, 1.0]),
        ("Neg Inf 2x2", dmatrix![f64::NEG_INFINITY, 0.0; 0.0, 1.0]),
        ("All NaN 2x2", dmatrix![f64::NAN, f64::NAN; f64::NAN, f64::NAN]),
        ("All zero 2x2", dmatrix![0.0, 0.0; 0.0, 0.0]),
    ];

    for (name, m) in &matrices {
        let chol = m.clone().cholesky();
        println!("{}: cholesky = {:?}", name, chol.is_some());

        let svd = m.clone().svd(true, true);
        let pinv = svd.pseudo_inverse(1e-9);
        println!("{}: pseudo_inverse = {:?}", name, pinv.is_ok());

        if pinv.is_ok() {
            println!("  result diagonal: {} {}", pinv.as_ref().unwrap()[(0,0)],
                if pinv.as_ref().unwrap().nrows() > 1 { pinv.as_ref().unwrap()[(1,1)] } else { 0.0 });
        }
    }

    // Also try empty 0x0 matrix
    let empty = nalgebra::DMatrix::<f64>::zeros(0, 0);
    println!("\nEmpty 0x0: cholesky = {:?}", empty.clone().cholesky().is_some());
    let svd_e = empty.clone().svd(true, true);
    println!("Empty 0x0: pseudo_inverse = {:?}", svd_e.pseudo_inverse(1e-9).is_ok());
}
