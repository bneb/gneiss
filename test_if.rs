fn main() {
    let f1 = 1575.42e6;
    let f2 = 1227.60e6;
    let gamma = (f1 * f1) / (f2 * f2);
    let alpha = gamma / (gamma - 1.0);
    let beta = -1.0 / (gamma - 1.0);
    println!("gamma: {}, alpha: {}, beta: {}", gamma, alpha, beta);
}
