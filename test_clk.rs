fn main() {
    let line = "AS G01  2020 12 24  0  0  0.000000  1    0.791467054363E-03";
    let bias_str = line[40..59].replace("D", "e");
    println!("bias_str: '{}'", bias_str);
    let bias = bias_str.trim().parse::<f64>().unwrap();
    println!("bias: {}", bias);
}
