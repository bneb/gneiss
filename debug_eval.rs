use std::fs;
fn main() {
    let sol_data = fs::read_to_string("benchmarks/rtklib_comparison/gneiss_Shinjuku_u-blox_spp.pos").unwrap();
    let mut sol_is_llh = sol_data.lines().any(|l| l.contains("latitude(deg) longitude(deg)"));
    println!("sol_is_llh: {}", sol_is_llh);
    let line = sol_data.lines().nth(1).unwrap();
    let parts: Vec<&str> = line.split_whitespace().collect();
    let x = parts[2].parse::<f64>().unwrap();
    let y = parts[3].parse::<f64>().unwrap();
    let z = parts[4].parse::<f64>().unwrap();
    println!("sol_ecef: {} {} {}", x, y, z);
    
    let truth_data = fs::read_to_string("datasets/urbannav/tokyo/Tokyo_Data/Shinjuku/reference.csv").unwrap();
    let line2 = truth_data.lines().nth(1).unwrap();
    let parts: Vec<&str> = line2.split(',').collect();
    let true_x = parts[5].trim().parse::<f64>().unwrap();
    let true_y = parts[6].trim().parse::<f64>().unwrap();
    let true_z = parts[7].trim().parse::<f64>().unwrap();
    println!("true_ecef: {} {} {}", true_x, true_y, true_z);
    
    println!("diff y: {}", y - true_y);
}
