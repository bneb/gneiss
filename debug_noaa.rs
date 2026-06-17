use std::env;
#[tokio::main]
async fn main() {
    let client = reqwest::Client::new();
    let text = client.get("https://geodesy.noaa.gov/corsdata/coord/coord_14/p222_14.coord.txt").send().await.unwrap().text().await.unwrap();
    let mut in_itrf_block = false;
    for line in text.lines() {
        if line.contains("ITRF2014 POSITION") {
            in_itrf_block = true;
        } else if line.contains("NAD_83") || line.contains("VELOCITY") || line.contains("L1 Phase Center") {
            in_itrf_block = false;
        }
        if in_itrf_block {
            if line.contains("X =") {
                println!("Found X line: {}", line);
                if let Some(val_str) = line.split("X =").nth(1).and_then(|s| s.split('m').next()) {
                    println!("Parsed X string: '{}'", val_str);
                    println!("Parsed X value: {:?}", val_str.trim().parse::<f64>());
                }
            }
        }
    }
}
