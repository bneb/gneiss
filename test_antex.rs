use std::fs::File;
use std::io::{BufRead, BufReader};
fn main() {
    let line = "BLOCK IIA           G01                 G032      1992-079A TYPE / SERIAL NO    ";
    let type_str = line[0..20].trim();
    let serial = line[20..40].trim();
    println!("type: '{}', serial: '{}'", type_str, serial);
}
