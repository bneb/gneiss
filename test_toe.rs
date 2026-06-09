use std::io::BufReader;

fn main() {
    let file_contents = "     3.03           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
R 6 2020 12 24 21 15  0  .189751386642E-03  .000000000000E+00  .422910000000E+06
     -.740158740234E+04 -.212037086487E+00  .000000000000E+00  .000000000000E+00
     -.206682856445E+05 -.176755714417E+01  .931322574615E-09 -.400000000000E+01
      .129489067383E+05 -.294115734100E+01 -.186264514923E-08  .000000000000E+00
";
    let mut reader = BufReader::new(file_contents.as_bytes());
    let eph = gneiss_parsers::rinex::parse_rinex_nav(&mut reader).unwrap();
    println!("toe.tow: {}", eph[0].toe().tow);
}
