use std::fs;

fn main() {
    let data = fs::read("datasets/rtkexplorer/sample_1/base.rtcm3").unwrap();
    let mut i = 0;
    while i < data.len() {
        if data[i] != 0xD3 {
            i += 1;
            continue;
        }
        if i + 3 > data.len() { break; }
        let len = ((data[i+1] as usize & 0x03) << 8) | (data[i+2] as usize);
        if i + 3 + len + 3 > data.len() { break; }
        let payload = &data[i+3..i+3+len];
        let msg_num = (payload[0] as u16) << 4 | (payload[1] as u16) >> 4;
        println!("Msg type: {}", msg_num);
        i += 3 + len + 3;
    }
}
