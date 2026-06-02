use std::io::{self, Read};
use gneiss_parsers::rtcm3::{parse_rtcm3_frame, RtcmParseError};

fn main() -> io::Result<()> {
    let mut stdin = io::stdin().lock();
    let mut buffer = Vec::new();
    let mut chunk = [0u8; 1024];

    // Simple unix filter loop
    loop {
        let n = stdin.read(&mut chunk)?;
        if n == 0 {
            break;
        }
        buffer.extend_from_slice(&chunk[..n]);

        loop {
            if buffer.is_empty() {
                break;
            }

            match parse_rtcm3_frame(&buffer) {
                Ok((remaining, frame)) => {
                    // In a full implementation, we would extract the message type from the payload.
                    // For RTCM3, the first 12 bits of the payload contain the message number.
                    if frame.payload.len() >= 2 {
                        let msg_type = (u16::from_be_bytes([frame.payload[0], frame.payload[1]]) >> 4) & 0x0FFF;
                        println!("{{\"rtcm_frame\": {{\"msg_type\": {}, \"payload_bytes\": {}}}}}", msg_type, frame.payload.len());
                    } else {
                        println!("{{\"rtcm_frame\": {{\"payload_bytes\": {}}}}}", frame.payload.len());
                    }
                    
                    let rem_len = remaining.len();
                    let buf_len = buffer.len();
                    buffer.drain(..(buf_len - rem_len));
                }
                Err(RtcmParseError::Incomplete) => {
                    break;
                }
                Err(RtcmParseError::InvalidPreamble) => {
                    if let Some(pos) = buffer.iter().position(|&b| b == 0xD3) {
                        buffer.drain(..pos);
                    } else {
                        buffer.clear();
                        break;
                    }
                }
                Err(RtcmParseError::CrcMismatch) => {
                    eprintln!("{{\"error\": \"crc_mismatch\"}}");
                    buffer.remove(0);
                }
                Err(RtcmParseError::UnsupportedMsmType) => {
                    eprintln!("{{\"error\": \"unsupported_msm_type\"}}");
                    buffer.remove(0);
                }
            }
        }
    }

    Ok(())
}
