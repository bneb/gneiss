use gneiss_parsers::rtcm3::{parse_rtcm3_frame, RtcmParseError, parse_msm_message, parse_1019};
use std::fs;
use std::path::PathBuf;

#[test]
fn test_real_world_rtcm3_parsing() {
    // Read the actual 30-second live capture from the Centipede RTK network
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("../../base_35JF.rtcm3");
    
    if !path.exists() {
        println!("Test skipped: Real world RTCM3 dataset not found.");
        return;
    }

    let mut buffer = fs::read(path).expect("Failed to read base_35JF.rtcm3");
    
    let mut total_frames = 0;
    let mut msm_frames = 0;
    let mut eph_frames = 0;
    
    loop {
        if buffer.is_empty() {
            break;
        }
        
        match parse_rtcm3_frame(&buffer) {
            Ok((remaining, frame)) => {
                total_frames += 1;
                
                // Attempt to parse the inner payload
                // The first 2 bytes of the payload contain the 12-bit message number
                if frame.payload.len() >= 2 {
                    let msg_num = (u16::from_be_bytes([frame.payload[0], frame.payload[1]]) >> 4) & 0x0FFF;
                    
                    if msg_num == 1019 {
                        if let Ok(_eph) = parse_1019(frame.payload) {
                            eph_frames += 1;
                        }
                    } else if let Ok(_msm) = parse_msm_message(frame.payload) {
                        msm_frames += 1;
                    }
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
            Err(RtcmParseError::CrcMismatch) | Err(RtcmParseError::UnsupportedMsmType) => {
                buffer.remove(0);
            }
        }
    }
    
    println!("Total Frames Parsed: {}", total_frames);
    println!("MSM Frames Decoded: {}", msm_frames);
    println!("GPS Ephemeris Frames Decoded: {}", eph_frames);
    
    assert!(total_frames > 10, "Should have parsed multiple frames from a 30s capture");
    // Some frames might be other types like 1005 (Station ARP), so msm_frames might be slightly less than total
    assert!(msm_frames > 0, "Should have parsed at least some MSM frames");
}
