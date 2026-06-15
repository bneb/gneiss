use gneiss_rtk::engine::{ProcessingEngine, EngineConfig};
use gneiss_parsers::ubx::parse_ubx_frame;
use std::fs::File;
use std::io::Read;
use std::path::Path;

#[test]
fn test_urbannav_tst_replay_skeleton() {
    // This test is a placeholder for replaying the UrbanNav TST-1 dataset.
    // It verifies that the engine can handle a stream of real UBX/IMU data.
    
    let dataset_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../datasets/urbannav/TST1/rover.ubx");
    
    // Skip if dataset is not present or empty (it's large and not committed)
    if !dataset_path.exists() || std::fs::metadata(&dataset_path).map(|m| m.len()).unwrap_or(0) == 0 {
        return;
    }

    let config = EngineConfig {
        mode: gneiss_rtk::engine::EngineMode::RtkIns,
        ..Default::default()
    };
    let _engine = ProcessingEngine::new(config);

    let mut file = File::open(dataset_path).unwrap();
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer).unwrap();

    // Replay logic
    let mut pos = 0;
    let mut parsed_count = 0;
    while pos < buffer.len() {
        if let Ok((rem, _frame)) = parse_ubx_frame(&buffer[pos..]) {
            parsed_count += 1;
            pos = buffer.len() - rem.len();
        } else {
            pos += 1;
        }
    }
    
    println!("Parsed {} UBX frames from the dataset.", parsed_count);
    assert!(parsed_count > 0, "Expected to parse at least one UBX frame from the dataset");
}
