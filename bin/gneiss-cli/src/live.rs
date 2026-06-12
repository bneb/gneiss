use tokio::io::AsyncReadExt;
use tokio_serial::SerialPortBuilderExt;
use tracing::{info, error};
use gneiss_rtk::engine::{ProcessingEngine, EngineConfig};
use gneiss_ntrip::client::{NtripClient, NtripConfig};
use gneiss_core::time::GpsTime;
use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode, ObsType, SignalCode};
use gneiss_core::sat::{SatelliteId, Constellation};
use gneiss_parsers::ubx::UbxRxmRawx;

pub struct LiveConfig {
    pub port: String,
    pub baud: u32,
    pub ntrip_url: Option<String>,
    pub ntrip_mount: Option<String>,
    pub ntrip_user: Option<String>,
    pub ntrip_pass: Option<String>,
    pub _output: Option<String>,
}

pub async fn run_live(
    live_cfg: LiveConfig,
    config: EngineConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    info!("Starting Live Real-Time GNSS Engine on port {} @ {}", live_cfg.port, live_cfg.baud);
    
    let (rtcm_tx, mut rtcm_rx) = tokio::sync::mpsc::channel::<Vec<u8>>(100);
    
    if let (Some(url), Some(mount)) = (live_cfg.ntrip_url, live_cfg.ntrip_mount) {
        let ntrip_config = NtripConfig {
            server_url: url,
            mountpoint: mount,
            username: live_cfg.ntrip_user,
            password: live_cfg.ntrip_pass,
        };
        let client = NtripClient::new(ntrip_config);
        info!("Connecting to NTRIP...");
        let mut ntrip_stream = client.connect().await?;
        
        tokio::spawn(async move {
            let mut buffer = Vec::new();
            while let Some(bytes) = ntrip_stream.recv().await {
                buffer.extend_from_slice(&bytes);
                loop {
                    match gneiss_parsers::rtcm3::parse_rtcm3_frame(&buffer) {
                        Ok((rem, frame)) => {
                            let _ = rtcm_tx.send(frame.payload.to_vec()).await;
                            buffer = rem.to_vec();
                        }
                        Err(gneiss_parsers::rtcm3::RtcmParseError::Incomplete) => break,
                        Err(_) => {
                            // Invalid data, drop a byte to resync
                            if !buffer.is_empty() {
                                buffer.remove(0);
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
        });
    }

    let mut serial_port = tokio_serial::new(live_cfg.port, live_cfg.baud).open_native_async()?;
    let mut buffer = Vec::new();
    let mut read_buf = [0u8; 4096];
    
    let mut engine = ProcessingEngine::new(config);
    let mut base_meas_cache: Option<EpochObs> = None;
    
    loop {
        tokio::select! {
            Some(rtcm_payload) = rtcm_rx.recv() => {
                if rtcm_payload.len() < 2 { continue; }
                let msg_num = u16::from_be_bytes([rtcm_payload[0], rtcm_payload[1]]) >> 4;
                match msg_num {
                    1077 | 1087 | 1097 | 1127 => {
                        if let Ok(msm) = gneiss_parsers::rtcm3::msm::parse_msm_message(&rtcm_payload) {
                            base_meas_cache = Some(msm.into_epoch_obs());
                        }
                    }
                    _ => {}
                }
            }
            res = serial_port.read(&mut read_buf) => {
                match res {
                    Ok(0) => break,
                    Ok(n) => {
                        buffer.extend_from_slice(&read_buf[..n]);
                        loop {
                            match gneiss_parsers::ubx::parse_ubx_frame(&buffer) {
                                Ok((rem, frame)) => {
                                    if frame.class == 0x02 && frame.id == 0x15 {
                                        if let Ok(rawx) = gneiss_parsers::ubx::parse_rxm_rawx(frame.payload) {
                                            let rinex = rawx_to_epoch(&rawx);
                                            let _ = engine.process_epoch(&rinex, base_meas_cache.as_ref());
                                            
                                            // Real-time output handling could go here
                                            if let Some(ref state) = engine.current_state {
                                                info!("Live Pos: {:.4}, {:.4}, {:.4} | Fix: {}", 
                                                      state.position.vector.x, state.position.vector.y, state.position.vector.z, state.is_fixed);
                                            }
                                        }
                                    }
                                    buffer = rem.to_vec();
                                }
                                Err(gneiss_parsers::ubx::UbxParseError::Incomplete) => break,
                                Err(_) => {
                                    if !buffer.is_empty() {
                                        buffer.remove(0);
                                    } else {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    Err(e) => {
                        error!("Serial read error: {}", e);
                        break;
                    }
                }
            }
        }
    }
    
    Ok(())
}

fn rawx_to_epoch(rawx: &UbxRxmRawx) -> EpochObs {
    let time = GpsTime::new(rawx.week.into(), rawx.rcv_tow);
    let mut sats = Vec::new();
    
    for meas in &rawx.measurements {
        let constel = match meas.gnss_id {
            0 => Constellation::Gps,
            2 => Constellation::Galileo,
            3 => Constellation::Beidou,
            6 => Constellation::Glonass,
            _ => continue, // Unknown/unsupported
        };
        let sat_id = SatelliteId { constellation: constel, prn: meas.sv_id };
        
        // Simple mapping, UBX sigId mapping is complex, assume L1/E1/B1 for now
        let sig = SignalCode { freq_band: 1, attribute: 'C' };
        
        let mut obs = Vec::new();
        if meas.pr_valid {
            obs.push(Observation { code: ObsCode { obs_type: ObsType::Pseudorange, signal: sig }, value: meas.pr_mes, lock_time: None });
        }
        if meas.cp_valid {
            obs.push(Observation { code: ObsCode { obs_type: ObsType::CarrierPhase, signal: sig }, value: meas.cp_mes, lock_time: Some(meas.locktime) });
        }
        obs.push(Observation { code: ObsCode { obs_type: ObsType::Doppler, signal: sig }, value: meas.do_mes as f64, lock_time: None });
        obs.push(Observation { code: ObsCode { obs_type: ObsType::Snr, signal: sig }, value: meas.cno as f64, lock_time: None });
        
        sats.push(SatObs { sat: sat_id, observations: obs });
    }
    
    EpochObs {
        time,
        satellites: sats,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_parsers::ubx::RxmRawxMeas;

    #[test]
    fn test_rawx_to_epoch() {
        let mut rawx = UbxRxmRawx {
            rcv_tow: 100000.5,
            week: 2100,
            leap_s: 18,
            num_meas: 2,
            rec_stat: 0x01,
            version: 0x01,
            measurements: vec![],
        };

        // Add GPS Sat
        rawx.measurements.push(RxmRawxMeas {
            pr_mes: 20000000.5,
            cp_mes: 1000000.5,
            do_mes: 123.4,
            gnss_id: 0,
            sv_id: 12,
            sig_id: 0,
            freq_id: 0,
            locktime: 5000,
            cno: 45,
            pr_stdev: 0.04,
            cp_stdev: 0.032,
            do_stdev: 0.004,
            pr_valid: true,
            cp_valid: true,
            half_cycle_valid: false,
            sub_half_cycle: false,
        });

        // Add Galileo Sat
        rawx.measurements.push(RxmRawxMeas {
            pr_mes: 25000000.5,
            cp_mes: 1200000.5,
            do_mes: -50.2,
            gnss_id: 2,
            sv_id: 5,
            sig_id: 0,
            freq_id: 0,
            locktime: 4000,
            cno: 42,
            pr_stdev: 0.04,
            cp_stdev: 0.032,
            do_stdev: 0.004,
            pr_valid: true,
            cp_valid: false, // Intentionally missing phase
            half_cycle_valid: false,
            sub_half_cycle: false,
        });

        let epoch = rawx_to_epoch(&rawx);

        assert_eq!(epoch.time.week, 2100);
        assert_eq!(epoch.time.tow, 100000.5);
        assert_eq!(epoch.satellites.len(), 2);

        // Check GPS 12
        let sat1 = &epoch.satellites[0];
        assert_eq!(sat1.sat.constellation, Constellation::Gps);
        assert_eq!(sat1.sat.prn, 12);
        assert_eq!(sat1.observations.len(), 4); // pr, cp, do, snr
        assert_eq!(sat1.observations[0].code.obs_type, ObsType::Pseudorange);
        assert_eq!(sat1.observations[0].value, 20000000.5);
        assert_eq!(sat1.observations[1].code.obs_type, ObsType::CarrierPhase);
        assert_eq!(sat1.observations[1].lock_time, Some(5000));

        // Check Galileo 5
        let sat2 = &epoch.satellites[1];
        assert_eq!(sat2.sat.constellation, Constellation::Galileo);
        assert_eq!(sat2.sat.prn, 5);
        assert_eq!(sat2.observations.len(), 3); // pr, do, snr (no cp)
    }
}
