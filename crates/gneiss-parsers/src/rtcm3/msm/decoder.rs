//! MSM message to EpochObs decoding.

use super::signals::{msm_signal_rinex_code, split_rinex_code, P2_10, P2_24, P2_29, P2_31, RANGE_MS};
use super::{MsmMessage, MsmType};
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

impl MsmMessage {
    pub fn into_epoch_obs(&self) -> EpochObs {
        let time = GpsTime::new(0, self.header.epoch_time as f64 / 1000.0);

        let constellation = match self.header.message_number / 10 {
            107 => Constellation::Gps,
            108 => Constellation::Glonass,
            109 => Constellation::Galileo,
            110 => Constellation::Sbas,
            111 => Constellation::Qzss,
            112 => Constellation::Beidou,
            _ => Constellation::Gps,
        };

        let active_sats: Vec<u8> = (0..64u8)
            .filter(|&i| (self.masks.satellite_mask & (1 << (63 - i))) != 0)
            .map(|i| i + 1)
            .collect();
        let active_sigs: Vec<u8> = (0..32u8)
            .filter(|&j| (self.masks.signal_mask & (1 << (31 - j))) != 0)
            .map(|j| j + 1)
            .collect();
        let n_sig = active_sigs.len();

        let (pr_scale, ph_scale) = match self.msm_type {
            MsmType::Msm4 | MsmType::Msm5 => (P2_24, P2_29),
            MsmType::Msm6 | MsmType::Msm7 => (P2_29, P2_31),
        };
        let pr_bits = match self.msm_type {
            MsmType::Msm4 | MsmType::Msm5 => 15,
            _ => 20,
        };
        let ph_bits = match self.msm_type {
            MsmType::Msm4 | MsmType::Msm5 => 22,
            _ => 24,
        };
        let pr_sentinel = -(1i64 << (pr_bits - 1));
        let ph_sentinel = -(1i64 << (ph_bits - 1));

        let mut satellites: Vec<SatObs> = Vec::with_capacity(active_sats.len());
        let mut cell_idx = 0usize;
        for (sat_pos, &prn) in active_sats.iter().enumerate() {
            let sat = SatelliteId {
                constellation,
                prn,
            };
            let mut observations = Vec::new();

            let rough_m = self
                .satellite_data
                .rough_range_int_ms
                .get(sat_pos)
                .copied()
                .and_then(|int_ms| {
                    if int_ms == 255 {
                        return None;
                    }
                    let modulo = *self.satellite_data.rough_ranges.get(sat_pos)?;
                    Some(int_ms as f64 * RANGE_MS + modulo as f64 * P2_10 * RANGE_MS)
                });

            for (sig_pos, &sig_id) in active_sigs.iter().enumerate() {
                let mask_pos = sat_pos * n_sig + sig_pos;
                if !self.masks.cell_mask.get(mask_pos).copied().unwrap_or(false) {
                    continue;
                }
                let k = cell_idx;
                cell_idx += 1;

                let Some(rough_m) = rough_m else { continue };
                let Some((band, attribute)) = msm_signal_rinex_code(constellation, sig_id)
                    .and_then(split_rinex_code)
                else {
                    continue;
                };

                push_pseudorange_obs(
                    &mut observations,
                    self.signal_data.fine_pseudoranges.get(k).copied(),
                    pr_sentinel,
                    rough_m,
                    pr_scale,
                    band,
                    attribute,
                );
                push_carrier_phase_obs(
                    &mut observations,
                    self.signal_data.fine_phase_ranges.get(k).copied(),
                    self.signal_data.lock_time_indicators.get(k).copied(),
                    ph_sentinel,
                    rough_m,
                    ph_scale,
                    band,
                    attribute,
                    constellation,
                );
            }

            satellites.push(SatObs { sat, observations });
        }

        EpochObs { time, satellites }
    }
}

#[allow(clippy::too_many_arguments)]
fn push_pseudorange_obs(
    observations: &mut Vec<Observation>,
    fine_pr: Option<i32>,
    pr_sentinel: i64,
    rough_m: f64,
    pr_scale: f64,
    band: u8,
    attribute: char,
) {
    let Some(fine_pr) = fine_pr else { return };
    if fine_pr as i64 == pr_sentinel {
        return;
    }
    observations.push(Observation {
        code: ObsCode {
            obs_type: ObsType::Pseudorange,
            signal: SignalCode {
                freq_band: band,
                attribute,
            },
        },
        value: rough_m + fine_pr as f64 * pr_scale * RANGE_MS,
        lock_time: None,
        lli: None,
    });
}

#[allow(clippy::too_many_arguments)]
fn push_carrier_phase_obs(
    observations: &mut Vec<Observation>,
    fine_ph: Option<i32>,
    lock_time: Option<u16>,
    ph_sentinel: i64,
    rough_m: f64,
    ph_scale: f64,
    band: u8,
    attribute: char,
    constellation: Constellation,
) {
    if constellation == Constellation::Glonass {
        return;
    }
    let Some(fine_ph) = fine_ph else { return };
    if fine_ph as i64 == ph_sentinel {
        return;
    }
    let range_m = rough_m + fine_ph as f64 * ph_scale * RANGE_MS;
    let freq_hz = gneiss_core::frequencies::track_c_frequency(constellation, band, 0);
    observations.push(Observation {
        code: ObsCode {
            obs_type: ObsType::CarrierPhase,
            signal: SignalCode {
                freq_band: band,
                attribute,
            },
        },
        value: range_m * freq_hz / SPEED_OF_LIGHT_M_S,
        lock_time,
        lli: None,
    });
}
