//! MSM message to EpochObs decoding.

use super::signals::{msm_signal_rinex_code, split_rinex_code, P2_10, P2_24, P2_29, P2_31, RANGE_MS};
use super::{MsmMessage, MsmType};
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::{GpsTime, BDT_OFFSET_SECONDS};

/// Seconds in a day, used for the GLONASS day wrap.
const SECONDS_PER_DAY: f64 = 86_400.0;
/// Seconds in a week.
const SECONDS_PER_WEEK: f64 = 604_800.0;
/// GLONASS broadcast time counts UTC(SU) plus three hours. RTKLIB `adjday_glot`
/// (`rtcm3.c`) adds 10800 s when entering GLONASS time and subtracts it again
/// when converting the decoded epoch back to GPST.
const GLONASS_TIME_OFFSET_S: f64 = 10_800.0;
/// Unix timestamp of the GPS epoch, 1980-01-06 00:00:00 UTC.
const GPS_EPOCH_UNIX_S: f64 = 315_964_800.0;

/// Approximate current GPST, used to resolve an MSM epoch that carries only a
/// time of week. Within ~37 s (the accumulated leap seconds); irrelevant for
/// week resolution, which is why the current time is only an anchor.
fn current_gps_time() -> GpsTime {
    let unix = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs_f64())
        .unwrap_or(GPS_EPOCH_UNIX_S);
    let since = unix - GPS_EPOCH_UNIX_S;
    let week = (since / SECONDS_PER_WEEK).floor().max(0.0);
    GpsTime::new(week as u32, since - week * SECONDS_PER_WEEK)
}

/// Places a bare time of week on the week that puts it closest to `reference`.
fn week_from_tow(tow: f64, reference: GpsTime) -> GpsTime {
    let target = reference.week as f64 * SECONDS_PER_WEEK + reference.tow;
    let week = ((target - tow) / SECONDS_PER_WEEK).round().max(0.0);
    GpsTime::new(week as u32, tow)
}

/// Resolves a GLONASS MSM epoch.
///
/// DF393 is a 3-bit day-of-week and DF392 a 27-bit time of day in milliseconds;
/// the time of day is expressed in GLONASS time (UTC(SU) + 3 h). This mirrors
/// RTKLIB `adjday_glot`, which anchors on the UTC day containing the reference
/// epoch, wraps the message time into +/-12 h of that anchor, and removes the
/// three-hour offset to return to GPST.
fn glonass_epoch(epoch_field: u32, reference: GpsTime) -> GpsTime {
    let tod = (epoch_field & 0x07FF_FFFF) as f64 / 1000.0;
    let glonass_tow = reference.week as f64 * SECONDS_PER_WEEK + reference.tow + GLONASS_TIME_OFFSET_S;
    let day_start = (glonass_tow / SECONDS_PER_DAY).floor() * SECONDS_PER_DAY;
    let tod_anchor = glonass_tow - day_start;
    let tod = if tod < tod_anchor - SECONDS_PER_DAY / 2.0 {
        tod + SECONDS_PER_DAY
    } else if tod > tod_anchor + SECONDS_PER_DAY / 2.0 {
        tod - SECONDS_PER_DAY
    } else {
        tod
    };
    let gpst = day_start + tod - GLONASS_TIME_OFFSET_S;
    let week = (gpst / SECONDS_PER_WEEK).floor().max(0.0);
    GpsTime::new(week as u32, gpst - week * SECONDS_PER_WEEK)
}

impl MsmMessage {
    /// Decodes this message into an epoch of observations, anchoring the MSM
    /// time-of-week on the current date.
    pub fn into_epoch_obs(&self) -> EpochObs {
        self.into_epoch_obs_at(current_gps_time())
    }

    /// Decodes this message into an epoch of observations, anchoring the MSM
    /// time-of-week on `reference`.
    ///
    /// An MSM epoch carries a time of week but no week number, so it cannot be
    /// made absolute without an anchor; without one every observation would be
    /// stamped in GPS week 0 (1980-01-06).
    pub fn into_epoch_obs_at(&self, reference: GpsTime) -> EpochObs {
        let constellation = match self.header.message_number / 10 {
            107 => Constellation::Gps,
            108 => Constellation::Glonass,
            109 => Constellation::Galileo,
            110 => Constellation::Sbas,
            111 => Constellation::Qzss,
            112 => Constellation::Beidou,
            _ => Constellation::Gps,
        };

        let time = match constellation {
            Constellation::Glonass => glonass_epoch(self.header.epoch_time, reference),
            // BeiDou transmits the epoch in BDT; BDT leads GPST by 14 s.
            Constellation::Beidou => week_from_tow(
                self.header.epoch_time as f64 / 1000.0 + BDT_OFFSET_SECONDS,
                reference,
            ),
            _ => week_from_tow(self.header.epoch_time as f64 / 1000.0, reference),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rtcm3::msm::parse_msm_message;

    /// Packs (width, value) pairs big-endian into a byte buffer.
    fn pack(pairs: &[(usize, u64)]) -> Vec<u8> {
        let total: usize = pairs.iter().map(|(b, _)| *b).sum();
        let mut bytes = vec![0u8; total.div_ceil(8)];
        let mut pos = 0;
        for &(bits, val) in pairs {
            for i in 0..bits {
                if (val >> (bits - 1 - i)) & 1 != 0 {
                    bytes[pos / 8] |= 1 << (7 - (pos % 8));
                }
                pos += 1;
            }
        }
        bytes
    }

    /// MSM4 with one satellite and one signal; `epoch` is DF004 (ms of week)
    /// and `fine_pr` / `fine_ph` are DF399 / DF400 in raw units.
    fn msm4(epoch: u32, rough_ms: u8, rate: u16, fine_pr: i32, fine_ph: i32) -> Vec<u8> {
        let sign = |v: i32, bits: usize| (v as u32 as u64) & ((1u64 << bits) - 1);
        pack(&[
            (12, 1074), (12, 1), (30, epoch as u64), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0),
            (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, rough_ms as u64), (10, rate as u64),
            (15, sign(fine_pr, 15)), (22, sign(fine_ph, 22)), (4, 3), (1, 0), (6, 35),
        ])
    }

    /// Same as [`msm4`] but with an explicit message number (1074 GPS,
    /// 1084 GLONASS, 1124 BeiDou) and epoch field.
    fn msm_n(msg_num: u16, epoch: u32) -> Vec<u8> {
        let sign = |v: i32, bits: usize| (v as u32 as u64) & ((1u64 << bits) - 1);
        pack(&[
            (12, msg_num as u64), (12, 1), (30, epoch as u64), (1, 0), (3, 2), (7, 0), (2, 0),
            (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (10, 100),
            (15, sign(50, 15)), (22, sign(500, 22)), (4, 3), (1, 0), (6, 35),
        ])
    }

    /// An MSM epoch carries a time of week but no week number. Without an
    /// anchor it used to land in GPS week 0, i.e. 1980-01-06.
    #[test]
    fn epoch_week_is_anchored_on_the_reference() {
        let msg = parse_msm_message(&msm4(345_600_000, 2, 100, 50, 500)).unwrap();
        // Reference at the same time of week: the mapping is then unambiguous.
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
        assert_eq!(epoch.time.week, 2110);
        assert!((epoch.time.tow - 345_600.0).abs() < 1e-9);
    }

    /// The convenience entry point must never fall back to week 0.
    #[test]
    fn default_epoch_is_not_week_zero() {
        let msg = parse_msm_message(&msm4(345_600_000, 2, 100, 50, 500)).unwrap();
        let epoch = msg.into_epoch_obs();
        assert!(
            epoch.time.week > 1000,
            "epoch week {} is implausible for any live GNSS stream",
            epoch.time.week
        );
    }

    /// BeiDou MSM epochs are in BDT, which leads GPST by exactly 14 s
    /// (`gneiss_core::time::BDT_OFFSET_SECONDS`, asserted to be 14.0 in core).
    #[test]
    fn beidou_epoch_carries_the_bdt_gpst_offset() {
        let payload = msm_n(1124, 345_600_000);
        let msg = parse_msm_message(&payload).unwrap();
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
        assert_eq!(epoch.time.week, 2110);
        assert!(
            (epoch.time.tow - (345_600.0 + 14.0)).abs() < 1e-9,
            "BeiDou epoch must be shifted by the 14 s BDT-GPST offset"
        );
        assert_eq!(epoch.satellites.len(), 1);
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Beidou);
    }

    /// GLONASS MSM epochs are DF393 (3 bits day-of-week) + DF392 (27 bits time
    /// of day in ms), expressed in GLONASS time = UTC(SU) + 3 h.
    ///
    /// With a reference at GPST week 2110, tow 345600 (Wednesday 00:00) the
    /// GLONASS day anchor is 3 h = 10800 s into that day, so a message time of
    /// day of 10800 s must convert back to GPST 345600, and 49200 s
    /// (13 h 40 min GLONASS) must convert back to GPST 38400 s into the same
    /// day, i.e. tow 345600 + 38400 = 384000.
    #[test]
    fn glonass_epoch_is_glonass_time_minus_three_hours() {
        for (tod_ms, expected_tow) in [(10_800_000u32, 345_600.0), (49_200_000, 384_000.0)] {
            let payload = msm_n(1084, (5u32 << 27) | tod_ms);
            let msg = parse_msm_message(&payload).unwrap();
            let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
            assert_eq!(epoch.time.week, 2110);
            assert!(
                (epoch.time.tow - expected_tow).abs() < 1e-6,
                "GLONASS TOD {tod_ms} ms gave tow {} instead of {expected_tow}",
                epoch.time.tow
            );
            assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Glonass);
        }
    }

    /// Rough range = DF010 (integer ms) + DF398 * 2^-10 ms, in metres.
    /// DF010 = 200, DF398 = 512 -> 200 + 0.5 = 200.5 ms.
    #[test]
    fn rough_range_assembles_integer_ms_plus_fraction() {
        let msg = parse_msm_message(&msm4(345_600_000, 200, 512, 0, 0)).unwrap();
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
        let pr = epoch.satellites[0]
            .observations
            .iter()
            .find(|o| o.code.obs_type == ObsType::Pseudorange)
            .expect("pseudorange present");
        let expected = (200.5 / 1000.0) * SPEED_OF_LIGHT_M_S;
        assert!((pr.value - expected).abs() < 1e-6, "{} vs {expected}", pr.value);
    }

    /// DF010 = 255 marks the rough range as invalid; the satellite must produce
    /// no observations rather than a pseudorange near zero.
    #[test]
    fn invalid_rough_range_yields_no_observations() {
        let msg = parse_msm_message(&msm4(345_600_000, 255, 100, 50, 500)).unwrap();
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 0.0));
        assert_eq!(epoch.satellites.len(), 1);
        assert!(
            epoch.satellites[0].observations.is_empty(),
            "an invalid rough range must not yield observations"
        );
    }

    /// MSM4 fine-field sentinels are -(2^(bits-1)): DF399 = -16384 and
    /// DF400 = -2097152. They mean "no data" and must be dropped.
    #[test]
    fn fine_field_sentinels_are_dropped() {
        let msg = parse_msm_message(&msm4(345_600_000, 2, 100, -16384, -2097152)).unwrap();
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
        assert!(epoch.satellites[0].observations.is_empty());

        // The largest non-sentinel values must survive.
        let msg = parse_msm_message(&msm4(345_600_000, 2, 100, 16383, 2097151)).unwrap();
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
        assert_eq!(epoch.satellites[0].observations.len(), 2);
    }

    /// The cell mask is satellite-major (DF400: for each satellite, for each
    /// signal). With 2 satellites x 2 signals the fine pseudoranges must land
    /// on cells 0,1 for satellite 1 and cells 2,3 for satellite 2.
    #[test]
    fn cell_mask_is_satellite_major() {
        let payload = pack(&[
            (12, 1074), (12, 9), (30, 345_600_000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0),
            (1, 0), (3, 0),
            (64, (1u64 << 63) | (1u64 << 62)),
            (32, (1u64 << 31) | (1u64 << 30)),
            (1, 1), (1, 1), (1, 1), (1, 1),
            (8, 2), (8, 3), (10, 0), (10, 0),
            (15, 11), (15, 22), (15, 33), (15, 44),
            (22, 0), (22, 0), (22, 0), (22, 0),
            (4, 1), (4, 1), (4, 1), (4, 1),
            (1, 0), (1, 0), (1, 0), (1, 0),
            (6, 30), (6, 31), (6, 32), (6, 33),
        ]);
        let msg = parse_msm_message(&payload).unwrap();
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
        assert_eq!(epoch.satellites.len(), 2);
        let codes: Vec<Vec<f64>> = epoch
            .satellites
            .iter()
            .map(|s| {
                s.observations
                    .iter()
                    .filter(|o| o.code.obs_type == ObsType::Pseudorange)
                    .map(|o| o.value)
                    .collect()
            })
            .collect();
        assert_eq!(codes.len(), 2);
        assert!(codes[0].windows(2).all(|w| w[1] > w[0]), "sat 1 cells 11 then 22");
        assert!(codes[1].windows(2).all(|w| w[1] > w[0]), "sat 2 cells 33 then 44");
        assert!(codes[0][0] < codes[1][0], "satellite 1 fine values precede satellite 2");
    }

    /// GLONASS carrier phase is currently dropped: the wavelength depends on
    /// the FDMA channel in DF419 (extended satellite info), which the decoder
    /// parses but does not apply. Silently emitting a wrong-channel phase would
    /// be worse than emitting none, so the drop is asserted here as the current
    /// contract rather than left implicit.
    #[test]
    fn glonass_carrier_phase_is_not_emitted_without_the_fdma_channel() {
        let payload = pack(&[
            (12, 1084), (12, 9), (30, 49_200_000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0),
            (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (10, 0),
            (15, 100), (22, 500), (4, 3), (1, 0), (6, 35),
        ]);
        let msg = parse_msm_message(&payload).unwrap();
        let epoch = msg.into_epoch_obs_at(GpsTime::new(2110, 345_600.0));
        let has_phase = epoch.satellites[0]
            .observations
            .iter()
            .any(|o| o.code.obs_type == ObsType::CarrierPhase);
        assert!(!has_phase, "GLONASS phase must not be emitted without DF419");
    }
}
