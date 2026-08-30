//! Per-satellite double-differenced pseudorange accumulator (O(1) ring buffer)
//! for fast code averaging across epochs in RTK mode.

use std::collections::{HashMap, VecDeque};

/// Key identifying a double-differenced satellite pair and frequency band.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DdSatKey {
    pub sat: u16,
    pub ref_sat: u16,
    pub frequency: u8,
}

/// Accumulator entry storing recent DD pseudorange residuals in an O(1) ring buffer.
#[derive(Debug, Clone)]
pub struct DdAccumulatorEntry {
    pub history: VecDeque<f64>,
    pub max_samples: usize,
    running_sum: f64,
}

impl DdAccumulatorEntry {
    pub fn new(max_samples: usize) -> Self {
        Self {
            history: VecDeque::with_capacity(max_samples),
            max_samples,
            running_sum: 0.0,
        }
    }

    #[inline]
    pub fn add(&mut self, sample: f64) {
        if self.history.len() >= self.max_samples {
            if let Some(old) = self.history.pop_front() {
                self.running_sum -= old;
            }
        }
        self.running_sum += sample;
        self.history.push_back(sample);
    }

    #[inline]
    pub fn mean(&self) -> Option<f64> {
        if self.history.is_empty() {
            None
        } else {
            Some(self.running_sum / self.history.len() as f64)
        }
    }

    pub fn std_dev(&self) -> Option<f64> {
        if self.history.len() < 2 {
            return None;
        }
        let mean = self.mean()?;
        let variance = self.history.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (self.history.len() - 1) as f64;
        Some(variance.sqrt())
    }

    #[inline]
    pub fn count(&self) -> usize {
        self.history.len()
    }
}

/// Multi-epoch DD pseudorange accumulator.
#[derive(Debug, Clone)]
pub struct DdPseudorangeAccumulator {
    entries: HashMap<DdSatKey, DdAccumulatorEntry>,
    max_samples: usize,
}

impl DdPseudorangeAccumulator {
    pub fn new(max_samples: usize) -> Self {
        Self {
            entries: HashMap::new(),
            max_samples,
        }
    }

    pub fn add_observation(&mut self, key: DdSatKey, value: f64) {
        self.entries
            .entry(key)
            .or_insert_with(|| DdAccumulatorEntry::new(self.max_samples))
            .add(value);
    }

    pub fn get_averaged(&self, key: &DdSatKey) -> Option<(f64, usize)> {
        let entry = self.entries.get(key)?;
        let mean = entry.mean()?;
        Some((mean, entry.count()))
    }

    pub fn get_stats(&self, key: &DdSatKey) -> Option<(f64, f64, usize)> {
        let entry = self.entries.get(key)?;
        let mean = entry.mean()?;
        let std_dev = entry.std_dev()?;
        Some((mean, std_dev, entry.count()))
    }

    pub fn clear(&mut self) {
        self.entries.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_accumulator_ring_buffer() {
        let mut acc = DdPseudorangeAccumulator::new(5);
        let key = DdSatKey { sat: 1, ref_sat: 2, frequency: 1 };

        for i in 1..=5 {
            acc.add_observation(key, i as f64);
        }
        let (mean, count) = acc.get_averaged(&key).unwrap();
        assert_eq!(count, 5);
        assert!((mean - 3.0).abs() < 1e-12); // (1+2+3+4+5)/5 = 3.0

        // Push 6th item — 1 drops out
        acc.add_observation(key, 6.0);
        let (mean2, count2) = acc.get_averaged(&key).unwrap();
        assert_eq!(count2, 5);
        assert!((mean2 - 4.0).abs() < 1e-12); // (2+3+4+5+6)/5 = 4.0
    }
}
