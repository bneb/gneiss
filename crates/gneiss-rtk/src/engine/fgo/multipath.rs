use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct MultipathConfig {
    pub threshold: f64,
    pub max_history: usize,
    pub eval_window: usize,
}

pub const DEFAULT_MULTIPATH_THRESHOLD: f64 = 5.0;
pub const DEFAULT_MAX_HISTORY: usize = 10;
pub const DEFAULT_EVAL_WINDOW: usize = 5;

impl Default for MultipathConfig {
    fn default() -> Self {
        Self {
            threshold: DEFAULT_MULTIPATH_THRESHOLD,
            max_history: DEFAULT_MAX_HISTORY,
            eval_window: DEFAULT_EVAL_WINDOW,
        }
    }
}

pub struct MultipathEvaluator {
    pub history: HashMap<u32, Vec<f64>>,
    pub config: MultipathConfig,
}

impl Default for MultipathEvaluator {
    fn default() -> Self {
        Self::new(MultipathConfig::default())
    }
}

impl MultipathEvaluator {
    pub fn new(config: MultipathConfig) -> Self {
        Self {
            history: HashMap::new(),
            config,
        }
    }

    /// Evaluates if a satellite should be rejected due to multipath based on recent residuals.
    ///
    /// This uses a median-based filtering technique over a sliding window. It takes the
    /// `eval_window` most recent residuals and calculates their median. If the median
    /// exceeds the configured `threshold`, the satellite is considered to be experiencing
    /// sustained multipath and is rejected. This prevents isolated outliers from
    /// incorrectly rejecting the satellite.
    pub fn evaluate_and_reject(&mut self, sat_id: u32, new_residual: f64) -> bool {
        let hist = self.history.entry(sat_id).or_default();
        hist.push(new_residual);
        if hist.len() > self.config.max_history {
            hist.remove(0);
        }

        if hist.len() < self.config.eval_window {
            return false;
        }

        let start_idx = hist.len() - self.config.eval_window;
        let mut recent_window = hist[start_idx..].to_vec();
        recent_window.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        recent_window[self.config.eval_window / 2] > self.config.threshold
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_config() -> MultipathConfig {
        MultipathConfig {
            threshold: 3.0,
            max_history: 10,
            eval_window: 5,
        }
    }

    #[test]
    fn test_single_outlier_not_rejected() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 1;

        assert!(!evaluator.evaluate_and_reject(sat_id, 1.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 1.5));
        assert!(!evaluator.evaluate_and_reject(sat_id, 1.2));
        assert!(!evaluator.evaluate_and_reject(sat_id, 1.1));

        // Outlier
        assert!(!evaluator.evaluate_and_reject(sat_id, 10.0));
    }

    #[test]
    fn test_sustained_multipath_rejected() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 2;

        assert!(!evaluator.evaluate_and_reject(sat_id, 4.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 4.5));
        assert!(!evaluator.evaluate_and_reject(sat_id, 3.8));
        assert!(!evaluator.evaluate_and_reject(sat_id, 4.1));

        // 5th consecutive high residual -> rejection
        assert!(evaluator.evaluate_and_reject(sat_id, 5.0));
    }

    #[test]
    fn test_history_limit() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 3;

        for _ in 0..15 {
            let _ = evaluator.evaluate_and_reject(sat_id, 1.0);
        }

        if let Some(hist) = evaluator.history.get(&sat_id) {
            // Must be exactly max_history, not max_history + 1
            assert_eq!(hist.len(), evaluator.config.max_history);
        }
    }

    #[test]
    fn test_history_less_than_eval_window() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 4;

        // Even if residuals are very high, not rejected if we don't have enough history
        assert!(!evaluator.evaluate_and_reject(sat_id, 10.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 10.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 10.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 10.0));
    }

    #[test]
    fn test_nan_values_handled() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 5;

        assert!(!evaluator.evaluate_and_reject(sat_id, f64::NAN));
        assert!(!evaluator.evaluate_and_reject(sat_id, f64::NAN));
        assert!(!evaluator.evaluate_and_reject(sat_id, f64::NAN));
        assert!(!evaluator.evaluate_and_reject(sat_id, f64::NAN));
        // Median might be NaN, comparison with threshold > 3.0 returns false
        assert!(!evaluator.evaluate_and_reject(sat_id, f64::NAN));
    }

    #[test]
    fn test_boundary_conditions() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 6;

        // Exact threshold should NOT trigger rejection (it is strictly > threshold)
        assert!(!evaluator.evaluate_and_reject(sat_id, 3.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 3.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 3.0));
        assert!(!evaluator.evaluate_and_reject(sat_id, 3.0));

        // Slightly above threshold should trigger it, but we need at least 3 values
        // to change the median of a window of 5!
        assert!(!evaluator.evaluate_and_reject(sat_id, 3.000001));
        assert!(!evaluator.evaluate_and_reject(sat_id, 3.000001));
        assert!(evaluator.evaluate_and_reject(sat_id, 3.000001));
    }

    #[test]
    fn test_max_history_boundary() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 7;

        // Push max_history elements
        for i in 0..10 {
            evaluator.evaluate_and_reject(sat_id, i as f64);
            assert_eq!(evaluator.history.get(&sat_id).unwrap().len(), i + 1);
        }

        // Push 11th element, length should remain 10
        evaluator.evaluate_and_reject(sat_id, 10.0);
        assert_eq!(evaluator.history.get(&sat_id).unwrap().len(), 10);
        // The first element (0.0) should be removed
        assert_eq!(evaluator.history.get(&sat_id).unwrap().first(), Some(&1.0));
    }

    #[test]
    fn test_catch_minus_to_div_mutant() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 8;

        // Push 7 elements to make hist.len() = 7.
        evaluator.evaluate_and_reject(sat_id, 1.0); // hist[0]
        evaluator.evaluate_and_reject(sat_id, 1.0); // hist[1]

        evaluator.evaluate_and_reject(sat_id, 10.0); // hist[2]
        evaluator.evaluate_and_reject(sat_id, 10.0); // hist[3]
        evaluator.evaluate_and_reject(sat_id, 10.0); // hist[4]
        evaluator.evaluate_and_reject(sat_id, 1.0); // hist[5]
        let rejected = evaluator.evaluate_and_reject(sat_id, 1.0); // hist[6]

        assert!(rejected);
    }

    #[test]
    fn test_catch_div_to_mod_mutant() {
        let mut evaluator = MultipathEvaluator::new(test_config());
        let sat_id = 9;

        evaluator.evaluate_and_reject(sat_id, 1.0);
        evaluator.evaluate_and_reject(sat_id, 2.0);
        evaluator.evaluate_and_reject(sat_id, 4.0);
        evaluator.evaluate_and_reject(sat_id, 5.0);
        let rejected = evaluator.evaluate_and_reject(sat_id, 6.0);

        assert!(rejected);
    }
}
