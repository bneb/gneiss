use libm::pow;

/// Calculates the Fixed Failure-rate Ratio Test (FF-RT) threshold.
/// Based on Hou et al. (2016) "An Efficient Implementation of Fixed Failure-Rate Ratio Test for GNSS Ambiguity Resolution"
/// 
/// `n` is the number of ambiguities.
/// `pf` is the target failure rate (e.g., 0.01 or 0.001).
pub fn calculate_threshold(n: usize, pf: f64) -> f64 {
    let (p1, p2, _p3) = if pf <= 0.001 {
        (1.45, 0.25, 0.0)
    } else { // Fallback for e.g. 0.01
        (1.15, 0.25, 0.0)
    };

    let n_f = n as f64;
    p1 * pow(n_f, p2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ffrt_threshold_pf_001() {
        // Pf = 0.001 -> p1 = 1.45, p2 = 0.25, p3 = 0.0
        assert!((calculate_threshold(1, 0.001) - 1.45).abs() < 1e-3);
        assert!((calculate_threshold(5, 0.001) - 2.168).abs() < 1e-3);
        assert!((calculate_threshold(10, 0.001) - 2.579).abs() < 1e-3);
        assert!((calculate_threshold(20, 0.001) - 3.066).abs() < 1e-3);
    }

    #[test]
    fn test_ffrt_threshold_pf_01() {
        // Pf = 0.01 -> p1 = 1.15, p2 = 0.25, p3 = 0.0
        assert!((calculate_threshold(1, 0.01) - 1.15).abs() < 1e-3);
        assert!((calculate_threshold(5, 0.01) - 1.719).abs() < 1e-3);
        assert!((calculate_threshold(10, 0.01) - 2.045).abs() < 1e-3);
        assert!((calculate_threshold(20, 0.01) - 2.431).abs() < 1e-3);
    }
}
