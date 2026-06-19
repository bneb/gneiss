#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LossFunction {
    Huber(f64),
    Cauchy(f64),
}

impl LossFunction {
    /// Computes the scaling weight for robust optimization given a squared Mahalanobis distance.
    ///
    /// For negative squared distances or NaNs (which are theoretically impossible but could
    /// arise from numerical instability), this function gracefully returns 1.0 (no downweighting).
    ///
    /// # Huber Loss Weight
    /// The Huber weight $w(x)$ for squared Mahalanobis distance $x$ and threshold $k^2$ is:
    /// $$
    /// w(x) = \begin{cases}
    /// 1.0 & \text{if } x \le k^2 \\
    /// \sqrt{\frac{k^2}{x}} & \text{if } x > k^2
    /// \end{cases}
    /// $$
    ///
    /// # Cauchy Loss Weight
    /// The Cauchy weight $w(x)$ for threshold $k^2$ is:
    /// $$
    /// w(x) = \frac{1}{1 + \frac{x}{k^2}}
    /// $$
    pub fn weight(&self, squared_mahalanobis: f64) -> f64 {
        if squared_mahalanobis.is_sign_negative() || squared_mahalanobis.is_nan() {
            return 1.0;
        }
        match self {
            LossFunction::Huber(k2) => {
                if squared_mahalanobis <= *k2 {
                    1.0
                } else {
                    (*k2 / squared_mahalanobis).sqrt()
                }
            }
            LossFunction::Cauchy(k2) => 1.0 / (1.0 + squared_mahalanobis / *k2),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_huber_weight_below_threshold() {
        let loss = LossFunction::Huber(4.0);
        assert!((loss.weight(2.0) - 1.0).abs() < 1e-9);
        assert!((loss.weight(4.0) - 1.0).abs() < 1e-9);
        assert!((loss.weight(0.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_huber_weight_above_threshold() {
        let loss = LossFunction::Huber(4.0);
        assert!((loss.weight(16.0) - 0.5).abs() < 1e-9);
        assert!((loss.weight(100.0) - 0.2).abs() < 1e-9);
        // Test strict boundary > 4.0
        assert!((loss.weight(4.000001) - 0.999999875).abs() < 1e-9);
    }

    #[test]
    fn test_cauchy_weight() {
        let loss = LossFunction::Cauchy(4.0);
        assert!((loss.weight(0.0) - 1.0).abs() < 1e-9);
        assert!((loss.weight(4.0) - 0.5).abs() < 1e-9);
        assert!((loss.weight(12.0) - 0.25).abs() < 1e-9);
    }

    #[test]
    fn test_edge_cases() {
        let huber = LossFunction::Huber(1.0);
        let cauchy = LossFunction::Cauchy(1.0);

        assert!(huber.weight(1e10) < 1e-4);
        assert!(cauchy.weight(1e10) < 1e-9);
    }

    #[test]
    fn test_negative_distance() {
        let huber = LossFunction::Huber(1.0);
        let cauchy = LossFunction::Cauchy(1.0);

        assert_eq!(huber.weight(-1.0), 1.0);
        assert_eq!(cauchy.weight(-5.0), 1.0);
        // Test strict boundary < 0.0
        assert_eq!(huber.weight(-0.000001), 1.0);
        assert_eq!(cauchy.weight(-0.000001), 1.0);
    }

    #[test]
    fn test_nan_distance() {
        let huber = LossFunction::Huber(1.0);
        let cauchy = LossFunction::Cauchy(1.0);

        assert_eq!(huber.weight(f64::NAN), 1.0);
        assert_eq!(cauchy.weight(f64::NAN), 1.0);
    }

    #[test]
    fn test_infinity_distance() {
        let huber = LossFunction::Huber(1.0);
        let cauchy = LossFunction::Cauchy(1.0);

        assert_eq!(huber.weight(f64::INFINITY), 0.0);
        assert_eq!(cauchy.weight(f64::INFINITY), 0.0);
    }
}
