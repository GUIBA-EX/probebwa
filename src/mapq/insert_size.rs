//! Insert-size distribution (online learning).

use crate::types::InsertSizeDistribution;

impl Default for InsertSizeDistribution {
    fn default() -> Self {
        Self::new()
    }
}

impl InsertSizeDistribution {
    pub fn new() -> Self {
        Self::with_params(250.0, 60.0)
    }

    /// Start from a caller-supplied initial mean/SD (e.g. from CLI options);
    /// these are the priors that online learning then refines from observed
    /// proper pairs.
    pub fn with_params(mean: f64, sd: f64) -> Self {
        Self { mean, sd, n_observed: 0 }
    }

    pub fn log_likelihood(&self, insert_size: f64) -> f64 {
        if self.sd <= 0.0 {
            return -1e300;
        }
        let diff = insert_size - self.mean;
        -0.5 * (diff * diff) / (self.sd * self.sd)
            - 0.5 * (2.0 * std::f64::consts::PI * self.sd * self.sd).ln()
    }

    pub fn update(&mut self, observed: f64) {
        self.n_observed += 1;
        let n = self.n_observed as f64;
        let delta = observed - self.mean;
        self.mean += delta / n;
        let delta2 = observed - self.mean;
        let variance = if n > 1.0 {
            ((self.sd * self.sd) * (n - 1.0) + delta * delta2) / n
        } else {
            self.sd * self.sd
        };
        self.sd = variance.sqrt().max(1.0);
    }
}
