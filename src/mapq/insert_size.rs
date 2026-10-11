//! Insert-size distribution (online learning).

use crate::types::InsertSizeDistribution;

const PRIOR_COUNT: f64 = 10.0;
const SD_CORRECTION: f64 = 0.9707;

impl Default for InsertSizeDistribution {
    fn default() -> Self {
        Self::new()
    }
}

impl InsertSizeDistribution {
    pub fn new() -> Self {
        Self::with_params(250.0, 80.0)
    }

    /// Start from a caller-supplied initial mean/SD (e.g. from CLI options);
    /// these are the priors that online learning then refines from observed
    /// proper pairs.
    pub fn with_params(mean: f64, sd: f64) -> Self {
        Self {
            mean,
            sd: sd / SD_CORRECTION,
            n_observed: 0,
            distance_sum: PRIOR_COUNT * mean,
            distance_square_sum: PRIOR_COUNT * (sd * sd + mean * mean),
        }
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
        // Stampy counts each mate as half an observation, with a ten-pair
        // prior. Both reports use the distribution at the start of the pair.
        if !(self.mean - 3.0 * self.sd < observed && observed < self.mean + 3.0 * self.sd) {
            return;
        }
        for _ in 0..2 {
            self.distance_sum += observed * 0.5;
            self.distance_square_sum += observed * observed * 0.5;
        }
        self.n_observed += 1;
        let count = PRIOR_COUNT + self.n_observed as f64;
        self.mean = self.distance_sum / count;
        self.sd = (self.distance_square_sum / count - self.mean * self.mean)
            .max(0.0)
            .sqrt()
            / SD_CORRECTION;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_stampy_insert_updates() {
        let mut model = InsertSizeDistribution::with_params(250.0, 60.0);
        for line in include_str!("../../tests/fixtures/stampy_insert.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let v: Vec<f64> = line.split('\t').map(|x| x.parse().unwrap()).collect();
            model.update(v[0]);
            assert!((model.mean - v[1]).abs() < 1e-10, "{line}");
            assert!((model.sd - v[2]).abs() < 1e-10, "{line}");
            assert_eq!(model.n_observed, v[3] as usize);
        }
    }
}
