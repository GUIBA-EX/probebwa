//! Candidate posterior, missing-locus statistics, entropy and BAQ evidence.
//! Alignment likelihoods come from the native-compatible coarse/sensitive
//! engines. Repeat-mask generation remains specific to this mapper.

use crate::types::*;

/// Phred-per-nat conversion constant, `0.1 * ln(0.1)`. Negative; dividing a
/// natural log of a probability (`<= 0`) by this yields a non-negative Phred
/// value.
const DB: f64 = -0.1 * std::f64::consts::LN_10; // = 0.1 * ln(0.1)
/// Floor under which a "probability of being wrong" is treated as exactly
/// this small, avoiding `ln(0)`.
fn posterior_floor() -> f64 {
    (99.0 * DB).exp()
}

/// Phred-space "probabilistic OR": approximates
/// `-10*log10(10^(-a/10) + 10^(-b/10))` via a coarse lookup table rather
/// than the exact (and slower) log-sum-exp. Used throughout to combine
/// independent phred-scale "this is wrong" probabilities.
pub(crate) fn logadd(a: f64, b: f64) -> f64 {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    if b >= a + 10.0 {
        a
    } else if b >= a + 4.0 {
        (a - 1.0).max(0.0)
    } else if b >= a + 2.0 {
        (a - 2.0).max(0.0)
    } else {
        (a - 3.0).max(0.0)
    }
}

/// The paired entropy combination must retain negative evidence: clamping
/// it here can turn a rejected random alignment into a mapped record.
pub(crate) fn logadd_unclamped(a: f64, b: f64) -> f64 {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    if b >= a + 10.0 { a } else if b >= a + 4.0 { a - 1.0 }
    else if b >= a + 2.0 { a - 2.0 } else { a - 3.0 }
}

/// Python's `int(0.5 + x)`: round-half-up via truncation-toward-zero after
/// adding 0.5 (not `floor`, which disagrees with `int()` for negative `x` —
/// and the `logentropy` terms below are frequently large and negative).
fn py_int_round(x: f64) -> f64 {
    (x + 0.5).trunc()
}

pub struct MapqCalculator {
    /// `--substitutionrate` expressed as a Phred value (port of
    /// `self.substitutionrate_phred`), used to inflate each base's
    /// effective error probability via `logadd`.
    substitution_rate_phred: f64,
    /// Total reference size, used by `logentropy` as the "random sequence"
    /// search space (port of `self.genome.getgenomesize()`); a bigger
    /// genome needs a more surprising alignment to be called non-random.
    genome_size: f64,
}

impl MapqCalculator {
    pub fn new(sub_rate: f64, _seed_k: usize, genome_size: usize) -> Self {
        let substitution_rate_phred = py_int_round(sub_rate.max(1e-30).ln() / DB);
        Self {
            substitution_rate_phred,
            genome_size: (genome_size.max(1)) as f64,
        }
    }

    pub fn adjusted_qualities(&self, qualities: &[u8]) -> Vec<u8> {
        qualities.iter().map(|&q| logadd(q as f64, self.substitution_rate_phred) as u8).collect()
    }

    // ------------------------------------------------------------------

    // ------------------------------------------------------------------
    // 2. Posterior over candidates (softmax) — unchanged in method, now
    //    exposed separately so the winning candidate's own CIGAR can be
    //    used for `logentropy` afterwards (pick the best candidate first,
    //    then compute entropy just for that one).
    // ------------------------------------------------------------------

    /// Softmax over candidate log-likelihoods. Returns `(best_idx, sum_exp)`
    /// where `sum_exp` is the *unnormalized* softmax denominator (the best
    /// candidate's own term is `exp(0) = 1`, so `1/sum_exp` is exactly its
    /// normalized posterior — this is the quantity `wrongmapping` below is
    /// built from).
    pub fn select_best(&self, log_lks: &[f64]) -> (usize, f64) {
        let max_lk = *log_lks.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let sum_exp: f64 = log_lks.iter().map(|lk| (lk - max_lk).exp()).sum();
        // Preserve the candidate order on ties; Iterator::max_by alone
        // returns the last maximum and would override the chosen priority.
        let best_idx = log_lks.iter().enumerate()
            .max_by(|(i, a), (j, b)| a.partial_cmp(b).unwrap().then_with(|| j.cmp(i)))
            .map(|(i, _)| i)
            .unwrap_or(0);
        (best_idx, sum_exp)
    }

    /// Per-candidate normalized posterior (softmax), rather than just the
    /// winner + denominator `select_best` returns. Used by
    /// `mapper::paired` to build the candidate shortlist: the locations
    /// that together constitute 99.9% of the single-read posterior mapping
    /// probability.
    pub fn posterior_fractions(&self, log_lks: &[f64]) -> Vec<f64> {
        if log_lks.is_empty() {
            return Vec::new();
        }
        let max_lk = log_lks.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exps: Vec<f64> = log_lks.iter().map(|lk| (lk - max_lk).exp()).collect();
        let sum: f64 = exps.iter().sum();
        exps.iter().map(|e| e / sum).collect()
    }

    /// `P(some other candidate is actually correct)`, i.e. the probability
    /// mass softmax puts everywhere except the best candidate. Exposed
    /// separately (as well as folded into `finalize_mapq`'s `wrongmapping`
    /// term) so `mapper::paired` can check it against the fast-path
    /// thresholds (both single reads mapping sufficiently uniquely: a
    /// posterior probability of less than 1% of having been mapped
    /// incorrectly due to near-repetitiveness) without having to run a full
    /// `finalize_mapq` (which also needs a concrete alignment/CIGAR this
    /// check happens before committing to).
    pub fn wrong_mapping_probability(&self, sum_exp: f64) -> f64 {
        let best_posterior = 1.0 / sum_exp.max(1.0);
        (1.0 - best_posterior).max(posterior_floor())
    }

    /// `missing_locus_phred`, exposed for the same reason as
    /// `wrong_mapping_probability` above (the fast path also requires the
    /// estimated probability of not having found the correct candidate to
    /// be sufficiently low for both mates).
    pub fn missing_locus_probability_estimate(&self, read: &Read, repeat_mask: &[f64]) -> f64 {
        10f64.powf(-self.missing_locus_phred(read, repeat_mask) / 10.0)
    }

    /// This is a rejection decision, separate from MAPQ: a repetitive but
    /// non-random alignment can legitimately be mapped with MAPQ zero.
    pub fn passes_random_match_test(&self, cigar: &[CigarOp], bases: &[u8], qualities: &[u8], reference: &[u8]) -> bool {
        self.log_entropy(cigar, bases, qualities, reference) > 3.0
    }

    pub(crate) fn mapping_evidence(&self, sum_exp: f64, read: &Read, repeat_mask: &[f64]) -> (f64, f64, f64) {
        let (missing, expected) = self.missing_locus_statistics(read, repeat_mask);
        if sum_exp == 0.0 { return (99.0, 0.0, expected); }
        (py_int_round(self.wrong_mapping_probability(sum_exp).ln() / DB), missing, expected)
    }

    pub(crate) fn pair_entropy(&self, cigar: &[CigarOp], bases: &[u8], qualities: &[u8], reference: &[u8], adjusted: bool) -> f64 {
        if adjusted { self.log_entropy(cigar, bases, qualities, reference) }
        else { f64::from(super::stampy::alignment_entropy(cigar, bases, qualities, reference,
            self.genome_size as usize, self.substitution_rate_phred.min(10.0) as u8)) }
    }

    /// Final MAPQ for the winning candidate, combining candidate-posterior
    /// confidence (`wrongmapping`, from `sum_exp`), missing-locus risk
    /// (`notfoundprob`), and random-sequence entropy (`entropy`, from the
    /// winning alignment's own CIGAR) via `logadd`:
    /// `posterior_phred = min(99, logadd(entropy, logadd(wrongmapping,
    /// notfoundprob)))`, with base-alignment confidence folded into wrongmapping.
    /// Returns `(mapq,
    /// posterior_probability)`, or `(0, 0.0)` if `entropy <= 0` (not a
    /// significantly better explanation than a random reference sequence of
    /// this size).
    #[allow(clippy::too_many_arguments)]
    pub fn finalize_mapq(
        &self,
        sum_exp: f64,
        read: &Read,
        winning_cigar: &[CigarOp],
        query_bases: &[u8],
        query_qual: &[u8],
        ref_window: &[u8],
        repeat_mask: &[f64],
        base_alignment_phred: f64,
    ) -> (u8, f64) {
        // --- wrongmapping: Phred(P(some other candidate is actually correct)) ---
        let p_wrong_posterior = self.wrong_mapping_probability(sum_exp);
        let wrongmapping = py_int_round(p_wrong_posterior.ln() / DB);

        // --- entropy: combinatorial surprise-vs-random-genome estimate ---
        let entropy = self.log_entropy(winning_cigar, query_bases, query_qual, ref_window) - 3.0;
        if entropy <= 0.0 {
            return (0, 0.0);
        }

        // --- notfoundprob: missing-locus probability (see
        // `missing_locus_phred`), converted to Phred ---
        let notfound_phred = self.missing_locus_phred(read, repeat_mask);

        let combined = logadd(logadd(base_alignment_phred, wrongmapping), notfound_phred);
        let posterior_phred = logadd(entropy, combined).clamp(0.0, 99.0);

        let mapq = posterior_phred.round() as u8;
        let posterior_prob = 10f64.powf(-posterior_phred / 10.0);
        (mapq, 1.0 - posterior_prob)
    }

    /// Stampy's random-match statistic, using its effective qualities,
    /// fixed alphabet parameter, indel boundaries and integer rounding.
    fn log_entropy(&self, cigar: &[CigarOp], query_bases: &[u8], query_qual: &[u8], ref_window: &[u8]) -> f64 {
        let qualities: Vec<u8> = query_qual.iter().map(|&q|
            logadd(f64::from(q), self.substitution_rate_phred) as u8).collect();
        f64::from(super::stampy::alignment_entropy(cigar, query_bases, &qualities, ref_window,
            self.genome_size as usize, self.substitution_rate_phred.min(10.0) as u8))
    }

    /// Stampy's length-dependent missing-candidate model. The current seed
    /// engine supplies repeat probabilities; convert these to its Phred mask
    /// convention. Exact seeding/mask parity remains a separate requirement.
    fn missing_locus_phred(&self, read: &Read, repeat_mask: &[f64]) -> f64 {
        self.missing_locus_statistics(read, repeat_mask).0
    }

    fn missing_locus_statistics(&self, read: &Read, repeat_mask: &[f64]) -> (f64, f64) {
        let qualities: Vec<u8> = read.qualities.scores.iter().map(|&q|
            logadd(f64::from(q), self.substitution_rate_phred) as u8).collect();
        let mask: Vec<u16> = (0..qualities.len()).map(|i| {
            let p = repeat_mask.get(i).copied().unwrap_or(0.0);
            if p <= 0.0 { 1000 } else { (-10.0 * p.min(1.0).log10()).round() as u16 }
        }).collect();
        let (missing, expected) = super::stampy::mutation_statistics(&qualities, &mask, -1);
        (f64::from(missing), f64::from(expected))
    }
}
