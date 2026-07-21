//! Bayesian MAPQ and probabilistic scoring.
//!
//! This module combines several independent sources of evidence about a
//! candidate alignment's quality:
//!
//!   1. Per-base likelihood from quality values, with sequencing error and
//!      biological substitution rate combined via `logadd` (a Phred-space
//!      probabilistic-OR approximation) rather than an ad hoc multiplicative
//!      model.
//!   2. Posterior over candidate loci via softmax.
//!   3. `wrongmapping`: Phred of `1 - best_posterior`, the probability that
//!      some other candidate is actually the correct one.
//!   4. `entropy`: `logentropy`, a combinatorial (not just `-10log10(p)`)
//!      estimate of how surprising this good an alignment would be against
//!      a random sequence of the actual reference size — this doubles as
//!      the rejection gate (entropy <= 0 means "not better than chance")
//!      and as a MAPQ component.
//!   5. `notfoundprob`: this module's own simplified missing-locus model
//!      (see `missing_locus_probability`) — an approximation, not a full
//!      symmetric-polynomial treatment over per-base error probabilities.
//!   6. Final MAPQ = `min(99, logadd(entropy, logadd(wrongmapping,
//!      notfoundprob)))`, capped at 99 (not 255) — chosen after finding a
//!      real weakness by comparing against another tool on the same input:
//!      an earlier, simpler formula here saturated to a 255-equivalent cap
//!      on almost every accepted alignment, including a clearly spurious
//!      repeat-region hit, giving it the same "maximum confidence" score as
//!      a clean unique match. A base-alignment-quality (BAQ) term is
//!      omitted: it needs per-base alignment-quality output this module's
//!      aligner doesn't produce, so it isn't in the chain.

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
fn logadd(a: f64, b: f64) -> f64 {
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

/// Python's `int(0.5 + x)`: round-half-up via truncation-toward-zero after
/// adding 0.5 (not `floor`, which disagrees with `int()` for negative `x` —
/// and the `logentropy` terms below are frequently large and negative).
fn py_int_round(x: f64) -> f64 {
    (x + 0.5).trunc()
}

/// Phred approximation of `log(C(n, m))` (log of the binomial coefficient),
/// used by `logentropy` to penalize the number of ways `m` variant
/// positions could have landed among `n` candidate positions — i.e. the
/// "this could just as easily have been a different subset" correction.
/// Port of `phrednchoosem` (`pyx/*/maptools.pyx`).
fn phred_n_choose_m(n: i64, m: i64) -> f64 {
    let mut m = if m > n / 2 { n - m } else { m };
    if m <= 0 {
        return 0.0;
    }
    let mut n = n;
    let mut phred = 0.0f64;
    while m > 1 {
        phred += (n as f64).ln() - (m as f64).ln();
        n -= 1;
        m -= 1;
    }
    phred += (n as f64).ln();
    py_int_round(10.0 * phred / std::f64::consts::LN_10)
}

pub struct MapqCalculator {
    /// `--substitutionrate` expressed as a Phred value (port of
    /// `self.substitutionrate_phred`), used to inflate each base's
    /// effective error probability via `logadd`.
    substitution_rate_phred: f64,
    /// Seed k-mer length, needed by the missing-locus probability model to
    /// reason about how many independent exact-match windows a read offers.
    seed_k: usize,
    /// Total reference size, used by `logentropy` as the "random sequence"
    /// search space (port of `self.genome.getgenomesize()`); a bigger
    /// genome needs a more surprising alignment to be called non-random.
    genome_size: f64,
}

impl MapqCalculator {
    pub fn new(sub_rate: f64, seed_k: usize, genome_size: usize) -> Self {
        let substitution_rate_phred = py_int_round(sub_rate.max(1e-30).ln() / DB);
        Self {
            substitution_rate_phred,
            seed_k,
            genome_size: (genome_size.max(1)) as f64,
        }
    }

    // ------------------------------------------------------------------
    // 1. Alignment likelihood  Pr(read | ref, pos)
    // ------------------------------------------------------------------

    /// Log-likelihood of observing `read` given that it truly aligns to
    /// `ref_seq` (starting at index 0) via `cigar`. See
    /// `alignment_log_likelihood_bases` for the per-base model.
    pub fn alignment_log_likelihood(&self, read: &Read, ref_seq: &[u8], cigar: &[CigarOp]) -> f64 {
        self.alignment_log_likelihood_bases(&read.sequence.bases, &read.qualities.scores, ref_seq, cigar)
    }

    /// Log-likelihood of observing `bases`/`qual` given that they truly align
    /// to `ref_seq` (starting at index 0) via `cigar`.
    ///
    /// Walking the CIGAR (rather than assuming an ungapped 1:1 correspondence
    /// between read and reference positions) matters as soon as there's an
    /// indel: past the indel, a naive linear walk compares every downstream
    /// base against the wrong reference position, tanking the likelihood of
    /// an otherwise-correct alignment.
    ///
    /// Each matched position's effective error probability combines
    /// sequencing error (from the Phred score) and true biological
    /// substitution (`substitution_rate`) via `logadd` in Phred space —
    /// `effective_q = logadd(base_quality, substitution_rate_phred)` —
    /// rather than this module's previous ad hoc
    /// `(1-e)*(1-substitution_rate)` multiplication. Inserted bases (present
    /// in the read but not the reference) have no reference base to compare
    /// against, so they're scored as a uniform-random draw; deleted bases
    /// (present in the reference but not the read) consume reference only
    /// and contribute no likelihood term.
    pub fn alignment_log_likelihood_bases(&self, bases: &[u8], qual: &[u8], ref_seq: &[u8], cigar: &[CigarOp]) -> f64 {
        let mut log_lk = 0.0;
        let mut qi = 0usize;
        let mut ri = 0usize;

        for op in cigar {
            match *op {
                CigarOp::Match(n) => {
                    for _ in 0..n {
                        if qi >= bases.len() || ri >= ref_seq.len() { break; }
                        let raw_q = qual.get(qi).copied().unwrap_or(30) as f64;
                        let eff_q = logadd(raw_q, self.substitution_rate_phred);
                        let e = 10f64.powf(-eff_q / 10.0);
                        let p_match = 1.0 - e;
                        let p = if bases[qi] == ref_seq[ri] { p_match } else { e / 3.0 };
                        log_lk += p.max(1e-300).ln();
                        qi += 1;
                        ri += 1;
                    }
                }
                CigarOp::Insertion(n) => {
                    for _ in 0..n {
                        if qi >= bases.len() { break; }
                        log_lk += (0.25f64).ln();
                        qi += 1;
                    }
                }
                CigarOp::Deletion(n) => {
                    ri += n as usize;
                }
                CigarOp::SoftClip(n) => {
                    qi += n as usize;
                }
                CigarOp::HardClip(_) => {}
            }
        }
        log_lk
    }

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
        let best_idx = log_lks.iter().enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
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

    /// `missing_locus_probability`, exposed for the same reason as
    /// `wrong_mapping_probability` above (the fast path also requires the
    /// estimated probability of not having found the correct candidate to
    /// be sufficiently low for both mates).
    pub fn missing_locus_probability_estimate(&self, read: &Read, repeat_mask: &[f64]) -> f64 {
        self.missing_locus_probability(read, repeat_mask)
    }

    /// Final MAPQ for the winning candidate, combining candidate-posterior
    /// confidence (`wrongmapping`, from `sum_exp`), missing-locus risk
    /// (`notfoundprob`), and random-sequence entropy (`entropy`, from the
    /// winning alignment's own CIGAR) via `logadd`:
    /// `posterior_phred = min(99, logadd(entropy, logadd(wrongmapping,
    /// notfoundprob)))` (see module docs for why a base-alignment-quality
    /// term isn't part of this chain). Returns `(mapq,
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
        // `missing_locus_probability`), converted to Phred ---
        let p_missing = self.missing_locus_probability(read, repeat_mask).max(1e-30);
        let notfound_phred = -10.0 * p_missing.log10();

        let combined = logadd(wrongmapping, notfound_phred);
        let posterior_phred = logadd(entropy, combined).clamp(0.0, 99.0);

        let mapq = posterior_phred.round() as u8;
        let posterior_prob = 10f64.powf(-posterior_phred / 10.0);
        (mapq, 1.0 - posterior_prob)
    }

    /// Combinatorial (not simple `-10log10(p)`) Phred estimate of how
    /// surprising this alignment would be against a *random* reference of
    /// size `self.genome_size`, following a random-match hypothesis test.
    ///
    /// This models the size of the space of "equally good" random
    /// alignments as
    ///
    ///   |S| ~= 2g * (b-1)^m * C(n,m) * C(l, 2i+d) * 30^d
    ///
    /// (`g`=genome size, factor 2 for strand; `n`=aligned positions, `m`=
    /// mismatches among them; `l`=read length; `i`,`d`=insertion/deletion
    /// *run* counts, each run consuming a length choice from 1..30; `C(l,
    /// 2i+d)` accounts for choosing where each indel run starts/ends within
    /// the read) and the probability of a random match this good as `P ~=
    /// |S| * b^-n`, giving `entropy_phred = -10log10(P) = 10n*log10(b) -
    /// 10log10(2g) - 10m*log10(b-1) - 10log10(C(n,m)) - 10log10(C(l,2i+d)) -
    /// 10d*log10(30)`. `b` is a GC-content-adjusted effective alphabet size
    /// (b=4 for uniform composition; b<4 for GC/AT-biased sequence, via
    /// `b = exp(-f*ln(f/2) - (1-f)*ln((1-f)/2))`, `f`=GC fraction), computed
    /// here from the reference window rather than a fixed constant.
    /// `min_q=10` is this module's default threshold for which aligned
    /// positions count towards `n`/`m`.
    fn log_entropy(&self, cigar: &[CigarOp], query_bases: &[u8], query_qual: &[u8], ref_window: &[u8]) -> f64 {
        const MIN_Q: u8 = 10;
        const MAX_INDEL_LEN: f64 = 30.0;

        let mut n: i64 = 0; // aligned (matched or mismatched), high-quality positions
        let mut m: i64 = 0; // mismatches among those
        let mut readlen: usize = 0;
        let mut refpos: usize = 0;
        let mut indel_runs: i64 = 0; // i + d: number of separate indel runs
        let mut deletion_runs: i64 = 0; // d
        let mut curtype: u8 = 0; // 0=match, 1=ins, 2=del

        for op in cigar {
            match *op {
                CigarOp::Insertion(count) => {
                    for _ in 0..count {
                        if curtype != 1 { indel_runs += 1; }
                        curtype = 1;
                        readlen += 1;
                    }
                }
                CigarOp::Deletion(count) => {
                    for _ in 0..count {
                        if curtype != 2 {
                            deletion_runs += 1;
                            indel_runs += 1;
                        }
                        curtype = 2;
                        refpos += 1;
                    }
                }
                CigarOp::Match(count) => {
                    for _ in 0..count {
                        if curtype != 0 { indel_runs += 1; }
                        curtype = 0;
                        if readlen >= query_bases.len() || refpos >= ref_window.len() { break; }
                        let q = query_qual.get(readlen).copied().unwrap_or(0);
                        let is_hq = q >= MIN_Q;
                        if is_hq { n += 1; }
                        let a = query_bases[readlen];
                        let r = ref_window[refpos];
                        if a != r && a != b'N' && r != b'N' && is_hq {
                            m += 1;
                        }
                        readlen += 1;
                        refpos += 1;
                    }
                }
                _ => {}
            }
        }
        let insertion_runs = indel_runs - deletion_runs;

        let gc_frac = gc_content(ref_window);
        let b = effective_alphabet_size(gc_frac);

        let two_g = 2.0 * self.genome_size;
        let match_term = py_int_round(10.0 * (n as f64 * b.ln() - two_g.ln()) / std::f64::consts::LN_10);
        let mismatch_term = py_int_round(10.0 * m as f64 * (b - 1.0).max(1e-6).ln() / std::f64::consts::LN_10);
        let position_choice_term = phred_n_choose_m(readlen as i64, 2 * insertion_runs + deletion_runs);
        let deletion_length_term = if deletion_runs > 0 {
            py_int_round(deletion_runs as f64 * 10.0 * MAX_INDEL_LEN.ln() / std::f64::consts::LN_10)
        } else {
            0.0
        };

        match_term - mismatch_term - phred_n_choose_m(n, m) - position_choice_term - deletion_length_term
    }

    /// P(the true locus produced no seed hit anywhere in the read): the
    /// probability that *every* independent candidate 15-mer window has 2 or
    /// more mutations (sequencing errors, combined with the biological
    /// substitution rate) — 0 or 1 mutations would have been caught by the
    /// exact-match or single-mismatch-neighbour search respectively (see
    /// `mapper::seeds`). Models reads whose 15-mers all get neighbour-probed:
    /// since the hash table only indexes every 5th genomic position, a
    /// true-locus window only has a chance of an exact hash hit once every 5
    /// read positions, and non-overlapping windows 15bp apart are treated as
    /// independent for tractability. For a single 15-mer window with
    /// per-base mutation probabilities `q_1..q_15`, writing `z =
    /// prod(1-q_i)` and `p1 = sum(q_i/(1-q_i))`, the probability of >= 2
    /// mutations is `1 - z*(1+p1)`. The true read/genome phase offset (which
    /// of the 5 possible frame-shifts aligns a window to an indexed
    /// position) is unknown, so the result is Bayesian-averaged over all 5
    /// offsets with a uniform prior.
    ///
    /// `repeat_mask[read_pos]` (from `mapper::seeds::find_seed_hits`, empty
    /// when unavailable — e.g. for mate-rescue alignments, which skip
    /// seeding entirely) is this project's probability estimate that the
    /// true reference 15-mer at that window is repetitive: 1.0 when the
    /// window's own exact 15-mer looked repetitive (a spurious hit in a
    /// repetitive region shouldn't be trusted even when variant-free), or,
    /// in the more refined case, the probability that a probed 1-mismatch
    /// neighbour really is the true sequence (derived from the read's own
    /// base quality there) when that neighbour turned out to be repetitive
    /// instead. Combined into each window's own miss probability via a
    /// probabilistic OR (`1 - (1-a)(1-b)`), since a window can miss either
    /// because of too many real mutations *or* because the true locus is
    /// repetitive, independently of each other.
    fn missing_locus_probability(&self, read: &Read, repeat_mask: &[f64]) -> f64 {
        let k = self.seed_k;
        let len = read.qualities.scores.len();
        if k == 0 || len < k {
            return 0.5;
        }

        let q: Vec<f64> = read.qualities.scores.iter()
            .map(|&raw_q| {
                let eff_q = logadd(raw_q as f64, self.substitution_rate_phred);
                10f64.powf(-eff_q / 10.0).min(1.0 - 1e-12)
            })
            .collect();

        let mut offset_means = 0.0;
        for phase in 0..k.min(5) {
            let mut log_p_all_windows_miss = 0.0;
            let mut start = phase;
            while start + k <= len {
                let repeat_p = repeat_mask.get(start).copied().unwrap_or(0.0);
                let window = &q[start..start + k];
                let z: f64 = window.iter().map(|&qi| 1.0 - qi).product();
                let p1: f64 = window.iter().map(|&qi| qi / (1.0 - qi)).sum();
                let p_two_or_more = (1.0 - z * (1.0 + p1)).clamp(0.0, 1.0);
                // Probabilistic OR: miss if too many real mutations, OR the
                // true locus is repetitive -- independent causes.
                let p_miss = (1.0 - (1.0 - p_two_or_more) * (1.0 - repeat_p)).clamp(1e-300, 1.0);
                log_p_all_windows_miss += p_miss.ln();
                start += k; // next non-overlapping window
            }
            offset_means += log_p_all_windows_miss.exp();
        }

        (offset_means / k.min(5) as f64).min(0.99)
    }
}

/// Fraction of G/C bases in `seq` (N and other ambiguity codes excluded from
/// both numerator and denominator), used to derive the GC-adjusted effective
/// alphabet size `b` for the random-match model (`log_entropy`).
fn gc_content(seq: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut total = 0usize;
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => { gc += 1; total += 1; }
            b'A' | b'T' => { total += 1; }
            _ => {}
        }
    }
    if total == 0 { 0.5 } else { gc as f64 / total as f64 }
}

/// GC-content-adjusted effective alphabet size,
/// `b = exp(-f*ln(f/2) - (1-f)*ln((1-f)/2))`.
/// `b=4` at `f=0.5` (uniform composition); GC- or AT-biased composition
/// (`f` away from 0.5) lowers `b`, since a biased random sequence is more
/// likely to coincidentally match by chance.
fn effective_alphabet_size(gc_frac: f64) -> f64 {
    let f = gc_frac.clamp(1e-6, 1.0 - 1e-6);
    (-f * (f / 2.0).ln() - (1.0 - f) * ((1.0 - f) / 2.0).ln()).exp()
}
