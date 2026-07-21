//! Paired-end read mapping.
//!
//! Each mate is seeded and clustered into its own candidate list
//! (`SingleEndMapper::map_candidates`) independently. From there, this
//! module:
//!
//! - Build a **shortlist** of candidate loci for each mate covering 99.9% of
//!   its own single-read posterior mass (clamped to 3..=20 locations), and
//!   for every shortlisted locus, align the *other* mate directly against
//!   the window the library insert-size distribution implies (reusing
//!   `SingleEndMapper::align_at_window`, the same primitive plain mate
//!   rescue uses).
//! - Score every resulting candidate **pair** with the joint likelihood
//!   `L(r1,r2,x,y) = Pr(r1|x)·Pr(r2|y)·Pd(y-x)`, where `Pd` is the
//!   insert-size Gaussian floored by a structural-variant prior
//!   (`--svprior`) for pairs whose orientation/separation looks nothing like
//!   a normal library fragment — this is what lets a read pair spanning a
//!   real deletion or a repeat-vs-repeat mismatch still get scored sensibly
//!   instead of just being penalized to nothing.
//! - Pick the best-scoring pair by softmax over all candidate pairs, and
//!   assign MAPQ accordingly: if the winning pair is simply each mate's own
//!   independent best pick, report each mate's own single-end MAPQ; if
//!   pairing was needed to pick out a locus that wasn't already the winner
//!   on its own (e.g. a mate landing in a repeat where its true copy is only
//!   distinguishable via its partner's position), both mates report the
//!   confidence of whichever one *was* already unambiguous on its own — see
//!   `mapq_for_winning_pair` for the exact interpretation, including the
//!   case where neither side's own best candidate survived.
//!
//! When one mate has no candidates of its own at all, this degenerates
//! naturally to plain mate rescue (a one-candidate "shortlist" built from the
//! other mate's own best pick).

use crate::types::*;
use crate::align::cigar_reference_span;
use crate::index::{GenomeIndex, HashTable};
use crate::mapper::single::{SingleEndMapper, CandidateAlignment, unmapped_record};
use crate::mapq::MapqCalculator;
use anyhow::Result;

/// Number of standard deviations from the mean a pair's separation must fall
/// within to be flagged "proper", and the search radius used both for plain
/// mate rescue and for shortlist cross-alignment. Five sigma is a
/// conventional statistical threshold for "this is not plausibly due to
/// chance" (the same convention used for discovery-level significance in
/// physics), used here as a generously permissive cutoff so that legitimate
/// biological variance in fragment size isn't mistaken for a discordant
/// pair.
const PROPER_PAIR_SD: f64 = 5.0;

/// How close (in standard deviations of the current insert-size estimate)
/// the naive best-best pair must be for the fast path below to even
/// consider it: the best locations for the single reads are "close
/// together" when their implied insert size is within 4 standard deviations
/// of the mean. Not the same constant as `PROPER_PAIR_SD` (5.0): that one
/// gates the "proper pair" SAM flag and the mate-rescue/shortlist search
/// radius, this one gates a different, narrower question.
const FAST_PATH_CLOSE_SD: f64 = 4.0;

/// Posterior-probability-of-being-wrong threshold both mates' independent
/// best pick must clear (in *both* the wrongmapping and missing-locus
/// senses) for the fast path to skip shortlist construction: a posterior
/// probability of less than 1% of having mapped incorrectly, and again less
/// than 1% of having missed the correct candidate.
const FAST_PATH_MAX_WRONG_PROB: f64 = 0.01;

/// Cumulative single-read posterior mass a mate's candidate shortlist must
/// cover, and the `[min, max]` candidate count it's clamped to: the
/// locations that together constitute 99.9% of the single-read posterior
/// mapping probability are extracted, up to a maximum of 20 locations, and
/// subject to a minimum of 3.
const SHORTLIST_POSTERIOR_MASS: f64 = 0.999;
const SHORTLIST_MIN: usize = 3;
const SHORTLIST_MAX: usize = 20;

pub struct PairedEndMapper<'a> {
    single_mapper: SingleEndMapper<'a>,
    insert_model: InsertSizeDistribution,
    /// Natural-log structural-variant prior probability (`--svprior`, Phred
    /// scale on input; default is Phred 55, i.e. probability 3e-6).
    sv_prior_log_prob: f64,
}

impl<'a> PairedEndMapper<'a> {
    pub fn new(
        genome: &'a GenomeIndex,
        hash_table: &'a HashTable,
        options: &'a MapperOptions,
    ) -> Self {
        let single = SingleEndMapper::new(genome, hash_table, options);
        // Seed the insert-size prior from the CLI options; online learning
        // then refines it from observed proper pairs.
        let insert = InsertSizeDistribution::with_params(
            options.insert_size_mean,
            options.insert_size_sd,
        );
        let sv_prior_log_prob = -options.sv_prior_phred / 10.0 * std::f64::consts::LN_10;
        Self { single_mapper: single, insert_model: insert, sv_prior_log_prob }
    }

    pub fn map_pair(&mut self, pair: &ReadPair) -> Result<(AlignmentRecord, AlignmentRecord)> {
        let insert_snapshot = self.insert_model.clone();
        let (r1, r2, observed_insert_size) = self.map_pair_readonly(pair, &insert_snapshot);
        if let Some(insert_size) = observed_insert_size {
            self.insert_model.update(insert_size);
        }
        Ok((r1, r2))
    }

    /// Batch version of `map_pair`: processes every pair in `pairs` against
    /// one frozen snapshot of the insert-size model, in parallel across
    /// `rayon` worker threads, then serially folds whatever proper,
    /// independently-confident observations came out of the batch into the
    /// model afterwards.
    ///
    /// This is what actually lets paired-end mapping use more than one core:
    /// `map_pair` itself can't be parallelized directly across pairs because
    /// each call both *reads* the current insert-size estimate (to size
    /// search windows and score candidate pairs) and potentially *writes* to
    /// it, and those writes need to seralize somewhere. Freezing the model
    /// for the batch sidesteps that entirely -- no shared mutable state
    /// during the parallel phase, so no locking, at the cost of the model
    /// being up to one batch "stale" (bounded by `batch.len()`, and
    /// negligible once the model has converged from a few thousand
    /// observations, which is the common case well before a batch boundary
    /// matters).
    pub fn map_pairs_batch(&mut self, pairs: &[ReadPair]) -> Vec<(AlignmentRecord, AlignmentRecord)> {
        use rayon::prelude::*;

        let insert_snapshot = self.insert_model.clone();
        let results: Vec<((AlignmentRecord, AlignmentRecord), Option<f64>)> = pairs
            .par_iter()
            .map(|pair| {
                let (r1, r2, observed) = self.map_pair_readonly(pair, &insert_snapshot);
                ((r1, r2), observed)
            })
            .collect();

        for (_, observed) in &results {
            if let Some(insert_size) = observed {
                self.insert_model.update(*insert_size);
            }
        }

        results.into_iter().map(|(rec, _)| rec).collect()
    }

    /// The actual mapping logic, parameterized by an explicit insert-size
    /// snapshot instead of reading `self.insert_model` directly, and
    /// returning (rather than immediately applying) the insert-size
    /// observation this pair would contribute to online learning --
    /// `&self`, not `&mut self`, so `map_pairs_batch` can call it from
    /// multiple `rayon` threads concurrently.
    fn map_pair_readonly(&self, pair: &ReadPair, insert_model: &InsertSizeDistribution) -> (AlignmentRecord, AlignmentRecord, Option<f64>) {
        // The two mates' seed+cluster+align passes are fully independent of
        // each other -- profiling against real E. coli data found this step
        // alone at ~775us/pair, so running the two sides on separate
        // `rayon` worker threads is free parallelism, on top of whatever
        // parallelism the caller (`map_pairs_batch`) is already getting
        // across pairs.
        let (cands1, cands2) = rayon::join(
            || self.single_mapper.map_candidates(&pair.read1),
            || self.single_mapper.map_candidates(&pair.read2),
        );
        let (cands1, repeat_mask1) = cands1;
        let (cands2, repeat_mask2) = cands2;

        let (mut r1, mut r2, independently_confident) = if cands1.is_empty() && cands2.is_empty() {
            (unmapped_record(&pair.read1), unmapped_record(&pair.read2), false)
        } else {
            self.joint_map(pair, &cands1, &repeat_mask1, &cands2, &repeat_mask2, insert_model)
        };

        r1.is_paired = true;
        r2.is_paired = true;
        r1.read_number = Some(1);
        r2.read_number = Some(2);

        let r1_mapped = r1.position != 0;
        let r2_mapped = r2.position != 0;
        r1.mate_unmapped = !r2_mapped;
        r2.mate_unmapped = !r1_mapped;
        r1.mate_strand = if r2_mapped { Some(r2.strand) } else { None };
        r2.mate_strand = if r1_mapped { Some(r1.strand) } else { None };

        let mut observed_insert_size = None;

        if r1_mapped && r2_mapped && r1.contig_name == r2.contig_name {
            // TLEN (SAM spec): distance from the leftmost mapped base of the
            // upstream mate to the rightmost mapped base of the downstream
            // mate -- NOT just the difference between their start
            // positions, which silently drops the downstream mate's own
            // reference-consumed span. Found by comparing against real E.
            // coli data: the old start-to-start calculation undercounted by
            // very close to one read length (100bp reads, ~100bp gap between
            // the previously-learned insert-size mean and the corrected one
            // on the same reads -- 330bp vs. 430bp), which then fed a
            // systematically-too-short prior into every
            // subsequent proper-pair check and rescue/shortlist search
            // window.
            let (upstream, downstream) = if r1.position <= r2.position { (&r1, &r2) } else { (&r2, &r1) };
            let downstream_end = downstream.position + cigar_reference_span(&downstream.cigar) - 1;
            let insert_size = (downstream_end as i64 - upstream.position as i64 + 1).max(0);

            let within_distance = (insert_size as f64 - insert_model.mean).abs()
                < PROPER_PAIR_SD * insert_model.sd;
            let correctly_oriented = is_fr_orientation(&r1, &r2);
            let is_proper = within_distance && correctly_oriented;

            r1.is_proper_pair = is_proper;
            r2.is_proper_pair = is_proper;
            r1.mate_contig = Some(r2.contig_name.clone());
            r2.mate_contig = Some(r1.contig_name.clone());
            r1.mate_position = Some(r2.position);
            r2.mate_position = Some(r1.position);
            // Signed per SAM convention: positive for the upstream mate,
            // negative for the downstream one (matches `r1`/`r2` regardless
            // of which turned out to be upstream).
            r1.insert_size = Some(if r1.position <= r2.position { insert_size } else { -insert_size });
            r2.insert_size = Some(if r2.position <= r1.position { insert_size } else { -insert_size });

            // Only learn from pairs neither mate needed cross-alignment help
            // to place: feeding a rescued/shortlist-disambiguated mate's
            // (possibly still slightly off) position back into the model
            // that *drives* that same rescue's search window is a feedback
            // loop that could in principle drag the mean away from the true
            // value over a run, independent of the TLEN bug above.
            if is_proper && independently_confident {
                observed_insert_size = Some(insert_size as f64);
            }
        }

        (r1, r2, observed_insert_size)
    }

    /// The core shortlist + cross-alignment + joint-posterior procedure.
    /// Precondition: `cands1` and `cands2` are not *both* empty (the caller
    /// short-circuits that case before calling this).
    fn joint_map(
        &self,
        pair: &ReadPair,
        cands1: &[CandidateAlignment],
        repeat_mask1: &[f64],
        cands2: &[CandidateAlignment],
        repeat_mask2: &[f64],
        insert_model: &InsertSizeDistribution,
    ) -> (AlignmentRecord, AlignmentRecord, bool) {
        let mq = self.single_mapper.mapq_calc();
        let lk1: Vec<f64> = cands1.iter().map(|c| c.log_lk).collect();
        let lk2: Vec<f64> = cands2.iter().map(|c| c.log_lk).collect();
        let best1_idx = argmax_index(&lk1);
        let best2_idx = argmax_index(&lk2);

        // Fast path: when both mates' own independent best picks are already
        // close together, each individually unambiguous, and each
        // individually unlikely to have missed the true candidate, report
        // that pair directly rather than building a shortlist at all -- the
        // expensive cross-alignment machinery below exists specifically for
        // the cases that *don't* clear this bar.
        // Skipping it here isn't a speed/fidelity tradeoff: running it
        // unconditionally (as this project did before this fast path was
        // added) was itself the fidelity gap, confirmed by comparing timing
        // against real E. coli data, where it accounted for the large
        // majority of paired-mapping wall time despite most real pairs
        // being exactly the unambiguous case this fast path is for.
        if let (Some(i1), Some(i2)) = (best1_idx, best2_idx) {
            let (sum_exp1, sum_exp2) = (mq.select_best(&lk1).1, mq.select_best(&lk2).1);
            let wrong1 = mq.wrong_mapping_probability(sum_exp1);
            let wrong2 = mq.wrong_mapping_probability(sum_exp2);
            if wrong1 < FAST_PATH_MAX_WRONG_PROB && wrong2 < FAST_PATH_MAX_WRONG_PROB {
                let missing1 = mq.missing_locus_probability_estimate(&pair.read1, repeat_mask1);
                let missing2 = mq.missing_locus_probability_estimate(&pair.read2, repeat_mask2);
                if missing1 < FAST_PATH_MAX_WRONG_PROB && missing2 < FAST_PATH_MAX_WRONG_PROB {
                    let (bx, by) = (&cands1[i1], &cands2[i2]);
                    let gap = implied_insert_size(bx, by) as f64;
                    let close = bx.contig_id == by.contig_id
                        && is_fr_pair(bx, by)
                        && (gap - insert_model.mean).abs() <= FAST_PATH_CLOSE_SD * insert_model.sd;
                    if close {
                        let (mapq1, post1, mapq2, post2) = self.mapq_for_winning_pair(
                            pair, bx, by, 1.0, &lk1, &lk2, repeat_mask1, repeat_mask2, true, true,
                        );
                        return (
                            self.single_mapper.to_record(&pair.read1, bx, mapq1, post1),
                            self.single_mapper.to_record(&pair.read2, by, mapq2, post2),
                            true,
                        );
                    }
                }
            }
        }

        let mut raw_pairs: Vec<(CandidateAlignment, CandidateAlignment)> = Vec::new();

        match (cands1.is_empty(), cands2.is_empty()) {
            (false, false) => {
                // Both sides have candidates: the "obvious" independent
                // best-best pick, plus shortlist-driven cross-alignment
                // extended from each side in turn.
                if let (Some(i1), Some(i2)) = (best1_idx, best2_idx) {
                    raw_pairs.push((cands1[i1].clone(), cands2[i2].clone()));
                }
                let post1 = mq.posterior_fractions(&lk1);
                for i in shortlist_indices(&post1) {
                    if let Some(y) = self.align_mate_near(&pair.read2, &cands1[i], insert_model) {
                        raw_pairs.push((cands1[i].clone(), y));
                    }
                }
                let post2 = mq.posterior_fractions(&lk2);
                for j in shortlist_indices(&post2) {
                    if let Some(x) = self.align_mate_near(&pair.read1, &cands2[j], insert_model) {
                        raw_pairs.push((x, cands2[j].clone()));
                    }
                }
            }
            (true, false) => {
                // Plain mate rescue: read1 has nothing of its own, anchor on
                // read2's own best candidate.
                let j = best2_idx.expect("cands2 non-empty implies a best index");
                if let Some(x) = self.align_mate_near(&pair.read1, &cands2[j], insert_model) {
                    raw_pairs.push((x, cands2[j].clone()));
                }
            }
            (false, true) => {
                let i = best1_idx.expect("cands1 non-empty implies a best index");
                if let Some(y) = self.align_mate_near(&pair.read2, &cands1[i], insert_model) {
                    raw_pairs.push((cands1[i].clone(), y));
                }
            }
            (true, true) => unreachable!("caller guarantees not both candidate lists are empty"),
        }

        if raw_pairs.is_empty() {
            // Nothing could be cross-aligned (e.g. the only anchor's implied
            // window ran off the end of its contig) -- fall back to each
            // side's independent single-end result.
            return (
                independent_record(&self.single_mapper, &pair.read1, cands1, repeat_mask1),
                independent_record(&self.single_mapper, &pair.read2, cands2, repeat_mask2),
                true,
            );
        }

        // Dedup identical (locus, locus) pairs -- the independent best-best
        // pick very often reappears verbatim from a shortlist extension.
        let mut pairs: Vec<(CandidateAlignment, CandidateAlignment)> = Vec::new();
        for (x, y) in raw_pairs {
            if !pairs.iter().any(|(px, py)| same_locus(px, &x) && same_locus(py, &y)) {
                pairs.push((x, y));
            }
        }

        let joint_lks: Vec<f64> = pairs.iter().map(|(x, y)| self.pair_joint_log_lk(x, y, insert_model)).collect();
        let (best_idx, sum_exp) = mq.select_best(&joint_lks);
        let (bx, by) = pairs[best_idx].clone();

        let bx_is_own_best = best1_idx.is_some_and(|i| same_locus(&bx, &cands1[i]));
        let by_is_own_best = best2_idx.is_some_and(|j| same_locus(&by, &cands2[j]));

        let (mapq1, post1, mapq2, post2) = self.mapq_for_winning_pair(
            pair, &bx, &by, sum_exp,
            &lk1, &lk2, repeat_mask1, repeat_mask2,
            bx_is_own_best, by_is_own_best,
        );

        (
            self.single_mapper.to_record(&pair.read1, &bx, mapq1, post1),
            self.single_mapper.to_record(&pair.read2, &by, mapq2, post2),
            bx_is_own_best && by_is_own_best,
        )
    }

    /// Assigns final MAPQ/posterior to the winning candidate pair. The
    /// intended combination rule: the posterior mapping quality is the
    /// product of the single-end mapping qualities when the top-scoring
    /// single-end hits are what got selected as the pair, or the single-end
    /// posterior of the anchoring read in other cases. This project's
    /// concrete implementation of that rule:
    /// - **Both mates' own independent best pick won** (`bx_is_own_best &&
    ///   by_is_own_best`): each mate reports its *own* single-end MAPQ
    ///   (computed from its own candidate set) -- pairing here only
    ///   confirmed the two picks are consistent, it didn't have to
    ///   disambiguate anything, so "product of the single-end qualities"
    ///   reduces to each read just keeping its own.
    /// - **Exactly one mate's own best pick won**: that mate is the anchor
    ///   the other one's true locus was picked out *via* (the other mate's
    ///   own top pick lost to a cross-aligned candidate). Both mates report
    ///   the anchor's own single-end confidence, since the rescued mate's
    ///   placement is only as trustworthy as the anchor it was found relative
    ///   to.
    /// - **Neither mate's own best pick won** (both were out-competed by a
    ///   cross-aligned candidate at a different locus): there's no single
    ///   obviously-correct rule for this case. This project's choice is to
    ///   score each mate against its own candidate set using the winning
    ///   *pair's* joint `sum_exp` (so the pairing evidence that changed the
    ///   outcome still counts).
    #[allow(clippy::too_many_arguments)]
    fn mapq_for_winning_pair(
        &self,
        pair: &ReadPair,
        bx: &CandidateAlignment,
        by: &CandidateAlignment,
        pair_sum_exp: f64,
        lk1: &[f64],
        lk2: &[f64],
        repeat_mask1: &[f64],
        repeat_mask2: &[f64],
        bx_is_own_best: bool,
        by_is_own_best: bool,
    ) -> (u8, f64, u8, f64) {
        let mq = self.single_mapper.mapq_calc();

        if bx_is_own_best && by_is_own_best {
            let (_, sum_exp1) = mq.select_best(lk1);
            let (_, sum_exp2) = mq.select_best(lk2);
            let (mapq1, post1) = mq.finalize_mapq(sum_exp1, &pair.read1, &bx.cigar, &bx.out_seq, &bx.out_qual, &bx.ref_window, repeat_mask1);
            let (mapq2, post2) = mq.finalize_mapq(sum_exp2, &pair.read2, &by.cigar, &by.out_seq, &by.out_qual, &by.ref_window, repeat_mask2);
            (mapq1, post1, mapq2, post2)
        } else if bx_is_own_best != by_is_own_best {
            let (anchor_sum_exp, anchor_read, anchor_cand, anchor_mask) = if bx_is_own_best {
                (mq.select_best(lk1).1, &pair.read1, bx, repeat_mask1)
            } else {
                (mq.select_best(lk2).1, &pair.read2, by, repeat_mask2)
            };
            let (anchor_mapq, anchor_post) = mq.finalize_mapq(
                anchor_sum_exp, anchor_read, &anchor_cand.cigar, &anchor_cand.out_seq,
                &anchor_cand.out_qual, &anchor_cand.ref_window, anchor_mask,
            );
            (anchor_mapq, anchor_post, anchor_mapq, anchor_post)
        } else {
            let (mapq1, post1) = mq.finalize_mapq(pair_sum_exp, &pair.read1, &bx.cigar, &bx.out_seq, &bx.out_qual, &bx.ref_window, repeat_mask1);
            let (mapq2, post2) = mq.finalize_mapq(pair_sum_exp, &pair.read2, &by.cigar, &by.out_seq, &by.out_qual, &by.ref_window, repeat_mask2);
            (mapq1, post1, mapq2, post2)
        }
    }

    /// Align `mate_read` against the window `anchor`'s position and the
    /// current insert-size estimate imply, assuming FR orientation (a
    /// forward-strand anchor implies the mate should be found downstream on
    /// the reverse strand, and vice versa). This is the one primitive both
    /// plain mate rescue and shortlist cross-alignment reduce to.
    ///
    /// `insert_model.mean` is a TLEN-style outer-span estimate (the model is
    /// trained on `map_pair_readonly`'s TLEN calculation, `downstream_end -
    /// upstream.position + 1`) -- so predicting the *mate's own start*
    /// position needs to account for the mate's own length, not just add
    /// `mean` directly to the anchor's start (which would place the window
    /// about one mate-length too far downstream/not far enough upstream; see
    /// `implied_insert_size`'s doc comment, found by code review).
    fn align_mate_near(&self, mate_read: &Read, anchor: &CandidateAlignment, insert_model: &InsertSizeDistribution) -> Option<CandidateAlignment> {
        let mean = insert_model.mean.max(0.0) as usize;
        let mate_len = mate_read.sequence.len();
        let search_radius = (PROPER_PAIR_SD * insert_model.sd).ceil() as usize + mate_len;

        let (expected_position, mate_strand) = if anchor.strand.is_forward() {
            // Anchor is upstream; the downstream mate's own start is the
            // implied fragment end minus its own length.
            (anchor.position + mean.saturating_sub(mate_len), Strand::Reverse)
        } else {
            // Anchor is downstream; the upstream mate's start is the
            // anchor's own rightmost reference base minus the implied
            // fragment length.
            let anchor_end = anchor.position + cigar_reference_span(&anchor.cigar);
            (anchor_end.saturating_sub(mean), Strand::Forward)
        };

        self.single_mapper.align_at_window(mate_read, anchor.contig_id, expected_position, search_radius, mate_strand)
    }

    /// Joint likelihood `L(r1,r2,x,y) = Pr(r1|x)·Pr(r2|y)·Pd(y-x)` (section
    /// 1.12 / main text) in log space: the two mates' own alignment
    /// log-likelihoods plus the insert-size term from `pair_distance_log_lk`.
    fn pair_joint_log_lk(&self, r1: &CandidateAlignment, r2: &CandidateAlignment, insert_model: &InsertSizeDistribution) -> f64 {
        r1.log_lk + r2.log_lk + self.pair_distance_log_lk(r1, r2, insert_model)
    }

    /// `Pd`: the insert-size Gaussian log-likelihood for an FR-oriented pair
    /// on the same contig, floored by the structural-variant prior
    /// (`--svprior`) — "close" pairs use whichever term is larger, i.e. the
    /// crossover point where the Gaussian's likelihood drops below the
    /// prior; implemented here simply as
    /// `max(gaussian_log_lk, sv_prior_log_prob)`, which has exactly that
    /// crossover behavior. Any other-contig or non-FR-oriented pair gets the
    /// SV prior outright, since the Gaussian model doesn't describe those
    /// configurations at all.
    fn pair_distance_log_lk(&self, r1: &CandidateAlignment, r2: &CandidateAlignment, insert_model: &InsertSizeDistribution) -> f64 {
        if r1.contig_id != r2.contig_id || !is_fr_pair(r1, r2) {
            return self.sv_prior_log_prob;
        }
        // `implied_insert_size` is always >= 0 by construction (the
        // downstream candidate's own end is never before its start, and its
        // start is never before the upstream candidate's start), so no
        // clamping is needed here.
        let gap = implied_insert_size(r1, r2) as f64;
        insert_model.log_likelihood(gap).max(self.sv_prior_log_prob)
    }
}

fn independent_record(single_mapper: &SingleEndMapper, read: &Read, cands: &[CandidateAlignment], repeat_mask: &[f64]) -> AlignmentRecord {
    if cands.is_empty() {
        return unmapped_record(read);
    }
    let mq: &MapqCalculator = single_mapper.mapq_calc();
    let lk: Vec<f64> = cands.iter().map(|c| c.log_lk).collect();
    let (best_idx, sum_exp) = mq.select_best(&lk);
    let best = &cands[best_idx];
    let (mapq, post) = mq.finalize_mapq(sum_exp, read, &best.cigar, &best.out_seq, &best.out_qual, &best.ref_window, repeat_mask);
    single_mapper.to_record(read, best, mapq, post)
}

fn argmax_index(lk: &[f64]) -> Option<usize> {
    lk.iter().enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
}

fn same_locus(a: &CandidateAlignment, b: &CandidateAlignment) -> bool {
    a.contig_id == b.contig_id && a.position == b.position && a.strand == b.strand
}

/// TLEN-style implied insert size between two candidate alignments: distance
/// from the leftmost mapped base of the upstream one to the rightmost mapped
/// base of the downstream one (via `cigar_reference_span`, so a downstream
/// candidate with an indel is measured correctly), matching the SAM TLEN
/// convention `map_pair_readonly` reports and the insert-size model is
/// trained on. Must be used everywhere a candidate pair's separation is
/// compared against `insert_model.mean`/`.sd` -- comparing against a plain
/// `position` difference instead systematically understates the true
/// separation by about one mate's read length (the downstream mate's own
/// span), since `position` is only its leftmost coordinate. Found by code
/// review: `pair_distance_log_lk`, `align_mate_near`, and the fast-path
/// "close" check all used to compute a plain start-to-start gap while the
/// model they compared it against was trained on this TLEN-style value.
fn implied_insert_size(a: &CandidateAlignment, b: &CandidateAlignment) -> i64 {
    let (upstream, downstream) = if a.position <= b.position { (a, b) } else { (b, a) };
    let downstream_end = downstream.position + cigar_reference_span(&downstream.cigar);
    downstream_end as i64 - upstream.position as i64
}

/// FR (forward-reverse) orientation for a candidate pair: the one at the
/// lower coordinate must be on the forward strand and the one at the higher
/// coordinate on the reverse strand (standard Illumina paired-end library
/// layout).
fn is_fr_pair(a: &CandidateAlignment, b: &CandidateAlignment) -> bool {
    let (upstream, downstream) = if a.position <= b.position { (a, b) } else { (b, a) };
    upstream.strand.is_forward() && !downstream.strand.is_forward()
}

/// Given a caller-supplied posterior distribution (already normalized,
/// summing to 1), return the indices (best-first) of the smallest prefix
/// that both covers `SHORTLIST_POSTERIOR_MASS` and contains at least
/// `SHORTLIST_MIN` entries, capped at `SHORTLIST_MAX`.
fn shortlist_indices(posterior: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..posterior.len()).collect();
    order.sort_by(|&a, &b| posterior[b].partial_cmp(&posterior[a]).unwrap());

    let mut cum = 0.0;
    let mut take = 0;
    for &i in &order {
        if take >= SHORTLIST_MAX { break; }
        if take >= SHORTLIST_MIN && cum >= SHORTLIST_POSTERIOR_MASS { break; }
        cum += posterior[i];
        take += 1;
    }
    order.into_iter().take(take).collect()
}

/// FR (forward-reverse) orientation: the mate mapped to the lower coordinate
/// must be on the forward strand and the mate mapped to the higher
/// coordinate must be on the reverse strand (standard Illumina paired-end
/// library layout). FF/RR/RF layouts are not "proper" here.
fn is_fr_orientation(r1: &AlignmentRecord, r2: &AlignmentRecord) -> bool {
    let (upstream, downstream) = if r1.position <= r2.position { (r1, r2) } else { (r2, r1) };
    upstream.strand.is_forward() && !downstream.strand.is_forward()
}
