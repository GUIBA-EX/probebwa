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
//! - Combine single-read and pair hypotheses in Phred space. Paired SAM
//!   records use pair confidence when the pair passes the random-match test.
//!   The probability primitives are oracle-tested; candidate search and
//!   slow-path rescue still differ from Stampy.
//!
//! When one mate has no candidates of its own at all, this degenerates
//! naturally to plain mate rescue (a one-candidate "shortlist" built from the
//! other mate's own best pick).

use crate::types::*;
use crate::align::cigar_reference_span;
use crate::index::{GenomeIndex, HashTable};
use crate::mapper::single::{SingleEndMapper, CandidateAlignment, unmapped_record};
use crate::mapq::MapqCalculator;
use crate::mapq::bayesian::{logadd, logadd_unclamped};
use anyhow::Result;

/// Stampy marks proper pairs strictly inside three standard deviations.
const PROPER_PAIR_SD: f64 = 3.0;
/// Search radius is independent of the proper-pair reporting threshold.
const RESCUE_SEARCH_SD: f64 = 5.0;

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
        let single = SingleEndMapper::new_paired(genome, hash_table, options);
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

    /// Seeding does not depend on the library model. Run it in parallel,
    /// then make pair decisions and update that model in input order. A
    /// frozen model per batch changes the output at arbitrary chunk borders.
    pub fn map_pairs_batch(&mut self, pairs: &[ReadPair]) -> Vec<(AlignmentRecord, AlignmentRecord)> {
        use rayon::prelude::*;
        let candidates: Vec<_> = pairs.par_iter().map(|pair| (
            self.single_mapper.map_candidates(&pair.read1),
            self.single_mapper.map_candidates(&pair.read2),
        )).collect();
        pairs.iter().zip(candidates).map(|(pair, (first, second))| {
            let (a, b, observed) = self.finish_pair(pair, &self.insert_model, first, second);
            if let Some(distance) = observed { self.insert_model.update(distance); }
            (a, b)
        }).collect()
    }

    fn map_pair_readonly(&self, pair: &ReadPair, insert_model: &InsertSizeDistribution) -> (AlignmentRecord, AlignmentRecord, Option<f64>) {
        let (first, second) = rayon::join(
            || self.single_mapper.map_candidates(&pair.read1),
            || self.single_mapper.map_candidates(&pair.read2),
        );
        self.finish_pair(pair, insert_model, first, second)
    }

    fn finish_pair(&self, pair: &ReadPair, insert_model: &InsertSizeDistribution,
        first: (Vec<CandidateAlignment>, Vec<f64>), second: (Vec<CandidateAlignment>, Vec<f64>),
    ) -> (AlignmentRecord, AlignmentRecord, Option<f64>) {
        let (cands1, repeat_mask1) = first;
        let (cands2, repeat_mask2) = second;

        let (mut r1, mut r2, observed_insert_size) = if cands1.is_empty() && cands2.is_empty() {
            (unmapped_record(&pair.read1), unmapped_record(&pair.read2), None)
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

        // Mate coordinates are independent of contig equality and proper-pair status.
        if r2_mapped {
            r1.mate_contig = Some(r2.contig_name.clone());
            r1.mate_position = Some(r2.position);
        }
        if r1_mapped {
            r2.mate_contig = Some(r1.contig_name.clone());
            r2.mate_position = Some(r1.position);
        }

        if r1_mapped && !r2_mapped {
            r1.mate_contig=Some(r1.contig_name.clone()); r1.mate_position=Some(r1.position);
            r2.sam_tags.extend([format!("MQ:i:{}",r1.mapq),format!("XQ:i:{}",(-10.0*r1.alignment_score/std::f64::consts::LN_10).round() as i32)]);
        } else if r2_mapped && !r1_mapped {
            r2.mate_contig=Some(r2.contig_name.clone()); r2.mate_position=Some(r2.position);
            r1.sam_tags.extend([format!("MQ:i:{}",r2.mapq),format!("XQ:i:{}",(-10.0*r2.alignment_score/std::f64::consts::LN_10).round() as i32)]);
        }


        if r1_mapped && r2_mapped && r1.contig_name == r2.contig_name {
            // Native SAM uses each mate's strand-specific outer endpoint.
            // Overlapping/RF pairs cannot be ordered merely by their starts.
            let endpoint = |r: &AlignmentRecord| if r.strand == Strand::Forward { r.position as i64 }
                else { (r.position + cigar_reference_span(&r.cigar)) as i64 };
            let signed_size = endpoint(&r2) - endpoint(&r1);
            let insert_size = signed_size.abs();
            let correctly_oriented = is_fr_orientation(&r1,&r2);
            let within_distance = (insert_size as f64 - 1.0 - insert_model.mean).abs() < PROPER_PAIR_SD * insert_model.sd;
            r1.is_proper_pair=within_distance && correctly_oriented;
            r2.is_proper_pair=r1.is_proper_pair;
            r1.insert_size=Some(signed_size); r2.insert_size=Some(-signed_size);


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
    ) -> (AlignmentRecord, AlignmentRecord, Option<f64>) {
        let mq = self.single_mapper.mapq_calc();
        let lk1: Vec<f64> = cands1.iter().map(|c| c.coarse_score as f64 * (-std::f64::consts::LN_10 / 10.0)).collect();
        let lk2: Vec<f64> = cands2.iter().map(|c| c.coarse_score as f64 * (-std::f64::consts::LN_10 / 10.0)).collect();
        let best1_idx = argmax_index(&lk1);
        let best2_idx = argmax_index(&lk2);

        if let (Some(i1),Some(i2)) = (best1_idx,best2_idx) {
            let sum1=mq.select_best(&lk1).1; let sum2=mq.select_best(&lk2).1;
            let (_,missing1,_)=mq.mapping_evidence(sum1,&pair.read1,repeat_mask1);
            let (_,missing2,_)=mq.mapping_evidence(sum2,&pair.read2,repeat_mask2);
            let (bx,by)=(&cands1[i1],&cands2[i2]);
            let competing = cands1.iter().enumerate().any(|(i,c)| i != i1 && c.coarse_score < bx.coarse_score+30
                && c.contig_id == by.contig_id && implied_insert_size(c,by) < 10000)
                || cands2.iter().enumerate().any(|(i,c)| i != i2 && c.coarse_score < by.coarse_score+30
                && c.contig_id == bx.contig_id && implied_insert_size(bx,c) < 10000);
            if self.distance_phred(bx,by,insert_model) < 40 && sum1+sum2 < 2.02 && (missing1 > 40.0 || missing2 > 40.0) && !competing {
                let evidence=self.mapq_for_winning_pair(pair,bx,by,1.0,&lk1,&lk2,repeat_mask1,repeat_mask2,true,true,insert_model,true);
                let (a,b)=self.pair_records(pair,bx,by,evidence,insert_model);
                return (a,b,self.observed_distance(bx,by));
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
                None,
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

        // Native pair_candidates.sort() compares full candidate tuples on
        // equal pair scores, rather than favoring the last rescue attempted.
        pairs.sort_by(|(a, b), (c, d)| self.pair_joint_phred(a, b, insert_model)
            .cmp(&self.pair_joint_phred(c, d, insert_model))
            .then_with(|| a.native_cmp(c)).then_with(|| b.native_cmp(d)));
        let joint_lks: Vec<f64> = pairs.iter().map(|(x, y)| self.pair_joint_log_lk(x, y, insert_model)).collect();
        let (best_idx, sum_exp) = mq.select_best(&joint_lks);
        let (bx, by) = pairs[best_idx].clone();

        let bx_is_own_best = best1_idx.is_some_and(|i| same_original_locus(&bx, &cands1[i]));
        let by_is_own_best = best2_idx.is_some_and(|j| same_original_locus(&by, &cands2[j]));

        let evidence = self.mapq_for_winning_pair(
            pair, &bx, &by, sum_exp,
            &lk1, &lk2, repeat_mask1, repeat_mask2,
            bx_is_own_best, by_is_own_best, insert_model, false,
        );

        let (a,b) = self.pair_records(pair,&bx,&by,evidence,insert_model);
        (a,b,self.observed_distance(&bx,&by))
    }

    // Native report() updates the library from the chosen pair before SAM
    // rejection. Low-confidence pairs still carry distance observations.
    fn observed_distance(&self, a: &CandidateAlignment, b: &CandidateAlignment) -> Option<f64> {
        if a.contig_id != b.contig_id || a.strand == b.strand { return None; }
        Some(if a.strand == Strand::Forward {
            (b.position + cigar_reference_span(&b.cigar)) as f64 - a.position as f64 - 1.0
        } else { (a.position + cigar_reference_span(&a.cigar)) as f64 - b.position as f64 - 1.0 })
    }

    fn pair_records(&self, pair: &ReadPair, x: &CandidateAlignment, y: &CandidateAlignment,
        e: (u8,f64,u8,f64,u8,u8), model: &InsertSizeDistribution) -> (AlignmentRecord,AlignmentRecord) {
        let mut a=self.single_mapper.to_record(&pair.read1,x,e.0,e.1);
        let mut b=self.single_mapper.to_record(&pair.read2,y,e.2,e.3);
        if a.position != 0 && b.position != 0 {
            let pq = x.coarse_score + y.coarse_score + self.distance_phred(x,y,model);
            a.sam_tags.extend([format!("PQ:i:{pq}"),format!("SM:i:{}",e.4),format!("MQ:i:{}",e.5),format!("XQ:i:{}",(-10.0*y.log_lk/std::f64::consts::LN_10).round() as i32)]);
            b.sam_tags.extend([format!("PQ:i:{pq}"),format!("SM:i:{}",e.5),format!("MQ:i:{}",e.4),format!("XQ:i:{}",(-10.0*x.log_lk/std::f64::consts::LN_10).round() as i32)]);
        }
        (a,b)
    }

    // Pair selection and PQ use candidate anchors before final refinement.
    // The reported fragment length and model updates use observed_distance.
    fn distance_phred(&self, a: &CandidateAlignment, b: &CandidateAlignment, model: &InsertSizeDistribution) -> i32 {
        if a.contig_id != b.contig_id || a.strand == b.strand {
            return (-10.0*self.sv_prior_log_prob/std::f64::consts::LN_10).round() as i32
                + (10.0*(self.single_mapper.genome().total_length as f64).log10()) as i32;
        }
        let distance = if a.strand == Strand::Forward {
            (b.coarse_position + b.out_seq.len()) as f64 - a.coarse_position as f64 - 1.0
        } else { (a.coarse_position + a.out_seq.len()) as f64 - b.coarse_position as f64 - 1.0 };
        if distance.abs() >= 10000.0 {
            return (-10.0*self.sv_prior_log_prob/std::f64::consts::LN_10).round() as i32
                + (10.0*(self.single_mapper.genome().total_length as f64).log10()) as i32;
        }
        let z=(distance-model.mean)/model.sd;
        (18 + (0.5+z*z*5.0/std::f64::consts::LN_10) as i32).min(40)
    }

    /// Combine the single and pair hypotheses in Stampy's Phred arithmetic.
    /// Includes native base-alignment confidence; slow-path candidate search
    /// and repeat-mask inputs still require further fidelity work.
    #[allow(clippy::too_many_arguments)]
    fn mapq_for_winning_pair(
        &self, pair: &ReadPair, bx: &CandidateAlignment, by: &CandidateAlignment,
        pair_sum_exp: f64, lk1: &[f64], lk2: &[f64],
        repeat_mask1: &[f64], repeat_mask2: &[f64],
        bx_is_own_best: bool, by_is_own_best: bool,
        insert_model: &InsertSizeDistribution, fast: bool,
    ) -> (u8, f64, u8, f64, u8, u8) {
        let mq = self.single_mapper.mapq_calc();
        let sum1 = if lk1.is_empty() { 0.0 } else { mq.select_best(lk1).1 };
        let sum2 = if lk2.is_empty() { 0.0 } else { mq.select_best(lk2).1 };
        let (mut wrong1, mut missing1, expected1) = mq.mapping_evidence(sum1, &pair.read1, repeat_mask1);
        let (mut wrong2, mut missing2, expected2) = mq.mapping_evidence(sum2, &pair.read2, repeat_mask2);
        wrong1 = logadd(bx.base_alignment_phred, wrong1);
        wrong2 = logadd(by.base_alignment_phred, wrong2);
        let db = 0.1 * 0.1_f64.ln();
        let mut pair_wrong = if fast { wrong1.max(wrong2) } else {
            (0.5 + (1.0 - 1.0 / pair_sum_exp.max(1.0000000001) + (99.0 * db).exp()).ln() / db).trunc()
        };
        let mut pair_missing = missing1 + missing2;
        let distance = implied_insert_size(bx, by) as f64;
        let z = (distance - insert_model.mean) / insert_model.sd;
        let genome_size = self.single_mapper.genome().total_length as f64;
        if fast {
            let no_repeats = |mask: &[f64], len: usize| {
                let sum: f64 = (0..len).map(|i| {
                    let p = mask.get(i).copied().unwrap_or(0.0);
                    if p <= 0.0 { 1000.0 } else { (-10.0 * p.log10()).round() }
                }).sum();
                sum > 666.0 * len as f64
            };
            if no_repeats(repeat_mask1, pair.read1.sequence.len())
                && no_repeats(repeat_mask2, pair.read2.sequence.len())
                && (bx.coarse_score as f64) < expected1 + 3.0 * expected1.sqrt()
                && (by.coarse_score as f64) < expected2 + 3.0 * expected2.sqrt() {
                pair_missing += 40.0;
            }
        } else if !bx_is_own_best {
            wrong1 = 0.0;
            pair_missing = missing2;
        } else if !by_is_own_best {
            wrong2 = 0.0;
            pair_missing = missing1;
        } else {
            let prior = if bx.contig_id != by.contig_id || distance.abs() > 10000.0 { 0.0 }
                else if z.abs() > 4.0 { (10.0 * (genome_size / 20000.0).log10()).trunc().max(0.0) }
                else { (10.0 * (genome_size / 100.0).log10()).trunc().max(0.0) };
            pair_wrong = logadd(pair_wrong, prior + logadd(wrong1, wrong2));
            pair_missing = logadd(pair_missing, prior + logadd(missing1, missing2));
        }
        let entropy1 = mq.pair_entropy(&bx.cigar, &bx.out_seq, &bx.out_qual, &bx.ref_window, fast);
        let entropy2 = mq.pair_entropy(&by.cigar, &by.out_seq, &by.out_qual, &by.ref_window, fast);
        missing1 = logadd_unclamped(missing1, entropy1);
        missing2 = logadd_unclamped(missing2, entropy2);
        pair_missing = logadd_unclamped(pair_missing, logadd_unclamped(entropy1, entropy2)) + 10.0;
        // Native remove_double_reads rejects the worse fine alignment at
        // the same signed locus; equal scores retain mate 2. SAM uses the
        // minimum pair confidence, so either rejection forces single-end
        // confidence for both records.
        if same_reported_locus(bx, by) {
            if bx.log_lk > by.log_lk { missing2 = -99.0; }
            else { missing1 = -99.0; }
            pair_missing = -99.0;
        }
        let pair_quality = logadd(pair_wrong, pair_missing).min(99.0);
        let single1 = logadd(wrong1, missing1).min(99.0);
        let single2 = logadd(wrong2, missing2).min(99.0);
        let confidence = |q: f64| {
            if q < 0.0 { (0, -1.0) } else { (q as u8, 1.0 - (q * db).exp()) }
        };
        let (q1, p1) = confidence(if pair_quality >= 0.0 { pair_quality } else { single1 });
        let (q2, p2) = confidence(if pair_quality >= 0.0 { pair_quality } else { single2 });
        (q1, p1, q2, p2, single1.max(0.0) as u8, single2.max(0.0) as u8)
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
        let search_radius = (RESCUE_SEARCH_SD * insert_model.sd).ceil() as usize + mate_len;

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

    /// Sum integer Phred terms before conversion so mathematically equal
    /// pair scores remain exact ties, independent of their decomposition.
    fn pair_joint_phred(&self, r1: &CandidateAlignment, r2: &CandidateAlignment, insert_model: &InsertSizeDistribution) -> i32 {
        r1.coarse_score + r2.coarse_score + self.distance_phred(r1, r2, insert_model)
    }

    fn pair_joint_log_lk(&self, r1: &CandidateAlignment, r2: &CandidateAlignment, insert_model: &InsertSizeDistribution) -> f64 {
        self.pair_joint_phred(r1, r2, insert_model) as f64 * (-std::f64::consts::LN_10 / 10.0)
    }

}

fn independent_record(single_mapper: &SingleEndMapper, read: &Read, cands: &[CandidateAlignment], repeat_mask: &[f64]) -> AlignmentRecord {
    if cands.is_empty() {
        return unmapped_record(read);
    }
    let mq: &MapqCalculator = single_mapper.mapq_calc();
    let lk: Vec<f64> = cands.iter().map(|c| c.coarse_score as f64 * (-std::f64::consts::LN_10 / 10.0)).collect();
    let (best_idx, sum_exp) = mq.select_best(&lk);
    let best = &cands[best_idx];
    if !mq.passes_random_match_test(&best.cigar, &best.out_seq, &best.out_qual, &best.ref_window) {
        return unmapped_record(read);
    }
    let (mapq, post) = mq.finalize_mapq(sum_exp, read, &best.cigar, &best.out_seq, &best.out_qual, &best.ref_window, repeat_mask, best.base_alignment_phred);
    single_mapper.to_record(read, best, mapq, post)
}

fn argmax_index(lk: &[f64]) -> Option<usize> {
    lk.iter().enumerate()
        .max_by(|(i, a), (j, b)| a.partial_cmp(b).unwrap().then_with(|| j.cmp(i)))
        .map(|(i, _)| i)
}

fn same_reported_locus(a: &CandidateAlignment, b: &CandidateAlignment) -> bool {
    let endpoint = |c: &CandidateAlignment| if c.strand == Strand::Forward { c.position }
        else { c.position + cigar_reference_span(&c.cigar) };
    a.contig_id == b.contig_id && a.strand == b.strand && endpoint(a) == endpoint(b)
}

// Native surewrongmap compares the original coarse anchor with the final
// signed locus, not two already-refined SAM starts.
fn same_original_locus(refined: &CandidateAlignment, original: &CandidateAlignment) -> bool {
    let anchor = if refined.strand == Strand::Forward { refined.position }
        else { (refined.position + cigar_reference_span(&refined.cigar)).saturating_sub(refined.out_seq.len()) };
    refined.contig_id == original.contig_id && refined.strand == original.strand
        && anchor == original.coarse_position
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
/// review: `distance_phred`, `align_mate_near`, and the fast-path
/// "close" check all used to compute a plain start-to-start gap while the
/// model they compared it against was trained on this TLEN-style value.
fn implied_insert_size(a: &CandidateAlignment, b: &CandidateAlignment) -> i64 {
    let (upstream, downstream) = if a.position <= b.position { (a, b) } else { (b, a) };
    let downstream_end = downstream.position + cigar_reference_span(&downstream.cigar);
    downstream_end as i64 - upstream.position as i64
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
