//! Single-end read mapping.

use crate::types::*;
use crate::index::{GenomeIndex, HashTable};
use crate::align::cigar_reference_span;
use crate::align::sensitive::{align_sensitive, SensitiveParameters};
use crate::align::stampy_band::coarse_align;
use crate::mapq::MapqCalculator;
use crate::mapper::{find_seed_hits, cluster_hits};
use anyhow::Result;

/// One candidate mapping location for a single read: its own alignment and
/// single-read log-likelihood, plus everything `MapqCalculator::finalize_mapq`
/// needs. `mapper::paired`'s shortlist/joint-posterior procedure combines
/// these across mates rather than picking one immediately the way
/// `SingleEndMapper::map` does.
#[derive(Clone)]
pub struct CandidateAlignment {
    pub contig_id: usize,
    /// 0-based true alignment start (add 1 for a SAM position).
    pub position: usize,
    pub strand: Strand,
    pub cigar: Vec<CigarOp>,
    pub log_lk: f64,
    pub coarse_score: i32,
    /// Coarse anchor and traceback in sequenced-read orientation. Refinement
    /// can move the SAM start; it must not rewrite candidate tie ordering.
    pub coarse_position: usize,
    pub coarse_read_alignment: Vec<u8>,
    pub coarse_reference_alignment: Vec<u8>,
    pub base_alignment_phred: f64,
    /// The reference span the CIGAR actually covers, decoded once so callers
    /// (both `finalize_mapq`'s entropy term and paired-end scoring) don't
    /// have to re-slice the genome.
    pub ref_window: Vec<u8>,
    /// Read sequence/qualities in the *aligned* orientation (i.e.
    /// reverse-complemented already, for a `Strand::Reverse` candidate).
    pub out_seq: Vec<u8>,
    pub out_qual: Vec<u8>,
}

impl CandidateAlignment {
    /// Python candidate tuples compare score, read alignment, reference
    /// alignment, then signed genomic locus (reverse loci precede forward).
    pub(crate) fn native_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.coarse_score.cmp(&other.coarse_score)
            .then_with(|| self.coarse_read_alignment.cmp(&other.coarse_read_alignment))
            .then_with(|| self.coarse_reference_alignment.cmp(&other.coarse_reference_alignment))
            .then_with(|| locus_key(self.contig_id, self.coarse_position, self.strand)
                .cmp(&locus_key(other.contig_id, other.coarse_position, other.strand)))
    }
}

fn locus_key(contig: usize, position: usize, strand: Strand) -> (bool, usize, usize) {
    if strand == Strand::Reverse { (false, usize::MAX - contig, usize::MAX - position) }
    else { (true, contig, position) }
}

pub struct SingleEndMapper<'a> {
    genome: &'a GenomeIndex,
    hash_table: &'a HashTable,
    mapq_calc: MapqCalculator,
    /// Max indel length the alignment window/band must accommodate; taken
    /// from `MapperOptions::max_indel_len`.
    max_indel_len: usize,
    paired: bool,
    /// Native coarse and sensitive gap-cost parameters.
    gap_open_phred: u32,
    gap_extend_phred: u32,
}

impl<'a> SingleEndMapper<'a> {
    pub fn new(
        genome: &'a GenomeIndex,
        hash_table: &'a HashTable,
        options: &'a MapperOptions,
    ) -> Self {
        // Band the DP wide enough to reach the largest indel the caller asked
        // for (previously hardcoded to 15, ignoring the option).
        let mq = MapqCalculator::new(options.substitution_rate, hash_table.k, genome.total_length);
        Self {
            genome,
            hash_table,
            mapq_calc: mq,
            max_indel_len: options.max_indel_len,
            paired: false,
            gap_open_phred: options.gap_open_phred,
            gap_extend_phred: options.gap_extend_phred,
        }
    }

    pub(crate) fn new_paired(genome: &'a GenomeIndex, hash_table: &'a HashTable, options: &'a MapperOptions) -> Self {
        let mut mapper=Self::new(genome,hash_table,options); mapper.paired=true; mapper
    }

    pub fn genome(&self) -> &GenomeIndex {
        self.genome
    }

    pub fn mapq_calc(&self) -> &MapqCalculator {
        &self.mapq_calc
    }

    /// Seed, cluster, and align every candidate locus for `read`, returning
    /// each as a `CandidateAlignment`, ordered by coarse likelihood.
    /// Also returns the per-window repeat mask (see
    /// `mapper::seeds::find_seed_hits`).
    pub fn map_candidates(&self, read: &Read) -> (Vec<CandidateAlignment>, Vec<f64>) {
        // 1. Dense seeding (exact + 1-mismatch neighbours)
        let (hits, repeat_mask) = find_seed_hits(read, self.hash_table, self.genome);

        // 2. Cluster into candidate loci
        let candidates = cluster_hits(&hits, read.sequence.len());
        for (rank, candidate) in candidates.iter().enumerate() {
            log::trace!(target: "probebwa::fidelity", "cluster\t{}\t{}\t{}\t{}\t{:?}\t{}",
                read.id, rank, self.genome.contigs[candidate.contig_id].name,
                candidate.position, candidate.strand, candidate.score);
        }

        // Rank actual coarse alignments before limiting candidate count. Seed
        // support alone is not an alignment likelihood. Keep 200 candidates
        // within 200 Phred of the best, as in the native candidate budget.
        // Qualities and score scaling depend on the read, not the candidate.
        // Reuse the two orientations across all coarse reference windows.
        let adjusted = self.mapq_calc.adjusted_qualities(&read.qualities.scores);
        let divisor = 1 + (adjusted.iter().map(|&q| q as usize).sum::<usize>() + 256) / 0x3a00;
        let scaled: Vec<u8> = adjusted.iter().map(|&q| (q as usize / divisor) as u8).collect();
        let scale = |cost: u32| cost.div_ceil(divisor as u32) as i32;
        let mut ranked = Vec::new();
        for cand in &candidates {
            let contig = &self.genome.contigs[cand.contig_id];
            // Native coarse traceback always consumes the sequenced read.
            // Reversing the query instead changes gap and equal-window ties.
            let seq = &read.sequence.bases;
            let starts = if cand.strand == Strand::Forward {
                [cand.position.saturating_sub(15), cand.position]
            } else {
                [cand.position, cand.position.saturating_sub(15)]
            };
            let mut best_coarse = None;
            for start in starts {
                let end = (start + seq.len() + 15).min(contig.length);
                let reference = contig.slice(start, end);
                let reference = if cand.strand == Strand::Reverse {
                    DnaSeq::from_ascii(&reference).reverse_complement().bases
                } else { reference };
                if let Some(a) = coarse_align(seq, &scaled, &reference, scale(self.gap_open_phred), scale(self.gap_extend_phred), scale(5)) {
                    let score = a.score * divisor as i32;
                    let position = if cand.strand == Strand::Forward { start + a.offset }
                        else { end.saturating_sub(a.offset + seq.len()) };
                    let span = a.reference.iter().filter(|&&b| b != b'-').count();
                    // Stampy's second window wins equal scores in read orientation.
                    if best_coarse.as_ref().is_none_or(|(best, _, _, _, _)| score <= *best) {
                        best_coarse = Some((score, position, span, a.read, a.reference));
                    }
                }
            }
            let (score, position, span, read_aln, ref_aln) = best_coarse.unwrap_or_else(||
                (0, cand.position, read.sequence.len(), seq.clone(), Vec::new()));
            ranked.push((score, cand.contig_id, position, cand.strand, span, read_aln, ref_aln));
        }
        ranked.sort_by_key(|a| (a.0, locus_key(a.1, a.2, a.3)));
        let mut seen = std::collections::HashSet::new();
        ranked.retain(|item| seen.insert((item.1,item.2,item.3)));
        let best = ranked.first().map(|a| a.0).unwrap_or(0);
        let mut alignments = Vec::new();
        for (score, cid, position, strand, span, coarse_read, coarse_ref) in ranked.into_iter().filter(|a| a.0 <= best + 200).take(200) {
            let padding = 20 + span.saturating_sub(read.sequence.len());
            let padding = if self.max_indel_len > 30 { padding.max(self.max_indel_len) } else { padding };
            if let Some(mut alignment) = self.align_at_window(read, cid, position, padding, strand) {
                let contig = &self.genome.contigs[cid];
                let left = position.saturating_sub(padding);
                let right = (position + read.sequence.len() + padding).min(contig.length);
                let edge = (alignment.position - left).min(right.saturating_sub(alignment.position + cigar_reference_span(&alignment.cigar)));
                if edge < 10 {
                    if let Some(expanded) = self.align_at_window_banded(read, cid, position, padding + read.sequence.len(), strand, (4 * self.max_indel_len.max(30)) as i32) {
                        alignment = expanded;
                    }
                }
                log::trace!(target: "probebwa::fidelity", "candidate\t{}\t{}\t{}\t{:?}\t{:?}\t{}\t{}",
                    read.id, contig.name, alignment.position, strand, alignment.cigar, score, alignment.log_lk);
                alignment.coarse_score = score;
                alignment.coarse_position = position;
                alignment.coarse_read_alignment = coarse_read;
                alignment.coarse_reference_alignment = coarse_ref;
                alignments.push(alignment);
            }
        }
        // Different seed anchors can converge on the same refined mapping.
        // Count that hypothesis once in the posterior denominator.
        alignments.sort_by_key(|a| (a.coarse_score, locus_key(a.contig_id, a.coarse_position, a.strand)));
        let mut unique = std::collections::HashSet::new();
        alignments.retain(|a| unique.insert((a.contig_id, a.position, a.strand)));

        (alignments, repeat_mask)
    }

    pub fn map(&self, read: &Read) -> Result<AlignmentRecord> {
        let (alignments, repeat_mask) = self.map_candidates(read);

        if alignments.is_empty() {
            return Ok(unmapped_record(read));
        }

        // Probabilistic scoring: pick the best candidate by posterior, then
        // compute its final MAPQ (candidate-ambiguity + missing-locus +
        // random-sequence-entropy, combined in Phred space).
        let (best_idx, sum_exp) = self.mapq_calc.select_best(
            &alignments.iter().map(|a| a.coarse_score as f64 * (-std::f64::consts::LN_10 / 10.0)).collect::<Vec<_>>(),
        );
        let best = &alignments[best_idx];

        if !self.mapq_calc.passes_random_match_test(&best.cigar, &best.out_seq, &best.out_qual, &best.ref_window) {
            return Ok(unmapped_record(read));
        }

        let (mapq, posterior) = self.mapq_calc.finalize_mapq(
            sum_exp,
            read,
            &best.cigar,
            &best.out_seq,
            &best.out_qual,
            &best.ref_window,
            &repeat_mask,
            best.base_alignment_phred,
        );

        Ok(self.to_record(read, best, mapq, posterior))
    }

    /// Align `read` directly against a window around `expected_position` on
    /// `contig_id` (on the given `strand`), skipping the normal
    /// seed-and-cluster pipeline entirely — the alignment-only half of what
    /// used to be `rescue`, now also reused by `mapper::paired`'s shortlist
    /// cross-alignment: the mate is aligned against the reference around the
    /// location implied by the library insert size distribution. Returns
    /// `None` only if the window falls outside the contig; unlike `rescue`,
    /// this does *not* gate on MAPQ — callers that need that gate (plain
    /// mate rescue) do it themselves via `MapqCalculator::finalize_mapq`.
    pub fn align_at_window(&self, read: &Read, contig_id: usize, expected_position: usize, search_radius: usize, strand: Strand) -> Option<CandidateAlignment> {
        self.align_at_window_banded(read, contig_id, expected_position, search_radius, strand, (2 * self.max_indel_len.max(30)) as i32)
    }

    // Native realign_single doubles the band as well as expanding the window
    // on an edge retry. BAQ must use the same band as the selected alignment.
    fn align_at_window_banded(&self, read: &Read, contig_id: usize, expected_position: usize, search_radius: usize, strand: Strand, band: i32) -> Option<CandidateAlignment> {
        let contig = &self.genome.contigs[contig_id];
        let ref_start = expected_position.saturating_sub(search_radius);
        let ref_end = (expected_position + read.sequence.len() + search_radius).min(contig.length);
        if ref_start >= ref_end { return None; }
        let ref_seg = contig.slice(ref_start, ref_end);

        // Only build the reverse-complement when this call is actually for
        // the reverse strand -- shortlist cross-alignment tries both
        // orientations across its candidates, so roughly half of all calls
        // used to pay for a `reverse_complement()` + quality-reverse alloc
        // whose result was then discarded unused.
        let (out_seq, out_qual): (Vec<u8>, Vec<u8>) = match strand {
            Strand::Forward => (read.sequence.bases.clone(), read.qualities.scores.clone()),
            Strand::Reverse => {
                let rc_seq = read.sequence.reverse_complement();
                let mut rc_qual = read.qualities.scores.clone();
                rc_qual.reverse();
                (rc_seq.bases, rc_qual)
            }
        };
        // Stampy aligns the sequenced read against the reference in its
        // strand orientation, then reverses the alignment for SAM reporting.
        let oriented_ref = if strand == Strand::Reverse {
            DnaSeq::from_ascii(&ref_seg).reverse_complement().bases
        } else { ref_seg.clone() };
        let (query_bases, query_qual): (&[u8], &[u8]) = (&read.sequence.bases, &read.qualities.scores);

        let adjusted = if self.paired { query_qual.to_vec() } else { self.mapq_calc.adjusted_qualities(query_qual) };
        let extend = self.gap_extend_phred.max(1) as f64;
        let norm = (0.5 + 10.0 * (1.0 - (-std::f64::consts::LN_10 * extend / 10.0).exp()
            / -std::f64::consts::LN_10).ln()) as i32;
        let parameters = SensitiveParameters {
            insertion_prior: 5.0,
            gap_open: self.gap_open_phred as f64 + norm as f64,
            gap_extend: extend,
            start_mean: 0,
            start_sd: 10000,
        };
        let result = align_sensitive(query_bases, &adjusted, &oriented_ref, parameters, band, 999999)?;
        let has_indel = result.primary.read.contains(&b'-') || result.primary.reference.contains(&b'-');
        let base_alignment_phred = if has_indel {
            let evidence = crate::align::confidence::sensitive_confidence(query_bases, &adjusted, &oriented_ref, parameters, band, &result)?;
            let probability = evidence[0].probability;
            if probability >= 1.0 { 999.0 } else { (0.5 - 10.0 * (1.0-probability).log10()).trunc().max(0.0) }
        } else { 999.0 };
        let offset = result.primary.offset;
        let mut path = if strand == Strand::Forward { result.alternate } else { result.primary };
        let coarse_read_alignment = path.read.clone();
        let coarse_reference_alignment = path.reference.clone();
        if strand == Strand::Reverse {
            let reverse = |s: &mut Vec<u8>| {
                s.reverse();
                for b in s { *b = match *b { b'A'=>b'T',b'T'=>b'A',b'C'=>b'G',b'G'=>b'C',x=>x }; }
            };
            reverse(&mut path.read); reverse(&mut path.reference);
        }
        let cigar = crate::align::cigar_from_alignment(&path.read, &path.reference)?;
        let ref_span = cigar_reference_span(&cigar);
        let true_start = if strand == Strand::Forward { ref_start + offset } else { ref_end.checked_sub(offset + ref_span)? };
        if ref_span == 0 || true_start + ref_span > contig.length { return None; }
        let ref_window = contig.slice(true_start, true_start + ref_span);
        let log_lk = result.score as f64 * (-std::f64::consts::LN_10 / 10.0);

        Some(CandidateAlignment {
            contig_id,
            position: true_start,
            strand,
            cigar,
            log_lk,
            coarse_score: result.score,
            coarse_position: if strand == Strand::Forward { true_start }
                else { (true_start + ref_span).saturating_sub(read.sequence.len()) },
            coarse_read_alignment,
            coarse_reference_alignment,
            base_alignment_phred,
            ref_window,
            out_seq,
            out_qual,
        })
    }

    /// Mate rescue: `align_at_window` plus the same MAPQ gate normal mapping
    /// uses. Used when a read's mate mapped confidently but the read itself
    /// produced no seed hits of its own — the mapped mate's position plus
    /// the current insert-size estimate gives a strong-enough prior for
    /// where this read should be to go straight to alignment rather than
    /// relying on this read's own (by definition unsuccessful) seeding.
    /// Returns `None` if the window falls outside the contig or the
    /// resulting alignment doesn't clear `finalize_mapq`'s entropy check (so
    /// a rescue "window" that's actually unrelated sequence still gets
    /// rejected rather than reported).
    pub fn rescue(&self, read: &Read, contig_id: usize, expected_position: usize, search_radius: usize, strand: Strand) -> Option<AlignmentRecord> {
        let cand = self.align_at_window(read, contig_id, expected_position, search_radius, strand)?;

        // A rescue only ever considers this one window (no competing
        // candidates), so the softmax denominator is trivially 1 — the
        // `finalize_mapq` gate that actually matters here is `logentropy`
        // (is this alignment meaningfully better than a random sequence of
        // the genome's size?), not candidate-vs-candidate ambiguity.
        // No `repeat_mask`: rescue skips seeding entirely (that's the whole
        // point — this read produced no seed hits of its own), so there's no
        // per-window repetitiveness signal to pass; an empty slice makes
        // `missing_locus_probability` fall back to its quality-only formula.
        let (mapq, posterior) = self.mapq_calc.finalize_mapq(
            1.0,
            read,
            &cand.cigar,
            &cand.out_seq,
            &cand.out_qual,
            &cand.ref_window,
            &[],
            cand.base_alignment_phred,
        );
        if !self.mapq_calc.passes_random_match_test(&cand.cigar, &cand.out_seq, &cand.out_qual, &cand.ref_window) {
            return None;
        }

        Some(self.to_record(read, &cand, mapq, posterior))
    }

    /// Build the SAM-ready parts of an `AlignmentRecord` from a
    /// `CandidateAlignment` and its final MAPQ/posterior — shared by
    /// `mapper::paired`'s joint-posterior path so it doesn't have to
    /// duplicate this field wiring.
    fn alignment_tags(&self, cand: &CandidateAlignment) -> Vec<String> {
        let (mut q, mut r, mut nm) = (0, 0, 0);
        for op in &cand.cigar {
            match *op {
                CigarOp::Match(n) => for _ in 0..n {
                    let (a,b)=(cand.out_seq[q],cand.ref_window[r]);
                    if a != b && a != b'N' && b != b'N' { nm += 1; }
                    q+=1; r+=1;
                },
                CigarOp::Insertion(n) => { q+=n as usize; nm+=n; },
                CigarOp::Deletion(n) => { r+=n as usize; nm+=n; },
                CigarOp::SoftClip(n) => q+=n as usize,
                CigarOp::HardClip(_) => (),
            }
        }
        vec![format!("UQ:i:{}", (cand.log_lk * (-10.0 / std::f64::consts::LN_10)).round() as i32), format!("NM:i:{nm}")]
    }

    pub fn to_record(&self, read: &Read, cand: &CandidateAlignment, mapq: u8, posterior: f64) -> AlignmentRecord {
        if posterior < 0.0 {
            return unmapped_record(read);
        }
        AlignmentRecord {
            read_id: read.id.clone(),
            contig_name: self.genome.contigs[cand.contig_id].name.clone(),
            position: cand.position + 1,
            mapq,
            cigar: cand.cigar.clone(),
            strand: cand.strand,
            read_sequence: cand.out_seq.clone(),
            read_qualities: cand.out_qual.clone(),
            is_paired: false,
            is_proper_pair: false,
            mate_contig: None,
            mate_position: None,
            insert_size: None,
            alignment_score: cand.log_lk,
            posterior_prob: posterior,
            read_number: None,
            mate_strand: None,
            mate_unmapped: false,
            sam_tags: self.alignment_tags(cand),
        }
    }
}

pub(crate) fn unmapped_record(read: &Read) -> AlignmentRecord {
    AlignmentRecord {
        read_id: read.id.clone(),
        contig_name: "*".to_string(),
        position: 0,
        mapq: 0,
        cigar: Vec::new(),
        strand: Strand::Forward,
        read_sequence: read.sequence.bases.clone(),
        read_qualities: read.qualities.scores.clone(),
        is_paired: false,
        is_proper_pair: false,
        mate_contig: None,
        mate_position: None,
        insert_size: None,
        alignment_score: 0.0,
        posterior_prob: 0.0,
        read_number: None,
        mate_strand: None,
        mate_unmapped: false,
        sam_tags: Vec::new(),
    }
}

#[cfg(test)]
mod selection_tests {
    use super::*;
    use std::io::Read as _;

    fn native_candidate(fields: &[&str]) -> CandidateAlignment {
        let locus: i64 = fields[3].parse().unwrap();
        CandidateAlignment {
            contig_id: 0, position: 0,
            strand: if locus < 0 { Strand::Reverse } else { Strand::Forward },
            cigar: vec![], log_lk: 0.0,
            coarse_score: fields[0].parse().unwrap(),
            coarse_position: locus.unsigned_abs() as usize,
            coarse_read_alignment: fields[1].as_bytes().to_vec(),
            coarse_reference_alignment: fields[2].as_bytes().to_vec(),
            base_alignment_phred: 999.0,
            ref_window: vec![], out_seq: vec![], out_qual: vec![],
        }
    }

    #[test]
    fn native_pair_tuple_tie_order() {
        let mut source = String::new();
        flate2::read::GzDecoder::new(&include_bytes!("../../tests/fixtures/stampy_pair_ties.tsv.gz")[..])
            .read_to_string(&mut source).unwrap();
        let mut cases = std::collections::BTreeMap::<usize, Vec<_>>::new();
        for line in source.lines().filter(|line| !line.starts_with('#')) {
            let fields: Vec<_> = line.split('\t').collect();
            cases.entry(fields[0].parse().unwrap()).or_default()
                .push((native_candidate(&fields[2..6]), native_candidate(&fields[6..10])));
        }
        assert_eq!(cases.len(), 16);
        assert_eq!(cases.values().map(Vec::len).sum::<usize>(), 93);
        for (case, expected) in cases {
            let mut actual = expected.clone();
            actual.reverse();
            actual.sort_by(|(a, b), (c, d)| a.native_cmp(c).then_with(|| b.native_cmp(d)));
            for ((a, b), (c, d)) in actual.iter().zip(&expected) {
                assert_eq!(a.native_cmp(c), std::cmp::Ordering::Equal, "case {case}, mate 1");
                assert_eq!(b.native_cmp(d), std::cmp::Ordering::Equal, "case {case}, mate 2");
            }
        }
    }
}
