//! Single-end read mapping.

use crate::types::*;
use crate::index::{GenomeIndex, HashTable};
use crate::align::{SmithWatermanAligner, cigar_reference_span};
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
    /// The reference span the CIGAR actually covers, decoded once so callers
    /// (both `finalize_mapq`'s entropy term and paired-end scoring) don't
    /// have to re-slice the genome.
    pub ref_window: Vec<u8>,
    /// Read sequence/qualities in the *aligned* orientation (i.e.
    /// reverse-complemented already, for a `Strand::Reverse` candidate).
    pub out_seq: Vec<u8>,
    pub out_qual: Vec<u8>,
}

pub struct SingleEndMapper<'a> {
    genome: &'a GenomeIndex,
    hash_table: &'a HashTable,
    sw_aligner: SmithWatermanAligner,
    mapq_calc: MapqCalculator,
    /// Max indel length the alignment window/band must accommodate; taken
    /// from `MapperOptions::max_indel_len`.
    max_indel_len: usize,
    /// Gap-cost inputs kept around (not just baked into `sw_aligner`) so
    /// `align_at_window` can build a second aligner with a much wider band —
    /// normal mapping only ever needs to search `max_indel_len` around a
    /// seed hit, but a rescue/shortlist window's positional uncertainty is
    /// the whole insert-size search radius, which can be an order of
    /// magnitude wider.
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
        let sw = SmithWatermanAligner::new(options.max_indel_len, options.gap_open_phred, options.gap_extend_phred);
        let mq = MapqCalculator::new(options.substitution_rate, hash_table.k, genome.total_length);
        Self {
            genome,
            hash_table,
            sw_aligner: sw,
            mapq_calc: mq,
            max_indel_len: options.max_indel_len,
            gap_open_phred: options.gap_open_phred,
            gap_extend_phred: options.gap_extend_phred,
        }
    }

    pub fn genome(&self) -> &GenomeIndex {
        self.genome
    }

    pub fn mapq_calc(&self) -> &MapqCalculator {
        &self.mapq_calc
    }

    /// Seed, cluster, and align every candidate locus for `read`, returning
    /// each as a `CandidateAlignment` (sorted by score-derived cluster order,
    /// same as the underlying candidate list — *not* re-sorted by
    /// likelihood, so callers that need best-first should sort themselves).
    /// Also returns the per-window repeat mask (see
    /// `mapper::seeds::find_seed_hits`).
    pub fn map_candidates(&self, read: &Read) -> (Vec<CandidateAlignment>, Vec<f64>) {
        // 1. Dense seeding (exact + 1-mismatch neighbours)
        let (hits, repeat_mask) = find_seed_hits(read, self.hash_table, self.genome);

        // 2. Cluster into candidate loci
        let candidates = cluster_hits(&hits, read.sequence.len());

        // Reverse-strand candidates must be aligned against the
        // reverse-complement of the read: the hash table's reverse-strand
        // entries mean the read (as sequenced) matches the genome's minus
        // strand, so ref_seg (always forward-strand genome sequence) has to
        // be compared to rc(read), not read itself. SAM output for a
        // reverse-strand record also uses the RC'd sequence and reversed
        // quality string, so we build them once up front and reuse for
        // every reverse candidate.
        let rc_seq = read.sequence.reverse_complement();
        let mut rc_qual = read.qualities.scores.clone();
        rc_qual.reverse();

        // Reference window padding: allow the alignment to start slightly
        // before the seed anchor (LEFT_PAD) and to extend past the read end by
        // up to a full max-length deletion (RIGHT_PAD), so the DP window can
        // actually contain the largest indel the band is set up to find.
        let left_pad = 10usize;
        let right_pad = self.max_indel_len.max(10);

        let mut alignments = Vec::new();

        for cand in candidates.iter().take(20) {
            let contig = &self.genome.contigs[cand.contig_id];
            let ref_start = cand.position.saturating_sub(left_pad);
            let ref_end = (cand.position + read.sequence.len() + right_pad).min(contig.length);
            if ref_start >= ref_end { continue; }
            // Decode only the small window needed for this candidate — the
            // genome itself stays 2-bit-packed (see types::Contig).
            let ref_seg = contig.slice(ref_start, ref_end);

            let (query_bases, query_qual): (&[u8], &[u8]) = match cand.strand {
                Strand::Forward => (&read.sequence.bases, &read.qualities.scores),
                Strand::Reverse => (&rc_seq.bases, &rc_qual),
            };

            let (_sw_score, cigar, ref_offset) = self.sw_aligner.align(
                query_bases,
                query_qual,
                &ref_seg,
                left_pad as i64,
            );

            let true_start = ref_start + ref_offset;
            // The reference span the CIGAR actually consumes (matches +
            // deletions) can exceed the read length when there's a net
            // deletion, so size the likelihood window off the CIGAR, not
            // off `query_bases.len()`.
            let ref_span = cigar_reference_span(&cigar);
            let ref_window = contig.slice(true_start, (true_start + ref_span).min(contig.length));
            let log_lk = self.mapq_calc.alignment_log_likelihood_bases(query_bases, query_qual, &ref_window, &cigar);

            let (out_seq, out_qual) = match cand.strand {
                Strand::Forward => (read.sequence.bases.clone(), read.qualities.scores.clone()),
                Strand::Reverse => (rc_seq.bases.clone(), rc_qual.clone()),
            };

            alignments.push(CandidateAlignment {
                contig_id: cand.contig_id,
                position: true_start,
                strand: cand.strand,
                cigar,
                log_lk,
                ref_window,
                out_seq,
                out_qual,
            });
        }

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
            &alignments.iter().map(|a| a.log_lk).collect::<Vec<_>>(),
        );
        let best = &alignments[best_idx];

        let (mapq, posterior) = self.mapq_calc.finalize_mapq(
            sum_exp,
            read,
            &best.cigar,
            &best.out_seq,
            &best.out_qual,
            &best.ref_window,
            &repeat_mask,
        );

        Ok(AlignmentRecord {
            read_id: read.id.clone(),
            contig_name: self.genome.contigs[best.contig_id].name.clone(),
            position: best.position + 1, // SAM is 1-based
            mapq,
            cigar: best.cigar.clone(),
            strand: best.strand,
            read_sequence: best.out_seq.clone(),
            read_qualities: best.out_qual.clone(),
            is_paired: false,
            is_proper_pair: false,
            mate_contig: None,
            mate_position: None,
            insert_size: None,
            alignment_score: best.log_lk,
            posterior_prob: posterior,
            read_number: None,
            mate_strand: None,
            mate_unmapped: false,
        })
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
        let (query_bases, query_qual): (&[u8], &[u8]) = (&out_seq, &out_qual);

        // Unlike normal seeded mapping, this window's positional uncertainty
        // is the whole search radius, not just `max_indel_len` —
        // `self.sw_aligner`'s band is far too narrow to reach the true
        // diagonal here (it's centered on a seed hit that, by construction,
        // this call doesn't have). Build a one-off aligner banded to the
        // actual window width instead, and center the diagonal on
        // `expected_position` (where in `ref_seg` this read is predicted to
        // start), not on the search radius itself.
        let window_aligner = SmithWatermanAligner::new(ref_seg.len(), self.gap_open_phred, self.gap_extend_phred);
        let ref_offset_hint = expected_position.saturating_sub(ref_start) as i64;
        let (_score, cigar, ref_offset) = window_aligner.align(
            query_bases,
            query_qual,
            &ref_seg,
            ref_offset_hint,
        );

        let true_start = ref_start + ref_offset;
        let ref_span = cigar_reference_span(&cigar);
        let ref_window = contig.slice(true_start, (true_start + ref_span).min(contig.length));
        let log_lk = self.mapq_calc.alignment_log_likelihood_bases(query_bases, query_qual, &ref_window, &cigar);

        Some(CandidateAlignment {
            contig_id,
            position: true_start,
            strand,
            cigar,
            log_lk,
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
        );
        if mapq == 0 {
            return None;
        }

        Some(AlignmentRecord {
            read_id: read.id.clone(),
            contig_name: self.genome.contigs[cand.contig_id].name.clone(),
            position: cand.position + 1,
            mapq,
            cigar: cand.cigar,
            strand,
            read_sequence: cand.out_seq,
            read_qualities: cand.out_qual,
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
        })
    }

    /// Build the SAM-ready parts of an `AlignmentRecord` from a
    /// `CandidateAlignment` and its final MAPQ/posterior — shared by
    /// `mapper::paired`'s joint-posterior path so it doesn't have to
    /// duplicate this field wiring.
    pub fn to_record(&self, read: &Read, cand: &CandidateAlignment, mapq: u8, posterior: f64) -> AlignmentRecord {
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
        }
    }
}

pub(crate) fn unmapped_record(read: &Read) -> AlignmentRecord {
    AlignmentRecord {
        read_id: read.id.clone(),
        contig_name: "*".to_string(),
        position: 0,
        mapq: 0,
        cigar: vec![CigarOp::SoftClip(read.sequence.len() as u32)],
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
    }
}
