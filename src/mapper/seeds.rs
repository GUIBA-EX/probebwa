//! Dense seeding — the sensitivity engine of this mapper.
//!
//! For every overlapping 15-mer in the read, we always query the exact
//! 15-mer. Single-mismatch neighbours (15 positions × 3 alternative bases)
//! are additionally queried, but only at a length-dependent subset of
//! positions, rather than at every position.
//!
//! The fraction of positions that get neighbour-probed is length-dependent:
//! for reads longer than 34bp, single-mismatch neighbours are considered for
//! a reduced fraction of initial 15-mers -- half of them for reads up to
//! 49bp, a third for reads of 50bp and above. This is designed around the
//! hash table only indexing every 5th genomic position (see
//! `index::hashtable::DEFAULT_STRIDE`) -- a read only needs one
//! neighbour-probed window per 5-position block to have a chance of hitting
//! an indexed locus. The one-third case uses a block of 5 consecutive
//! positions probed, then 10 skipped, repeating; the one-half case is the
//! analogous block-of-5-probed/block-of-5-skipped pattern (5 out of every
//! 10), the natural generalisation that yields a 1/2 fraction.
//!
//! Probing every mismatch variant at every seed position is a strict
//! superset that's more sensitive but multiplies the number of hash probes
//! by ~46x; restricting those probes to the appropriate fraction — while
//! still doing exact matching everywhere — trades a little sensitivity for a
//! large speedup.
//!
//! Performance notes: mismatch-variant probing is the hot inner loop here
//! (up to k×3 hash lookups per sampled position), so it's written to do zero
//! heap allocation — `probe_1mm_neighbours` mutates one stack buffer in place
//! rather than materializing a fresh `Vec<u8>` per variant. Variant probing
//! is also skipped entirely at positions whose exact-match hit count already
//! looks repetitive: those loci are already ambiguous from the exact match
//! alone, so spending ~45 more probes there buys little (this repeat-skip
//! isn't itself part of the published fraction rule, but is a natural
//! extension of it).

use crate::types::*;
use crate::index::{GenomeIndex, HashTable};

/// Max k-mer length supported by the stack-buffer neighbour probe below.
const MAX_K: usize = 32;

/// If a position's exact-match lookup already returns at least this many
/// hits, treat it as a repetitive locus and skip mismatch-variant probing
/// there.
const REPEAT_SKIP_THRESHOLD: usize = 64;

/// Number of non-overlapping 4-mers making up a similarity fingerprint
/// (three 4-nucleotide words).
const FINGERPRINT_WORD_LEN: usize = 4;
const FINGERPRINT_WORDS: usize = 3;
const FINGERPRINT_SPAN: usize = FINGERPRINT_WORD_LEN * FINGERPRINT_WORDS; // 12
/// How many of the 12 flanking bases are drawn from after vs. before the
/// window when both sides have room (2 words after, 1 word before) — the
/// exact placement isn't load-bearing (the words just need to fall close to
/// but not overlapping the 15-mer, within the read), so this split is this
/// module's own reasonable choice.
const FINGERPRINT_AFTER_WANT: usize = 8;
const FINGERPRINT_BEFORE_WANT: usize = 4;

/// Similarity-filter acceptance threshold on `fingerprint_distance`, and how
/// much lower it is for a 1-mismatch ("neighbour") 15-mer: only locations
/// whose fingerprint match value doesn't exceed this threshold are
/// considered, reduced by 2 for one-mismatch 15-mers. This module picks a
/// permissive default based on the statistic's own property that it
/// increases by at most 2 with every single-nucleotide change — loose
/// enough to tolerate several real substitutions in the flanking sequence
/// without discarding a true candidate, while still rejecting flanking
/// sequence that looks nothing alike.
const FINGERPRINT_BASE_THRESHOLD: u32 = 10;
const FINGERPRINT_NEIGHBOUR_DISCOUNT: u32 = 2;

/// Find all seed hits for a read: exact matches at every position, plus
/// single-mismatch neighbours at a length-dependent, block-sampled subset of
/// non-repetitive positions (see module docs).
///
/// Also returns a per-window repeat-mask *probability* (not just a bool):
/// `repeat_prob[read_pos]` is this project's estimate of "the true reference
/// 15-mer here is repetitive", used by `mapq::bayesian`'s missing-locus model
/// to avoid treating a repetitive window as a reliable candidate (see
/// `mapq::bayesian::MapqCalculator::missing_locus_probability`): 1.0 when
/// the window's own exact 15-mer is repetitive; otherwise, when a probed
/// 1-mismatch neighbour is repetitive, the probability that *that*
/// neighbour is the true sequence — i.e. that this exact base really was a
/// sequencing error/variant — derived from the read's own quality at the
/// mutated position (`10^(-q/10)`), taking the maximum such probability
/// across all repetitive neighbours found (a "worst case wins" rule); 0.0
/// otherwise. Only checked at positions that actually get neighbour-probed
/// (`neighbour_probed_position`) — the same length-dependent sampling
/// seeding itself uses, since checking every neighbour everywhere would
/// undo the whole point of that sampling.
pub fn find_seed_hits(read: &Read, hash_table: &HashTable, genome: &GenomeIndex) -> (Vec<SeedHit>, Vec<f64>) {
    let k = hash_table.k;
    let mut hits = Vec::new();

    let read_len = read.sequence.bases.len();
    if read_len < k || k > MAX_K {
        return (hits, Vec::new());
    }
    let mut repeat_prob = vec![0.0f64; read_len - k + 1];

    // Reverse-strand hits are filtered in the same orientation the aligner
    // ultimately uses (rc(read) against the forward-strand genome — see
    // `mapper::single`), so the RC is computed once up front rather than
    // per hit.
    let rc_read = read.sequence.reverse_complement();

    for (read_pos, window) in read.sequence.bases.windows(k).enumerate() {
        // A window containing N (or another non-ACGT ambiguity code) is
        // skipped rather than hashed: `HashTable::hash_kmer` treats N the
        // same as 'A', so hashing it here would spuriously match whatever
        // fabricated 'A'-run entries (or genuine poly-A loci) share that
        // hash, and the reference side already skips N-containing windows
        // when building the table (see `HashTable::is_acgt_only`).
        if !HashTable::is_acgt_only(window) { continue; }

        // 1. Exact match — tried at every position, unconditionally.
        let exact_hits = hash_table.lookup(window);
        for &packed in exact_hits {
            let (cid, rpos, strand) = HashTable::unpack_position(packed);
            if passes_similarity_filter(read, &rc_read, read_pos, k, genome, cid, rpos, strand, false) {
                hits.push(SeedHit { read_pos, ref_contig: cid, ref_pos: rpos, k, hit_strand: strand });
            }
        }

        let looks_repetitive = exact_hits.len() >= REPEAT_SKIP_THRESHOLD;
        if looks_repetitive {
            repeat_prob[read_pos] = 1.0;
        }

        // 2. Single-mismatch neighbours (up to k × 3), only at the
        // length-appropriate fraction of non-repetitive positions.
        if neighbour_probed_position(read_pos, read_len) && !looks_repetitive {
            let qual = read.qualities.scores.get(read_pos).copied().unwrap_or(30);
            let neighbour_repeat_prob = 10f64.powf(-(qual as f64) / 10.0);
            probe_1mm_neighbours(window, |neighbour| {
                let neighbour_hits = hash_table.lookup(neighbour);
                if neighbour_hits.len() >= REPEAT_SKIP_THRESHOLD {
                    repeat_prob[read_pos] = repeat_prob[read_pos].max(neighbour_repeat_prob);
                }
                for &packed in neighbour_hits {
                    let (cid, rpos, strand) = HashTable::unpack_position(packed);
                    if passes_similarity_filter(read, &rc_read, read_pos, k, genome, cid, rpos, strand, true) {
                        hits.push(SeedHit { read_pos, ref_contig: cid, ref_pos: rpos, k, hit_strand: strand });
                    }
                }
            });
        }
    }

    (hits, repeat_prob)
}

/// Similarity filter: a cheap base-composition comparison of the sequence
/// flanking the seed window in the read vs. the candidate genomic position,
/// applied before a hit is even added to the candidate-clustering pool —
/// screening out obviously-unrelated near-repeat hash collisions well before
/// the much more expensive full alignment stage sees them. See
/// `fingerprint_distance` for the actual statistic and
/// `FINGERPRINT_BASE_THRESHOLD` for why the acceptance threshold's value is
/// this module's own choice.
#[allow(clippy::too_many_arguments)]
fn passes_similarity_filter(
    read: &Read,
    rc_read: &DnaSeq,
    read_pos: usize,
    k: usize,
    genome: &GenomeIndex,
    contig_id: usize,
    ref_pos: usize,
    strand: Strand,
    is_neighbour: bool,
) -> bool {
    let contig = &genome.contigs[contig_id];

    // Orient both sides identically: for a reverse-strand hit, `ref_pos` is
    // the forward-genome coordinate where rc(read)'s window matches (see
    // `index::hashtable` build docs and `mapper::candidates`'s diagonal
    // derivation), so using rc(read) here — rather than the read as
    // sequenced — makes the "window start" and "before/after" directions
    // agree with the forward-strand genome coordinates used below.
    let (read_local, window_start) = match strand {
        Strand::Forward => (&read.sequence.bases, read_pos),
        Strand::Reverse => (&rc_read.bases, read.sequence.len() - read_pos - k),
    };

    let (read_after, read_before) = clip_flanks(read_local.len(), window_start, k, FINGERPRINT_AFTER_WANT, FINGERPRINT_BEFORE_WANT);

    // Fetch only the small genomic region actually needed (the genome stays
    // 2-bit-packed otherwise — see `types::Contig`), sized generously enough
    // to always contain the full wanted flank even before clipping.
    let region_start = ref_pos.saturating_sub(FINGERPRINT_BEFORE_WANT);
    let region_end = (ref_pos + k + FINGERPRINT_AFTER_WANT).min(contig.length);
    if region_start >= region_end {
        return true; // no genomic context available — don't filter
    }
    let genome_local = contig.slice(region_start, region_end);
    let genome_window_start = ref_pos - region_start;
    let (genome_after, genome_before) = clip_flanks(genome_local.len(), genome_window_start, k, FINGERPRINT_AFTER_WANT, FINGERPRINT_BEFORE_WANT);

    // Use whichever side (read or genome) has less room, so both
    // fingerprints are computed over exactly the same span — a sequence edge
    // on either side just shrinks the compared region rather than
    // misaligning the two sides.
    let after_len = read_after.len().min(genome_after.len());
    let before_len = read_before.len().min(genome_before.len());
    if after_len + before_len < FINGERPRINT_SPAN {
        return true; // not enough flanking sequence on either side to judge
    }

    let mut read_fp = compute_fingerprint(&read_local[read_after.start..read_after.start + after_len]);
    read_fp.merge(compute_fingerprint(&read_local[read_before.end - before_len..read_before.end]));

    let mut genome_fp = compute_fingerprint(&genome_local[genome_after.start..genome_after.start + after_len]);
    genome_fp.merge(compute_fingerprint(&genome_local[genome_before.end - before_len..genome_before.end]));

    let distance = read_fp.distance(&genome_fp);
    let threshold = FINGERPRINT_BASE_THRESHOLD.saturating_sub(if is_neighbour { FINGERPRINT_NEIGHBOUR_DISCOUNT } else { 0 });
    distance <= threshold
}

/// Up to `after_want` bases immediately following `[window_start,
/// window_start+k)`, and up to `before_want` bases immediately preceding it,
/// both clipped to `[0, seq_len)`.
fn clip_flanks(seq_len: usize, window_start: usize, k: usize, after_want: usize, before_want: usize) -> (std::ops::Range<usize>, std::ops::Range<usize>) {
    let after_start = (window_start + k).min(seq_len);
    let after_end = (after_start + after_want).min(seq_len);

    let before_end = window_start.min(seq_len);
    let before_start = before_end.saturating_sub(before_want);

    (after_start..after_end, before_start..before_end)
}

/// Base-composition counts over a flanking region, used by the similarity
/// filter. Counts all 4 bases directly (rather than inferring T by
/// subtraction from a fixed total) — equivalent when every position is a
/// called ACGT base, and more robust if an N slips into a flanking region.
#[derive(Clone, Copy, Default)]
struct Fingerprint { a: u32, c: u32, g: u32, t: u32 }

fn compute_fingerprint(seq: &[u8]) -> Fingerprint {
    let mut fp = Fingerprint::default();
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'A' => fp.a += 1,
            b'C' => fp.c += 1,
            b'G' => fp.g += 1,
            b'T' => fp.t += 1,
            _ => {}
        }
    }
    fp
}

impl Fingerprint {
    fn merge(&mut self, other: Fingerprint) {
        self.a += other.a;
        self.c += other.c;
        self.g += other.g;
        self.t += other.t;
    }

    /// Sum of absolute per-base-count differences (the similarity-filter
    /// statistic).
    fn distance(&self, other: &Fingerprint) -> u32 {
        (self.a as i64 - other.a as i64).unsigned_abs() as u32
            + (self.c as i64 - other.c as i64).unsigned_abs() as u32
            + (self.g as i64 - other.g as i64).unsigned_abs() as u32
            + (self.t as i64 - other.t as i64).unsigned_abs() as u32
    }
}

/// Whether `read_pos` is one of the 15-mer start positions that gets
/// single-mismatch neighbour probing, per the length-dependent fraction and
/// block pattern described in the module docs.
fn neighbour_probed_position(read_pos: usize, read_len: usize) -> bool {
    if read_len <= 34 {
        true
    } else if read_len <= 49 {
        // One half: alternating blocks of 5 probed / 5 skipped.
        (read_pos / 5).is_multiple_of(2)
    } else {
        // One third: a block of 5 probed, then 10 skipped.
        (read_pos % 15) < 5
    }
}

/// Visit exactly (k × 3) single-mismatch neighbours of `kmer`, without
/// allocating: a fixed-size stack buffer is mutated one base at a time and
/// handed to `visit`, then restored before moving to the next position.
fn probe_1mm_neighbours<F: FnMut(&[u8])>(kmer: &[u8], mut visit: F) {
    const BASES: [u8; 4] = *b"ACGT";
    let k = kmer.len();
    debug_assert!(k <= MAX_K);

    let mut buf = [0u8; MAX_K];
    buf[..k].copy_from_slice(kmer);

    for i in 0..k {
        let original = buf[i];
        for &base in &BASES {
            if base != original {
                buf[i] = base;
                visit(&buf[..k]);
            }
        }
        buf[i] = original;
    }
}
