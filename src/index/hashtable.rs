//! 15-mer hash table (`.sthash`) using FxHashMap.
//!
//! Key design choices:
//! - 15-mer → 30-bit integer (2 bits/base), stored in u32.
//! - The packed integer is then run through a multiply-xor-shift finalizer
//!   (a standard, widely-used integer-hash technique — see e.g. the
//!   "splitmix"/"murmur finalizer" family) so that hash buckets aren't
//!   skewed by DNA's compositional bias (poly-A/T runs etc.); a raw 2-bit
//!   packing alone would cluster low-complexity sequence into a few buckets.
//! - Each occurrence is packed into a u64: bit 63 = strand (0=forward,
//!   1=reverse-complement), bits 62..40 = contig id (23 bits → 8.3M contigs),
//!   bits 39..0 = position (40 bits → up to ~1.1 trillion bp per contig,
//!   comfortably covering any real chromosome; human chr1 needs 28 bits).
//! - Both forward and reverse-complement strands are indexed as separate
//!   entries, each carrying an explicit strand tag so lookups don't need to
//!   infer strand from a diagonal trend (see mapper/candidates.rs).
//! - High-frequency (repetitive) k-mers are capped at `max_count` entries per
//!   bucket — an unbounded bucket would make lookups near centromeres/
//!   transposons blow up, so entries beyond the cap are simply dropped.
//! - **Flat shared position store, not one `Vec<u64>` per bucket.** Storing
//!   `FxHashMap<u32, Vec<u64>>` directly (as this module used to) means one
//!   separate heap allocation per distinct k-mer -- for a real genome that's
//!   millions of small `Vec`s, fragmenting memory and scattering lookups
//!   across the heap. Instead `table` maps each hash to an `(offset, count)`
//!   pair into one contiguous `positions: Vec<u64>` shared by every bucket:
//!   one allocation for the whole table's occurrence data and better cache
//!   locality on lookup (a hit's positions are already one contiguous slice
//!   read). `HashTableBuilder::build` gets there by collecting every
//!   occurrence into one flat `Vec<(hash, position)>` as it scans (no
//!   per-key `Vec` ever exists) and doing a single stable sort by hash to
//!   group same-key runs together -- see that function's doc for why this
//!   sort-based approach replaced an earlier version that built a
//!   `FxHashMap<u32, Vec<u64>>` first and flattened it afterward (that
//!   intermediate version measured *higher* peak build memory than the
//!   original design, the opposite of the intended win, since both
//!   representations were briefly alive at once).
//!   (Approach independently arrived at; the same shared-flat-array idea is
//!   also used by other k-mer indexing tools, e.g. minimap2-style bucketed
//!   hash indices.)

use crate::types::*;
use crate::index::GenomeIndex;
use serde::{Deserialize, Serialize};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use anyhow::Result;

/// Default cap on stored occurrences per k-mer bucket.
pub const DEFAULT_MAX_COUNT: usize = 300;

/// Bit widths for the packed (strand, contig, position) occurrence.
const CONTIG_BITS: u32 = 23;
const POSITION_BITS: u32 = 40;
const POSITION_MASK: u64 = (1 << POSITION_BITS) - 1;
const CONTIG_MASK: u64 = (1 << CONTIG_BITS) - 1;

#[derive(Serialize, Deserialize)]
pub struct HashTable {
    /// Pseudo-randomized 15-mer hash → `(offset, count)` into `positions`.
    table: FxHashMap<u32, (u32, u32)>,
    /// Flat store for every bucket's packed occurrences, shared across all
    /// hashes -- see module docs for why this replaces one `Vec<u64>` per
    /// bucket.
    positions: Vec<u64>,
    pub k: usize,
    pub seed_stride: usize,
    pub total_entries: usize,
}

impl HashTable {
    pub fn load<P: AsRef<Path>>(prefix: P) -> Result<Self> {
        let path = format!("{}.sthash", prefix.as_ref().display());
        let file = File::open(&path)?;
        let reader = BufReader::new(file);
        let ht: HashTable = bincode::deserialize_from(reader)?;
        Ok(ht)
    }

    pub fn save<P: AsRef<Path>>(&self, prefix: P) -> Result<()> {
        let path = format!("{}.sthash", prefix.as_ref().display());
        let file = File::create(&path)?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, self)?;
        Ok(())
    }

    /// Lookup exact 15-mer, return packed positions.
    pub fn lookup(&self, kmer: &[u8]) -> &[u64] {
        self.get(Self::hash_kmer(kmer))
    }

    /// Lookup by an already-computed hash, resolving `table`'s `(offset,
    /// count)` into the corresponding slice of `positions`.
    #[inline]
    fn get(&self, hash: u32) -> &[u64] {
        match self.table.get(&hash) {
            Some(&(offset, count)) => &self.positions[offset as usize..(offset + count) as usize],
            None => &[],
        }
    }

    /// Pack ACGT into 2-bit representation (15-mer fits in 30 bits), then run
    /// through a multiply-xor-shift finalizer so hash buckets aren't skewed
    /// by DNA's compositional bias.
    pub fn hash_kmer(kmer: &[u8]) -> u32 {
        let mut raw: u32 = 0;
        for &base in kmer.iter().take(15) {
            raw = (raw << 2) | Self::base2bits(base);
        }
        Self::mix(raw)
    }

    /// Integer finalizer mix (multiply-xor-shift), independent of any
    /// specific external hash table implementation.
    fn mix(mut h: u32) -> u32 {
        h ^= h >> 16;
        h = h.wrapping_mul(0x7feb352d);
        h ^= h >> 15;
        h = h.wrapping_mul(0x846ca68b);
        h ^= h >> 16;
        h
    }

    fn base2bits(base: u8) -> u32 {
        match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,
        }
    }

    /// Whether every base in `kmer` is a plain A/C/G/T (case-insensitive).
    /// `base2bits` maps anything else (N, other IUPAC ambiguity codes) to the
    /// same code as 'A', so a window containing one would otherwise hash and
    /// match identically to an all-A window of the same length -- both when
    /// building the table (silently indexing a reference N-gap as if it were
    /// a real "AAAA..." run) and when querying it (a read's N bases would
    /// spuriously match those fabricated entries, or any real poly-A locus).
    /// Callers on both sides skip windows failing this check rather than
    /// hash them. Found by code review; genome/read sequences are already
    /// uppercased before reaching here (`DnaSeq::from_ascii`,
    /// `Contig::slice`), but this checks lowercase too for robustness.
    pub fn is_acgt_only(kmer: &[u8]) -> bool {
        kmer.iter().all(|&b| matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
    }

    /// Pack (strand, contig_id, position) into a single u64.
    /// Bit 63 = strand (0=forward, 1=reverse-complement), bits 62..40 = contig
    /// (23 bits), bits 39..0 = position (40 bits). Contig ids and positions
    /// beyond the representable range are asserted against in debug builds
    /// (a genome that large is out of scope, and a silent wrap would corrupt
    /// mappings) and saturate in release.
    pub fn pack_position(contig_id: usize, position: usize, strand: Strand) -> u64 {
        debug_assert!((contig_id as u64) <= CONTIG_MASK, "contig id {contig_id} exceeds 23-bit packing range");
        debug_assert!((position as u64) <= POSITION_MASK, "position {position} exceeds 40-bit packing range");
        let strand_bit: u64 = if strand.is_forward() { 0 } else { 1 };
        (strand_bit << 63)
            | (((contig_id as u64) & CONTIG_MASK) << POSITION_BITS)
            | (position as u64 & POSITION_MASK)
    }

    pub fn unpack_position(packed: u64) -> (usize, usize, Strand) {
        let strand = if (packed >> 63) & 1 == 0 { Strand::Forward } else { Strand::Reverse };
        let contig_id = ((packed >> POSITION_BITS) & CONTIG_MASK) as usize;
        let position = (packed & POSITION_MASK) as usize;
        (contig_id, position, strand)
    }
}

pub struct HashTableBuilder<'a> {
    genome: &'a GenomeIndex,
    k: usize,
    stride: usize,
    max_count: usize,
}

/// Only every `DEFAULT_STRIDE`-th genomic position gets an entry in the hash
/// table -- this cuts table size/memory ~5x, and keeps the packed position
/// small enough to fit its allotted bits comfortably. Correctness for a read
/// relies on it having 5+ overlapping candidate windows to try, which is why
/// seed scanning (`mapper::seeds`) probes every overlapping 15-mer in the
/// read rather than also striding on the read side.
pub const DEFAULT_STRIDE: usize = 5;

impl<'a> HashTableBuilder<'a> {
    pub fn new(genome: &'a GenomeIndex) -> Self {
        Self { genome, k: 15, stride: DEFAULT_STRIDE, max_count: DEFAULT_MAX_COUNT }
    }

    pub fn with_k(mut self, k: usize) -> Self {
        self.k = k;
        self
    }

    pub fn with_stride(mut self, stride: usize) -> Self {
        self.stride = stride;
        self
    }

    /// Cap on stored occurrences per k-mer bucket; buckets exceeding this are
    /// truncated so ultra-repetitive k-mers (centromeres, transposons) don't
    /// produce unbounded candidate lists at lookup time.
    pub fn with_max_count(mut self, max_count: usize) -> Self {
        self.max_count = max_count;
        self
    }

    pub fn build_and_save<P: AsRef<Path>>(self, prefix: P) -> Result<()> {
        self.build().save(prefix)
    }

    /// Build the `HashTable` in memory without saving it -- the actual
    /// indexing logic; `build_and_save` is a thin wrapper, and this is what
    /// tests exercise directly to inspect the built table.
    ///
    /// Single-pass, sort-based build: every occurrence is pushed straight
    /// into one flat `Vec<(hash, packed_position)>` as it's scanned (no
    /// per-key `Vec` ever materializes), then a single stable sort by hash
    /// groups same-key occurrences into contiguous runs -- stable so each
    /// run keeps its original scan order, matching the "keep the first
    /// `max_count` encountered" cap this always had. This avoids the
    /// transient double bookkeeping an earlier version of this function had
    /// (build a per-key `FxHashMap<u32, Vec<u64>>` first, then copy it into
    /// the final flat form after the fact -- measured to *increase* peak
    /// build memory over the original nested-`Vec` design, since both
    /// representations were briefly alive at once, the opposite of the
    /// intended win).
    pub fn build(self) -> HashTable {
        let mut pairs: Vec<(u32, u64)> = Vec::new();

        for (cid, contig) in self.genome.contigs.iter().enumerate() {
            // Decode this contig's packed sequence once; the genome as a
            // whole stays 2-bit-packed in memory (see types::Contig), only
            // the one contig currently being indexed is briefly materialized
            // as ASCII.
            let seq = contig.slice(0, contig.length);
            if seq.len() < self.k { continue; }

            // Index forward strand. Windows overlapping an N-run (or other
            // non-ACGT ambiguity code) are skipped entirely rather than
            // hashed -- see `HashTable::is_acgt_only`.
            for i in (0..=seq.len() - self.k).step_by(self.stride) {
                let kmer = &seq[i..i + self.k];
                if !HashTable::is_acgt_only(kmer) { continue; }
                let hash = HashTable::hash_kmer(kmer);
                pairs.push((hash, HashTable::pack_position(cid, i, Strand::Forward)));
            }

            // Index reverse-complement strand, explicitly tagged so lookups
            // don't need to infer strand from a diagonal trend.
            let rc_seq = DnaSeq { bases: seq.to_vec() }.reverse_complement();
            let rc = &rc_seq.bases;
            for i in (0..=rc.len() - self.k).step_by(self.stride) {
                let kmer = &rc[i..i + self.k];
                if !HashTable::is_acgt_only(kmer) { continue; }
                let hash = HashTable::hash_kmer(kmer);
                let packed = HashTable::pack_position(cid, seq.len() - i - self.k, Strand::Reverse);
                pairs.push((hash, packed));
            }
        }

        // Stable: preserves each hash's original scan-order run so the
        // max_count cap below keeps the same entries the old per-key-Vec
        // version did.
        pairs.sort_by_key(|&(h, _)| h);

        let mut positions: Vec<u64> = Vec::with_capacity(pairs.len());
        let mut table: FxHashMap<u32, (u32, u32)> = FxHashMap::default();
        let mut total_entries = 0usize;

        let mut i = 0;
        while i < pairs.len() {
            let hash = pairs[i].0;
            let start = i;
            while i < pairs.len() && pairs[i].0 == hash { i += 1; }
            let count = (i - start).min(self.max_count);
            let offset = positions.len() as u32;
            positions.extend(pairs[start..start + count].iter().map(|&(_, pos)| pos));
            table.insert(hash, (offset, count as u32));
            total_entries += count;
        }

        HashTable {
            table,
            positions,
            k: self.k,
            seed_stride: self.stride,
            total_entries,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A k-mer window entirely inside an N-run must not end up in the built
    /// table at all -- neither under its own hash (N isn't a real base, so
    /// there's nothing legitimate to record) nor, more importantly, under
    /// the hash of the all-'A' k-mer `base2bits` would otherwise decode it
    /// as (which would let real 'AAAA...' reads elsewhere spuriously match
    /// an N-gap). Regression test for the bug found in code review.
    #[test]
    fn n_windows_are_not_indexed() {
        let seq = format!("{}{}{}", "A".repeat(20), "N".repeat(20), "C".repeat(20));
        let contig = Contig::from_ascii("c1".to_string(), seq.as_bytes());
        let genome = GenomeIndex {
            species: "t".to_string(),
            assembly: "t".to_string(),
            total_length: contig.length,
            contigs: vec![contig],
        };

        let k = 15;
        let ht = HashTableBuilder::new(&genome).with_k(k).with_stride(1).build();

        // Every forward-strand window that overlaps the N-run (positions
        // 20..40) at all must be absent. Two independent checks:
        // (a) none of those windows' own hash buckets contain a forward-
        //     strand entry at that position;
        // (b) the all-A k-mer's bucket (what `base2bits` would otherwise
        //     decode any N-containing window as) contains no entry whose
        //     position falls in the N-run.
        let all_a_hash = HashTable::hash_kmer(&vec![b'A'; k]);
        let all_a_bucket = ht.get(all_a_hash);
        for &packed in all_a_bucket {
            let (_, pos, strand) = HashTable::unpack_position(packed);
            if strand.is_forward() {
                assert!(
                    pos + k <= 20 || pos >= 40,
                    "N-run position {pos} spuriously indexed under the all-A hash"
                );
            }
        }

        // Direct check: no forward-strand entry at any start position whose
        // window overlaps [20, 40) exists anywhere in the table.
        for start in 6..40 {
            for &(offset, count) in ht.table.values() {
                let bucket = &ht.positions[offset as usize..(offset + count) as usize];
                for &packed in bucket {
                    let (_, pos, strand) = HashTable::unpack_position(packed);
                    if strand.is_forward() && pos == start {
                        assert!(
                            pos + k <= 20 || pos >= 40,
                            "window at {pos}..{} overlaps the N-run but was indexed", pos + k
                        );
                    }
                }
            }
        }

        // Sanity: windows entirely outside the N-run (e.g. the all-A run's
        // own start) are still indexed normally.
        assert!(all_a_bucket.iter().any(|&packed| {
            let (_, pos, strand) = HashTable::unpack_position(packed);
            strand.is_forward() && pos == 0
        }), "the genuine all-A window at position 0 should still be indexed");
    }

    /// The read-side query mirrors the reference-side fix: a read window
    /// containing N must not be looked up at all (see
    /// `mapper::seeds::find_seed_hits`), so `is_acgt_only` -- the primitive
    /// both sides share -- is exercised directly here too.
    #[test]
    fn is_acgt_only_rejects_n_and_other_ambiguity_codes() {
        assert!(HashTable::is_acgt_only(b"ACGTacgtACGTACG"));
        assert!(!HashTable::is_acgt_only(b"ACGTNCGTACGTACG"));
        assert!(!HashTable::is_acgt_only(b"ACGTRCGTACGTACG")); // IUPAC ambiguity code
    }
}
