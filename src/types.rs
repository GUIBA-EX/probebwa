//! Core data structures.

use serde::{Deserialize, Serialize};

/// A nucleotide sequence (read or reference segment).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DnaSeq {
    pub bases: Vec<u8>,
}

impl DnaSeq {
    pub fn from_ascii(s: &[u8]) -> Self {
        let bases = s.iter().map(|&b| b.to_ascii_uppercase()).collect();
        Self { bases }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// Extract all overlapping k-mers of length k.
    pub fn kmers(&self, k: usize) -> Vec<&[u8]> {
        if self.len() < k {
            return vec![];
        }
        self.bases.windows(k).collect()
    }

    /// Generate the reverse-complement of this sequence.
    pub fn reverse_complement(&self) -> Self {
        let mut rc = Vec::with_capacity(self.bases.len());
        for &base in self.bases.iter().rev() {
            rc.push(match base {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                b'N' => b'N',
                _ => b'N',
            });
        }
        Self { bases: rc }
    }
}

/// Quality scores (Phred scale).
#[derive(Clone, Debug)]
pub struct QualityScores {
    pub scores: Vec<u8>,
}

impl QualityScores {
    pub fn from_phred_ascii(ascii: &[u8], offset: u8) -> Self {
        let scores = ascii.iter().map(|&c| c.saturating_sub(offset)).collect();
        Self { scores }
    }

    /// Error probability at position i: P(error) = 10^(-Q/10).
    pub fn error_probability(&self, pos: usize) -> f64 {
        if pos >= self.scores.len() {
            return 0.01; // default when quality unavailable
        }
        let q = self.scores[pos] as f64;
        10f64.powf(-q / 10.0)
    }

    /// Probability that the base is *correct*.
    pub fn correct_probability(&self, pos: usize) -> f64 {
        1.0 - self.error_probability(pos)
    }
}

/// A sequencing read (single-end or one mate).
#[derive(Clone, Debug)]
pub struct Read {
    pub id: String,
    pub sequence: DnaSeq,
    pub qualities: QualityScores,
    pub is_reverse: bool,
}

/// A paired-end read pair.
#[derive(Clone, Debug)]
pub struct ReadPair {
    pub read1: Read,
    pub read2: Read,
}

/// Reference contig / chromosome.
///
/// Sequence is stored 2-bit-packed (4 bases/byte) rather than as raw ASCII —
/// a 4x memory reduction and better cache locality for the random-access
/// window fetches mapping does, mirroring how production aligners store the
/// reference (e.g. minimap2's packed `mm_idx_t::S`, minibwa's `l2bit`
/// format). Positions that weren't plain A/C/G/T (N, ambiguity codes, gaps)
/// can't be represented in 2 bits, so their packed bits are a placeholder and
/// the true runs are kept separately in `n_runs`, restored to 'N' on decode
/// (see `Contig::slice` in index::genome). This is lossy for anything other
/// than N — case (soft-masking) and specific IUPAC ambiguity codes are not
/// round-tripped — which is acceptable here since neither is used elsewhere
/// in this codebase (sequences are uppercased on ingest regardless).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Contig {
    pub name: String,
    pub packed: Vec<u8>,
    pub n_runs: Vec<(u32, u32)>,
    pub length: usize,
}

/// Genome index (`.stidx`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GenomeIndex {
    pub species: String,
    pub assembly: String,
    pub contigs: Vec<Contig>,
    pub total_length: usize,
}

/// A seed hit: k-mer match between read and reference.
///
/// `hit_strand` comes directly from the hash table entry (see
/// `index::HashTable::unpack_position`): `Forward` means the read's k-mer, as
/// given, equals the genome's forward-strand k-mer at `ref_pos`; `Reverse`
/// means it equals the genome's reverse-complement k-mer there, i.e. the read
/// aligns to the minus strand at this locus.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SeedHit {
    pub read_pos: usize,
    pub ref_contig: usize,
    pub ref_pos: usize,
    pub k: usize,
    pub hit_strand: Strand,
}

/// Candidate mapping location.
#[derive(Clone, Debug)]
pub struct Candidate {
    pub contig_id: usize,
    pub position: usize,
    pub strand: Strand,
    pub score: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn is_forward(&self) -> bool {
        matches!(self, Strand::Forward)
    }

    pub fn complement(&self) -> Self {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
        }
    }
}

/// CIGAR operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),      // M
    Insertion(u32),  // I
    Deletion(u32),   // D
    SoftClip(u32),   // S
    HardClip(u32),   // H
}

// Display impl for CigarOp lives in align::cigar.

/// SAM-like alignment record.
#[derive(Clone, Debug)]
pub struct AlignmentRecord {
    pub read_id: String,
    pub contig_name: String,
    pub position: usize, // 1-based
    pub mapq: u8,
    pub cigar: Vec<CigarOp>,
    pub strand: Strand,
    pub read_sequence: Vec<u8>,
    pub read_qualities: Vec<u8>,
    pub is_paired: bool,
    pub is_proper_pair: bool,
    pub mate_contig: Option<String>,
    pub mate_position: Option<usize>,
    pub insert_size: Option<i64>,
    pub alignment_score: f64,
    pub posterior_prob: f64,
    /// 1 or 2 for paired-end reads (SAM FLAG 0x40/0x80); None for single-end.
    pub read_number: Option<u8>,
    /// Strand the mate aligned to, if the mate is mapped (SAM FLAG 0x20).
    pub mate_strand: Option<Strand>,
    /// True if this read is paired but its mate did not map (SAM FLAG 0x8).
    pub mate_unmapped: bool,
}

/// Mapper runtime options.
#[derive(Clone, Debug)]
pub struct MapperOptions {
    pub substitution_rate: f64,
    pub insert_size_mean: f64,
    pub insert_size_sd: f64,
    pub max_indel_len: usize,
    pub threads: usize,
    pub output_format: OutputFormat,
    /// Gap-open penalty, Phred scale (`--gapopen`).
    pub gap_open_phred: u32,
    /// Gap-extend penalty, Phred scale (`--gapextend`).
    pub gap_extend_phred: u32,
    /// Prior probability (Phred scale) that a discordant/anomalous pair
    /// separation reflects a real structural variant rather than a mapping
    /// error, used as a floor under the Gaussian insert-size likelihood when
    /// scoring paired-end candidates (`--svprior`; default Phred 55, i.e.
    /// probability 3e-6).
    pub sv_prior_phred: f64,
    /// Accept `/1`/`/2`-suffixed pre-CASAVA-1.8 read IDs (`--casava8`) --
    /// see `io::fastx::ReadPreprocessing` for why this doesn't actually gate
    /// different parsing behavior (the suffix is always stripped).
    pub casava8: bool,
    /// Adapter sequence to trim from each read's 3' end before mapping
    /// (`--adapter-strip=SEQ`); see `io::fastx::strip_adapter`.
    pub adapter_strip: Option<Vec<u8>>,
}

impl Default for MapperOptions {
    fn default() -> Self {
        Self {
            substitution_rate: 0.001,
            insert_size_mean: 250.0,
            insert_size_sd: 60.0,
            max_indel_len: 30,
            threads: 1,
            output_format: OutputFormat::Sam,
            gap_open_phred: 40,
            gap_extend_phred: 3,
            sv_prior_phred: 55.0,
            casava8: false,
            adapter_strip: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    Sam,
}

/// Read input format, selected via `--inputformat=fasta|fastq`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InputFormat {
    Fastq,
    Fasta,
}

/// Insert-size distribution (learned online from proper pairs).
/// Impl (`new`/`log_likelihood`/`update`) lives in mapq::insert_size.
#[derive(Clone, Debug)]
pub struct InsertSizeDistribution {
    pub mean: f64,
    pub sd: f64,
    pub n_observed: usize,
}
