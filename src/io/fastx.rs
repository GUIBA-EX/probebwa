//! FASTA / FASTQ parsing (supports .gz).
//!
//! Quality encoding: FASTQ historically shipped in two incompatible
//! encodings — Sanger/Illumina 1.8+ (Phred+33) and the older Solexa/Illumina
//! 1.0-1.3 (Phred+64). Defaulting to Phred+33 while silently ignoring the
//! possibility of Phred+64 input would misparse older data, so the offset is
//! a parameter threaded from the CLI down to the parser instead of a
//! hardcoded constant.

use crate::types::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use flate2::read::MultiGzDecoder;
use anyhow::Result;

/// Phred+33 (Sanger/Illumina 1.8+), the modern FASTQ default.
pub const PHRED33_OFFSET: u8 = 33;
/// Phred+64 (Solexa/Illumina 1.0-1.3), selected via `--phred64`.
pub const PHRED64_OFFSET: u8 = 64;

/// Read one or more input files (FASTQ or FASTA) into `Read`s.
pub fn read_reads_files<P: AsRef<Path>>(paths: &[P], format: InputFormat, quality_offset: u8) -> Result<Vec<Read>> {
    read_reads_files_with_preprocessing(paths, format, quality_offset, &ReadPreprocessing::default())
}

/// Read-ID and read-content preprocessing applied as each record is parsed
/// (`--casava8`, `--adapter-strip`).
#[derive(Clone, Default)]
pub struct ReadPreprocessing {
    /// CASAVA 1.8+ headers (`@INST:RUN:FLOWCELL:LANE:TILE:X:Y
    /// READNUM:FILTER:CONTROL:INDEX`) put the mate number in a
    /// space-separated field rather than a `/1`/`/2` suffix on the ID
    /// itself, which this parser's plain whitespace-split ID extraction
    /// already handles correctly either way. The one thing that differs by
    /// format is the *older* Illumina convention's `/1`/`/2` suffix directly
    /// on the ID — harmless to strip unconditionally (CASAVA 1.8+ IDs never
    /// have it), so `--casava8` doesn't gate any different behavior here;
    /// it's accepted as a familiar flag name for users coming from other
    /// short-read mapping tools.
    pub casava8: bool,
    /// Adapter sequence to trim from the read's 3' end (`--adapter-strip`),
    /// if the read's tail overlaps the adapter's start with at most one
    /// mismatch per 10 bases of overlap and at least `MIN_ADAPTER_OVERLAP`
    /// bases of overlap.
    pub adapter: Option<Vec<u8>>,
}

/// Read one or more input files (FASTQ or FASTA) into `Read`s, applying
/// `--casava8`/`--adapter-strip` preprocessing to each record.
pub fn read_reads_files_with_preprocessing<P: AsRef<Path>>(
    paths: &[P],
    format: InputFormat,
    quality_offset: u8,
    preprocessing: &ReadPreprocessing,
) -> Result<Vec<Read>> {
    let mut reads = Vec::new();
    for path in paths {
        match format {
            InputFormat::Fastq => reads.extend(read_fastq(path, quality_offset)?),
            InputFormat::Fasta => reads.extend(read_fasta_as_reads(path)?),
        }
    }
    for read in &mut reads {
        strip_mate_suffix(&mut read.id);
        if let Some(adapter) = &preprocessing.adapter {
            strip_adapter(read, adapter);
        }
    }
    Ok(reads)
}

/// Strip a trailing `/1` or `/2` (pre-CASAVA-1.8 Illumina mate-number
/// convention) from a read ID — modern SAM QNAME convention doesn't carry
/// it (mate number is FLAG bits 0x40/0x80 instead), and leaving it in place
/// makes the same physical read pair show up under two different QNAMEs.
fn strip_mate_suffix(id: &mut String) {
    if let Some(stripped) = id.strip_suffix("/1").or_else(|| id.strip_suffix("/2")) {
        id.truncate(stripped.len());
    }
}

/// Minimum overlap (bases) between a read's 3' end and the adapter's start
/// for `strip_adapter` to trim it — shorter apparent overlaps are treated as
/// coincidental rather than real adapter read-through.
const MIN_ADAPTER_OVERLAP: usize = 6;

/// Trim `adapter` from `read`'s 3' end in place, if found: scans overlap
/// lengths from the read's full length down to `MIN_ADAPTER_OVERLAP`,
/// comparing the read's last `overlap` bases against the adapter's first
/// `overlap` bases, allowing up to 1 mismatch per 10 bases (rounded down,
/// minimum 0) — a real sequencing-error-tolerant adapter match, not an exact
/// one. The longest qualifying overlap wins (checked longest-first, so a
/// genuine long read-through isn't missed in favor of a shorter incidental
/// match). Trims both the sequence and quality track to keep them the same
/// length.
fn strip_adapter(read: &mut Read, adapter: &[u8]) {
    let seq = &read.sequence.bases;
    let read_len = seq.len();
    let max_overlap = read_len.min(adapter.len());
    if max_overlap < MIN_ADAPTER_OVERLAP {
        return;
    }
    for overlap in (MIN_ADAPTER_OVERLAP..=max_overlap).rev() {
        let read_tail = &seq[read_len - overlap..];
        let adapter_head = &adapter[..overlap];
        let mismatches = read_tail.iter().zip(adapter_head)
            .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
            .count();
        if mismatches <= overlap / 10 {
            let keep = read_len - overlap;
            read.sequence.bases.truncate(keep);
            read.qualities.scores.truncate(keep);
            return;
        }
    }
}

/// Read one or more FASTQ files (Phred+33).
pub fn read_fastq_files<P: AsRef<Path>>(paths: &[P]) -> Result<Vec<Read>> {
    let mut reads = Vec::new();
    for path in paths {
        reads.extend(read_fastq(path, PHRED33_OFFSET)?);
    }
    Ok(reads)
}

fn open_maybe_gz<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>> {
    let file = File::open(&path)?;
    let is_gz = path.as_ref().extension().and_then(|s| s.to_str()) == Some("gz");
    Ok(if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

/// Read a single FASTQ file (plain or gzipped) with the given quality offset.
pub fn read_fastq<P: AsRef<Path>>(path: P, quality_offset: u8) -> Result<Vec<Read>> {
    let reader = open_maybe_gz(path)?;

    let mut parser = FastqParser::new(reader);
    let mut reads = Vec::new();
    while let Some(record) = parser.next() {
        let rec = record?;
        let id = String::from_utf8_lossy(&rec.id).to_string();
        let seq = DnaSeq::from_ascii(&rec.seq);
        let qual = QualityScores::from_phred_ascii(&rec.qual, quality_offset);
        reads.push(Read { id, sequence: seq, qualities: qual, is_reverse: false });
    }
    Ok(reads)
}

/// Default quality assigned to FASTA-derived reads, which carry no quality
/// track of their own — a flat high-confidence value is the standard
/// convention wherever a tool accepts FASTA as read input.
const FASTA_DEFAULT_QUAL: u8 = 40;

/// Read a single FASTA file as reads (no quality track — filled with a flat
/// high-confidence default).
pub fn read_fasta_as_reads<P: AsRef<Path>>(path: P) -> Result<Vec<Read>> {
    let reader = open_maybe_gz(path)?;

    let mut parser = FastaAsReadsParser::new(reader);
    let mut reads = Vec::new();
    while let Some(record) = parser.next() {
        let rec = record?;
        let id = String::from_utf8_lossy(&rec.id).to_string();
        let seq = DnaSeq::from_ascii(&rec.seq);
        let qual = QualityScores { scores: vec![FASTA_DEFAULT_QUAL; seq.len()] };
        reads.push(Read { id, sequence: seq, qualities: qual, is_reverse: false });
    }
    Ok(reads)
}

// ---------- Minimal FASTA-as-reads parser ----------
struct FastaAsReadsParser<R: BufRead> {
    reader: R,
    line: String,
    /// A header line read while scanning the *previous* record's sequence
    /// (i.e. the line that ended it) gets carried over here instead of
    /// discarded, so every record is actually returned. Losing that line
    /// used to silently drop every other record in a multi-record FASTA
    /// file — confirmed by comparing against another tool on the same
    /// multi-record test input.
    pending_header: Option<String>,
}

struct FastaRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
}

impl<R: BufRead> FastaAsReadsParser<R> {
    fn new(reader: R) -> Self {
        Self { reader, line: String::new(), pending_header: None }
    }

    fn next(&mut self) -> Option<Result<FastaRecord>> {
        let header_line = if let Some(h) = self.pending_header.take() {
            h
        } else {
            loop {
                self.line.clear();
                match self.reader.read_line(&mut self.line) {
                    Ok(0) => return None,
                    Ok(_) => {
                        if self.line.starts_with('>') { break; }
                    }
                    Err(e) => return Some(Err(e.into())),
                }
            }
            std::mem::take(&mut self.line)
        };
        let id = header_line.trim_start_matches('>')
            .split_whitespace().next().unwrap_or("").as_bytes().to_vec();

        let mut seq = Vec::new();
        loop {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => break,
                Ok(_) => {
                    let line = self.line.trim();
                    if line.starts_with('>') {
                        self.pending_header = Some(self.line.clone());
                        break;
                    }
                    seq.extend(line.bytes());
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
        Some(Ok(FastaRecord { id, seq }))
    }
}

// ---------- Minimal FASTQ parser ----------
struct FastqParser<R: BufRead> {
    reader: R,
    buf: String,
}

struct FastqRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl<R: BufRead> FastqParser<R> {
    fn new(reader: R) -> Self {
        Self { reader, buf: String::new() }
    }

    fn next(&mut self) -> Option<Result<FastqRecord>> {
        // ID line (@...)
        self.buf.clear();
        match self.reader.read_line(&mut self.buf) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e.into())),
        }
        let id = self.buf.trim_start_matches('@')
            .split_whitespace().next().unwrap_or("").as_bytes().to_vec();

        // Sequence line(s)
        self.buf.clear();
        let mut seq = Vec::new();
        loop {
            match self.reader.read_line(&mut self.buf) {
                Ok(0) => break,
                Ok(_) => {
                    let line = self.buf.trim();
                    if line.starts_with('+') {
                        break;
                    }
                    seq.extend(line.bytes());
                    self.buf.clear();
                }
                Err(e) => return Some(Err(e.into())),
            }
        }

        // Quality line(s) — must match sequence length. Read at least one
        // line unconditionally (checking the length only *after*), so a
        // record with a zero-length sequence still consumes its (blank)
        // quality line instead of leaving it in the stream to be misread as
        // the next record's ID line, desyncing every record after it.
        self.buf.clear();
        let mut qual = Vec::new();
        loop {
            match self.reader.read_line(&mut self.buf) {
                Ok(0) => break,
                Ok(_) => {
                    qual.extend(self.buf.trim().bytes());
                    self.buf.clear();
                    if qual.len() >= seq.len() { break; }
                }
                Err(e) => return Some(Err(e.into())),
            }
        }

        Some(Ok(FastqRecord { id, seq, qual }))
    }
}
