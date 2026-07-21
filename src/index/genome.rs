//! Genome index (`.stidx`) — reference sequence storage.

use crate::types::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use anyhow::Result;
use flate2::read::MultiGzDecoder;

/// 2-bit code → base, inverse of the packing done in `Contig::from_ascii`.
const BITS2BASE: [u8; 4] = *b"ACGT";

impl Contig {
    /// Pack raw ASCII sequence into 2 bits/base, recording any non-ACGT runs
    /// (N, ambiguity codes) separately so they can be restored on decode.
    pub fn from_ascii(name: String, seq: &[u8]) -> Self {
        let length = seq.len();
        let mut packed = vec![0u8; length.div_ceil(4)];
        let mut n_runs = Vec::new();
        let mut run_start: Option<usize> = None;

        for (pos, &b) in seq.iter().enumerate() {
            let (bits, is_acgt) = match b {
                b'A' | b'a' => (0u8, true),
                b'C' | b'c' => (1u8, true),
                b'G' | b'g' => (2u8, true),
                b'T' | b't' => (3u8, true),
                _ => (0u8, false),
            };
            packed[pos / 4] |= bits << ((pos % 4) * 2);
            if is_acgt {
                if let Some(s) = run_start.take() {
                    n_runs.push((s as u32, pos as u32));
                }
            } else if run_start.is_none() {
                run_start = Some(pos);
            }
        }
        if let Some(s) = run_start.take() {
            n_runs.push((s as u32, length as u32));
        }

        Self { name, packed, n_runs, length }
    }

    /// Decode the half-open range `[start, end)` back to ASCII bases.
    pub fn slice(&self, start: usize, end: usize) -> Vec<u8> {
        let start = start.min(self.length);
        let end = end.min(self.length);
        if start >= end {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(end - start);
        for pos in start..end {
            let byte = self.packed[pos / 4];
            let bits = (byte >> ((pos % 4) * 2)) & 0b11;
            out.push(BITS2BASE[bits as usize]);
        }
        for &(s, e) in &self.n_runs {
            let os = (s as usize).max(start);
            let oe = (e as usize).min(end);
            for p in os..oe {
                out[p - start] = b'N';
            }
        }
        out
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GenomeIndex {
    pub species: String,
    pub assembly: String,
    pub contigs: Vec<Contig>,
    pub total_length: usize,
}

impl GenomeIndex {
    pub fn load<P: AsRef<Path>>(prefix: P) -> Result<Self> {
        let path = format!("{}.stidx", prefix.as_ref().display());
        let file = File::open(&path)?;
        let reader = BufReader::new(file);
        let index: GenomeIndex = bincode::deserialize_from(reader)?;
        Ok(index)
    }

    pub fn save<P: AsRef<Path>>(&self, prefix: P) -> Result<()> {
        let path = format!("{}.stidx", prefix.as_ref().display());
        let file = File::create(&path)?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, self)?;
        Ok(())
    }

    pub fn get_sequence(&self, contig_id: usize, start: usize, end: usize) -> Vec<u8> {
        self.contigs[contig_id].slice(start, end)
    }
}

pub struct GenomeIndexBuilder {
    species: String,
    assembly: String,
    contigs: Vec<Contig>,
}

impl GenomeIndexBuilder {
    pub fn new(species: &str, assembly: &str) -> Self {
        Self {
            species: species.to_string(),
            assembly: assembly.to_string(),
            contigs: Vec::new(),
        }
    }

    pub fn add_fasta<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let file = File::open(&path)?;
        let reader: Box<dyn BufRead> = if path.as_ref().extension()
            .and_then(|s| s.to_str()) == Some("gz")
        {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut fasta_reader = FastaReader::new(reader);
        while let Some(record) = fasta_reader.next() {
            let record = record?;
            let name = String::from_utf8_lossy(record.id())
                .split_whitespace().next().unwrap_or("").to_string();
            self.contigs.push(Contig::from_ascii(name, record.seq()));
        }
        Ok(())
    }

    pub fn build_and_save<P: AsRef<Path>>(self, prefix: P) -> Result<()> {
        let total_length = self.contigs.iter().map(|c| c.length).sum();
        let index = GenomeIndex {
            species: self.species,
            assembly: self.assembly,
            contigs: self.contigs,
            total_length,
        };
        index.save(prefix)?;
        Ok(())
    }
}

// ---------- Minimal FASTA parser ----------
struct FastaReader<R: BufRead> {
    reader: R,
    line: String,
    /// A header line read while scanning the *previous* record's sequence
    /// (i.e. the line that ended it) gets carried over here instead of
    /// discarded, so every record — every contig — actually gets indexed.
    /// Losing that line used to silently drop every other contig in a
    /// multi-contig reference FASTA.
    pending_header: Option<String>,
}

struct FastaRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
}

impl FastaRecord {
    fn id(&self) -> &[u8] { &self.id }
    fn seq(&self) -> &[u8] { &self.seq }
}

impl<R: BufRead> FastaReader<R> {
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
