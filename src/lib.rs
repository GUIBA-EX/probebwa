//! probebwa — Rust mapper for aligning UCE (ultra-conserved element)
//! bait/probe sequences to a reference genome.
//!
//! Implements a dense seed-and-extend mapping algorithm:
//!   per-position exact + stride-sampled single-mismatch k-mer seeding,
//!   candidate clustering, banded affine-gap alignment, and Bayesian
//!   posterior MAPQ scoring.
//! See README.md for usage and current limitations.

pub mod index;
pub mod align;
pub mod io;
pub mod mapq;
pub mod mapper;
pub mod types;

pub use types::*;

use std::path::Path;
use anyhow::Result;

/// Build a genome index (`.stidx`) from FASTA reference files.
pub fn build_genome_index<P: AsRef<Path>>(
    species: &str,
    assembly: &str,
    output_prefix: P,
    fasta_files: &[P],
) -> Result<()> {
    let mut builder = index::GenomeIndexBuilder::new(species, assembly);
    for file in fasta_files {
        builder.add_fasta(file)?;
    }
    builder.build_and_save(output_prefix)?;
    Ok(())
}

/// Build a 15-mer hash table (`.sthash`) from a genome index.
pub fn build_hash_table<P: AsRef<Path>>(
    genome_prefix: P,
    output_prefix: P,
) -> Result<()> {
    let genome = index::GenomeIndex::load(genome_prefix)?;
    let builder = index::HashTableBuilder::new(&genome);
    builder.build_and_save(output_prefix)?;
    Ok(())
}

/// Map reads to the reference genome.
pub fn map_reads<P: AsRef<Path>>(
    genome_prefix: P,
    hash_prefix: P,
    read_files: &[P],
    options: MapperOptions,
) -> Result<Vec<AlignmentRecord>> {
    map_reads_with_format(genome_prefix, hash_prefix, read_files, InputFormat::Fastq, io::PHRED33_OFFSET, options)
}

/// Map reads to the reference genome, choosing input format and quality
/// encoding explicitly (`--inputformat`, `--phred64`).
///
/// Follows the common short-read aligner convention of "two read files means
/// paired-end": exactly two read files triggers paired-end mapping, with
/// mate 1 and mate 2 read from the two files in lockstep (paired by stream
/// order, not by matching read IDs); any other count is single-ended.
pub fn map_reads_with_format<P: AsRef<Path>>(
    genome_prefix: P,
    hash_prefix: P,
    read_files: &[P],
    input_format: InputFormat,
    quality_offset: u8,
    options: MapperOptions,
) -> Result<Vec<AlignmentRecord>> {
    let genome = index::GenomeIndex::load(genome_prefix)?;
    let hash_table = index::HashTable::load(hash_prefix)?;

    if read_files.len() == 2 {
        return map_paired(&genome, &hash_table, &read_files[0], &read_files[1], input_format, quality_offset, &options);
    }

    let preprocessing = io::ReadPreprocessing { casava8: options.casava8, adapter: options.adapter_strip.clone() };
    let reads = io::read_reads_files_with_preprocessing(read_files, input_format, quality_offset, &preprocessing)?;

    use rayon::prelude::*;
    let single_mapper = mapper::SingleEndMapper::new(&genome, &hash_table, &options);

    let map_all = || -> Vec<AlignmentRecord> {
        reads.par_iter()
            .map(|read| single_mapper.map(read).unwrap_or_else(|_| mapper::single::unmapped_record(read)))
            .collect()
    };

    // `--threads` controls parallelism: run inside a scoped pool sized to
    // `options.threads` rather than leaking a global-pool configuration.
    // threads >= 1 pins the pool to exactly that many workers (so `-t 1` is
    // genuinely serial, as the CLI's default implies); threads == 0 means
    // "let rayon use all available cores".
    let records = if options.threads >= 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(options.threads)
            .build()?;
        pool.install(map_all)
    } else {
        map_all()
    };

    Ok(records)
}

/// How many read pairs `map_paired` processes per parallel batch (see
/// `mapper::PairedEndMapper::map_pairs_batch`). Each batch freezes the
/// insert-size model at its current estimate for the whole batch, so this is
/// a tradeoff: larger batches mean more parallelism per `rayon` dispatch
/// (less scheduling overhead relative to work), but let the model go longer
/// between updates. 2000 is small relative to the number of observations
/// (typically thousands) it usually takes for the model to have converged
/// and stabilized in the first place, so the staleness this introduces is
/// negligible in practice for any run long enough for batching to matter --
/// *after* the warm-up prefix below, which is what actually keeps the first
/// batch from being scored against a completely unconverged model.
const PAIR_BATCH_SIZE: usize = 2000;

/// How many pairs `map_paired` processes strictly serially (via `map_pair`,
/// updating the insert-size model after every pair as it goes) before
/// switching to parallel batching. Comparing batched output against the old
/// fully-serial baseline on real E. coli data found that without this, the
/// *entire first batch* is scored against the model's untouched initial
/// prior (`--insertsize`/`--insertsd`, default 250/60) -- unlike every
/// subsequent batch, which at least starts from a model already refined by
/// every prior batch, the first batch has nothing to inherit from, so
/// batching's usual "stale by at most one batch" framing understates its
/// actual cost there. A serial warm-up prefix bounds that cold-start cost to
/// a small fixed prefix instead of a full batch, at a small fixed cost (500
/// pairs processed without the cross-pair parallelism batching provides,
/// negligible against a real run's total size).
const PAIR_WARMUP_SIZE: usize = 500;

/// Paired-end mapping path. A short prefix is processed serially to warm up
/// the insert-size model (see `PAIR_WARMUP_SIZE`); the remainder runs in
/// parallel batches (see `mapper::PairedEndMapper::map_pairs_batch`) -- the
/// one thing that can't be parallelized across an unbounded number of pairs
/// is the insert-size model's online learning
/// (`InsertSizeDistribution::update`), which needs a consistent,
/// serially-updated shared state; batching (after warm-up) bounds that to
/// "stale by at most one batch" instead of serializing the whole run.
fn map_paired<P: AsRef<Path>>(
    genome: &index::GenomeIndex,
    hash_table: &index::HashTable,
    file1: &P,
    file2: &P,
    input_format: InputFormat,
    quality_offset: u8,
    options: &MapperOptions,
) -> Result<Vec<AlignmentRecord>> {
    let preprocessing = io::ReadPreprocessing { casava8: options.casava8, adapter: options.adapter_strip.clone() };
    let reads1 = io::read_reads_files_with_preprocessing(std::slice::from_ref(file1), input_format, quality_offset, &preprocessing)?;
    let reads2 = io::read_reads_files_with_preprocessing(std::slice::from_ref(file2), input_format, quality_offset, &preprocessing)?;
    anyhow::ensure!(
        reads1.len() == reads2.len(),
        "paired-end input files have different read counts ({} vs {})",
        reads1.len(), reads2.len()
    );

    let mut paired_mapper = mapper::PairedEndMapper::new(genome, hash_table, options);
    let mut records = Vec::with_capacity(reads1.len() * 2);

    let pairs: Vec<ReadPair> = reads1.into_iter().zip(reads2)
        .map(|(read1, read2)| ReadPair { read1, read2 })
        .collect();

    let warmup_len = pairs.len().min(PAIR_WARMUP_SIZE);
    for pair in &pairs[..warmup_len] {
        let (r1, r2) = paired_mapper.map_pair(pair)?;
        records.push(r1);
        records.push(r2);
    }

    for batch in pairs[warmup_len..].chunks(PAIR_BATCH_SIZE) {
        for (r1, r2) in paired_mapper.map_pairs_batch(batch) {
            records.push(r1);
            records.push(r2);
        }
    }

    Ok(records)
}
