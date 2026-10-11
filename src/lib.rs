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
/// order, with matching normalized IDs required); any other count is single-ended.
pub fn map_reads_with_format<P: AsRef<Path>>(
    genome_prefix: P,
    hash_prefix: P,
    read_files: &[P],
    input_format: InputFormat,
    quality_offset: u8,
    options: MapperOptions,
) -> Result<Vec<AlignmentRecord>> {
    anyhow::ensure!(!read_files.is_empty(), "at least one read file is required");
    anyhow::ensure!(options.substitution_rate.is_finite() && (0.0..=1.0).contains(&options.substitution_rate),
        "substitution rate must be finite and between 0 and 1 (inclusive)");
    anyhow::ensure!(options.insert_size_mean.is_finite(), "insert size mean must be finite");
    anyhow::ensure!(options.insert_size_sd.is_finite() && options.insert_size_sd > 0.0,
        "insert size SD must be finite and positive");
    let genome = index::GenomeIndex::load(genome_prefix)?;
    let hash_table = index::HashTable::load(hash_prefix)?;

    let read_files: Vec<&Path> = read_files.iter().map(AsRef::as_ref).collect();
    let map_all = || -> Result<Vec<AlignmentRecord>> {
        if read_files.len() == 2 {
            return map_paired(&genome, &hash_table, &read_files[0], &read_files[1], input_format, quality_offset, &options);
        }

        let preprocessing = io::ReadPreprocessing { casava8: options.casava8, adapter: options.adapter_strip.clone() };
        let reads = io::read_reads_files_with_preprocessing(&read_files, input_format, quality_offset, &preprocessing)?;

        use rayon::prelude::*;
        let single_mapper = mapper::SingleEndMapper::new(&genome, &hash_table, &options);
        Ok(reads.par_iter()
            .map(|read| single_mapper.map(read).unwrap_or_else(|_| mapper::single::unmapped_record(read)))
            .collect())
    };

    // `--threads` controls parallelism: run inside a scoped pool sized to
    // `options.threads` rather than leaking a global-pool configuration.
    // threads >= 1 pins the pool to exactly that many workers (so `-t 1` is
    // genuinely serial, as the CLI's default implies); threads == 0 means
    // "let rayon use all available cores".
    if options.threads >= 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(options.threads)
            .build()?;
        pool.install(map_all)
    } else {
        map_all()
    }
}

/// Bound candidate memory while seeding pairs in parallel. Pair decisions
/// and library updates remain in input order, independent of chunk size.
const PAIR_BATCH_SIZE: usize = 2000;

/// Paired mapping with parallel seeding and sequential model updates.
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

    for (i, (read1, read2)) in reads1.iter().zip(&reads2).enumerate() {
        anyhow::ensure!(read1.id == read2.id,
            "paired-end read IDs differ at pair {}: '{}' vs '{}'", i + 1, read1.id, read2.id);
    }

    let mut paired_mapper = mapper::PairedEndMapper::new(genome, hash_table, options);
    let mut records = Vec::with_capacity(reads1.len() * 2);

    let pairs: Vec<ReadPair> = reads1.into_iter().zip(reads2)
        .map(|(read1, read2)| ReadPair { read1, read2 })
        .collect();

    for batch in pairs.chunks(PAIR_BATCH_SIZE) {
        for (r1, r2) in paired_mapper.map_pairs_batch(batch) {
            records.push(r1);
            records.push(r2);
        }
    }

    Ok(records)
}
