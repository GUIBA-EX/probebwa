//! Scratch profiling harness (not part of the test suite): times
//! `SingleEndMapper::map_candidates` for both mates separately from the full
//! `PairedEndMapper::map_pair` call, over a real FASTQ pair, to isolate how
//! much of paired-mapping wall time is the shared seed+align work (also paid
//! by single-end mapping) vs. the shortlist/cross-alignment machinery added
//! on top of it this session.

use probebwa::index::{GenomeIndex, HashTable};
use probebwa::io::read_reads_files;
use probebwa::mapper::{PairedEndMapper, SingleEndMapper};
use probebwa::{InputFormat, MapperOptions, ReadPair};
use std::time::{Duration, Instant};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let genome_prefix = &args[1];
    let hash_prefix = &args[2];
    let fq1 = &args[3];
    let fq2 = &args[4];

    let genome = GenomeIndex::load(genome_prefix).unwrap();
    let hash_table = HashTable::load(hash_prefix).unwrap();
    let opts = MapperOptions::default();

    let reads1 = read_reads_files(&[fq1], InputFormat::Fastq, probebwa::io::PHRED33_OFFSET).unwrap();
    let reads2 = read_reads_files(&[fq2], InputFormat::Fastq, probebwa::io::PHRED33_OFFSET).unwrap();

    let single_mapper = SingleEndMapper::new(&genome, &hash_table, &opts);
    let mut paired_mapper = PairedEndMapper::new(&genome, &hash_table, &opts);

    let mut candidates_time = Duration::ZERO;
    let mut full_pair_time = Duration::ZERO;
    let mut n = 0usize;

    for (r1, r2) in reads1.iter().zip(reads2.iter()) {
        let t0 = Instant::now();
        let _c1 = single_mapper.map_candidates(r1);
        let _c2 = single_mapper.map_candidates(r2);
        candidates_time += t0.elapsed();

        let pair = ReadPair { read1: r1.clone(), read2: r2.clone() };
        let t1 = Instant::now();
        let _ = paired_mapper.map_pair(&pair);
        full_pair_time += t1.elapsed();

        n += 1;
    }

    println!("pairs: {n}");
    println!("map_candidates (both mates, seed+cluster+align, no pairing logic): {:?} total, {:?}/pair", candidates_time, candidates_time / n as u32);
    println!("full map_pair (includes shortlist+cross-alignment+scoring):        {:?} total, {:?}/pair", full_pair_time, full_pair_time / n as u32);
    let overhead = full_pair_time.saturating_sub(candidates_time);
    println!("implied pairing-specific overhead (full - candidates, note: candidates work is *redone* inside map_pair too, so this is a lower bound on total, not a clean subtraction): {:?} total, {:?}/pair", overhead, overhead / n as u32);
}
