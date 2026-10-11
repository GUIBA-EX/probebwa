//! End-to-end mapping tests.
//!
//! These build a genome index + hash table from a synthetic reference on
//! disk (via tempfile) and map synthetic reads through the same public API
//! the CLI uses (`build_genome_index`, `build_hash_table`, `map_reads`),
//! then assert on the resulting `AlignmentRecord` fields. This exercises the
//! full seed -> cluster -> align -> score pipeline, not just individual
//! functions.

use probebwa::*;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn write_fasta(dir: &TempDir, name: &str, contig_name: &str, seq: &str) -> PathBuf {
    let path = dir.path().join(name);
    fs::write(&path, format!(">{}\n{}\n", contig_name, seq)).unwrap();
    path
}

fn write_fastq(dir: &TempDir, name: &str, records: &[(String, String)]) -> PathBuf {
    let path = dir.path().join(name);
    let mut content = String::new();
    for (id, seq) in records {
        content.push_str(&format!("@{}\n{}\n+\n{}\n", id, seq, "I".repeat(seq.len())));
    }
    fs::write(&path, content).unwrap();
    path
}

fn revcomp(s: &str) -> String {
    s.chars().rev().map(|c| match c {
        'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N',
    }).collect()
}

/// Build a `.stidx`/`.sthash` pair from a single-contig FASTA and return the
/// shared prefix.
fn build_index(dir: &TempDir, fasta: &PathBuf, prefix: &str) -> PathBuf {
    let genome_prefix = dir.path().join(prefix);
    build_genome_index("test", "t1", &genome_prefix, &[fasta]).unwrap();
    build_hash_table(&genome_prefix, &genome_prefix).unwrap();
    genome_prefix
}

/// Small deterministic PRNG so tests are reproducible without adding a `rand`
/// dependency.
struct Xorshift64(u64);
impl Xorshift64 {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn next_base(&mut self) -> char {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        BASES[(self.next_u64() % 4) as usize]
    }
    fn gen_seq(&mut self, len: usize) -> String {
        (0..len).map(|_| self.next_base()).collect()
    }
}

#[test]
fn test_forward_and_reverse_strand_mapping() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(42).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "ref");

    let fwd_read = genome_seq[40..80].to_string();
    let rev_read = revcomp(&genome_seq[200..240]);
    let fastq = write_fastq(&dir, "reads.fq", &[
        ("fwd".to_string(), fwd_read),
        ("rev".to_string(), rev_read),
    ]);

    let records = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 2);

    let fwd = &records[0];
    assert_eq!(fwd.contig_name, "chr1");
    assert_eq!(fwd.position, 41);
    assert_eq!(fwd.strand, Strand::Forward);
    assert_eq!(fwd.cigar, vec![CigarOp::Match(40)]);
    // Native Stampy calcmutprobs([30; 40], [1000; 40], -1) returns
    // Phred 87. This short-read case uses its overlapping-window model.
    assert_eq!(fwd.mapq, 87);

    let rev = &records[1];
    assert_eq!(rev.contig_name, "chr1");
    assert_eq!(rev.position, 201);
    assert_eq!(rev.strand, Strand::Reverse);
    assert_eq!(rev.cigar, vec![CigarOp::Match(40)]);
    assert_eq!(rev.mapq, 87);
}

#[test]
fn test_mismatch_deletion_and_unmapped() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(7).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "ref");

    let base_read = &genome_seq[40..80];

    let mut mismatch_chars: Vec<char> = base_read.chars().collect();
    mismatch_chars[10] = if mismatch_chars[10] == 'A' { 'C' } else { 'A' };
    let mismatch_read: String = mismatch_chars.into_iter().collect();

    // Drop reference bases [60,62) relative to genome (positions 20-21 of
    // the read window) to create a clean 2bp deletion.
    let del_read = format!("{}{}", &base_read[0..20], &base_read[22..]);

    let garbage = "T".repeat(40);

    let fastq = write_fastq(&dir, "reads.fq", &[
        ("mismatch".to_string(), mismatch_read),
        ("del".to_string(), del_read),
        ("garbage".to_string(), garbage),
    ]);

    let records = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 3);

    let mismatch = &records[0];
    assert_eq!(mismatch.position, 41);
    assert_eq!(mismatch.cigar, vec![CigarOp::Match(40)]);
    // See the comment on the equivalent assertion in
    // test_forward_and_reverse_strand_mapping for why a 40bp read caps
    // at 87 rather than 99.
    assert_eq!(mismatch.mapq, 87);

    // The two dropped bases can legitimately land in the CIGAR as one 2bp
    // deletion or two 1bp deletions a couple of matches apart (both score
    // identically here since the synthetic genome has no homopolymer bias
    // to break the tie), so check aggregate op totals rather than the exact
    // operation sequence.
    let del = &records[1];
    assert_eq!(del.position, 41);
    let (m, d, i) = del.cigar.iter().fold((0u32, 0u32, 0u32), |(m, d, i), op| match op {
        CigarOp::Match(n) => (m + n, d, i),
        CigarOp::Deletion(n) => (m, d + n, i),
        CigarOp::Insertion(n) => (m, d, i + n),
        _ => (m, d, i),
    });
    assert_eq!((m, d, i), (38, 2, 0), "expected 38 matched bases and a net 2bp deletion, got cigar {:?}", del.cigar);
    assert!(del.mapq > 0, "a clean 2bp-deletion alignment should not be rejected by the LRT");

    let garbage_rec = &records[2];
    assert_eq!(garbage_rec.position, 0);
    assert_eq!(garbage_rec.mapq, 0);
}

/// Match the native alignment choice, not the synthetic generating CIGAR.
#[test]
fn test_native_short_read_long_deletion_choice() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(555).gen_seq(500);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "ref");

    let window = &genome_seq[100..160]; // 60bp true window
    let deletion_len = 25;
    let read = format!("{}{}", &window[..25], &window[25 + deletion_len..]);
    let fastq = write_fastq(&dir, "reads.fq", &[("bigdel".to_string(), read)]);

    let records = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 1);

    let rec = &records[0];
    assert_eq!(rec.position, 101);
    // Native Stampy 1.0.32 on this exact input chooses a terminal insertion,
    // rather than the synthetic generating deletion: UQ=122, MAPQ=65.
    assert_eq!(rec.cigar, vec![CigarOp::Match(25), CigarOp::Insertion(10)]);
    assert_eq!(rec.mapq, 65);
    assert!(rec.sam_tags.contains(&"UQ:i:122".to_owned()));
}

#[test]
fn test_paired_end_proper_and_improper() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(99).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "ref");

    // Proper (FR) pair: mate1 forward at 41, mate2 the RC of the
    // forward-strand window at 201..241, so mate2 lands downstream on the
    // reverse strand.
    let m1 = genome_seq[40..80].to_string();
    let m2 = revcomp(&genome_seq[200..240]);
    let fq1 = write_fastq(&dir, "m1.fq", &[("pair1".to_string(), m1)]);
    let fq2 = write_fastq(&dir, "m2.fq", &[("pair1".to_string(), m2)]);

    let records = map_reads(&prefix, &prefix, &[&fq1, &fq2], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 2);
    let (r1, r2) = (&records[0], &records[1]);
    assert_eq!(r1.read_number, Some(1));
    assert_eq!(r2.read_number, Some(2));
    assert!(r1.is_proper_pair && r2.is_proper_pair, "FR pair within insert-size range should be proper");
    assert_eq!(r1.strand, Strand::Forward);
    assert_eq!(r2.strand, Strand::Reverse);
    // TLEN (SAM convention): from mate1's start (41) to mate2's end
    // (201 + 40 - 1 = 240) inclusive = 200 -- not the 160bp start-to-start
    // gap, which undercounts by mate2's own 40bp span.
    assert_eq!(r1.insert_size, Some(200));
    assert_eq!(r2.insert_size, Some(-200));

    // Improper (FF) pair: both mates forward at the same coordinates.
    let m1b = genome_seq[40..80].to_string();
    let m2b = genome_seq[200..240].to_string();
    let fq1b = write_fastq(&dir, "m1b.fq", &[("ffpair".to_string(), m1b)]);
    let fq2b = write_fastq(&dir, "m2b.fq", &[("ffpair".to_string(), m2b)]);

    let records_ff = map_reads(&prefix, &prefix, &[&fq1b, &fq2b], MapperOptions::default()).unwrap();
    assert_eq!(records_ff.len(), 2);
    assert!(!records_ff[0].is_proper_pair && !records_ff[1].is_proper_pair, "FF orientation must not be flagged proper");
}

/// Regression test for position packing: map reads whose true coordinate is
/// well past 65,535 bp. With the old 16-bit position field these wrapped and
/// mapped to the wrong place (or not at all); the whole synthetic-scale test
/// below stays under 20kb, so this is the one case that actually exercises a
/// chromosome-scale coordinate.
#[test]
fn test_large_contig_positions() {
    let dir = TempDir::new().unwrap();
    let genome_len = 80_000usize; // > 65_535, i.e. past a u16
    let genome_seq = Xorshift64(2024).gen_seq(genome_len);
    let fasta = write_fasta(&dir, "big.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "big");

    // A read sourced from a coordinate that overflows 16 bits.
    let start = 70_000usize;
    let read = genome_seq[start..start + 50].to_string();
    let fastq = write_fastq(&dir, "reads.fq", &[("deep".to_string(), read)]);

    let records = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].contig_name, "chr1");
    assert_eq!(records[0].position, start + 1, "position past 65535 must not wrap");
    // Stampy 1.0.32 on this same generated reference/read reports MAPQ 96,
    // not 99: combining the single-read posterior floor with missing-locus
    // confidence reduces the final single-end score.
    assert_eq!(records[0].mapq, 96);
}

#[test]
fn test_alignment_band_tracks_clipped_left_padding() {
    let dir = TempDir::new().unwrap();
    let genome = Xorshift64(7052).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome);
    let prefix = build_index(&dir, &fasta, "ref");
    let inputs: Vec<_> = [0, 1, 9, 10].into_iter().map(|start|
        (format!("edge{start}"), genome[start..start + 60].to_string())).collect();
    let fq = write_fastq(&dir, "reads.fq", &inputs);
    let opts = MapperOptions { max_indel_len: 0, ..Default::default() };
    let results = map_reads(&prefix, &prefix, &[&fq], opts).unwrap();
    for (record, start) in results.iter().zip([0, 1, 9, 10]) {
        assert_eq!(record.position, start + 1, "{}", record.read_id);
        assert_eq!(record.cigar, vec![CigarOp::Match(60)]);
    }
}

#[test]
fn test_pair_batching_preserves_sequential_library_updates() {
    use probebwa::mapper::PairedEndMapper;
    let dir = TempDir::new().unwrap();
    let sequence = Xorshift64(1049).gen_seq(3000);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &sequence);
    let prefix = build_index(&dir, &fasta, "ref");
    let genome = index::GenomeIndex::load(&prefix).unwrap();
    let hash = index::HashTable::load(&prefix).unwrap();
    let options = MapperOptions::default();
    let read = |id: String, bases: String| Read {
        id, sequence: DnaSeq { bases: bases.into_bytes() },
        qualities: QualityScores { scores: vec![30; 80] }, is_reverse: false,
    };
    let pairs: Vec<_> = (0..30).map(|i| {
        let start = 100 + i * 50;
        let length = 200 + (i * 17) % 200;
        ReadPair {
            read1: read(format!("p{i}"), sequence[start..start + 80].into()),
            read2: read(format!("p{i}"), revcomp(&sequence[start + length - 80..start + length])),
        }
    }).collect();
    let mut sequential = PairedEndMapper::new(&genome, &hash, &options);
    let expected: Vec<_> = pairs.iter().flat_map(|p| {
        let (a, b) = sequential.map_pair(p).unwrap();
        [io::sam::format_sam(&a, None), io::sam::format_sam(&b, None)]
    }).collect();
    for threads in [1, 4] {
        for chunk_size in [1, 7, 30] {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let mut batched = PairedEndMapper::new(&genome, &hash, &options);
            let actual: Vec<_> = pool.install(|| pairs.chunks(chunk_size).flat_map(|batch| {
                batched.map_pairs_batch(batch).into_iter().flat_map(|(a, b)|
                    [io::sam::format_sam(&a, None), io::sam::format_sam(&b, None)])
            }).collect());
            assert_eq!(actual, expected, "threads={threads}, chunk_size={chunk_size}");
        }
    }
}

#[test]
fn test_random_match_rejection_is_distinct_from_ambiguous_mapping() {
    let dir = TempDir::new().unwrap();
    let segment = Xorshift64(9876).gen_seq(400);
    let reference = segment.repeat(12);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &reference);
    let prefix = build_index(&dir, &fasta, "ref");
    let bases = &segment[50..150];
    let fastq = dir.path().join("reads.fq");
    fs::write(&fastq, format!("@zero\n{bases}\n+\n{}\n@repeat\n{bases}\n+\n{}\n",
        "!".repeat(100), "I".repeat(100))).unwrap();
    let records = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
    assert_eq!(records[0].position, 0);
    assert!(records[0].cigar.is_empty());
    assert_eq!(records[1].mapq, 0);
    assert_ne!(records[1].position, 0, "MAPQ zero alone must not unmap a non-random repeated hit");
}

/// Native rejection fixture: this short, highly divergent mate cannot seed.
#[test]
fn test_mate_rescue_rejects_native_low_confidence_case() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(2468).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "ref");

    let m1 = genome_seq[40..80].to_string();
    let true_m2 = revcomp(&genome_seq[200..240]);
    let mut m2_chars: Vec<char> = true_m2.chars().collect();
    let mut pos = 0;
    while pos < m2_chars.len() {
        let orig = m2_chars[pos];
        m2_chars[pos] = ['A', 'C', 'G', 'T'].into_iter().find(|&b| b != orig).unwrap();
        pos += 6;
    }
    let m2: String = m2_chars.into_iter().collect();

    // Confirm mate2 truly can't seed on its own.
    let fq2_alone = write_fastq(&dir, "m2_alone.fq", &[("pair1".to_string(), m2.clone())]);
    let alone = map_reads(&prefix, &prefix, &[&fq2_alone], MapperOptions::default()).unwrap();
    assert_eq!(alone.len(), 1);
    assert_eq!(alone[0].position, 0, "test setup: mate2 should NOT map on its own without rescue");

    // Native paired mapping keeps mate1 and rejects mate2.
    let fq1 = write_fastq(&dir, "m1.fq", &[("pair1".to_string(), m1)]);
    let fq2 = write_fastq(&dir, "m2.fq", &[("pair1".to_string(), m2)]);
    let records = map_reads(&prefix, &prefix, &[&fq1, &fq2], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 2);
    let (r1, r2) = (&records[0], &records[1]);
    assert_eq!(r1.position, 41);
    // Stampy 1.0.32 rejects this 40 bp / seven-substitution mate on the
    // same synthetic reference (native FLAG 133, CIGAR *, MAPQ 0).
    // The old custom scoring recovered it by assigning no deletion cost.
    assert_eq!(r2.position, 0);
    assert_eq!(r2.mapq, 0);
    assert!(r2.cigar.is_empty());

    // Its SAM coordinate is copied from the mapped mate while FLAG 4 remains
    // set. Such records must sort at their emitted position for BAM indexing.
    let sam = io::format_sam(r2, None);
    let fields: Vec<_> = sam.split('\t').collect();
    assert_eq!(fields[2], "chr1");
    assert_eq!(fields[3], "41");
    assert_eq!(fields[1].parse::<u16>().unwrap() & 4, 4);
    let mut late = r1.clone(); late.position=301; late.read_id="late".into();
    let mut sorted=vec![late,r1.clone(),r2.clone()];
    let contigs=vec![("chr1".into(),400)];
    io::sort_for_coordinate_order(&mut sorted,&contigs);
    assert_eq!(sorted[2].read_id,"late");
    let sam=io::format_sam_records_sorted(&sorted,&contigs,"test",None,true);
    io::write_bam_indexed(&sam,&dir.path().join("copied-coordinate.bam")).unwrap();
}

#[test]
fn test_n_run_round_trip() {
    let dir = TempDir::new().unwrap();
    let mut rng = Xorshift64(7);
    let left = rng.gen_seq(100);
    let right = rng.gen_seq(100);
    let genome_seq = format!("{}{}{}", left, "N".repeat(20), right);
    let fasta = write_fasta(&dir, "refN.fa", "chrN", &genome_seq);
    let prefix = build_index(&dir, &fasta, "refN");

    let in_gap = genome_seq[105..135].to_string();
    let right_unique = genome_seq[150..190].to_string();
    let fastq = write_fastq(&dir, "reads.fq", &[
        ("in_gap".to_string(), in_gap),
        ("right_unique".to_string(), right_unique),
    ]);

    let records = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 2);

    let in_gap_rec = &records[0];
    assert_eq!(in_gap_rec.position, 106);
    assert_eq!(&in_gap_rec.read_sequence[0..15], b"NNNNNNNNNNNNNNN");

    let right_rec = &records[1];
    assert_eq!(right_rec.position, 151);
    // See the comment on the equivalent assertion in
    // test_forward_and_reverse_strand_mapping for why a 40bp read caps
    // at 87 rather than 99.
    assert_eq!(right_rec.mapq, 87);
}

/// Larger-scale synthetic accuracy check: a 20kb random genome and 300
/// simulated reads (some mismatched, some reverse-complemented), asserting a
/// concrete mapping-accuracy threshold and reporting elapsed time. This is
/// still synthetic data, not real sequencing data or a real genome, but it's
/// a meaningfully larger and more systematic check than the small
/// hand-constructed cases above.
#[test]
fn test_synthetic_scale_mapping_accuracy() {
    let dir = TempDir::new().unwrap();
    let mut rng = Xorshift64(12345);
    let genome_len = 20_000usize;
    let genome_seq = rng.gen_seq(genome_len);
    let fasta = write_fasta(&dir, "big.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "big");

    let read_len = 60usize;
    let n_reads = 300usize;
    let mut inputs = Vec::with_capacity(n_reads);
    let mut truth = Vec::with_capacity(n_reads);

    for i in 0..n_reads {
        let start = (rng.next_u64() as usize) % (genome_len - read_len);
        let mut bases: Vec<char> = genome_seq[start..start + read_len].chars().collect();
        if rng.next_u64() % 10 < 3 {
            let pos = (rng.next_u64() as usize) % read_len;
            let orig = bases[pos];
            bases[pos] = ['A', 'C', 'G', 'T'].into_iter().find(|&b| b != orig).unwrap();
        }
        let mut seq: String = bases.into_iter().collect();
        let reverse = rng.next_u64().is_multiple_of(2);
        if reverse {
            seq = revcomp(&seq);
        }
        inputs.push((format!("r{i}"), seq));
        truth.push((start + 1, if reverse { Strand::Reverse } else { Strand::Forward }));
    }

    let fastq = write_fastq(&dir, "many.fq", &inputs);

    let started = std::time::Instant::now();
    let records = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
    let elapsed = started.elapsed();

    let mut correct = 0;
    for (rec, (expected_pos, expected_strand)) in records.iter().zip(truth.iter()) {
        if rec.contig_name == "chr1"
            && rec.strand == *expected_strand
            && (rec.position as i64 - *expected_pos as i64).abs() <= 2
        {
            correct += 1;
        }
    }
    let accuracy = correct as f64 / n_reads as f64;
    eprintln!(
        "synthetic mapping accuracy: {correct}/{n_reads} ({:.1}%) in {elapsed:?}",
        accuracy * 100.0
    );
    assert!(accuracy >= 0.95, "mapping accuracy too low: {:.1}%", accuracy * 100.0);
}

/// Regression/validation test for the paired-end candidate-shortlist +
/// joint-posterior procedure (`mapper::paired`): a read whose mate anchors
/// it uniquely should get correctly redirected to
/// the true copy of a repeat, even when its *own* single-end best pick lands
/// on a different (better-matching, but biologically wrong) copy of that
/// repeat far away.
///
/// Genome layout on one contig:
///   [0, 100)      unique region (mate1's source)
///   [100, 210)    unique spacer
///   [210, 310)    "true" repeat copy (one base different from the repeat
///                 unit mate2 is drawn from, at a low-quality position)
///   [310, 10310)  unique filler (10kb -- puts the decoy copy far outside
///                 any plausible insert-size window)
///   [10310, 10410) "decoy" repeat copy (an exact copy of the repeat unit)
///
/// mate2 is drawn to exactly match the decoy copy and differ from the true
/// copy by one low-quality base, so standalone single-end mapping prefers
/// the (wrong) decoy: a perfect match beats a match with one low-quality
/// mismatch, with no insert-size information to say otherwise. Paired
/// mapping has that information -- the true copy sits ~210bp from mate1
/// (within the default insert-size model), while the decoy sits ~10060bp
/// away (many SDs out, floored to the structural-variant prior) -- so the
/// joint posterior should favor the true copy despite its slightly worse
/// standalone likelihood.
#[test]
fn test_paired_disambiguates_repeat_copy_via_anchor() {
    use probebwa::mapper::{PairedEndMapper, SingleEndMapper};

    let mut rng = Xorshift64(2026);
    let unique_a = rng.gen_seq(100);
    let spacer = rng.gen_seq(110);
    let repeat_unit = rng.gen_seq(100);
    let filler = rng.gen_seq(10_000);

    // The true copy differs from the repeat unit at one position (index 20,
    // comfortably inside the first 40bp mate2 is drawn from).
    const DIFF_POS: usize = 20;
    let mut true_copy: Vec<u8> = repeat_unit.bytes().collect();
    true_copy[DIFF_POS] = match true_copy[DIFF_POS] {
        b'A' => b'C', b'C' => b'A', b'G' => b'T', _ => b'G',
    };
    let true_copy = String::from_utf8(true_copy).unwrap();
    let decoy_copy = repeat_unit.clone();

    let genome_seq = format!("{unique_a}{spacer}{true_copy}{filler}{decoy_copy}");
    let true_copy_pos = unique_a.len() + spacer.len(); // 210
    let decoy_copy_pos = genome_seq.len() - decoy_copy.len(); // 10310

    let dir = TempDir::new().unwrap();
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = build_index(&dir, &fasta, "ref");
    let genome = index::GenomeIndex::load(&prefix).unwrap();
    let hash_table = index::HashTable::load(&prefix).unwrap();
    let options = MapperOptions::default();

    // mate1: a clean, unique 40bp read from the unique region -- the anchor.
    let read1 = Read {
        id: "pair1".to_string(),
        sequence: DnaSeq::from_ascii(&unique_a.as_bytes()[0..40]),
        qualities: QualityScores { scores: vec![35u8; 40] },
        is_reverse: false,
    };

    // mate2: reverse-complement of the repeat unit's first 40bp (so it maps
    // reverse-strand, matching FR orientation with mate1's forward mapping),
    // with a single low-quality base at the position that differs between
    // the true and decoy copies -- low quality keeps that one mismatch's
    // cost small relative to the insert-size evidence, but still large
    // enough that a perfect match (the decoy) wins on single-end likelihood
    // alone.
    let mate2_fwd_window = &repeat_unit[0..40];
    let mate2_seq = revcomp(mate2_fwd_window);
    let mut mate2_qual = vec![35u8; 40];
    mate2_qual[39 - DIFF_POS] = 10; // index flips under reverse-complement
    let read2 = Read {
        id: "pair1".to_string(),
        sequence: DnaSeq::from_ascii(mate2_seq.as_bytes()),
        qualities: QualityScores { scores: mate2_qual },
        is_reverse: false,
    };

    // Sanity check: mapped alone (no pairing information), mate2 should
    // indeed prefer the decoy copy -- confirming this scenario actually
    // tests disambiguation, not something seeding/alignment would have
    // gotten right anyway.
    let single_mapper = SingleEndMapper::new(&genome, &hash_table, &options);
    let standalone = single_mapper.map(&read2).unwrap();
    assert_eq!(
        standalone.position, decoy_copy_pos + 1,
        "sanity check failed: standalone mate2 should map to the decoy copy \
         (position {}), got {} -- scenario doesn't actually exercise \
         disambiguation", decoy_copy_pos + 1, standalone.position
    );

    // Paired mapping should instead find the true copy, using mate1 as the
    // anchor.
    let mut paired_mapper = PairedEndMapper::new(&genome, &hash_table, &options);
    let pair = ReadPair { read1, read2 };
    let (r1, r2) = paired_mapper.map_pair(&pair).unwrap();

    assert_eq!(r1.position, 1, "mate1 should map at the start of the unique region");
    assert_eq!(
        r2.position, true_copy_pos + 1,
        "mate2 should be redirected to the true copy (position {}) via its \
         anchor, not left at its own standalone best pick (the decoy, {})",
        true_copy_pos + 1, decoy_copy_pos + 1
    );
    assert_eq!(r2.strand, Strand::Reverse);
    assert!(r1.is_proper_pair && r2.is_proper_pair, "should be flagged as a proper pair once correctly placed");
    assert!(r2.mapq > 0, "disambiguated mate should still pass the entropy gate");
}

#[test]
fn test_paired_names_must_match_after_normalization() {
    let dir = TempDir::new().unwrap();
    let genome = Xorshift64(123).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome);
    let prefix = build_index(&dir, &fasta, "ref");
    let first = write_fastq(&dir, "r1.fq", &[("sampleA/1".into(), genome[40..80].into())]);
    let second = write_fastq(&dir, "r2.fq", &[("sampleB/2".into(), revcomp(&genome[200..240]))]);
    let error = map_reads(&prefix, &prefix, &[&first, &second], MapperOptions::default()).unwrap_err();
    assert!(error.to_string().contains("pair 1: 'sampleA' vs 'sampleB'"));
    let second = write_fastq(&dir, "r2.fq", &[("sampleA/2".into(), revcomp(&genome[200..240]))]);
    let records = map_reads(&prefix, &prefix, &[&first, &second], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 2);
    assert!(records.iter().all(|r| r.read_id == "sampleA"));
}

#[test]
fn test_cross_contig_pairs_preserve_mate_coordinates() {
    let dir = TempDir::new().unwrap();
    let mut rng = Xorshift64(137);
    let chr1 = rng.gen_seq(2000);
    let chr2 = rng.gen_seq(2000);
    let fasta = dir.path().join("ref.fa");
    fs::write(&fasta, format!(">chr1\n{chr1}\n>chr2\n{chr2}\n")).unwrap();
    let prefix = build_index(&dir, &fasta, "ref");
    let fq1 = write_fastq(&dir, "r1.fq", &[("cross/1".into(), chr1[100..220].into())]);
    let fq2 = write_fastq(&dir, "r2.fq", &[("cross/2".into(), revcomp(&chr2[300..420]))]);
    let records = map_reads(&prefix, &prefix, &[&fq1, &fq2], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 2);
    for (record, contig, pos, mate_contig, mate_pos) in [
        (&records[0], "chr1", 101, "chr2", 301),
        (&records[1], "chr2", 301, "chr1", 101),
    ] {
        assert_eq!(record.contig_name, contig);
        assert_eq!(record.position, pos);
        assert_eq!(record.mate_contig.as_deref(), Some(mate_contig));
        assert_eq!(record.mate_position, Some(mate_pos));
        assert!(!record.mate_unmapped);
        assert!(!record.is_proper_pair);
        let sam = io::format_sam(record, None);
        let fields: Vec<_> = sam.split('\t').collect();
        assert_eq!(fields[6], mate_contig);
        assert_eq!(fields[7], mate_pos.to_string());
        assert_eq!(fields[8], "0");
    }
}

#[test]
fn test_tlen_signs_for_equal_starts_and_swapped_mates() {
    let dir = TempDir::new().unwrap();
    let genome = Xorshift64(321).gen_seq(2000);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome);
    let prefix = build_index(&dir, &fasta, "ref");
    for (start1, start2, length) in [(100, 100, 120), (100, 300, 320)] {
        let first = write_fastq(&dir, "r1.fq", &[("pair".into(), genome[start1..start1 + 120].into())]);
        let second = write_fastq(&dir, "r2.fq", &[("pair".into(), revcomp(&genome[start2..start2 + 120]))]);
        for (files, pos1, pos2) in [([&first, &second], start1, start2), ([&second, &first], start2, start1)] {
            let records = map_reads(&prefix, &prefix, &files, MapperOptions::default()).unwrap();
            assert_eq!(records[0].position, pos1 + 1);
            assert_eq!(records[1].position, pos2 + 1);
            let signed = if records[0].strand == Strand::Forward { length } else { -length };
            assert_eq!(records[0].insert_size, Some(signed));
            assert_eq!(records[1].insert_size, Some(-signed));
        }
    }
}

#[test]
fn test_mapping_validates_options_but_accepts_empty_files_and_rate_boundaries() {
    let dir = TempDir::new().unwrap();
    let missing = dir.path().join("missing");
    let error = map_reads(&missing, &missing, &[], MapperOptions::default()).unwrap_err();
    assert!(error.to_string().contains("at least one read file"));
    for rate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 1.1] {
        let options = MapperOptions { substitution_rate: rate, ..Default::default() };
        let error = map_reads(&missing, &missing, &[&missing], options).unwrap_err();
        assert!(error.to_string().contains("substitution rate must be finite"));
    }
    let genome = Xorshift64(322).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome);
    let prefix = build_index(&dir, &fasta, "ref");
    let empty = write_fastq(&dir, "empty.fq", &[]);
    let reads = write_fastq(&dir, "reads.fq", &[("read".into(), genome[40..160].into())]);
    for rate in [0.0, 1.0] {
        let options = MapperOptions { substitution_rate: rate, ..Default::default() };
        assert!(map_reads(&prefix, &prefix, &[&empty], options.clone()).unwrap().is_empty());
        assert_eq!(map_reads(&prefix, &prefix, &[&reads], options).unwrap().len(), 1);
    }
}

#[test]
fn native_single_and_paired_use_different_refinement_qualities() {
    let dir=TempDir::new().unwrap();
    let genome=Xorshift64(88291).gen_seq(1000);
    let fasta=write_fasta(&dir,"quality.fa","chr1",&genome);
    let prefix=build_index(&dir,&fasta,"quality");
    let mut read=genome.as_bytes()[100..160].to_vec();
    read[10]=if read[10]==b'A' { b'C' } else { b'A' };
    let first=write_fastq(&dir,"quality1.fq",&[("quality".into(),String::from_utf8(read).unwrap())]);
    let second=write_fastq(&dir,"quality2.fq",&[("quality".into(),revcomp(&genome[300..360]))]);
    let single=map_reads(&prefix,&prefix,&[&first],MapperOptions::default()).unwrap();
    let paired=map_reads(&prefix,&prefix,&[&first,&second],MapperOptions::default()).unwrap();
    assert_eq!(single[0].position,101);
    assert_eq!(paired[0].position,101);
    assert!(single[0].sam_tags.contains(&"UQ:i:30".into()));
    assert!(paired[0].sam_tags.contains(&"UQ:i:40".into()));
    assert!(paired[1].sam_tags.contains(&"XQ:i:40".into()));
}

#[test]
fn nearby_seed_anchors_do_not_chain_across_distinct_loci() {
    let hits: Vec<_>=(100..104).map(|pos| SeedHit { read_pos:0,ref_contig:0,ref_pos:pos,k:15,hit_strand:Strand::Forward }).collect();
    let candidates=mapper::cluster_hits(&hits,40);
    let mut positions:Vec<_>=candidates.iter().map(|c|c.position).collect(); positions.sort();
    assert_eq!(positions,vec![100,102]);
}

#[test]
fn native_slow_pair_posterior_floor() {
    let dir=TempDir::new().unwrap();
    let genome=Xorshift64(99721).gen_seq(1000);
    let fasta=write_fasta(&dir,"posterior.fa","chr1",&genome);
    let prefix=build_index(&dir,&fasta,"posterior");
    let first=write_fastq(&dir,"posterior1.fq",&[("short-pair".into(),genome[100..200].into())]);
    let second=write_fastq(&dir,"posterior2.fq",&[("short-pair".into(),revcomp(&genome[110..210]))]);
    let records=map_reads(&prefix,&prefix,&[&first,&second],MapperOptions { insert_size_sd:30.0, ..MapperOptions::default() }).unwrap();
    assert_eq!(records[0].position,101); assert_eq!(records[1].position,111);
    // Native realign_pairs clamps its normalization sum to 1.0000000001.
    // The Phred-99 error floor first gives 96. On this tiny reference the
    // distance-prior correction then combines two Phred-96 terms into 93.
    assert_eq!(records[0].mapq,93); assert_eq!(records[1].mapq,93);
}

#[test]
fn native_contig_edge_retry_keeps_terminal_insertion() {
    // Native Stampy independently verified on this reduced reference. The
    // read at the contig edge was previously rejected with an empty CIGAR.
    let fields: Vec<_> = include_str!("fixtures/stampy_edge_retry.tsv")
        .lines().find(|line| !line.starts_with('#')).unwrap().split('\t').collect();
    let dir = TempDir::new().unwrap();
    let fasta = write_fasta(&dir, "edge.fa", "edge", fields[0]);
    let prefix = build_index(&dir, &fasta, "edge");
    let first = dir.path().join("edge1.fq");
    let second = dir.path().join("edge2.fq");
    fs::write(&first, format!("@edge/1\n{}\n+\n{}\n", fields[1], fields[2])).unwrap();
    fs::write(&second, format!("@edge/2\n{}\n+\n{}\n", fields[3], fields[4])).unwrap();
    let records = map_reads(&prefix, &prefix, &[&first, &second], MapperOptions {
        insert_size_mean: 250.0, insert_size_sd: 60.0, ..MapperOptions::default()
    }).unwrap();
    assert_eq!(records[0].position, 902);
    assert_eq!(records[0].strand, Strand::Reverse);
    assert_eq!(records[0].cigar, vec![CigarOp::Match(99), CigarOp::Insertion(1)]);
    assert_eq!(records[1].position, 561);
    assert_eq!(records[1].strand, Strand::Forward);
    assert_eq!(records[1].cigar, vec![CigarOp::Match(100)]);
}

#[test]
fn reverse_seed_clusters_follow_signed_locus_order() {
    // Native signed loci run in the opposite order to forward coordinates.
    // Three adjacent anchors expose the grouping difference hidden by an
    // even-length run: (102,101), (100), rather than (100,101), (102).
    let hits: Vec<_> = (100..103).map(|pos| SeedHit {
        read_pos: 0, ref_contig: 0, ref_pos: pos, k: 15, hit_strand: Strand::Reverse,
    }).collect();
    let candidates = mapper::cluster_hits(&hits, 40);
    let mut positions: Vec<_> = candidates.iter().map(|c| c.position).collect();
    positions.sort();
    assert_eq!(positions, vec![75, 77]);
}

#[test]
fn equal_likelihood_preserves_first_candidate() {
    let calculator = mapq::MapqCalculator::new(0.001, 15, 1000);
    let (index, sum) = calculator.select_best(&[0.0, -1.0, 0.0]);
    assert_eq!(index, 0);
    assert!((sum - (2.0 + (-1.0_f64).exp())).abs() < 1e-12);
}

#[test]
fn native_pair_tie_selects_reverse_repeat_copy() {
    let dir = TempDir::new().unwrap();
    let sequence = Xorshift64(891273).gen_seq(700);
    let block = &sequence[100..500];
    let reference = format!("{}{}{}{}{}", &sequence[..100], block, &sequence[500..600], revcomp(block), &sequence[600..]);
    let fasta = write_fasta(&dir, "repeat.fa", "repeat", &reference);
    let prefix = build_index(&dir, &fasta, "repeat");
    let first = write_fastq(&dir, "repeat1.fq", &[("repeat/1".into(), block[..100].into())]);
    let second = write_fastq(&dir, "repeat2.fq", &[("repeat/2".into(), revcomp(&block[300..]))]);
    let records = map_reads(&prefix, &prefix, &[&first, &second], MapperOptions {
        insert_size_mean: 250.0, insert_size_sd: 60.0, ..MapperOptions::default()
    }).unwrap();
    // Native Stampy 1.0.32: read 1 at 901 reverse, read 2 at 601 forward.
    // Both repeat copies have identical pair scores and alignment strings.
    assert_eq!((records[0].position, records[0].strand), (901, Strand::Reverse));
    assert_eq!((records[1].position, records[1].strand), (601, Strand::Forward));
    assert_eq!(records[0].read_sequence, revcomp(&block[..100]).as_bytes());
    assert_eq!(records[1].read_sequence, revcomp(&block[300..]).as_bytes());
    assert!(records.iter().all(|r| r.cigar == vec![CigarOp::Match(100)]));
}

#[test]
fn native_same_signed_locus_falls_back_to_single_confidence() {
    let dir = TempDir::new().unwrap();
    let genome = Xorshift64(99721).gen_seq(1000);
    let fasta = write_fasta(&dir, "twins.fa", "chr1", &genome);
    let prefix = build_index(&dir, &fasta, "twins");
    for reverse in [false, true] {
        // Reverse reads have different left edges but the same sequenced
        // 5' endpoint. Native drops mate 1 on equal fine alignment scores.
        let (s1, s2, flags) = if reverse {
            (revcomp(&genome[120..200]), revcomp(&genome[100..200]), [101, 153])
        } else { (genome[100..200].into(), genome[100..200].into(), [69, 137]) };
        let first = write_fastq(&dir, "twin1.fq", &[("twin/1".into(), s1)]);
        let second = write_fastq(&dir, "twin2.fq", &[("twin/2".into(), s2)]);
        let records = map_reads(&prefix, &prefix, &[&first, &second], MapperOptions {
            insert_size_mean: 250.0, insert_size_sd: 60.0, ..MapperOptions::default()
        }).unwrap();
        assert_eq!(records[0].position, 0);
        assert_eq!(records[1].position, 101);
        assert_eq!(records[1].mapq, 96);
        assert!(records[0].sam_tags.contains(&"MQ:i:96".into()));
        for (record, flag) in records.iter().zip(flags) {
            let sam = io::format_sam(record, None);
            let fields: Vec<_> = sam.split('\t').collect();
            assert_eq!(fields[1], flag.to_string());
            assert_eq!(fields[3], "101");
            assert!(!record.sam_tags.iter().any(|tag| tag.starts_with("PQ:")));
        }
    }
}

#[test]
fn native_refinement_shift_invalidates_original_single_locus() {
    let fields: Vec<_> = include_str!("fixtures/stampy_refined_anchor.tsv")
        .lines().find(|line| !line.starts_with('#')).unwrap().split('\t').collect();
    let dir = TempDir::new().unwrap();
    let fasta = write_fasta(&dir, "anchor.fa", "anchor", fields[0]);
    let prefix = build_index(&dir, &fasta, "anchor");
    let first = dir.path().join("anchor1.fq");
    let second = dir.path().join("anchor2.fq");
    fs::write(&first, format!("@anchor/1\n{}\n+\n{}\n", fields[1], fields[2])).unwrap();
    fs::write(&second, format!("@anchor/2\n{}\n+\n{}\n", fields[3], fields[4])).unwrap();
    let records = map_reads(&prefix, &prefix, &[&first, &second], MapperOptions {
        insert_size_mean: 425.0, insert_size_sd: 40.0, ..MapperOptions::default()
    }).unwrap();
    assert_eq!(records[0].position, 607);
    assert_eq!(records[1].position, 424);
    assert!(records[0].sam_tags.contains(&"SM:i:96".into()));
    assert!(records[0].sam_tags.contains(&"MQ:i:0".into()));
    assert!(records[1].sam_tags.contains(&"SM:i:0".into()));
    assert_eq!(records[1].cigar, vec![CigarOp::Insertion(8), CigarOp::Match(92)]);
}
