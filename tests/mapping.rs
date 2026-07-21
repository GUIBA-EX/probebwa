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
    // Not 99: a 40bp read only offers 2 non-overlapping 15-mer candidate
    // windows per read/genome phase offset (see
    // `mapq::bayesian::MapqCalculator::missing_locus_probability`), so the
    // default 0.1% substitution-rate prior alone caps confidence around 80
    // even for a perfect, uniquely-mapping match — this is the corrected
    // (paper-faithful) missing-locus model, not a mapping defect.
    assert_eq!(fwd.mapq, 80);

    let rev = &records[1];
    assert_eq!(rev.contig_name, "chr1");
    assert_eq!(rev.position, 201);
    assert_eq!(rev.strand, Strand::Reverse);
    assert_eq!(rev.cigar, vec![CigarOp::Match(40)]);
    assert_eq!(rev.mapq, 80);
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
    // around 80 rather than 99.
    assert_eq!(mismatch.mapq, 80);

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

/// Regression test: a large deletion (25bp, well within the default
/// `max_indel_len` of 30) used to be silently unfindable — the DP would
/// prefer an alternative alignment that avoided the gap's steep cost
/// entirely, producing a CIGAR with 0-4 net deleted bases instead of 25 (see
/// the gap-cost recalibration this test locks in: gap_extend went from -2 to
/// -1). The recovered alignment can still legitimately fragment the deletion
/// into several smaller pieces (ties broken by incidental sequence
/// composition), so this checks the aggregate deleted-base total rather than
/// the exact op sequence.
#[test]
fn test_large_deletion_is_found() {
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
    let total_deleted: u32 = rec.cigar.iter().map(|op| match op {
        CigarOp::Deletion(n) => *n,
        _ => 0,
    }).sum();
    assert_eq!(total_deleted, deletion_len as u32, "expected a net {deletion_len}bp deletion, got cigar {:?}", rec.cigar);
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
    assert_eq!(records[0].mapq, 99);
}

/// Regression test for mate rescue: mate2 is diverged heavily enough (a
/// mismatch every 6bp) that *every* overlapping 15-mer window contains at
/// least 2 mismatches — defeating the seeding pass's single-mismatch
/// tolerance everywhere, so mate2 produces zero candidates on its own
/// (confirmed below by mapping it alone). Mapped as a pair, mate1 anchors a
/// direct-alignment search window for mate2 via the insert-size estimate,
/// and should recover it anyway.
#[test]
fn test_mate_rescue_recovers_unseedable_mate() {
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

    // Now map as a pair: mate1 anchors, mate2 should be rescued.
    let fq1 = write_fastq(&dir, "m1.fq", &[("pair1".to_string(), m1)]);
    let fq2 = write_fastq(&dir, "m2.fq", &[("pair1".to_string(), m2)]);
    let records = map_reads(&prefix, &prefix, &[&fq1, &fq2], MapperOptions::default()).unwrap();
    assert_eq!(records.len(), 2);
    let (r1, r2) = (&records[0], &records[1]);
    assert_eq!(r1.position, 41);
    assert_ne!(r2.position, 0, "mate2 should have been rescued using mate1's anchor");
    assert_eq!(r2.contig_name, "chr1");
    assert_eq!(r2.strand, Strand::Reverse);
    assert!((r2.position as i64 - 201).abs() <= 5, "rescued position should land near the true locus (201), got {}", r2.position);
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
    // around 80 rather than 99.
    assert_eq!(right_rec.mapq, 80);
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
