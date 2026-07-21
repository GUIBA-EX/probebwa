//! Integration tests.

use probebwa::*;
use probebwa::align::SmithWatermanAligner;

#[test]
fn test_dna_seq_kmers() {
    let seq = DnaSeq::from_ascii(b"ACGTACGT");
    let kmers = seq.kmers(4);
    assert_eq!(kmers.len(), 5);
    assert_eq!(kmers[0], b"ACGT");
}

#[test]
fn test_reverse_complement() {
    let seq = DnaSeq::from_ascii(b"ACGT");
    let rc = seq.reverse_complement();
    assert_eq!(rc.bases, b"ACGT"); // palindrome
    let seq2 = DnaSeq::from_ascii(b"AACCGGTT");
    assert_eq!(seq2.reverse_complement().bases, b"AACCGGTT");
}

#[test]
fn test_quality_scores() {
    let qs = QualityScores::from_phred_ascii(b"!\"#", 33);
    assert_eq!(qs.scores, vec![0, 1, 2]);
    assert!((qs.error_probability(0) - 1.0).abs() < 1e-10);
}

#[test]
fn test_hash_kmer_15mer() {
    // hash_kmer runs the packed k-mer through a pseudo-randomizing finalizer,
    // so we check its contract rather than a raw 2-bit value: deterministic,
    // and distinct for distinct k-mers.
    let h1 = index::HashTable::hash_kmer(b"ACGTACGTACGTACG");
    let h2 = index::HashTable::hash_kmer(b"ACGTACGTACGTACG");
    assert_eq!(h1, h2);
    let h3 = index::HashTable::hash_kmer(b"TTTTTTTTTTTTTTT");
    assert_ne!(h1, h3);
}

#[test]
fn test_pack_unpack_position() {
    let packed = index::HashTable::pack_position(5, 1234, Strand::Forward);
    let (c, p, s) = index::HashTable::unpack_position(packed);
    assert_eq!(c, 5);
    assert_eq!(p, 1234);
    assert_eq!(s, Strand::Forward);

    let packed_rev = index::HashTable::pack_position(5, 1234, Strand::Reverse);
    let (c2, p2, s2) = index::HashTable::unpack_position(packed_rev);
    assert_eq!(c2, 5);
    assert_eq!(p2, 1234);
    assert_eq!(s2, Strand::Reverse);

    // Positions well beyond 16 bits must round-trip exactly — this is the
    // whole point of u64 packing (a chromosome-scale coordinate).
    let big_pos = 200_000_000usize; // ~human chr1 scale
    let packed_big = index::HashTable::pack_position(42, big_pos, Strand::Forward);
    let (c3, p3, s3) = index::HashTable::unpack_position(packed_big);
    assert_eq!(c3, 42);
    assert_eq!(p3, big_pos);
    assert_eq!(s3, Strand::Forward);
}

/// Regression test: a record with a zero-length sequence used to desync the
/// FASTQ parser — its (blank) quality line was never consumed, so the next
/// record's `@id` line got misread as that quality line, and everything
/// after that was misparsed too. Here a malformed empty-sequence record is
/// followed by a normal one; both must parse cleanly.
#[test]
fn test_fastq_empty_sequence_record_does_not_desync_parser() {
    let dir = tempfile::TempDir::new().unwrap();
    let path = dir.path().join("empty_seq.fq");
    std::fs::write(&path, "@empty\n\n+\n\n@normal\nACGT\n+\nIIII\n").unwrap();

    let reads = io::read_fastq(&path, io::PHRED33_OFFSET).unwrap();
    assert_eq!(reads.len(), 2);
    assert_eq!(reads[0].id, "empty");
    assert_eq!(reads[0].sequence.bases, b"");
    assert_eq!(reads[1].id, "normal");
    assert_eq!(reads[1].sequence.bases, b"ACGT");
    assert_eq!(reads[1].qualities.scores, vec![40, 40, 40, 40]);
}

/// Regression test: both FASTA parsers (reads via `--inputformat=fasta`, and
/// genome-building) used to discard the header line that terminates a
/// record's sequence scan instead of carrying it forward to the next call,
/// silently dropping every other record in a multi-record file. Found by
/// comparing against another tool on the same multi-record FASTA input,
/// where probebwa quietly processed only 11 of 22 reads.
#[test]
fn test_fasta_multi_record_does_not_drop_every_other_entry() {
    let dir = tempfile::TempDir::new().unwrap();

    // FASTA-as-reads path.
    let reads_path = dir.path().join("multi.fa");
    std::fs::write(&reads_path, ">r1\nACGT\n>r2\nTTTT\n>r3\nGGGG\n>r4\nCCCC\n").unwrap();
    let reads = io::read_fasta_as_reads(&reads_path).unwrap();
    let ids: Vec<&str> = reads.iter().map(|r| r.id.as_str()).collect();
    assert_eq!(ids, vec!["r1", "r2", "r3", "r4"]);

    // Genome-building path.
    let genome_path = dir.path().join("multi_genome.fa");
    std::fs::write(&genome_path, ">c1\nAAAA\n>c2\nCCCC\n>c3\nGGGG\n>c4\nTTTT\n").unwrap();
    let mut builder = index::GenomeIndexBuilder::new("test", "t1");
    builder.add_fasta(&genome_path).unwrap();
    let prefix = dir.path().join("genome");
    builder.build_and_save(&prefix).unwrap();
    let loaded = index::GenomeIndex::load(&prefix).unwrap();
    let contig_names: Vec<&str> = loaded.contigs.iter().map(|c| c.name.as_str()).collect();
    assert_eq!(contig_names, vec!["c1", "c2", "c3", "c4"]);
}

/// Regression test found by benchmarking against real data (E. coli K-12
/// MG1655, NCBI NC_000913.3, vs. real paired-end reads from ENA run
/// DRR002055): a wide alignment window (as `mapper::single::align_at_window`
/// uses for mate rescue / paired shortlist cross-alignment) containing a
/// locally repeat-like region used to let the DP "trade" several real
/// mismatches for a handful of spurious short indels that happen to land on
/// a cleaner-looking but biologically wrong local register — because, under
/// the old cost balance, a few gap-opens were cheap enough to be worth it for
/// the resulting run of extra matches, even though nothing here is actually
/// insertion/deletion variation.
///
/// This exact 900bp window and read are real E. coli sequence/data (not
/// synthetic) — the true placement is position 3607868 (1-based) with
/// CIGAR `100M` (13 mismatches, no indels); this test's window is centered
/// the same way `align_at_window` would build it for that real read.
#[test]
fn test_divergent_region_does_not_invent_spurious_indels() {
    const WINDOW: &[u8] = b"TCCCTGCTGGTGCCACCAGCGTAGACCGTCTGGTGACGTTGGAAGTGCTGTCAGAACCGGGAGCCAGCGCCATTGACCGGATTCTGAAACTGATCGAAGAAGCCGAAGAGCGTCGCGCTCCCATTGAGCGGTTTATCGACCGTTTCAGCCGTATCTATACGCCCGCGATTATGGCCGTCGCTCTGCTGGTGACGCTGGTGCCACCGCTGCTGTTTGCCGCCAGCTGGCAGGAGTGGATTTATAAAGGGCTGACGCTGCTGCTGATTGGCTGCCCGTGTGCGTTAGTTATCTCAACGCCTGCGGCGATTACCTCCGGGCTGGCGGCGGCAGCGCGTCGTGGGGCGTTGATTAAAGGCGGAGCGGCGCTGGAACAGCTGGGTCGTGTTACTCAGGTGGCGTTTGATAAAACCGGTACGCTGACCGTCGGTAAACCGCGCGTTACCGCGATTCATCCGGCAACGGGTATTAGTGAATCTGAACTGCTGACACTGGCGGCGGCGGTCGAGCAAGGCGCGACGCATCCACTGGCGCAAGCCATCGTACGCGAAGCACAGGTTGCTGAACTCGCCATTCCCACCGCCGAATCACAGCGGGCGCTGGTCGGGTCTGGCATTGAAGCGCAGGTTAACGGTGAGCGCGTATTGATTTGCGCTGCCGGGAAACATCCCGCTGATGCATTTACTGGTTTAATTAACGAACTGGAAAGCGCCGGGCAAACGGTAGTGCTGGTAGTACGTAACGATGACGTGCTTGGTGTCATTGCGTTACAGGATACCCTGCGCGCCGATGCTGCAACTGCCATCAGTGAACTGAACGCGCTGGGCGTCAAAGGGGTGATCCTCACCGGCGATAATCCACGCGCAGCGGCGGCAATTGCCGGGGAGCTGGGGCTGGAGTTTA";
    const READ: &[u8] = b"CGGCGCGGCCCCTGGCGGAGGCGGTCGCCCGCGAAGCACAGGTTGCTCACCTCGCCCTTCCCACCGCCGAATCACAGCGGGCGCTGGTCGGGTCTGGCAT";
    // Real Illumina quality (mostly Q30, "?" in Phred+33) round-tripped from
    // the actual read.
    let qual = vec![30u8; READ.len()];

    // Mirrors how `align_at_window` builds its aligner: banded to the full
    // window width, diagonal centered where the read is expected to start
    // (513 here — the true, correct offset).
    let aligner = SmithWatermanAligner::new(WINDOW.len(), 40, 3);
    let (_score, cigar, ref_start_offset) = aligner.align(READ, &qual, WINDOW, 513);

    let indel_bases: u32 = cigar.iter().map(|op| match op {
        CigarOp::Insertion(n) | CigarOp::Deletion(n) => *n,
        _ => 0,
    }).sum();
    let indel_ops = cigar.iter().filter(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_))).count();

    assert_eq!(ref_start_offset, 513, "should land on the true, collinear offset, not a shifted register");
    assert_eq!(indel_ops, 0, "a genuinely indel-free (mismatch-only) region shouldn't get spurious indels invented; got CIGAR {cigar:?}");
    assert_eq!(indel_bases, 0);
}
