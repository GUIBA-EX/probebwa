//! CLI-level tests for flags that live in `main.rs` rather than the library
//! (`--outputformat=bam`, `--labelfilter`), so they run the compiled binary
//! directly instead of calling library functions.

use probebwa::*;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
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

struct Xorshift64(u64);
impl Xorshift64 {
    fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn gen_seq(&mut self, len: usize) -> String {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        (0..len).map(|_| BASES[(self.next_u64() % 4) as usize]).collect()
    }
}

fn probebwa_bin() -> &'static str {
    env!("CARGO_BIN_EXE_probebwa")
}

#[test]
fn test_bam_output_is_valid_bgzf_bam() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(99).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = dir.path().join("ref");
    build_genome_index("t", "t1", prefix.to_str().unwrap(), &[fasta.to_str().unwrap()]).unwrap();
    build_hash_table(prefix.to_str().unwrap(), prefix.to_str().unwrap()).unwrap();

    let read = genome_seq[40..80].to_string();
    let fastq = write_fastq(&dir, "reads.fq", &[("r1".to_string(), read)]);
    let bam_path = dir.path().join("out.bam");

    let status = Command::new(probebwa_bin())
        .args([
            "map",
            "--genome", prefix.to_str().unwrap(),
            "--hash", prefix.to_str().unwrap(),
            "-M", fastq.to_str().unwrap(),
            "--outputformat", "bam",
            "--output", bam_path.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let compressed = fs::read(&bam_path).unwrap();
    // BGZF blocks are valid gzip members; the standard `gzip` crate decodes
    // (and, per the gzip spec, stops at) the first member, which is enough
    // to check the BAM magic without a full BGZF multi-block reader.
    let mut decoder = flate2::read::GzDecoder::new(&compressed[..]);
    let mut decompressed = Vec::new();
    std::io::Read::read_to_end(&mut decoder, &mut decompressed).unwrap();
    assert_eq!(&decompressed[0..4], b"BAM\x01", "output should start with the BAM magic after BGZF decompression");
}

#[test]
fn test_labelfilter_restricts_output_to_matching_read_ids() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(100).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = dir.path().join("ref");
    build_genome_index("t", "t1", prefix.to_str().unwrap(), &[fasta.to_str().unwrap()]).unwrap();
    build_hash_table(prefix.to_str().unwrap(), prefix.to_str().unwrap()).unwrap();

    let fastq = write_fastq(&dir, "reads.fq", &[
        ("keep_1".to_string(), genome_seq[40..80].to_string()),
        ("drop_1".to_string(), genome_seq[200..240].to_string()),
    ]);

    let output = Command::new(probebwa_bin())
        .args([
            "map",
            "--genome", prefix.to_str().unwrap(),
            "--hash", prefix.to_str().unwrap(),
            "-M", fastq.to_str().unwrap(),
            "--labelfilter", "keep_",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let sam = String::from_utf8(output.stdout).unwrap();

    assert!(sam.lines().any(|l| l.starts_with("keep_1\t")), "matching read should be present:\n{sam}");
    assert!(!sam.lines().any(|l| l.starts_with("drop_1\t")), "non-matching read should be filtered out:\n{sam}");
}

#[test]
fn test_adapter_strip_trims_3prime_readthrough() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(101).gen_seq(400);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = dir.path().join("ref");
    build_genome_index("t", "t1", prefix.to_str().unwrap(), &[fasta.to_str().unwrap()]).unwrap();
    build_hash_table(prefix.to_str().unwrap(), prefix.to_str().unwrap()).unwrap();

    // A read that reads through into 20bp of adapter sequence past the true
    // 40bp genomic insert -- without stripping, the adapter tail is
    // non-genomic sequence that should prevent a clean full-length match.
    let adapter = "AGATCGGAAGAGCACACGTC";
    let genomic = &genome_seq[40..80];
    let read_with_adapter = format!("{genomic}{adapter}");
    let fastq = write_fastq(&dir, "reads.fq", &[("r1".to_string(), read_with_adapter)]);

    let output = Command::new(probebwa_bin())
        .args([
            "map",
            "--genome", prefix.to_str().unwrap(),
            "--hash", prefix.to_str().unwrap(),
            "-M", fastq.to_str().unwrap(),
            "--adapter-strip", adapter,
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let sam = String::from_utf8(output.stdout).unwrap();
    let record_line = sam.lines().find(|l| l.starts_with("r1\t")).expect("read should map");
    let fields: Vec<&str> = record_line.split('\t').collect();
    assert_eq!(fields[3], "41", "should map at the true genomic start");
    assert_eq!(fields[5], "40M", "adapter tail should be trimmed, leaving a clean 40M");
    assert_eq!(fields[9].len(), 40, "sequence field should be trimmed to the genomic portion");
}

#[test]
fn test_index_produces_a_valid_bai_alongside_the_bam() {
    let dir = TempDir::new().unwrap();
    let genome_seq = Xorshift64(102).gen_seq(2000);
    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome_seq);
    let prefix = dir.path().join("ref");
    build_genome_index("t", "t1", prefix.to_str().unwrap(), &[fasta.to_str().unwrap()]).unwrap();
    build_hash_table(prefix.to_str().unwrap(), prefix.to_str().unwrap()).unwrap();

    // Reads deliberately out of coordinate order in the input, so a
    // non-trivial sort has to happen for the index to be meaningful.
    let fastq = write_fastq(&dir, "reads.fq", &[
        ("late".to_string(), genome_seq[1000..1040].to_string()),
        ("early".to_string(), genome_seq[40..80].to_string()),
    ]);
    let bam_path = dir.path().join("out.bam");
    let bai_path = dir.path().join("out.bam.bai");

    let status = Command::new(probebwa_bin())
        .args([
            "map",
            "--genome", prefix.to_str().unwrap(),
            "--hash", prefix.to_str().unwrap(),
            "-M", fastq.to_str().unwrap(),
            "--outputformat", "bam",
            "--output", bam_path.to_str().unwrap(),
            "--index",
        ])
        .status()
        .unwrap();
    assert!(status.success());
    assert!(bai_path.exists(), "expected a .bai file alongside the .bam");

    // A real BAI starts with the "BAI\1" magic.
    let bai_bytes = fs::read(&bai_path).unwrap();
    assert_eq!(&bai_bytes[0..4], b"BAI\x01");

    // The BAM itself should be coordinate-sorted now (the "early" read's
    // record should come before "late"'s).
    let compressed = fs::read(&bam_path).unwrap();
    let mut decoder = flate2::read::GzDecoder::new(&compressed[..]);
    let mut decompressed = Vec::new();
    std::io::Read::read_to_end(&mut decoder, &mut decompressed).unwrap();
    let early_pos = decompressed.windows(5).position(|w| w == b"early");
    let late_pos = decompressed.windows(4).position(|w| w == b"late");
    assert!(early_pos.is_some() && late_pos.is_some());
    assert!(early_pos < late_pos, "coordinate-sorted BAM should have the lower-position read first");
}
