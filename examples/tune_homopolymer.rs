//! Scratch tuning harness (not part of the test suite): sweeps candidate
//! homopolymer gap-open discount slopes and mismatch-penalty clamp floors
//! against a battery of synthetic deletion/insertion scenarios placed
//! inside vs. outside homopolymer runs, to empirically pick values that
//! maximize correct-CIGAR recall. Parameters are tuned against this
//! project's own accuracy metric.

use probebwa::*;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

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
    fn base(&mut self) -> char {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        BASES[(self.next_u64() % 4) as usize]
    }
}

fn write_fasta(dir: &TempDir, name: &str, contig_name: &str, seq: &str) -> PathBuf {
    let path = dir.path().join(name);
    fs::write(&path, format!(">{}\n{}\n", contig_name, seq)).unwrap();
    path
}
/// Build a genome containing, at a fixed offset, a homopolymer run of length
/// `run_len` (or ordinary random sequence when `run_len <= 1`), then a read
/// that spans that region with a `del_len`bp deletion placed at the run
/// (or, for the non-homopolymer case, at the equivalent offset). `n_errors`
/// low-quality single-base substitution errors are additionally scattered
/// into the read (simulating realistic sequencing noise alongside the
/// indel, which is what actually stresses the mismatch-penalty/gap-open
/// interaction rather than a pristine indel-only read).
fn scenario(seed: u64, run_len: usize, del_len: usize, n_errors: usize) -> (String, String, usize, Vec<usize>) {
    let mut rng = Xorshift64(seed);
    let prefix = rng.gen_seq(200);
    let run_base = rng.base();
    let run: String = if run_len > 1 { std::iter::repeat_n(run_base, run_len).collect() } else { rng.gen_seq(1) };
    let suffix = rng.gen_seq(200);
    let genome = format!("{prefix}{run}{suffix}");

    // Read = genome with `del_len` bases removed starting at the run's start
    // (so the deletion sits inside the homopolymer run when run_len > del_len,
    // matching how real indels concentrate in/near homopolymer runs).
    let del_start = prefix.len();
    let read_start = del_start.saturating_sub(37);
    let read_end = (del_start + 40).min(genome.len());
    let mut read = String::new();
    read.push_str(&genome[read_start..del_start]);
    read.push_str(&genome[(del_start + del_len).min(genome.len())..read_end.max(del_start + del_len)]);

    let mut read: Vec<u8> = read.into_bytes();
    const BASES: [u8; 4] = *b"ACGT";
    let mut error_positions = Vec::new();
    for _ in 0..n_errors {
        let i = (rng.next_u64() as usize) % read.len();
        let mut b = BASES[(rng.next_u64() % 4) as usize];
        while b == read[i] { b = BASES[(rng.next_u64() % 4) as usize]; }
        read[i] = b;
        error_positions.push(i);
    }
    (genome, String::from_utf8(read).unwrap(), read_start + 1, error_positions)
}

/// Like `scenario`, but with per-base quality scores: `n_errors` positions
/// get a low quality score (Q10, i.e. plausibly-a-sequencing-error) and the
/// rest get Q35 (high confidence) — a FASTQ round-trip through
/// `write_fastq_with_quals` instead of the flat-Q40 default, so the
/// mismatch-penalty clamp actually sees varied input.
fn write_fastq_with_quals(dir: &TempDir, name: &str, id: &str, seq: &[u8], quals: &[u8]) -> PathBuf {
    let path = dir.path().join(name);
    let qual_str: String = quals.iter().map(|&q| (q + 33) as char).collect();
    fs::write(&path, format!("@{id}\n{}\n+\n{qual_str}\n", String::from_utf8_lossy(seq))).unwrap();
    path
}

fn evaluate(homopolymer_slope: f64, mismatch_floor: f64) -> (usize, usize) {
    // NOTE: this harness can only *measure* against the aligner's actual
    // compiled-in constants (there's no runtime knob), so it's used as a
    // print-and-inspect tool between manual edits to smith_waterman.rs
    // rather than a fully automated sweep. Slopes/floors are echoed so the
    // operator knows which build they're looking at.
    eprintln!("(evaluating build compiled with homopolymer_slope~{homopolymer_slope}, mismatch_floor~{mismatch_floor})");

    let mut correct = 0usize;
    let mut total = 0usize;
    for &run_len in &[1usize, 3, 5, 8, 12] {
        for &del_len in &[1usize, 2, 3, 5, 8] {
            if del_len >= run_len && run_len > 1 { continue; } // deletion must fit inside the run to test "inside homopolymer"
            for &n_errors in &[0usize, 2, 4] {
                for seed in 0..5u64 {
                    total += 1;
                    let variant_seed = seed * 1000 + run_len as u64 * 10 + del_len as u64 + n_errors as u64 * 100_000;
                    let (genome, read, true_pos, error_positions) = scenario(variant_seed, run_len, del_len, n_errors);
                    let dir = TempDir::new().unwrap();
                    let fasta = write_fasta(&dir, "ref.fa", "chr1", &genome);
                    let prefix = dir.path().join("ref");
                    build_genome_index("t", "t1", prefix.to_str().unwrap(), &[fasta.to_str().unwrap()]).unwrap();
                    build_hash_table(prefix.to_str().unwrap(), prefix.to_str().unwrap()).unwrap();
                    // Flat Q10 throughout (not just the injected-error
                    // positions): a flat, sharply-differentiated quality
                    // scheme (e.g. Q35 correct / Q8 error) gives the DP
                    // enough signal to trivially tell real errors from the
                    // indel, which defeats the point of this stress test —
                    // flat low quality is what actually makes the DP have to
                    // rely on the gap-open/mismatch cost *balance* (rather
                    // than the quality signal) to find the right explanation,
                    // and Q10 (scale 0.33) sits inside the floor values being
                    // swept below, so the floor's value actually matters here
                    // (unlike Q20, where scale 0.67 sat above every floor
                    // tried and the clamp never engaged).
                    let _ = &error_positions;
                    let quals = vec![10u8; read.len()];
                    let fastq = write_fastq_with_quals(&dir, "r.fq", "r", read.as_bytes(), &quals);
                    let recs = map_reads(&prefix, &prefix, &[&fastq], MapperOptions::default()).unwrap();
                    let rec = &recs[0];
                    let deletion_total: u32 = rec.cigar.iter().map(|op| match op {
                        CigarOp::Deletion(n) => *n,
                        _ => 0,
                    }).sum();
                    let close_enough = (rec.position as i64 - true_pos as i64).abs() <= 2;
                    if close_enough && deletion_total as usize == del_len {
                        correct += 1;
                    } else {
                        eprintln!("  MISS run_len={run_len} del_len={del_len} n_errors={n_errors} seed={seed}: pos={} (true {true_pos}) mapq={} cigar={:?}", rec.position, rec.mapq, rec.cigar);
                    }
                }
            }
        }
    }
    (correct, total)
}

fn main() {
    let (correct, total) = evaluate(1.0, 0.2);
    println!("recall: {correct}/{total} ({:.1}%)", 100.0 * correct as f64 / total as f64);
}
