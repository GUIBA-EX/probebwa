//! Candidate mapping position clustering.
//!
//! Seed hits are grouped by (contig, strand, approximate diagonal) to infer
//! candidate loci. Strand comes directly from the hash table entry (see
//! `index::HashTable::unpack_position`) rather than being inferred from a
//! trend in `ref_pos`, since that heuristic silently fails whenever a true
//! reverse-strand locus doesn't happen to produce hits in read-position order.
//!
//! The diagonal invariant differs by strand:
//!   - Forward: `ref_pos - read_pos` is constant for a true forward mapping
//!     (it's the genome start position of the alignment).
//!   - Reverse: for a reverse-complement hit, `ref_pos + read_pos` is constant
//!     instead (derivable from how `HashTableBuilder` stores RC k-mer
//!     positions: `ref_pos = contig_len - i - k` where `i` is the RC-strand
//!     index that coincides with `read_pos` at a true match, so
//!     `ref_pos + read_pos = genome_start + read_len - k`).
//!
//! Within each (contig, strand) group, hits are sorted by diagonal and merged
//! with a sliding window (consecutive hits within `MERGE_DISTANCE` join the
//! same cluster) rather than binned into a fixed-width grid — a fixed grid
//! can split hits from one true locus across a bin boundary.

use crate::types::*;
use rustc_hash::FxHashMap;

/// Max gap between consecutive sorted diagonals to still merge into one cluster.
const MERGE_DISTANCE: i64 = 50;

/// Cluster seed hits into candidate mapping positions.
pub fn cluster_hits(hits: &[SeedHit], read_len: usize) -> Vec<Candidate> {
    // Group by (contig_id, strand) first; diagonal is only meaningful within
    // a fixed strand (mixing strands would fold two unrelated invariants
    // together).
    let mut groups: FxHashMap<(usize, Strand), Vec<(i64, &SeedHit)>> = FxHashMap::default();

    for hit in hits {
        let diagonal = match hit.hit_strand {
            Strand::Forward => hit.ref_pos as i64 - hit.read_pos as i64,
            Strand::Reverse => hit.ref_pos as i64 + hit.read_pos as i64,
        };
        groups.entry((hit.ref_contig, hit.hit_strand)).or_default().push((diagonal, hit));
    }

    let mut candidates: Vec<Candidate> = Vec::new();

    for ((cid, strand), mut diag_hits) in groups {
        diag_hits.sort_by_key(|(d, _)| *d);

        let mut cluster_start = 0usize;
        for idx in 1..=diag_hits.len() {
            let boundary = idx == diag_hits.len()
                || diag_hits[idx].0 - diag_hits[idx - 1].0 > MERGE_DISTANCE;
            if boundary {
                let group = &diag_hits[cluster_start..idx];
                let avg_diag = group.iter().map(|(d, _)| *d).sum::<i64>() / group.len() as i64;

                // Forward: diagonal already equals the genome start position.
                // Reverse: genome start = diagonal - read_len + k.
                let position = match strand {
                    Strand::Forward => avg_diag,
                    Strand::Reverse => avg_diag - read_len as i64 + group[0].1.k as i64,
                };

                candidates.push(Candidate {
                    contig_id: cid,
                    position: position.max(0) as usize,
                    strand,
                    score: group.len() as f64,
                });

                cluster_start = idx;
            }
        }
    }

    candidates.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
    candidates
}
