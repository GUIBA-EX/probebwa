//! Banded affine-gap alignment.
//!
//! NOTE: This is a scalar placeholder. For production, replace with:
//!   - `block-aligner` crate (SIMD, fast)
//!   - `bio::alignment::pairwise` (rust-bio, well-tested)
//!   - FFI to KSW2 (minimap2's aligner)
//!
//! Design notes:
//!
//! - **Semi-global in the read**: a read that was placed here by seeding
//!   should, in the overwhelming majority of cases, align end-to-end against
//!   a correctly windowed reference segment — soft-clipping is better
//!   decided by candidate/window selection upstream than by letting the DP
//!   freely truncate either end. We therefore only look for the best score
//!   in the row where the whole read has been consumed (`i == read.len()`),
//!   instead of taking the global max over the whole matrix (which is what a
//!   textbook local Smith-Waterman does).
//! - **Quality-aware mismatch penalty, with gap costs tuned for this
//!   architecture rather than scored as raw linear Phred**: an early version
//!   of this module scored mismatch as the read base's Phred quality
//!   directly (linear, unclamped) and match as exactly 0, with
//!   `--gapopen`/`--gapextend` (defaults 40/3) scaled the same linear way.
//!   Those magnitudes were tried and reverted after testing: transplanted
//!   directly into this module's DP (which searches a generously *padded*
//!   reference window and allows either the M or I state as a valid
//!   endpoint), they let the DP find a *higher-scoring* alternative that
//!   "gives up" on a genuine trailing match run by reporting it as
//!   free-floating Insertions instead of paying for the Deletion that
//!   actually explains it. This wasn't a bug — the alternative really did
//!   score higher under those exact numbers, in this architecture — a
//!   10-seed regression sweep with a 25bp deletion found 6/10 silently
//!   dropping the deletion. `gap_open_phred`/`gap_extend_phred` are still
//!   accepted as configurable inputs (using familiar `--gapopen`/
//!   `--gapextend` flag names/defaults) but are scaled relative to this
//!   module's own validated baseline, not applied as a raw quarter-Phred
//!   magnitude, and mismatch keeps this module's own clamped-quality-ratio
//!   model rather than a linear/unclamped one.
//! - **Gap-open recalibrated after a real-data disagreement:** benchmarking
//!   against real E. coli sequencing data (ENA run DRR002055 vs. the K-12
//!   MG1655 reference, cross-checked against another mapper's output on the
//!   same data) turned up reads where a locally repeat-like window let this
//!   DP trade several genuine mismatches for a few short spurious indels
//!   that happened to land on a cleaner-looking but biologically wrong
//!   register — at the original baseline (-50), 2-3 gap-opens were cheap
//!   enough to be worth it for the resulting extra matches, even with zero
//!   real indel variation present. Raised to -90 (empirically: -80 still
//!   reproduces the bug, -85 already fixes it; -90 keeps a safety margin)
//!   fixes this
//!   (`tests/integration.rs::test_divergent_region_does_not_invent_spurious_indels`,
//!   built from that exact real window/read) while still passing the
//!   original 25bp-deletion regression and the full homopolymer/indel
//!   battery (`examples/tune_homopolymer.rs`) at 100%.
//! - **Homopolymer-aware gap-open discount**: indels are far more common
//!   inside homopolymer runs, so opening a gap there should be cheaper. The
//!   discount curve used here is this module's own approximation of that
//!   idea, tuned against this project's own synthetic accuracy harness
//!   (`examples/tune_homopolymer.rs`), not derived from any external table.
//! - **Banding is offset, not just symmetric**: the reference window handed
//!   in is typically padded on the left (see `mapper/single.rs`), so the true
//!   alignment diagonal is `j ≈ i + ref_offset_hint`, not `j == i`.
//! - **Separate traceback per DP state (M/I/D)**: a single shared traceback
//!   array is not valid for affine-gap DP. The M-state's arrival pointer
//!   ("did the best way to *match* here come from a match, insertion, or
//!   deletion predecessor?") only tells you how to continue if the optimal
//!   path is currently *in* the M state at that cell — reusing it to
//!   reconstruct a path that's actually in the I or D state produces a CIGAR
//!   inconsistent with the score that was actually found. This was a real,
//!   long-standing bug here: confirmed by constructing a read with a clean
//!   25bp deletion where the DP correctly computed the optimal score (which
//!   matches a single clean 25bp deletion exactly) but the single-traceback
//!   version reconstructed a wrong, fragmented CIGAR with real mismatches
//!   hidden inside it — the reported score and the reported CIGAR simply
//!   didn't agree. Each state now gets its own traceback array.

use crate::types::*;

/// Match bonus. Large enough relative to this module's own (deliberately
/// gentle) gap costs to make explaining a read's tail via a genuine
/// reference position always beat giving up on it as free-floating
/// insertions — see module docs.
const MATCH_SCORE: i32 = 20;

/// This module's validated baseline gap costs, reached at the default
/// `--gapopen=40 --gapextend=3`. Non-default values scale proportionally
/// from here rather than using a raw `phred*4` magnitude directly (see
/// module docs for why).
const BASELINE_GAP_OPEN_PHRED: u32 = 40;
const BASELINE_GAP_OPEN_SCORE: f64 = -90.0;
const BASELINE_GAP_EXTEND_PHRED: u32 = 3;
const BASELINE_GAP_EXTEND_SCORE: f64 = -1.0;

/// Default slope of the homopolymer gap-open discount (see
/// `homopolymer_gap_open`) and floor of the mismatch-penalty quality clamp
/// (see `mismatch_penalty`) — both empirically tuned against this project's
/// own synthetic accuracy harness (`examples/tune_homopolymer.rs`), not
/// derived from or matched to any external source.
const DEFAULT_HOMOPOLYMER_SLOPE: f64 = 1.0;
const DEFAULT_MISMATCH_FLOOR: f64 = 0.2;

pub struct SmithWatermanAligner {
    gap_open: i32,
    gap_extend: i32,
    band_width: i64,
    homopolymer_slope: f64,
    mismatch_floor: f64,
}

impl SmithWatermanAligner {
    /// `gap_open_phred`/`gap_extend_phred` are the `--gapopen`/`--gapextend`
    /// CLI flags (defaults 40/3), scaled relative to this module's own
    /// validated baseline (see module docs), not applied as a raw
    /// quarter-Phred magnitude.
    pub fn new(band: usize, gap_open_phred: u32, gap_extend_phred: u32) -> Self {
        let gap_open = (BASELINE_GAP_OPEN_SCORE * gap_open_phred as f64 / BASELINE_GAP_OPEN_PHRED as f64)
            .round() as i32;
        let gap_extend = (BASELINE_GAP_EXTEND_SCORE * gap_extend_phred as f64 / BASELINE_GAP_EXTEND_PHRED as f64)
            .round() as i32;
        Self {
            gap_open: gap_open.min(-1),
            // Keep a non-zero extend cost even at very low gap_extend_phred:
            // the small-but-nonzero granularity is what makes "one long gap"
            // and "several short gaps interspersed with mismatches" resolve
            // consistently instead of hitting exact integer ties.
            gap_extend: gap_extend.min(-1),
            band_width: band as i64,
            homopolymer_slope: DEFAULT_HOMOPOLYMER_SLOPE,
            mismatch_floor: DEFAULT_MISMATCH_FLOOR,
        }
    }

    /// Align `read` (fully consumed) against a window of `ref_seg`, banded
    /// around the diagonal `j = i + ref_offset_hint`.
    ///
    /// Returns `(score, cigar, ref_start_offset)` where `ref_start_offset` is
    /// the 0-based offset into `ref_seg` where the alignment begins (callers
    /// must add this to the window's genome start to get the true SAM
    /// position — the DP is free to start anywhere in the padded window).
    ///
    /// Performance: the DP rows and the traceback matrices are each
    /// allocated once up front and reused for the whole alignment (the row
    /// buffers are swapped rather than reallocated every row, and each
    /// traceback matrix is one flat row-major buffer rather than one `Vec`
    /// per row) — a naive per-row-allocation DP does several heap
    /// allocations per read base, which dominates runtime at these read
    /// lengths; production aligners (ksw2 in minimap2/minibwa, and others)
    /// all avoid this the same way, with a rolling pair of row buffers
    /// instead of a freshly allocated matrix.
    ///
    /// A `thread_local!`-scratch-buffer version of this (reusing the DP rows
    /// and traceback matrices across calls instead of allocating fresh ones
    /// each time) was tried and reverted: it measured consistently *slower*
    /// against real E. coli data (~48s vs. ~33s for the same 20000-pair
    /// benchmark, reproduced across repeated runs, so not measurement noise)
    /// despite removing real per-call allocations, most likely because
    /// routing the DP's hot inner loop through `RefCell`-borrowed struct
    /// fields (vs. plain independent stack-local `Vec`s) gave LLVM less
    /// certainty about alias-freedom between the buffers, costing more in
    /// lost vectorization/hoisting than the allocations it saved — a
    /// plausible explanation, not confirmed with a disassembly-level look,
    /// but the *measurement* (reproduced, and reverted back to this simpler
    /// version to confirm the timing returned to baseline) is what actually
    /// justifies keeping the simpler version here rather than the theory.
    pub fn align(&self, read: &[u8], qual: &[u8], ref_seg: &[u8], ref_offset_hint: i64) -> (i32, Vec<CigarOp>, usize) {
        let m = read.len();
        let n = ref_seg.len();
        if m == 0 || n == 0 {
            return (0, vec![CigarOp::SoftClip(m as u32)], 0);
        }

        let homopolymer_run = precompute_homopolymer_runs(ref_seg);
        let width = n + 1;

        const NEG_INF: i32 = i32::MIN / 2;
        let mut dp_m = vec![0i32; width];
        let mut dp_i = vec![NEG_INF; width];
        let mut dp_d = vec![NEG_INF; width];
        let mut dp_m_next = vec![NEG_INF; width];
        let mut dp_i_next = vec![NEG_INF; width];
        let mut dp_d_next = vec![NEG_INF; width];

        // One flat row-major traceback matrix per state: cell (i, j) lives at
        // i*width + j. tb_m encodes which predecessor STATE (0=M, 1=I, 2=D)
        // fed the match/mismatch step; tb_i/tb_d each encode just whether
        // that state was reached by opening a fresh gap (0, predecessor is
        // M) or extending an existing one (1, predecessor is the same state).
        let mut tb_m = vec![0u8; (m + 1) * width];
        let mut tb_i = vec![0u8; (m + 1) * width];
        let mut tb_d = vec![0u8; (m + 1) * width];

        for i in 1..=m {
            let center = i as i64 + ref_offset_hint;
            let j_start = (center - self.band_width).max(1) as usize;
            let j_end = ((center + self.band_width).min(n as i64)).max(0) as usize;
            if j_start > j_end { continue; }

            dp_m_next.iter_mut().for_each(|v| *v = NEG_INF);
            dp_i_next.iter_mut().for_each(|v| *v = NEG_INF);
            dp_d_next.iter_mut().for_each(|v| *v = NEG_INF);

            let q = *qual.get(i - 1).unwrap_or(&30);
            let row_base = i * width;

            for j in j_start..=j_end {
                let match_sc = if read[i - 1] == ref_seg[j - 1] {
                    MATCH_SCORE
                } else {
                    self.mismatch_penalty(q)
                };

                let from_m = dp_m[j - 1] + match_sc;
                let from_i = dp_i[j - 1] + match_sc;
                let from_d = dp_d[j - 1] + match_sc;
                let m_score = from_m.max(from_i).max(from_d);

                let gap_open_here = self.homopolymer_gap_open(homopolymer_run[j - 1]);

                // Insertion (consume read, stay in ref): open from M at (i-1,j), extend from I at (i-1,j).
                let ins_open = dp_m[j] + gap_open_here;
                let ins_ext = dp_i[j] + self.gap_extend;
                let i_score = ins_open.max(ins_ext);

                // Deletion (consume ref, stay in read): open from M at (i,j-1), extend from D at (i,j-1).
                let del_open = dp_m_next[j - 1] + gap_open_here;
                let del_ext = dp_d_next[j - 1] + self.gap_extend;
                let d_score = del_open.max(del_ext);

                dp_m_next[j] = m_score;
                dp_i_next[j] = i_score;
                dp_d_next[j] = d_score;

                tb_m[row_base + j] = if m_score == from_m { 0 } else if m_score == from_i { 1 } else { 2 };
                tb_i[row_base + j] = if i_score == ins_open { 0 } else { 1 };
                tb_d[row_base + j] = if d_score == del_open { 0 } else { 1 };
            }

            std::mem::swap(&mut dp_m, &mut dp_m_next);
            std::mem::swap(&mut dp_i, &mut dp_i_next);
            std::mem::swap(&mut dp_d, &mut dp_d_next);
        }

        // Semi-global in the read: the whole read must be consumed, so the
        // best endpoint is searched only in row i == m (ref end is free).
        // Only the M and I states are eligible endpoints — both represent
        // "the read's last base was placed here" (matched/mismatched, or
        // inserted); the D state never consumes a read base, so ending in D
        // would mean padding the alignment with reference-only bases the
        // read says nothing about, which is never meaningful to report.
        let (best_j, best_score, best_state) = (0..=n)
            .flat_map(|j| [(j, dp_m[j], 0u8), (j, dp_i[j], 1u8)])
            .max_by_key(|&(_, score, _)| score)
            .unwrap_or((n, 0, 0));

        let (cigar, ref_start_offset) = self.traceback(&tb_m, &tb_i, &tb_d, width, m, best_j, best_state);
        (best_score, cigar, ref_start_offset)
    }

    /// Mismatch penalty scaled by the read base's own Phred quality, capped
    /// to a `[DEFAULT_MISMATCH_FLOOR, 1.0]` fraction of the base penalty.
    /// A linear, unclamped quality-based penalty was tried and reverted (see
    /// module docs): the much larger, unclamped magnitude is part of what
    /// let the DP "give up" on genuine deletions in this architecture's
    /// padded-window search. A low-quality mismatch is cheap (more likely
    /// sequencing error); a high-quality mismatch costs close to the full
    /// base penalty.
    ///
    /// `examples/tune_homopolymer.rs` swept this floor from 0.0 to 1.0
    /// against a synthetic battery of homopolymer/non-homopolymer indels
    /// with injected low-quality sequencing errors: recall is flat (284/285)
    /// for floor <= ~0.3, degrades monotonically above that, and drops to
    /// 274/285 at floor = 1.0 (i.e. the unclamped model) — this project's
    /// original choice of 0.2 was already on the good plateau.
    fn mismatch_penalty(&self, qual: u8) -> i32 {
        let scale = (qual as f64 / 30.0).clamp(self.mismatch_floor, 1.0);
        (-30.0 * scale).round() as i32
    }

    /// Indels are far more common inside homopolymer runs; discount the
    /// gap-open cost by `DEFAULT_HOMOPOLYMER_SLOPE` per extra base of run
    /// length (never past -1, so it never becomes a net gain). The discount
    /// curve is this module's own approximation, tuned against this
    /// project's own synthetic accuracy harness rather than derived from any
    /// external source.
    ///
    /// `examples/tune_homopolymer.rs` swept this slope from 0.5 to 3.0
    /// against the same synthetic battery: recall is flat (284/285) for
    /// slope in roughly [0.5, 2.0], and degrades above that (282/285 at
    /// 3.0) — the original choice of 1.0 was already on the good plateau.
    fn homopolymer_gap_open(&self, run_len: u16) -> i32 {
        let discount = (run_len.saturating_sub(1) as f64 * self.homopolymer_slope).round() as i32;
        (self.gap_open + discount).min(-1)
    }

    /// Traceback from `(read_len, best_j, best_state)` until the read is
    /// fully consumed (`i == 0`). Each of `tb_m`/`tb_i`/`tb_d` is a flat
    /// row-major matrix (`width = n+1` columns per row) holding that state's
    /// own predecessor encoding (see field docs above). Returns `(cigar,
    /// ref_start_offset)`.
    #[allow(clippy::too_many_arguments)]
    fn traceback(&self, tb_m: &[u8], tb_i: &[u8], tb_d: &[u8], width: usize, read_len: usize, best_j: usize, best_state: u8) -> (Vec<CigarOp>, usize) {
        let mut ops = Vec::new();
        let mut i = read_len;
        let mut j = best_j;
        let mut state = best_state;

        while i > 0 {
            if j == 0 {
                // Ran out of reference inside the band; fall back to treating
                // the remainder as insertions rather than indexing out of bounds.
                ops.push(CigarOp::Insertion(1));
                i -= 1;
                continue;
            }
            match state {
                0 => {
                    // M state: consult tb_m to see which state fed the
                    // match/mismatch step, then move diagonally.
                    let parent = tb_m[i * width + j];
                    ops.push(CigarOp::Match(1));
                    i -= 1;
                    j -= 1;
                    state = parent;
                }
                1 => {
                    // I state: tb_i says whether this was a fresh gap-open
                    // (predecessor M) or an extension (stay in I).
                    let opened = tb_i[i * width + j] == 0;
                    ops.push(CigarOp::Insertion(1));
                    i -= 1;
                    state = if opened { 0 } else { 1 };
                }
                _ => {
                    // D state: tb_d says whether this was a fresh gap-open
                    // (predecessor M) or an extension (stay in D). Deletion
                    // doesn't consume a read base, so `i` is untouched.
                    let opened = tb_d[i * width + j] == 0;
                    ops.push(CigarOp::Deletion(1));
                    j -= 1;
                    state = if opened { 0 } else { 2 };
                }
            }
        }
        ops.reverse();
        (collapse_cigar(&ops), j)
    }
}

/// For each position, the length of the homopolymer run it belongs to.
fn precompute_homopolymer_runs(seg: &[u8]) -> Vec<u16> {
    let n = seg.len();
    let mut runs = vec![1u16; n];
    if n == 0 { return runs; }
    let mut start = 0;
    for idx in 1..=n {
        if idx == n || seg[idx] != seg[start] {
            let len = (idx - start).min(u16::MAX as usize) as u16;
            for r in runs.iter_mut().take(idx).skip(start) {
                *r = len;
            }
            start = idx;
        }
    }
    runs
}

fn collapse_cigar(ops: &[CigarOp]) -> Vec<CigarOp> {
    let mut result = Vec::new();
    for op in ops {
        if let Some(last) = result.last_mut() {
            match (last, op) {
                (CigarOp::Match(n), CigarOp::Match(1)) => *n += 1,
                (CigarOp::Insertion(n), CigarOp::Insertion(1)) => *n += 1,
                (CigarOp::Deletion(n), CigarOp::Deletion(1)) => *n += 1,
                (CigarOp::SoftClip(n), CigarOp::SoftClip(m)) => *n += m,
                _ => result.push(*op),
            }
        } else {
            result.push(*op);
        }
    }
    result
}
