//! Stampy forward/backward base-alignment evidence. Coordinates intentionally
//! follow frontend.cc's pre-emission posterior lookup, including the first base.
use super::sensitive::{SensitiveParameters, SensitivePath};

#[derive(Debug)]
pub struct AlignmentConfidence {
    pub quality_string: Vec<u8>,
    pub baq_offsets: Vec<u8>,
    pub probability: f64,
}

const M: usize = 0;
const I: usize = 1;
const D: usize = 2;
const ZERO: f64 = -1e300;
fn add(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    if hi - lo > 36.7 {
        hi
    } else {
        hi + (1.0 + (lo - hi).exp()).ln()
    }
}

/// Calculate both ambiguity paths' evidence inside the already selected window.
/// `band=None` is unbanded; otherwise this is the native fixed band, including
/// possible degenerate widths returned by coarse window selection.
pub fn confidence(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    p: SensitiveParameters,
    band: Option<i32>,
    paths: [&SensitivePath; 2],
) -> Option<[AlignmentConfidence; 2]> {
    let n = read.len();
    let len = reference.len();
    let w = len + 1;
    if n == 0 || qualities.len() != n || p.start_sd <= 0 {
        return None;
    }
    let mut fw = vec![[ZERO; 3]; (n + 1) * w];
    let mut bw = fw.clone();
    let db = -std::f64::consts::LN_10 / 10.0;
    let ge = p.gap_extend * db;
    let go = p.gap_open * db + (1.0 - ge.exp()).ln();
    let emission = |q: usize, r: usize| {
        if read[q] == reference[r] || read[q] == b'N' || reference[r] == b'N' {
            0.0
        } else {
            qualities[q] as f64 * db
        }
    };
    let ins = |q: usize| {
        if read[q] == b'N' {
            0.0
        } else {
            p.insertion_prior.min(qualities[q] as f64 / 2.0) * db
        }
    };
    let prior =
        |r: i64| -((r - p.start_mean as i64) as f64).powi(2) / (2.0 * (p.start_sd as f64).powi(2));
    let bounds = |r: usize, forward: bool| -> (usize, usize) {
        let Some(b) = band else {
            return (1, n);
        };
        let b = b as i64;
        let diag = |r: i64| (2 * r - len as i64 + n as i64 + 1) / 2;
        let lower = if forward && r == 0 {
            1
        } else {
            (diag(r as i64) - b / 2).max(1) as usize
        };
        let upper = if !forward && r == len {
            n
        } else {
            ((b - 1) / 2 + (diag(r as i64 + 1) - 1).max(diag(r as i64)))
                .min(n as i64)
                .max(1) as usize
        };
        if forward {
            (lower, upper.max(lower))
        } else {
            (lower.min(upper), upper)
        }
    };
    for r in 0..w {
        let (lo, hi) = bounds(r, true);
        for q in lo..=hi.min(n) {
            let mut v = [ZERO; 3];
            if r > 0 {
                let previous = fw[(q - 1) * w + r - 1];
                v[M] = if q == 1 {
                    prior(r as i64 - 1)
                } else {
                    add(add(previous[I], previous[D]), previous[M])
                } + emission(q - 1, r - 1);
                let prev = fw[q * w + r - 1];
                v[D] = add(add(prev[I] + go, prev[D] + ge), prev[M] + go);
            }
            v[I] = if q == 1 {
                prior(r as i64) + go
            } else {
                let prev = fw[(q - 1) * w + r];
                add(prev[I] + ge, prev[M] + go)
            } + ins(q - 1);
            fw[q * w + r] = v;
        }
    }
    for r in (0..w).rev() {
        let (lo, hi) = bounds(r, false);
        for q in (lo..=hi.min(n)).rev() {
            let mut v = [ZERO; 3];
            if q < n && r < len {
                let value = emission(q, r) + bw[(q + 1) * w + r + 1][M];
                v = [value; 3];
            }
            if r < len {
                let next = bw[q * w + r + 1][D];
                v[M] = add(v[M], go + next);
                v[I] = add(v[I], go + next);
                v[D] = add(v[D], ge + next);
            }
            if q < n {
                let next = ins(q) + bw[(q + 1) * w + r][I];
                v[M] = add(v[M], go + next);
                v[I] = add(v[I], ge + next);
            }
            if q == n {
                v[M] = add(v[M], 0.0);
                v[I] = add(v[I], 0.0);
            }
            bw[q * w + r] = v;
        }
    }
    let mut total = ZERO;
    // Same summation order as the backward padding recurrence.
    for r in (0..w).rev() {
        let mut v = if r < len {
            prior(r as i64 - 1) + emission(0, r) + bw[w + r + 1][M]
        } else {
            ZERO
        };
        if r < len {
            v = add(v, total);
        }
        total = add(v, go + prior(r as i64) + ins(0) + bw[w + r][I]);
    }
    let project = |path: &SensitivePath| -> Option<AlignmentConfidence> {
        let (mut q, mut r) = (0, path.offset);
        let mut quality_string = Vec::new();
        let mut baq_offsets = vec![64; n];
        let mut minima = std::collections::VecDeque::<(usize, f64)>::new();
        let mut maximum: f64 = 0.0;
        for (&a, &b) in path.read.iter().zip(&path.reference) {
            if q >= n || r > len {
                return None;
            }
            let state = if a == b'-' {
                D
            } else if b == b'-' {
                I
            } else {
                M
            };
            let probability = (fw[q * w + r][state] + bw[q * w + r][state] - total).exp();
            quality_string.push(b'a' + ((0.5 + 25.0 * probability) as i32).clamp(0, 25) as u8);
            let original = qualities[q];
            let adjusted = if probability <= 0.0 {
                0
            } else if probability < 1.0 {
                original.min((-4.34294481903252 * (1.0 - probability).ln()).min(94.0) as u8)
            } else {
                original
            };
            baq_offsets[q] = 64 + original - adjusted;
            let value = if b == b'-' { 0.0 } else { probability };
            while minima.back().is_some_and(|&(_, p)| p >= value) {
                minima.pop_back();
            }
            minima.push_back((q, value));
            while minima.front().is_some_and(|&(x, _)| x + 15 < q) {
                minima.pop_front();
            }
            if q >= 15 {
                maximum = maximum.max(minima.front()?.1);
            }
            if a != b'-' {
                q += 1;
            }
            if b != b'-' {
                r += 1;
            }
        }
        Some(AlignmentConfidence {
            quality_string,
            baq_offsets,
            probability: maximum,
        })
    };
    Some([project(paths[0])?, project(paths[1])?])
}

/// Evidence for an alignment returned by `align_sensitive`, using exactly its
/// cropped reference and shifted location prior.
pub fn sensitive_confidence(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    mut p: SensitiveParameters,
    band: i32,
    alignment: &super::sensitive::SensitiveAlignment,
) -> Option<[AlignmentConfidence; 2]> {
    if band < 0 {
        return confidence(
            read,
            qualities,
            reference,
            p,
            None,
            [&alignment.primary, &alignment.alternate],
        );
    }
    let band = band.max(reference.len() as i32 - read.len() as i32);
    if band <= 15 {
        return confidence(
            read,
            qualities,
            reference,
            p,
            Some(band),
            [&alignment.primary, &alignment.alternate],
        );
    }
    let window = super::stampy_band::find_window(
        read,
        qualities,
        reference,
        p.gap_open as i32,
        p.gap_extend as i32,
        p.insertion_prior as i32,
        band,
    )?;
    p.start_mean -= window.start as i32;
    let mut first = alignment.primary.clone();
    let mut second = alignment.alternate.clone();
    first.offset = first.offset.checked_sub(window.start)?;
    second.offset = second.offset.checked_sub(window.start)?;
    confidence(
        read,
        qualities,
        &reference[window.start..window.start + window.length],
        p,
        Some(window.width),
        [&first, &second],
    )
}
