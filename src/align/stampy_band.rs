//! Exact replay of the original 8-lane coarse aligner used to choose a
//! sensitive-alignment window. Integer labels and wrapping match Stampy 1.0.32
//! c/alignutils.c; replacing this with generic affine DP changes tie placement.

#[derive(Debug)]
pub struct CoarseAlignment {
    pub score: i32,
    pub offset: usize,
    pub read: Vec<u8>,
    pub reference: Vec<u8>,
}

pub fn coarse_align(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    open: i32,
    extend: i32,
    prior: i32,
) -> Option<CoarseAlignment> {
    if read.is_empty() || qualities.len() != read.len() || reference.len() != read.len() + 15 {
        return None;
    }
    #[cfg(target_arch = "x86_64")]
    // SSE2 is part of the x86-64 baseline. Other targets use the same scalar
    // recurrence; neither path changes integer wrapping or traceback ties.
    let (best, best_s, back) =
        unsafe { coarse_forward_sse2(read, qualities, reference, open, extend, prior) };
    #[cfg(not(target_arch = "x86_64"))]
    let (best, best_s, back) =
        coarse_forward_scalar(read, qualities, reference, open, extend, prior);
    let mut s = best_s as isize;
    let mut lane = s / 2 - read.len() as isize;
    let mut y = read.len();
    let mut x = best_s - y;
    let mut state = back[s as usize][lane as usize] & 3;
    s -= 2;
    let mut ar = Vec::new();
    let mut ag = Vec::new();
    while y > 0 {
        let newstate = (back
            .get(usize::try_from(s).ok()?)?
            .get(usize::try_from(lane).ok()?)?
            >> (2 * state))
            & 3;
        match state {
            0 => {
                s -= 2;
                x = x.checked_sub(1)?;
                y -= 1;
                ag.push(*reference.get(x)?);
                ar.push(read[y]);
            }
            1 => {
                lane += s & 1;
                s -= 1;
                y -= 1;
                ag.push(b'-');
                ar.push(read[y]);
            }
            3 => {
                s -= 1;
                lane -= s & 1;
                x = x.checked_sub(1)?;
                ag.push(*reference.get(x)?);
                ar.push(b'-');
            }
            _ => return None,
        }
        state = newstate;
    }
    ar.reverse();
    ag.reverse();
    Some(CoarseAlignment {
        score: (best as i32 + 0x8000) >> 2,
        offset: x,
        read: ar,
        reference: ag,
    })
}

#[cfg(any(test, not(target_arch = "x86_64")))]
fn coarse_forward_scalar(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    open: i32,
    extend: i32,
    prior: i32,
) -> (i16, usize, Vec<[i16; 8]>) {
    const INF: i16 = 0x7800;
    let go = (open * 4) as i16;
    let ge = ((extend - open) * 4) as i16;
    let gi = ((open + prior) * 4) as i16;
    let mut m1 = [INF; 8];
    let mut m2 = m1;
    let mut i1 = m1;
    let mut i2 = m1;
    let mut d1 = m1;
    let mut d2 = m1;
    let mut back = vec![[0_i16; 8]; 2 * (read.len() + 8)];
    let mut best = INF;
    let mut best_s = 0;
    let label =
        |m: &mut [i16; 8], ins: &mut [i16; 8], del: &mut [i16; 8], pointers: &mut [i16; 8]| {
            for k in 0..8 {
                pointers[k] = (m[k] & 3) | ((ins[k] & 3) << 2) | ((del[k] & 3) << 6);
                m[k] &= !3;
                ins[k] = (ins[k] & !3) | 1;
                del[k] = (del[k] & !3) | 3;
            }
        };
    for t in 0..read.len() + 8 {
        if t < 8 {
            m1[t] = i16::MIN;
            m2[t] = i16::MIN;
        }
        for k in 0..8 {
            m1[k] = m1[k].min(i1[k]).min(d1[k]);
        }
        if t >= read.len() && m1[t - read.len()] < best {
            best = m1[t - read.len()];
            best_s = 2 * t;
        }
        for k in 0..8 {
            let q = t.checked_sub(k);
            let (base, quality) = match q {
                Some(q) if q < read.len() => (read[q] as i16, qualities[q] as i16 * 4),
                Some(_) => (b'0' as i16, 256),
                None => (INF, 256),
            };
            let refbase = reference.get(t + k).copied().unwrap_or(b'G') as i16;
            m1[k] = m1[k].wrapping_add(if base == refbase { 0 } else { quality });
            d1[k] = if k == 0 {
                INF
            } else {
                d2[k - 1].wrapping_add(ge).min(m2[k - 1]).wrapping_add(go)
            };
            i1[k] = i2[k].wrapping_add(ge).min(m2[k]).wrapping_add(gi);
        }
        label(&mut m1, &mut i1, &mut d1, &mut back[2 * t]);
        for k in 0..8 {
            m2[k] = m2[k].min(i2[k]).min(d2[k]);
        }
        if t >= read.len() && m2[t - read.len()] < best {
            best = m2[t - read.len()];
            best_s = 2 * t + 1;
        }
        for k in 0..8 {
            let q = t.checked_sub(k);
            let (base, quality) = match q {
                Some(q) if q < read.len() => (read[q] as i16, qualities[q] as i16 * 4),
                Some(_) => (b'0' as i16, 256),
                None => (INF, 256),
            };
            let refbase = reference.get(t + k + 1).copied().unwrap_or(b'G') as i16;
            m2[k] = m2[k].wrapping_add(if base == refbase { 0 } else { quality });
            d2[k] = d1[k].wrapping_add(ge).min(m1[k]).wrapping_add(go);
            i2[k] = if k == 7 {
                INF
            } else {
                i1[k + 1].wrapping_add(ge).min(m1[k + 1]).wrapping_add(gi)
            };
        }
        label(&mut m2, &mut i2, &mut d2, &mut back[2 * t + 1]);
    }
    (best, best_s, back)
}

#[derive(Debug, PartialEq, Eq)]
pub struct SensitiveWindow {
    pub start: usize,
    pub length: usize,
    // Native findband can return -1 for a degenerate cropped window.
    pub width: i32,
}

/// Original coarse pre-alignment and exact-match-run window selection.
/// `band` is the requested width after expanding to the reference/read length
/// difference. Native callers bypass this routine for widths <= 15.
pub fn find_window(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    open: i32,
    extend: i32,
    prior: i32,
    band: i32,
) -> Option<SensitiveWindow> {
    if read.is_empty() || qualities.len() != read.len() || band <= 15 {
        return None;
    }
    let n = i32::try_from(read.len()).ok()?;
    let len_ref = i32::try_from(reference.len()).ok()?;
    let divisor = 1 + (qualities.iter().map(|&q| q as i32).sum::<i32>() + 256) / 0x3a00;
    let q: Vec<u8> = qualities
        .iter()
        .map(|&q| (q as i32 / divisor) as u8)
        .collect();
    let scale = |v: i32| if v == 0 { 0 } else { 1 + (v - 1) / divisor };
    let (open, extend, prior) = (scale(open), scale(extend), scale(prior));
    let mut padded = vec![b'G'; 16];
    padded.extend_from_slice(reference);
    padded.extend_from_slice(&[b'G'; 16]);
    let mut min_band = band;
    let mut max_band = -band;
    let mut min_width = 25;
    let mut spine = 0;
    let mut max_band_len = 0;
    let mut local_len = vec![0; read.len()];
    let mut local_band = vec![0; read.len()];
    let mut band_len = vec![0; band as usize + 32];
    let max_index = (band / 2 + 14 - (-band / 2) + 14) / 15;
    for half in 0..max_index - max_index / 2 {
        if min_width <= 1 {
            break;
        }
        for side in 0..2 {
            if side == 1 && (max_index / 2 - 1 - half < 0 || min_width <= 1) {
                break;
            }
            let mut i = -band / 2
                + if side == 0 {
                    max_index / 2 + half
                } else {
                    max_index / 2 - 1 - half
                } * 15;
            if i >= band / 2 {
                i = band / 2 - 1;
            }
            let mut ypos = (len_ref - n) / 2 + i;
            let mut xmin = 0;
            let mut xmax = n;
            let mut match_size = 8;
            if ypos < -15 {
                xmin = -15 - ypos;
                ypos = -15;
                match_size = 10;
            }
            if ypos + n >= len_ref {
                xmax = len_ref - ypos;
                match_size = 10;
            }
            if xmin >= xmax {
                return None;
            }
            let start = (16 + ypos) as usize;
            let length = (xmax - xmin) as usize;
            let a = coarse_align(
                &read[xmin as usize..xmax as usize],
                &q[xmin as usize..xmax as usize],
                padded.get(start..start + length + 15)?,
                open,
                extend,
                prior,
            )?;
            if a.score < open {
                min_width = 1;
            }
            let (mut x, mut y, mut matches) = (0_i32, 0_i32, 0_i32);
            for j in 0..a.read.len() {
                let equal = a.read[j] == a.reference[j];
                let last = j + 1 == a.read.len();
                if equal {
                    matches += 1;
                }
                if !equal || last {
                    if matches > match_size {
                        // Preserve the native terminal-run convention, including
                        // its coordinate update when the last column mismatches.
                        if last {
                            x += 1;
                            y += 1;
                        }
                        let band_y = a.offset as i32 + y - x + i;
                        let accum = band_len.get_mut(usize::try_from(band_y + band / 2).ok()?)?;
                        *accum += matches;
                        if *accum > max_band_len {
                            max_band_len = *accum;
                            spine = band_y;
                        }
                        for k in 1..=matches {
                            let at = usize::try_from(x - k).ok()?;
                            if *local_len.get(at)? < matches {
                                local_len[at] = matches;
                                local_band[at] = band_y;
                            }
                        }
                    }
                    matches = 0;
                }
                if a.read[j] != b'-' {
                    x += 1;
                }
                if a.reference[j] != b'-' {
                    y += 1;
                }
            }
        }
    }
    for (&length, &b) in local_len.iter().zip(&local_band) {
        if length > 0 && spine - 75 < b && b < spine + 75 {
            min_band = min_band.min(b);
            max_band = max_band.max(b);
        }
    }
    if min_band == band {
        min_band = (-band / 2).max(-min_width / 2);
        max_band = (band / 2).min(min_width / 2);
    }
    if max_band - min_band < min_width {
        let extra = min_width - (max_band - min_band);
        min_band = (min_band - extra / 2).max(-band / 2);
        max_band = (max_band + extra / 2).min(band / 2);
    }
    let diagonal = (len_ref - n) / 2;
    let start = (diagonal + min_band).max(0);
    let end = (diagonal + max_band + n).min(len_ref);
    let length = end - start;
    let new_diagonal = (length - n) / 2 + start;
    let width = (-2 * (diagonal + min_band - new_diagonal) + 1)
        .max(2 * (diagonal + max_band - new_diagonal) + 1);
    Some(SensitiveWindow {
        start: usize::try_from(start).ok()?,
        length: usize::try_from(length).ok()?,
        width,
    })
}

/// # Safety
/// SSE2 must be available. Inputs must satisfy coarse_align's nonempty read,
/// matching quality length, and read-length-plus-15 reference contract.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn coarse_forward_sse2(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    open: i32,
    extend: i32,
    prior: i32,
) -> (i16, usize, Vec<[i16; 8]>) {
    use std::arch::x86_64::*;
    const INF: i16 = 0x7800;
    let n = read.len();
    // Each load follows query positions t,t-1,...,t-7. Prefix and suffix
    // sentinels match the scalar recurrence, including short reads.
    let mut bases = vec![INF; n + 16];
    let mut quals = vec![256_i16; n + 16];
    bases[..8].fill(b'0' as i16);
    for q in 0..n {
        bases[n + 7 - q] = read[q] as i16;
        quals[n + 7 - q] = qualities[q] as i16 * 4;
    }
    let mut reference = reference.to_vec();
    reference.push(b'G');
    let inf = _mm_set1_epi16(INF);
    let go = _mm_set1_epi16((open * 4) as i16);
    let ge = _mm_set1_epi16(((extend - open) * 4) as i16);
    let gi = _mm_set1_epi16(((open + prior) * 4) as i16);
    let zero = _mm_setzero_si128();
    let lanes = _mm_set_epi16(7, 6, 5, 4, 3, 2, 1, 0);
    let clear = _mm_set1_epi16(!3);
    let low = _mm_set1_epi16(3);
    let one = _mm_set1_epi16(1);
    let low_inf = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, INF);
    let high_inf = _mm_set_epi16(INF, 0, 0, 0, 0, 0, 0, 0);
    let (mut m1, mut m2, mut i1, mut i2, mut d1, mut d2) = (inf, inf, inf, inf, inf, inf);
    let mut back = vec![[0_i16; 8]; 2 * (n + 8)];
    let (mut best, mut best_s) = (INF, 0);
    let mut scores = [0_i16; 8];
    for t in 0..n + 8 {
        if t < 8 {
            let mask = _mm_cmpeq_epi16(lanes, _mm_set1_epi16(t as i16));
            let start = _mm_and_si128(mask, _mm_set1_epi16(i16::MIN));
            m1 = _mm_or_si128(_mm_andnot_si128(mask, m1), start);
            m2 = _mm_or_si128(_mm_andnot_si128(mask, m2), start);
        }
        m1 = _mm_min_epi16(_mm_min_epi16(m1, i1), d1);
        if t >= n {
            _mm_storeu_si128(scores.as_mut_ptr().cast(), m1);
            if scores[t - n] < best {
                best = scores[t - n];
                best_s = 2 * t;
            }
        }
        // All 16-byte query loads and 8-byte reference loads fit their
        // padded buffers. Backpointer stores fit the allocated row exactly.
        let base = _mm_loadu_si128(bases.as_ptr().add(n + 7 - t).cast());
        let quality = _mm_loadu_si128(quals.as_ptr().add(n + 7 - t).cast());
        let reference1 = _mm_unpacklo_epi8(_mm_loadl_epi64(reference.as_ptr().add(t).cast()), zero);
        m1 = _mm_add_epi16(
            m1,
            _mm_andnot_si128(_mm_cmpeq_epi16(base, reference1), quality),
        );
        d1 = _mm_or_si128(
            _mm_slli_si128::<2>(_mm_add_epi16(_mm_min_epi16(_mm_add_epi16(d2, ge), m2), go)),
            low_inf,
        );
        i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(i2, ge), m2), gi);
        let pointers = _mm_or_si128(
            _mm_and_si128(m1, low),
            _mm_or_si128(
                _mm_slli_epi16::<2>(_mm_and_si128(i1, low)),
                _mm_slli_epi16::<6>(_mm_and_si128(d1, low)),
            ),
        );
        _mm_storeu_si128(back[2 * t].as_mut_ptr().cast(), pointers);
        m1 = _mm_and_si128(m1, clear);
        i1 = _mm_or_si128(_mm_and_si128(i1, clear), one);
        d1 = _mm_or_si128(_mm_and_si128(d1, clear), low);

        m2 = _mm_min_epi16(_mm_min_epi16(m2, i2), d2);
        if t >= n {
            _mm_storeu_si128(scores.as_mut_ptr().cast(), m2);
            if scores[t - n] < best {
                best = scores[t - n];
                best_s = 2 * t + 1;
            }
        }
        let reference2 =
            _mm_unpacklo_epi8(_mm_loadl_epi64(reference.as_ptr().add(t + 1).cast()), zero);
        m2 = _mm_add_epi16(
            m2,
            _mm_andnot_si128(_mm_cmpeq_epi16(base, reference2), quality),
        );
        d2 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(d1, ge), m1), go);
        i2 = _mm_or_si128(
            _mm_srli_si128::<2>(_mm_add_epi16(_mm_min_epi16(_mm_add_epi16(i1, ge), m1), gi)),
            high_inf,
        );
        let pointers = _mm_or_si128(
            _mm_and_si128(m2, low),
            _mm_or_si128(
                _mm_slli_epi16::<2>(_mm_and_si128(i2, low)),
                _mm_slli_epi16::<6>(_mm_and_si128(d2, low)),
            ),
        );
        _mm_storeu_si128(back[2 * t + 1].as_mut_ptr().cast(), pointers);
        m2 = _mm_and_si128(m2, clear);
        i2 = _mm_or_si128(_mm_and_si128(i2, clear), one);
        d2 = _mm_or_si128(_mm_and_si128(d2, clear), low);
    }
    (best, best_s, back)
}

#[cfg(all(test, target_arch = "x86_64"))]
mod simd_tests {
    use super::*;

    #[test]
    fn simd_matches_scalar_scores_and_every_backpointer() {
        let mut state = 918273_u64;
        let mut random = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as usize
        };
        for n in [1, 2, 7, 8, 9, 15, 16, 31, 100, 251, 600, 4096] {
            for (open, extend, prior) in [(40, 3, 5), (0, 0, 0), (5, 20, 1), (10000, 1, 19)] {
                let read: Vec<_> = (0..n).map(|_| b"ACGTN"[random() % 5]).collect();
                let qualities: Vec<_> = (0..n).map(|_| (random() % 256) as u8).collect();
                let reference: Vec<_> = (0..n + 15).map(|_| b"ACGTN"[random() % 5]).collect();
                let scalar =
                    coarse_forward_scalar(&read, &qualities, &reference, open, extend, prior);
                // x86-64 always provides SSE2; buffer sizes satisfy the kernel.
                let simd = unsafe {
                    coarse_forward_sse2(&read, &qualities, &reference, open, extend, prior)
                };
                assert_eq!(simd, scalar, "length {n}, costs {open}/{extend}/{prior}");
            }
        }
    }
}
