//! Stampy 1.0.32 statistical primitives, checked against its native extension.
//! Inputs are already substitution-adjusted Phred qualities. The repeat mask
//! is also in Phred units, one entry per read base (not per seed).

use super::stampy_tables::{MISS_CACHE, PHRED_PROB};
use crate::types::CigarOp;

/// Random-match statistic for a fixed alignment. Stampy uses a fixed 3.5
/// alphabet parameter and one genome-sized search space, not a GC estimate.
pub(super) fn alignment_entropy(
    cigar: &[CigarOp],
    read: &[u8],
    quals: &[u8],
    reference: &[u8],
    searchspace: usize,
    min_q: u8,
) -> i32 {
    let (mut qi, mut ri) = (0usize, 0usize);
    let (mut matches, mut variants, mut boundaries, mut deletions) = (0i32, 0i32, 0i32, 0i32);
    let mut previous = 0;
    for &op in cigar {
        match op {
            CigarOp::Insertion(n) if n > 0 => {
                if previous != 1 {
                    boundaries += 1;
                }
                previous = 1;
                qi += n as usize;
            }
            CigarOp::Deletion(n) if n > 0 => {
                if previous != 2 {
                    boundaries += 1;
                    deletions += 1;
                }
                previous = 2;
                ri += n as usize;
            }
            CigarOp::Match(n) if n > 0 => {
                if previous != 0 {
                    boundaries += 1;
                }
                previous = 0;
                for _ in 0..n {
                    if qi >= read.len() || ri >= reference.len() {
                        return i32::MIN / 2;
                    }
                    if quals.get(qi).copied().unwrap_or(0) >= min_q {
                        matches += 1;
                        if read[qi] != reference[ri] && read[qi] != b'N' && reference[ri] != b'N' {
                            variants += 1;
                        }
                    }
                    qi += 1;
                    ri += 1;
                }
            }
            CigarOp::SoftClip(n) => qi += n as usize,
            _ => {}
        }
    }
    let initial = (0.5
        + 10.0 * (f64::from(matches - variants) * 3.5_f64.ln() - (searchspace.max(1) as f64).ln())
            / 10.0_f64.ln()) as i32;
    initial
        - phred_choose(matches, variants)
        - phred_choose(qi as i32, boundaries - deletions)
        - if deletions > 0 {
            (0.5 + f64::from(deletions) * 10.0 * 30.0_f64.ln() / 10.0_f64.ln()) as i32
        } else {
            0
        }
}

fn phred_choose(mut n: i32, mut m: i32) -> i32 {
    if m > n / 2 {
        m = n - m;
    }
    if m == 0 {
        return 0;
    }
    let mut result = 0.0;
    while m > 1 {
        result += f64::from(n).ln() - f64::from(m).ln();
        n -= 1;
        m -= 1;
    }
    result += f64::from(n).ln();
    (0.5 + 10.0 * result / 10.0_f64.ln()) as i32
}

pub(super) fn mutation_statistics(
    qualities: &[u8],
    repeat_phred: &[u16],
    variants: i32,
) -> (i32, i32) {
    let len = qualities.len();
    let mode = if variants == -1 {
        if len <= 34 {
            1
        } else if len <= 49 {
            2
        } else {
            3
        }
    } else {
        variants
    };
    let mut expected = 0.0;
    let mut unmasked = Vec::with_capacity(len);
    let mut masked = Vec::with_capacity(len);
    let mut odds = Vec::with_capacity(len);
    for (i, &q) in qualities.iter().enumerate() {
        let p = if q == 0 {
            0.8
        } else {
            phred_probability(q as u16)
        };
        expected += f64::from(q) * p;
        let repeat = phred_probability(repeat_phred.get(i).copied().unwrap_or(1000));
        let adjusted = (p + (1.0 - p) * repeat).min(0.9);
        unmasked.push(p);
        masked.push(adjusted);
        odds.push(adjusted / (1.0 - adjusted));
    }

    let probability = if mode == 3 {
        let mut average = 0.0;
        for phase in 0..5 {
            let mut all_miss = 1.0;
            for start in (phase..len.saturating_sub(14)).step_by(15) {
                let mut zero = (1.0 - masked[start]) * (1.0 - masked[start + 1]);
                let mut one = odds[start] + odds[start + 1];
                for &p in &unmasked[start + 2..start + 15] {
                    zero *= 1.0 - p;
                    one += p / (1.0 - p);
                }
                all_miss *= 1.0 - zero * (1.0 + one);
            }
            average += all_miss / 5.0;
        }
        average
    } else {
        let mut minimum = 1.0;
        let mut best_length = 0;
        for end in 15..=len {
            let (mut s1, mut s2, mut s3, mut s4) = (0.0, 0.0, 0.0, 0.0);
            let mut zero = 1.0;
            for start in (end.saturating_sub(52)..end).rev() {
                zero *= 1.0 - masked[start];
                let x = odds[start];
                let x2 = x * x;
                let x3 = x2 * x;
                let x4 = x3 * x;
                s1 += x;
                s2 += x2;
                s3 += x3;
                s4 += x4;
                // Preserve evaluation order and cancellation of the native code.
                let e2 = (s1 * s1 - s2) / 2.0;
                let e3 = (s1 * s1 * s1 - 3.0 * s1 * s2 + 2.0 * s3) / 6.0;
                let e4 = (s1 * s1 * s1 * s1 + 8.0 * s3 * s1 + 3.0 * s2 * s2
                    - 6.0 * s2 * s1 * s1
                    - 6.0 * s4)
                    / 24.0;
                let five_plus = 1.0 - zero * (1.0 + s1 + e2 + e3 + e4);
                if five_plus > 0.01 {
                    break;
                }
                let length = end - start;
                if length < 15 {
                    continue;
                }
                let cache = &MISS_CACHE[5 * (length - 15)..];
                let mut miss = zero * cache[0];
                miss += (1.0 - miss).max(0.0) * zero * cache[1] * s1;
                miss += (1.0 - miss).max(0.0) * zero * cache[2] * e2;
                miss += (1.0 - miss).max(0.0) * zero * cache[3] * e3;
                miss += (1.0 - miss).max(0.0) * zero * cache[4] * e4;
                miss += (1.0 - miss).max(0.0) * five_plus;
                if miss < minimum || (miss < minimum + 1e-8 && length > best_length) {
                    minimum = miss;
                    best_length = length;
                }
            }
        }
        minimum
    };
    (
        (10.0 * probability.max(1e-10).ln() / 0.1_f64.ln()) as i32,
        (expected + 0.5) as i32,
    )
}

fn phred_probability(q: u16) -> f64 {
    PHRED_PROB.get(q as usize).copied().unwrap_or(1e-8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_stampy_alignment_entropy() {
        let mut count = 0;
        for line in include_str!("../../tests/fixtures/stampy_entropy.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<_> = line.split('\t').collect();
            let mut read = Vec::new();
            let mut reference = Vec::new();
            let mut cigar = Vec::new();
            for (r, g) in cols[1].bytes().zip(cols[2].bytes()) {
                if r != b'-' {
                    read.push(r);
                }
                if g != b'-' {
                    reference.push(g);
                }
                cigar.push(if r == b'-' {
                    CigarOp::Deletion(1)
                } else if g == b'-' {
                    CigarOp::Insertion(1)
                } else {
                    CigarOp::Match(1)
                });
            }
            let q: Vec<u8> = cols[3].split(',').map(|v| v.parse().unwrap()).collect();
            let actual = alignment_entropy(
                &cigar,
                &read,
                &q,
                &reference,
                cols[4].parse().unwrap(),
                cols[5].parse().unwrap(),
            );
            assert_eq!(actual, cols[6].parse::<i32>().unwrap(), "{}", cols[0]);
            count += 1;
        }
        assert_eq!(count, 400);
    }

    #[test]
    fn native_stampy_mutation_statistics() {
        let mut count = 0;
        for line in include_str!("../../tests/fixtures/stampy_mutation.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<_> = line.split('\t').collect();
            let q: Vec<u8> = cols[1].split(',').map(|v| v.parse().unwrap()).collect();
            let mask: Vec<u16> = cols[2].split(',').map(|v| v.parse().unwrap()).collect();
            let actual = mutation_statistics(&q, &mask, cols[3].parse().unwrap());
            assert_eq!(
                actual,
                (cols[4].parse().unwrap(), cols[5].parse().unwrap()),
                "{}",
                cols[0]
            );
            count += 1;
        }
        assert!(count >= 100);
    }
}
