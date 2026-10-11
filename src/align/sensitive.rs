//! Scalar reference implementation of Stampy's sensitive Viterbi model.
//!
//! Used by the production mapper and independently verified against native calls.
//! Includes coarse window selection, band traversal, ambiguity and early bailout.
//! Forward/backward confidence is implemented in the sibling confidence module.
//! Model/trace order: Stampy 1.0.32 ext/readalign/readalign.cc and frontend.cc.

#[derive(Clone, Copy, Debug)]
pub struct SensitiveParameters {
    pub insertion_prior: f64,
    pub gap_open: f64,
    pub gap_extend: f64,
    pub start_mean: i32,
    pub start_sd: i32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SensitivePath {
    pub read: Vec<u8>,
    pub reference: Vec<u8>,
    pub offset: usize,
}

#[derive(Debug)]
pub struct SensitiveAlignment {
    pub score: i32,
    /// Native ambiguity bias 1.001; the native reported offset belongs here
    /// even when its Python wrapper exchanges the primary and alternate paths.
    pub primary: SensitivePath,
    /// Native ambiguity bias 0.999.
    pub alternate: SensitivePath,
}

const M: usize = 0;
const I: usize = 1;
const D: usize = 2;
const PAD: usize = 3;

/// Consume the whole read, allowing free reference prefix/suffix padding.
/// Qualities are numeric Phred scores, not FASTQ characters.
pub fn align_unbanded(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    parameters: SensitiveParameters,
) -> Option<SensitiveAlignment> {
    align_model(read, qualities, reference, parameters, None, None)
}

/// The fixed diagonal band used by the native sensitive core. This does not
/// include `findband` (the separate coarse pre-alignment/window selection).
pub fn align_fixed_band(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    parameters: SensitiveParameters,
    band: i32,
) -> Option<SensitiveAlignment> {
    align_model(
        read,
        qualities,
        reference,
        parameters,
        Some(band as i64),
        None,
    )
}

/// Native window selection, Viterbi score/path and score-cutoff behavior.
/// Returned paths retain native trace order; wrapper modes 3/5 exchange their
/// strings (but retain the primary offset). BAQ is not implemented here.
pub fn align_sensitive(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    mut parameters: SensitiveParameters,
    band: i32,
    max_score: i32,
) -> Option<SensitiveAlignment> {
    if band < 0 {
        return align_model(
            read,
            qualities,
            reference,
            parameters,
            None,
            Some(max_score),
        );
    }
    let band = band.max(i32::try_from(reference.len()).ok()? - i32::try_from(read.len()).ok()?);
    if band <= 15 {
        return align_model(
            read,
            qualities,
            reference,
            parameters,
            Some(band as i64),
            Some(max_score),
        );
    }
    let window = super::stampy_band::find_window(
        read,
        qualities,
        reference,
        parameters.gap_open as i32,
        parameters.gap_extend as i32,
        parameters.insertion_prior as i32,
        band,
    )?;
    parameters.start_mean = parameters
        .start_mean
        .checked_sub(i32::try_from(window.start).ok()?)?;
    let mut result = align_model(
        read,
        qualities,
        reference.get(window.start..window.start + window.length)?,
        parameters,
        Some(window.width as i64),
        Some(max_score),
    )?;
    result.primary.offset += window.start;
    result.alternate.offset += window.start;
    Some(result)
}

fn align_model(
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    parameters: SensitiveParameters,
    band: Option<i64>,
    max_score: Option<i32>,
) -> Option<SensitiveAlignment> {
    if read.is_empty()
        || qualities.len() != read.len()
        || parameters.start_sd <= 0
        || !parameters.insertion_prior.is_finite()
        || parameters.insertion_prior < 0.0
        || !parameters.gap_open.is_finite()
        || parameters.gap_open < 0.0
        || !parameters.gap_extend.is_finite()
        || parameters.gap_extend <= 0.0
    {
        return None;
    }
    let width = reference.len().checked_add(1)?;
    // Native Logspace::clear is finite (-1e300), not log(0). This affects
    // degenerate bands: their traceback still chooses the first transition.
    const CLEARED: f64 = 10.0 / std::f64::consts::LN_10 * 1e300;
    // Reproduce the backward iterator, including its endpoint and C++ integer
    // division toward zero. Its first/last columns are not a symmetric mask.
    let bounds = |r: usize| -> (usize, usize) {
        if let Some(b) = band {
            let diagonal = |r: i64| (2 * r - reference.len() as i64 + read.len() as i64 + 1) / 2;
            let lower = (diagonal(r as i64) - b / 2).max(1) as usize;
            let upper = if r == reference.len() {
                read.len()
            } else {
                ((b - 1) / 2 + (diagonal(r as i64 + 1) - 1).max(diagonal(r as i64)))
                    .min(read.len() as i64)
                    .max(1) as usize
            };
            (lower.min(upper), upper)
        } else {
            (1, read.len())
        }
    };
    // Store only rows visited by each reference column. All other cells
    // retain the native finite clear value through `at`; in particular, an
    // unvisited terminal row must not be mistaken for the free suffix.
    let mut columns = Vec::with_capacity(width);
    let mut size = 0usize;
    for r in 0..width {
        let (lower, upper) = bounds(r);
        columns.push((lower, upper, size));
        size = size.checked_add(upper - lower + 1)?;
    }
    let mut dp = vec![[CLEARED; 3]; size];
    let at = |dp: &[[f64; 3]], q: usize, r: usize, state: usize| {
        let (lower, upper, offset) = columns[r];
        if q < lower || q > upper {
            CLEARED
        } else {
            dp[offset + q - lower][state]
        }
    };
    let open =
        parameters.gap_open - 10.0 * (1.0 - 10_f64.powf(-parameters.gap_extend / 10.0)).log10();
    let insertion: Vec<f64> = read
        .iter()
        .zip(qualities)
        .map(|(&base, &q)| {
            if base == b'N' {
                0.0
            } else {
                parameters.insertion_prior.min(q as f64 / 2.0)
            }
        })
        .collect();
    let emission = |q: usize, r: usize| {
        if read[q] == reference[r] || read[q] == b'N' || reference[r] == b'N' {
            0.0
        } else {
            qualities[q] as f64
        }
    };
    let location = |r: usize, matching: bool| {
        let delta = r as f64 - f64::from(matching) - parameters.start_mean as f64;
        10.0 / std::f64::consts::LN_10 * delta * delta
            / (2.0 * (parameters.start_sd as f64).powi(2))
    };
    for r in (0..width).rev() {
        let (lower, upper, offset) = columns[r];
        if upper == read.len() {
            dp[offset + read.len() - lower] = [0.0, 0.0, f64::INFINITY];
        }
        for q in (lower..=upper.min(read.len() - 1)).rev() {
            let matching = if r < reference.len() {
                emission(q, r) + at(&dp, q + 1, r + 1, M)
            } else {
                f64::INFINITY
            };
            let deleting = if r < reference.len() {
                at(&dp, q, r + 1, D)
            } else {
                f64::INFINITY
            };
            let inserting = insertion[q] + at(&dp, q + 1, r, I);
            dp[offset + q - lower] = [
                matching.min(open + deleting).min(open + inserting),
                matching
                    .min(open + deleting)
                    .min(parameters.gap_extend + inserting),
                matching.min(parameters.gap_extend + deleting),
            ];
        }
    }
    let mut pad = vec![f64::INFINITY; width];
    for r in (0..width).rev() {
        let matching = if r < reference.len() {
            location(r, true) + emission(0, r) + at(&dp, 1, r + 1, M)
        } else {
            f64::INFINITY
        };
        let skipping = if r < reference.len() {
            pad[r + 1]
        } else {
            f64::INFINITY
        };
        pad[r] = matching
            .min(skipping)
            .min(open + location(r, false) + insertion[0] + at(&dp, 1, r, I));
    }
    // The reference x86 build converts an out-of-range double to INT_MIN;
    // frontend.cc clamps it to zero only after checking the bailout threshold.
    // Define that observed behavior explicitly rather than relying on a Rust
    // saturating cast (or reproducing C++'s undefined conversion).
    let raw_score = if pad[0] + 0.5 > i32::MAX as f64 {
        i32::MIN
    } else {
        (0.5 + pad[0]) as i32
    };
    let score = raw_score.max(0);
    if max_score.is_some_and(|max| raw_score > max) {
        let path = SensitivePath {
            read: read.to_vec(),
            reference: reference.get(..read.len())?.to_vec(),
            offset: 0,
        };
        return Some(SensitiveAlignment {
            score: score.clamp(0, 65535),
            primary: path.clone(),
            alternate: path,
        });
    }
    let trace = |bias: f64| -> Option<SensitivePath> {
        let tolerance = 10.0 * bias.log10();
        let (mut q, mut r, mut state) = (0, 0, PAD);
        let mut offset = 0;
        let mut aligned_read = Vec::new();
        let mut aligned_ref = Vec::new();
        let mut padding = true;
        while q < read.len() {
            let mut best = f64::INFINITY;
            let mut next = None;
            // Preserve native transition visitation order. Bias applies only
            // to the challenger, not to the stored incumbent's probability.
            let mut consider = |cost: f64, target: usize| {
                if cost - tolerance < best {
                    best = cost;
                    next = Some(target);
                }
            };
            if state == PAD && r < reference.len() {
                consider(pad[r + 1], PAD);
            }
            if r < reference.len() {
                consider(
                    emission(q, r)
                        + at(&dp, q + 1, r + 1, M)
                        + if state == PAD { location(r, true) } else { 0.0 },
                    M,
                );
            }
            if state != PAD && r < reference.len() {
                consider(
                    at(&dp, q, r + 1, D)
                        + if state == D {
                            parameters.gap_extend
                        } else {
                            open
                        },
                    D,
                );
            }
            if state != D {
                consider(
                    insertion[q]
                        + at(&dp, q + 1, r, I)
                        + if state == I {
                            parameters.gap_extend
                        } else {
                            open
                        }
                        + if state == PAD {
                            location(r, false)
                        } else {
                            0.0
                        },
                    I,
                );
            }
            let next = next?;
            if state == PAD && next != PAD {
                offset = r;
            }
            match next {
                PAD => r += 1,
                M => {
                    padding = false;
                    aligned_read.push(read[q]);
                    aligned_ref.push(reference[r]);
                    q += 1;
                    r += 1;
                }
                I => {
                    aligned_read.push(read[q]);
                    aligned_ref.push(b'-');
                    q += 1;
                }
                D => {
                    if band.is_some() && padding {
                        offset += 1;
                    } else {
                        aligned_read.push(b'-');
                        aligned_ref.push(reference[r]);
                    }
                    r += 1;
                }
                _ => unreachable!(),
            }
            state = next;
        }
        Some(SensitivePath {
            read: aligned_read,
            reference: aligned_ref,
            offset,
        })
    };
    Some(SensitiveAlignment {
        score,
        primary: trace(1.001)?,
        alternate: trace(0.999)?,
    })
}
