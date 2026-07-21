//! CIGAR utilities.

use crate::types::CigarOp;
use std::fmt;

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CigarOp::Match(n) => write!(f, "{}M", n),
            CigarOp::Insertion(n) => write!(f, "{}I", n),
            CigarOp::Deletion(n) => write!(f, "{}D", n),
            CigarOp::SoftClip(n) => write!(f, "{}S", n),
            CigarOp::HardClip(n) => write!(f, "{}H", n),
        }
    }
}

/// Reference bases a CIGAR consumes (Match + Deletion; Insertion/SoftClip/
/// HardClip consume the read, not the reference). Used to find a mapped
/// read's rightmost reference coordinate (e.g. for TLEN, or for locating a
/// mate relative to an anchor whose own alignment may include indels), since
/// `position` alone is only the leftmost.
pub fn cigar_reference_span(cigar: &[CigarOp]) -> usize {
    cigar.iter().map(|op| match op {
        CigarOp::Match(n) | CigarOp::Deletion(n) => *n as usize,
        _ => 0,
    }).sum()
}

pub fn parse_cigar(s: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num = 0u32;
    for ch in s.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + ch.to_digit(10).unwrap();
        } else {
            let op = match ch {
                'M' => CigarOp::Match(num),
                'I' => CigarOp::Insertion(num),
                'D' => CigarOp::Deletion(num),
                'S' => CigarOp::SoftClip(num),
                'H' => CigarOp::HardClip(num),
                _ => CigarOp::Match(num),
            };
            ops.push(op);
            num = 0;
        }
    }
    ops
}
