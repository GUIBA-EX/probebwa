//! BAM output.
//!
//! Rather than re-deriving flags/CIGAR/tag encoding a second time against
//! `noodles`'s binary record types, this reuses the existing SAM text
//! formatting (`format_sam_records`) as the single source of truth and lets
//! `noodles_sam::io::Reader` parse it back into typed records, which
//! `noodles_bam::io::Writer` then re-encodes as BGZF-compressed BAM. This
//! keeps the two output formats from silently drifting apart, at the cost of
//! one extra text round-trip — a non-issue at the read volumes this mapper
//! targets.

use anyhow::Result;
use noodles_sam::alignment::io::Write as _;
use std::io::{Cursor, Write};
use std::path::Path;

/// Write `sam_text` (a complete SAM stream: header + records, as produced by
/// `format_sam_records`) to `out` as BAM.
pub fn write_bam(sam_text: &str, out: impl Write) -> Result<()> {
    let mut reader = noodles_sam::io::Reader::new(Cursor::new(sam_text.as_bytes()));
    let header = reader.read_header()?;

    let mut writer = noodles_bam::io::Writer::new(out);
    writer.write_header(&header)?;

    let mut record = noodles_sam::alignment::RecordBuf::default();
    while reader.read_record_buf(&header, &mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;
    }
    writer.try_finish()?;

    Ok(())
}

/// Write `sam_text` as BAM to `bam_path`, then build and write a BAI index
/// alongside it at `{bam_path}.bai` (`--index`).
///
/// `sam_text` must already be coordinate-sorted with its `@HD` line
/// declaring `SO:coordinate` (see `io::sam::sort_for_coordinate_order` /
/// `format_sam_records_sorted`) -- `noodles_bam::fs::index` refuses to index
/// a BAM that isn't marked that way, and the resulting index would be
/// meaningless (and unusable by tools like samtools/IGV) against records
/// that aren't actually in coordinate order regardless.
pub fn write_bam_indexed(sam_text: &str, bam_path: &Path) -> Result<()> {
    let file = std::fs::File::create(bam_path)?;
    write_bam(sam_text, file)?;

    let index = noodles_bam::fs::index(bam_path)?;
    let bai_path = {
        let mut p = bam_path.as_os_str().to_owned();
        p.push(".bai");
        p
    };
    noodles_bam::bai::fs::write(bai_path, &index)?;

    Ok(())
}
