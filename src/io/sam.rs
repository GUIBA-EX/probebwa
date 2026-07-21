//! SAM output formatting.

use crate::types::*;
use std::fmt::Write;

/// A parsed `--readgroup=ID:id,tag:value,...` specification: an `@RG` SAM
/// header line plus the group id every record's `RG:Z:` tag should carry.
pub struct ReadGroup {
    pub id: String,
    header_line: String,
}

impl ReadGroup {
    /// Parse `"ID:id,tag:value,..."` (SAM-format read-group tags,
    /// comma-separated). Requires at least an `ID` tag, matching the SAM
    /// spec's requirement that every `@RG` line have one.
    pub fn parse(spec: &str) -> Result<Self, String> {
        let mut id = None;
        let mut fields = Vec::new();
        for pair in spec.split(',') {
            let mut it = pair.splitn(2, ':');
            let tag = it.next().unwrap_or("");
            let value = it.next()
                .ok_or_else(|| format!("--readgroup: expected TAG:VALUE, got '{pair}'"))?;
            if tag == "ID" {
                id = Some(value.to_string());
            }
            fields.push(format!("{tag}:{value}"));
        }
        let id = id.ok_or_else(|| "--readgroup must set at least an ID tag".to_string())?;
        Ok(Self { id, header_line: format!("@RG\t{}", fields.join("\t")) })
    }
}

/// Format a single record as SAM, optionally tagged with a read-group id
/// (`RG:Z:<id>`).
pub fn format_sam(rec: &AlignmentRecord, read_group_id: Option<&str>) -> String {
    let cigar: String = rec.cigar.iter().map(|c| c.to_string()).collect();
    let flag = sam_flag(rec);
    let mut line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.read_id,
        flag,
        rec.contig_name,
        rec.position,
        rec.mapq,
        if cigar.is_empty() { "*".to_string() } else { cigar },
        rec.mate_contig.as_deref().unwrap_or("*"),
        rec.mate_position.unwrap_or(0),
        rec.insert_size.unwrap_or(0),
        String::from_utf8_lossy(&rec.read_sequence),
        String::from_utf8_lossy(&rec.read_qualities),
    );
    if let Some(id) = read_group_id {
        write!(line, "\tRG:Z:{id}").unwrap();
    }
    line
}

/// Write a full SAM stream: `@HD`, one `@SQ` per reference contig, an `@PG`
/// line identifying this program/version/invocation, an optional `@RG`
/// header line, then every record (tagged with the read group's id when one
/// is given).
///
/// `@SQ` and `@PG` are required by the SAM spec for a well-formed header;
/// without them, tools like samtools/Picard/IGV that validate or index the
/// output will reject or warn on it even though the alignment records
/// themselves are well-formed.
pub fn format_sam_records(
    records: &[AlignmentRecord],
    contigs: &[(String, usize)],
    command_line: &str,
    read_group: Option<&ReadGroup>,
) -> String {
    format_sam_records_sorted(records, contigs, command_line, read_group, false)
}

/// Like `format_sam_records`, but lets the caller declare the stream
/// coordinate-sorted (`SO:coordinate` instead of `SO:unsorted`) -- used for
/// `--index`, which requires both an actually-sorted BAM *and* that header
/// declaration (`noodles_bam::fs::index` refuses to index a BAM that isn't
/// marked `SO:coordinate`, even if the records happen to be in order).
pub fn format_sam_records_sorted(
    records: &[AlignmentRecord],
    contigs: &[(String, usize)],
    command_line: &str,
    read_group: Option<&ReadGroup>,
    coordinate_sorted: bool,
) -> String {
    let mut out = String::new();
    let sort_order = if coordinate_sorted { "coordinate" } else { "unsorted" };
    writeln!(out, "@HD\tVN:1.6\tSO:{sort_order}").unwrap();
    for (name, length) in contigs {
        writeln!(out, "@SQ\tSN:{name}\tLN:{length}").unwrap();
    }
    writeln!(
        out,
        "@PG\tID:probebwa\tPN:probebwa\tVN:{}\tCL:{command_line}",
        env!("CARGO_PKG_VERSION"),
    ).unwrap();
    if let Some(rg) = read_group {
        writeln!(out, "{}", rg.header_line).unwrap();
    }
    for rec in records {
        writeln!(out, "{}", format_sam(rec, read_group.map(|rg| rg.id.as_str()))).unwrap();
    }
    out
}

/// Sort `records` into coordinate order (reference index per `contigs`'
/// order, then position), the order a BAM index requires. Unmapped records
/// (`position == 0`) sort last, per SAM/BAM convention (they have no
/// meaningful coordinate to place in the index).
pub fn sort_for_coordinate_order(records: &mut [AlignmentRecord], contigs: &[(String, usize)]) {
    let contig_index: std::collections::HashMap<&str, usize> = contigs.iter()
        .enumerate().map(|(i, (name, _))| (name.as_str(), i)).collect();
    records.sort_by_key(|r| {
        if r.position == 0 {
            (usize::MAX, usize::MAX)
        } else {
            (contig_index.get(r.contig_name.as_str()).copied().unwrap_or(usize::MAX), r.position)
        }
    });
}

fn sam_flag(rec: &AlignmentRecord) -> u16 {
    let mut flag = 0u16;
    let unmapped = rec.position == 0;
    if rec.is_paired { flag |= 0x01; }
    if rec.is_proper_pair { flag |= 0x02; }
    if unmapped { flag |= 0x04; }
    if rec.mate_unmapped { flag |= 0x08; }
    if !unmapped && !rec.strand.is_forward() { flag |= 0x10; }
    if let Some(mate_strand) = rec.mate_strand {
        if !mate_strand.is_forward() { flag |= 0x20; }
    }
    match rec.read_number {
        Some(1) => flag |= 0x40,
        Some(2) => flag |= 0x80,
        _ => {}
    }
    flag
}
