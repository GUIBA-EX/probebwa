//! CLI entry point.

use mimalloc::MiMalloc;

/// This mapper's hot paths (banded-DP alignment, seed/candidate lists) do a
/// lot of small, short-lived `Vec` allocations per read; mimalloc is
/// consistently faster than the system allocator for that pattern, and
/// swapping it in is a pure performance change with no effect on program
/// behavior (only the binary opts in -- the library itself doesn't impose an
/// allocator choice on downstream consumers).
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use clap::{Parser, Subcommand};
use probebwa::{build_genome_index, build_hash_table, map_reads_with_format, MapperOptions, InputFormat};
use probebwa::index::GenomeIndex;
use probebwa::io::{format_sam_records, format_sam_records_sorted, sort_for_coordinate_order, write_bam, write_bam_indexed, ReadGroup, PHRED33_OFFSET, PHRED64_OFFSET};
use anyhow::{anyhow, Result};
use std::env;
use std::fs::File;
use std::io::{self, Write};

#[derive(Parser)]
#[command(name = "probebwa")]
#[command(about = "Rust mapper for aligning UCE bait/probe sequences to a reference genome")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build genome index (.stidx)
    BuildGenome {
        #[arg(short, long)]
        species: String,
        #[arg(short, long)]
        assembly: String,
        #[arg(short = 'G')]
        output: String,
        files: Vec<String>,
    },
    /// Build 15-mer hash table (.sthash)
    BuildHash {
        #[arg(short, long)]
        genome: String,
        #[arg(short = 'H')]
        output: String,
    },
    /// Map reads
    Map {
        #[arg(short, long)]
        genome: String,
        #[arg(short, long)]
        hash: String,
        #[arg(short, long, default_value_t = 0.001)]
        substitution_rate: f64,
        #[arg(short = 't', long, default_value_t = 1)]
        threads: usize,
        /// Gap open penalty, Phred scale.
        #[arg(long, default_value_t = 40)]
        gapopen: u32,
        /// Gap extension penalty, Phred scale.
        #[arg(long, default_value_t = 3)]
        gapextend: u32,
        /// Prior probability (Phred scale) that a discordant paired-end
        /// separation reflects a real structural variant rather than a
        /// mapping error; floors the insert-size likelihood when scoring
        /// paired-end candidates.
        #[arg(long, default_value_t = 55.0)]
        svprior: f64,
        /// Input quality encoding is Solexa/Illumina 1.0-1.3 (Phred+64) instead of Sanger (Phred+33).
        #[arg(long, default_value_t = false)]
        phred64: bool,
        /// Input file format: "fastq" (default) or "fasta".
        #[arg(long, default_value = "fastq")]
        inputformat: String,
        /// Set read-group tags (SAM format), e.g. "ID:rg1,SM:sample1,PL:illumina".
        /// Requires at least an ID tag; adds an @RG header line and tags
        /// every record with RG:Z:<id>.
        #[arg(long)]
        readgroup: Option<String>,
        /// Output format: "sam" (default) or "bam" (BGZF-compressed binary,
        /// via the `noodles` crates).
        #[arg(short = 'f', long, default_value = "sam")]
        outputformat: String,
        /// Write output to FILE instead of stdout. Required for
        /// `--outputformat=bam` (binary output shouldn't go to a terminal).
        #[arg(short, long)]
        output: Option<String>,
        /// Only report reads whose ID starts with this prefix.
        #[arg(long)]
        labelfilter: Option<String>,
        /// Accept pre-CASAVA-1.8 `/1`/`/2`-suffixed read IDs (modern CASAVA
        /// 1.8+ headers are always handled correctly regardless of this
        /// flag; see `io::fastx::ReadPreprocessing`).
        #[arg(long, default_value_t = false)]
        casava8: bool,
        /// Trim this adapter sequence from each read's 3' end before
        /// mapping (sequencing-error-tolerant: up to 1 mismatch per 10bp of
        /// overlap).
        #[arg(long)]
        adapter_strip: Option<String>,
        /// Coordinate-sort the output and write a BAM index (`{output}.bai`)
        /// alongside it. Requires `--outputformat=bam`.
        #[arg(long, default_value_t = false)]
        index: bool,
        /// Read file(s) to map: one file for single-ended, two (mate1 mate2) for paired-end.
        #[arg(short = 'M', num_args = 1..=2)]
        reads: Vec<String>,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::BuildGenome { species, assembly, output, files } => {
            println!("Building genome index for {} {}...", species, assembly);
            let files: Vec<&String> = files.iter().collect();
            build_genome_index(&species, &assembly, &output, &files)?;
            println!("Done. Index written to {}.stidx", output);
        }
        Commands::BuildHash { genome, output } => {
            println!("Building 15-mer hash table...");
            build_hash_table(&genome, &output)?;
            println!("Done. Hash table written to {}.sthash", output);
        }
        Commands::Map { genome, hash, substitution_rate, threads, gapopen, gapextend, svprior, phred64, inputformat, readgroup, outputformat, output, labelfilter, casava8, adapter_strip, index, reads } => {
            let opts = MapperOptions {
                substitution_rate,
                threads,
                gap_open_phred: gapopen,
                gap_extend_phred: gapextend,
                sv_prior_phred: svprior,
                casava8,
                adapter_strip: adapter_strip.map(|s| s.into_bytes()),
                ..Default::default()
            };
            let quality_offset = if phred64 { PHRED64_OFFSET } else { PHRED33_OFFSET };
            let input_format = match inputformat.as_str() {
                "fasta" => InputFormat::Fasta,
                _ => InputFormat::Fastq,
            };
            let read_group = readgroup.as_deref().map(ReadGroup::parse).transpose().map_err(|e| anyhow!(e))?;
            let contigs: Vec<(String, usize)> = GenomeIndex::load(&genome)?
                .contigs.iter().map(|c| (c.name.clone(), c.length)).collect();
            let command_line = env::args().collect::<Vec<_>>().join(" ");
            let reads: Vec<&String> = reads.iter().collect();
            let mut records = map_reads_with_format(&genome, &hash, &reads, input_format, quality_offset, opts)?;

            if let Some(prefix) = &labelfilter {
                records.retain(|r| r.read_id.starts_with(prefix.as_str()));
            }

            if index {
                if outputformat != "bam" {
                    return Err(anyhow!("--index requires --outputformat=bam"));
                }
                if output.is_none() {
                    return Err(anyhow!("--index requires --output FILE (an index needs a real file to index)"));
                }
                sort_for_coordinate_order(&mut records, &contigs);
            }

            let sam_text = if index {
                format_sam_records_sorted(&records, &contigs, &command_line, read_group.as_ref(), true)
            } else {
                format_sam_records(&records, &contigs, &command_line, read_group.as_ref())
            };

            match outputformat.as_str() {
                "bam" => {
                    let Some(path) = &output else {
                        return Err(anyhow!("--outputformat=bam requires --output FILE (binary output shouldn't go to a terminal)"));
                    };
                    if index {
                        write_bam_indexed(&sam_text, std::path::Path::new(path))?;
                    } else {
                        write_bam(&sam_text, File::create(path)?)?;
                    }
                }
                _ => {
                    let mut sink: Box<dyn Write> = match &output {
                        Some(path) => Box::new(File::create(path)?),
                        None => Box::new(io::stdout()),
                    };
                    sink.write_all(sam_text.as_bytes())?;
                }
            }
        }
    }
    Ok(())
}
