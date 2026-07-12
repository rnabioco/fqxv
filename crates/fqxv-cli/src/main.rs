//! `fqxv` command-line interface — a thin front-end over the [`fqxv`] library.

use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

/// Reference-free FASTQ archiver for short-read data.
#[derive(Debug, Parser)]
#[command(name = "fqxv", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Number of worker threads (0 = all available cores).
    #[arg(long, global = true, default_value_t = 0)]
    threads: usize,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Compress a FASTQ file to `.fqxv`.
    Compress {
        /// Input FASTQ (optionally gzipped).
        input: PathBuf,
        /// Output `.fqxv` path.
        #[arg(short, long)]
        output: PathBuf,
        /// Compression effort level.
        #[arg(short, long, default_value_t = 5)]
        level: u8,
        /// Reorder reads to exploit depth redundancy (may not preserve order).
        #[arg(long)]
        reorder: bool,
        /// Preserve exact input read order even when reordering.
        #[arg(long)]
        keep_order: bool,
        /// Opt-in lossy quality binning.
        #[arg(long, value_enum, default_value_t = QualityBin::Lossless)]
        quality_bin: QualityBin,
    },
    /// Decompress a `.fqxv` file back to FASTQ.
    Decompress {
        /// Input `.fqxv` file.
        input: PathBuf,
        /// Output FASTQ path.
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Print the container header and per-stream statistics.
    Inspect {
        /// Input `.fqxv` file.
        input: PathBuf,
    },
}

/// Lossy quality quantization choices exposed on the CLI.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum QualityBin {
    /// Fully lossless (default).
    Lossless,
    /// Illumina 8-level binning (lossy).
    Bin8,
    /// Illumina 4-level binning (lossy).
    Bin4,
    /// 2-level binning (lossy).
    Bin2,
}

fn main() -> anyhow::Result<()> {
    // Container/format version is owned by the library crate.
    let _ = fqxv::FORMAT_VERSION;

    let cli = Cli::parse();
    match cli.command {
        Command::Compress { .. } => {
            anyhow::bail!("`fqxv compress` is not implemented yet (M5)");
        }
        Command::Decompress { .. } => {
            anyhow::bail!("`fqxv decompress` is not implemented yet (M5)");
        }
        Command::Inspect { .. } => {
            anyhow::bail!("`fqxv inspect` is not implemented yet (M5)");
        }
    }
}
