//! `fqxv` command-line interface — a thin front-end over the [`fqxv`] library.

use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Parser, Subcommand, ValueEnum};
use flate2::read::MultiGzDecoder;

/// Terminal color scheme for `--help` (shared with the rnabioco tooling look):
/// yellow-bold headings, green-bold literals, cyan value placeholders.
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Yellow.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Yellow.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default());

/// Reference-free FASTQ archiver for short-read data.
#[derive(Debug, Parser)]
#[command(name = "fqxv", version, about, long_about = None, styles = STYLES)]
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
        /// Input FASTQ (plain or gzipped; gzip is detected automatically).
        input: PathBuf,
        /// Output `.fqxv` path.
        #[arg(short, long)]
        output: PathBuf,
        /// Compression effort level (1-9); higher raises the sequence order.
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

impl From<QualityBin> for fqxv::QualityBinning {
    fn from(b: QualityBin) -> Self {
        match b {
            QualityBin::Lossless => fqxv::QualityBinning::Lossless,
            QualityBin::Bin8 => fqxv::QualityBinning::Bin8,
            QualityBin::Bin4 => fqxv::QualityBinning::Bin4,
            QualityBin::Bin2 => fqxv::QualityBinning::Bin2,
        }
    }
}

/// Map a 1-9 effort level to a sequence context order (higher = better ratio).
fn level_to_order(level: u8) -> u8 {
    (level as usize + 6).clamp(1, 11) as u8
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Compress {
            input,
            output,
            level,
            reorder,
            keep_order,
            quality_bin,
        } => {
            if reorder {
                anyhow::bail!("--reorder (read reordering) is not implemented yet (M4)");
            }
            let _ = keep_order;
            let params = fqxv::Params {
                seq_order: level_to_order(level),
                quality_binning: quality_bin.into(),
                threads: cli.threads,
            };
            let in_size = std::fs::metadata(&input).map(|m| m.len()).unwrap_or(0);
            let reader = open_input(&input)?;
            let out = File::create(&output)?;

            let t0 = Instant::now();
            let stats = fqxv::compress(reader, out, params)?;
            let secs = t0.elapsed().as_secs_f64();

            let ratio = if stats.out_bytes > 0 {
                in_size as f64 / stats.out_bytes as f64
            } else {
                0.0
            };
            eprintln!(
                "compressed {} reads in {} block(s) -> {} bytes ({:.2}x vs {} input) in {:.1}s",
                stats.reads, stats.blocks, stats.out_bytes, ratio, in_size, secs
            );
        }
        Command::Decompress { input, output } => {
            let f = File::open(&input)?;
            let out = File::create(&output)?;
            let t0 = Instant::now();
            let stats = fqxv::decompress(f, out, cli.threads)?;
            eprintln!(
                "decompressed {} reads from {} block(s) in {:.1}s",
                stats.reads,
                stats.blocks,
                t0.elapsed().as_secs_f64()
            );
        }
        Command::Inspect { input } => {
            let info = fqxv::inspect(File::open(&input)?)?;
            let total = info.names_bytes + info.seq_bytes + info.qual_bytes;
            let pct = |x: u64| {
                if total > 0 {
                    100.0 * x as f64 / total as f64
                } else {
                    0.0
                }
            };
            println!("fqxv container");
            println!("  reads          {}", info.reads);
            println!("  blocks         {}", info.blocks);
            println!("  sequence order {}", info.seq_order);
            println!("  quality bin    {}", info.quality_binning);
            println!(
                "  plus line      {}",
                if info.plus_normalized {
                    "normalized"
                } else {
                    "verbatim"
                }
            );
            println!(
                "  names  {:>12} bytes ({:.1}%)",
                info.names_bytes,
                pct(info.names_bytes)
            );
            println!(
                "  seq    {:>12} bytes ({:.1}%)",
                info.seq_bytes,
                pct(info.seq_bytes)
            );
            println!(
                "  qual   {:>12} bytes ({:.1}%)",
                info.qual_bytes,
                pct(info.qual_bytes)
            );
            if info.reads > 0 {
                println!(
                    "  total  {:>12} bytes ({:.2} bytes/read)",
                    total,
                    total as f64 / info.reads as f64
                );
            }
        }
    }
    Ok(())
}

/// Open a FASTQ input, transparently decoding gzip (detected by magic bytes).
fn open_input(path: &Path) -> anyhow::Result<Box<dyn Read>> {
    let mut f = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = read_up_to(&mut f, &mut magic)?;
    let head = io::Cursor::new(magic[..n].to_vec());
    let chained = head.chain(f);
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(MultiGzDecoder::new(chained)))
    } else {
        Ok(Box::new(chained))
    }
}

fn read_up_to<R: Read>(r: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..])? {
            0 => break,
            k => n += k,
        }
    }
    Ok(n)
}
