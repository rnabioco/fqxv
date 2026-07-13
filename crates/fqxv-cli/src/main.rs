//! `fqxv` command-line interface — a thin front-end over the [`fqxv`] library.

use std::fs::File;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::Context;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Parser, Subcommand, ValueEnum};
use flate2::read::MultiGzDecoder;
use tracing::info;
use tracing_subscriber::{fmt, EnvFilter};

/// Terminal color scheme for `--help` (shared with the rnabioco tooling look):
/// yellow-bold headings, green-bold literals, cyan value placeholders.
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Yellow.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Yellow.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default());

/// Worked examples, appended under the top-level `--help`.
const EXAMPLES: &str = "\
Examples:
  # Compress a single-end FASTQ (gzip is auto-detected by magic bytes)
  fqxv compress reads.fastq.gz -o reads.fqxv

  # Paired-end: mates are interleaved per spot into one archive
  fqxv compress R1.fastq.gz R2.fastq.gz -o sample.fqxv

  # Squeeze harder: top level plus read reordering (does not preserve order)
  fqxv compress reads.fastq.gz -o reads.fqxv -l 9 --reorder

  # Decompress to stdout and pipe straight into an aligner
  fqxv decompress sample.fqxv | bwa mem -p ref.fa -

  # Restore the separate paired files
  fqxv decompress sample.fqxv --split sample

  # Inspect an archive without decompressing it
  fqxv info sample.fqxv

  # Download from SRA with sracha (rnabioco/sracha-rs) and archive in one pass,
  # nothing hitting disk: -Z streams interleaved FASTQ to stdout, '-' reads it,
  # and paired data is auto-detected as interleaved (override with --interleaved N)
  sracha get -Z --split interleaved SRR2584863 | fqxv compress - -o SRR2584863.fqxv
";

/// Reference-free FASTQ archiver for short-read data.
#[derive(Debug, Parser)]
#[command(name = "fqxv", version, about, long_about = None, styles = STYLES, after_help = EXAMPLES)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Number of worker threads, capped at available cores (0 = all cores).
    #[arg(long, global = true, default_value_t = 16)]
    threads: usize,

    /// Increase log verbosity: -v debug, -vv trace, -vvv trace with targets,
    /// thread ids, and span timing. Overridden by RUST_LOG if set.
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Silence all output except warnings and errors (suppresses the summary).
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,
}

/// Install the `tracing` subscriber. Verbosity comes from `-v/-vv/-vvv` and
/// `--quiet`; a set `RUST_LOG` overrides the computed level entirely. All output
/// goes to **stderr** so decompressed FASTQ on stdout stays uncontaminated.
fn init_tracing(verbose: u8, quiet: bool) {
    let level = if quiet {
        "warn"
    } else {
        match verbose {
            0 => "info",
            1 => "debug",
            _ => "trace",
        }
    };
    // Scope to our own crates so dependency noise (rayon, noodles) stays out even
    // at trace; RUST_LOG, when set, takes over completely.
    let filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new(format!("fqxv={level},fqxv_cli={level}")));
    let builder = fmt()
        .with_writer(io::stderr)
        .with_env_filter(filter)
        .with_level(true);
    // -vvv adds the noisy-but-useful detail; lower levels stay terse.
    if verbose >= 3 {
        builder
            .with_target(true)
            .with_thread_ids(true)
            .with_span_events(fmt::format::FmtSpan::CLOSE)
            .init();
    } else {
        builder.with_target(false).without_time().init();
    }
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Compress FASTQ to `.fqxv`. Give multiple inputs to interleave per-spot
    /// files (paired mates, or single-cell R1/R2/I1/I2) into one archive.
    Compress {
        /// Input FASTQ file(s), plain or gzipped; `-` reads one stream from
        /// stdin. One file = single-end, 2 = paired, 3-4 = single-cell; multiple
        /// files are interleaved per spot. Order is preserved for `--split`.
        #[arg(num_args = 1..)]
        inputs: Vec<PathBuf>,
        /// Output `.fqxv` path.
        #[arg(short, long)]
        output: PathBuf,
        /// Compression effort level (1-9); higher raises the sequence order.
        #[arg(short, long, default_value_t = 5)]
        level: u8,

        /// Interleaving of a single input, in members per spot. Auto-detected
        /// from read names by default; pass to force (1 = single-end, 2 = paired
        /// as from `sracha get -Z`). Ignored with multiple inputs.
        #[arg(long, value_name = "N", help_heading = "Advanced")]
        interleaved: Option<u8>,
        /// Reorder reads to exploit depth redundancy (may not preserve order).
        #[arg(long, help_heading = "Advanced")]
        reorder: bool,
        /// With `--reorder`, store a permutation so the original order is
        /// restored on decompress.
        #[arg(long, requires = "reorder", help_heading = "Advanced")]
        keep_order: bool,
        /// Opt-in lossy quality binning (changes the data; default is lossless).
        #[arg(long, value_enum, default_value_t = QualityBin::Lossless, help_heading = "Advanced")]
        quality_bin: QualityBin,
    },
    /// Decompress a `.fqxv` file to FASTQ.
    Decompress {
        /// Input `.fqxv` file.
        input: PathBuf,
        /// Interleaved FASTQ output; omit or use `-` for stdout (pipe to an
        /// aligner: `fqxv decompress x.fqxv | bwa mem -p ref -`).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Restore separate per-spot files as `<prefix>_1.fastq … _G.fastq`.
        #[arg(long, conflicts_with = "output")]
        split: Option<PathBuf>,
    },
    /// Print `.fqxv` container metadata and per-stream sizes.
    Info {
        /// Input `.fqxv` file.
        input: PathBuf,
        /// Emit a single machine-readable TSV line instead of the human
        /// report (stable columns for the benchmark harness / scripts).
        #[arg(long)]
        tsv: bool,
    },
}

/// Lossy quality quantization choices exposed on the CLI.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum QualityBin {
    /// Fully lossless (default).
    Lossless,
    /// Illumina standard 8-level binning (lossy).
    Bin8,
    /// Illumina documented 4-level binning (NovaSeq X / RTA4; lossy).
    Bin4,
    /// Custom 2-level binning (lossy).
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

/// Map a 1-9 effort level to reads per block. Larger blocks train the sequence
/// model on more reads (better ratio) at the cost of parallelism and memory.
fn level_to_block(level: u8) -> usize {
    match level {
        0..=2 => 128 << 10,
        3..=4 => 256 << 10,
        5..=6 => 1 << 20,
        7..=8 => 2 << 20,
        _ => 4 << 20,
    }
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    init_tracing(cli.verbose, cli.quiet);
    match cli.command {
        Command::Compress {
            inputs,
            output,
            level,
            interleaved,
            reorder,
            keep_order,
            quality_bin,
        } => {
            if inputs.is_empty() {
                anyhow::bail!("at least one input FASTQ is required");
            }
            if reorder && inputs.len() > 1 {
                anyhow::bail!("--reorder is single-end only (would break spot grouping)");
            }
            let interleaved = interleaved.filter(|_| inputs.len() == 1);
            if reorder && interleaved.is_some_and(|g| g > 1) {
                anyhow::bail!("--reorder cannot combine with --interleaved (would break spots)");
            }
            let params = fqxv::Params {
                seq_order: level_to_order(level),
                block_reads: level_to_block(level),
                quality_binning: quality_bin.into(),
                reorder,
                keep_order,
                threads: cli.threads,
            };
            let in_size: u64 = inputs
                .iter()
                .filter_map(|p| std::fs::metadata(p).ok())
                .map(|m| m.len())
                .sum();
            let out = File::create(&output)
                .with_context(|| format!("creating output {}", output.display()))?;

            let t0 = Instant::now();
            let stats = if inputs.len() == 1 {
                match interleaved {
                    None => fqxv::compress_auto(open_input(&inputs[0])?, out, params)?,
                    Some(g) if g > 1 => {
                        fqxv::compress_interleaved(open_input(&inputs[0])?, out, params, g)?
                    }
                    Some(_) => fqxv::compress(open_input(&inputs[0])?, out, params)?,
                }
            } else {
                let readers: Vec<Box<dyn Read + Send>> = inputs
                    .iter()
                    .map(|p| open_input(p))
                    .collect::<anyhow::Result<_>>()?;
                fqxv::compress_multi(readers, out, params)?
            };
            let secs = t0.elapsed().as_secs_f64();
            // Ratio needs the input size; a stdin stream (`-`) has none, so omit it.
            let ratio = if in_size > 0 && stats.out_bytes > 0 {
                format!(" ({:.2}x)", in_size as f64 / stats.out_bytes as f64)
            } else {
                String::new()
            };
            let layout = match stats.group_size {
                0 | 1 => "single-end".to_string(),
                2 => "paired".to_string(),
                g => format!("grouped x{g}"),
            };
            info!(
                reads = stats.reads,
                inputs = inputs.len(),
                blocks = stats.blocks,
                out_bytes = stats.out_bytes,
                secs = format_args!("{secs:.1}"),
                "compressed {layout}{ratio}"
            );
        }
        Command::Decompress {
            input,
            output,
            split,
        } => {
            let open_in =
                || File::open(&input).with_context(|| format!("opening input {}", input.display()));
            let t0 = Instant::now();
            let stats = if let Some(prefix) = split {
                let g = fqxv::peek(open_in()?)?.group_size as usize;
                let mut files: Vec<File> = (1..=g)
                    .map(|i| {
                        let path = format!("{}_{}.fastq", prefix.display(), i);
                        File::create(&path).with_context(|| format!("creating output {path}"))
                    })
                    .collect::<anyhow::Result<_>>()?;
                fqxv::decompress_split(open_in()?, &mut files, cli.threads)?
            } else {
                fqxv::decompress(open_in()?, open_output(output.as_deref())?, cli.threads)?
            };
            info!(
                reads = stats.reads,
                blocks = stats.blocks,
                secs = format_args!("{:.1}", t0.elapsed().as_secs_f64()),
                "decompressed"
            );
        }
        Command::Info { input, tsv } => print_info(&input, tsv)?,
    }
    Ok(())
}

fn print_info(path: &Path, tsv: bool) -> anyhow::Result<()> {
    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let info = fqxv::inspect(
        File::open(path).with_context(|| format!("opening input {}", path.display()))?,
    )?;
    let total = info.names_bytes + info.seq_bytes + info.qual_bytes;
    if tsv {
        // Stable, tab-separated columns for scripts. Keep field order fixed;
        // append new fields at the end so existing parsers don't break.
        println!(
            "file_size\treads\tblocks\tgroup_size\tseq_order\tquality_binning\treordered\tnames_bytes\tseq_bytes\tqual_bytes"
        );
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            file_size,
            info.reads,
            info.blocks,
            info.group_size,
            info.seq_order,
            info.quality_binning,
            info.reordered as u8,
            info.names_bytes,
            info.seq_bytes,
            info.qual_bytes,
        );
        return Ok(());
    }
    let pct = |x: u64| {
        if total > 0 {
            100.0 * x as f64 / total as f64
        } else {
            0.0
        }
    };
    let layout = match info.group_size {
        1 => "single-end".to_string(),
        2 => "paired".to_string(),
        g => format!("grouped x{g} (single-cell)"),
    };
    let quality = match info.quality_binning {
        0 => "lossless",
        1 => "lossy (Illumina 8-bin)",
        2 => "lossy (Illumina 4-bin, RTA4)",
        3 => "lossy (2-bin, custom)",
        _ => "unknown",
    };

    println!("{}", path.display());
    println!("  layout         {layout} (group size {})", info.group_size);
    println!("  reads          {}", info.reads);
    if info.group_size > 1 {
        println!(
            "  spots          {}",
            info.reads / info.group_size.max(1) as u64
        );
    }
    println!("  blocks         {}", info.blocks);
    println!("  sequence order {}", info.seq_order);
    println!("  quality        {quality}");
    println!(
        "  reordered      {}",
        if info.reordered { "yes" } else { "no" }
    );
    println!(
        "  plus line      {}",
        if info.plus_normalized {
            "normalized"
        } else {
            "verbatim"
        }
    );
    println!("  file size      {file_size} bytes");
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
            "  streams total  {total} bytes ({:.2} bytes/read)",
            total as f64 / info.reads as f64
        );
    }
    Ok(())
}

/// Open a FASTQ input, transparently decoding gzip (detected by magic bytes).
/// A path of `-` reads from stdin, so a downloader can pipe straight in.
fn open_input(path: &Path) -> anyhow::Result<Box<dyn Read + Send>> {
    let mut src: Box<dyn Read + Send> = if path.as_os_str() == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening input {}", path.display()))?)
    };
    let mut magic = [0u8; 2];
    let n = read_up_to(&mut src, &mut magic)?;
    let head = io::Cursor::new(magic[..n].to_vec());
    let chained = head.chain(src);
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(MultiGzDecoder::new(chained)))
    } else {
        Ok(Box::new(chained))
    }
}

/// Open an output sink: a file, or stdout when the path is absent or `-`.
fn open_output(path: Option<&Path>) -> anyhow::Result<Box<dyn Write>> {
    match path {
        Some(p) if p.as_os_str() != "-" => {
            Ok(Box::new(File::create(p).with_context(|| {
                format!("creating output {}", p.display())
            })?))
        }
        _ => Ok(Box::new(io::stdout())),
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
