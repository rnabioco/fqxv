//! `fqxv` command-line interface — a thin front-end over the [`fqxv`] library.

use std::fs::File;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::Context;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Parser, Subcommand, ValueEnum};
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use tabled::builder::Builder as TableBuilder;
use tabled::settings::object::Columns;
use tabled::settings::{Alignment, Modify, Style};
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

  # Squeeze harder: top level plus read reordering (single-end order may change)
  fqxv compress reads.fastq.gz -o reads.fqxv -l 9 --order any

  # Decompress to stdout and pipe straight into an aligner
  fqxv decompress sample.fqxv | bwa mem -p ref.fa -

  # Restore the separate paired files
  fqxv decompress sample.fqxv --split sample

  # Inspect an archive without decompressing it (add --json or --tsv for scripts)
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
        /// Read-order guarantee. `preserve` (default) restores the original order
        /// on decompress. `any` allows read reordering to exploit depth redundancy
        /// for a better ratio — single-end reads may come back in a different
        /// order; paired/grouped input still round-trips in order (the mate
        /// interleaving requires it), it just costs a stored permutation.
        #[arg(long, value_enum, default_value_t = ReadOrder::Preserve, help_heading = "Advanced")]
        order: ReadOrder,
        /// With `--order any`, use the literal-rescue sequence codec: keep every
        /// contig alive and re-attach would-be literals to any contig they
        /// overlap. Smaller sequence stream on deep data, at a higher encode
        /// cost. No effect without `--order any`.
        #[arg(long, help_heading = "Advanced")]
        rescue: bool,
        /// With `--order any`, force original read order to be restored on
        /// decompress (store a permutation, code names/quality in original order).
        /// By default single-end reorder picks this automatically when it makes
        /// the archive smaller — counter-style names (e.g. an incrementing
        /// `.N N`) delta-code to almost nothing in original order, so the
        /// permutation beats a scrambled-name stream. Pass this to force it on.
        #[arg(long, help_heading = "Advanced")]
        keep_order: bool,
        /// Opt-in lossy quality binning (changes the data; default is lossless).
        #[arg(long, value_enum, default_value_t = QualityBin::Lossless, help_heading = "Advanced")]
        quality_bin: QualityBin,
        /// Sequencing platform to record in the archive. Auto-detected from the
        /// read names by default; pass to force it (e.g. for unusual name
        /// conventions the detector doesn't recognize).
        #[arg(long, value_enum, help_heading = "Advanced")]
        platform: Option<Platform>,
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
        /// Best-effort recovery of a corrupted archive: skip blocks that fail
        /// their CRC and decode the rest, reporting what was lost. Interleaved
        /// output only (plain, non-reordered archives).
        #[arg(long, conflicts_with = "split")]
        recover: bool,
    },
    /// Print `.fqxv` container metadata and per-stream sizes.
    Info {
        /// Input `.fqxv` file.
        input: PathBuf,
        /// Emit a single machine-readable TSV line instead of the human
        /// report (stable columns for the benchmark harness / scripts).
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit a JSON object instead of the human report.
        #[arg(long)]
        json: bool,
    },
    /// Verify archive integrity (CRC checks) without writing any output.
    /// Prints a table of per-check results; exits non-zero if the archive is
    /// corrupt.
    Verify {
        /// Input `.fqxv` file.
        input: PathBuf,
        /// Faster, weaker check: verify each block's stored CRC via the footer
        /// index (parallel positioned reads) instead of the whole-file digest.
        /// Skips the header, footer, and inter-block framing bytes.
        #[arg(long)]
        quick: bool,
        /// Emit tab-separated per-check rows instead of the human table.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit a JSON object instead of the human table.
        #[arg(long)]
        json: bool,
    },
}

/// Read-order guarantee exposed on the CLI. Reordering (`Any`) is a compression
/// technique, not a user-facing knob — the user picks the property they care
/// about (does my order survive?), and the library chooses the mechanism.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum ReadOrder {
    /// Original read order is restored on decompress (default).
    Preserve,
    /// Allow reordering for a better ratio; single-end order may change.
    Any,
    /// Like `any`, but if the read names are purely positional (a counter, e.g.
    /// SRA `@RUN.N N`), discard the original order entirely: renumber the reads
    /// and regenerate the names, dropping the permutation and the name stream for
    /// the smallest archive. **Reorder-lossy — reads are renumbered** (sequence
    /// and quality preserved exactly). Falls back to `any` when names aren't a
    /// counter. Single-end only.
    Shuffle,
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

/// Sequencing platform override for `compress --platform` (absent = auto-detect).
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Platform {
    /// Illumina (colon-delimited instrument names).
    Illumina,
    /// Oxford Nanopore (UUID read names).
    Nanopore,
    /// PacBio (movie/zmw read names).
    Pacbio,
    /// MGI / BGI (V/E/DP-prefixed read names).
    Mgi,
}

impl From<Platform> for fqxv::Platform {
    fn from(p: Platform) -> Self {
        match p {
            Platform::Illumina => fqxv::Platform::Illumina,
            Platform::Nanopore => fqxv::Platform::Nanopore,
            Platform::Pacbio => fqxv::Platform::PacBio,
            Platform::Mgi => fqxv::Platform::MgiBgi,
        }
    }
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

/// Resolve the `--threads` budget to a concrete worker count: 0 means all
/// available cores, and any explicit request is clamped to what physically
/// exists so we never oversubscribe. Mirrors the library's compression pool
/// sizing so decode and compress use the same budget.
fn resolve_threads(threads: usize) -> usize {
    let available = std::thread::available_parallelism().map_or(1, |n| n.get());
    if threads == 0 {
        available
    } else {
        threads.min(available)
    }
}

/// Size rayon's global thread pool to the `--threads` budget. The parallel BGZF
/// decoder ([`noodles_bgzf::io::MultithreadedReader`]) runs on this global pool,
/// so this makes `--threads` govern decode as well as compress. Must run before
/// any rayon use; a re-init error would only occur if the pool were already
/// built, which never happens this early.
fn configure_global_pool(threads: usize) -> Result<(), rayon::ThreadPoolBuildError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(resolve_threads(threads))
        .build_global()
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    init_tracing(cli.verbose, cli.quiet);
    configure_global_pool(cli.threads).context("configuring the global thread pool")?;
    match cli.command {
        Command::Compress {
            inputs,
            output,
            level,
            interleaved,
            order,
            rescue,
            keep_order,
            quality_bin,
            platform,
        } => {
            if inputs.is_empty() {
                anyhow::bail!("at least one input FASTQ is required");
            }
            let interleaved = interleaved.filter(|_| inputs.len() == 1);
            // `--order any`/`shuffle` turn on reordering; the library forces the
            // permutation (keep_order) back on for grouped input, so paired
            // archives still round-trip in order regardless. `shuffle` additionally
            // opts into name regeneration (reorder-lossy) for single-end input.
            let reorders = order != ReadOrder::Preserve;
            let params = fqxv::Params {
                seq_order: level_to_order(level),
                block_reads: level_to_block(level),
                quality_binning: quality_bin.into(),
                reorder: reorders,
                keep_order: keep_order && reorders,
                rescue: rescue && reorders,
                regenerate_names: order == ReadOrder::Shuffle,
                threads: cli.threads,
                platform: platform.map(Into::into),
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
            recover,
        } => {
            let open_in =
                || File::open(&input).with_context(|| format!("opening input {}", input.display()));
            let t0 = Instant::now();
            if recover {
                let rec = fqxv::decompress_recover(
                    open_in()?,
                    open_output(output.as_deref())?,
                    cli.threads,
                )?;
                if rec.blocks_skipped > 0 {
                    // User-facing summary on stderr so it shows even when the
                    // decoded FASTQ is piped from stdout.
                    eprintln!(
                        "warning: recovered {} block(s), skipped {} corrupt block(s) — {} read(s) lost",
                        rec.blocks_recovered, rec.blocks_skipped, rec.reads_lost
                    );
                }
                info!(
                    reads = rec.stats.reads,
                    blocks = rec.blocks_recovered,
                    skipped = rec.blocks_skipped,
                    reads_lost = rec.reads_lost,
                    secs = format_args!("{:.1}", t0.elapsed().as_secs_f64()),
                    "recovered"
                );
            } else {
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
        }
        Command::Info { input, tsv, json } => print_info(&input, tsv, json)?,
        Command::Verify {
            input,
            quick,
            tsv,
            json,
        } => print_verify(&input, quick, tsv, json)?,
    }
    Ok(())
}

/// A per-stream size row in the machine-readable [`InfoReport`].
#[derive(Debug, Serialize)]
struct StreamJson {
    bytes: u64,
    /// Share of the three compressed streams, in percent (0 when empty).
    pct: f64,
}

/// JSON shape emitted by `fqxv info --json`. Mirrors the human report and adds
/// derived fields (labels, percentages, bytes/read) so consumers don't have to
/// recompute them. Unlike the TSV line this is free to grow.
#[derive(Debug, Serialize)]
struct InfoReport {
    file: String,
    file_size: u64,
    /// Human-readable file size (e.g. `"1.06 MB"`).
    file_size_human: String,
    platform: String,
    reads: u64,
    /// Present only for grouped/paired archives (`group_size > 1`).
    #[serde(skip_serializing_if = "Option::is_none")]
    spots: Option<u64>,
    blocks: u64,
    layout: String,
    group_size: u8,
    sequence_order: u8,
    quality: String,
    quality_binning: u8,
    reordered: bool,
    read_order_preserved: bool,
    plus_normalized: bool,
    streams: InfoStreams,
    /// Compressed stream bytes divided by read count (null when there are none).
    #[serde(skip_serializing_if = "Option::is_none")]
    bytes_per_read: Option<f64>,
}

/// The `streams` object nested in [`InfoReport`].
#[derive(Debug, Serialize)]
struct InfoStreams {
    names: StreamJson,
    sequence: StreamJson,
    quality: StreamJson,
    total: StreamJson,
}

fn print_info(path: &Path, tsv: bool, json: bool) -> anyhow::Result<()> {
    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let info = fqxv::inspect(
        File::open(path).with_context(|| format!("opening input {}", path.display()))?,
    )?;
    let total = info.names_bytes + info.seq_bytes + info.qual_bytes;

    if tsv {
        // Stable, tab-separated columns for scripts. Keep field order fixed;
        // append new fields at the end so existing parsers don't break.
        println!(
            "file_size\treads\tblocks\tgroup_size\tseq_order\tquality_binning\treordered\tnames_bytes\tseq_bytes\tqual_bytes\tplatform"
        );
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            info.platform.token(),
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
    let spots = (info.group_size > 1).then(|| info.reads / info.group_size.max(1) as u64);
    let bytes_per_read = (info.reads > 0).then(|| total as f64 / info.reads as f64);

    if json {
        let report = InfoReport {
            file: path.display().to_string(),
            file_size,
            file_size_human: human_bytes(file_size),
            platform: info.platform.token().to_string(),
            reads: info.reads,
            spots,
            blocks: info.blocks,
            layout,
            group_size: info.group_size,
            sequence_order: info.seq_order,
            quality: quality.to_string(),
            quality_binning: info.quality_binning,
            reordered: info.reordered,
            read_order_preserved: info.keep_order,
            plus_normalized: info.plus_normalized,
            streams: InfoStreams {
                names: StreamJson {
                    bytes: info.names_bytes,
                    pct: pct(info.names_bytes),
                },
                sequence: StreamJson {
                    bytes: info.seq_bytes,
                    pct: pct(info.seq_bytes),
                },
                quality: StreamJson {
                    bytes: info.qual_bytes,
                    pct: pct(info.qual_bytes),
                },
                total: StreamJson {
                    bytes: total,
                    pct: pct(total),
                },
            },
            bytes_per_read,
        };
        println!("{}", serde_json::to_string_pretty(&report)?);
        return Ok(());
    }

    // Human report: a metadata table plus a per-stream size table.
    let read_order = if info.keep_order {
        "preserved (permutation stored)"
    } else if info.regenerated_names {
        "discarded (reads renumbered, names regenerated)"
    } else {
        "clustered (not preserved)"
    };
    let mut meta = TableBuilder::default();
    let mut meta_row = |k: &str, v: String| meta.push_record([k.to_string(), v]);
    meta_row("property", "value".to_string());
    meta_row(
        "layout",
        format!("{layout} (group size {})", info.group_size),
    );
    meta_row("reads", group_digits(info.reads));
    if let Some(spots) = spots {
        meta_row("spots", group_digits(spots));
    }
    meta_row("blocks", group_digits(info.blocks));
    meta_row("platform", info.platform.label().to_string());
    meta_row("sequence order", info.seq_order.to_string());
    meta_row("quality", quality.to_string());
    meta_row("reordered", bool_word(info.reordered, "yes", "no"));
    if info.reordered {
        meta_row("read order", read_order.to_string());
    }
    meta_row(
        "plus line",
        bool_word(info.plus_normalized, "normalized", "verbatim"),
    );
    meta_row(
        "file size",
        format!(
            "{} ({} bytes)",
            human_bytes(file_size),
            group_digits(file_size)
        ),
    );
    let mut meta = meta.build();
    meta.with(Style::rounded());

    let mut streams = TableBuilder::default();
    streams.push_record([
        "stream".to_string(),
        "bytes".to_string(),
        "share".to_string(),
    ]);
    for (label, bytes) in [
        ("names", info.names_bytes),
        ("sequence", info.seq_bytes),
        ("quality", info.qual_bytes),
        ("total", total),
    ] {
        streams.push_record([
            label.to_string(),
            group_digits(bytes),
            format!("{:.1}%", pct(bytes)),
        ]);
    }
    let mut streams = streams.build();
    streams
        .with(Style::rounded())
        .with(Modify::new(Columns::new(1..)).with(Alignment::right()));

    println!("{}", path.display());
    println!("{meta}");
    println!("{streams}");
    if let Some(bpr) = bytes_per_read {
        println!("{bpr:.2} bytes/read");
    }
    Ok(())
}

/// One check row in the machine-readable [`VerifyJson`].
#[derive(Debug, Serialize)]
struct VerifyCheckJson {
    name: String,
    ok: bool,
    detail: String,
}

/// JSON shape emitted by `fqxv verify --json`.
#[derive(Debug, Serialize)]
struct VerifyJson {
    file: String,
    passed: bool,
    checks: Vec<VerifyCheckJson>,
    /// Indices of blocks whose CRC failed (omitted when intact).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    failed_blocks: Vec<u64>,
}

/// Run `fqxv::verify_report` and render it as a table (default), TSV, or JSON.
/// Exits the process non-zero when the archive is corrupt so scripts can branch
/// on `$?` regardless of the chosen format.
fn print_verify(path: &Path, quick: bool, tsv: bool, json: bool) -> anyhow::Result<()> {
    let file = File::open(path).with_context(|| format!("opening input {}", path.display()))?;
    let report = fqxv::verify_report(&file, quick)
        .with_context(|| format!("{} is not a readable fqxv archive", path.display()))?;
    let passed = report.passed();

    if json {
        let out = VerifyJson {
            file: path.display().to_string(),
            passed,
            checks: report
                .checks
                .iter()
                .map(|c| VerifyCheckJson {
                    name: c.name.clone(),
                    ok: c.ok,
                    detail: c.detail.clone(),
                })
                .collect(),
            failed_blocks: report.failed_blocks.clone(),
        };
        println!("{}", serde_json::to_string_pretty(&out)?);
    } else if tsv {
        println!("check\tresult\tdetail");
        for c in &report.checks {
            println!(
                "{}\t{}\t{}",
                c.name,
                if c.ok { "ok" } else { "fail" },
                c.detail
            );
        }
    } else {
        let mut table = TableBuilder::default();
        table.push_record([
            "check".to_string(),
            "result".to_string(),
            "detail".to_string(),
        ]);
        for c in &report.checks {
            table.push_record([
                c.name.clone(),
                if c.ok { "ok" } else { "FAIL" }.to_string(),
                c.detail.clone(),
            ]);
        }
        let mut table = table.build();
        table.with(Style::rounded());
        println!("{}", path.display());
        println!("{table}");
        println!(
            "{}: {}",
            path.display(),
            if passed { "OK" } else { "CORRUPT" }
        );
    }

    if !passed {
        std::process::exit(1);
    }
    Ok(())
}

/// Pick one of two words by a flag and own it (small helper so the metadata
/// table rows are all `String`).
fn bool_word(flag: bool, yes: &str, no: &str) -> String {
    if flag { yes } else { no }.to_string()
}

/// Format a byte count with a binary-scaled unit (KB/MB/GB/TB/PB, 1024-based like
/// `du -h`). Bytes stay integral; larger units show two decimals.
fn human_bytes(n: u64) -> String {
    const UNITS: [&str; 6] = ["bytes", "KB", "MB", "GB", "TB", "PB"];
    if n < 1024 {
        return format!("{n} bytes");
    }
    let mut value = n as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit < UNITS.len() - 1 {
        value /= 1024.0;
        unit += 1;
    }
    format!("{value:.2} {}", UNITS[unit])
}

/// Format an integer with `,` thousands separators (e.g. `1234567` -> `1,234,567`)
/// so large byte and read counts stay legible in the human report.
fn group_digits(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    let first = digits.len() % 3;
    for (i, ch) in digits.chars().enumerate() {
        if i != 0 && i >= first && (i - first).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// True if `hdr` begins with a BGZF block header: gzip magic (`1f 8b`), deflate
/// method, the `FEXTRA` flag, and the mandatory `BC` extra subfield (`SI1='B'`,
/// `SI2='C'`) that BGZF places first in the header. BGZF is the block-gzip
/// variant emitted by `bgzip`/`samtools`; unlike plain gzip its blocks are
/// independently inflatable, so decode can be parallelized. Matches htslib's
/// fixed-offset check (a spec-conformant BGZF header is at least 18 bytes).
fn is_bgzf(hdr: &[u8]) -> bool {
    hdr.len() >= 18
        && hdr[0] == 0x1f
        && hdr[1] == 0x8b
        && hdr[2] == 0x08 // CM = deflate
        && (hdr[3] & 0x04) != 0 // FLG.FEXTRA
        && hdr[12] == b'B'
        && hdr[13] == b'C'
}

/// Open a FASTQ input, transparently decoding gzip (detected by magic bytes).
/// A path of `-` reads from stdin, so a downloader can pipe straight in.
///
/// BGZF (block-gzip) input is inflated in parallel: its blocks are independent,
/// so [`noodles_bgzf::io::MultithreadedReader`] decodes them across rayon's
/// global thread pool, sized to `--threads` by [`configure_global_pool`]. Plain
/// gzip is a single DEFLATE stream and stays serial. The decoded byte stream —
/// and therefore the archive — is identical regardless of how the input was
/// compressed or how many decode workers run.
fn open_input(path: &Path) -> anyhow::Result<Box<dyn Read + Send>> {
    let mut src: Box<dyn Read + Send> = if path.as_os_str() == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening input {}", path.display()))?)
    };
    // Peek enough to classify: 2 bytes for gzip magic, 18 for a full BGZF header.
    let mut magic = [0u8; 18];
    let n = read_up_to(&mut src, &mut magic)?;
    let head = io::Cursor::new(magic[..n].to_vec());
    let chained = head.chain(src);
    if is_bgzf(&magic[..n]) {
        Ok(Box::new(noodles_bgzf::io::MultithreadedReader::new(
            chained,
        )))
    } else if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A real 18-byte BGZF block header (from `bgzip`): gzip magic, deflate,
    /// FEXTRA, then the `BC` subfield carrying BSIZE.
    const BGZF_HEADER: [u8; 18] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
        0x00, 0x5f, 0x00,
    ];

    #[test]
    fn detects_bgzf_header() {
        assert!(is_bgzf(&BGZF_HEADER));
    }

    #[test]
    fn plain_gzip_is_not_bgzf() {
        // gzip magic + deflate, FNAME flag (0x08) instead of FEXTRA — plain gzip.
        let hdr = [
            0x1f, 0x8b, 0x08, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00,
        ];
        assert!(!is_bgzf(&hdr));
    }

    #[test]
    fn non_gzip_and_short_inputs_are_not_bgzf() {
        assert!(!is_bgzf(b"@read\n")); // plain FASTQ
        assert!(!is_bgzf(&[0x1f, 0x8b])); // gzip magic but too short for a BGZF header
        assert!(!is_bgzf(&[])); // empty
                                // FEXTRA present but the extra subfield is not `BC`.
        let mut hdr = BGZF_HEADER;
        hdr[12] = b'X';
        assert!(!is_bgzf(&hdr));
    }

    #[test]
    fn resolve_threads_zero_is_all_cores_and_explicit_is_clamped() {
        let available = std::thread::available_parallelism().map_or(1, |n| n.get());
        assert_eq!(resolve_threads(0), available); // 0 = all cores
        assert_eq!(resolve_threads(1), 1);
        assert_eq!(resolve_threads(usize::MAX), available); // never oversubscribe
    }
}
