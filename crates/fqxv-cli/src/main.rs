//! `fqxv` command-line interface — a thin front-end over the [`fqxv`] library.

mod progress;
mod report;

use std::fs::File;
use std::io::{self, IsTerminal, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Instant;

use anyhow::Context;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Parser, Subcommand, ValueEnum};
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use tabled::builder::Builder as TableBuilder;
use tabled::settings::object::Columns;
use tabled::settings::{Alignment, Modify, Style};
use tracing::warn;
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
  # Compress a single-end FASTQ (gzip auto-detected; -o defaults to reads.fqxv)
  fqxv compress reads.fastq.gz

  # Paired-end: mates are interleaved per spot into one archive
  fqxv compress sample_R1.fastq.gz sample_R2.fastq.gz -o sample.fqxv

  # Squeeze harder: top level plus read reordering (single-end order may change)
  fqxv compress reads.fastq.gz -l 9 --order any

  # Stream to stdout (-Z, always raw) and pipe straight into an aligner
  fqxv decompress sample.fqxv -Z | bwa mem -p ref.fa -

  # Restore the separate mate files -> sample_R1.fastq.gz, sample_R2.fastq.gz
  fqxv decompress sample.fqxv --split sample
  # ...or plain, numbered: sample_1.fastq, sample_2.fastq
  fqxv decompress sample.fqxv --split sample --no-gzip --mate-style num

  # Inspect an archive without decompressing it (add --json or --tsv for scripts)
  fqxv info sample.fqxv
  # ...with content stats (read-length spread, GC%, quality dist); decodes fully
  fqxv info sample.fqxv --stats

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
        /// Output `.fqxv` path. Defaults to the first input's name with the
        /// FASTQ/gzip extension replaced by `.fqxv` (`reads.fastq.gz` ->
        /// `reads.fqxv`), written alongside the input. Required when the input is
        /// stdin (`-`).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Overwrite the output archive if it already exists. Without this, an
        /// existing output path is left untouched and the command errors, so a
        /// stray `compress` can't silently clobber an archive.
        #[arg(short = 'f', long)]
        force: bool,
        /// Maximum-compression preset: the deepest sequence context plus read
        /// reordering *where it helps*. Overrides `--level`/`--order`. Reordering
        /// is applied to short-read data and automatically skipped for long reads
        /// (nanopore/PacBio), where it costs ~10x the time and memory for no ratio
        /// gain — so `--max` adapts to the input rather than forcing one fixed
        /// setting. Single-end short reads may come back reordered; names,
        /// sequence, and quality are preserved exactly (still lossless).
        #[arg(long)]
        max: bool,
        /// Compression effort level (1-9); higher raises the sequence order.
        #[arg(short, long, default_value_t = 5, help_heading = "Advanced")]
        level: u8,
        /// Reads per row group (block). Overrides the block size `--level` would
        /// pick, decoupling granularity from effort: smaller groups give finer
        /// random access and more parallelism at some ratio cost (the order-k
        /// sequence model trains on fewer reads), larger groups the reverse. Useful
        /// when archiving to object storage where you want to fetch small ranges.
        /// Sequence order still follows `--level`. Ignored by the reorder path
        /// (`--order any`/`--max`), which clusters globally.
        #[arg(long, value_name = "N", help_heading = "Advanced")]
        block_reads: Option<usize>,

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
        /// With `--order any`, disable the adaptive assembly codecs (block-local
        /// literal-rescue and the whole-file global reference) and use the faster
        /// single-contig sequence codec only. By default reorder codes each block
        /// with every codec and keeps the smaller — plus, when a shared global
        /// reference nets a whole-file win, stores it once and codes reads as
        /// positions on it (never worse than block-local) — recovering reads the
        /// single-contig codec would strand as literals at a higher encode cost.
        /// No effect without `--order any`.
        #[arg(long = "no-rescue", help_heading = "Advanced")]
        no_rescue: bool,
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
        /// Estimate the compression ratio and archive size from a sample of the
        /// input, then exit without writing anything. Codes the leading reads
        /// with the real codecs at the chosen `--level`/`--quality-bin` and
        /// projects the whole-file archive. Reordering (`--order any`/`--max`) is
        /// not modeled, so with it the estimate is a conservative lower bound —
        /// the real archive comes out this size or smaller. Pass `tsv` for a
        /// two-line machine-readable table (`file`, `input_bytes`,
        /// `est_fqxv_bytes`, `ratio`) instead of the human report.
        #[arg(long, value_enum, num_args = 0..=1, default_missing_value = "human", value_name = "FORMAT")]
        estimate: Option<EstimateFormat>,
        /// After writing the archive, re-read and fully decode it to confirm it
        /// round-trips before reporting success — recommended before deleting the
        /// source FASTQ. Catches a codec or memory error that produced a
        /// CRC-valid-but-wrong archive. Roughly doubles wall time (adds a decode
        /// pass; the reorder/`--max` path is the heaviest). On failure the archive
        /// is left in place and the command exits non-zero.
        #[arg(long, conflicts_with = "estimate")]
        verify: bool,
    },
    /// Decompress a `.fqxv` file to FASTQ.
    Decompress {
        /// Input `.fqxv` file.
        input: PathBuf,
        /// Interleaved FASTQ output file. A `.gz` extension writes block-gzip
        /// (BGZF); any other extension writes plain FASTQ. Use `-Z/--stdout` (or
        /// `-o -`) to stream to stdout instead.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Stream interleaved FASTQ (always raw) to stdout, e.g. to pipe into an
        /// aligner: `fqxv decompress x.fqxv -Z | bwa mem -p ref -`. Required to
        /// write to stdout — a bare `decompress` with no `-o`/`--split`/`-Z` errors
        /// rather than flooding the terminal with reads.
        #[arg(short = 'Z', long, conflicts_with_all = ["output", "split"])]
        stdout: bool,
        /// Restore separate per-spot files as `<prefix>_R1.fastq.gz …
        /// _R<G>.fastq.gz` (block-gzip by default; see `--mate-style`/`--no-gzip`).
        #[arg(long, conflicts_with = "output")]
        split: Option<PathBuf>,
        /// Labeling for `--split` outputs: `r` gives `_R1`,`_R2`,… (Illumina
        /// convention, the default); `num` gives `_1`,`_2`,….
        #[arg(long, value_enum, default_value_t = MateStyle::R, requires = "split")]
        mate_style: MateStyle,
        /// Write plain `.fastq` for `--split` instead of the default block-gzip
        /// `.fastq.gz`. (For `-o FILE`, compression follows the file extension.)
        #[arg(long, requires = "split")]
        no_gzip: bool,
        /// Best-effort recovery of a corrupted archive: skip blocks that fail
        /// their CRC and decode the rest, reporting what was lost. Interleaved
        /// output only (plain, non-reordered archives).
        #[arg(long, conflicts_with = "split")]
        recover: bool,
        /// Overwrite output FASTQ file(s) if they already exist. Without this, an
        /// existing output (`-o FILE` or any `--split` mate file) is left
        /// untouched and the command errors before decoding. Ignored when writing
        /// to stdout (`-Z`/`-o -`).
        #[arg(short = 'f', long)]
        force: bool,
    },
    /// Print `.fqxv` container metadata and per-stream sizes. Give several files
    /// or a directory (scanned recursively for `*.fqxv`) to report a batch — one
    /// TSV/JSON entry per archive, keyed by filename.
    Info {
        /// Input `.fqxv` file(s), or a directory scanned recursively for `*.fqxv`.
        #[arg(num_args = 1..)]
        inputs: Vec<PathBuf>,
        /// Emit machine-readable TSV instead of the human report. A single file
        /// keeps the stable one-row columns (benchmark harness / scripts); a batch
        /// prepends a `file` column and prints one row per archive.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit a JSON object instead of the human report.
        #[arg(long)]
        json: bool,
        /// Also report content statistics — read-length spread, base
        /// composition, GC%, and the quality distribution. This decodes the whole
        /// archive (unlike the default metadata-only report, which just reads the
        /// header and footer index), so it costs a full decompress.
        #[arg(long, short = 's')]
        stats: bool,
    },
    /// Verify archive integrity (CRC checks) without writing any output.
    /// Prints a table of per-check results; exits non-zero if any archive is
    /// corrupt. Give several files or a directory (scanned recursively for
    /// `*.fqxv`) to verify a batch.
    Verify {
        /// Input `.fqxv` file(s), or a directory scanned recursively for `*.fqxv`.
        #[arg(num_args = 1..)]
        inputs: Vec<PathBuf>,
        /// Faster, weaker check: verify each block's stored CRC via the footer
        /// index (parallel positioned reads) instead of the whole-file digest.
        /// Skips the header, footer, and inter-block framing bytes.
        #[arg(long)]
        quick: bool,
        /// Emit tab-separated per-check rows instead of the human table. A batch
        /// prepends a `file` column so rows from different archives stay distinct.
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
    /// Discard the original read order AND names entirely for the smallest
    /// archive: renumber the reads and regenerate the names, dropping both the
    /// order permutation and the name stream. Names that are a positional counter
    /// (e.g. SRA `@RUN.N N`) are reproduced exactly; otherwise the reads are
    /// renumbered with a fresh 1..n counter. **Reorder-lossy — original names and
    /// order are not recoverable** (sequence and quality preserved exactly).
    /// Single-end only.
    Shuffle,
}

/// Output format for `compress --estimate`. Bare `--estimate` selects `Human`
/// (via `default_missing_value`); `--estimate tsv` selects the scriptable table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum)]
enum EstimateFormat {
    /// Human-readable projection report (default).
    #[default]
    Human,
    /// Two lines — a header then one data row: `file`, `input_bytes`,
    /// `est_fqxv_bytes`, `ratio` (current on-disk size / estimated archive size).
    Tsv,
}

/// Labeling for the per-mate files produced by `decompress --split`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum MateStyle {
    /// `_R1`, `_R2`, … — the Illumina convention (default).
    R,
    /// `_1`, `_2`, … — bare member numbers.
    Num,
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
    /// Oxford Nanopore 4-level binning (lossy; CoLoRd ONT cutpoints).
    Ont,
    /// PacBio HiFi 5-level binning (lossy; CoLoRd HiFi cutpoints, Q93 kept).
    Hifi,
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
            QualityBin::Ont => fqxv::QualityBinning::BinOnt,
            QualityBin::Hifi => fqxv::QualityBinning::BinHifi,
        }
    }
}

/// Map a 1-9 effort level to a sequence context order (higher = better ratio).
fn level_to_order(level: u8) -> u8 {
    (level as usize + 6).clamp(1, 11) as u8
}

/// Map a 1-9 effort level to the hashed high-order sequence tier `(order, bits)`
/// (`1 << bits` table slots). Off below level 8; the table costs ~`16 << bits`
/// bytes per active block, so it is gated to the top levels where the user has
/// opted into maximum compression. `0` order disables it. Applies to the
/// non-reorder (`--order preserve`) path only.
fn level_to_hash(level: u8) -> (u8, u8) {
    match level {
        0..=7 => (0, 0),
        _ => (13, 25), // level 8 and up
    }
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
    // Color the user-facing summaries only on an interactive stderr with color
    // permitted; honor the `NO_COLOR` convention. Quiet mode prints no summary,
    // so color is moot there.
    report::set_color(
        !cli.quiet && io::stderr().is_terminal() && std::env::var_os("NO_COLOR").is_none(),
    );
    configure_global_pool(cli.threads).context("configuring the global thread pool")?;
    match cli.command {
        Command::Compress {
            inputs,
            output,
            force,
            level,
            block_reads,
            max,
            interleaved,
            order,
            no_rescue,
            keep_order,
            quality_bin,
            platform,
            estimate,
            verify,
        } => {
            if inputs.is_empty() {
                anyhow::bail!("at least one input FASTQ is required");
            }
            // `--max` is the best-ratio preset so users don't have to know the
            // knobs: top effort level plus a request to reorder. The library then
            // decides *dynamically* whether reordering actually pays off for the
            // input — it is applied to short reads and skipped for long reads —
            // so `--max` tracks "the best available for this data" rather than a
            // fixed pair of flags. Overrides `--level`/`--order`.
            let level = if max { 9 } else { level };
            let order = if max { ReadOrder::Any } else { order };
            let interleaved = interleaved.filter(|_| inputs.len() == 1);
            // `--order any`/`shuffle` turn on reordering; the library forces the
            // permutation (keep_order) back on for grouped input, so paired
            // archives still round-trip in order regardless. `shuffle` additionally
            // opts into name regeneration (reorder-lossy) for single-end input.
            let reorders = order != ReadOrder::Preserve;
            let params = fqxv::Params {
                seq_order: level_to_order(level),
                seq_hash_order: level_to_hash(level).0,
                seq_hash_bits: level_to_hash(level).1,
                // An explicit --block-reads overrides the level's block size,
                // decoupling random-access granularity from effort. A zero is
                // meaningless; treat it as "unset" and fall back to the level. Tiny
                // blocks no longer bomb: the sequence codec caps its per-block model
                // to the block's own base count (see `fqxv_seq::encode`).
                block_reads: block_reads
                    .filter(|&n| n > 0)
                    .unwrap_or(level_to_block(level)),
                quality_binning: quality_bin.into(),
                reorder: reorders,
                keep_order: keep_order && reorders,
                rescue: !no_rescue && reorders,
                regenerate_names: order == ReadOrder::Shuffle,
                threads: cli.threads,
                platform: platform.map(Into::into),
            };
            warn_redundant_binning(&inputs, params.quality_binning);
            // `--estimate` samples the input and reports a projected ratio/size
            // without writing an archive, so it short-circuits the whole compress
            // path (no output file is resolved or created).
            if let Some(fmt) = estimate {
                return run_estimate(&inputs, params, reorders, fmt);
            }
            let output = match output {
                Some(p) => p,
                None => default_archive_name(&inputs[0])?,
            };
            let in_size: u64 = inputs
                .iter()
                .filter_map(|p| std::fs::metadata(p).ok())
                .map(|m| m.len())
                .sum();
            let out = create_output(&output, force)?;

            // The reorder path opens with a long single-threaded prelude, so name
            // the spinner after the slow work to set expectations. On a `?` bail
            // the spinner's Drop clears the line before the error prints.
            let spin_msg = if params.reorder {
                "compressing (reorder — this can take a while)…"
            } else {
                "compressing…"
            };
            let sp = progress::Spinner::start(spin_msg, cli.quiet);

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
            sp.abandon();

            // Optional read-after-write verification: fully decode the archive we
            // just wrote and confirm it round-trips before reporting success. The
            // output handle was flushed and dropped by the compress call above, so
            // reopen the path. On any failure the archive is left in place (for
            // inspection) and we bail non-zero via `?` — a bad archive is never
            // reported as done. The spinner's Drop clears its line on that bail.
            if verify {
                let sp = progress::Spinner::start("verifying…", cli.quiet);
                let decoded = File::open(&output)
                    .with_context(|| format!("reopening {} to verify", output.display()))
                    .and_then(|f| {
                        fqxv::verify_roundtrip(f)
                            .with_context(|| format!("verifying {}", output.display()))
                    })?;
                if decoded != stats.reads {
                    anyhow::bail!(
                        "verification of {} decoded {decoded} reads but {} were written",
                        output.display(),
                        stats.reads
                    );
                }
                sp.abandon();
            }

            // User-facing summary: printed to stderr, gated only on `--quiet` and
            // deliberately independent of the tracing log level. The structured
            // line below is diagnostics only (visible at `-v`).
            if !cli.quiet {
                let names: Vec<String> = inputs.iter().map(|p| p.display().to_string()).collect();
                eprintln!(
                    "{}",
                    report::compress_summary(&names, &output, in_size, &stats, secs)
                );
                if verify {
                    eprintln!("verified: archive fully decodes ({} reads)", stats.reads);
                }
            }
            tracing::debug!(
                reads = stats.reads,
                inputs = inputs.len(),
                blocks = stats.blocks,
                out_bytes = stats.out_bytes,
                secs = format_args!("{secs:.1}"),
                "compressed"
            );
        }
        Command::Decompress {
            input,
            output,
            stdout,
            split,
            mate_style,
            no_gzip,
            recover,
            force,
        } => {
            let open_in =
                || File::open(&input).with_context(|| format!("opening input {}", input.display()));
            // Where the decoded FASTQ is headed, for the summary line. Computed
            // from borrows before `split`/`output` are consumed below.
            let dest = if let Some(prefix) = split.as_ref() {
                format!("{}_*", prefix.display())
            } else if stdout || output.as_deref().is_some_and(|p| p.as_os_str() == "-") {
                "stdout".to_string()
            } else if let Some(o) = output.as_ref() {
                o.display().to_string()
            } else {
                "stdout".to_string()
            };
            let t0 = Instant::now();
            if recover {
                // Recovery deliberately tolerates missing/corrupt blocks, so the
                // completeness check below does not apply here.
                let sp = progress::Spinner::start("recovering…", cli.quiet);
                let mut sink = open_sink(output.as_deref(), stdout, force)?;
                let rec = fqxv::decompress_recover(open_in()?, &mut sink, cli.threads)?;
                sink.finish()?;
                let secs = t0.elapsed().as_secs_f64();
                sp.abandon();
                if rec.blocks_skipped > 0 {
                    // On stderr so it shows even when decoded FASTQ is piped out.
                    eprintln!(
                        "warning: recovered {} block(s), skipped {} corrupt block(s) — {} read(s) lost",
                        rec.blocks_recovered, rec.blocks_skipped, rec.reads_lost
                    );
                }
                if !cli.quiet {
                    eprintln!(
                        "{}",
                        report::decompress_summary(&input, &dest, &rec.stats, secs)
                    );
                }
                tracing::debug!(
                    reads = rec.stats.reads,
                    blocks = rec.blocks_recovered,
                    skipped = rec.blocks_skipped,
                    reads_lost = rec.reads_lost,
                    secs = format_args!("{secs:.1}"),
                    "recovered"
                );
            } else {
                let sp = progress::Spinner::start("decompressing…", cli.quiet);
                // Read the footer's authoritative read count first. For a truncated
                // archive (which also lost its footer) this errors before any output
                // is written; otherwise we compare it against what we actually
                // decode so a silently short file is caught rather than trusted.
                let expected = fqxv::expected_reads(open_in()?)?;
                let stats = if let Some(prefix) = split {
                    let g = fqxv::peek(open_in()?)?.group_size as usize;
                    let mut sinks: Vec<FastqSink> = (1..=g)
                        .map(|i| {
                            open_fastq_file(&split_path(&prefix, i, mate_style, !no_gzip), force)
                        })
                        .collect::<anyhow::Result<_>>()?;
                    let stats = fqxv::decompress_split(open_in()?, &mut sinks, cli.threads)?;
                    for sink in sinks {
                        sink.finish()?;
                    }
                    stats
                } else {
                    let mut sink = open_sink(output.as_deref(), stdout, force)?;
                    let stats = fqxv::decompress(open_in()?, &mut sink, cli.threads)?;
                    sink.finish()?;
                    stats
                };
                if let Some(expected) = expected {
                    if stats.reads != expected {
                        anyhow::bail!(
                            "archive is truncated: footer declares {expected} reads but only \
                             {} were decoded — the output is incomplete",
                            stats.reads
                        );
                    }
                }
                let secs = t0.elapsed().as_secs_f64();
                sp.abandon();
                if !cli.quiet {
                    eprintln!(
                        "{}",
                        report::decompress_summary(&input, &dest, &stats, secs)
                    );
                }
                tracing::debug!(
                    reads = stats.reads,
                    blocks = stats.blocks,
                    secs = format_args!("{secs:.1}"),
                    "decompressed"
                );
            }
        }
        Command::Info {
            inputs,
            tsv,
            json,
            stats,
        } => print_info(&inputs, tsv, json, stats, cli.threads)?,
        Command::Verify {
            inputs,
            quick,
            tsv,
            json,
        } => print_verify(&inputs, quick, tsv, json)?,
    }
    Ok(())
}

/// A per-stream size row in the machine-readable [`InfoReport`].
#[derive(Debug, Serialize)]
struct StreamJson {
    bytes: u64,
    /// Share of the three compressed streams, in percent (0 when empty).
    pct: f64,
    /// This stream's compressed bytes divided by read count (null when there are
    /// no reads).
    #[serde(skip_serializing_if = "Option::is_none")]
    per_read: Option<f64>,
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
    /// On-disk container format version.
    format_version: u16,
    /// Stored whole-file CRC-32C as lowercase hex (absent for the footer-less
    /// whole-file-reorder layout and truncated archives). The value `verify`
    /// recomputes.
    #[serde(skip_serializing_if = "Option::is_none")]
    whole_file_crc: Option<String>,
    /// Content statistics from a full decode; present only with `--stats`.
    #[serde(skip_serializing_if = "Option::is_none")]
    stats: Option<StatsJson>,
}

/// Content statistics nested in [`InfoReport`] under `--stats` — mirrors
/// [`fqxv::ContentStats`] with derived percentages precomputed.
#[derive(Debug, Serialize)]
struct StatsJson {
    reads: u64,
    bases: u64,
    min_length: u32,
    max_length: u32,
    /// Null when there are no reads.
    #[serde(skip_serializing_if = "Option::is_none")]
    mean_length: Option<f64>,
    /// True when every read shares one length.
    fixed_length: bool,
    /// `(G + C) / (A + C + G + T)`; null for an all-N or empty archive.
    #[serde(skip_serializing_if = "Option::is_none")]
    gc_fraction: Option<f64>,
    /// Mean Phred quality over every base; null when there are no bases.
    #[serde(skip_serializing_if = "Option::is_none")]
    mean_quality: Option<f64>,
    base_composition: BaseCompositionJson,
    /// Per-Phred quality counts, only for values that occur (ascending).
    quality_histogram: Vec<QualBucketJson>,
}

/// Absolute base counts nested in [`StatsJson`].
#[derive(Debug, Serialize)]
struct BaseCompositionJson {
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
    other: u64,
}

/// One occupied Phred bucket in [`StatsJson::quality_histogram`].
#[derive(Debug, Serialize)]
struct QualBucketJson {
    phred: u8,
    count: u64,
}

/// The `streams` object nested in [`InfoReport`].
#[derive(Debug, Serialize)]
struct InfoStreams {
    names: StreamJson,
    sequence: StreamJson,
    quality: StreamJson,
    total: StreamJson,
}

/// Resolve `info`/`verify` inputs to a concrete `.fqxv` file list. A file path is
/// taken as-is; a directory is walked recursively for `*.fqxv` (sorted for a
/// deterministic order). Returns the file list and whether this is a *batch* — more
/// than one input path, or any directory — which selects the per-file
/// (filename-keyed) output. A single explicit file keeps the original
/// single-archive output, so scripts and the benchmark harness are unaffected.
fn resolve_fqxv_inputs(inputs: &[PathBuf]) -> anyhow::Result<(Vec<PathBuf>, bool)> {
    let mut files = Vec::new();
    let mut batch = inputs.len() > 1;
    for input in inputs {
        let meta =
            std::fs::metadata(input).with_context(|| format!("reading {}", input.display()))?;
        if meta.is_dir() {
            batch = true;
            let mut dir_files = Vec::new();
            collect_fqxv(input, &mut dir_files)?;
            // Sort within each directory subtree for a stable, reproducible order;
            // explicit file arguments keep the order the user gave them.
            dir_files.sort();
            files.extend(dir_files);
        } else {
            files.push(input.clone());
        }
    }
    if files.is_empty() {
        anyhow::bail!("no .fqxv files found in the given input(s)");
    }
    Ok((files, batch))
}

/// Recursively collect files with a `.fqxv` extension (case-insensitive) under
/// `dir` into `out`.
fn collect_fqxv(dir: &Path, out: &mut Vec<PathBuf>) -> anyhow::Result<()> {
    let entries =
        std::fs::read_dir(dir).with_context(|| format!("reading directory {}", dir.display()))?;
    for entry in entries {
        let path = entry
            .with_context(|| format!("reading directory {}", dir.display()))?
            .path();
        if path.is_dir() {
            collect_fqxv(&path, out)?;
        } else if path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("fqxv"))
        {
            out.push(path);
        }
    }
    Ok(())
}

/// One archive's inspected metadata plus the optional `--stats` content decode,
/// gathered once by [`gather_info`] and shared across the TSV/JSON/human renderers.
struct FileInfo {
    path: PathBuf,
    file_size: u64,
    info: fqxv::Info,
    content: Option<fqxv::ContentStats>,
    /// Whole-file CRC-32C as lowercase hex, or `None` for the footer-less layouts.
    crc_hex: Option<String>,
}

/// Inspect one archive's header/footer, and — with `stats` — decode it fully for
/// content statistics. Shared by every `info` output format.
fn gather_info(path: &Path, stats: bool, threads: usize) -> anyhow::Result<FileInfo> {
    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let info = fqxv::inspect(
        File::open(path).with_context(|| format!("opening input {}", path.display()))?,
    )?;
    // `--stats` requires a full decode (unlike the metadata above, which comes
    // from the header + footer index); do it once and share across formats.
    let content = if stats {
        Some(
            fqxv::content_stats(
                File::open(path).with_context(|| format!("opening input {}", path.display()))?,
                threads,
            )
            .with_context(|| format!("decoding {} for --stats", path.display()))?,
        )
    } else {
        None
    };
    let crc_hex = info.whole_file_crc.map(|c| format!("{c:08x}"));
    Ok(FileInfo {
        path: path.to_path_buf(),
        file_size,
        info,
        content,
        crc_hex,
    })
}

fn print_info(
    inputs: &[PathBuf],
    tsv: bool,
    json: bool,
    stats: bool,
    threads: usize,
) -> anyhow::Result<()> {
    let (files, batch) = resolve_fqxv_inputs(inputs)?;

    if tsv {
        // Header printed once. A batch prepends a `file` column and one row per
        // archive; a single file keeps the stable columns the benchmark harness
        // parses by position (see the `Info` doc comment).
        let mut header = String::new();
        if batch {
            header.push_str("file\t");
        }
        header.push_str(&info_tsv_header(stats));
        println!("{header}");
        for path in &files {
            let fi = gather_info(path, stats, threads)?;
            if batch {
                println!("{}\t{}", path.display(), info_tsv_row(&fi));
            } else {
                println!("{}", info_tsv_row(&fi));
            }
        }
        return Ok(());
    }

    if json {
        // A batch emits a JSON array; a single file stays a bare object so
        // existing single-archive consumers keep parsing one object.
        if batch {
            let reports = files
                .iter()
                .map(|p| gather_info(p, stats, threads).map(|fi| info_json_report(&fi)))
                .collect::<anyhow::Result<Vec<_>>>()?;
            println!("{}", serde_json::to_string_pretty(&reports)?);
        } else {
            let fi = gather_info(&files[0], stats, threads)?;
            println!("{}", serde_json::to_string_pretty(&info_json_report(&fi))?);
        }
        return Ok(());
    }

    for (i, path) in files.iter().enumerate() {
        if i > 0 {
            println!();
        }
        let fi = gather_info(path, stats, threads)?;
        print_info_human(&fi);
    }
    Ok(())
}

/// The stable `info --tsv` header columns (no leading `file` column — a batch
/// prepends that). `--stats` appends its columns only when set.
fn info_tsv_header(stats: bool) -> String {
    let mut header = String::from(
        "file_size\treads\tblocks\tgroup_size\tseq_order\tquality_binning\treordered\tnames_bytes\tseq_bytes\tqual_bytes\tplatform\tformat_version\twhole_file_crc",
    );
    if stats {
        header.push_str("\tbases\tmin_len\tmax_len\tgc_fraction\tmean_quality");
    }
    header
}

/// One archive's `info --tsv` data row (no leading `file` column). Keep field
/// order fixed and append new fields at the end so existing parsers don't break.
fn info_tsv_row(fi: &FileInfo) -> String {
    let info = &fi.info;
    let mut row = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        fi.file_size,
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
        info.format_version,
        fi.crc_hex.as_deref().unwrap_or(""),
    );
    if let Some(cs) = &fi.content {
        row.push_str(&format!(
            "\t{}\t{}\t{}\t{}\t{}",
            cs.bases,
            cs.min_len,
            cs.max_len,
            cs.gc_fraction()
                .map(|x| format!("{x:.6}"))
                .unwrap_or_default(),
            cs.mean_quality()
                .map(|x| format!("{x:.4}"))
                .unwrap_or_default(),
        ));
    }
    row
}

/// Build the [`InfoReport`] JSON shape for one archive.
fn info_json_report(fi: &FileInfo) -> InfoReport {
    let info = &fi.info;
    let content = &fi.content;
    let total = info.names_bytes + info.seq_bytes + info.qual_bytes;
    let pct = |x: u64| {
        if total > 0 {
            100.0 * x as f64 / total as f64
        } else {
            0.0
        }
    };
    let per_read = |x: u64| (info.reads > 0).then(|| x as f64 / info.reads as f64);
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
        4 => "lossy (Nanopore 4-bin)",
        5 => "lossy (PacBio HiFi 5-bin)",
        _ => "unknown",
    };
    let spots = (info.group_size > 1).then(|| info.reads / info.group_size.max(1) as u64);
    let bytes_per_read = (info.reads > 0).then(|| total as f64 / info.reads as f64);

    let stream = |bytes: u64| StreamJson {
        bytes,
        pct: pct(bytes),
        per_read: per_read(bytes),
    };
    InfoReport {
        file: fi.path.display().to_string(),
        file_size: fi.file_size,
        file_size_human: human_bytes(fi.file_size),
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
            names: stream(info.names_bytes),
            sequence: stream(info.seq_bytes),
            quality: stream(info.qual_bytes),
            total: stream(total),
        },
        bytes_per_read,
        format_version: info.format_version,
        whole_file_crc: fi.crc_hex.clone(),
        stats: content.as_ref().map(|cs| StatsJson {
            reads: cs.reads,
            bases: cs.bases,
            min_length: cs.min_len,
            max_length: cs.max_len,
            mean_length: cs.mean_len(),
            fixed_length: cs.fixed_length(),
            gc_fraction: cs.gc_fraction(),
            mean_quality: cs.mean_quality(),
            base_composition: BaseCompositionJson {
                a: cs.a,
                c: cs.c,
                g: cs.g,
                t: cs.t,
                n: cs.n,
                other: cs.other,
            },
            quality_histogram: cs
                .qual_hist
                .iter()
                .enumerate()
                .filter(|(_, &count)| count > 0)
                .map(|(phred, &count)| QualBucketJson {
                    phred: phred as u8,
                    count,
                })
                .collect(),
        }),
    }
}

/// Print the human `info` report for one archive: a metadata table, a per-stream
/// size table, and — with `--stats` — content statistics and a quality histogram.
fn print_info_human(fi: &FileInfo) {
    let info = &fi.info;
    let content = &fi.content;
    let total = info.names_bytes + info.seq_bytes + info.qual_bytes;
    let pct = |x: u64| {
        if total > 0 {
            100.0 * x as f64 / total as f64
        } else {
            0.0
        }
    };
    let per_read = |x: u64| (info.reads > 0).then(|| x as f64 / info.reads as f64);
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
        4 => "lossy (Nanopore 4-bin)",
        5 => "lossy (PacBio HiFi 5-bin)",
        _ => "unknown",
    };
    let spots = (info.group_size > 1).then(|| info.reads / info.group_size.max(1) as u64);
    let bytes_per_read = (info.reads > 0).then(|| total as f64 / info.reads as f64);

    // Human report: a metadata table, a per-stream size table, and — with
    // `--stats` — a content-statistics table plus a quality histogram.
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
    // Blocks, annotated with the average reads per block when both are known.
    let blocks = if info.blocks > 0 && info.reads > 0 {
        format!(
            "{} (avg {} reads)",
            group_digits(info.blocks),
            group_digits(info.reads / info.blocks)
        )
    } else {
        group_digits(info.blocks)
    };
    meta_row("blocks", blocks);
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
        "format",
        format!(
            "v{}.{}",
            info.format_version >> 8,
            info.format_version & 0xff
        ),
    );
    meta_row(
        "whole-file crc",
        fi.crc_hex.clone().unwrap_or_else(|| "—".to_string()),
    );
    meta_row(
        "file size",
        format!(
            "{} ({} bytes)",
            human_bytes(fi.file_size),
            group_digits(fi.file_size)
        ),
    );
    let mut meta = meta.build();
    meta.with(Style::rounded());

    let mut streams = TableBuilder::default();
    streams.push_record([
        "stream".to_string(),
        "bytes".to_string(),
        "share".to_string(),
        "bytes/read".to_string(),
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
            per_read(bytes).map_or_else(|| "—".to_string(), |x| format!("{x:.3}")),
        ]);
    }
    let mut streams = streams.build();
    streams
        .with(Style::rounded())
        .with(Modify::new(Columns::new(1..)).with(Alignment::right()));

    println!("{}", fi.path.display());
    println!("{meta}");
    println!("{streams}");
    if let Some(bpr) = bytes_per_read {
        println!("{bpr:.2} bytes/read");
    }
    if let Some(cs) = content {
        print!("{}", render_content_stats(cs));
    }
}

/// Render the `--stats` content-statistics table and quality histogram as a
/// printable block (trailing newline included). Kept separate from
/// [`print_info`] so the metadata report stays readable.
fn render_content_stats(cs: &fqxv::ContentStats) -> String {
    use std::fmt::Write as _;

    let frac_of_bases = |x: u64| {
        if cs.bases > 0 {
            100.0 * x as f64 / cs.bases as f64
        } else {
            0.0
        }
    };
    let mut t = TableBuilder::default();
    t.push_record(["content".to_string(), "value".to_string()]);
    let mut row = |k: &str, v: String| t.push_record([k.to_string(), v]);
    row("reads", group_digits(cs.reads));
    row("bases", group_digits(cs.bases));
    let length = if cs.reads == 0 {
        "—".to_string()
    } else if cs.fixed_length() {
        format!("{} (fixed)", cs.min_len)
    } else {
        format!(
            "{}–{} (mean {:.1})",
            cs.min_len,
            cs.max_len,
            cs.mean_len().unwrap_or(0.0)
        )
    };
    row("read length", length);
    row(
        "GC content",
        cs.gc_fraction()
            .map_or_else(|| "—".to_string(), |g| format!("{:.1}%", 100.0 * g)),
    );
    row(
        "A / C / G / T",
        format!(
            "{:.1}% / {:.1}% / {:.1}% / {:.1}%",
            frac_of_bases(cs.a),
            frac_of_bases(cs.c),
            frac_of_bases(cs.g),
            frac_of_bases(cs.t),
        ),
    );
    row(
        "N bases",
        format!("{} ({:.2}%)", group_digits(cs.n), frac_of_bases(cs.n)),
    );
    if cs.other > 0 {
        row(
            "other bases",
            format!(
                "{} ({:.2}%)",
                group_digits(cs.other),
                frac_of_bases(cs.other)
            ),
        );
    }
    row(
        "mean quality",
        cs.mean_quality()
            .map_or_else(|| "—".to_string(), |q| format!("Q{q:.1}")),
    );
    let mut table = t.build();
    table.with(Style::rounded());

    let mut out = format!("{table}\n");
    // Quality distribution, bucketed into width-5 Phred ranges; only occupied
    // ranges are shown, each with a bar scaled to the busiest range.
    let mut buckets: Vec<(usize, u64)> = Vec::new();
    let mut lo = 0;
    while lo < fqxv::QUAL_MAX {
        let hi = (lo + 5).min(fqxv::QUAL_MAX);
        let count: u64 = cs.qual_hist[lo..hi].iter().sum();
        if count > 0 {
            buckets.push((lo, count));
        }
        lo = hi;
    }
    // Buckets are only pushed when non-empty, so `max` is > 0 whenever the loop
    // runs — no division guard needed.
    let max = buckets.iter().map(|&(_, c)| c).max().unwrap_or(0);
    if !buckets.is_empty() {
        out.push_str("quality distribution\n");
        for (lo, count) in buckets {
            let bar = (24 * count / max) as usize;
            let _ = writeln!(
                out,
                "  Q{:>2}–{:<2} {:<24} {:>5.1}%",
                lo,
                lo + 4,
                "█".repeat(bar),
                frac_of_bases(count),
            );
        }
    }
    out
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

/// Verify one or more archives and render each as a table (default), TSV, or
/// JSON. A single file keeps the original single-archive output; a batch (several
/// files or a directory) prepends a `file` column (TSV) or emits a JSON array,
/// and stays resilient — an unreadable archive becomes a failing entry instead of
/// aborting the run. Exits the process non-zero when *any* archive is corrupt so
/// scripts can branch on `$?` regardless of the chosen format.
fn print_verify(inputs: &[PathBuf], quick: bool, tsv: bool, json: bool) -> anyhow::Result<()> {
    let (files, batch) = resolve_fqxv_inputs(inputs)?;

    // Verify each archive. In single-file mode a bad archive is a hard error (as
    // before); in a batch it becomes a failing entry so one corrupt file doesn't
    // abort the rest.
    let mut results: Vec<(PathBuf, Result<fqxv::VerifyReport, String>)> =
        Vec::with_capacity(files.len());
    for path in &files {
        let r = File::open(path)
            .map_err(|e| format!("opening input {}: {e}", path.display()))
            .and_then(|f| {
                fqxv::verify_report(&f, quick)
                    .map_err(|e| format!("{} is not a readable fqxv archive: {e}", path.display()))
            });
        if !batch {
            // Single file: surface the error with context, exactly as before.
            results.push((path.clone(), Ok(r.map_err(anyhow::Error::msg)?)));
        } else {
            results.push((path.clone(), r));
        }
    }

    if json {
        let entries: Vec<VerifyJson> = results.iter().map(|(p, r)| verify_json_of(p, r)).collect();
        if batch {
            println!("{}", serde_json::to_string_pretty(&entries)?);
        } else {
            println!("{}", serde_json::to_string_pretty(&entries[0])?);
        }
    } else if tsv {
        let mut header = String::new();
        if batch {
            header.push_str("file\t");
        }
        header.push_str("check\tresult\tdetail");
        println!("{header}");
        for (p, r) in &results {
            let prefix = if batch {
                format!("{}\t", p.display())
            } else {
                String::new()
            };
            match r {
                Ok(rep) => {
                    for c in &rep.checks {
                        println!(
                            "{prefix}{}\t{}\t{}",
                            c.name,
                            if c.ok { "ok" } else { "fail" },
                            c.detail
                        );
                    }
                }
                Err(e) => println!("{prefix}readable\tfail\t{e}"),
            }
        }
    } else {
        for (i, (p, r)) in results.iter().enumerate() {
            if i > 0 {
                println!();
            }
            match r {
                Ok(rep) => print_verify_human_one(p, rep),
                Err(e) => println!("{}: CORRUPT ({e})", p.display()),
            }
        }
    }

    let all_passed = results
        .iter()
        .all(|(_, r)| r.as_ref().map(|rep| rep.passed()).unwrap_or(false));
    if !all_passed {
        std::process::exit(1);
    }
    Ok(())
}

/// Build the [`VerifyJson`] for one archive, mapping an unreadable archive to a
/// single failing `readable` check.
fn verify_json_of(path: &Path, r: &Result<fqxv::VerifyReport, String>) -> VerifyJson {
    match r {
        Ok(rep) => VerifyJson {
            file: path.display().to_string(),
            passed: rep.passed(),
            checks: rep
                .checks
                .iter()
                .map(|c| VerifyCheckJson {
                    name: c.name.clone(),
                    ok: c.ok,
                    detail: c.detail.clone(),
                })
                .collect(),
            failed_blocks: rep.failed_blocks.clone(),
        },
        Err(e) => VerifyJson {
            file: path.display().to_string(),
            passed: false,
            checks: vec![VerifyCheckJson {
                name: "readable".to_string(),
                ok: false,
                detail: e.clone(),
            }],
            failed_blocks: Vec::new(),
        },
    }
}

/// Print the human `verify` report for one archive: the per-check table and a
/// final `path: OK/CORRUPT` line.
fn print_verify_human_one(path: &Path, report: &fqxv::VerifyReport) {
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
        if report.passed() { "OK" } else { "CORRUPT" }
    );
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
/// Warn when the requested lossy quality binning can't actually shrink the
/// input, so it would only add distortion. Binning to N levels is pointless when
/// the quality stream already has no more than N distinct levels — it can't
/// reduce cardinality, and if its representatives don't line up with the native
/// ones it re-quantizes nearly every base off-grid for no size gain (the classic
/// case: `--quality-bin bin4` on already-4-level NovaSeq quality). Peeks the
/// first input's quality; silent on stdin (can't rewind) or when inconclusive.
fn warn_redundant_binning(inputs: &[PathBuf], binning: fqxv::QualityBinning) {
    let (name, target_levels) = match binning {
        fqxv::QualityBinning::Lossless => return,
        fqxv::QualityBinning::Bin8 => ("bin8", 8usize),
        fqxv::QualityBinning::Bin4 => ("bin4", 4),
        fqxv::QualityBinning::Bin2 => ("bin2", 2),
        fqxv::QualityBinning::BinOnt => ("ont", 4),
        fqxv::QualityBinning::BinHifi => ("hifi", 5),
    };
    let Some(path) = inputs.first() else { return };
    if path.as_os_str() == "-" {
        return; // stdin: peeking would consume the stream the compressor needs.
    }
    let Ok(mut rdr) = open_input(path) else {
        return;
    };
    let mut buf = vec![0u8; 1 << 20];
    let mut filled = 0;
    while filled < buf.len() {
        match rdr.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(_) => return,
        }
    }
    buf.truncate(filled);
    // Quality is line 4 of each 4-line FASTQ record; a trailing partial line is
    // ignored. (fqxv targets single-line-per-field FASTQ.)
    let mut seen = [false; 256];
    let (mut distinct, mut sampled, mut changed) = (0usize, 0u64, 0u64);
    for (i, line) in buf.split(|&b| b == b'\n').enumerate() {
        if i % 4 != 3 || line.is_empty() {
            continue;
        }
        for &q in line {
            if !std::mem::replace(&mut seen[q as usize], true) {
                distinct += 1;
            }
            sampled += 1;
            changed += u64::from(binning.apply(q) != q);
        }
    }
    if sampled == 0 || distinct > target_levels {
        return;
    }
    let pct = 100.0 * changed as f64 / sampled as f64;
    let hint = if target_levels > 2 {
        " Use lossless, or a coarser bin (e.g. bin2) if you want a real size reduction."
    } else {
        " Use lossless instead."
    };
    warn!(
        "input quality already has {distinct} distinct level(s); --quality-bin {name} \
         targets {target_levels} and re-quantizes {pct:.0}% of sampled bases — this adds \
         distortion for little or no size gain.{hint}"
    );
}

fn open_input(path: &Path) -> anyhow::Result<Box<dyn Read + Send>> {
    decode_input(raw_input(path)?)
}

/// The raw (still-compressed) input source: stdin for `-`, else the file.
fn raw_input(path: &Path) -> anyhow::Result<Box<dyn Read + Send>> {
    if path.as_os_str() == "-" {
        Ok(Box::new(io::stdin()))
    } else {
        Ok(Box::new(File::open(path).with_context(|| {
            format!("opening input {}", path.display())
        })?))
    }
}

/// Wrap a raw source in the right decoder, sniffing gzip/BGZF from the leading
/// magic bytes (identical decoded stream either way — see the module note above).
fn decode_input(mut src: Box<dyn Read + Send>) -> anyhow::Result<Box<dyn Read + Send>> {
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

/// Like [`open_input`], but tallies bytes pulled from the *raw* (pre-decode)
/// source into the returned counter. `--estimate` uses this to learn how much
/// on-disk input a sample consumed, so it can project the whole-file archive
/// size — for a gzip input the count is on-disk (compressed) bytes.
fn open_input_counted(path: &Path) -> anyhow::Result<(Box<dyn Read + Send>, Arc<AtomicU64>)> {
    let count = Arc::new(AtomicU64::new(0));
    let counted: Box<dyn Read + Send> = Box::new(CountingReader {
        inner: raw_input(path)?,
        count: Arc::clone(&count),
    });
    Ok((decode_input(counted)?, count))
}

/// A pass-through reader that adds every byte it yields to a shared counter.
struct CountingReader<R> {
    inner: R,
    count: Arc<AtomicU64>,
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.count.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

/// Cap on the number of reads coded for `--estimate`. One default-level block is
/// `1 << 20` reads; a sample this size codes in well under a second and warms the
/// context models the same way a real block does, so the sample's ratio tracks
/// the full run. Split across inputs so grouped input still samples ~one block.
const ESTIMATE_SAMPLE_READS: usize = 1 << 20;

/// One input's estimate: the coded sample sizes, the on-disk bytes the sample
/// consumed (compressed, for a gzip input), and the input's total on-disk size
/// (`None` for stdin, which has no length).
struct InputEstimate {
    est: fqxv::Estimate,
    consumed: u64,
    disk: Option<u64>,
}

/// Run `--estimate`: sample each input, code it with the real codecs, and print a
/// projected ratio and archive size without writing an archive. `reorders` is the
/// resolved `--order any`/`--max` flag; when set the estimate is a lower bound
/// (reordering is not modeled — see the note in the report and `fqxv::estimate`).
fn run_estimate(
    inputs: &[PathBuf],
    params: fqxv::Params,
    reorders: bool,
    fmt: EstimateFormat,
) -> anyhow::Result<()> {
    let g = inputs.len().max(1);
    // Keep the aggregate sample ~one block: split the cap across grouped inputs.
    let per_input = (ESTIMATE_SAMPLE_READS / g).max(1);

    // Sample every input, tracking whether a whole-file size projection is
    // possible (it is unless a non-exhausted stdin stream hides the total).
    let mut parts = Vec::with_capacity(inputs.len());
    for path in inputs {
        let (reader, counter) = open_input_counted(path)?;
        let est = fqxv::estimate(reader, params, per_input)
            .with_context(|| format!("estimating {}", path.display()))?;
        let consumed = counter.load(Ordering::Relaxed);
        let disk = if path.as_os_str() == "-" {
            None
        } else {
            std::fs::metadata(path).ok().map(|m| m.len())
        };
        parts.push(InputEstimate {
            est,
            consumed,
            disk,
        });
    }

    // Aggregate the sample, and project each input's whole-file contribution by
    // the fraction of its bytes the sample consumed (an exhausted input needs no
    // projection — its sample *is* the whole file).
    let mut reads = 0u64;
    let mut bases = 0u64;
    let (mut names, mut seq, mut qual) = (0u64, 0u64, 0u64);
    let mut sample_raw = 0u64;
    let mut sample_archive = 0u64;
    let (mut proj_raw, mut proj_archive, mut disk_in) = (0f64, 0f64, 0u64);
    let mut projectable = true;
    for p in &parts {
        let e = &p.est;
        reads += e.sample_reads;
        bases += e.sample_bases;
        names += e.names_bytes;
        seq += e.seq_bytes;
        qual += e.qual_bytes;
        sample_raw += e.raw_bytes;
        sample_archive += e.archive_bytes;
        if e.exhausted {
            // Whole input coded: exact, and its on-disk size is what we consumed.
            proj_raw += e.raw_bytes as f64;
            proj_archive += e.archive_bytes as f64;
            disk_in += p.disk.unwrap_or(p.consumed);
        } else if let Some(d) = p.disk.filter(|&d| d > 0 && p.consumed > 0) {
            let frac = (p.consumed as f64 / d as f64).min(1.0);
            proj_raw += e.raw_bytes as f64 / frac;
            proj_archive += e.archive_bytes as f64 / frac;
            disk_in += d;
        } else {
            // A non-exhausted stream with no known length: report the sample only.
            projectable = false;
        }
    }
    // Machine-readable table: the input file(s), current on-disk size, projected
    // archive size, and their ratio — a header line then one data row (a compress
    // estimate is one archive, even when several inputs interleave into it). When
    // the whole-file total is known (files, or a fully-consumed stream) it reports
    // the projection; a streaming input of unknown length falls back to the
    // sample's own figures.
    if fmt == EstimateFormat::Tsv {
        let (cur, est, ratio) = if projectable && proj_archive > 0.0 && disk_in > 0 {
            (disk_in, proj_archive as u64, disk_in as f64 / proj_archive)
        } else if sample_archive > 0 {
            (
                sample_raw,
                sample_archive,
                sample_raw as f64 / sample_archive as f64,
            )
        } else {
            (sample_raw, sample_archive, 0.0)
        };
        // Interleaved inputs all feed one archive, so join their names into the
        // single `file` field (`stdin` for a `-` stream).
        let file = inputs
            .iter()
            .map(|p| {
                if p.as_os_str() == "-" {
                    "stdin".to_string()
                } else {
                    p.display().to_string()
                }
            })
            .collect::<Vec<_>>()
            .join(",");
        println!("file\tinput_bytes\test_fqxv_bytes\tratio");
        println!("{file}\t{cur}\t{est}\t{ratio:.4}");
        return Ok(());
    }

    let streams = names + seq + qual;
    let name_of = |p: &Path| {
        if p.as_os_str() == "-" {
            "stdin".to_string()
        } else {
            p.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| p.display().to_string())
        }
    };

    let mut out = String::new();
    use std::fmt::Write as _;

    // Lead with the punchline: input size(s) -> estimated fqxv size, % smaller.
    if projectable && proj_archive > 0.0 && disk_in > 0 {
        let ratio_disk = disk_in as f64 / proj_archive;
        let pct = (1.0 - proj_archive / disk_in as f64) * 100.0;
        let archive = human_bytes(proj_archive as u64);
        if inputs.len() == 1 {
            let _ = writeln!(
                out,
                "{} ({})  \u{2192}  estimated fqxv ~{}  ({:.0}% smaller, ~{:.2}x)",
                name_of(&inputs[0]),
                human_bytes(disk_in),
                archive,
                pct,
                ratio_disk,
            );
        } else {
            let _ = writeln!(
                out,
                "{} FASTQ files ({} total)  \u{2192}  estimated fqxv ~{}  ({:.0}% smaller, ~{:.2}x)",
                inputs.len(),
                human_bytes(disk_in),
                archive,
                pct,
                ratio_disk,
            );
            for (path, part) in inputs.iter().zip(&parts) {
                let d = part.disk.unwrap_or(part.consumed);
                let _ = writeln!(out, "    {:<28} {:>11}", name_of(path), human_bytes(d));
            }
        }
    } else {
        // Streaming input of unknown length: no on-disk total to reduce against,
        // so lead with the sample's own reduction (vs uncompressed FASTQ).
        let sample_ratio = if sample_archive == 0 {
            0.0
        } else {
            sample_raw as f64 / sample_archive as f64
        };
        let pct = if sample_raw > 0 {
            (1.0 - sample_archive as f64 / sample_raw as f64) * 100.0
        } else {
            0.0
        };
        let _ = writeln!(
            out,
            "sample: {} uncompressed FASTQ  \u{2192}  ~{} fqxv  ({:.0}% smaller, ~{:.2}x)",
            human_bytes(sample_raw),
            human_bytes(sample_archive),
            pct,
            sample_ratio,
        );
        let _ = writeln!(
            out,
            "(streaming input has no on-disk size; give a file to project the whole-file total)"
        );
    }

    // Per-stream breakdown from the sample: compressed size, share of the coded
    // streams, and the natural rate for each (names amortize over reads; seq/qual
    // over bases).
    let share = |b: u64| {
        if streams == 0 {
            0.0
        } else {
            b as f64 / streams as f64 * 100.0
        }
    };
    let per_read = |b: u64| {
        if reads == 0 {
            0.0
        } else {
            b as f64 * 8.0 / reads as f64
        }
    };
    let per_base = |b: u64| {
        if bases == 0 {
            0.0
        } else {
            b as f64 * 8.0 / bases as f64
        }
    };
    let _ = writeln!(
        out,
        "\nEstimated from a {}-read sample ({} uncompressed FASTQ):",
        group_digits(reads),
        human_bytes(sample_raw),
    );
    let _ = writeln!(out, "  stream      compressed     share   rate");
    let _ = writeln!(
        out,
        "  names     {:>11}   {:>5.1}%   {:.2} bits/read",
        human_bytes(names),
        share(names),
        per_read(names),
    );
    let _ = writeln!(
        out,
        "  sequence  {:>11}   {:>5.1}%   {:.3} bits/base",
        human_bytes(seq),
        share(seq),
        per_base(seq),
    );
    let _ = writeln!(
        out,
        "  quality   {:>11}   {:>5.1}%   {:.3} bits/base",
        human_bytes(qual),
        share(qual),
        per_base(qual),
    );

    // Secondary figure: the vs-uncompressed-FASTQ ratio — the data-level
    // compression, distinct from the vs-on-disk reduction in the headline.
    if projectable && proj_archive > 0.0 {
        let _ = writeln!(
            out,
            "  vs uncompressed FASTQ (~{}):  ~{:.2}x",
            human_bytes(proj_raw as u64),
            proj_raw / proj_archive,
        );
    }

    if reorders {
        let _ = writeln!(
            out,
            "\nnote: --order/--max reordering is not modeled; the real archive will\n\
             be this size or smaller."
        );
    }

    print!("{out}");
    Ok(())
}

/// A FASTQ output sink: raw bytes, or BGZF (block-gzip) compression.
///
/// BGZF is selected by a `.gz` filename ([`is_gzip_path`]); stdout is always raw so
/// a decode piped straight into an aligner isn't gzip-wrapped. The multithreaded
/// BGZF writer compresses blocks on rayon's global pool (sized by `--threads`), so
/// gzipping the output doesn't serialize the otherwise-parallel decode.
enum FastqSink {
    Raw(Box<dyn Write + Send>),
    Bgzf(noodles_bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>),
}

impl Write for FastqSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            FastqSink::Raw(w) => w.write(buf),
            FastqSink::Bgzf(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            FastqSink::Raw(w) => w.flush(),
            FastqSink::Bgzf(w) => w.flush(),
        }
    }
}

impl FastqSink {
    /// Wrap `inner` in a BGZF encoder when `gzip`, else pass it through raw.
    fn new(inner: Box<dyn Write + Send>, gzip: bool) -> Self {
        if gzip {
            FastqSink::Bgzf(noodles_bgzf::io::MultithreadedWriter::new(inner))
        } else {
            FastqSink::Raw(inner)
        }
    }

    /// Flush and, for BGZF, append the mandatory EOF marker block. The library
    /// decode paths only flush their writers, so finalization lives here; a dropped
    /// `MultithreadedWriter` also writes the EOF block, but calling `finish`
    /// surfaces a write error (a full disk, say) instead of swallowing it in `Drop`.
    fn finish(self) -> io::Result<()> {
        match self {
            FastqSink::Raw(mut w) => w.flush(),
            FastqSink::Bgzf(mut w) => w.finish().and_then(|mut inner| inner.flush()),
        }
    }
}

/// True when `path` ends in a `.gz` extension (case-insensitive), i.e. the caller
/// asked for gzip-compressed output.
fn is_gzip_path(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
}

/// Resolve the interleaved-output sink from `-o`/`-Z`. An explicit path (or `-o -`)
/// wins; `-Z/--stdout` streams raw to stdout; with neither we refuse rather than
/// dump a whole FASTQ to the terminal. Stdout is always raw so a piped decode isn't
/// gzip-wrapped; a file's `.gz` extension selects BGZF.
fn open_sink(output: Option<&Path>, to_stdout: bool, force: bool) -> anyhow::Result<FastqSink> {
    match output {
        Some(p) if p.as_os_str() == "-" => Ok(FastqSink::Raw(Box::new(io::stdout()))),
        Some(p) => open_fastq_file(p, force),
        None if to_stdout => Ok(FastqSink::Raw(Box::new(io::stdout()))),
        None => anyhow::bail!(
            "no output specified: pass -o FILE (\".gz\" for BGZF), --split PREFIX for separate \
             mate files, or -Z/--stdout to stream interleaved FASTQ to stdout"
        ),
    }
}

/// Create a FASTQ output file, choosing BGZF compression from its `.gz` extension.
fn open_fastq_file(path: &Path, force: bool) -> anyhow::Result<FastqSink> {
    let file = create_output(path, force)?;
    Ok(FastqSink::new(Box::new(file), is_gzip_path(path)))
}

/// Create `path` for writing, refusing to clobber an existing file unless `force`.
/// Without `force` the create is atomic (`create_new`), so it fails cleanly if the
/// path appears between a check and the open rather than truncating a file another
/// process just wrote. With `force` this truncates any existing file, as before.
fn create_output(path: &Path, force: bool) -> anyhow::Result<File> {
    if force {
        return File::create(path).with_context(|| format!("creating output {}", path.display()));
    }
    match File::options().write(true).create_new(true).open(path) {
        Ok(f) => Ok(f),
        Err(e) if e.kind() == io::ErrorKind::AlreadyExists => anyhow::bail!(
            "output {} already exists; pass -f/--force to overwrite",
            path.display()
        ),
        Err(e) => Err(e).with_context(|| format!("creating output {}", path.display())),
    }
}

/// Build the path for member `i` (1-based) of a `--split` decode: `<prefix>` plus a
/// mate label (`_R1` or `_1`, per `style`) and a `.fastq`/`.fastq.gz` extension.
fn split_path(prefix: &Path, i: usize, style: MateStyle, gzip: bool) -> PathBuf {
    let label = match style {
        MateStyle::R => format!("R{i}"),
        MateStyle::Num => i.to_string(),
    };
    let ext = if gzip { "fastq.gz" } else { "fastq" };
    PathBuf::from(format!("{}_{}.{}", prefix.display(), label, ext))
}

/// Derive the default archive name from the first input: strip a gzip suffix and a
/// FASTQ extension, append `.fqxv`, and keep it in the input's directory
/// (`path/to/reads.fastq.gz` -> `path/to/reads.fqxv`). Errors on stdin (`-`), which
/// has no name to derive from.
fn default_archive_name(first: &Path) -> anyhow::Result<PathBuf> {
    if first.as_os_str() == "-" {
        anyhow::bail!("reading FASTQ from stdin (-) has no name to derive from; pass -o/--output");
    }
    let name = first
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let base = name.strip_suffix(".gz").unwrap_or(&name);
    let base = base
        .strip_suffix(".fastq")
        .or_else(|| base.strip_suffix(".fq"))
        .unwrap_or(base);
    let base = if base.is_empty() { "out" } else { base };
    Ok(first.with_file_name(format!("{base}.fqxv")))
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
    fn create_output_refuses_existing_unless_forced() {
        let dir = std::env::temp_dir().join(format!("fqxv-create-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("archive.fqxv");

        // A fresh path is created regardless of `force`.
        create_output(&path, false).unwrap();
        assert!(path.exists());
        std::fs::write(&path, b"existing contents").unwrap();

        // Without force, an existing path errors and is left untouched.
        let err = create_output(&path, false).unwrap_err().to_string();
        assert!(err.contains("already exists"), "unexpected error: {err}");
        assert_eq!(std::fs::read(&path).unwrap(), b"existing contents");

        // With force, the file is truncated for rewriting.
        create_output(&path, true).unwrap();
        assert!(std::fs::read(&path).unwrap().is_empty());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn resolve_threads_zero_is_all_cores_and_explicit_is_clamped() {
        let available = std::thread::available_parallelism().map_or(1, |n| n.get());
        assert_eq!(resolve_threads(0), available); // 0 = all cores
        assert_eq!(resolve_threads(1), 1);
        assert_eq!(resolve_threads(usize::MAX), available); // never oversubscribe
    }
}
