//! `fqxv` command-line interface — a thin front-end over the [`fqxv`] library.

mod output;
mod progress;
mod report;

use std::fs::File;
use std::io::{self, IsTerminal, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
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
use tracing_subscriber::{EnvFilter, fmt};

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

/// Version reported by `--version`: the crate version alone for a clean
/// release-tag build, or with the git description appended for anything else
/// (`0.3.0 (v0.3.0-7-gab12cd34-dirty)`). Assembled in `build.rs`, which also
/// explains the release/development distinction.
const VERSION: &str = env!("FQXV_VERSION");

/// Reference-free FASTQ archiver for short-read data.
#[derive(Debug, Parser)]
#[command(name = "fqxv", version = VERSION, about, long_about = None, styles = STYLES, after_help = EXAMPLES)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Number of worker threads (0 = all cores).
    ///
    /// Capped at the number of cores that physically exist, so an explicit
    /// value never oversubscribes. Governs both compression and the parallel
    /// BGZF decode path.
    #[arg(long, global = true, default_value_t = 16)]
    threads: usize,

    /// Increase log verbosity (repeatable: -v, -vv, -vvv).
    ///
    /// -v adds info, -vv adds debug, -vvv adds trace with targets, thread ids,
    /// and span timing. All logs go to stderr, so piped FASTQ stays clean.
    /// Overridden by RUST_LOG if set.
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Silence all output except warnings and errors.
    ///
    /// Suppresses the progress indicator and the end-of-run summary; only
    /// warnings and errors reach stderr.
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,
}

/// Install the `tracing` subscriber. Verbosity comes from `-v/-vv/-vvv` and
/// `--quiet`; a set `RUST_LOG` overrides the computed level entirely. All output
/// goes to **stderr** so decompressed FASTQ on stdout stays uncontaminated.
///
/// The default (no `-v`) is `warn`: the user-facing surface is the progress
/// indicator plus the end-of-run summary, so routine `info` diagnostics — which
/// otherwise interleave with and smear the live spinner/bar — are held back to
/// `-v`. `-vv` adds `debug`, `-vvv` adds `trace` with targets and timing.
fn init_tracing(verbose: u8, quiet: bool) {
    let level = if quiet {
        "error"
    } else {
        match verbose {
            0 => "warn",
            1 => "info",
            2 => "debug",
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
    /// Compress FASTQ to `.fqxv`.
    ///
    /// Give multiple inputs to interleave per-spot files (paired mates, or
    /// single-cell R1/R2/I1/I2) into one archive.
    Compress {
        /// Input FASTQ file(s), or `-` for stdin.
        ///
        /// Plain or gzipped; `-` reads one stream from stdin. One file =
        /// single-end, 2 = paired, 3-4 = single-cell; multiple files are
        /// interleaved per spot. Order is preserved for `--split`.
        #[arg(num_args = 1..)]
        inputs: Vec<PathBuf>,
        /// Output `.fqxv` path (default: alongside the input).
        ///
        /// Defaults to the first input's name with the FASTQ/gzip extension
        /// replaced by `.fqxv` (`reads.fastq.gz` -> `reads.fqxv`), written
        /// alongside the input. Required when the input is stdin (`-`).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Overwrite the output archive if it already exists.
        ///
        /// Without this, an existing output path is left untouched and the
        /// command errors, so a stray `compress` can't silently clobber an
        /// archive.
        #[arg(short = 'f', long)]
        force: bool,
        /// Maximum-compression preset (overrides `--level`/`--order`).
        ///
        /// The deepest sequence context plus read reordering *where it helps*.
        /// Reordering is applied to short-read data and automatically skipped
        /// for long reads (nanopore/PacBio), where it costs ~10x the time and
        /// memory for no ratio gain — so `--max` adapts to the input rather than
        /// forcing one fixed setting. Single-end short reads may come back
        /// reordered; names, sequence, and quality are preserved exactly (still
        /// lossless).
        #[arg(long)]
        max: bool,
        /// Re-decode the archive after writing to confirm it round-trips.
        ///
        /// Recommended before deleting the source FASTQ. Catches a codec or
        /// memory error that produced a CRC-valid-but-wrong archive. Roughly
        /// doubles wall time (adds a decode pass; the reorder/`--max` path is
        /// the heaviest). On failure the archive is left in place and the
        /// command exits non-zero.
        #[arg(long, conflicts_with = "estimate")]
        verify: bool,
        /// Estimate the ratio and archive size from a sample, then exit.
        ///
        /// Writes nothing. Codes the leading reads with the real codecs at the
        /// chosen `--level`/`--quality-bin` and projects the whole-file archive.
        /// Reordering (`--order any`/`--max`) is not modeled, so with it the
        /// estimate is a conservative lower bound — the real archive comes out
        /// this size or smaller. Pass `tsv` for a two-line machine-readable
        /// table (`file`, `input_bytes`, `est_fqxv_bytes`, `ratio`) instead of
        /// the human report.
        #[arg(long, value_enum, num_args = 0..=1, default_missing_value = "human", value_name = "FORMAT")]
        estimate: Option<EstimateFormat>,
        /// Compression effort level (1-9); higher raises the sequence order.
        ///
        /// The primary effort knob when you aren't using `--max`. Higher levels
        /// raise the sequence context order (and, at the top levels, enable a
        /// hashed high-order model) for a better ratio at more time and memory.
        #[arg(short, long, default_value_t = 5, help_heading = "Advanced")]
        level: u8,
        /// Opt-in lossy quality binning (default: lossless).
        ///
        /// Changes the data. Quantizes quality scores to a fixed set of levels
        /// (Illumina 8/4/2-bin, or the Nanopore/HiFi CoLoRd cutpoints) for a
        /// smaller quality stream at the cost of exact fidelity.
        #[arg(long, value_enum, default_value_t = QualityBin::Lossless, help_heading = "Advanced")]
        quality_bin: QualityBin,
        /// Read-order guarantee: `preserve` (default) or `any`.
        ///
        /// `preserve` restores the original order on decompress. `any` allows
        /// read reordering to exploit depth redundancy for a better ratio —
        /// single-end reads may come back in a different order; paired/grouped
        /// input still round-trips in order (the mate interleaving requires it),
        /// it just costs a stored permutation.
        #[arg(long, value_enum, default_value_t = ReadOrder::Preserve, help_heading = "Advanced")]
        order: ReadOrder,
        /// With `--order any`, force original read order to be restored.
        ///
        /// Stores a permutation and codes names/quality in original order. By
        /// default single-end reorder picks this automatically when it makes the
        /// archive smaller — counter-style names (e.g. an incrementing `.N N`)
        /// delta-code to almost nothing in original order, so the permutation
        /// beats a scrambled-name stream. Pass this to force it on.
        #[arg(long, help_heading = "Advanced")]
        keep_order: bool,
        /// With `--order any`, disable the adaptive assembly codecs.
        ///
        /// Uses the faster single-contig sequence codec only. By default reorder
        /// codes each block with every codec and keeps the smaller — plus, when
        /// a shared global reference nets a whole-file win, stores it once and
        /// codes reads as positions on it (never worse than block-local) —
        /// recovering reads the single-contig codec would strand as literals at
        /// a higher encode cost. No effect without `--order any`.
        #[arg(long = "no-rescue", help_heading = "Advanced")]
        no_rescue: bool,
        /// Reads per row group (block); overrides the level's block size.
        ///
        /// Decouples random-access granularity from effort: smaller groups give
        /// finer random access and more parallelism at some ratio cost (the
        /// order-k sequence model trains on fewer reads), larger groups the
        /// reverse. Useful when archiving to object storage where you want to
        /// fetch small ranges. Sequence order still follows `--level`. Ignored by
        /// the reorder path (`--order any`/`--max`), which clusters globally.
        #[arg(long, value_name = "N", help_heading = "Advanced")]
        block_reads: Option<usize>,
        /// Force single-input interleaving, in members per spot.
        ///
        /// Auto-detected from read names by default; pass to force (1 =
        /// single-end, 2 = paired as from `sracha get -Z`). Ignored with
        /// multiple inputs. A group size of 0 is meaningless and rejected.
        #[arg(
            long,
            value_name = "N",
            help_heading = "Advanced",
            value_parser = clap::value_parser!(u8).range(1..)
        )]
        interleaved: Option<u8>,
        /// Sequencing platform to record (default: auto-detect).
        ///
        /// Auto-detected from the read names by default; pass to force it (e.g.
        /// for unusual name conventions the detector doesn't recognize).
        #[arg(long, value_enum, help_heading = "Advanced")]
        platform: Option<Platform>,
    },
    /// Decompress a `.fqxv` file to FASTQ.
    ///
    /// Writes interleaved FASTQ by default; `--split` restores separate
    /// per-spot mate files and `-Z` streams to stdout.
    ///
    /// The input may be `-` for stdin: the archive is streamed and decoded on the
    /// fly, so a remote archive is read by piping a transfer tool in — e.g.
    /// `aws s3 cp s3://bkt/reads.fqxv - | fqxv decompress - -Z | bwa mem ref.fa -`
    /// (or `curl`/`gsutil cat`) aligns reads as they arrive without staging the file
    /// to disk, and the transfer tool handles all the auth, retries, and resume.
    /// (`--recover` and `--split` still need a seekable file — a stream can't be
    /// rewound.)
    Decompress {
        /// Input `.fqxv` file, or `-` to stream from stdin.
        input: PathBuf,
        /// Interleaved FASTQ output file (`.gz` => BGZF).
        ///
        /// A `.gz` extension writes block-gzip (BGZF); any other extension
        /// writes plain FASTQ. Use `-Z/--stdout` (or `-o -`) to stream to stdout
        /// instead.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Stream interleaved FASTQ (always raw) to stdout.
        ///
        /// E.g. to pipe into an aligner:
        /// `fqxv decompress x.fqxv -Z | bwa mem -p ref -`. Required to write to
        /// stdout — a bare `decompress` with no `-o`/`--split`/`-Z` errors rather
        /// than flooding the terminal with reads.
        #[arg(short = 'Z', long, conflicts_with_all = ["output", "split"])]
        stdout: bool,
        /// Restore separate per-spot mate files under `<prefix>`.
        ///
        /// Writes `<prefix>_R1.fastq.gz … _R<G>.fastq.gz` (block-gzip by default;
        /// see `--mate-style`/`--no-gzip`).
        #[arg(long, conflicts_with = "output")]
        split: Option<PathBuf>,
        /// Overwrite output FASTQ file(s) if they already exist.
        ///
        /// Without this, an existing output (`-o FILE` or any `--split` mate
        /// file) is left untouched and the command errors before decoding.
        /// Ignored when writing to stdout (`-Z`/`-o -`).
        #[arg(short = 'f', long)]
        force: bool,
        /// Labeling for `--split` outputs: `r` (`_R1`) or `num` (`_1`).
        ///
        /// `r` gives `_R1`,`_R2`,… (Illumina convention, the default); `num`
        /// gives `_1`,`_2`,….
        #[arg(long, value_enum, default_value_t = MateStyle::R, requires = "split", help_heading = "Advanced")]
        mate_style: MateStyle,
        /// Write plain `.fastq` for `--split` (default is `.fastq.gz`).
        ///
        /// For `-o FILE`, compression follows the file extension instead.
        #[arg(long, requires = "split", help_heading = "Advanced")]
        no_gzip: bool,
        /// Best-effort recovery of a corrupted archive.
        ///
        /// Skip blocks that fail their CRC and decode the rest, reporting what
        /// was lost. Interleaved output only (plain, non-reordered archives).
        #[arg(long, conflicts_with = "split", help_heading = "Advanced")]
        recover: bool,
    },
    /// Print `.fqxv` container metadata and per-stream sizes.
    ///
    /// Give several files or a directory (scanned recursively for `*.fqxv`) to
    /// report a batch — one TSV/JSON entry per archive, keyed by filename.
    Info {
        /// Input `.fqxv` file(s), or a directory scanned for `*.fqxv`.
        #[arg(num_args = 1..)]
        inputs: Vec<PathBuf>,
        /// Also report content statistics (decodes the whole archive).
        ///
        /// Read-length spread, base composition, GC%, and the quality
        /// distribution. Unlike the default metadata-only report (which just
        /// reads the header and footer index), this costs a full decompress.
        #[arg(long, short = 's')]
        stats: bool,
        /// Emit machine-readable TSV instead of the human report.
        ///
        /// A single file keeps the stable one-row columns (benchmark harness /
        /// scripts); a batch prepends a `file` column and prints one row per
        /// archive.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit a JSON object instead of the human report.
        #[arg(long)]
        json: bool,
    },
    /// Verify archive integrity (CRC checks) without writing output.
    ///
    /// Prints a table of per-check results; exits non-zero if any archive is
    /// corrupt. Give several files or a directory (scanned recursively for
    /// `*.fqxv`) to verify a batch.
    Verify {
        /// Input `.fqxv` file(s), or a directory scanned for `*.fqxv`.
        #[arg(num_args = 1..)]
        inputs: Vec<PathBuf>,
        /// Faster, weaker check: per-block CRC via the footer index.
        ///
        /// Verifies each block's stored CRC through the footer index (parallel
        /// positioned reads) instead of the whole-file digest. Skips the header,
        /// footer, and inter-block framing bytes.
        #[arg(long)]
        quick: bool,
        /// Emit tab-separated per-check rows instead of the human table.
        ///
        /// A batch prepends a `file` column so rows from different archives stay
        /// distinct.
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

/// Map a 1-9 effort level to the multi-reference tiler's `(band, max_refs)`
/// (Nanopore long reads only — the only path that runs the tiler). Best-of-N
/// reference selection is the dominant ONT sequence-ratio lever, so it ramps with
/// effort: the default levels keep the cheap single-reference cover (byte-identical
/// to what ships without this knob), the upper levels turn on best-of-N, and
/// `--max` (level 9) reaches the CoLoRd-parity operating point — measured
/// band 768 / best-of-4 = ~0.89 b/base on `ecoli_ont`, below CoLoRd's ~0.91, at
/// ~6× the tiler's alignment time. Above AVX2 the alignment is already vectorised,
/// so the cost is purely the extra candidate alignments the user opted into.
fn level_to_tile(level: u8) -> (usize, usize) {
    match level {
        0..=6 => (256, 1), // codec default: single reference, current ONT speed
        7 => (256, 2),     // best-of-2: the cheap knee (~-14% for ~1.5× tiler time)
        8 => (256, 4),     // best-of-4 (~-19%)
        _ => (768, 4),     // level 9 / --max: wide band + best-of-4 → CoLoRd parity
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
    // Clean up an in-progress temp archive if the user Ctrl-C's mid-compress, so
    // an interrupted stream never orphans a partial file (the atomic rename
    // already guarantees no corrupt archive lands at the destination path).
    output::install_interrupt_handler();
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
            let (tile_band, tile_max_refs) = level_to_tile(level);
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
                tile_band,
                tile_max_refs,
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
            // Write to a sibling temp file and rename into place only once the
            // whole archive is on disk. An interrupted or failed run then never
            // leaves a footer-less, corrupt `.fqxv` at the destination path.
            let (archive, out) = output::AtomicOutput::create(&output, force)?;

            // The reorder path reads the whole input up front before its
            // single-threaded compute begins, so input-byte progress would stall
            // at 100% mid-run — give it an indeterminate readout and name it after
            // the slow work. The streaming path tracks input consumption, so it
            // gets a live bar (a percentage when the input size is known, a
            // bytes + rate readout for a stdin stream of unknown length).
            let spin_msg = if params.reorder {
                "compressing (reorder — this can take a while)…"
            } else {
                "compressing…"
            };

            let t0 = Instant::now();
            let stats = if inputs.len() == 1 {
                // Create the byte counter before opening the reader so the
                // gzip-sniffing read — which blocks on a not-yet-ready stdin — is
                // already counted and the indicator is live while we wait on the
                // upstream producer (a pause reads as `0 B`, not a hang).
                let counter = Arc::new(AtomicU64::new(0));
                let total = single_input_len(&inputs[0], params.reorder);
                let reader = open_input_with_counter(&inputs[0], Arc::clone(&counter))?;
                let bar = progress::Bar::start(spin_msg, cli.quiet, counter, total);
                let stats = match interleaved {
                    None => fqxv::compress_auto(reader, out, params)?,
                    Some(g) if g > 1 => fqxv::compress_interleaved(reader, out, params, g)?,
                    Some(_) => fqxv::compress(reader, out, params)?,
                };
                bar.abandon();
                stats
            } else {
                // Interleaving reads every input in lockstep, so one shared counter
                // over all of them against their summed length is the same measure
                // the single-input bar uses — paired mates get a real percentage
                // rather than the bare spinner this path used to fall back to.
                let counter = Arc::new(AtomicU64::new(0));
                let total = multi_input_len(&inputs, params.reorder);
                let readers: Vec<Box<dyn Read + Send>> = inputs
                    .iter()
                    .map(|p| open_input_with_counter(p, Arc::clone(&counter)))
                    .collect::<anyhow::Result<_>>()?;
                let bar = progress::Bar::start(spin_msg, cli.quiet, Arc::clone(&counter), total);
                let stats = fqxv::compress_multi(readers, out, params)?;
                bar.abandon();
                stats
            };
            let secs = t0.elapsed().as_secs_f64();

            // Optional read-after-write verification: fully decode the archive we
            // just wrote and confirm it round-trips *before* publishing it. We
            // verify the temp file and only commit (atomically rename into place)
            // on success, so `--verify` never leaves a bad archive at the
            // destination. On a mismatch or decode error the temp is kept off the
            // final path for inspection and we bail non-zero. On a `?` bail
            // elsewhere the spinner's Drop clears its line before the error prints.
            if verify {
                let sp = progress::Spinner::start("verifying…", cli.quiet);
                let verified = File::open(archive.temp_path())
                    .with_context(|| {
                        format!("reopening {} to verify", archive.temp_path().display())
                    })
                    .and_then(|f| {
                        fqxv::verify_roundtrip(f, cli.threads).context("verifying the archive")
                    });
                match verified {
                    Ok(decoded) if decoded == stats.reads => sp.abandon(),
                    Ok(decoded) => {
                        let kept = archive.keep_for_inspection();
                        anyhow::bail!(
                            "verification decoded {decoded} reads but {} were written; \
                             unverified output kept at {} (not published to {})",
                            stats.reads,
                            kept.display(),
                            output.display()
                        );
                    }
                    Err(e) => {
                        let kept = archive.keep_for_inspection();
                        return Err(e.context(format!(
                            "unverified output kept at {} (not published to {})",
                            kept.display(),
                            output.display()
                        )));
                    }
                }
            }
            // Publish: atomically move the finished, verified archive into place.
            archive.commit()?;

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
            // Input is a local path or stdin (`-`). stdin streams the archive body
            // straight into the decoder (`decompress` only needs `Read`) — the way
            // to read from S3/HTTP is to pipe a transfer tool in (see `is_stdin`).
            // A file is seekable, so `open_in` can reopen it (for the footer
            // pre-check and the `--split` header peek); stdin can only be consumed
            // once, so those paths are gated off below.
            let stdin_input = is_stdin(&input);
            let open_in = || -> anyhow::Result<Box<dyn Read + Send>> {
                if stdin_input {
                    Ok(Box::new(io::stdin()))
                } else {
                    Ok(Box::new(File::open(&input).with_context(|| {
                        format!("opening input {}", input.display())
                    })?))
                }
            };
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
                // Recovery rescans for block-sync markers, which needs to seek —
                // stdin can't be rewound. Write the stream to a file first if you
                // need `--recover`.
                if stdin_input {
                    anyhow::bail!(
                        "--recover needs a seekable file; stdin can't be rewound to rescan \
                         for block boundaries — write the stream to a file first"
                    );
                }
                // Recovery deliberately tolerates missing/corrupt blocks, so the
                // completeness check below does not apply here.
                let sp = progress::Spinner::start("recovering…", cli.quiet);
                let (pending, mut sink) = open_sink(output.as_deref(), stdout, force)?;
                let archive = File::open(&input)
                    .with_context(|| format!("opening input {}", input.display()))?;
                let rec = fqxv::decompress_recover(archive, &mut sink, cli.threads)?;
                sink.finish()?;
                // Recovery's output is deliberately incomplete, but it is still a
                // *finished* artifact — publish it only once the salvage is done,
                // so an interrupt mid-recover doesn't masquerade as one.
                if let Some(pending) = pending {
                    pending.commit()?;
                }
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
                // A `--split` decode needs the group size before decoding (to open
                // the output files), which means a header peek and then a second
                // read from the start — impossible on a single-pass stream. Write
                // the stream to a file, or decode interleaved (`-Z`) and split
                // downstream.
                if stdin_input && split.is_some() {
                    anyhow::bail!(
                        "--split can't read from stdin (it needs a second pass for the header); \
                         write the stream to a file first, or decode interleaved with -Z"
                    );
                }
                // Read the footer's authoritative read count first. For a truncated
                // archive (which also lost its footer) this errors before any output
                // is written; otherwise we compare it against what we actually
                // decode so a silently short file is caught rather than trusted.
                // stdin can't seek to the footer, so it skips this: the terminator
                // frame marks a clean end and per-block CRCs catch a truncated
                // stream (decode hits EOF before the terminator).
                let expected = if stdin_input {
                    None
                } else {
                    fqxv::expected_reads(
                        File::open(&input)
                            .with_context(|| format!("opening input {}", input.display()))?,
                    )?
                };
                // Reads written drive the readout, measured against the footer's
                // count. Archive bytes consumed would be the easier counter but a
                // misleading one: decode pulls blocks in batches of `--threads`, so
                // input consumption runs ahead of the work and pins at 100% while
                // the last batch is still decoding. Output records track what is
                // actually finished. A global-reorder archive has no footer count,
                // and degrades to a reads-so-far readout.
                let done = Arc::new(AtomicU64::new(0));
                let lines = Arc::new(AtomicU64::new(0));
                let bar = progress::Bar::start_with_unit(
                    "decompressing…",
                    cli.quiet,
                    Arc::clone(&done),
                    expected,
                    progress::Unit::Reads,
                );
                // Commits are deferred until after the completeness check below, so
                // a truncated archive leaves nothing at the destination names.
                let (stats, pending) = if let Some(prefix) = split {
                    let g = fqxv::peek(open_in()?)?.group_size as usize;
                    let (pending, sinks): (Vec<_>, Vec<FastqSink>) = (1..=g)
                        .map(|i| {
                            open_fastq_file(&split_path(&prefix, i, mate_style, !no_gzip), force)
                        })
                        .collect::<anyhow::Result<Vec<_>>>()?
                        .into_iter()
                        .unzip();
                    let mut sinks: Vec<CountingSink> = sinks
                        .into_iter()
                        .map(|s| CountingSink::new(s, Arc::clone(&lines), Arc::clone(&done)))
                        .collect();
                    let stats = fqxv::decompress_split(open_in()?, &mut sinks, cli.threads)?;
                    for sink in sinks {
                        sink.into_inner().finish()?;
                    }
                    (stats, pending)
                } else {
                    let (pending, sink) = open_sink(output.as_deref(), stdout, force)?;
                    let mut sink = CountingSink::new(sink, Arc::clone(&lines), Arc::clone(&done));
                    let stats = fqxv::decompress(open_in()?, &mut sink, cli.threads)?;
                    sink.into_inner().finish()?;
                    (stats, pending.into_iter().collect())
                };
                if let Some(expected) = expected
                    && stats.reads != expected
                {
                    anyhow::bail!(
                        "archive is truncated: footer declares {expected} reads but only \
                             {} were decoded — the output is incomplete",
                        stats.reads
                    );
                }
                for out in pending {
                    out.commit()?;
                }
                let secs = t0.elapsed().as_secs_f64();
                bar.abandon();
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
        } => print_verify(&inputs, quick, tsv, json, cli.threads)?,
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
                .filter(|&(_, &count)| count > 0)
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
fn print_verify(
    inputs: &[PathBuf],
    quick: bool,
    tsv: bool,
    json: bool,
    threads: usize,
) -> anyhow::Result<()> {
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
                fqxv::verify_report(&f, quick, threads)
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

/// Whether the decompress input is stdin (`-`). A stream can only be read once and
/// can't seek, so the caller streams it straight through the decoder. Remote inputs
/// are handled by piping a transfer tool into stdin — e.g.
/// `aws s3 cp s3://bucket/reads.fqxv - | fqxv decompress - -Z | bwa mem ref.fa -` —
/// which keeps all the auth/retry/resume in the tool built for it and adds no HTTP
/// dependency here.
fn is_stdin(input: &Path) -> bool {
    input.as_os_str() == "-"
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
/// source into `count`. The compress progress bar drives a live readout from
/// this counter; passing the counter in (rather than getting it back) lets the
/// caller start the indicator *before* the reader's gzip-sniffing first read,
/// which blocks on a not-yet-ready stdin. For a gzip input the count is on-disk
/// (compressed) bytes, matching the file length the bar measures against.
fn open_input_with_counter(
    path: &Path,
    count: Arc<AtomicU64>,
) -> anyhow::Result<Box<dyn Read + Send>> {
    let counted: Box<dyn Read + Send> = Box::new(CountingReader {
        inner: raw_input(path)?,
        count,
    });
    decode_input(counted)
}

/// [`open_input_with_counter`] with a fresh counter returned alongside. `--estimate`
/// uses this to learn how much on-disk input a sample consumed, so it can project
/// the whole-file archive size.
fn open_input_counted(path: &Path) -> anyhow::Result<(Box<dyn Read + Send>, Arc<AtomicU64>)> {
    let count = Arc::new(AtomicU64::new(0));
    let reader = open_input_with_counter(path, Arc::clone(&count))?;
    Ok((reader, count))
}

/// The on-disk byte length to drive a compress progress bar, or `None` for an
/// indeterminate readout. `None` when the input is stdin (`-`, no length) or
/// when `reorder` is set — the reorder codec reads the whole input up front, so
/// an input-byte bar would misleadingly stall at 100% during its compute.
fn single_input_len(path: &Path, reorder: bool) -> Option<u64> {
    if reorder || path.as_os_str() == "-" {
        return None;
    }
    std::fs::metadata(path).ok().map(|m| m.len())
}

/// [`single_input_len`] for an interleaved multi-input compress: the summed
/// on-disk length of every input, or `None` for an indeterminate readout. The
/// `reorder` and stdin carve-outs apply for the same reasons, and a single
/// unstattable input makes the whole total meaningless — a bar measuring
/// against a short total would run past 100%, so fall back rather than lie.
fn multi_input_len(paths: &[PathBuf], reorder: bool) -> Option<u64> {
    if reorder {
        return None;
    }
    paths
        .iter()
        .map(|p| single_input_len(p, false))
        .sum::<Option<u64>>()
}

/// A pass-through FASTQ sink that counts the records passing through it, so a
/// decompress bar can report reads finished against the footer's read count.
///
/// Counts newlines and divides by four rather than parsing: fqxv emits exactly
/// four lines per record (the `+` line is normalized to a bare `+`), so the
/// division is exact for what we write, and a byte scan costs far less than the
/// decode it is measuring.
///
/// Both counters are shared across every sink of a `--split` decode — the bar
/// wants one total, not one per mate. The running value is approximate between
/// updates (two sinks can interleave a fetch_add and a store), which is all a
/// progress readout needs; the authoritative count is `Stats::reads`.
struct CountingSink {
    inner: FastqSink,
    /// Newlines seen, the raw quantity this can count incrementally.
    lines: Arc<AtomicU64>,
    /// Records finished — `lines / 4`, kept separately because the bar reads a
    /// counter directly and cannot divide.
    reads: Arc<AtomicU64>,
}

impl CountingSink {
    fn new(inner: FastqSink, lines: Arc<AtomicU64>, reads: Arc<AtomicU64>) -> Self {
        Self {
            inner,
            lines,
            reads,
        }
    }

    fn into_inner(self) -> FastqSink {
        self.inner
    }
}

impl Write for CountingSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        let newlines = buf[..n].iter().filter(|&&b| b == b'\n').count() as u64;
        if newlines > 0 {
            let total = self.lines.fetch_add(newlines, Ordering::Relaxed) + newlines;
            self.reads.store(total / 4, Ordering::Relaxed);
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
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

    // Empty input: no reads to sample. Report it cleanly (0 reads / 0 bytes)
    // rather than projecting a "0x" ratio or mislabeling a file as a streaming
    // input. Mirrors `compress`, which accepts empty input and writes an empty
    // archive.
    if reads == 0 {
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
        match fmt {
            EstimateFormat::Tsv => {
                println!("file\tinput_bytes\test_fqxv_bytes\tratio");
                println!("{file}\t{disk_in}\t0\t0.0000");
            }
            EstimateFormat::Human => {
                println!("{file}: input has no reads — nothing to estimate");
            }
        }
        return Ok(());
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
/// Returns the sink together with a pending commit when the destination is a
/// file; stdout has no destination name to protect, so it commits to nothing.
fn open_sink(
    output: Option<&Path>,
    to_stdout: bool,
    force: bool,
) -> anyhow::Result<(Option<output::AtomicOutput>, FastqSink)> {
    match output {
        Some(p) if p.as_os_str() == "-" => Ok((None, FastqSink::Raw(Box::new(io::stdout())))),
        Some(p) => open_fastq_file(p, force).map(|(pending, sink)| (Some(pending), sink)),
        None if to_stdout => Ok((None, FastqSink::Raw(Box::new(io::stdout())))),
        None => anyhow::bail!(
            "no output specified: pass -o FILE (\".gz\" for BGZF), --split PREFIX for separate \
             mate files, or -Z/--stdout to stream interleaved FASTQ to stdout"
        ),
    }
}

/// Create a FASTQ output file, choosing BGZF compression from its `.gz` extension.
///
/// The bytes go to a sibling temp that reaches `path` only when the returned
/// [`output::AtomicOutput`] is committed, so an interrupted decode leaves no
/// plausible-looking partial FASTQ at the destination name. The caller must
/// finish the sink (flushing BGZF's EOF block) *before* committing.
fn open_fastq_file(path: &Path, force: bool) -> anyhow::Result<(output::AtomicOutput, FastqSink)> {
    let (pending, file) = output::AtomicOutput::create(path, force)?;
    Ok((pending, FastqSink::new(Box::new(file), is_gzip_path(path))))
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
    fn fastq_output_refuses_existing_unless_forced() {
        let dir = std::env::temp_dir().join(format!("fqxv-create-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("reads.fastq");

        // A fresh path is accepted regardless of `force`.
        let (pending, sink) = open_fastq_file(&path, false).unwrap();
        sink.finish().unwrap();
        pending.commit().unwrap();
        assert!(path.exists());
        std::fs::write(&path, b"existing contents").unwrap();

        // Without force, an existing path errors and is left untouched.
        // `unwrap_err` would need `Debug` on the sink, which wraps a `dyn Write`.
        let err = match open_fastq_file(&path, false) {
            Ok(_) => panic!("expected an error for an existing output"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("already exists"), "unexpected error: {err}");
        assert_eq!(std::fs::read(&path).unwrap(), b"existing contents");

        // With force, the destination keeps its old contents until the commit
        // swaps the finished output in — a decode that dies partway leaves the
        // previous file intact rather than a half-overwritten one.
        let (pending, mut sink) = open_fastq_file(&path, true).unwrap();
        sink.write_all(b"replacement").unwrap();
        sink.finish().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"existing contents");
        pending.commit().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"replacement");

        std::fs::remove_dir_all(&dir).ok();
    }

    /// The wart this guards: an interrupted `--split` decode used to leave
    /// mate files at their real names, and a truncated `.fastq.gz` is
    /// indistinguishable from a complete one. Dropping the pending output
    /// (what a `?` bail does; the signal handler does the same via the ACTIVE
    /// list) must leave *nothing* at the destination.
    #[test]
    fn abandoned_fastq_output_never_reaches_the_destination() {
        let dir = std::env::temp_dir().join(format!("fqxv-abandon-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("sample_R1.fastq.gz");

        let (pending, mut sink) = open_fastq_file(&path, false).unwrap();
        sink.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        sink.flush().unwrap();
        // No `finish`/`commit` — models an interrupt mid-decode.
        drop(sink);
        drop(pending);

        assert!(
            !path.exists(),
            "a partial mate file must never appear at the destination name"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    /// Whatever `build.rs` resolves from git, the crate version must lead —
    /// scripts and humans key off that prefix, and a git description alone
    /// (or an empty string from a failed command) would break them.
    #[test]
    fn version_always_leads_with_the_crate_version() {
        let pkg = env!("CARGO_PKG_VERSION");
        assert!(
            VERSION.starts_with(pkg),
            "version {VERSION:?} must start with the crate version {pkg:?}"
        );
        // Either the bare version (a clean release tag) or it plus a
        // parenthesized git description — never anything else.
        let suffix = &VERSION[pkg.len()..];
        assert!(
            suffix.is_empty() || (suffix.starts_with(" (") && suffix.ends_with(')')),
            "unexpected version suffix {suffix:?}"
        );
    }

    #[test]
    fn resolve_threads_zero_is_all_cores_and_explicit_is_clamped() {
        let available = std::thread::available_parallelism().map_or(1, |n| n.get());
        assert_eq!(resolve_threads(0), available); // 0 = all cores
        assert_eq!(resolve_threads(1), 1);
        assert_eq!(resolve_threads(usize::MAX), available); // never oversubscribe
    }
}
