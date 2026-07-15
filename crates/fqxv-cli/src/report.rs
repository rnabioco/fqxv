//! Human-facing run summaries and the formatting/styling they rely on.
//!
//! These are printed to **stderr** (stdout is reserved for piped FASTQ) directly
//! via `eprintln!`, gated only on `--quiet` — deliberately decoupled from the
//! `tracing` log level so `-v` adds diagnostics without reshaping the summary and
//! `--quiet` is the single switch that hides it. Color is applied only when
//! [`set_color`] was told stderr is an interactive terminal with color allowed.

use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use fqxv::Stats;

/// Whether summary output should be colorized. Set once at startup by
/// [`set_color`]; read by the `style::*` helpers.
static COLOR: AtomicBool = AtomicBool::new(false);

/// Enable or disable ANSI color in summaries. Call once at startup with the
/// result of "not quiet, stderr is a TTY, and `NO_COLOR` is unset".
pub(crate) fn set_color(on: bool) {
    COLOR.store(on, Ordering::Relaxed);
}

fn colored() -> bool {
    COLOR.load(Ordering::Relaxed)
}

/// Semantic styling for summary text — mirrors the palette used by the sibling
/// `sracha` tooling so the rnabioco CLIs read as one family.
mod style {
    use super::colored;
    use owo_colors::OwoColorize;
    use std::fmt::Display;

    /// Section headers and file names: bold.
    pub(crate) fn header<T: Display>(s: T) -> String {
        if colored() {
            format!("{}", s.bold())
        } else {
            format!("{s}")
        }
    }

    /// Important counts (reads, blocks): green.
    pub(crate) fn count<T: Display>(n: T) -> String {
        if colored() {
            format!("{}", n.green())
        } else {
            format!("{n}")
        }
    }

    /// Values in key/value context (sizes, ratios, rates): cyan.
    pub(crate) fn value<T: Display>(s: T) -> String {
        if colored() {
            format!("{}", s.cyan())
        } else {
            format!("{s}")
        }
    }

    /// File paths: cyan.
    pub(crate) fn path<T: Display>(s: T) -> String {
        value(s)
    }
}

/// Format a byte count as a human-readable string using decimal SI units
/// (e.g. `276.15 MB`), matching how disk/file sizes are conventionally reported
/// to users. `KB`/`MB`/`GB`/`TB` are powers of 1000.
pub(crate) fn format_size(bytes: u64) -> String {
    const KB: u64 = 1_000;
    const MB: u64 = 1_000_000;
    const GB: u64 = 1_000_000_000;
    const TB: u64 = 1_000_000_000_000;
    if bytes >= TB {
        format!("{:.2} TB", bytes as f64 / TB as f64)
    } else if bytes >= GB {
        format!("{:.2} GB", bytes as f64 / GB as f64)
    } else if bytes >= MB {
        format!("{:.2} MB", bytes as f64 / MB as f64)
    } else if bytes >= KB {
        format!("{:.2} KB", bytes as f64 / KB as f64)
    } else {
        format!("{bytes} B")
    }
}

/// Insert thousands separators into an integer (e.g. `1234567` -> `1,234,567`).
pub(crate) fn thousands(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

/// Throughput in MB/s (decimal) over `bytes` processed in `secs`, or `None` when
/// there is nothing to report (no bytes, or an immeasurably short run).
fn mb_per_s(bytes: u64, secs: f64) -> Option<f64> {
    if bytes == 0 || secs <= 0.0 {
        return None;
    }
    Some(bytes as f64 / 1_000_000.0 / secs)
}

/// Human label for a `group_size` interleaving.
fn layout_label(group_size: u8) -> String {
    match group_size {
        0 | 1 => "single-end".to_string(),
        2 => "paired".to_string(),
        g => format!("grouped x{g}"),
    }
}

/// Build the multi-line compress summary. `in_size` is 0 for a stdin (`-`) input,
/// where the input size isn't known up front — the before/ratio figures are then
/// omitted, but the resulting archive size is always shown.
pub(crate) fn compress_summary(
    inputs: &[String],
    output: &Path,
    in_size: u64,
    stats: &Stats,
    secs: f64,
) -> String {
    let src = inputs.join(&style::header(" + "));
    let mut out = format!(
        "{} {src} {} {}",
        style::header("compressed"),
        style::header("->"),
        style::path(output.display()),
    );
    // Size line: before -> after (ratio, % saved) whenever the input size is
    // known; otherwise just the archive size. This is the headline result, so it
    // comes before the read/block/timing detail.
    if in_size > 0 && stats.out_bytes > 0 {
        let ratio = in_size as f64 / stats.out_bytes as f64;
        let saved = 100.0 * (1.0 - stats.out_bytes as f64 / in_size as f64);
        out.push_str(&format!(
            "\n  {} {} {}  ({}, {} saved)",
            style::value(format_size(in_size)),
            style::header("->"),
            style::value(format_size(stats.out_bytes)),
            style::value(format!("{ratio:.2}x")),
            style::value(format!("{saved:.1}%")),
        ));
    } else if stats.out_bytes > 0 {
        out.push_str(&format!(
            "\n  archive {}",
            style::value(format_size(stats.out_bytes))
        ));
    }
    // Detail line: layout, counts, wall time, and throughput (measured against
    // the input when known, else the archive written).
    let mut detail = format!(
        "\n  {} · {} reads · {} block(s) in {}",
        layout_label(stats.group_size),
        style::count(thousands(stats.reads)),
        style::count(thousands(stats.blocks)),
        style::value(format!("{secs:.1}s")),
    );
    let rate_bytes = if in_size > 0 {
        in_size
    } else {
        stats.out_bytes
    };
    if let Some(rate) = mb_per_s(rate_bytes, secs) {
        detail.push_str(&format!("  ·  {}", style::value(format!("{rate:.1} MB/s"))));
    }
    out.push_str(&detail);
    out
}

/// Build the multi-line decompress summary. `dest` describes where the FASTQ went
/// (a path, a `--split` prefix, or `stdout`).
pub(crate) fn decompress_summary(input: &Path, dest: &str, stats: &Stats, secs: f64) -> String {
    let mut out = format!(
        "{} {} {} {}",
        style::header("decompressed"),
        style::path(input.display()),
        style::header("->"),
        style::path(dest),
    );
    let mut line = format!(
        "\n  {} reads · {} block(s) in {}",
        style::count(thousands(stats.reads)),
        style::count(thousands(stats.blocks)),
        style::value(format!("{secs:.1}s")),
    );
    if stats.out_bytes > 0 {
        line.push_str(&format!(
            "  ·  {}",
            style::value(format_size(stats.out_bytes))
        ));
        if let Some(rate) = mb_per_s(stats.out_bytes, secs) {
            line.push_str(&format!(" ({})", style::value(format!("{rate:.1} MB/s"))));
        }
    }
    out.push_str(&line);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_size_scales() {
        assert_eq!(format_size(0), "0 B");
        assert_eq!(format_size(999), "999 B");
        assert_eq!(format_size(1_000), "1.00 KB");
        assert_eq!(format_size(1_000_000), "1.00 MB");
        assert_eq!(format_size(1_000_000_000), "1.00 GB");
        assert_eq!(format_size(1_000_000_000_000), "1.00 TB");
    }

    #[test]
    fn thousands_basics() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(1_000), "1,000");
        assert_eq!(thousands(1_234_567), "1,234,567");
    }

    #[test]
    fn mb_per_s_guards_zero() {
        assert!(mb_per_s(0, 1.0).is_none());
        assert!(mb_per_s(1_000_000, 0.0).is_none());
        assert!(mb_per_s(1_000_000, 1.0).is_some());
    }

    #[test]
    fn compress_summary_always_shows_archive_size() {
        // stdin input (in_size == 0): no ratio, but the archive size is present.
        let stats = Stats {
            reads: 5,
            blocks: 1,
            out_bytes: 2_000_000,
            group_size: 1,
        };
        let s = compress_summary(&["-".to_string()], Path::new("out.fqxv"), 0, &stats, 1.0);
        assert!(s.contains("2.00 MB"));
        assert!(!s.contains("saved"));
    }

    #[test]
    fn compress_summary_omits_ratio_for_stdin() {
        let stats = Stats {
            reads: 10,
            blocks: 1,
            out_bytes: 100,
            group_size: 2,
        };
        let s = compress_summary(&["-".to_string()], Path::new("out.fqxv"), 0, &stats, 1.0);
        assert!(s.contains("paired"));
        assert!(!s.contains("saved"));
    }

    #[test]
    fn compress_summary_reports_ratio_with_input_size() {
        let stats = Stats {
            reads: 10,
            blocks: 1,
            out_bytes: 250,
            group_size: 1,
        };
        let s = compress_summary(
            &["reads.fastq".to_string()],
            Path::new("reads.fqxv"),
            1000,
            &stats,
            1.0,
        );
        assert!(s.contains("4.00x"));
        assert!(s.contains("saved"));
    }
}
