//! Terminal progress for the slow phases of compress/decompress.
//!
//! Two shapes, both on **stderr** (so decoded FASTQ on stdout stays clean):
//!
//! * [`Spinner`] — an indeterminate ticking line for phases with no measurable
//!   progress (verify, decompress, recover, or the reorder compress path, whose
//!   work isn't tracked by input consumption).
//! * [`Bar`] — a live readout of a compress driven by a bytes-read counter. When
//!   the input size is known (a file) it renders a real percentage bar; for a
//!   stream of unknown length (stdin) it degrades to a bytes + rate spinner, so
//!   a pause waiting on an upstream producer reads as `0 B` rather than a hang.
//!
//! On a non-TTY (CI logs, redirects) both print the start message as a plain
//! line and do no live drawing, so logs stay free of control characters. Under
//! `--quiet` both are fully inert. If a handle is dropped without an explicit
//! `finish`/`abandon` — e.g. a `?` bails mid-run — [`Drop`] clears the live line
//! so the error message isn't tangled with a half-drawn indicator.

use std::io::IsTerminal;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

/// Braille tick frames, matching the rnabioco tooling look.
const TICK_STRINGS: &[&str] = &["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"];

/// How often the monitor thread samples the byte counter and redraws.
const REFRESH: Duration = Duration::from_millis(120);

/// A single-line indeterminate status indicator (TTY-animated, non-TTY plain,
/// quiet-inert).
#[derive(Debug)]
pub(crate) struct Spinner {
    /// `Some` while a live indicatif bar is running; taken by
    /// finish/abandon/drop so exactly one of them acts.
    pb: Option<ProgressBar>,
}

impl Spinner {
    /// Start a spinner with `message`. `quiet` silences it entirely; otherwise it
    /// animates on an interactive stderr and degrades to a printed line on a
    /// non-TTY.
    pub(crate) fn start(message: impl Into<String>, quiet: bool) -> Self {
        if quiet {
            return Self { pb: None };
        }
        let message = message.into();
        if std::io::stderr().is_terminal() {
            let pb = ProgressBar::new_spinner();
            pb.set_style(
                ProgressStyle::with_template("{spinner:.cyan} {msg}")
                    .expect("valid spinner template")
                    .tick_strings(TICK_STRINGS),
            );
            pb.set_message(message);
            pb.enable_steady_tick(Duration::from_millis(100));
            Self { pb: Some(pb) }
        } else {
            eprintln!("{message}");
            Self { pb: None }
        }
    }

    /// Clear the indicator without printing a final line — the caller prints its
    /// own multi-line summary next, so the two don't collide.
    pub(crate) fn abandon(mut self) {
        if let Some(pb) = self.pb.take() {
            pb.finish_and_clear();
        }
    }
}

impl Drop for Spinner {
    fn drop(&mut self) {
        // Only fires on an early return before finish/abandon (they take `pb`).
        if let Some(pb) = self.pb.take() {
            pb.finish_and_clear();
        }
    }
}

/// A live compress readout driven by `counter` (bytes pulled from the input).
/// `total` is the input's byte length when known (a file) — then a real
/// percentage bar is drawn; `None` (stdin) draws a bytes + rate spinner instead.
#[derive(Debug)]
pub(crate) struct Bar {
    pb: Option<ProgressBar>,
}

impl Bar {
    pub(crate) fn start(
        message: impl Into<String>,
        quiet: bool,
        counter: Arc<AtomicU64>,
        total: Option<u64>,
    ) -> Self {
        if quiet {
            return Self { pb: None };
        }
        let message = message.into();
        if !std::io::stderr().is_terminal() {
            eprintln!("{message}");
            return Self { pb: None };
        }
        let pb = match total {
            Some(len) => {
                let pb = ProgressBar::new(len);
                pb.set_style(
                    ProgressStyle::with_template(
                        "{spinner:.cyan} {msg} [{bar:28.cyan/blue}] {bytes}/{total_bytes} ({percent}%)",
                    )
                    .expect("valid bar template")
                    .tick_strings(TICK_STRINGS)
                    .progress_chars("=>-"),
                );
                pb
            }
            None => {
                let pb = ProgressBar::new_spinner();
                pb.set_style(
                    ProgressStyle::with_template(
                        "{spinner:.cyan} {msg} {bytes} in ({bytes_per_sec})",
                    )
                    .expect("valid stream template")
                    .tick_strings(TICK_STRINGS),
                );
                pb
            }
        };
        pb.set_message(message);
        pb.enable_steady_tick(Duration::from_millis(100));
        // Feed the bar from the shared counter until it is finished. `is_finished`
        // flips on `finish_and_clear`, so the thread exits when the caller
        // abandons or drops the bar — no separate stop flag needed.
        let feeder = pb.clone();
        std::thread::spawn(move || {
            while !feeder.is_finished() {
                feeder.set_position(counter.load(Ordering::Relaxed));
                std::thread::sleep(REFRESH);
            }
        });
        Self { pb: Some(pb) }
    }

    /// Clear the indicator without printing a final line (the summary follows).
    pub(crate) fn abandon(mut self) {
        if let Some(pb) = self.pb.take() {
            pb.finish_and_clear();
        }
    }
}

impl Drop for Bar {
    fn drop(&mut self) {
        if let Some(pb) = self.pb.take() {
            pb.finish_and_clear();
        }
    }
}
