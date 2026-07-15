//! Terminal spinner for the slow phases of compress/decompress.
//!
//! On a TTY, [`Spinner::start`] renders a ticking spinner next to a message and
//! swaps it for a final line on [`Spinner::finish`]. On a non-TTY (CI logs,
//! redirects) the start message is printed as a plain line and the finish line
//! follows, so logs still record what happened without any control characters.
//! When `--quiet` is set the spinner is fully inert — nothing is drawn or
//! printed — because the user-facing summary is suppressed too.
//!
//! The spinner always renders on **stderr**, so it coexists with FASTQ streamed
//! to stdout (`decompress -Z`). If the handle is dropped without an explicit
//! [`Spinner::finish`]/[`Spinner::abandon`] — e.g. a `?` bails mid-run — the
//! [`Drop`] impl clears the live line so the error message isn't tangled with a
//! half-drawn spinner.

use std::io::IsTerminal;
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

/// Braille tick frames, matching the rnabioco tooling look.
const TICK_STRINGS: &[&str] = &["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"];

/// A single-line status indicator that adapts to TTY, non-TTY, and quiet modes.
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

    /// Clear the spinner without printing a final line — the caller prints its own
    /// multi-line summary next, so the two don't collide.
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
