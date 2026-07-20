//! Atomic file output, for both directions of the codec.
//!
//! Every named output file — the `.fqxv` an archive compresses to, and the
//! FASTQ a decompress writes — streams into a sibling temp file that is moved
//! to the destination only once the whole thing is on disk. An interrupted or
//! failed run therefore never leaves a partial file at the destination path.
//!
//! Both directions need this, for different reasons. A truncated, footer-less
//! `.fqxv` is at least *detectable* — [`fqxv::verify`] rightly calls it
//! CORRUPT. Truncated FASTQ is not. BGZF blocks are self-contained, so a
//! partial `.fastq.gz` inflates cleanly — `zcat` exits 0 with no warning — and
//! FASTQ carries no terminator or record count, so nothing downstream can tell
//! a half-written mate file from a complete one. Interrupting a `--split`
//! decode of a 600k-read pair six seconds in produced two files that gzip
//! called valid, held 524,100 reads apiece, and ended mid-record on a short
//! quality line. Silently losing 13% of a run that way is the failure this
//! guards against: the partial bytes must never reach the destination name.
//!
//! Cleanup is covered on both exits either direction can take: a `?` error
//! bail runs [`AtomicOutput`]'s [`Drop`], and a Ctrl-C/terminate runs the
//! signal handler installed by [`install_interrupt_handler`]. The handler
//! matters on its own — it calls [`std::process::exit`], which does not
//! unwind, so `Drop` alone would never run on a signal.

use std::fs::{self, File};
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use anyhow::Context;

/// Temp files currently being written, so the signal handler can delete them.
/// A list (not a single slot) because one process may write several archives,
/// and an interrupt must clear every in-progress one.
static ACTIVE: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());

/// Install a Ctrl-C / SIGTERM / SIGHUP handler that removes any in-progress
/// temp outputs and exits. Call once from `main`, before any compress or
/// decompress. A leftover temp would be harmless (it never sits at a
/// destination path), but this keeps the working tree clean and lets the exit
/// code carry the interrupt.
pub(crate) fn install_interrupt_handler() {
    // A failure here only means we lose the tidy-up-on-signal nicety; the
    // atomic rename still guarantees no corrupt archive is published, so it is
    // not worth aborting startup over.
    let _ = ctrlc::set_handler(|| {
        remove_all_active();
        // 128 + SIGINT(2): the shell convention for "terminated by Ctrl-C".
        std::process::exit(130);
    });
}

fn remove_all_active() {
    if let Ok(mut active) = ACTIVE.lock() {
        for path in active.drain(..) {
            let _ = fs::remove_file(&path);
        }
    }
}

fn register(path: &Path) {
    if let Ok(mut active) = ACTIVE.lock() {
        active.push(path.to_path_buf());
    }
}

fn unregister(path: &Path) {
    if let Ok(mut active) = ACTIVE.lock() {
        active.retain(|p| p != path);
    }
}

/// A pending output file: bytes stream into a temp file that becomes
/// `final_path` only on [`commit`](AtomicOutput::commit). Dropping without
/// committing removes the temp.
#[derive(Debug)]
pub(crate) struct AtomicOutput {
    tmp: PathBuf,
    final_path: PathBuf,
    /// Whether the caller passed `-f/--force`, i.e. may replace an existing
    /// destination. Recorded at create so `commit` can pick a publish that
    /// enforces the same rule the initial check did.
    force: bool,
    /// Set by `commit`/`keep_for_inspection` so `Drop` leaves the temp alone.
    settled: bool,
}

impl AtomicOutput {
    /// Create the temp file alongside `final_path` and return it with an open
    /// [`File`] to write the output into. Honors the same overwrite rule as a
    /// direct write against the *final* path: without `force`, an existing
    /// destination is refused.
    ///
    /// This check is the fail-fast one — it rejects a doomed run before minutes
    /// of coding work rather than after. It is not the only one: `commit`
    /// re-enforces the rule atomically, so a destination that appears while we
    /// are working is still not clobbered.
    pub(crate) fn create(final_path: &Path, force: bool) -> anyhow::Result<(Self, File)> {
        if !force && final_path.exists() {
            anyhow::bail!(
                "output {} already exists; pass -f/--force to overwrite",
                final_path.display()
            );
        }
        let tmp = temp_path(final_path);
        // Truncate-create: a temp left by a crashed run was never a published
        // output, so overwriting it is safe.
        let file = File::create(&tmp)
            .with_context(|| format!("creating temporary output {}", tmp.display()))?;
        register(&tmp);
        Ok((
            Self {
                tmp,
                final_path: final_path.to_path_buf(),
                force,
                settled: false,
            },
            file,
        ))
    }

    /// The temp path — for a read-after-write verify before committing.
    pub(crate) fn temp_path(&self) -> &Path {
        &self.tmp
    }

    /// Atomically move the finished temp into place. After this the output
    /// exists at `final_path` and the temp is gone.
    ///
    /// Without `force` this publishes via `hard_link`, which fails rather than
    /// replaces when the destination exists — the atomic create-if-absent that
    /// a direct `File::options().create_new(true)` open used to provide, kept
    /// intact now that the bytes route through a temp. `rename` alone would
    /// silently clobber a file that appeared while we were decoding. With
    /// `force` the caller has asked to replace, so a plain `rename` is right
    /// (and is the only one of the two that is a single atomic step).
    pub(crate) fn commit(mut self) -> anyhow::Result<()> {
        if self.force {
            fs::rename(&self.tmp, &self.final_path).with_context(|| {
                format!(
                    "finalizing {} (from {})",
                    self.final_path.display(),
                    self.tmp.display()
                )
            })?;
        } else {
            match fs::hard_link(&self.tmp, &self.final_path) {
                Ok(()) => {
                    // The destination now has the content; drop our temp name for it.
                    let _ = fs::remove_file(&self.tmp);
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {
                    anyhow::bail!(
                        "output {} already exists; pass -f/--force to overwrite",
                        self.final_path.display()
                    );
                }
                // Not every filesystem supports hard links. Where it is refused
                // outright we fall back to the rename, accepting its weaker
                // create-if-absent guarantee rather than failing the whole run —
                // the `create` check above already covered the common case.
                Err(e)
                    if matches!(
                        e.kind(),
                        io::ErrorKind::Unsupported | io::ErrorKind::PermissionDenied
                    ) =>
                {
                    fs::rename(&self.tmp, &self.final_path).with_context(|| {
                        format!(
                            "finalizing {} (from {})",
                            self.final_path.display(),
                            self.tmp.display()
                        )
                    })?;
                }
                Err(e) => {
                    return Err(e).with_context(|| {
                        format!(
                            "finalizing {} (from {})",
                            self.final_path.display(),
                            self.tmp.display()
                        )
                    });
                }
            }
        }
        unregister(&self.tmp);
        self.settled = true;
        Ok(())
    }

    /// Disarm cleanup and hand back the temp path, leaving the bytes on disk.
    /// Used when a `--verify` read-after-write check fails: we refuse to publish
    /// a bad archive to `final_path`, but keep the output so the failure can be
    /// investigated.
    pub(crate) fn keep_for_inspection(mut self) -> PathBuf {
        unregister(&self.tmp);
        self.settled = true;
        self.tmp.clone()
    }
}

impl Drop for AtomicOutput {
    fn drop(&mut self) {
        if !self.settled {
            let _ = fs::remove_file(&self.tmp);
            unregister(&self.tmp);
        }
    }
}

/// Sibling temp path for `final_path`: same directory (so the `rename` is atomic
/// — a cross-filesystem move would not be), dot-hidden, and pid-tagged so two
/// concurrent writers targeting the same archive don't stomp each other's temp.
fn temp_path(final_path: &Path) -> PathBuf {
    let name = final_path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "out.fqxv".to_string());
    final_path.with_file_name(format!(".{name}.{}.partial", std::process::id()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn scratch_dir(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("fqxv-atomic-{}-{}", tag, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn commit_publishes_and_removes_temp() {
        let dir = scratch_dir("commit");
        let out = dir.join("archive.fqxv");
        let (ao, mut f) = AtomicOutput::create(&out, false).unwrap();
        let tmp = ao.temp_path().to_path_buf();
        assert!(tmp.exists() && !out.exists());
        f.write_all(b"payload").unwrap();
        drop(f);
        ao.commit().unwrap();
        assert!(out.exists() && !tmp.exists());
        assert_eq!(std::fs::read(&out).unwrap(), b"payload");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn drop_without_commit_removes_temp_and_never_creates_final() {
        let dir = scratch_dir("drop");
        let out = dir.join("archive.fqxv");
        let tmp = {
            let (ao, mut f) = AtomicOutput::create(&out, false).unwrap();
            f.write_all(b"partial").unwrap();
            let tmp = ao.temp_path().to_path_buf();
            // `ao` dropped here without commit — models a `?` bail mid-compress.
            drop(ao);
            tmp
        };
        assert!(!tmp.exists(), "temp must be cleaned up on drop");
        assert!(!out.exists(), "final archive must never appear on failure");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn refuses_existing_final_unless_forced() {
        let dir = scratch_dir("exists");
        let out = dir.join("archive.fqxv");
        std::fs::write(&out, b"existing").unwrap();

        let err = AtomicOutput::create(&out, false).unwrap_err().to_string();
        assert!(err.contains("already exists"), "unexpected error: {err}");
        assert_eq!(std::fs::read(&out).unwrap(), b"existing");

        // With force the existing archive is left in place until commit swaps it.
        let (ao, mut f) = AtomicOutput::create(&out, true).unwrap();
        f.write_all(b"replacement").unwrap();
        drop(f);
        assert_eq!(std::fs::read(&out).unwrap(), b"existing");
        ao.commit().unwrap();
        assert_eq!(std::fs::read(&out).unwrap(), b"replacement");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn keep_for_inspection_preserves_temp_off_the_final_path() {
        let dir = scratch_dir("keep");
        let out = dir.join("archive.fqxv");
        let (ao, mut f) = AtomicOutput::create(&out, false).unwrap();
        f.write_all(b"suspect").unwrap();
        drop(f);
        let kept = ao.keep_for_inspection();
        assert!(kept.exists(), "kept temp must survive");
        assert!(
            !out.exists(),
            "a failed verify must not publish to the final path"
        );
        assert_eq!(std::fs::read(&kept).unwrap(), b"suspect");
        std::fs::remove_dir_all(&dir).ok();
    }
}
