//! Atomic archive output.
//!
//! Compression writes to a sibling temp file and `rename`s it into place only
//! once the whole archive — header, blocks, *and* footer trailer — is on disk.
//! An interrupted or failed run therefore never leaves a truncated,
//! footer-less `.fqxv` at the destination path (which [`fqxv::verify`] would
//! rightly call CORRUPT). Cleanup is covered on both exits a compress can take:
//! a `?` error bail runs [`AtomicOutput`]'s [`Drop`], and a Ctrl-C/terminate
//! runs the signal handler installed by [`install_interrupt_handler`].

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use anyhow::Context;

/// Temp files currently being written, so the signal handler can delete them.
/// A list (not a single slot) because one process may write several archives,
/// and an interrupt must clear every in-progress one.
static ACTIVE: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());

/// Install a Ctrl-C / SIGTERM / SIGHUP handler that removes any in-progress
/// temp outputs and exits. Call once from `main`, before any compression. A
/// leftover temp would be harmless (it is never a valid archive and never sits
/// at the destination path), but this keeps the working tree clean and lets the
/// exit code carry the interrupt.
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

/// A pending archive: bytes stream into a temp file that becomes `final_path`
/// only on [`commit`](AtomicOutput::commit). Dropping without committing removes
/// the temp.
#[derive(Debug)]
pub(crate) struct AtomicOutput {
    tmp: PathBuf,
    final_path: PathBuf,
    /// Set by `commit`/`keep_for_inspection` so `Drop` leaves the temp alone.
    settled: bool,
}

impl AtomicOutput {
    /// Create the temp file alongside `final_path` and return it with an open
    /// [`File`] to write the archive into. Honors the same overwrite rule as a
    /// direct write against the *final* path: without `force`, an existing
    /// archive is refused.
    pub(crate) fn create(final_path: &Path, force: bool) -> anyhow::Result<(Self, File)> {
        if !force && final_path.exists() {
            anyhow::bail!(
                "output {} already exists; pass -f/--force to overwrite",
                final_path.display()
            );
        }
        let tmp = temp_path(final_path);
        // Truncate-create: a temp left by a crashed run was never a real archive,
        // so overwriting it is safe.
        let file = File::create(&tmp)
            .with_context(|| format!("creating temporary output {}", tmp.display()))?;
        register(&tmp);
        Ok((
            Self {
                tmp,
                final_path: final_path.to_path_buf(),
                settled: false,
            },
            file,
        ))
    }

    /// The temp path — for a read-after-write verify before committing.
    pub(crate) fn temp_path(&self) -> &Path {
        &self.tmp
    }

    /// Atomically move the finished temp into place. After this the archive
    /// exists at `final_path` and the temp is gone.
    pub(crate) fn commit(mut self) -> anyhow::Result<()> {
        fs::rename(&self.tmp, &self.final_path).with_context(|| {
            format!(
                "finalizing {} (from {})",
                self.final_path.display(),
                self.tmp.display()
            )
        })?;
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
