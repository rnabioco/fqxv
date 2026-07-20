//! Bakes the git provenance of a *development* build into `--version`.
//!
//! A released binary should say exactly what it is — `fqxv 0.3.0` — while a
//! build from a working tree should say which commit it came from and whether
//! that tree had uncommitted changes, so a bug report from one is actionable.
//!
//! "Release" means git positively confirms HEAD is a clean checkout of the tag
//! matching this crate's version (`v0.3.0` for 0.3.0). Anything else — commits
//! past the tag, a dirty tree, an untagged branch — appends the description:
//!
//! ```text
//! fqxv 0.3.0                          # the v0.3.0 tag, clean
//! fqxv 0.3.0 (v0.3.0-7-gab12cd34)     # 7 commits past it
//! fqxv 0.3.0 (v0.3.0-7-gab12cd34-dirty)
//! fqxv 0.3.0 (ab12cd34-dirty)         # no tag reachable (shallow clone)
//! ```
//!
//! The bias is deliberate: an unconfirmed state is reported as a dev build, not
//! silently as a release. A released binary mislabeled with its own commit is
//! noise; a dev build claiming to be 0.3.0 is a support dead end.
//!
//! With no git at all — a crates.io tarball, an exported source tree — every
//! command below fails and the version stays the bare `CARGO_PKG_VERSION`,
//! which is the honest answer when there is no provenance to report.

use std::path::Path;
use std::process::Command;

fn main() {
    let pkg = std::env::var("CARGO_PKG_VERSION").unwrap_or_default();

    // Rebuild when the checkout moves, so the string does not go stale across
    // commits and branch switches.
    //
    // The `-dirty` flag tracks changes that could affect the binary, not every
    // change in the worktree. Editing a crate in the compile graph reruns this
    // script (measured), so a build whose code differs from the tag always says
    // so. Dirtying only files nothing compiles — docs, bench scripts, CI — does
    // not, and the version will still read clean. That is the right answer
    // about the *binary*, which is what a version string describes.
    if let Some(git_dir) = git_dir() {
        for f in ["HEAD", "index"] {
            let p = git_dir.join(f);
            if p.exists() {
                println!("cargo:rerun-if-changed={}", p.display());
            }
        }
    }
    println!("cargo:rerun-if-changed=build.rs");

    println!("cargo:rustc-env=FQXV_VERSION={}", version(&pkg));
}

/// Resolve the git directory, asking git rather than assuming `.git/` — in a
/// linked worktree `.git` is a *file* pointing elsewhere, and this repo is
/// worked on through worktrees.
fn git_dir() -> Option<std::path::PathBuf> {
    let out = git(&["rev-parse", "--absolute-git-dir"])?;
    let p = std::path::PathBuf::from(out);
    p.is_dir().then_some(p)
}

/// `CARGO_PKG_VERSION`, plus a git description when this is not a confirmed
/// clean release build.
fn version(pkg: &str) -> String {
    // One command covers all of it: nearest tag, commits since, abbreviated
    // hash, and a `-dirty` suffix when the worktree has changes. `--always`
    // falls back to a bare hash when no tag is reachable.
    let Some(desc) = git(&["describe", "--tags", "--always", "--dirty=-dirty"]) else {
        return pkg.to_string();
    };
    if desc == format!("v{pkg}") {
        // Exactly the release tag, and clean — `--dirty` would have suffixed it.
        pkg.to_string()
    } else {
        format!("{pkg} ({desc})")
    }
}

/// Run a git command in the manifest directory, returning trimmed stdout on a
/// clean exit. Any failure — git missing, not a repository, no commits yet — is
/// `None`, and the caller degrades to the plain version.
fn git(args: &[&str]) -> Option<String> {
    let dir = std::env::var("CARGO_MANIFEST_DIR").ok()?;
    let out = Command::new("git")
        .arg("-C")
        .arg(Path::new(&dir))
        .args(args)
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8(out.stdout).ok()?.trim().to_string();
    (!s.is_empty()).then_some(s)
}
