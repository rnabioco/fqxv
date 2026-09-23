//! Passphrase resolution for `--encrypt` / `--password-file` / `$FQXV_PASSWORD`.
//!
//! Precedence is the same for every command: `--password-file` (the file's raw
//! bytes) beats `$FQXV_PASSWORD` (its UTF-8 string bytes) beats an interactive
//! hidden-terminal prompt. A `str` source (env var, prompt) is always UTF-8; a
//! `--password-file` is used byte-for-byte with no encoding assumed — mixing
//! sources for the same archive only works if the bytes happen to agree, so
//! prefer one source consistently for a given archive.

use anyhow::Context;
use std::path::Path;

/// Env var carrying a passphrase for scripted/CI use — visible to other
/// processes of the same user (e.g. via `/proc/<pid>/environ`), but avoids
/// shell history and works without a TTY.
const FQXV_PASSWORD_ENV: &str = "FQXV_PASSWORD";

/// Resolve a passphrase for `compress --encrypt`.
///
/// Interactive fallback asks twice and requires the two entries to match: a
/// typo here is unrecoverable (there is no way to decrypt without the exact
/// passphrase), unlike a decode-time guess, which just fails cleanly and can
/// be retried. Also rejects an empty passphrase, which would provide
/// essentially no confidentiality.
pub(crate) fn resolve_passphrase_for_compress(
    password_file: Option<&Path>,
) -> anyhow::Result<Vec<u8>> {
    if let Some(path) = password_file {
        return read_password_file(path);
    }
    if let Ok(pw) = std::env::var(FQXV_PASSWORD_ENV) {
        return Ok(pw.into_bytes());
    }
    let first = rpassword::prompt_password("passphrase: ").context("reading passphrase")?;
    if first.is_empty() {
        anyhow::bail!("passphrase must not be empty");
    }
    let second =
        rpassword::prompt_password("confirm passphrase: ").context("reading passphrase")?;
    if first != second {
        anyhow::bail!("passphrases did not match");
    }
    Ok(first.into_bytes())
}

/// Resolve a passphrase for `decompress`/`verify`.
///
/// `prompt_if_needed` gates the interactive fallback: the caller passes `true`
/// only when it has already confirmed (via a cheap header peek) that the
/// archive is encrypted and a real terminal is available to prompt on —
/// `verify` never sets it, so a routine `fqxv verify *.fqxv` batch run never
/// blocks on a surprise TTY prompt; stdin input never sets it either, since a
/// stream can't pause mid-read for a prompt.
///
/// Returns `None` when nothing is configured and prompting isn't appropriate —
/// the library surfaces `fqxv::Error::PasswordRequired` if the archive turns
/// out to need one.
pub(crate) fn resolve_passphrase_for_decode(
    password_file: Option<&Path>,
    prompt_if_needed: bool,
) -> anyhow::Result<Option<Vec<u8>>> {
    if let Some(path) = password_file {
        return Ok(Some(read_password_file(path)?));
    }
    if let Ok(pw) = std::env::var(FQXV_PASSWORD_ENV) {
        return Ok(Some(pw.into_bytes()));
    }
    if prompt_if_needed {
        let pw = rpassword::prompt_password("passphrase: ").context("reading passphrase")?;
        return Ok(Some(pw.into_bytes()));
    }
    Ok(None)
}

/// Read a `--password-file`'s raw bytes as the passphrase. A single trailing
/// `\n` (or `\r\n`) is stripped, since a shell-authored file commonly ends in
/// one; every other byte is used verbatim — no encoding is assumed.
fn read_password_file(path: &Path) -> anyhow::Result<Vec<u8>> {
    let mut bytes = std::fs::read(path)
        .with_context(|| format!("reading passphrase file {}", path.display()))?;
    if bytes.last() == Some(&b'\n') {
        bytes.pop();
        if bytes.last() == Some(&b'\r') {
            bytes.pop();
        }
    }
    Ok(bytes)
}
