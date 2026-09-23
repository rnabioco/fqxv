//! End-to-end CLI tests for `--encrypt` / `--password-file` / `$FQXV_PASSWORD`.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

/// The fqxv binary under test (set by Cargo for integration tests).
const FQXV: &str = env!("CARGO_BIN_EXE_fqxv");

/// A private temp dir Cargo manages for this test binary.
fn tmp(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    p.push(name);
    p
}

const SAMPLE: &[u8] = b"\
@read.1 lane\n\
ACGTACGTACGT\n\
+\n\
IIIIFFFF####\n\
@read.2 lane\n\
NNGGCCTAGGCC\n\
+\n\
0:F,+.HHHIII\n";

fn run(args: &[&str]) -> std::process::Output {
    let out = Command::new(FQXV).args(args).output().expect("spawn fqxv");
    assert!(
        out.status.success(),
        "fqxv {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

/// A TSV `info` column's value by header name (single-archive, no `file` column).
fn tsv_col(stdout: &str, name: &str) -> String {
    let mut lines = stdout.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let data: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = header
        .iter()
        .position(|&h| h == name)
        .unwrap_or_else(|| panic!("no {name} column in {header:?}"));
    data[col].to_string()
}

#[test]
fn encrypt_password_file_roundtrip() {
    let in_path = tmp("enc_pf_in.fastq");
    let arc_path = tmp("enc_pf.fqxv");
    let pw_path = tmp("enc_pf.pw");
    let rt_path = tmp("enc_pf_rt.fastq");
    fs::write(&in_path, SAMPLE).unwrap();
    fs::write(&pw_path, b"correct horse battery staple\n").unwrap();

    run(&[
        "compress",
        in_path.to_str().unwrap(),
        "-o",
        arc_path.to_str().unwrap(),
        "--force",
        "--encrypt",
        "--password-file",
        pw_path.to_str().unwrap(),
        "--threads",
        "1",
    ]);

    // info --tsv reports encrypted without a password.
    let info = run(&["info", arc_path.to_str().unwrap(), "--tsv"]);
    let stdout = String::from_utf8(info.stdout).unwrap();
    assert_eq!(tsv_col(&stdout, "encrypted"), "1");
    assert!(tsv_col(&stdout, "encrypted_bytes").parse::<u64>().unwrap() > 0);

    run(&[
        "decompress",
        arc_path.to_str().unwrap(),
        "-o",
        rt_path.to_str().unwrap(),
        "--force",
        "--password-file",
        pw_path.to_str().unwrap(),
        "--threads",
        "1",
    ]);
    let rt = fs::read(&rt_path).unwrap();
    assert_eq!(
        rt, SAMPLE,
        "decrypted output must match the original FASTQ exactly"
    );
}

#[test]
fn encrypt_fqxv_password_env_roundtrip() {
    let in_path = tmp("enc_env_in.fastq");
    let arc_path = tmp("enc_env.fqxv");
    let rt_path = tmp("enc_env_rt.fastq");
    fs::write(&in_path, SAMPLE).unwrap();

    let compress_out = Command::new(FQXV)
        .args([
            "compress",
            in_path.to_str().unwrap(),
            "-o",
            arc_path.to_str().unwrap(),
            "--force",
            "--encrypt",
            "--threads",
            "1",
        ])
        .env("FQXV_PASSWORD", "env-var-passphrase")
        .output()
        .expect("spawn fqxv compress");
    assert!(
        compress_out.status.success(),
        "compress failed: {}",
        String::from_utf8_lossy(&compress_out.stderr)
    );

    let decompress_out = Command::new(FQXV)
        .args([
            "decompress",
            arc_path.to_str().unwrap(),
            "-o",
            rt_path.to_str().unwrap(),
            "--force",
            "--threads",
            "1",
        ])
        .env("FQXV_PASSWORD", "env-var-passphrase")
        .output()
        .expect("spawn fqxv decompress");
    assert!(
        decompress_out.status.success(),
        "decompress failed: {}",
        String::from_utf8_lossy(&decompress_out.stderr)
    );
    assert_eq!(fs::read(&rt_path).unwrap(), SAMPLE);
}

#[test]
fn wrong_password_fails_cleanly_not_a_panic() {
    let in_path = tmp("enc_wrong_in.fastq");
    let arc_path = tmp("enc_wrong.fqxv");
    let pw_path = tmp("enc_wrong.pw");
    let wrong_pw_path = tmp("enc_wrong.wrong.pw");
    fs::write(&in_path, SAMPLE).unwrap();
    fs::write(&pw_path, b"the real passphrase").unwrap();
    fs::write(&wrong_pw_path, b"a wrong guess").unwrap();

    run(&[
        "compress",
        in_path.to_str().unwrap(),
        "-o",
        arc_path.to_str().unwrap(),
        "--force",
        "--encrypt",
        "--password-file",
        pw_path.to_str().unwrap(),
        "--threads",
        "1",
    ]);

    let out = Command::new(FQXV)
        .args([
            "decompress",
            arc_path.to_str().unwrap(),
            "-o",
            tmp("enc_wrong_rt.fastq").to_str().unwrap(),
            "--force",
            "--password-file",
            wrong_pw_path.to_str().unwrap(),
            "--threads",
            "1",
        ])
        .output()
        .expect("spawn fqxv decompress");
    assert!(
        !out.status.success(),
        "decompress with the wrong password must fail"
    );
    // A clean process exit code, not a signal (crash/panic/abort).
    assert!(
        out.status.code().is_some(),
        "must exit cleanly, not crash: {out:?}"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("panicked"),
        "must not panic; stderr: {stderr}"
    );
}

#[test]
fn missing_password_reports_password_required() {
    let in_path = tmp("enc_missing_in.fastq");
    let arc_path = tmp("enc_missing.fqxv");
    let pw_path = tmp("enc_missing.pw");
    fs::write(&in_path, SAMPLE).unwrap();
    fs::write(&pw_path, b"needs a password").unwrap();

    run(&[
        "compress",
        in_path.to_str().unwrap(),
        "-o",
        arc_path.to_str().unwrap(),
        "--force",
        "--encrypt",
        "--password-file",
        pw_path.to_str().unwrap(),
        "--threads",
        "1",
    ]);

    // No --password-file, no FQXV_PASSWORD, and stdin/stderr are not a TTY under
    // the test harness, so this must not hang waiting on a prompt.
    let out = Command::new(FQXV)
        .args([
            "decompress",
            arc_path.to_str().unwrap(),
            "-o",
            tmp("enc_missing_rt.fastq").to_str().unwrap(),
            "--force",
            "--threads",
            "1",
        ])
        .env_remove("FQXV_PASSWORD")
        .output()
        .expect("spawn fqxv decompress");
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("passphrase") || stderr.contains("password"),
        "stderr should mention the missing passphrase, got: {stderr}"
    );
}

#[test]
fn verify_without_password_runs_crc_checks_with_password_runs_footer_tag() {
    let in_path = tmp("enc_verify_in.fastq");
    let arc_path = tmp("enc_verify.fqxv");
    let pw_path = tmp("enc_verify.pw");
    fs::write(&in_path, SAMPLE).unwrap();
    fs::write(&pw_path, b"verify passphrase").unwrap();

    run(&[
        "compress",
        in_path.to_str().unwrap(),
        "-o",
        arc_path.to_str().unwrap(),
        "--force",
        "--encrypt",
        "--password-file",
        pw_path.to_str().unwrap(),
        "--threads",
        "1",
    ]);

    // No password: the archive still verifies (keyless CRC checks only), and
    // the footer-tag row is explicitly reported as skipped, not silently
    // omitted.
    let no_pw = run(&["verify", arc_path.to_str().unwrap(), "--tsv"]);
    let no_pw_out = String::from_utf8(no_pw.stdout).unwrap();
    let no_pw_row = no_pw_out
        .lines()
        .find(|l| l.starts_with("footer auth tag"))
        .unwrap_or_else(|| panic!("no footer auth tag row in: {no_pw_out}"));
    assert!(no_pw_row.contains("skipped"), "row was: {no_pw_row}");

    // With the correct password: the footer auth tag check also runs and passes.
    let with_pw = run(&[
        "verify",
        arc_path.to_str().unwrap(),
        "--tsv",
        "--password-file",
        pw_path.to_str().unwrap(),
    ]);
    let with_pw_out = String::from_utf8(with_pw.stdout).unwrap();
    let with_pw_row = with_pw_out
        .lines()
        .find(|l| l.starts_with("footer auth tag"))
        .unwrap_or_else(|| panic!("no footer auth tag row in: {with_pw_out}"));
    assert!(with_pw_row.contains("\tok\t"), "row was: {with_pw_row}");
}

#[test]
fn encrypt_and_reorder_are_rejected_before_any_output_is_written() {
    let in_path = tmp("enc_reorder_in.fastq");
    let arc_path = tmp("enc_reorder.fqxv");
    fs::write(&in_path, SAMPLE).unwrap();
    // A stray leftover from a previous run would defeat the "no output" check.
    let _ = fs::remove_file(&arc_path);

    let out = Command::new(FQXV)
        .args([
            "compress",
            in_path.to_str().unwrap(),
            "-o",
            arc_path.to_str().unwrap(),
            "--force",
            "--encrypt",
            "--order",
            "any",
        ])
        .env("FQXV_PASSWORD", "unused")
        .output()
        .expect("spawn fqxv compress");
    assert!(
        !out.status.success(),
        "--encrypt --order any must be rejected"
    );
    assert!(
        !arc_path.exists(),
        "no output archive must be written when the combination is rejected"
    );
}
