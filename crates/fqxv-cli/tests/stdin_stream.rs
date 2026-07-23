//! `fqxv decompress -` streams an archive from stdin and must produce
//! byte-identical interleaved FASTQ to a local-file decompress. This is the path
//! behind remote reads: `aws s3 cp s3://… - | fqxv decompress - -Z | <aligner>`,
//! where the transfer tool owns the download and fqxv just decodes the stream.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

const FQXV: &str = env!("CARGO_BIN_EXE_fqxv");

fn tmp(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    p.push(name);
    p
}

fn write_fastq(path: &PathBuf) {
    let mut state: u64 = 0x1234_5678_9abc_def0;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut out = String::new();
    for i in 0..4000 {
        let seq: String = (0..120)
            .map(|_| b"ACGT"[(next() % 4) as usize] as char)
            .collect();
        out.push_str(&format!("@read{i} x={i}\n{seq}\n+\n{}\n", "I".repeat(120)));
    }
    std::fs::write(path, out).unwrap();
}

fn build_archive(tag: &str) -> PathBuf {
    let fastq = tmp(&format!("{tag}.fastq"));
    let archive = tmp(&format!("{tag}.fqxv"));
    write_fastq(&fastq);
    let ok = Command::new(FQXV)
        .args(["--quiet", "compress", "--force"])
        .arg(&fastq)
        .arg("-o")
        .arg(&archive)
        .status()
        .expect("spawn fqxv compress");
    assert!(ok.success(), "compress failed");
    archive
}

#[test]
fn stdin_decompress_matches_file() {
    let archive = build_archive("stdin_stream");

    // Reference: decompress the file to stdout.
    let local = Command::new(FQXV)
        .args(["--quiet", "decompress"])
        .arg(&archive)
        .arg("-Z")
        .output()
        .expect("spawn fqxv");
    assert!(local.status.success());

    // Stream the same bytes into `decompress -` over a pipe (what `aws s3 cp -`
    // does upstream).
    let bytes = std::fs::read(&archive).unwrap();
    let mut child = Command::new(FQXV)
        .args(["--quiet", "decompress", "-", "-Z"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn fqxv");
    child.stdin.take().unwrap().write_all(&bytes).unwrap();
    let streamed = child.wait_with_output().expect("wait fqxv");
    assert!(
        streamed.status.success(),
        "stdin decompress failed: {}",
        String::from_utf8_lossy(&streamed.stderr)
    );

    assert!(!local.stdout.is_empty());
    assert!(
        local.stdout == streamed.stdout,
        "stdin FASTQ differs from file decode"
    );
}

#[test]
fn stdin_split_is_rejected() {
    // --split needs a second pass for the header; refuse cleanly on a stream.
    let out = Command::new(FQXV)
        .args(["--quiet", "decompress", "-", "--split", "/tmp/x"])
        .stdin(Stdio::null())
        .output()
        .expect("spawn fqxv");
    assert!(!out.status.success(), "--split from stdin must fail");
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("split"),
        "expected a --split error, got: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}
