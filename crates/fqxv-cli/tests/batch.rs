//! End-to-end CLI tests for the batch (`info`/`verify` over multiple files or a
//! directory) and `compress --estimate tsv` machine-readable outputs.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The fqxv binary under test (set by Cargo for integration tests).
const FQXV: &str = env!("CARGO_BIN_EXE_fqxv");

/// A private, per-test temp dir Cargo manages for this test binary.
fn tmp(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    p.push(name);
    p
}

/// A tiny four-record FASTQ with `n`-influenced content so archives differ.
fn fastq(seed: u8) -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..4u8 {
        let base = b"ACGT"[((seed + i) % 4) as usize];
        out.extend_from_slice(format!("@r.{seed}.{i} lane\n").as_bytes());
        out.extend_from_slice(&[base; 8]);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(b"IIIIIIII\n");
    }
    out
}

fn run(args: &[&str]) -> std::process::Output {
    let out = Command::new(FQXV).args(args).output().expect("spawn fqxv");
    assert!(
        out.status.success(),
        "fqxv {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

/// Compress one small FASTQ to `arc` (single-thread, deterministic).
fn compress_to(arc: &Path, seed: u8) {
    let src = arc.with_extension("fastq");
    fs::write(&src, fastq(seed)).unwrap();
    run(&[
        "compress",
        src.to_str().unwrap(),
        "-o",
        arc.to_str().unwrap(),
        "--force",
        "--threads",
        "1",
    ]);
}

#[test]
fn info_tsv_single_is_legacy_shape() {
    let dir = tmp("info_single");
    fs::create_dir_all(&dir).unwrap();
    let arc = dir.join("a.fqxv");
    compress_to(&arc, 1);

    let out = run(&["info", arc.to_str().unwrap(), "--tsv"]);
    let stdout = String::from_utf8(out.stdout).unwrap();
    let mut lines = stdout.lines();
    let header = lines.next().unwrap();
    // A single file keeps the stable, no-`file`-column layout the bench harness
    // parses by position (names/seq/qual at indices 7/8/9).
    assert!(
        !header.starts_with("file\t"),
        "single file must not add a file column"
    );
    assert_eq!(header.split('\t').next(), Some("file_size"));
    let cols: Vec<&str> = header.split('\t').collect();
    assert_eq!(cols[7], "names_bytes");
    assert_eq!(cols[8], "seq_bytes");
    assert_eq!(cols[9], "qual_bytes");
    // Exactly one data row.
    assert_eq!(lines.count(), 1);
}

#[test]
fn info_tsv_batch_recurses_and_keys_by_file() {
    let dir = tmp("info_batch");
    fs::create_dir_all(dir.join("sub")).unwrap();
    compress_to(&dir.join("a.fqxv"), 1);
    compress_to(&dir.join("b.fqxv"), 2);
    compress_to(&dir.join("sub/c.fqxv"), 3);

    let out = run(&["info", dir.to_str().unwrap(), "--tsv"]);
    let stdout = String::from_utf8(out.stdout).unwrap();
    let mut lines = stdout.lines();
    let header = lines.next().unwrap();
    assert!(
        header.starts_with("file\t"),
        "batch must lead with a file column"
    );
    // One row per archive, found recursively (including sub/c.fqxv).
    let rows: Vec<&str> = lines.collect();
    assert_eq!(
        rows.len(),
        3,
        "recursive scan should find all three archives"
    );
    assert!(rows.iter().any(|r| r.contains("sub/c.fqxv")));
    // Each row's first field is the path.
    for r in &rows {
        assert!(r.split('\t').next().unwrap().ends_with(".fqxv"));
    }
}

#[test]
fn info_json_batch_is_array_single_is_object() {
    let dir = tmp("info_json");
    fs::create_dir_all(&dir).unwrap();
    compress_to(&dir.join("a.fqxv"), 1);
    compress_to(&dir.join("b.fqxv"), 2);

    let single = run(&["info", dir.join("a.fqxv").to_str().unwrap(), "--json"]);
    assert_eq!(
        String::from_utf8(single.stdout)
            .unwrap()
            .trim_start()
            .chars()
            .next(),
        Some('{')
    );

    let batch = run(&["info", dir.to_str().unwrap(), "--json"]);
    assert_eq!(
        String::from_utf8(batch.stdout)
            .unwrap()
            .trim_start()
            .chars()
            .next(),
        Some('[')
    );
}

#[test]
fn verify_batch_reports_all_and_fails_if_any_corrupt() {
    let dir = tmp("verify_batch");
    fs::create_dir_all(&dir).unwrap();
    compress_to(&dir.join("good.fqxv"), 1);
    compress_to(&dir.join("bad.fqxv"), 2);
    // A file that isn't a container at all still yields a failing entry.
    fs::write(dir.join("garbage.fqxv"), [0u8; 20]).unwrap();

    let out = Command::new(FQXV)
        .args(["verify", dir.to_str().unwrap(), "--tsv"])
        .output()
        .expect("spawn fqxv");
    assert!(
        !out.status.success(),
        "a batch with a bad archive must exit non-zero"
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let header = stdout.lines().next().unwrap();
    assert!(
        header.starts_with("file\t"),
        "batch verify leads with a file column"
    );
    // The good archive is still reported despite the garbage one.
    assert!(stdout.contains("good.fqxv"));
    // The unreadable file becomes a failing `readable` row rather than aborting.
    assert!(stdout.contains("garbage.fqxv\treadable\tfail"));
}

#[test]
fn estimate_tsv_reports_file_and_sizes() {
    let dir = tmp("estimate_tsv");
    fs::create_dir_all(&dir).unwrap();
    let src = dir.join("reads.fastq");
    // Enough records that the archive is non-empty and the ratio is finite.
    let mut fq = Vec::new();
    for i in 0..500 {
        fq.extend_from_slice(format!("@r.{i} lane\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").as_bytes());
    }
    fs::write(&src, fq).unwrap();

    let out = run(&["compress", src.to_str().unwrap(), "--estimate", "tsv"]);
    let stdout = String::from_utf8(out.stdout).unwrap();
    let mut lines = stdout.lines();
    assert_eq!(
        lines.next(),
        Some("file\tinput_bytes\test_fqxv_bytes\tratio")
    );
    let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(row.len(), 4);
    assert!(
        row[0].ends_with("reads.fastq"),
        "first column is the input file"
    );
    assert!(row[1].parse::<u64>().unwrap() > 0);
    assert!(row[2].parse::<u64>().unwrap() > 0);
    assert!(row[3].parse::<f64>().unwrap() > 0.0);
}

#[test]
fn estimate_on_empty_input_succeeds() {
    // An empty input has no reads, but `compress` accepts it, so `--estimate` must
    // report it cleanly (exit 0) rather than erroring. Covers both output forms.
    let dir = tmp("estimate_empty");
    fs::create_dir_all(&dir).unwrap();
    let src = dir.join("empty.fastq");
    fs::write(&src, b"").unwrap();

    let human = run(&["compress", src.to_str().unwrap(), "--estimate"]);
    assert!(
        String::from_utf8_lossy(&human.stdout).contains("no reads"),
        "human estimate should note there are no reads"
    );

    let tsv = run(&["compress", src.to_str().unwrap(), "--estimate", "tsv"]);
    let stdout = String::from_utf8(tsv.stdout).unwrap();
    let mut lines = stdout.lines();
    assert_eq!(
        lines.next(),
        Some("file\tinput_bytes\test_fqxv_bytes\tratio")
    );
    let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(row.len(), 4);
    assert_eq!(row[2], "0", "estimated archive bytes are 0 for empty input");
}

#[test]
fn compress_verify_confirms_roundtrip() {
    let dir = tmp("compress_verify");
    fs::create_dir_all(&dir).unwrap();
    let src = dir.join("reads.fastq");
    fs::write(&src, fastq(3)).unwrap();
    let arc = dir.join("reads.fqxv");

    let out = run(&[
        "compress",
        src.to_str().unwrap(),
        "-o",
        arc.to_str().unwrap(),
        "--force",
        "--verify",
        "--threads",
        "1",
    ]);
    assert!(arc.exists(), "verified compress still writes the archive");
    // The verification note is printed to stderr alongside the summary.
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(
        stderr.contains("verified"),
        "expected a verification note, got: {stderr}"
    );
}

#[test]
fn compress_verify_conflicts_with_estimate() {
    let dir = tmp("verify_estimate_conflict");
    fs::create_dir_all(&dir).unwrap();
    let src = dir.join("reads.fastq");
    fs::write(&src, fastq(4)).unwrap();

    let out = Command::new(FQXV)
        .args(["compress", src.to_str().unwrap(), "--verify", "--estimate"])
        .output()
        .expect("spawn fqxv");
    assert!(
        !out.status.success(),
        "--verify with --estimate must be rejected"
    );
}
