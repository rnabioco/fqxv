//! End-to-end CLI test for `fqxv extract`: a read-index range must yield exactly
//! the same reads as the corresponding slice of a full `decompress`, open-ended
//! and clamped ranges must behave, and reordered archives must be rejected.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

/// The fqxv binary under test (set by Cargo for integration tests).
const FQXV: &str = env!("CARGO_BIN_EXE_fqxv");

fn tmp(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    p.push(name);
    p
}

fn run(args: &[&str]) -> Vec<u8> {
    let out = Command::new(FQXV).args(args).output().expect("spawn fqxv");
    assert!(
        out.status.success(),
        "fqxv {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn run_fail(args: &[&str]) {
    let out = Command::new(FQXV).args(args).output().expect("spawn fqxv");
    assert!(
        !out.status.success(),
        "fqxv {args:?} unexpectedly succeeded"
    );
}

/// Split FASTQ into per-record byte blocks (4 lines each, trailing newline kept).
fn records(fastq: &[u8]) -> Vec<Vec<u8>> {
    let lines: Vec<&[u8]> = fastq.split(|&b| b == b'\n').collect();
    lines
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| {
            let mut v = Vec::new();
            for l in c {
                v.extend_from_slice(l);
                v.push(b'\n');
            }
            v
        })
        .collect()
}

#[test]
fn extract_range_matches_decompress_slice() {
    let mut input = Vec::new();
    for i in 0..30 {
        input.extend_from_slice(format!("@read.{i} lane\nACGTACGTAC\n+\nIIIIFFFF##\n").as_bytes());
    }
    let in_path = tmp("ex_in.fastq");
    let arc = tmp("ex.fqxv");
    fs::write(&in_path, &input).unwrap();
    run(&[
        "compress",
        in_path.to_str().unwrap(),
        "-o",
        arc.to_str().unwrap(),
        "--threads",
        "1",
    ]);

    let full = run(&["decompress", arc.to_str().unwrap(), "--threads", "1"]);
    let recs = records(&full);
    assert_eq!(recs.len(), 30);

    let arc = arc.to_str().unwrap();
    // A middle range, an open end, an open start, and a clamped end.
    for (spec, lo, hi) in [
        ("10..15", 10usize, 15usize),
        ("27..", 27, 30),
        ("..4", 0, 4),
        ("25..1000", 25, 30),
    ] {
        let got = run(&["extract", arc, "--range", spec]);
        let want: Vec<u8> = recs[lo..hi].concat();
        assert_eq!(got, want, "extract {spec}");
    }

    // Malformed ranges fail cleanly.
    run_fail(&["extract", arc, "--range", "15..10"]);
    run_fail(&["extract", arc, "--range", "nonsense"]);
}

#[test]
fn extract_rejects_reordered_archive() {
    let mut input = Vec::new();
    for i in 0..20 {
        input.extend_from_slice(
            format!("@read.{i}\nACGTTTGACCGATTGCAACGT\n+\nIIIIIIIIIIIIIIIIIIIII\n").as_bytes(),
        );
    }
    let in_path = tmp("ex_ro_in.fastq");
    let arc = tmp("ex_ro.fqxv");
    fs::write(&in_path, &input).unwrap();
    run(&[
        "compress",
        in_path.to_str().unwrap(),
        "-o",
        arc.to_str().unwrap(),
        "--reorder",
        "--keep-order",
        "--threads",
        "1",
    ]);
    run_fail(&["extract", arc.to_str().unwrap(), "--range", "0..5"]);
}
